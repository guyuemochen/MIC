# 用 Python 和 Java 建模

一次建模可以同时使用分配变量和有序列表。排课时，用两个分配变量分别表示一节课的时段和教室。配送时，用一个列表组表示所有车辆，一名客户在完整结果中恰好属于一辆车，列表内的位置就是访问顺序。

## 从仓库根目录运行

Python 3.8 及以上版本用标准库加载 Rust 动态库。构建和运行四个示例只需要两条命令。

```sh
cargo build --release --manifest-path planner-ffi/Cargo.toml
python3 micpython/examples/planner_cases.py
```

Java 8 及以上版本使用 JDK 自带的本地调用接口（JNI）。准备好 JDK、Rust 和 C 编译器后执行以下命令。构建脚本也可使用 `JAVA_HOME` 指定 JDK。Windows 的 MSVC 工具链从开发者命令行运行。

```sh
python3 planner-java/build.py
javac -encoding UTF-8 -cp planner-java/target/micplanner.jar -d planner-java/target planner-java/examples/PlannerCases.java
java -cp "planner-java/target:planner-java/target/micplanner.jar" PlannerCases
```

Windows 使用分号分隔 Java 类路径。部署时将 JAR、同目录的 JNI 库和 Rust 动态库一起复制。Python 的库路径可以用 `MIC_PLANNER_LIB` 指定；Java 的 JNI 库路径可以用 `MIC_PLANNER_JNI` 指定。

## 声明业务规则

下面给两个重叠班次选员工，要求两班由不同的人承担，并尽量满足第二个班次的人员偏好。

```python
from micpython.planner import Model

with Model("排班") as m:
    morning = m.addVar([0, 1])
    evening = m.addVar([0, 1])
    m.addConstr(morning != evening, "班次冲突")
    m.addConstr(evening == 1, "晚班偏好", weight=2, hard=False)
    result = m.optimize(seed=7, steps=200)
    if result.feasible:
        print(result.value(morning), result.value(evening), result.penalties)
```

Java 的模型和表达式有相同的含义，比较和运算使用方法调用。

```java
import org.micplanner.Model;

// 以下代码放在调用方的方法内
try (Model m = new Model()) {
    Model.Expr morning = m.addVar(0, 1);
    Model.Expr evening = m.addVar(0, 1);
    m.addConstr(morning.ne(evening), "班次冲突");
    m.addConstr(evening.eq(1), "晚班偏好", 2, false);
    Model.Result r = m.optimize(7, 200, -1, 1, "tabu", "ff");
    if (r.feasible()) System.out.println(r.value(morning) + ", " + r.value(evening));
}
```

- **规则和目标**。每条规则都有唯一名称。`addConstr` 在条件不成立时计一次罚分；`addPenalty` 把非负违反量乘以权重计入罚分。`setObjective` 指定要最小化的整数表达式。硬分为所有硬罚分之和，软分为目标值加所有软罚分，硬分优先。规则、目标和分数使用有符号 64 位整数，溢出会报告错误。
- **分组与查表**。用 `groupCount(keys, key)` 统计指定组的成员，用 `groupSum(keys, values, key)` 计算该组的合计。用 `count` 按条件计数，`sum` 或 `max` 聚合表达式，`choose` 按条件选取数值，`lookup` 按从零开始的整数索引查业务表。员工、教室和时段可通过业务代码映射为整数编号。
- **路线与生产**。`addListVar(owners, eligible, durations)` 声明所有者数量、每个元素允许的所有者及可选工时。列表的 `load` 和 `distance` 分别计算载重和有向往返总距离，距离表的最后一行和一列表示仓库。`adjacent(owner, before, after)` 可直接声明同一列表中两元素的有向相邻关系，例如禁止某两道工序连续加工。`precede(a,b)` 规定工艺先后关系，`finish(i)` 计算同时满足工艺和机器队列后的最早完工时刻。模型拒绝互相矛盾的工艺依赖，搜索会跳过引入循环等待的换序。结果通过 `sequence(group, owner)` 读取。
- **预算与取解**。默认种子 42、50,000 个局部搜索步、单线程、首次适配构造和禁忌搜索。算法可选 `tabu`、`sa`、`la`，构造可选 `ff`、`bf`。时间预算以毫秒计，包含构造和搜索，默认 -1 表示只限制步数。零步仍会构造初始解。达到上限时 `status` 为 `LimitReached`，`reason` 指明 `Steps` 或 `Time`。构造中触及时限可能返回部分方案，所以先检查 `complete` 和可行性。没有可移动的变量时可以提前返回 `Completed/NoMoves`。这些状态描述启发式搜索的停止原因。
- **重复调用**。固定模型、种子、线程数和步数，可在同一二进制上复现结果。每次求解都从构造阶段开始，返回独立快照；可以改模型后再次求解。原生求解报错后关闭该模型，使用错误信息修正输入，再创建新模型。不同业务模型可并行运行。

[Python 四个完整示例](../micpython/examples/planner_cases.py) · [Java 四个完整示例](../planner-java/examples/PlannerCases.java) · [Rust 模型定义](examples/business_models/cases.rs)

## 验证三种语言一致

```sh
cargo test --manifest-path planner/Cargo.toml --release --test business_acceptance
python3 planner/tests/check_bindings.py
```

第一条用穷举核对微型配送和生产排程的最佳值，并验证三种搜索策略。第二条重建动态库，让 Python 和 Java 加载同一份本地产物，逐项对比四个模型的分配、序列、罚分和停止原因。

较大算例另对比了同一时间预算下 MIC 与 OptaPlanner 的首个可行解和最终代价。[四类业务等时实测](../bench/planner_comparison/results/2026-11-24/report.md)

最佳适配构造（`bf`）会从硬约束中识别可以共同选值的变量对。例如一节课用「时段 × 教室数 + 教室」表示占用位置，构造时就能同时考虑时段偏好和教室冲突。变量对互不重叠且候选组合至多 1024 个时启用。若一组约束比较了组内所有配对，评分器会记录每个取值有多少项占用，移动时只更新涉及的计数。[排课优化与复测](../bench/planner_comparison/rounds/paired-scoring/README.md)

禁忌搜索会在最好可行解附近检查整对互换和重新选值，例如让两门课直接互换时段和教室。离开最好解探索一段时间后，搜索也可以回到已保存的方案再次尝试。[局部搜索诊断与排课复测](../bench/planner_comparison/rounds/local-search/README.md)

## 后续开发方向

业务方已有禁忌搜索（`tabu`）、模拟退火（`sa`）和迟到接受（`la`）可选。后续希望让 Python 和 Java 用户进一步控制以下三件事。

- **选择或组合移动算子**。按业务指定一次允许改哪些决定，例如整节课同时换时段和教室、两辆车交换一段路线。移动的预演、提交和撤销共用 Rust 评分器。
- **配置搜索阶段**。把候选移动、接受策略和每步候选数量分别配置，让同一模型可以比较禁忌、退火等策略，或分阶段切换策略。
- **扩展群体搜索**。给自定义阶段提供维护多份候选解的接口，再评估鸟群、粒子群等算法。离散的员工、教室和路线需要定义候选解如何组合以及怎样修复冲突。新增策略沿用同一时间预算、首解记录和独立验分方法，与内置策略做对照。
