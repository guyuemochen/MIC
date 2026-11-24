# OptaPlanner 对照程序

该程序链接 Maven Central 上的 `org.optaplanner:optaplanner-core-impl:8.44.0.Final`，使用 OptaPlanner 的规划变量、构造启发式、局部搜索和原生增量评分接口。构建需要 JDK 17 或更高版本及 Maven，依赖仅属于此基准程序。[OptaPlanner 8.44 官方手册](https://docs.optaplanner.org/8.44.0.Final/optaplanner-docs/html_single/)

从仓库根目录单独运行一例。以下命令适用于 macOS 和 Linux；正式试验由上级目录的统一脚本调度。

```sh
mvn -q -f bench/planner_comparison/opta/pom.xml compile dependency:build-classpath -Dmdep.outputFile=target/classpath
cp="bench/planner_comparison/opta/target/classes:$(cat bench/planner_comparison/opta/target/classpath)"
java -cp "$cp" benchmark.Worker bench/planner_comparison/instances/roster.txt 7 10000
```

Worker 无参数启动时会从标准输入逐行读取同样的三个参数，每次重新创建空解和求解器。每行输出的 `RESULT` 含实际返回时间、首个可行事件的时间与解、截止时的解和分数。时间单位为毫秒。外部文件解析、模型构建、SolverFactory 构建和结果序列化均在求解计时区外。

排班使用各班次的合格员工值域。排课为每节课使用独立的时段和教室变量。两个模型均按占用计数变化计算重复碰撞对。配送和生产均使用 OptaPlanner 原生列表变量。配送按改变区段的连接边和需求变化维护总路程及超载量。生产在构造阶段过滤不合格机器，局部搜索仅改变同一机器的顺序。工序计分在受影响的后继子图重算最早完工时间；成环状态的硬罚分为 1，从成环状态继续移动时重新检查整图。

各算例均从未分配状态开始。构造阶段采用框架的队列式构造和择优放置。局部搜索使用单求解线程、固定随机种子，标量变量采用实体禁忌（长度 8），列表采用值禁忌（长度 8），每步收集 32 个可接受候选。标量使用改值和交换，列表使用移位和交换。相同种子在墙钟预算下会因机器负载不同经过不同迭代数。

首解时间从调用 `solve()` 前一刻算起，在最佳解事件回调收到第一个完成初始化且硬分为零的解时记录。每次在截止前收到的最佳解均复制为独立快照。如果 `solve()` 恰好在期限内返回，也检查返回值。这样，尾部的额外运行时间不会改善归档成绩。
