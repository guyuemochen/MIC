# AGENTS.md

## 沟通规范（最重要）

**回答的时候尽量说人话：**

- 先用一句大白话给结论，再讲细节。别让读者翻三段才看到答案。
- 少用黑话和缩写。术语第一次出现时，用一句话解释它是什么（比如"BTRAN 就是把转置矩阵作用到向量上，用来算对偶价格"）。
- 能用一段自然语言说清楚的，就不要铺表格、列表和层层标题。列表只留给真正需要逐条对照的东西。
- 不确定就说不确定，别用含糊的措辞掩盖。报结果要如实：测试挂了就说挂了，附上输出；跳过了哪步就说跳过了。
- 解释代码时讲"为什么这么写、约束是什么"，不要复述一遍代码本身。

## 项目是什么

MIC（Mixed-Integer & Constraint Solvers）是一组互相独立的 Rust 求解器 crate，纯标准库、零第三方依赖、无 `unsafe`：

- `milp/` — MILP 精确求解器：两阶段有界单纯形（稠密 `simplex.rs` + 稀疏 `simplex_sparse.rs`）、DFS 分支定界、Gomory/聚合 MIR 割、presolve 与启发式、多线程。
- `planner/` — OptaPlanner 式启发式约束求解器（构造 + 禁忌/SA/迟到接受局部搜索），与 milp 完全不共享代码。
- `milp-ffi/` — C ABI cdylib 绑定（panic 不跨 FFI 边界）。
- `micpython/` — gurobipy 风格的 Python 包（ctypes，零依赖，Py3.8+）。

当前主要工作在 `feat/gurobi-large-scale` 分支（大规模场景优化），逐轮实测记录见 `docs/large-scale-report.md`。

## 常用命令

注意：仓库根**没有** Cargo.toml，一律用 `--manifest-path`。

```bash
# 构建 + 运行
cargo run --manifest-path milp/Cargo.toml --release -- milp/instances/knapsack-10.milp
cargo run --manifest-path planner/Cargo.toml --release -- planner/instances/roster-small.txt --steps 2000 --seed 7 --algo tabu --construct bf

# 测试（milp 和 planner 都要跑；验收测试必须 --release）
cargo test --manifest-path milp/Cargo.toml
cargo test --manifest-path planner/Cargo.toml
cargo test --manifest-path milp/Cargo.toml --release --test acceptance

# Python 侧（先出 DLL；仓库根即包根，不需要 pip install）
cargo build --release --manifest-path milp-ffi/Cargo.toml
python micpython/examples/knapsack.py
python -m unittest discover -s micpython/tests -t .

# 基准（会自动先构建 profile_stats）
python bench/run_timings.py --set gate --threads 1 --json bench/gate_runs.json
```

## 代码铁律（违反即回退，不商量）

- **零第三方依赖、无 `unsafe`**。所有 crate 都是纯标准库。
- **确定性**：`threads=1` 时结果完全确定，并列时一律取最小下标破键。改求解逻辑时不得破坏这一点。
- **诚实性**：只有证到最优才返回 `Optimal`；触发任何上限（时间/节点/迭代）一律 `LimitReached`。绝不把启发式解谎报成最优。
- **FFI 边界**：panic 不跨边界，错误粘滞（出错后句柄保持错误态）。ABI 版本恒为 1。
- **稀疏/稠密两套单纯形是同一 LP 语义的两个实现**：改任何一个的比值测试、反循环、两阶段纪律，都要检查另一个是否需要对齐。契约逐条写在 `milp/src/simplex_sparse.rs` 模块头注释里（前 60 行就是架构文档，改代码前先读）。

## 约定

- 注释、文档、commit 信息全部用**中文**。
- commit 前缀：`milp S2:`（大规模轮次）、`bench:`、`docs:` 这类"模块 + 轮次/主题"格式。
- 测试框架就是 Rust 内置 `#[test]` 和 Python `unittest`，不要引入别的测试依赖。
- 改完代码至少跑对应 crate 的 `cargo test`；动到求解核心时补跑 `--release --test acceptance`（debug 下验收测试会超时或数值不稳定）。

## 别踩的坑

- 根目录的 `nul` 文件是 Windows 设备名重定向的残留，已 gitignore，**勿删勿动**（在 Windows 上很多工具打不开它）。
- `.zcodeignore` 顶部注释之后的部分才是自定义区，上半部分与 `.gitignore` 同步维护，别乱动。
- `micpython/pyproject.toml` 只是元数据声明，安装路径不纳入门禁，别围绕它做发布逻辑。
- `bench/instances/` 里背包算例的 `.opt` 文件是 DP 算出的真值，基准对账靠它们，别覆盖。
