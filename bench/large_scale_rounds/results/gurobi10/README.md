# Gurobi 10.0.3 冻结 20 例外部对照

本机官方 Gurobi 10.0.3 C 动态库完成全部 20 次运行。每例单线程、30 秒、2000000 节点，
两种 gap 均设为 0。设置的整数容差为 `1e-6`，对偶预处理及其他可行和最优性容差沿用
Gurobi 默认。2 行列和 2001 行列的授权探针均完成，实际 20 例均获许可。工具冻结于 `855127c`。

Gurobi 证优 15 例，5 例返回 `TIME_LIMIT`。12 个旧模型全部证优。现代 8 例中证优的是
`neos17`、`drayage-100-23` 和 `drayage-25-23`。其余 `bnatt500`、`seymour`、
`neos-4954672-berkel`、`gen-ip054`、`icir97_tension` 均未在预算内完成证明。
Gurobi 的 `bnatt500` 记录有有限界 0、无可行解。已归档的官方结论为不可行。

| 模型 | Gurobi 证优目标 | 单次优化秒数 | R12 同模型证明 |
| --- | ---: | ---: | --- |
| `p0201` | 7615 | 0.0847 | 基础 0/3，候选 3/3 |
| `neos17` | 0.1500025774 | 0.4614 | 基础 0/3，候选 0/3 |
| `drayage-100-23` | 103333.874073 | 0.0467 | 基础 0/3，候选 0/3 |
| `drayage-25-23` | 101282.647018 | 0.3805 | 基础 0/3，候选 0/3 |

`independent-check.json` 核验了 272 条记录，涵盖归档的 HiGHS、SCIP 基线和同机 R12 正式配对。
233 个可行解和 126 条有外部最优证据的 MIC 界通过复核。20 个模型的 Gurobi 载入读回检查
均通过，原约束、目标、常数、列和界与归档数据核对。`comparison.json` 保存逐模型的同预算
R12 配对，并列出终态、目标、界、节点、迭代、割、耗时和峰值内存。求解器的节点与迭代
记账口径各自定义。耗时基于单次外部观测，峰值内存覆盖整个外部进程，包括模型准备和读回。

`original/` 保存原始结果、命令、许可探针及动态库和输入摘要。归档版 `gurobi10.jsonl` 只改
模型路径，引用 `../models/`。在含 R12 归档的仓库根目录复算无须运行求解器。

```sh
rounds=bench/large_scale_rounds
refs="$rounds/results"
out=$(mktemp -d)
python3.13 "$rounds/check_results.py" "$refs/external-baseline.jsonl" \
  "$refs/miplib2017-baseline.jsonl" "$refs/root-cut-recheck/paired/root_cut_recheck-paired.jsonl" \
  "$refs/gurobi10/gurobi10.jsonl" --public-truth "$refs/miplib2017-sources.json" \
  --output "$out/check.json"
```
