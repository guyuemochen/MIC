# 单线程自动选择线性引擎的对照测试

基线是 `e765d6a`。两版 MIC 都从默认求解入口运行，模型加载完成后开始计时，关闭自动提前停止。Gurobi 10.0.3 使用相同的冻结模型、单线程、200 万节点和零目标间隙。

[冻结百例的对照程序](compare.py) 先核验输入及导出模型的摘要，再比对两者的变量、目标与约束。每次运行都按原模型核验终解及每次回调，也核对可行目标和证明界。两版 MIC 按样例和轮次交替先后顺序。2 秒短测接近时间边界，代表样例另给 60 秒预算，重复计时取中位数。

[百例逐项对照](corpus.md) 列出最终短测每个模型的耗时、状态、可行目标、证明界及迭代数。[实测记录](timings.json.gz) 包含两轮独立的百例短测（`initial_corpus`、`corpus`）、九个样例三轮 MIC 和 Gurobi 配对（`paired`）、周期等式模型的一轮 60 秒对照（`periodic`）及多线程、单节点和计时波动复核。每组内保存程序 SHA-256、模型名、求解状态、解、界和工作量；解字段是已经按原变量向量校验的目标值。

[独立模型生成与校验](holdouts.py) 构造 11 个匹配、网络流、容量和作业冲突模型。匹配与流模型有已知可行解及同值的费用下界，容量模型有可达的加权约束上界，作业冲突使用独立动态规划求最优值。每个模型各运行三次，时限 12 秒。稠密对照使用每行 60 项的 240 列容量模型。归档分别保存[旧版默认](holdouts-baseline-auto.json)、[新版默认](holdouts-auto.json)和[显式选择两种引擎](holdouts.json)的结果。

另在同一 240 列、240 行的生成模型上分别使用每行 14、15、16 个非零项，覆盖密度门槛的两侧及边界。两版默认和显式引擎各跑三遍，54 次都证明了独立计算的最优值。`timings.json.gz` 内的 `density_old`、`density_new` 保存逐次结果。

[未纳入默认行为的试验](experiments.json.gz) 包含 `air03` 的原例、行列重排、缩放测试和周期搜索的费用、界及工作量。`within_batch_duplicate_cut` 的六个随机种子依次用 Python 的同一个随机数生成器重排列、重排行并按新列序更名。`periodic` 的 24 个生成模型经过单节点与正常预算测试。原例和所有等价模型的可行解都经过复核。

在仓库根目录可以运行以下命令。Gurobi 对照需要可用的本机许可，`GUROBI_LIBRARY` 指向其动态库。

```sh
work=$(mktemp -d)
git worktree add --detach "$work/baseline" e765d6a
cargo build --locked --release --manifest-path "$work/baseline/milp/Cargo.toml" --example round_probe
cargo build --locked --release --manifest-path milp/Cargo.toml --example round_probe
old="$work/baseline/milp/target/release/examples/round_probe"
new=milp/target/release/examples/round_probe
checks=bench/large_scale_rounds/results/auto-lp-engine
python3 "$checks/compare.py" "$old" "$new" "$work/corpus.json" --seconds 2
python3 "$checks/compare.py" "$old" "$new" "$work/paired.json" --seconds 60 --reps 3 \
  --case neos17 assign-40x40 lseu air03 setcover-400x800 p0201 lp-sparse-4000x800 p0282 p0548 \
  --gurobi-library "$GUROBI_LIBRARY"
python3 "$checks/compare.py" "$old" "$new" "$work/periodic.json" --seconds 60 \
  --case icir97_tension --gurobi-library "$GUROBI_LIBRARY"
python3 "$checks/holdouts.py" "$new" "$work/holdouts.json" --engines auto --reps 3
python3 "$checks/holdouts.py" "$new" "$work/density.json" --density-boundary --engines auto tableau revised
```
