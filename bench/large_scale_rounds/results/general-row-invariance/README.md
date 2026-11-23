# 等价行变换与通用预处理的独立核验

[核验程序](check.py) 在临时目录生成 35 个模型，依次从 MIC 的 `timed-solve` 入口和 Gurobi 10.0.3 求解。同一模型使用单线程、12 秒、200 万节点、零目标间隙。模型读取完成后开始计时，每个程序对每个模型运行一次。单次耗时用于观察实际开销，线性松弛迭代数、搜索节点数、删行数用于核对等价写法是否进入相同路径。MIC 的专用整数搜索与通用分支定界的节点定义不同。

- 六个覆盖写法来自两个种子，分别保留原行、只将第一行乘以 2、或逐行循环使用 3、−0.5、0.25、−5、2。另一例把缩放后的二元变量声明为上下界 0 和 1 的普通整数变量。程序重建了上轮发现的种子 817391 反例，另外一组使用新种子 905731。
- 两个覆盖反例只改一行。将 `x + y >= 1` 改成 `2x + 2y >= 2.5` 会要求两列均取一，最优费用由 1 变为 3。
- 三个新种子的稠密装载模型各有 194 列、9 行。每例分别检查原行和带正负与分数倍率的等价写法。容量分别取各行总重量的 1/28、1/28 和 1/6，以检验不同紧度下同一算法的获益与开销。
- 19 个两变量模型包含连续、混合、纯整数、等式及冲突约束。除了核对 Gurobi，检查器还用精确分数枚举连续多边形的顶点、整数切片的端点或全部整数点。近似平行反例的两行系数只相差 2⁻²¹，最优值为 16384，误删一行会变成 32768。
- 80 列的连续模型从 100 条原行派生 600 条约束，含相反方向的等价写法和右端较松的约束。3 个模型另经双线程入口求解，4 个模型另以单节点预算复核中断时的解和界。

另外对 140 个小模型测试了从 2⁻²⁸ 到 2²⁸ 的行倍率，并混用负方向、负系数与松紧不同的右端。两版 MIC 的每个终态和界分别与精确有理数枚举对账，原始模型序列的摘要记在 `exact_scale_sweep`。这一组的参照是独立枚举。

每个程序的输入都从同一 `.milp` 文本重建，MIC 的实际模型导出和 Gurobi 读回的变量、行、目标需与输入相同。验证器对所有终解和求解中回调的向量逐行验可行性并重算原目标。全二元模型再对取整后的每个向量作整数精确复核。所有已证优的目标与界须在 1e−6 内符合独立参照，所有触限界不得排除参照最优解。[逐例结果及程序指纹](checks.json)

在仓库根目录重放，需要 Python 3.10+、两个发布版 `round_probe` 和 Gurobi 10.0.3 的本地动态库。脚本的 `old` 对应 `0b47f2a`。下例的临时目录会保存模型及包含向量的原始运行记录。

```sh
work=$(mktemp -d)
git worktree add --detach "$work/old-src" 0b47f2a
cargo build --locked --release --manifest-path "$work/old-src/milp/Cargo.toml" --example round_probe --target-dir "$work/old-target"
cargo build --locked --release --manifest-path milp/Cargo.toml --example round_probe --target-dir "$work/new-target"
script=bench/large_scale_rounds/results/general-row-invariance/check.py
python3 "$script" old --binary "$work/old-target/release/examples/round_probe" --scratch "$work/runs"
python3 "$script" gurobi --binary /Library/gurobi1003/macos_universal2/lib/libgurobi100.dylib --scratch "$work/runs"
python3 "$script" new --binary "$work/new-target/release/examples/round_probe" --scratch "$work/runs"
python3 "$script" summarize --scratch "$work/runs" --output "$work/checks.json"
```

## 冻结百例回归

[百例逐项结果](corpus.json) 保存两个二进制各自的 2 秒短测、当前版 14 个触限模型从头以 60 秒预算重跑的结果，以及两个连续模型的三次短预算交错重测。每次求解均为单线程、200 万节点、零目标间隙。60 秒扩展按三组并行处理独立模型，计时仅用于确认预算。合并两档的已证明结论后有 88 例证优、2 例证不可行、1 例证无界，其余 9 例触限。程序逐个核对冻结模型摘要、最终向量及回调，并用已归档的有效解或已知结论核对界。下面的命令使用上一节构建的程序。

```sh
corpus=bench/large_scale_rounds/results/general-row-invariance/corpus.py
python3 "$corpus" "$work/old-target/release/examples/round_probe" "$work/old-2s.json"
python3 "$corpus" "$work/new-target/release/examples/round_probe" "$work/new-2s.json"
# 单独复查任意模型时指定其名字和时限。
python3 "$corpus" "$work/new-target/release/examples/round_probe" "$work/example-60s.json" --seconds 60 --instances s1234 circ10-3
```

对通用预处理另做了平行约束删行的隔离试验。100 个模型均通过相同的缩放和预处理入口，4 例有额外删行；在 `neos-4954672-berkel` 上，两种二进制各跑一次 60 秒，最终可行解及中途回调均通过原约束复算。[模型扫描与单项试验](parallel-rows.json) · [含向量的原始数据](parallel-runs.jsonl.gz)
