# 加权集合覆盖的独立核验

在仓库根目录构建发布版探针，再指定本机安装的 Gurobi 10.0.3 动态库。脚本只依赖 Python 3.11 标准库和仓库已有的输入解析器与 Gurobi 调用代码。每次求解都使用单线程、60 秒、零目标间隙，普通运行的节点上限为 200 万。

```sh
cargo build --locked --release --manifest-path milp/Cargo.toml --example round_probe
python3 bench/large_scale_rounds/results/weighted-setcover/validation/check.py \
  milp/target/release/examples/round_probe \
  --gurobi /Library/gurobi1003/macos_universal2/lib/libgurobi100.dylib \
  --output /tmp/weighted-setcover-validation.json
```

- 原题和新种子模型都从 MILP 文本读出费用与约束。对 MIC 和 Gurobi 的每个可行解回调和最终解，逐项检查所有变量在距 0 或 1 的 `1e-6` 范围内，再用整数计算每一行的覆盖和目标。结果保存每次最终选中的原列号（从 0 开始）和已核验的回调费用。
- 原题额外打乱行序、列序和列名，把交替的行整体乘以负一并翻转不等号。另一变体将所有费用乘七，新增一条只能由费用为 19 的新列覆盖的需求。两者的解都映射回原题复算。
- 每个小问题用十二列构造 32 条需求，再补到 256 列。新增列只覆盖某个原列的部分或全部需求，费用至少等于该原列；任一使用新增列的解都可用费用不增的原列替换。因此十二列的全部 4096 个子集足以给出各模型的精确最优值。十二个固定种子覆盖了同费重复列、受支配列及最高一百万的单列费用。
- 原题额外在一个节点、完整证优所需节点数减一、恰好达到完整节点数时运行。各小模型也检查一个节点的预算。触限时保留的解和下界都与各自的独立最优值比较。
- 两份新模型重用仓库生成器的抽样规则，规模均为 400 行、800 列，每项以 0.08 的概率入行，列费用从 1 至 100 均匀抽取。种子预先固定为 7341、7342。每份输入分别由 MIC 与 Gurobi 证优。
