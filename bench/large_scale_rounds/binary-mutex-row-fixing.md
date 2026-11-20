# 二元互斥结合短行固定的定点试验

候选特性 `binary_mutex_row_fixing` 默认关闭。先验证预处理能否用已认证的二元
互斥关系和同一缩减模型的变量界，在短行上排除某个二元变量取 1。两组从同一干净
提交构建，继承相同的生产默认，候选组只增加本特性。

正式定点只跑 `drayage-100-23` 和 `drayage-25-23`，两例各组各一次，按原固定
矩阵中这两例的相对顺序与组别轮换执行。沿用第七轮参数：单线程、2000000 节点、
30 秒、二元传播自动、关闭探测复用与提前停止。正式计时关闭诊断，记录原模型哈希、
构建特性、终态、可行目标、对偶界、节点、LP 迭代和峰值 RSS。

用原模型独立检验可行解、重新计算目标并核对界方向和参考答案。两例都须零正确性
错误、保留基础组证明，每次候选 RSS 至多 max(64 MiB, 2 倍同次基础组)。
基础组已取得的有限可行目标和对偶界，在原矩阵的数值核对容差内都须保留或改善。
若两组共同证优，则每例耗时比须不超过 1.25。在满足以上条件时，候选至少增加一例
完成的证明或新增一例有可行解，或使两例均出现新的有限对偶界或严格改善已有界，
定点才算通过。诊断可另开，只用于解释实际固定的数量和发生位置。

这两例来自同一模型家族。定点通过后候选继续保持默认关闭；若要推广，另按固定 20
例各一次 10 秒筛查及各三次 30 秒正式配对验证，门槛在开始计时前冻结。

## 活动度证据

在当前工作模型的变量界上检查短行。先复用静态二元行的精确格点认证，确认一组变量
中任意两列不能同时为 1。对于同组的两个自由二元列，依次假设一列为 1、另一列为 0，
使用外向舍入计算剩余各项可能达到的活动度区间。只有区间完全越过带容差的约束边界
才固定假设取 1 的列为 0。余量在当前工作行尺度上覆盖单纯形、取整复核和终检形式的
误差，并给贴界留出余量。完整求解后仍在原模型空间独立复核。

在干净源码树根目录构建两组库与正式探针。核查和定点脚本使用 Python 3.11 或更新版本。
以下命令将产物写到临时目录。正式定点的
输入使用现有的 `miplib2017-baseline.jsonl` 与官方来源清单，其中归档模型及原
`.milp` 路径须可读且哈希吻合。

```sh
audit=$(mktemp -d)
source=$(pwd)
rounds="$source/bench/large_scale_rounds"
for variant in base candidate; do
  set --
  if [ "$variant" = candidate ]; then
    set -- --features binary_mutex_row_fixing
  fi
  target="$audit/$variant-target"
  cargo build --manifest-path "$source/milp/Cargo.toml" --release --lib \
    --example round_probe --target-dir "$target" "$@"
  cp "$target/release/examples/round_probe" "$audit/round_probe_$variant"
  rustc --edition 2021 -O --extern "milp=$target/release/libmilp.rlib" \
    "$rounds/binary_mutex_presolve_probe.rs" -o "$audit/presolve_$variant"
done
cat > "$audit/binaries.json" <<'JSON'
{"base":{"binary":"round_probe_base","features":[]},"candidate":{"binary":"round_probe_candidate","features":["binary_mutex_row_fixing"]}}
JSON
```

将下面的 `model_dir` 设为两个原始 `.milp` 的目录。探针只运行生产缩放与预处理，
输出原列固定集。独立 Python 核查读旧版缩减模型；使用精确有理活动度复核短行的
固定集，另外模拟历史逐行盒传播到不动点。逐行模拟是机制对照，当前生产没有启用该
历史规则。

```sh
for case in drayage-100-23 drayage-25-23; do
  for variant in base candidate; do
    "$audit/presolve_$variant" "$model_dir/$case.milp" "$audit/$case-$variant.milp" \
      > "$audit/$case-$variant.tsv"
  done
  python3 "$rounds/binary_mutex_drayage_audit.py" "$model_dir/$case.milp" \
    "$audit/$case-base.milp" "$audit/$case-base.tsv" "$audit/$case-candidate.tsv"
done
python3.13 "$rounds/run_binary_mutex_point.py" --source "$source" \
  --binaries "$audit/binaries.json" \
  --reference "$rounds/results/miplib2017-baseline.jsonl" \
  --public-truth "$rounds/results/miplib2017-sources.json" \
  --output "$audit/point.jsonl" --gate "$audit/gate.json"
```

## 定点结果

两例的 30 秒对偶下界均严格改善，按冻结判定通过定点。四次求解均在根节点诚实
限时，尚无整数可行解。两例来自同一家族，特性继续默认关闭。原始字段、固定集、
独立检查与复现说明见 [定点归档](results/binary-mutex-row-fixing/README.md)。
