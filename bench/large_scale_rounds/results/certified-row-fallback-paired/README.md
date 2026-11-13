# 第十一轮正式配对归档

正式运行使用筛查原二进制和冻结源码 `3d136d6`。固定 20 例，单线程，每组 30 秒各三次。
共 120 条记录。执行、独立复核、汇总、推广门禁的退出码分别为 0、0、0、1。

两组各有 27 次证优、33 次 LimitReached（含数值失败映射），均有 9 例三次证优。无证明丢失、无新增稳定证明、
无 RSS 或单例耗时超标。8 例进入几何平均，耗时比为 0.9270517832354401。正确性复核
覆盖 252 条记录、211 个可行解、126 条有外部最优值可对照的界，无错误。唯一未过的是性能
门槛，保持默认关闭。

`original/` 内是未经改写的执行文件、原始 JSONL 和原报告。`archive-relocation.json`
登记原始和归档 JSONL 的 SHA-256，及独立重算的结果。仓库版 JSONL 只改了模型路径，
引用已归档的 `../models/`。重算的独立复核和汇总与原件一致，门禁仅输入路径不同。

从仓库根目录重算。仅读取归档，不启动求解器。

```sh
rounds=bench/large_scale_rounds
refs="$rounds/results"
archive="$refs/certified-row-fallback-paired"
raw="$archive/certified_row_fallback-paired.jsonl"
replay=$(mktemp -d)
python3.13 "$rounds/check_results.py" "$refs/external-baseline.jsonl" \
  "$refs/miplib2017-baseline.jsonl" "$raw" --public-truth "$refs/miplib2017-sources.json" \
  --output "$replay/check.json"
python3.13 "$rounds/summarize_feature_trial.py" "$raw" \
  --legacy-baseline "$refs/external-baseline.jsonl" \
  --modern-baseline "$refs/miplib2017-baseline.jsonl" \
  --public-truth "$refs/miplib2017-sources.json" --output "$replay/summary.json"
python3.13 "$rounds/gate_feature_trial.py" "$replay/summary.json" \
  --output "$replay/gates.json"
```

最后一条命令预期返回 1，并保存未达性能门槛的结果。
