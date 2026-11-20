# drayage 首解与证优耗时

MIC 和 Gurobi 10.0.3 在两例上各做三遍同机单线程求解，12 次均证优。以下为毫秒中位数，
每格来自三次完整运行。独立核验使用原模型及 MIPLIB 已归档的官方结论。

| 模型 | MIC 首解 | Gurobi 首解 | MIC 证优 | Gurobi 证优 | 最优目标 |
|---|---:|---:|---:|---:|---:|
| drayage-100-23 | 0.868 | 1.399 | 14.94 | 47.27 | 103333.874073 |
| drayage-25-23 | 0.848 | 1.368 | 87.94 | 375.80 | 101282.647018 |

首解为首次收到且事后通过原始行、界、整数性和目标核验的解向量。MIC 的首解在两例都用
37 条路线承运全部 65 项任务，目标 119252.214832。Gurobi 的首解分别承运 35 和 37 项，
目标为 574298.822022 和 543511.027279。两者最终解均承运全部任务，39 个解回调都通过核验。
证优时间从准备好输入后调用求解器起算，至求解器返回最优为止，包含预处理和回调成本。
两组均为 30 秒、200 万节点、零目标间隙。每遍轮换求解器或模型顺序。

第一遍 MIC 最优解已解码到 [optimal-routes.json](optimal-routes.json)。100-23 共有
32 条路线，25-23 有 30 条。每条路线的任务按执行顺序排列，`clock` 是原模型的时刻
变量名，`window` 和 `time` 为该变量的时间窗和解值，`inbound`、`outbound` 分别为
选中的进入、离开弧的二元变量名。解码时还用精确二进制分数复查了第一遍 MIC 的首解
和最终解，所有原始行、变量界和整数性均满足。

## 离线核验

在仓库根执行。这里只读取归档，无须 Gurobi 许可证。

```sh
rounds=bench/large_scale_rounds
run="$rounds/results/drayage-route-timing"
python3.13 "$rounds/compare_drayage_timing.py" verify "$run/timings.jsonl"
python3.13 "$run/decode_routes.py" "$run/timings.jsonl" mic > /tmp/mic-verified-routes.json
cmp /tmp/mic-verified-routes.json "$run/optimal-routes.json"
```

原始终态和每个回调的完整向量保存在 [timings.jsonl](timings.jsonl)，各遍首解、目标、
时限判定和独立核验在 [verified.json](verified.json)。[freeze.json](freeze.json) 记录
源码、探针和 Gurobi 动态库摘要、模型指纹、平台及执行命令。模型读取和导出一致性检查在
计时开始前完成。Gurobi 复用 `compare_gurobi10.py` 的参数并启用原生 MIPSOL 回调。
MIC 调用与普通求解相同的受控入口并记录接收解的时刻。[recognition.tsv](recognition.tsv)
逐例列出冻结 20 例的输入摘要和结构识别结果。

## 重新测量

冻结源码为 `3d81c18`。在独立的干净检出中构建默认 release 探针，提供冻结的两个
`.milp` 文件所在目录及本机 Gurobi 10.0.3 动态库。输出目录须尚未存在。

```sh
rounds=bench/large_scale_rounds
lib=/Library/gurobi1003/macos_universal2/lib/libgurobi100.dylib
cargo build --release --manifest-path milp/Cargo.toml --example round_probe
python3.13 "$rounds/compare_drayage_timing.py" run --library "$lib" \
  --mic milp/target/release/examples/round_probe --instances /tmp/mic-rounds-miplib2017 \
  --output /tmp/mic-drayage-new-run
python3.13 "$rounds/compare_drayage_timing.py" verify /tmp/mic-drayage-new-run/timings.jsonl
```

普通 CLI 也可直接求解，例如
`cargo run --release --manifest-path milp/Cargo.toml -- /tmp/mic-rounds-miplib2017/drayage-25-23.milp --threads 1 --time-limit-secs 30 --early-stop off`。
普通 CLI 无须 Gurobi 许可证。
