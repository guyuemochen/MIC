# 混合父域内存优化结果

本轮通过事前门槛，求解器改动已合入6f470a1。对照使用固定工作量，不宣称完成求解的加速比。
三个200000节点模型的RSS中位下降84.2%—84.5%，耗时几何平均比0.9650。

- `regression.jsonl`：8例×3配置×新旧版本，48条，50000节点。
- `memory.jsonl`：3例×新旧版本×3次，18条，200000节点，轮次交替顺序。
- `summary.json`：全部确定性字段相等，逐例耗时中位/范围、RSS中位及预注册门槛裁决。
- `check.json`：66条解的独立检查，使用上一轮HiGHS参考；参考限额结果不当作真最优。
- `metadata.json`：提交、二进制与探针SHA256、实际开关和平台。
- 三份测试日志：候选完整milp、release验收与planner；既有忽略原样保留。

复现（在仓库根执行，需macOS的`/usr/bin/time -l`取得本轮同口径RSS）：

```bash
git worktree add /tmp/mic-domain-old 22efddd
git worktree add /tmp/mic-domain-new 6f470a1
cp milp/examples/domain_stats.rs /tmp/mic-domain-old/milp/examples/domain_stats.rs
cargo build --release --manifest-path /tmp/mic-domain-old/milp/Cargo.toml --example domain_stats
cargo build --release --manifest-path /tmp/mic-domain-new/milp/Cargo.toml --example domain_stats
python3 bench/auto_validation/generate.py --output /tmp/mic-domain-cases
python3 bench/domain_memory/run.py --baseline /tmp/mic-domain-old/milp/target/release/examples/domain_stats --candidate /tmp/mic-domain-new/milp/target/release/examples/domain_stats --directory /tmp/mic-domain-cases --output /tmp/domain-regression.jsonl --plan regression
python3 bench/domain_memory/run.py --baseline /tmp/mic-domain-old/milp/target/release/examples/domain_stats --candidate /tmp/mic-domain-new/milp/target/release/examples/domain_stats --directory /tmp/mic-domain-cases --output /tmp/domain-memory.jsonl --plan memory
python3 bench/domain_memory/check.py /tmp/domain-regression.jsonl /tmp/domain-memory.jsonl
python3 bench/auto_validation/check.py --directory /tmp/mic-domain-cases --results /tmp/domain-regression.jsonl /tmp/domain-memory.jsonl --reference bench/auto_validation/results/reference.jsonl
```

全程无墙钟截断、不主动降级，探针保留默认启发式；每次子进程另有180秒看门狗，
若失败会记录driver_error并使核验失败。计时期间不要编译或运行其它求解任务。
本轮复用了旧验证模型，属于回归，不能据此宣称自动配置选择通过独立验收。
