# 自动配置实验结果（2026-11-12）

本轮未通过验收，选择器只保留在`codex/auto-config`，主分支没有`--config auto`。
详细裁决及后续工作见`docs/large-scale-report.md`第11节。首次测量前的规则在上级README，未事后修改。

## 原始数据口径

- `screen.jsonl`：最初42条筛查，旧测量工具未输出解向量，仅留作过程证据，不计入独立解检查。
- `verified-screen.jsonl`：完整21例×2模式×1次×5秒，42条，包含解向量。
- `paired.jsonl`：5个主要命中模型×4模式×3次×10秒，60条。
- `negative-paired.jsonl`：3个未命中负控和1个命中边界×2模式×3次×10秒，24条。
- `summary.json`：后两组按实例/模式统计中位数、范围、全部终态和目标、最高逐进程RSS。
- `compatibility.json`：原主分支CLI与实验CLI默认/显式existing/域外auto输出逐字节对照。
- `reference.jsonl`：外部SciPy-HiGHS独立交叉核对，15秒，threads=1、相对gap=0；不是性能对比。
- `check.json`：126条完整记录的独立解、目标、状态、界、gap与非限额重复确定性核对结果。

基线主分支5b98c70，实验实现94df771，完整解向量/四组测量工具03473a6。
模型生成规则78bdb20，manifest内SHA256校验通过。运行平台Apple M4/macOS arm64。
计时从结构选择前开始，包含选择与完整求解、不含文件读取与JSON输出，统一自适应启发式开启。
不同于墙钟达到时限，已有主动降级机制可能提前退出；保留auto_degraded，不将提前退出解释为加速。
RSS由macOS `/usr/bin/time -l`逐进程采集，单位字节；不存在并发求解或编译。
结果归档时补入已固定的time_limit_secs与threads字段，其余原值保留。
21例auto中7例命中、14例回退；小于0.1秒的模型不计百分比性能门槛。

## 复现

在独立worktree检出03473a6并构建测量程序：

```bash
git worktree add /tmp/mic-auto-reproduce 03473a6
cargo build --release --manifest-path /tmp/mic-auto-reproduce/milp/Cargo.toml --example auto_config_stats
python3 bench/auto_validation/generate.py --output /tmp/mic-auto-cases
python3 bench/auto_validation/run.py --binary /tmp/mic-auto-reproduce/milp/target/release/examples/auto_config_stats --directory /tmp/mic-auto-cases --output /tmp/auto-screen.jsonl --seconds 5
python3 bench/auto_validation/run.py --binary /tmp/mic-auto-reproduce/milp/target/release/examples/auto_config_stats --directory /tmp/mic-auto-cases --output /tmp/auto-paired.jsonl --seconds 10 --repeats 3 --names packing_many_rows_s12017 packing_many_rows_s12043 packing_correlated_s12017 packing_correlated_s12043 sparse_packing --modes existing dfs_controls hybrid auto
python3 bench/auto_validation/run.py --binary /tmp/mic-auto-reproduce/milp/target/release/examples/auto_config_stats --directory /tmp/mic-auto-cases --output /tmp/auto-negative.jsonl --seconds 10 --repeats 3 --names covering_120x200_s12109 packing_medium_s12017 negative_coefficient boundary_n200
# 外部核对脚本需独立Python环境安装SciPy；不属于求解器依赖。
python3 bench/auto_validation/reference.py --directory /tmp/mic-auto-cases --output /tmp/auto-reference.jsonl
python3 bench/auto_validation/check.py --directory /tmp/mic-auto-cases --results /tmp/auto-screen.jsonl /tmp/auto-paired.jsonl /tmp/auto-negative.jsonl --reference /tmp/auto-reference.jsonl
```

测量顺序按重复轮次翻转；同时不要编译或运行其它求解。外部参考的限额解不能当作真最优。
本轮没有新增公共留出模型，也没有完成完整30秒方案，不能宣称通过广泛性能验证。

实际外部环境SciPy1.18.1；21例中16例证优、5例限额。独立检查通过不代表所有模型都由外部证优。
实验实现完整测试及release验收日志、参数顺序与探针JSON定向测试日志一并归档。
完整测试对应94df771，两个定向测试对应03473a6；最后c7edd6f仅为实验分支补上失败裁决说明。
本轮主分支没有Rust实现改动，不重复运行planner/FFI/Python全套测试。
