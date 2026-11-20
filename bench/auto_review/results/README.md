# 自动选择器复查证据

结论：结构选择规则仍未通过验收，不合入、不默认启用。发现的二元域校验缺口已修复于主分支d49055c。
接口补强和组合测试保留在codex/auto-selector-review。详见docs/large-scale-report.md第13节。

- screen.jsonl：21例×existing/auto，5秒各1次，42条。
- paired.jsonl：5个命中模型×两模式×3次10秒，30条。
- summary.json：重复结果与旧门槛裁决。存在两次同预算失去证明的记录；只有1例满足完整速度比较条件。
- check.json：72条独立核验通过；HiGHS参考沿用上一轮，对参考已证优部分检查目标和界。
- selector-tests/selector-acceptance：测量版本完整milp及release验收日志。
- binary-bounds-before/after：新增回归测试先复现Optimal(2)错误，再确认修复后拒绝非法二元域。
- main-tests/main-acceptance/planner-tests：主分支二元域修复的完整验证。
- combined-bounds/combined-selector：实验分支同步主分支修复后的定向验证。

计时版本4692638包含内存优化和selector入口校验。主分支后来补上的Binary域校验没有回写历史测量；
所有计时输入合法，已用于规则失败裁决，未声称最终补强代码通过完整性能验收。
原始候选门槛未变，旧模型只作回归。macOS峰值RSS为逐进程字节值；编译测试与计时串行隔离。
自动降级保留LimitReached，不按提前退出时间计算加速。

复现计时版本：

```bash
git worktree add /tmp/mic-auto-review-reproduce 4692638
cargo build --release --manifest-path /tmp/mic-auto-review-reproduce/milp/Cargo.toml --example auto_config_stats
python3 bench/auto_validation/generate.py --output /tmp/mic-auto-review-cases
python3 bench/auto_validation/run.py --binary /tmp/mic-auto-review-reproduce/milp/target/release/examples/auto_config_stats --directory /tmp/mic-auto-review-cases --output /tmp/review-screen.jsonl --seconds 5
python3 bench/auto_validation/run.py --binary /tmp/mic-auto-review-reproduce/milp/target/release/examples/auto_config_stats --directory /tmp/mic-auto-review-cases --output /tmp/review-paired.jsonl --seconds 10 --repeats 3 --names packing_many_rows_s12017 packing_many_rows_s12043 packing_correlated_s12017 packing_correlated_s12043 sparse_packing
python3 bench/auto_review/summarize.py --screen /tmp/review-screen.jsonl --paired /tmp/review-paired.jsonl
python3 bench/auto_validation/check.py --directory /tmp/mic-auto-review-cases --results /tmp/review-screen.jsonl /tmp/review-paired.jsonl --reference bench/auto_validation/results/reference.jsonl
```

未补做新公共留出集、完整30秒方案或FFI/Python全套测试；没有改变ABI，生产CLI有效输入行为保持不变。
