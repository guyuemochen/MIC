# 无人工比值反例与宽度溢出实测

无人工Phase-2反例已确认，两个LP引擎均会在原约束被违反时误报Optimal；仅保护人工变量不足以覆盖此漏洞。

本次在主树c88c47e先执行release库构建并链接，再执行debug库构建并链接，仅运行公开LP API的微型harness。源码为 `primal_ratio_repro.rs`，逐条完整命令和库来源在各编号JSON及 `manifest.json`，原始标准输出／错误保留在03、06日志，结构化观测在 `observations.json`。没有运行全测或MILP搜索，没有修改生产代码。

模型为x,y>=0，max x+0.1*y，约束按顺序为x<=1000+d、x+0.5*y<=1000，其中d=2^-21。真实最优是1000，见证x=1000,y=0；因为目标等于(x+0.5*y)-0.4*y，所以不能超过1000。Min对照使用负目标，真实最优-1000。

危险行序的Dense/Revised × Min/Max × release/debug共8例全部返回Optimal，x=1000.0000004768372、y=0，第二条原约束残差精确为4.76837158203125e-7，目标绝对值1000.0000004768372。这超过默认feas=1e-9，也超过10*feas终检尺度。交换行序的8个对照全部返回正确最优且没有正约束残差。

另外验证x∈[-1e308,1e308]、max x／min(-x)、无约束的有限宽度溢出例。release的两引擎两方向4例均返回正确有限端点x=1e308；debug的4例均在do_flip中含Inf-Inf的断言处panic，catch_unwind捕获并保留原始panic日志。没有观察到falseUnbounded。harness退出0仅表示观测完成，不能表述为这4个LP调用通过。

完整安全比值设计一并保存为 `safe-ratio-design.md`，工作副本在 `/tmp/mic-safe-ratio-design.md`。尚未对新实现做验证，不据此宣称p0033恢复证明能力。
