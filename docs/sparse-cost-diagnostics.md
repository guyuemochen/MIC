# 稀疏 LP 成本诊断与列缓冲复用实验

本轮不修改定价、比值测试、反循环或数值容差。历史 S8 混合定价已产生负收益，本轮不重试该方向。

## 诊断口径

`sparse_profile` 是默认关闭的 Cargo 特性；关闭时所有计时入口均不参与编译。
开启后按当前线程累计互斥时间，嵌套调用暂停父类别，不重复加总：

- `refactor`：基列收集、排序、LU 分解及基本解重算。
- `pricing`：定价，包括检验数与候选扫描；嵌套基求解另计。
- `basis_solve`：FTRAN/BTRAN 及其 ETA 更新作用。
- `basis_update`：主元更新、ETA 追加、Devex 权重维护；嵌套重构和基求解另计。

四类不覆盖全部求解工作，不能把其和当作端到端耗时。类别调用次数包含同类嵌套入口。
统计不随求解器快照回滚，且不汇总其它线程。诊断时间只供观察，绝不参与求解决策。

```bash
cargo run --release --manifest-path milp/Cargo.toml --features sparse_profile --example sparse_cost -- bench/instances/lp-sparse-4000x800.milp
cargo run --release --manifest-path milp/Cargo.toml --features sparse_profile,sparse_reuse_columns --example sparse_cost -- bench/instances/lp-sparse-4000x800.milp
```

工具直接求解原始模型的 LP 松弛，不运行预处理、割或分支定界。可在路径后加单次 LP 迭代上限；限额终态照实输出。
正式 A/B 耗时应关闭 `sparse_profile`，只切换 `sparse_reuse_columns`，并对比目标、状态和迭代数。
模型构造计入工具总时间，文件加载不计入。

## 候选优化与边界

代码证据：`refactor_cols.clear()` 会销毁所有内层向量，紧接着 `resize_with` 重新创建，
与原有“复用容量”注释不符。实验特性 `sparse_reuse_columns` 跳过外层 clear，保留各位置容量；
随后每列仍 clear 并按原顺序收集，LU 输入及运算顺序保持不变。行数减小时照常释放尾部列。

风险是保留较大历史列容量导致内存上升；预期只减少分配，不能预先声称有明显总时间收益。
该特性默认关闭，须由同机串行基准决定是否保留或转为默认。
