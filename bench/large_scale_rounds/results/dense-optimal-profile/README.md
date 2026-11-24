# 稠密自然最优出口只读诊断

诊断源码提交 `d200934`，从冻结求解核心 `ad7b4df` 建立独立分支 `codex/dense-optimal-profile`。源码只修改自然 Optimal 返回前的只读钩子、独立模块及 round_probe 的显式采集；没有改求解算序、定价、缓存、状态、容差或恢复策略。源码及 diagnostic-source.patch 均不应合入生产，主线只移植本目录证据。

钩子位于现有终检降级完成之后、objective/y/duals 已按原路径计算之后。只观察真正将返回 Optimal 的 Dense LP，包括根、启发式与探测调用；全部计数在线程内，不随引擎快照回退。失败阶段快照和终检降级调用不计入本统计。仅 proof_profile 编译且探针显式 reset(1) 时启用；非单线程禁用。

每次出口重建未夹界 raw 结构向量，并借用真实 reported x、内部完整 y。按当前内部原行（含割）检查两份向量的界、行及保存的增广等式，按原静态列计算 reduced costs 与一般有界拉格朗日表达式 L(y)。固定变量的侧最优性不要求 r=0，基驻点一致性另列。Max 统一到内部 Min，并反转目标常数；记录的 internal_box_lower_bound 是内部 Min 口径，不能直接叫 Max 模型下界。

L(y) 只作浮点诊断，不是精确证书。某个负 reduced cost 指向无限上界时，L=-∞ 是箱域推导结果，单独统计，不假装有限，也不把它自动叫非有限算术错误。源列、y 或有限乘加出现 NaN/Inf 则显式记异常。人工在 Phase2 的有效界为 [0,0]，L 中贡献为0；只有全部人工实际精确零时才无需原列映射复核增广状态，否则明确 unavailable，绝不反推人工来源。

所有自然 Optimal 都汇总，只有前四个触发诊断阈值的出口最多保存完整状态；后续异常仍累计。每份序列化数据估算上限8MiB，总已保存序列化数据32MiB。前置checked估算包含实际 reported_duals/constraint_row/cut_row 长度，超额保留触发类别等元数据。它不是 allocator 或 RSS 的硬上限。完整数据全部为原始 f64 位模式，坐标明确为 normalized_current_lp；不存在 presolve 后内部模型的独立精确导出，不能仅根据 row_scale 反除、或变量数量相同，就与 outer MIP 坐标等同。

anomaly_categories_mask 的位序与 metrics 字段顺序一致：raw_bounds、reported_bounds、raw_rows、reported_rows、augmented_equations、basic_bounds、static_basis_stationarity、side_optimality、dual_gap_absolute、dual_bound_above_primal、reported_objective_consistency、raw_reported_delta、invalid_nonbasic_side、phase1_pending、nonfinite_dual_data、invalid_dual_bound_arithmetic。strict 与10倍阈值只用于分类，不改变生产终检；10倍分类不表示生产存在同样的所有行检查。

## 验证和唯一公开算例

4项最小 capture 测试全部通过，覆盖 Min/Max/offset 的实际出口输出与状态逐位不变、Clone 不回退统计、终检降级不计、人工污染仅用于采样覆盖及四份上限、实际长度/溢出预算估计。人工污染测试不是自然误最优的证据。

同源 OFF/ON 小背包采用 time-limit=0（无墙钟终止）、node-limit=100；23个非时间字段完全一致，两边 fresh export 哈希相同，原模型向量/目标检查无错误。ON统计11次真实Optimal调用。默认关闭与proof_profile两种二进制均由已提交干净源码构建；hash、特性和命令在 run-manifest.json。

只运行一次 p0201 根节点诊断：沿冻结默认命令，仅替换二进制并把 node-limit 改为1，time-limit仍30秒。实际1节点、6591 LP迭代、65割、LimitReached，根界7227.000000683597，没有MIC incumbent（objective/values均null）。fresh export 与冻结模型哈希相同。独立检查通过的是两份外部解的原模型复核，以及MIC根界没有排除外部已知最优值；不是对不存在的MIC解向量做了检查，也不能验证中间LP剪枝界。未运行整树、额外公开算例或追加重采样。

## 观察结果及限制

188次自然Optimal中，98次触发至少一项诊断阈值。raw 行/界没有超过严格门；reported 行仅序号104超过 feasibility，最大1.529708582026501e-9，未超过10倍门。增广等式98次超过10倍门，最大5.4912391078687506e-8；静态基驻点95次超过 optimality，最大1.7070078683900647e-9，均未超过10倍门。没有非有限原始数据、非法side或人工欠账；188次均可重建增广等式。169次浮点L有限，19次因箱无限端点为-∞。

有限浮点L中有2次高于报告内部目标，最大2.598986611701548e-6，发生在序号183–184；它们尚无精确复核。这些诊断计数不能合并称为98个错误最优，也不能直接把浮点L称为错误下界证书。

严格遵循“前四异常”预算，完整快照仅序号88–91，触发的都是增广等式与静态基驻点一致性；它们的L均为-∞。序号104及183–184被前四采样遮蔽，没有完整快照，不能用88–91代替它们作证明。后续94个异常仅保留汇总。原始完整探针输出保存在 [p0201-root.stdout.gz](p0201-root.stdout.gz)，轻量分类在 p0201-root-summary.json。[深度离线审计](../dense-optimal-offline/README.md)另行处理，不能据本目录预先断言自然误Optimal或错误剪枝。
