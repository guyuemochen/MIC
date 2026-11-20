# 根启发式工作步限额验证归档

本候选完成了定向回归、默认配置、新特性、组合特性与 release 验收；目前只有正确性验证，未证明性能收益，也没有新增求解证明。`root_heuristic_work_limit` 默认关闭，启用时自动依赖原有 `root_heuristic_work_budget`。本次归档只复制已有结果，没有再次构建、测试或求解。

源码提交为 `0a41ddece73ea9520dc21b7b8d8e68c4555d3aa5`，基线为 `07f4d2341e9ab86d1b78007467c8f2dbda81a44b`，工作分支为 `codex/heuristic-work-limit`。源码六文件差异保存在 [source.diff](source.diff)，SHA256 为 `134c84da215e935111bfcedd9c6747df74bebe1ce03a3591c4f6d7161e587595`，与 `/root/risk_refresh_validation` 二审记录一致；二审未发现阻断问题，审查者没有改文件或运行测试。[manifest.json](manifest.json) 记录完整命令、工作目录、退出码、每份日志校验和，以及六个源码文件的 Git blob 和 SHA256。

为保留校验和，原始日志末尾空行和 `source.diff` 的空上下文行前缀没有清理；对整份归档执行 `git diff --check` 会报告这些格式提示。README、manifest 和失败摘录单独通过空白检查，源码本身没有新增修改。

## 计费与恢复约束

新口径的一个“工作步”是一次实际执行的换基或变量上下界翻转。余额只借给当前 `resolve` 调用，并贯穿稠密、稀疏实现的两阶段求解、人工变量清理、对偶修复和冷重建；它不属于引擎状态，所以快照和重建不能重新领取额度。原调用者的迭代提示值、绝对期限保持不变。原有 LP 迭代统计，包括根 LP 的迭代统计，均未改口径，不能把工作步写成 LP 迭代数。

批量界翻转（BFRT，即先翻转若干变量的界，再执行一次换基）在修改前检查完整的 `k + 1` 步额度，额度不足则整批不开始。检查不预扣整批费用；每个真实动作开始前才扣一份，未执行的动作不扣，已执行后发生数值失败的动作不退。因此，翻转后、主元执行前失败，只收已执行翻转的费用。不能沿修复失败的回退出口绕过余额限制。

已获准开始的一次原子主元包含原有的周期重构或数值刷新，即使最后一步把余额扣到零，也会完成这项维护。零余额阻止的是下一轮、额外数值重试或回退重建，不能截断已经开始的动作。维护失败后仍保留已扣费用和损坏因子状态，下一次普通求解负责恢复。该规则不引入延迟维护状态，也不保证硬墙钟上限。

预算耗尽走独立限额出口；非最优返回的目标值为 NaN，不把启发式子问题的限额或数值失败当成原问题不可行、无界或最优的证明。第一阶段中断保留待完成标记。预算作用域结束后，普通求解和快照后的普通求解可以继续恢复；根界和根割的恢复、根割主循环、主树搜索均在这项启发式预算作用域之外。

只开启旧 `root_heuristic_work_budget` 时，仍在 `resolve` 边界按 LP 迭代差扣费；新特性才启用调用内工作步计费。原有 `max(1000, 4 × 初始根 LP 迭代数)` 额度和最多 1000 次调用的参数不变，根割前后及嵌套启发式共享额度。旧诊断仍使用 `remaining_iterations`；新特性使用 `remaining_work_steps`。没有按模型调整系数。

## 已有验证结果

下面列出的原始日志均对应最终源码 `0a41dde`，退出码均为 0。完整执行命令保存在 manifest；每条执行时都将标准输出和错误合并重定向到原 `/tmp/mic-work-limit-*.log`，归档副本未改内容。表中“忽略”是 Rust 测试框架原样报告的数量。

| 配置与范围 | 通过数与忽略数 | 原始日志副本 |
| --- | --- | --- |
| 新特性 `--lib work_limit` 专项 | 14 通过 | [targeted.log](targeted.log) |
| 默认关闭，全量 milp | lib 401 通过、4 忽略；bin 29 通过；acceptance 14 通过、2 忽略；warm/cold 1 通过 | [default.log](default.log) |
| 新特性，全量 milp | lib 414 通过、4 忽略；bin 29 通过；acceptance 14 通过、2 忽略；warm/cold 1 通过 | [feature.log](feature.log) |
| 新特性，release acceptance，含 ignored | 16 通过 | [feature-release-acceptance.log](feature-release-acceptance.log) |
| 组合特性，全量 milp | lib 439 通过、4 忽略；bin 29 通过；acceptance 14 通过、2 忽略；warm/cold 1 通过 | [combined.log](combined.log) |
| 组合特性，release acceptance，含 ignored | 16 通过 | [combined-release-acceptance.log](combined-release-acceptance.log) |
| 仅旧预算特性 `--lib work_budget` | 5 通过 | [legacy-budget.log](legacy-budget.log) |
| 全量 planner | lib 103、bin 13、acceptance 11 通过；另 2 个既有集成测试忽略 | [planner.log](planner.log) |

全量测试的文档测试均为 0 项。组合特性参数为 `root_heuristic_work_limit,binary_coefficient_strengthening,risk_pivot_refresh,lu_failure_diagnostics,proof_profile,root_snapshot_restore`；旧预算由依赖自动开启，没有额外开启 `sparse_profile` 或统计特性。

14 项专项覆盖微预算、跨阶段和清理共享余额、两层对偶修复的独立限额出口、批量翻转原子性、失败后已执行动作计费、最后一步周期重构失败、普通求解和快照恢复、调用者提示与期限不变、嵌套启发式共享预算、根恢复在预算外、无限额路径等价、纯界翻转，以及零余额不能提前制造无界或不可行证明。具体测试名称保留在专项日志中。

## 首次失败与证据边界

初次尝试的原始完整日志被后续对同一路径的重定向覆盖，没有保存该未提交源码的精确内容哈希。以下文件仅转录本任务保留下来的工具输出，明确标注为摘录，不能视为原始日志，也不能归因于最终已通过的 `0a41dde` 源码：

- [initial-compile-excerpt.txt](initial-compile-excerpt.txt)：首次编译报 E0560，并行根启发式镜像仍引用已更名的 `iterations` 字段。修正镜像字段与新旧诊断键后继续验证。
- [initial-targeted-excerpt.txt](initial-targeted-excerpt.txt)：首次可执行的 13 项专项为 12 过、1 失败。失败是测试错误假定初始求解后 `just_refactored` 必为真；改为测试准备阶段显式成功重构，建立故障注入的前置条件。随后补入额外反例，最终 14 项全部通过。没有为绕过断言改变生产数学语义。

本归档不包含 FFI/Python 整合结果，该部分由主代理后续验证。R8/R9 正式计时使用既有冻结基线二进制，不含本候选；不能用那些结果推断本候选性能。工作步限额也不承诺一定给主树留下固定墙钟时间。
