# 第十二轮根割候选复查归档

`root_cut_recheck` 通过冻结的默认推广门槛。筛查和正式配对均按固定 20 例逐对轮换，
单线程、2000000 节点，两组执行、独立复核、汇总、门禁的退出码都是 0/0/0/0。
根节点诊断另行启动，执行和复核退出码为 0/0。

- 10 秒筛查共 40 条实测记录。基础组 7 例证优，候选组 8 例证优，增加的是 `lseu`。
  172 条含参考记录、152 个可行解、66 条可对照的 MIC 界独立核验通过。6 例进入耗时
  几何平均，候选/基础组耗时比为 0.63046979。
- 30 秒正式配对共 120 条实测记录。基础组合计 27 次证优，候选组 30 次。分别有 9 例和
  10 例三遍证优，新增 `p0201`，没有逐次丢失证明。252 条含参考记录、214 个可行解及
  126 条可对照的 MIC 界核验通过。8 个有效模型的候选/基础中位耗时比几何平均为
  0.56345409，最大共同证优单例中位比为 1.00437160，全部 60 对峰值内存合格。
- `p0201` 基础组三遍均在 30 秒耗尽预算，已知可行目标 7615、界 7227.000000639117，
  保留 65 条根割。候选组三遍以 29.539–29.582 秒证得 7615，均为 22 条割、17147 个节点、
  1508731 次 LP 迭代。`lseu` 两组均证得 1120，中位耗时由 17.586 秒降到 2.875 秒，
  LP 迭代由 1560614 降到 568182，保留根割由 84 条降到 44 条。
- 现代 8 例两组共 48 次均为 LimitReached。两个 `drayage` 每次均只有 1 个节点、
  0 条割、无可行解，同模型界各自稳定为约 11714.142。全部 12 次 LP 迭代在 171662 到 173751 之间。

## 根节点诊断

`p0201` 与 `setcover-300x600` 分别执行两组，节点上限 1，时间预算 30 秒。探针中的每个
候选原序号和压缩后位置均按实际小批验账。只有候选组在小批返回 Optimal 且候选未用完时
出现复查，随后添加的行按余池稳定顺序施加。`p0201` 前两轮分别有 25、29 个候选，
初批优化后分别跳过 14、18 项，每轮施加 11 条。第三轮两组均回退失败的 10 条，
前缀分别保留 65、22 条，结束界分别约 7227 和 7215。`setcover-300x600` 两组在首轮均无
通过卫生检查的候选，根界同为 131.46305476904783，根诊断均有可行目标 145。归档保留
逐轮及逐项顺序。四条诊断和 108 条相关参考记录均通过独立原模型检查。

## 来源与复算

`original/` 保留执行命令、构建清单、测试日志、逐次 JSONL 原件及原始报告。`input-recovery.json` 记录
冻结 20 例的原模型与转换文件摘要。旧临时目录清空后，从两份冻结清单的 ZIB 官方 URL
重新下载 16 个公开模型，并用 HiGHS 1.15.1 和仓库已有转换器恢复 `.milp`。压缩文件、
解压 MPS 和转换文件逐级核对 SHA-256 后通过临时目录链接恢复原基线命令的路径。
两份基线 JSONL 及仓库原模型未改写。另行重验 132 条参考结果并保存报告。

两组正式探针与两组诊断探针都来自同一干净的 `5f22868` 工作树，MILP 子树为
`c6846f32888cc3b3b0c0307dbe8f075056569c4e`。两组共享 `sparse_row_products` 和
`binary_coefficient_strengthening`，正式候选仅多开 `root_cut_recheck`。诊断两组另开
`proof_profile`。`original/builds/` 记录有效特性、构建命令、工具链和探针摘要，筛查与正式配对
使用同一对正式探针。系统计时器预先验证能提供正值峰值内存。

`original/validation/` 保留测量前开发验证的日志和退出码。03 默认全测与 04 候选全测均从开发工作树
启动，03 后的最后编辑涉及注释与候选专属代码。其后记录包括 planner 全测、两组各 16 项 release
验收、候选与诊断组合的库测试，以及候选与根快照恢复组合的专项测试，矩阵退出码均为 0。
这些测试日志未独立记录构建时的源码树哈希，正式与诊断探针的树和二进制摘要见构建清单。

各阶段归档版 JSONL 逐行只改 `model` 字段，指向仓库 `../../models/`。`archive-relocation.json`
登记原件和归档副本摘要、退出码、以及独立复核、汇总和门禁重算是否一致。全量结果仍带原始求解命令。
以下命令从仓库根目录重算正式配对的独立复核和门禁，无需求解器二进制。

```sh
rounds=bench/large_scale_rounds
refs="$rounds/results"
archive="$refs/root-cut-recheck"
raw="$archive/paired/root_cut_recheck-paired.jsonl"
replay=$(mktemp -d)
python3.13 "$rounds/check_results.py" "$refs/external-baseline.jsonl" \
  "$refs/miplib2017-baseline.jsonl" "$raw" --public-truth "$refs/miplib2017-sources.json" \
  --output "$replay/check.json"
python3.13 "$rounds/summarize_feature_trial.py" "$raw" \
  --legacy-baseline "$refs/external-baseline.jsonl" --modern-baseline "$refs/miplib2017-baseline.jsonl" \
  --public-truth "$refs/miplib2017-sources.json" --output "$replay/summary.json"
python3.13 "$rounds/gate_feature_trial.py" "$replay/summary.json" --output "$replay/gates.json"
```

单独重放根诊断核验。

```sh
PYTHONPATH="$rounds" python3.13 "$archive/root-cut-diagnostics.py" summarize \
  --results "$archive/diagnostics/root_cut_recheck-diagnostics.jsonl" \
  --baseline "$refs/external-baseline.jsonl" --output "$replay/diagnostics.json"
```
