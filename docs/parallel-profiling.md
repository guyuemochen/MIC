# milp 并行插桩取证（P1，2026-10-13）：worker 时间分解实测

> 本轮为 R9 空集回退后的**证据轮**（docs/large-scale-report.md §R9：插桩做完
> 即删、数字未进汇报 ⇒ 空集回退）。交付物 = 证据本身：全部数字、由数字导出
> 的结论与下一轮机制提案。插桩代码用后完整移除，净树背书见 §7。

## 0. 摘要表（详细口径见 §2/§3）

**mknap-5x300（主对象），env-on 臂 4 rep 中位，Σworkers 口径占比：**

| 环节 | t2 占比 | t4 占比 |
|---|---|---|
| **LP resolve（节点 resolve + dive 回取）** | **98.20%** | **98.15%** |
| heur（round/dive/probe/branch_children） | 1.22% | 1.29% |
| 池锁等待 pool_wait | 0.04% | 0.04% |
| 池临界区持锁 pool_hold | 0.09% | 0.08% |
| condvar 等待 cv_wait | 0.01% | 0.03% |
| 窃取扫描 scan（⊂pool_hold，单列） | 0.03% | 0.02% |
| incumbent RwLock（读+写） | 0.04% | 0.04% |
| 控制器锁 ctrl（snapshot+decide/record） | 0.03% | 0.03% |
| **六项协调合计** | **0.24%** | **0.24%** |
| startup（Simplex 冷建+割重放） | 0.01% | 0.01% |
| residual（未计时段） | 0.34% | 0.32% |

**判定：既非锁受限、也非等待受限——worker 全速烧 LP（计算吞吐受限）；
并行无效率 = 节点膨胀（t2 ×1.88 / t4 ×2.68）× 每节点 LP 成本膨胀（×1.49/×1.55），
两者乘积恰好吞掉 2×/4× 线程收益（t4/t1 墙钟 ≈ 0.93×）。**
R9 遗留方向「锁竞争拆分」被数字否决：全部锁相关环节合计上限 0.24%，
理论上限收益 <1%，不值得任何复杂度。

**milp-mixed-250i40c（对照小算例）：并行 = 纯开销。** root 单线程段占
root+scope 的 ~42%（t2/t4 同）；t4 下 3/4 worker 饥饿（Σcv_wait 57.29%、
某 worker 中位节点数 = 0），lp_resolve 只占 23.50%。

## 1. 方法论

### 1.1 插桩设计（纯 std Instant，用后整体移除）

- 门控：环境变量 `MILP_PAR_PROF=1` 才计时与捕获；探针把分解打到
  **stderr**（`milp::parallel_bb::take_par_prof_text`），stdout 契约 JSON
  一字不动（bench/run_timings.py 只读 stdout 末行）。
- 每桶语义与重叠规则（实现于当轮临时补丁，随测量结束移除）：
  - `wall`：worker 线程函数入口→返回全寿命；
  - `startup`：Simplex::new + 根割清单重放（每 worker 一次性）；
  - `lp_resolve`：节点 resolve + dive 后记账性回取 resolve；
  - `heur`：round_and_check / dive / probe_reliability / branch_children；
  - `pool_wait`：池 Mutex lock() 阻塞段；`pool_hold`：临界区持锁段
    （wait_timeout 释放锁期间不计——held_since 纪律）；
  - `cv_wait`：wait_timeout 调用全程（含醒后重获取，与 pool_hold 不重叠）；
  - `scan`：take_local_or_steal 调用段（池锁内，pool_hold 子集，单列不重复计）；
  - `inc`：incumbent RwLock read / try_update 全程（worker 调用点）；
  - `ctrl`：控制器 snapshot+decide / record（含其 Mutex；snapshot 内嵌的一次
    incumbent 读计入 ctrl 而非 inc）；
  - `residual = wall − Σ(其余桶)`。
- 每 worker 独立累加（worker 私有结构，无共享写）；观测通道无额外锁竞争。

### 1.2 三臂 × 交错 4 rep

- **head**：HEAD 基线二进制（插桩前构建）——纯净参照；
- **off**：插桩二进制、env 未设（付 Instant 对固定开销、不累加）；
- **on**：插桩二进制、env=1（计时+捕获活跃）——分解数字来源。
- 交错：每 rep 内 14 个配置连续跑，偶数 rep 正序、奇数 rep 逆序；
  每配置 4 rep ≥ 3（纪律下限）。全部最优收敛（见 §1.4 正确性锚）。
- 调用形态复刻 bench：`profile_stats <inst>.milp --threads N`（默认选项：
  类自适应分支、root_cuts on、auto 启发式、无时限无 mipgap——
  maybe_gap/maybe_degrade 全程早退，无进度线程）。

### 1.3 环境与噪声地板（诚实声明）

- 机器：Intel Core Ultra 5 225U（**15W 笔记本 U 系列**，12C/14T）。
- 工具链：本轮进行中遭外部 rustup 升级（1.95.0→1.99.0，msvc 工具链 bin
  一度被清空后恢复）。**全部测量二进制统一 1.99.0 构建**；与 R8 时代
  （1.95.0）绝对墙钟不可比——t1 节点数逐位恒同（87,800）而墙钟从历史
  11.5-13.7s 升到本轮 31-45s（编译器版本与机器功率态二义，未分辨）。
  **占比（run 内比值）不受影响**，跨 run 绝对墙钟对比只在本轮三臂内部有效。
- 噪声地板实证：t1 head 臂同一二进制、同一节点数（87,800）4 rep 墙钟
  **31.2 / 31.4 / 35.8 / 44.5s（±18%）**——U 系列持续稠密 LP 负载下的
  功率/热节流。因此本轮墙钟中位只能分辨 >±20% 的效应；插桩观察者效应
  （结构上界 ~0.3-1%，§4.3）被噪声淹没，无法在墙钟上直接读出。
- 结论稳健性纪律：**分解占比对节流稳健**（各桶随频率同缩放）；每 rep 的
  lp_resolve 占比 t2 ∈ [98.1, 98.4]%、t4 ∈ [98.1, 98.3]%，rep 间极差 <0.3pp。

### 1.4 正确性锚

- mknap-5x300：t1/t2/t4 × 3 臂 × 4 rep = 36 次求解全部 optimal，
  **objective 恒 4762.0**；t1 节点数 87,800 与 R8 基线（final_runs.json /
  large-scale-report §3.4）逐位一致 ⇒ 同一搜索、同一最优。
- milp-mixed-250i40c：全部 24 次求解 objective 恒 1275.791752
  （节点数 53/85 两种搜索形态，均最优）。

## 2. mknap-5x300（主对象）全部数字

### 2.1 墙钟三臂（s；4 rep 明细 + 中位 [min,max]）

| 配置 | head | off | on |
|---|---|---|---|
| t1 | 35.832 / 31.447 / 44.548 / 31.233；med **33.639** [31.2,44.5] | —（顺序路径，无插桩） | — |
| t2 | 55.417 / 42.684 / 112.329 / 31.016；med **49.051** [31.0,112.3] | 64.995 / 45.447 / 104.468 / 27.421；med **55.221** [27.4,104.5] | 93.839 / 27.728 / 48.230 / 31.003；med **39.616** [27.7,93.8] |
| t4 | 30.862 / 22.575 / 38.484 / 31.778；med **31.320** [22.6,38.5] | 52.977 / 21.977 / 20.860 / 50.568；med **36.273** [20.9,53.0] | 54.946 / 39.630 / 24.671 / 40.044；med **39.837** [24.7,54.9] |

t4/t1（head 臂中位）= 31.320/33.639 = **0.93×**；t2/t1 = 1.46×
（t2 因节点膨胀慢于 t1，与 R6 结论一致）。三臂间中位差（−19%~+27%）
逐 rep 符号翻转（t4 on-vs-off 四 rep：+3.7%/+80%/+18%/−21%），
与 §1.3 噪声地板同量级 ⇒ 归因机器噪声，不构成插桩效应读数。

### 2.2 节点数与膨胀（head 臂；膨胀 = nodes/t1 节点数）

| 配置 | 每 rep 节点数 | 每节点成本（µs = wall×threads/nodes，per rep） | 中位膨胀 |
|---|---|---|---|
| t1 | 87,800 ×4（逐位恒同） | 408 / 358 / 507 / 356；med **383** | 1.00 |
| t2 | 171,665 / 171,195 / 159,141 / 152,243 | 646 / 499 / 1412 / 407；med **572**（×1.49） | **1.88** |
| t4 | 190,664 / 250,432 / 245,839 / 225,320 | 647 / 361 / 626 / 564；med **595**（×1.55） | **2.68** |

- 吞吐换算（head 臂中位口径）：t2 聚合吞吐 = 2×(383/572) = **1.34× t1**；
  t4 = 4×(383/595) = **2.58× t1**。墙钟比预测 t2 = 1.88/1.34 = 1.40×
  （实测 1.46×）、t4 = 2.68/2.58 = 1.04×（实测 0.93×，噪声内）——
  **膨胀 × 成本两因子乘积即墙钟**，无残余未解释项需要锁/等待假设。
- env-on 臂节点数（t2：178,856/172,326/168,884/168,134；t4：237,805/
  244,558/239,597/271,694）与 off/head 同量级——插桩不改变搜索形态。

### 2.3 t2 worker 分解（env-on，4 rep 中位，ms）

root_ms 中位 67.077；scope_ms 中位 39,548.027。

| worker | nodes | wall | startup | lp_resolve | heur | pool_wait | pool_hold | cv_wait | scan | inc | ctrl | residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 86,299 | 39,544.4 | 3.864 | 38,840.6 | 474.8 | 14.8 | 36.4 | 2.9 | 10.8 | 17.5 | 11.3 | 131.5 |
| 1 | 85,485 | 39,543.8 | 2.927 | 38,822.1 | 489.2 | 15.0 | 36.9 | 0.0 | 10.3 | 17.3 | 11.5 | 138.6 |
| Σ | 171,784 | 79,088.2 | 7.543 | 77,662.8 | 964.0 | 29.8 | 73.3 | 4.0 | 21.1 | 34.7 | 22.8 | 268.2 |

Σ占比：**lp_resolve 98.20%**、heur 1.22%、pool_wait 0.04%、pool_hold
0.09%、cv_wait 0.01%、scan 0.03%、inc 0.04%、ctrl 0.03%、startup 0.01%、
residual 0.34%（六项协调合计 **0.24%**）。

### 2.4 t4 worker 分解（env-on，4 rep 中位，ms）

root_ms 中位 34.422；scope_ms 中位 39,801.485。

| worker | nodes | wall | startup | lp_resolve | heur | pool_wait | pool_hold | cv_wait | scan | inc | ctrl | residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 60,414 | 39,799.6 | 3.365 | 39,099.7 | 470.6 | 16.3 | 34.0 | 14.0 | 8.6 | 16.4 | 10.3 | 126.4 |
| 1 | 62,318 | 39,799.2 | 3.355 | 39,048.4 | 514.2 | 15.7 | 32.1 | 17.3 | 8.1 | 15.8 | 11.6 | 132.5 |
| 2 | 60,114 | 39,800.6 | 3.195 | 39,055.5 | 531.8 | 15.8 | 33.3 | 7.8 | 9.6 | 17.5 | 11.1 | 115.0 |
| 3 | 60,900 | 39,798.5 | 3.026 | 39,041.5 | 533.5 | 17.1 | 34.5 | 0.7 | 9.6 | 16.1 | 10.5 | 132.0 |
| Σ | 243,746 | 159,197.9 | 12.852 | 156,245.1 | 2,050.1 | 64.96 | 133.9 | 41.3 | 35.9 | 65.8 | 43.5 | 504.5 |

Σ占比：**lp_resolve 98.15%**、heur 1.29%、pool_wait 0.04%、pool_hold
0.08%、cv_wait 0.03%、scan 0.02%、inc 0.04%、ctrl 0.03%、startup 0.01%、
residual 0.32%（六项协调合计 **0.24%**）。

### 2.5 负载均衡与同步尾

- 节点分派极均衡：t2 两 worker 中位 86,299/85,485（差 0.9%）；
  t4 四 worker 60,114-62,318（极差 3.5%）。
- 每 run 内各 worker wall 几乎相同（如 t2 rep3：30,941.468 vs 30,942.093ms，
  差 0.6ms；t4 rep1 四 worker 39,587.3/39,586.8/39,588.2/39,586.0ms）——
  无掉队者、join 尾部浪费可忽略。

### 2.6 每 rep 原始分解行（env-on，ms；全量落档）

```text
mknap-5x300 t2 rep0 w0 nodes=91219  wall=93730.6 startup=15.5 lp_resolve=92217.5 heur=923.9 pool_wait=99.0 pool_hold=100.0 cv_wait=0.0   scan=43.4 inc=31.6 ctrl=22.4
mknap-5x300 t2 rep0 w1 nodes=87865  wall=93728.8 startup=2.5  lp_resolve=92284.9 heur=925.3 pool_wait=66.5 pool_hold=72.0  cv_wait=2.4   scan=18.4 inc=38.1 ctrl=20.8
mknap-5x300 t2 rep1 w0 nodes=87943  wall=27699.8 startup=2.1  lp_resolve=27233.3 heur=304.4 pool_wait=10.8 pool_hold=23.5  cv_wait=16.7  scan=6.2  inc=12.4 ctrl=7.8
mknap-5x300 t2 rep1 w1 nodes=84685  wall=27699.6 startup=2.1  lp_resolve=27254.6 heur=299.5 pool_wait=10.1 pool_hold=23.5  cv_wait=0.0   scan=6.6  inc=11.8 ctrl=7.7
mknap-5x300 t2 rep2 w0 nodes=82865  wall=48147.3 startup=4.3  lp_resolve=47301.8 heur=560.2 pool_wait=17.9 pool_hold=44.6  cv_wait=5.6   scan=13.1 inc=21.4 ctrl=13.1
mknap-5x300 t2 rep2 w1 nodes=86285  wall=48145.5 startup=4.0  lp_resolve=47265.4 heur=591.0 pool_wait=18.6 pool_hold=45.7  cv_wait=0.0   scan=12.6 inc=21.6 ctrl=13.7
mknap-5x300 t2 rep3 w0 nodes=84655  wall=30941.5 startup=3.5  lp_resolve=30379.4 heur=389.4 pool_wait=11.7 pool_hold=28.2  cv_wait=0.2   scan=8.5  inc=13.5 ctrl=9.4
mknap-5x300 t2 rep3 w1 nodes=83783  wall=30942.1 startup=3.3  lp_resolve=30378.9 heur=387.5 pool_wait=11.4 pool_hold=28.1  cv_wait=0.0   scan=7.9  inc=12.9 ctrl=9.4
mknap-5x300 t4 rep0 w0 nodes=60732  wall=54872.3 startup=8.6  lp_resolve=53850.1 heur=688.1 pool_wait=24.4 pool_hold=55.3  cv_wait=22.8  scan=21.3 inc=18.5 ctrl=14.7
mknap-5x300 t4 rep0 w1 nodes=59130  wall=54872.9 startup=7.6  lp_resolve=53740.6 heur=848.9 pool_wait=23.0 pool_hold=49.3  cv_wait=1.3   scan=13.2 inc=18.2 ctrl=12.9
mknap-5x300 t4 rep0 w2 nodes=57724  wall=54870.0 startup=4.7  lp_resolve=54011.5 heur=565.4 pool_wait=28.2 pool_hold=39.8  cv_wait=24.5  scan=10.2 inc=19.7 ctrl=13.1
mknap-5x300 t4 rep0 w3 nodes=60624  wall=54868.5 startup=4.3  lp_resolve=53772.5 heur=766.4 pool_wait=37.0 pool_hold=50.3  cv_wait=2.0   scan=14.2 inc=30.8 ctrl=13.5
mknap-5x300 t4 rep1 w0 nodes=60096  wall=39587.3 startup=3.6  lp_resolve=38919.7 heur=433.3 pool_wait=14.2 pool_hold=32.2  cv_wait=22.0  scan=7.8  inc=16.3 ctrl=10.2
mknap-5x300 t4 rep1 w1 nodes=62550  wall=39586.8 startup=3.5  lp_resolve=38855.4 heur=494.8 pool_wait=15.3 pool_hold=32.4  cv_wait=22.9  scan=7.8  inc=15.9 ctrl=10.4
mknap-5x300 t4 rep1 w2 nodes=61083  wall=39588.2 startup=3.3  lp_resolve=38870.8 heur=508.0 pool_wait=15.5 pool_hold=33.8  cv_wait=3.9   scan=9.9  inc=15.3 ctrl=9.9
mknap-5x300 t4 rep1 w3 nodes=61177  wall=39586.0 startup=3.1  lp_resolve=38871.6 heur=497.6 pool_wait=14.6 pool_hold=34.5  cv_wait=0.0   scan=9.5  inc=16.4 ctrl=10.6
mknap-5x300 t4 rep2 w0 nodes=58136  wall=24650.6 startup=3.0  lp_resolve=24198.1 heur=328.9 pool_wait=8.6  pool_hold=18.6  cv_wait=0.9   scan=4.9  inc=9.1  ctrl=6.1
mknap-5x300 t4 rep2 w1 nodes=62085  wall=24650.3 startup=3.2  lp_resolve=24191.1 heur=318.9 pool_wait=9.9  pool_hold=21.3  cv_wait=11.8  scan=6.4  inc=9.9  ctrl=6.7
mknap-5x300 t4 rep2 w2 nodes=59146  wall=24649.2 startup=3.1  lp_resolve=24222.0 heur=286.0 pool_wait=9.9  pool_hold=20.3  cv_wait=11.8  scan=6.3  inc=9.3  ctrl=6.8
mknap-5x300 t4 rep2 w3 nodes=60513  wall=24649.7 startup=2.9  lp_resolve=24168.5 heur=348.4 pool_wait=8.7  pool_hold=21.3  cv_wait=0.0   scan=6.4  inc=9.7  ctrl=6.4
mknap-5x300 t4 rep3 w0 nodes=66769  wall=40012.0 startup=3.1  lp_resolve=39279.7 heur=507.8 pool_wait=18.5 pool_hold=35.8  cv_wait=5.9   scan=9.4  inc=16.5 ctrl=10.4
mknap-5x300 t4 rep3 w1 nodes=68877  wall=40011.7 startup=3.0  lp_resolve=39241.5 heur=533.6 pool_wait=16.1 pool_hold=31.7  cv_wait=23.9  scan=8.3  inc=15.8 ctrl=12.9
mknap-5x300 t4 rep3 w2 nodes=67306  wall=40013.0 startup=2.9  lp_resolve=39240.1 heur=555.7 pool_wait=16.2 pool_hold=32.8  cv_wait=2.6   scan=9.3  inc=19.8 ctrl=12.4
mknap-5x300 t4 rep3 w3 nodes=69144  wall=40011.0 startup=2.9  lp_resolve=39211.5 heur=569.4 pool_wait=19.5 pool_hold=34.5  cv_wait=1.3   scan=9.7  inc=15.8 ctrl=10.3
```

## 3. milp-mixed-250i40c（对照小算例）全部数字

### 3.1 墙钟（s；head 臂 med：t1 0.027 / t2 0.028 / t4 0.026）

并行无收益（t2/t1 ≈ 1.0×、t4/t1 ≈ 0.97×，量级 20-30ms 全在启动/root 段）。
三臂全配置 objective 恒 1275.791752 optimal；节点数两种形态（53/85）。

### 3.2 t2 分解（env-on，4 rep 中位，ms）

root_ms 中位 10.886；scope_ms 中位 14.518 ⇒ root 占 root+scope 的 **42.9%**。

| worker | nodes | wall | startup | lp_resolve | heur | pool_wait | pool_hold | cv_wait | scan | inc | ctrl | residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 40 | 13.970 | 0.822 | 12.532 | 0.092 | 0.005 | 0.025 | 0.253 | 0.004 | 0.004 | 0.004 | 0.231 |
| 1 | 40 | 13.764 | 0.903 | 11.697 | 0.054 | 0.005 | 0.017 | 0.238 | 0.004 | 0.003 | 0.004 | 0.839 |
| Σ | 80 | 27.733 | 1.725 | 23.346 | 0.167 | 0.010 | 0.043 | 0.687 | 0.007 | 0.006 | 0.007 | 1.734 |

Σ占比：lp_resolve 84.18%、startup 6.22%、cv_wait 2.48%、residual 6.25%、
协调合计 2.76%。

### 3.3 t4 分解（env-on，4 rep 中位，ms）

root_ms 中位 9.228；scope_ms 中位 12.784 ⇒ root 占 **41.9%**。

| worker | nodes | wall | startup | lp_resolve | heur | pool_wait | pool_hold | cv_wait | scan | inc | ctrl | residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 16 | 12.002 | 1.124 | 5.505 | 0.036 | 0.003 | 0.006 | 4.362 | 0.001 | 0.002 | 0.002 | 0.963 |
| 1 | **0** | 11.992 | 1.089 | 0.000 | 0.000 | 0.000 | 0.001 | 9.328 | 0.001 | 0.000 | 0.000 | 1.574 |
| 2 | 6 | 11.627 | 1.022 | 2.627 | 0.012 | 0.001 | 0.021 | 8.542 | 0.001 | 0.001 | 0.001 | −0.601 |
| 3 | 4 | 11.488 | 0.909 | 2.554 | 0.000 | 0.001 | 0.003 | 8.204 | 0.002 | 0.001 | 0.000 | −0.185 |
| Σ | 26 | 47.109 | 4.144 | 11.070 | 0.078 | 0.005 | 0.030 | 26.991 | 0.005 | 0.004 | 0.004 | 4.779 |

Σ占比：**cv_wait 57.29%**、lp_resolve 23.50%、startup 8.80%、residual
10.14%。worker1 中位节点数 = 0（≥2 rep 整段饥饿）；residual 个别 rep 为
负（中位滤波下各桶中位独立聚合的算术效应，|值| < 1.6ms，不影响结论）。

## 4. 结论

### 4.1 受限模式判定（由数字导出）

1. **锁不限**：池锁等待 0.04% + 持锁 0.09%/0.08% + incumbent 0.04% +
   控制器 0.03%——单锁池在 4 线程下竞争可忽略（mknap t2/t4 均如此）。
   **R9 遗留方向「锁竞争拆分/分锁」的收益上限 <1%，数字否决。**
2. **等待不限（大算例）**：cv_wait 0.01%/0.03%、窃取扫描 0.02-0.03%；
   负载均衡 ±3.5% 内、worker 同进同出——协调协议本身工作良好。
3. **计算吞吐受限**：lp_resolve 占 98%。并行无效率拆成两个乘性因子
   （head 臂中位口径）：
   - **节点膨胀**：t2 ×1.88、t4 ×2.68（并行子树交错延迟 incumbent/
     剪枝信息——R6 机理判断的直接量化）；
   - **每节点 LP 成本膨胀**：t2 ×1.49（383→572µs）、t4 ×1.55（→595µs）。
   乘积 t2 = 2.80× vs 2 线程 ⇒ 墙钟 1.46×；t4 = 4.15× vs 4 线程 ⇒ 墙钟
   0.93×（噪声内 ≈1.0×）。两因子解释全部观测，无需第三因子。
4. **小算例（对照）**：并行 = 纯开销。root 单线程段 ~42% + 饥饿等待
   （t4 Σcv_wait 57%）+ 每 worker 冷启动；树太小（≤85 节点）填不满
   2/4 线程。gate 集快算例 t4/t1 ≈ 1× 的历史现象由此解释。

### 4.2 每节点成本膨胀的机制二义（诚实声明）

×1.49-1.55 的每节点多线程劣化有两个候选，本机墙钟数据**无法分离**：
- **功率/频率**：15W U 系列全核 turbo 显著低于单核（仅此一项即可达 ×1.4）；
- **共享内存子系统**：每 worker 全尺寸稠密表格流式读写，DRAM 带宽/L3
  容量竞争。
判别实验（列为 M2a，未在本轮执行）：固定高性能电源模式 + 核绑定
（`start /affinity`：同物理核 HT 对 vs 异物理核）各 4 rep——若同核 HT 对
几乎无额外劣化而异核劣化同现 ⇒ 频率主导；若随物理核数单调劣化 ⇒ 带宽主导。

### 4.3 观察者效应（插桩扭曲）评估

- 结构上界：每节点约 10 个计时段 × 2 次 QPC ≈ 0.5-1µs，t2 on 臂
  178k 节点 ≈ 0.1-0.2s ≈ **scope 的 0.3-0.5%**；累加/捕获开销更小。
- 实测三臂中位差（off−head：+12.6%/+15.8%；on−head：−19.2%/+27.2%）
  逐 rep 符号翻转、幅值与 §1.3 噪声地板（同二进制同节点 ±18%）同量级
  ⇒ 墙钟无法在本机分辨 <1% 的插桩效应；以结构上界为准（≤~0.5%）。
- **占比结论不受影响**：观察者效应按节点均匀分布，进入 lp_resolve 与
  residual 的份额不改变各桶比值 >0.3pp 的任何结论。

## 5. 下一轮机制提案（按杠杆排序）

1. **M1（最高杠杆）：压节点膨胀**。膨胀 ×1.88/×2.68 是最大损失项。
   现成未测开关：`BbOptions::global_work_stack`（R6-retry 全局 DFS 栈，
   parallel_bb.rs 已实现、探针未接线）——提案：给 profile_stats 加
   `--work-stack global|local` 旗标，t2/t4 A/B 膨胀 vs 亲和损失
   （预期：全局栈压缩膨胀但削弱热启动亲和，净效应需实测）。
2. **M2a（判别实验先行）：每节点成本 ×1.5 劣化的机制分离**（§4.2）。
   频率主导 ⇒ 软件无解，收益路径转向 M1/M3；带宽主导 ⇒ simplex 表格
   内存布局（块化/降低足迹）进入备选（工作量大于 M1，先判别再投入）。
3. **M3（小算例调度门）**：solve 层按问题规模/类检测预测小树（或首个
   节点响应时间）直接走 t1 顺序路径——对照算例证明并行对 20-30ms 量级
   求解是纯开销（root 42% + 饥饿 57%），门可白拿 gate 集快算例的
   t4≈t1×1 变成 t4=t1（免并行开销）。
4. **反提案（明确不做，数字背书）**：锁拆分/分锁、无锁池、窃取协议
   调优——六项协调合计 0.24% 封顶全部此类收益；共享割池（R6 已判先决
   不成立）与全局伪代价（R6 实证负收益）维持搁置。

## 6. 复现口径

- 二进制：`cargo build --release --example profile_stats`（rustc 1.99.0，
  2026-10-13）；head 臂 = 插桩前同源构建。
- 命令形态：`MILP_PAR_PROF=1 profile_stats.exe bench/instances/<inst>.milp
  --threads N`（stdout 契约 JSON；stderr parprof 行）。
- 原始数据：本轮测量落盘于 bench/parprof_scratch/（.out/.err × 56 run +
  run.log），**插桩移除时一并删除**——全量数字已内嵌本文 §2/§3。

## 7. 插桩移除与净树验证（已执行）

- [x] `milp/src/parallel_bb.rs` 与 `milp/examples/profile_stats.rs`
      恢复 HEAD 原样（git restore）；bench/parprof_scratch/（56 个 run 的
      .out/.err + 驱动脚本 + 基线二进制副本）全部删除，全量数字已内嵌
      本文 §2/§3。
- [x] 净树：`git status` 仅 `?? docs/parallel-profiling.md`，`git diff`
      为空。注：本轮进行中，先于本轮存在的 bench 三文件改动
      （run_timings.py 的 cargo 回退定位修复 + 编排方 gate 运行产物）被
      **外部行为体还原**（无新提交、无 stash，14:20-14:30 窗口）——如实
      记录，非本轮动作。该修复被还原不构成回归：外部 rustup 修复已重建
      `~/.cargo/bin` 垫片（cargo 1.99.0），miniconda Python 下
      `shutil.which('cargo')` 重新可解析，HEAD 版 run_timings.py 的
      `subprocess.run(["cargo", ...])` 不再复现 WinError 2（本轮末实测）。
- [x] 移除后重构建（release，3.48s）+ threads=1 确定性抽检：knapsack-2000
      t1 双跑，行为字段（objective/nodes/lp_iterations/cuts/bound/
      root_bound/status）**逐字节相同**，仅 wall_secs/stats_time_ms 两个
      计时字段不同（探针 JSON 按设计包含墙钟，非行为差异）。
- [x] 插桩存在期间 bench 契约未破坏：插桩版 stdout 末行 JSON 字段集与
      HEAD 完全一致（MILP_PAR_PROF 未设时 stderr 亦空，t1/t2 双形态
      实测）；分解数据只走 stderr。
