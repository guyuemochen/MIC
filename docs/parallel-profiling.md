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

## 8. P2 机制轮（2026-10-13 同日续窗）：work-stack A/B 与裁决

> **提交信息勘误（独立复核指认）**：提交 c4fdc89 标题沿用任务名「严格按取证
> 结论实施拆锁与 best-bound 协调，t4/t1 冲击 1.5 倍以上」，与实际内容相反——
> 拆锁与 best-bound 协调**均未实施**（被 §4/§5 证据否决，见 §8.1）；该提交
> 实际内容为探针 `--work-stack` 旗标接线 + 本节文档。标题不可据此考古。

### 8.1 选型（严格按 §4/§5 证据）

- 候选① 拆锁：**P1 数字否决**（六项协调合计 0.24%，收益上限 <1%），
  未实施。
- 候选② best-bound 协调：incumbent 已全局共享 + est 弹出剪枝已存在
  （integer 目标门控，parallel_bb 主循环）；新增协调读只加池锁负载，
  且任何偏离 DFS 序的取件都会重演 global 栈对亲和的破坏（见 8.3 的
  knapsack-500 证据）——无证据支持，未实施。
- 实施 = §5 M1：`--work-stack local|global` 接线进 profile_stats 探针
  （旗标 → `BbOptions.global_work_stack`，默认 local 不变），引擎零改动。

### 8.2 mknap-5x300 同窗紧交错（5 rep 轮转序；单二进制，旗标区分）

| 配置 | 每 rep 墙钟（s） | 中位 | 每 rep 节点 | 中位 | 加速比（t1/tN，>1 为并行收益） |
|---|---|---|---|---|---|
| t1 | 49.5/13.1/29.1/34.0/18.2 | **29.06** | 87,800 ×5（恒同） | 87,800 | 1× |
| t2 local | 85.7/39.5/34.7/27.7/26.0 | 34.66 | 174,251/173,497/170,575/158,289/173,909 | 173,909 | 0.84× |
| t2 global | 40.7/34.6/20.1/20.7/24.5 | **24.46** | 161,153/160,705/146,626/141,009/114,576 | 146,626 | 1.19× |
| t4 local | 18.6/37.0/31.2/20.2/27.4 | 27.43 | 187,111/255,731/192,959/197,085/191,866 | 192,959 | 1.06× |
| t4 global | 24.8/18.7/33.2/33.2/13.3 | **24.81** | 201,751/165,309/163,026/171,959/147,057 | 165,309 | **1.17×** |

- global vs local：t2 墙钟 **−29.4%**、节点 5/5 rep 一致更少（−7.4%~−34.1%,
  中位 −15.6%）；t4 墙钟 −9.6%（中位口径，rep 间噪声 ±40% 内）、节点
  4/5 更少（中位 −14.3%，rep0 +7.8% 例外）。
- **验收线未达：t4/t1 = 1.17×（global）/ 1.06×（local）< 1.5×。**

### 8.3 类别画像（t4，global vs local，各 3 rep）

| 算例 | local 墙钟/节点 | global 墙钟/节点 | 裁决 |
|---|---|---|---|
| mknap-3x100 | 0.074-0.094s / 4,994-7,240 | **0.052-0.069s / 2,444-3,079** | global 显著优（−26%~−45%） |
| knapsack-2000 | 0.438-0.491s / 109-261 | 0.395-0.492s / 112-333 | 持平略优 |
| knapsack-500 | **0.028-0.041s / 665-863** | 0.042-0.109s / 1,031-3,856 | **global 回退 ×1.2-×3.0** |
| knapsack-200 | 0.003-0.004s / 1（根闭合） | 同 | 中性 |
| setcover-100x150 | 0.292-0.325s / 13 | 0.333-0.368s / 13-15 | global 轻微回退（+8%~+26%） |
| milp-mixed-250i40c | 0.012-0.020s | 0.013-0.018s | 持平（ms 级） |

全部 47 次求解（mknap 25 + 类别画像 18 + 探索 4）均 optimal 且目标值
恒同已知最优（4762.0 / 2331.0 / 8975.0 / 35328.0 / 3409.0 / 276.0 /
1275.791752）。

### 8.4 裁决

- **默认保持 local**：global 是按问题类别分化的杠杆（mknap 家族受益、
  knapsack-500 回退至 3×、setcover 轻微回退）——翻默认会以 gate 集
  knapsack-500 类算例的确定性回退换 mknap 的不确定收益，证据不支持。
- **保留增量 = 探针 `--work-stack local|global` 旗标接线**（默认 local，
  全部既有调用行为零变化：t1 三算例行为字段与改动前二进制逐字节一致
  实测；parallel_bb 单测 29/29、acceptance 13/13 通过）。引擎
  （parallel_bb.rs）零改动——P1 证据（协调 0.24%）+ 本轮证据（唯一
  现成杠杆类别分化）共同指向：t4/t1 的 1.5× 在当前「稠密表格 LP +
  共享内存子系统」架构下没有协调层解，剩余路径为 §5 M2a（每节点成本
  ×1.55 的频率/带宽判别）与更强单节点 LP（数据布局），或按类选栈
  （需先有可靠的类检测器——R3 类自适应基础设施可复用，另立轮次）。

## 9. P3 收尾轮（2026-10-13 同日续窗）：global 转正评估 + 类检测收紧

### 9.1 ① `--work-stack global` 转正：**不转正，保持显式开关**

转正前提（「P2 使 t4/t1 达标且 global 参与达成」）不成立：P2 同窗紧交错
5 rep 中位 t4/t1 = **1.17× < 1.5×**（§8.2），且 global 是类别分化杠杆
（§8.3：knapsack-500 回退 ×1.2-3.0、setcover-100x150 轻微回退、mknap 族
受益）——翻默认会以 gate 集算例的确定性回退换噪声内收益。CLI 默认
local 不变（main.rs `--work-stack` 帮助文本与 USAGE 已如实描述）；
milp/README.md 已有 `--work-stack` 条目（参数表）且默认 local 未变、内容仍
准确；根 README 无 work-stack 条目，无需同步。多线程默认口径零变化。

### 9.2 ② detect_pseudocost_favorable 行数下限 1→2：**收紧（A/B 定夺）**

同窗交错 6 rep（探针 `--branch-rule mostfrac` 显式覆盖 ≡ 收紧后单行类
的解析结果；节点/迭代逐 rep 完全确定，墙钟给中位）：

| 算例（类） | def=Pseudocost+probes8 | mf=MostFractional | 裁决 |
|---|---|---|---|
| knapsack-500（单行） | 536 节点/10,485 lpi/~44ms | **360/4,889/~23ms（−48%）** | 单行纯付探测开销 |
| knapsack-2000（单行） | 39/16,994/~273ms | 38/16,787/~269ms | 持平（无收益） |
| knapsack-5000（单行守门） | 52/53,042/~2.34s | **37/38,242/~1.74s（−26%）** | 守门通过（无反向风险） |
| mknap-3x100（3 行） | **2,120/23,978/~88ms** | 7,793/57,892/~268ms | **多行收益 −67% 由 Pseudocost 贡献，保留** |

- 判据（ask 原文）：mknap 族收益主要由多行类贡献 ⇒ 收紧。实测：单行
  族零收益纯开销、多行族收益真实且巨大——收紧成立。
- 实施：`branch_bound.rs::detect_pseudocost_favorable` 行数下限
  `is_empty()`（≥1）→ `len() < 2`（≥2）；三处文档注释同步（枚举/字段/
  函数）。**mknap-5x300 默认行为逐位不变**（t1 87,800 节点/1,747,982
  lpi/4762.0 optimal 实测复核）。
- 测试重钉 3 处（理由入注释，objective 均为已知最优不变）：
  1. `detect_pseudocost_favorable_classification`：命中载体 1 行→2 行，
     新增单行排除向量（原单行正例变排除项）；
  2. `resolve_branch_rule_defaults_paths`：命中载体 p 补第 2 行
     （单行不再命中，三路径语义「命中⇒pseudo+probes≥8 / 未命中⇒
     mostfrac / 显式覆盖优先」不变）；
  3. `acceptance::acceptance_cli_heuristic_on_identical`：hard-knapsack-12
     （单行）nodes 196→159（回落 MostFractional 即 R3 前旧值；objective
     523 = 2^12 穷举真值恒同；显式 `--branch-rule pseudocost` 实测复现
     196，已记入注释）。
- 收紧后默认实测 = A/B 预测值逐位一致：knapsack-500 360/4,889、
  knapsack-2000 38/16,787、knapsack-5000 37/38,242（墙钟 2.34→1.66s）、
  mknap-3x100 2,120/23,978 不变；knapsack-200 根闭合（nodes=1）不变。
- 自检（本轮实跑）：`cargo test --release --lib branch_bound`（39/39）、
  `--lib parallel_bb`（29/29）、`--test acceptance`（13/13，重钉后）；
  t1 确定性双跑（knapsack-5000）行为字段逐字节相同。

## 10. W1 实现轮（2026-10-13 同日续窗）：按类选栈 auto 默认 + 小算例调度门

> 实现 §5 M1（类别分化杠杆的类检测接驾）与 §5 M3（小算例调度门）。
> 引擎（parallel_bb.rs）**零代码改动**——选栈在派发层解析、门在派发层
> 前置；本节全部数字为同窗交错多 rep 中位口径（§1.3 噪声地板 ±18-20%
> 依旧，跨窗绝对墙钟不可比，同窗 A/B 有效）。

### 10.1 实现

1. **按类选栈（`--work-stack auto|local|global`，多线程默认 auto）**：
   `BbOptions` 新增 `work_stack_auto: bool`（默认 true）；解析在
   `branch_bound.rs::parallel_entry`（纯并行概念，threads=1 顺序路径不读
   该字段 ⇒ t1 逐字节不变）：auto ⇒ `detect_pseudocost_favorable`（与 R3
   分支规则**同一检测器**：全 Binary + 正系数 Le 行 + 行数 2..=64）命中
   ⇒ `global_work_stack = true`（mknap 族），未命中（单行背包/等式行
   marketsplit/Ge 行 setcover/连续列 milp-mixed）⇒ local；显式
   local/global 覆盖优先。CLI/探针/USAGE/README 同步（main.rs、
   profile_stats.rs 均接受三态 token，默认 auto）。
2. **小算例调度门（M3）**：threads ≥ 2 ∧ LP 单位 = vars×(rows+1) ≤ 65,536
   ∧（node_limit=0 或 > 512）∧ **处于 auto 决策链**（显式 local/global
   点名并行引擎，不进门——P2 式栈形态 A/B 保持可测量）时，先以 512
   节点预算跑顺序路径探测：预算内出终态 ⇒ 门命中直返（**t4 输出与 t1
   逐字节相同**，同一确定性代码路径）；预算耗尽（LimitReached ∧ nodes
   ≥ 512 ∧ 非 gap/降级）⇒ 树不小，探测作废、并行从零（既有行为 + 有界
   前缀：≤512 个小 LP 节点 + 根阶段重跑）。
   - **512 取值**：两簇空档——并行纯亏簇 t1 节点 ≤ 360（milp-mixed
     53/85、knapsack-500 360、knapsack-2000 38、knapsack-5000 37、
     setcover-100x150 13、hard-knapsack-12 159），并行有收益/大树簇
     ≥ 2,120（mknap-3x100；mknap-5x300 87,800）。
   - **65,536 取值**：门内最大 = lp-wide-200x300（60,200）与
     setcover-100x150（15,150）；门外最小 = assign-40x40（129,600）与
     setcover-300x600（180,600）——以单节点 LP 成本为界钉住探测前缀
     墙钟上限。
3. **测试调整（断言均不变，逐条记因）**：
   - `threads4_matches_default_on_acceptance_instances` /
     `threads4_hard_knapsack_12_repeated_100`（branch_bound）：显式
     node_limit=512 旁路门——这些小算例 W1 起默认被门接走顺序路径，
     而两测试意图是受控入口的**并行**正确性/压测；512 ≥ 完备需求。
   - `acceptance::acceptance_cli_progress_log_format` 的 hard-knapsack-12
     变体：补 `--work-stack local`——「极小周期下并行协调线程必出周期行」
     的覆盖需要并行路径（门命中下亚毫秒顺序搜索可能只有首行）。
   - 新增 4 测：`resolve_work_stack_mode_paths`（解析纯函数全路径）、
     `small_tree_gate_applies_boundaries`（阈/预算边界向量）、
     `small_tree_gate_hit_returns_t1_identical_result`（门命中 t4=t1
     行为字段全等）、`small_tree_gate_miss_dispatches_parallel`
    （mknap-3x100 门未命中 → 并行完备最优 2331.0×3 rep）。
   - `resolve_branch_rule_defaults` 保持原名原义（工作栈解析移至
     parallel_entry，函数只管分支规则；其测试回原断言）。

### 10.2 正确性锚（全部实测，本窗）

- **t1 逐字节不变**：HEAD vs NEW 二进制 t1 行为字段（status/objective/
  nodes/lp_iterations/cuts/presolve_reductions/bound/gap_percent/
  root_bound）**15 算例全等**（knapsack-200/500/2000/5000、
  milp-mixed-120i30c/250i40c、mknap-3x100、setcover-60x100/100x150、
  marketsplit-2x20/4x30/4x40/5x50、assign-16x16、lp-wide-200x300）+
  **mknap-5x300**（87,800 节点 / 1,747,982 lpi / 4762.0 optimal —— R8
  锚逐位一致）；t1 双跑确定性复核（knapsack-5000）逐字节相同。
- **门命中 t4 = t1 逐字节**：12 快算例（上列中 t1 ≤ 512 节点者 +
  lp-wide-200x300 纯 LP）t4 输出与 t1 行为字段全等，含 knapsack-500
  （360）、milp-mixed-250i40c（53）、setcover-100x150（13）、
  marketsplit-2x20/4x30（1/1）、knapsack-2000（38）、knapsack-5000（37）。
- **全部目标值恒为已知最优**：本窗全部求解 objective —— mknap-5x300
  4762.0（20 次：t1×4 + t4 A/B×15 + 复核×1，全部无时限口径）、mknap-3x100 2331.0（13）、mknap-10x200 1791.0（10）、
  knapsack-200/500/2000/5000 = 3409/8975/35328/88526（.opt DP 真值）、
  milp-mixed 1275.791752/980.853268、setcover 276.0/243.0、marketsplit
  10/15/19.0/25.0、assign 165/176、lp-wide 112.668615/112.600268——
  无一偏离。marketsplit-5x50 t4（HEAD 与 NEW 同）2M 节点上限
  limit_reached 携带 incumbent 25.0（钉值不降）。

### 10.3 t4 A/B（同窗交错 5 rep 中位 [min,max]；行为见 §10.2）

| 算例（臂） | NEW auto | 对照臂 | 裁决 |
|---|---|---|---|
| mknap-5x300 t4 无时限 | nodes **122,432** [100,839,175,863]；wall 16.72s [10.7,33.9] | NEW local：225,905 [194,909,253,041] / 37.19s；HEAD local 默认：245,641 / 29.57s | **节点中位 −45.8%**（auto<local 逐 rep 5/5，区间不相交）；墙钟中位 −55%/−43%（区间重叠=噪声，中位与逐 rep 符号一致） |
| mknap-10x200 t4（首次进 global 栈画像） | nodes **92,137** [72,143,101,912]；wall 5.84s [5.2,7.2] | local：111,840 [81,729,114,984] / 5.44s [5.0,9.7] | 节点中位 **−17.6%**；墙钟区间重叠 ⇒ **持平**（无回退） |
| mknap-3x100 t4 | nodes **3,509** [2,437,6,057]；wall 0.12s | local：4,738 [2,927,4,926] / 0.12s | 节点中位 **−26%**；墙钟持平（与 §8.3 global 画像一致） |
| knapsack-500 t4 | **360 节点**（=t1，门命中）；wall 0.075s [0.075,0.081] | HEAD 并行 local：572 [559,595] / 0.110s [0.095,0.123] | 区间不相交 ⇒ **改善 −32%**；节点确定性回归 t1。**不回退 ✓** |
| milp-mixed-250i40c t4 | 53 节点（=t1）；wall 0.046s [0.043,0.064] | HEAD：51-85 节点 / 0.042s [0.041,0.047] | 区间重叠 ⇒ **持平**（验收线达成）；节点数终局确定化 |
| setcover-100x150 t4 | 13 节点（=t1）；wall 0.319s [0.284,0.322] | HEAD：13 节点 / 0.388s [0.373,0.498] | 区间不相交 ⇒ **改善 −18%**（global 回退类 ⇒ 保持 local + 门） |
| assign-40x40 t4（门外，local 不变） | 1 节点；wall 0.893s [0.442,0.903] | HEAD：1 节点 / 0.707s [0.596,1.068] | 区间重叠 ⇒ **持平** ✓ |
| marketsplit-4x40 t4 | 2M 节点 limit_reached obj 20.0；wall 21.6-29.0s | HEAD：同终态同 obj / 25.9-27.6s | **门不缓解**（结构不可能：t1 = 854,952 节点 ≫ 512 探测预算 ⇒ 门未命中走并行；诊断 1 rep `--work-stack global` 同样 2M 上限 obj 20.0——延迟发现最优是树内机理，两种栈都未解，另立轮次） |
| marketsplit-5x50 t4 | 2M limit_reached obj 25.0 | HEAD 同 | 持平（incumbent 钉值不降） |
| lp-wide-1000x2000 t4 | 纯 LP 根终局，HEAD vs NEW 行为字段逐字节相同 | — | 持平 ✓ |

- **显式栈旁路门**（可测量性验证）：knapsack-500 `--work-stack global`
  t4 = 1,092 节点、`local` = 787（均在 §8.3 画像区间 ⇒ 并行确实执行），
  默认 auto = 360（门命中）——显式旗标按契约不进门。
- **门未命中前缀实测**：mknap-5x300 auto t4 无时限 wall 中位 16.7s vs
  local 37.2s——512 节点探测前缀（≈0.2s 量级）被 global 栈收益完全
  吸收，无净损耗。

### 10.4 bench 基线翻转预期（供编排方重钉；30s 时限 t4 对照 HEAD 各 1 rep + 补测）

- knapsack-200 / setcover-100x150 / lp-sparse / setcover-400x800：输出
  无变化（前两者门命中=t1 且 HEAD 恰同值；后两者门外/纯 LP）。
- knapsack-2000 t4：nodes 40（并行）→ 38（=t1，门命中），objective 不变。
- **mknap-5x300 t4 30s**：终态在 optimal 4762 与 auto_degraded
  4746-4755 之间随机翻转——HEAD 4 rep = 3×optimal/1×degraded、
  NEW 4 rep = 2×optimal/2×degraded（机器噪声主导 30s 预算下的降级
  时机，与工作区既录 t1 30s 钉值 auto_degraded 4754 同现象，W3 漂移
  轮的正题）；NEW 节点中位仍低于 HEAD（127-199k vs 197-242k）。
  无时限口径恒 4762.0 optimal（20/20 次，见 §10.2 计数）。
- t4 全快算例 nodes 基线将确定化（= t1 值），并行形态的 t4 基线比较
  需显式 `--work-stack local|global`。

### 10.5 遗留与观察

1. marketsplit-4x40 的 t4 延迟发现最优（t1 19.0@854,952 节点 vs t4
   2M 上限未证优、incumbent 20.0）：调度门结构上不可能缓解（树太大），
   global 栈也不缓解（2M 上限 obj 20.0 诊断 rep）——属树内机理（根界
   弱 + 等式行类），候选路径为根割强化或 presolve（W4/W5 排期）。
2. 门命中时 `--log-secs` 用户在门未命中路径会看到重复首行（probe 与
   并行各发一次）——纯外观，文档已记。
3. mknap-10x200 在 global 栈下节点 −17.6% 但墙钟持平——节点收益被
   每节点成本膨胀吃掉（P1 §4.2 二义未分辨，M2a 仍缺位）。


## 11. W2 实现轮（2026-10-14）：纯 LP 路径内部时限检查（节点边界盲区消除）

> 现状（本轮前）：时限只在节点边界检查，纯 LP 单节点跑满全程——
> lp-sparse-4000x800 gate 名义 30s 实测 39.8-44.0s optimal（W1 窗）/
> 41.3s（工作区既录）。实现：Simplex 迭代循环内周期性条件读钟，超限诚实
> 截断。**红线达成：时限关闭时求解路径零墙钟读点（t1 逐字节不变，17 算例
> 对账），时限开启时 30s 名义值实测 30.1-30.3s 内诚实终态。**

### 11.1 实现

1. **`LpStatus::TimeLimit`（独立状态码，刻意与 IterationLimit 分离）**：
   R8 割循环回滚只对迭代上限触发——墙钟超限若复用同一状态码，「回滚
   重建 + 重放整段根 LP」会把超限点再拖过 deadline（冷审风险 4 的场景）。
   TimeLimit 在割循环出口走 `!= Optimal` 分支如实传播（整树诚实
   LimitReached 终态，不做回滚）；迭代上限触发的回滚若在重放中撞上墙钟
   ⇒ 重建 Simplex **透传原 deadline**（cuts.rs:1760-1762），TimeLimit
   诚实传播、不冒充 Optimal。
2. **条件读钟（红线）**：`Simplex::time_deadline: Option<Instant>`（默认
   None）+ `set_time_deadline`/`time_deadline`；run_phase 循环顶部仅当
   `Some(deadline)` 时递增计数器、每 `TIME_CHECK_INTERVAL = 64` 次主元读
   一次钟（`start_call` 预置计数器 ⇒ 首迭代即检查——重建后已超限的调用
   零主元即返回）。None 态每迭代只有 Option 判别式检查，**零墙钟读点**。
   读钟审计：simplex.rs 全文 `Instant::now()` 恰 1 处（run_phase 检查点，
   `if let Some(deadline)` 门内）——`grep -n "Instant::now" milp/src/
   simplex.rs` 实测仅 simplex.rs:1996（测试内 2709 不计）。
3. **间隔取值（实测定档）**：256 ⇒ lp-sparse 30s 三 rep 30.41/31.20/
   30.55s（rep1 超出 30-31s 验收窗——节流态每主元 ~4ms × 256 ≈ 1s 粒度）；
   收紧到 64 ⇒ 30.17/30.20/30.34s 全入窗。64 次主元读一次钟（~20-50ns）
   对 µs-ms 级主元成本的摊薄份额 <0.001%（「足够稀疏」约束裕量充足），
   且计数器只在注入态递增。
4. **deadline 接线（`branch_bound::bb_time_deadline`，pub(crate) 共用）**：
   `time_limit_ms > 0` ⇒ 求解层起点 + 预算的绝对时刻；锚定沿用既有约定
   （各引擎层自身 start——sequential_bb / parallel 根 / worker 各自锚定；
   W1 门未命中时并行从零重锚，探测前缀不占并行预算、前缀本身 ≤ 512 小
   LP 节点有界）。注入点：sequential_bb 初始 + 预割重建 + R8 回滚重建（3
   处）、parallel 根同 3 处、worker 1 处、cuts.rs 回滚重建透传 1 处。
5. **全匹配点语义**（编译器驱动的 7 处穷尽匹配 + 相等比较审计）：
   - 纯 LP 捷径（顺序/并行根）：TimeLimit ⇒ LimitReached（不携带 best——
     纯 LP 无 incumbent 概念，实测 knapsack-5000 1s 预算 obj=None）；
   - 根终局映射（root_terminal / root_terminal_status）：LimitReached；
   - MILP 根 IterationLimit 终态的 incumbent 提取（R8 ③）扩展到
     TimeLimit（中断点同为原始可行——round_and_check 复核后随
     LimitReached 携带）；
   - 树节点（顺序 classify_node → 新 `NodeAction::TruncateTime`；并行
     worker match）：与 PruneLimit 同处置（剪枝 + limit_hit 粘滞），随后
     节点边界检查立即停树（并行侧 wall_exceeded 置 STOP_TIME）；
   - dive 内部修复失败臂（heuristics.rs）：退化为取整复核（同
     IterationLimit）；dive 回取（dive_back_action）：Discard；
   - probe_one：`== Optimal` 才记账——TimeLimit 自然跳过（探测只读界并
     恢复，无正确性影响）。
6. **文件清单对账**：solver.rs **零改动**（grep 证实无 LpStatus 引用，
   SolveStatus 映射全在 branch_bound）；main.rs **零改动**（输出契约已
   覆盖 limit_reached + best 行）；README 三处口径同步（参数表 / 终态
   诚实判定节 / lp-sparse 算例行）。

### 11.2 验收数字（本轮实跑）

- **lp-sparse-4000x800 30s（目标算例）**：t1 ×3 = 30.17/30.20/30.34s、
  t4 ×1 = 30.11s，全部 `limit_reached` 诚实终态（对照本轮前 39.8-44.0s
  optimal 越窗）——**验收线达成**。
- **无时限路径逐字节不变**：W2 二进制 vs pre-W1 基线二进制 t1 行为字段
  （status/objective/nodes/lp_iterations/cuts/presolve_reductions/bound/
  gap_percent/root_bound）**17 算例全等**——含 lp-sparse（1 节点 / 7,774
  lpi / 5893.705757）、mknap-5x300（87,800 / 1,747,982 / 4762.0）、
  marketsplit-4x40（854,952 节点）、marketsplit-5x50（2M 上限）等全部
  大件。此对账横跨 W1+W2 两轮（pre-W1 基线传递覆盖两轮红线）。
- **R8 交互（冷审调整 3 的用例）**：mknap-5x300 `--time-limit-secs 3` ⇒
  auto_degraded 4727.0@10,962 节点 **cuts=70（割集完整保留）**；
  `--time-limit-secs 10` ⇒ auto_degraded 4754.0@27,176 cuts=70；
  knapsack-5000 `1s` ⇒ limit_reached、wall=1.00s（贴 deadline 截断、
  obj=None 不冒充）；setcover-100x150 `1s` ⇒ optimal 276.0（0.53s 内
  完成，不受扰）。全部终态诚实，无错误状态。
- **gate 集翻转表（t1 + 30s，pre-W2 基线 vs W2 各 1 rep）**：
  knapsack-200/setcover-100x150/knapsack-2000 输出不变；**lp-sparse
  39.8s optimal → 30.1s limit_reached（编排方需重钉的既录基线翻转）**；
  setcover-400x800 两轮同为 limit_reached@nodes=1（该例终态由 24k 迭代
  上限先于/后于 deadline 决定：快 rep 18-25s 迭代上限收场、节流 rep 此后
  由内部墙钟封在 ~30s——本轮前节流态曾越窗至 65s）；mknap-5x300 本窗
  12.8-13.1s optimal@87,800（机器快态未触截断，与 pre-W2 逐位同）。
- **测试**：lib 248/248（新增 4：simplex 条件读钟/截断/清除复原、
  classify_node 两状态码分离向量、setcover-100x150 1ms 预算诚实截断、
  纯 LP 截断不冒充 objective；扩展 1：dive_back_action 补 TimeLimit ⇒
  Discard 断言）、bin milp 21/21、acceptance 13/13、warm_cold_parity
  1/1。全量套件由编排脚本统一跑。

### 11.3 已知边界与交接

1. 超窗残余 = 预处理/收尾开销 + ≤64 主元检查粒度（lp-sparse 实测
   +0.1~0.3s）；deadline 锚在各引擎层自身 start（既有约定），wall 口径
   差异为 solve_controlled 包裹的 presolve/postsolve。
2. 纯 LP 截断不携带 best（无 incumbent 概念）；MILP 根截断经 round_and_check
   复核可携带诚实 incumbent（R8 ③ 口径扩展）。
3. **W3 交接（冷审调整 2/3 的契约点）**：budget 模式（时限换算为节点/
   迭代预算）下，求解路径零墙钟读点的实现口径 = **不注入 deadline**
   （`bb_time_deadline` 是唯一注入点，budget 化时返回 None 即可）+ 既有
   节点边界 `time_exceeded` 的抑制/换算；grep 审计口径见 11.1.2。
4. t4 并行侧 TimeLimit 的终态原因码经 limit_hit/STOP_TIME 复合映射为
   LimitReached（BbResult 不区分 Node/Time 触发原因——CLI 语义同为
   limit_reached，与既有口径一致）。
