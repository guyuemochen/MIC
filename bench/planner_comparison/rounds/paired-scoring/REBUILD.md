# 复现旧版 MIC 对照

在仓库根目录运行以下命令。`before.patch` 将本轮修改过的 Rust 源文件还原到上一轮与 OptaPlanner 对照时的版本。脚本在临时副本应用补丁，逐文件核对上一轮 `results/2026-11-24/environment.json` 中所有 Rust 源文件和基准入口的 SHA-256，再用 release 配置编译。还原后的入口没有引用本轮新增的模块。

```sh
round=bench/planner_comparison/rounds/paired-scoring
workdir=$(mktemp -d)
python3 "$round/rebuild_before.py" --out "$workdir/before"
cargo build --manifest-path planner/Cargo.toml --release --example compare_worker
```

在本轮 Apple M4 上实测，19 个 Rust 文件的指纹全部吻合，临时目录重新编译的可执行文件与正式计时使用的旧版逐字节一致，SHA-256 为 `c9a580c715835bbf2bedd0b46cfdd0fe34550a2b8b39732d54ceb2a27934e5ef`。换用其他工具链或平台时，使用脚本核验源码后打印的二进制指纹运行三方对照，原始记录写入新目录。

```sh
before_sha=$(python3 -c 'import hashlib, pathlib, sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$workdir/before")
python3 "$round/ab.py" run --with-opta --before "$workdir/before" --before-sha256 "$before_sha" \
    --after planner/target/release/examples/compare_worker --out "$workdir/results"
python3 "$round/ab.py" verify --out "$workdir/results"
```

默认排课每个引擎预热两次，每次 2 秒，随后五个种子各给 10 秒。排班、配送、生产排程用两版 MIC 各跑一次 500 毫秒作为回归。快速试跑可加 `--seeds 7 --budget-ms 350 --warmup-ms 45 --regression-ms 80`。Java 比较使用与上一轮相同的 Maven 工程和 JVM 参数。每轮完成后，用同一条 `verify` 命令可以核对算例与模型指纹、运行顺序、首个可行解及完整终解分数。
