# 复测局部搜索

上一轮发布版的排课首解代价为 109。还原脚本在临时副本应用本轮的逆补丁，并用上一轮归档中的 21 个 SHA-256 逐一核对源码，再编译同一个发布版基准入口。在仓库根目录执行以下命令。

本轮在 Apple M4 上从定稿源码还原的 109 版本，其二进制 SHA-256 与冻结程序一致，为 `734f54a8959b417996bf62ccf5ad228aa7db42a2d3956cc4e0c2bf87ac948828`。继续还原上一轮开始前的版本，19 个 Rust 文件和发布版程序也都与首轮归档吻合。

```sh
runner=bench/planner_comparison/rounds/paired-scoring
workdir=$(mktemp -d)
python3 "$runner/rebuild_before.py" --round local-search --out "$workdir/before"
cargo build --manifest-path planner/Cargo.toml --release --example compare_worker
before_sha=$(python3 -c 'import hashlib, pathlib, sys; print(hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest())' "$workdir/before")
python3 "$runner/ab.py" run --with-opta --before "$workdir/before" --before-sha256 "$before_sha" \
    --after planner/target/release/examples/compare_worker --out "$workdir/results"
python3 "$runner/ab.py" verify --out "$workdir/results"
```

三方使用同一算例和种子 7、19、43、73、101，每次给 10 秒，按种子轮换出场顺序。首解和终解由独立的业务规则核对。要复查已归档的结果，直接运行 `python3 "$runner/ab.py" verify --out bench/planner_comparison/rounds/local-search/results`。

需要复现更早的版本时，沿用上一轮的命令 `python3 "$runner/rebuild_before.py" --out "$workdir/earlier"`，脚本会依次还原两轮改动。[上轮复现说明](../paired-scoring/REBUILD.md)

要单独重放旧版每个种子的前 10,000 步，跟踪脚本会在临时目录中还原源码并核验指纹。独立的邻域枚举直接读取上一轮保存的可行解。以下两条命令分别输出归档中的 `trace.json` 和 `landscape.json`。

```sh
round=bench/planner_comparison/rounds/local-search
python3 "$round/trace_before.py"
python3 "$round/diagnose.py"
```
