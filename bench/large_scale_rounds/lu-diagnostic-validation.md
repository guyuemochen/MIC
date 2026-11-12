# LU 首次失败诊断验证

诊断默认关闭。对照复用 `round_probe`，开启组在同一线程调用 `reset()`，求解后取出
首条失败。关闭组沿用默认特性。单线程、稀疏单纯形、节点上限 500、无求解时限，
对照全部非耗时字段和完整解。输出记录每份输入和完整解的 SHA-256。

```sh
cargo test --manifest-path milp/Cargo.toml --features lu_failure_diagnostics --lib lu_diagnostic_
cargo test --manifest-path milp/Cargo.toml --features lu_failure_diagnostics
cargo test --manifest-path milp/Cargo.toml --release --features lu_failure_diagnostics --test acceptance -- --include-ignored
cargo test --manifest-path planner/Cargo.toml

python3 bench/large_scale_rounds/run_lu_diagnostic_parity.py \
  /tmp/mic-rounds-public/{p0033,air03,p0282,lseu}.milp \
  --output bench/large_scale_rounds/results/lu-diagnostic-parity.json
```

用原生产路径采集 `bnatt500`，不在这次求解运行额外的基分解或历史重放。
用保存的旧首条失败和历史证据核对完整失败对象，其中包括基、矩阵列和换基记录。

```sh
cargo build --manifest-path milp/Cargo.toml --release --features lu_failure_diagnostics --example lu_failure_probe
milp/target/release/examples/lu_failure_probe \
  /tmp/mic-rounds-miplib2017/bnatt500.milp /tmp/mic-bnatt500-lu-final-capture.json
shasum -a 256 milp/target/release/examples/lu_failure_probe \
  /tmp/mic-rounds-miplib2017/bnatt500.milp \
  /tmp/mic-bnatt500-lu-{final-capture,first-failure,history}.json \
  /tmp/mic-bnatt500-independent-lu.json
python3 - <<'PY'
import json
from pathlib import Path

root = Path("/tmp")
current = json.loads((root / "mic-bnatt500-lu-final-capture.json").read_text())
for name in ("mic-bnatt500-lu-first-failure.json", "mic-bnatt500-lu-history.json"):
    previous = json.loads((root / name).read_text())
    actual, expected = current["failure"], previous["failure"]
    assert actual is not None and expected is not None, name
    fields = actual.keys() | expected.keys()
    changed = sorted(
        key for key in fields
        if key not in actual or key not in expected or actual[key] != expected[key]
    )
    assert not changed, (name, changed)
    print(name, "全部失败字段一致")
PY
```
