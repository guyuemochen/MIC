#!/usr/bin/env python3
"""用正式编码器和搜索器重放证明；在临时 crate 中增加只写轨迹的测试入口。"""
import argparse
import gzip
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PROBE = """
#[test]
fn export_trace() {
    use std::{env, io::{BufWriter, Write}, path::PathBuf, fs::File};
    let out = PathBuf::from(env::var_os("MIC_BOOLEAN_AUDIT_DIR").unwrap());
    let model = crate::instance::load_instance(&out.join("model.milp")).unwrap();
    let net = encoding::Encoding::detect(&model, Default::default(), None).unwrap();
    let mut cnf = BufWriter::new(File::create(out.join("network.cnf")).unwrap());
    for (i, j) in net.columns.iter().enumerate() {
        writeln!(cnf, "c var {} {}", i + 1, j).unwrap();
    }
    writeln!(cnf, "p cnf {} {}", net.columns.len(), net.clauses.len()).unwrap();
    for clause in &net.clauses {
        for lit in clause { write!(cnf, "{lit} ").unwrap(); }
        writeln!(cnf, "0").unwrap();
    }
    cnf.flush().unwrap();
    let mut trace = BufWriter::new(File::create(out.join("full.rup")).unwrap());
    writeln!(trace, "c assumption 0").unwrap();
    let mut nodes = 0;
    let status = search::solve(net.columns.len(), &net.clauses, None, 2_000_000, None, &mut nodes,
        &mut |clause| {
            for lit in clause { write!(trace, "{lit} ").unwrap(); }
            writeln!(trace, "0").unwrap();
        });
    trace.flush().unwrap();
    assert_eq!(status, search::Outcome::Unsat);
    println!("完整公式无解。主列 {}，子句 {}，搜索节点 {nodes}。", net.columns.len(), net.clauses.len());
}
"""


def capture(output):
    with tempfile.TemporaryDirectory(prefix='mic-boolean-capture-') as tmp:
        work = Path(tmp)
        crate = work / 'milp'
        crate.mkdir()
        for name in ('Cargo.toml', 'Cargo.lock'):
            shutil.copyfile(ROOT / 'milp' / name, crate / name)
        for name in ('src', 'examples', 'tests'):
            shutil.copytree(ROOT / 'milp' / name, crate / name)
        with (crate / 'src/boolean_network.rs').open('a') as module:
            module.write(PROBE)
        frozen = HERE.parent / 'all-gurobi10-remaining-60s/inputs/bnatt500.milp.gz'
        (work / 'model.milp').write_bytes(gzip.decompress(frozen.read_bytes()))
        subprocess.run(['cargo', 'test', '--locked', '--release', '--lib', 'boolean_network::export_trace',
                        '--', '--exact', '--nocapture'], check=True, cwd=crate,
                       env={**os.environ, 'MIC_BOOLEAN_AUDIT_DIR': str(work)})
        output.mkdir(parents=True, exist_ok=True)
        for name in ('network.cnf', 'full.rup'):
            (output / (name + '.gz')).write_bytes(gzip.compress((work / name).read_bytes(), mtime=0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='存放两个压缩证明文件的目录')
    capture(parser.parse_args().output)
