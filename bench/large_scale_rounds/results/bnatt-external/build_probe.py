#!/usr/bin/env python3
"""复制求解器到临时目录，用同一个识别器观察命中情况，按需导出较优目标的证明。"""
import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]

PROBE = """
pub fn inspect_external(model: &Problem, proof_dir: Option<&std::path::Path>) -> std::io::Result<()> {
    use std::{fs::File, io::{BufWriter, Write}};
    let Some(net) = encoding::Encoding::detect(model, Default::default(), None) else {
        println!("{{\\"eligible\\":false}}");
        return Ok(());
    };
    println!("{{\\"eligible\\":true,\\"primary_variables\\":{},\\"clauses\\":{},\\"preferred\\":{}}}",
        net.columns.len(), net.clauses.len(),
        net.preferred.map_or_else(|| "null".to_string(), |lit| lit.to_string()));
    if let Some(dir) = proof_dir {
        let mut cnf = BufWriter::new(File::create(dir.join("network.cnf"))?);
        for (i, j) in net.columns.iter().enumerate() {
            writeln!(cnf, "c var {} {}", i + 1, j)?;
        }
        writeln!(cnf, "p cnf {} {}", net.columns.len(), net.clauses.len())?;
        for clause in &net.clauses {
            for lit in clause { write!(cnf, "{lit} ")?; }
            writeln!(cnf, "0")?;
        }
        cnf.flush()?;
        let assumption = net.preferred.expect("需有单列目标");
        let mut proof = BufWriter::new(File::create(dir.join("better.rup"))?);
        writeln!(proof, "c assumption {assumption}")?;
        let mut nodes = 0;
        let status = search::solve(net.columns.len(), &net.clauses, Some(assumption),
            2_000_000, None, &mut nodes, &mut |clause| {
                for lit in clause { write!(proof, "{lit} ").unwrap(); }
                writeln!(proof, "0").unwrap();
            });
        proof.flush()?;
        assert_eq!(status, search::Outcome::Unsat, "较优取值未证无解");
    }
    Ok(())
}
"""
MAIN = """
use std::{env, error::Error, path::PathBuf};
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<PathBuf> = env::args_os().skip(1).map(PathBuf::from).collect();
    if !(1..=2).contains(&args.len()) { return Err("用法：probe 模型.milp [证明目录]".into()); }
    let model = milp::instance::load_instance(&args[0])?;
    milp::inspect_external(&model, args.get(1).map(|p| p.as_path()))?;
    Ok(())
}
"""


def build(output):
    with tempfile.TemporaryDirectory(prefix='mic-external-probe-') as tmp:
        crate = Path(tmp)
        for name in ('Cargo.toml', 'Cargo.lock'):
            shutil.copyfile(ROOT / 'milp' / name, crate / name)
        for name in ('src', 'examples', 'tests'):
            shutil.copytree(ROOT / 'milp' / name, crate / name)
        with (crate / 'src/boolean_network.rs').open('a') as f:
            f.write(PROBE)
        with (crate / 'src/lib.rs').open('a') as f:
            f.write('\npub use boolean_network::inspect_external;\n')
        (crate / 'examples/external_gate.rs').write_text(MAIN)
        subprocess.run(['cargo', 'build', '--locked', '--release', '--example', 'external_gate'],
                       cwd=crate, check=True)
        shutil.copy2(crate / 'target/release/examples/external_gate', output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    build(parser.parse_args().output)
