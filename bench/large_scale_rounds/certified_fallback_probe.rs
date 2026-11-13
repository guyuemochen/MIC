//! 复用生产预处理和认证函数，核对旧行回退的载荷与输出；不创建 LP 或搜索树。
#[allow(dead_code)]
#[path = "../../milp/src/binary_propagation.rs"]
mod binary_propagation;
#[allow(dead_code)]
#[path = "../../milp/src/presolve.rs"]
mod presolve;

use milp::{error, model, simplex};
use std::{error::Error, path::Path};

fn identical_problem(a: &model::Problem, b: &model::Problem) {
    // 结构比较覆盖名称、类型、列号和关系；位比较另覆盖全部浮点字段。
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    let bits = |p: &model::Problem| -> Vec<u64> {
        p.vars.iter().flat_map(|v| [v.lb.to_bits(), v.ub.to_bits()])
            .chain(std::iter::once(p.obj_offset.to_bits()))
            .chain(p.obj.terms.iter().map(|&(_, value)| value.to_bits()))
            .chain(p.constraints.iter().flat_map(|row| std::iter::once(row.rhs.to_bits())
                .chain(row.expr.terms.iter().map(|&(_, value)| value.to_bits()))))
            .collect()
    };
    assert_eq!(bits(a), bits(b));
}

fn identical_map(a: &presolve::PostsolveMap, b: &presolve::PostsolveMap) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    let bits = |map: &presolve::PostsolveMap| -> Vec<u64> {
        map.fixed.iter().map(|&(_, value)| value.to_bits())
            .chain(map.aggregations.iter().filter_map(|(_, _, value)| value.map(f64::to_bits)))
            .collect()
    };
    assert_eq!(bits(a), bits(b));
}

fn main() -> Result<(), Box<dyn Error>> {
    assert!(cfg!(feature = "certified_row_fallback"), "须显式开启候选特性编译探针");
    let tol = model::Tolerances::default();
    for path in std::env::args_os().skip(1) {
        let path = Path::new(&path);
        let original = milp::instance::load_instance(path)?;
        original.validate()?;
        let (scaled, _) = milp::scaling::scale_problem(&original);
        let (plain, plain_map) = presolve::presolve(&scaled, &tol)?;
        let (candidate, candidate_map, preserved) = presolve::presolve_for_solve(&scaled, &tol, true)?;
        identical_problem(&plain, &candidate);
        identical_map(&plain_map, &candidate_map);
        for value in [0.0, 1.0] {
            let reduced_values = vec![value; plain.vars.len()];
            let a = presolve::postsolve(&plain_map, &reduced_values, original.vars.len(), &tol);
            let b = presolve::postsolve(&candidate_map, &reduced_values, original.vars.len(), &tol);
            assert_eq!(a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                b.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        }
        let plain_stats = binary_propagation::BinaryPropagation::new(&plain, tol).stats;
        let candidate_stats = binary_propagation::BinaryPropagation::with_preserved(&candidate, tol, preserved)
            .expect("生产载荷须通过逐位见证核验").stats;
        assert!(candidate_stats.accepted_rows >= plain_stats.accepted_rows);
        assert_eq!(candidate_stats.accepted_rows - plain_stats.accepted_rows,
            plain_stats.skipped_rows - candidate_stats.skipped_rows);
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\ttrue\ttrue", path.display(),
            candidate.vars.len(), candidate.constraints.len(), plain_stats.accepted_rows,
            plain_stats.skipped_rows, candidate_stats.accepted_rows, candidate_stats.skipped_rows);
    }
    Ok(())
}
