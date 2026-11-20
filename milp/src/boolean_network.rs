//! 从原约束识别二层布尔网络，消去可回填的列。先判定整个网络是否有解；
//! 单列二元目标最多有两个取值，只有首解取了较差值才需要第二次搜索。

use crate::branch_bound::BbStats;
use crate::error::{MilpError, Result};
use crate::model::{Problem, Sense, SolveStatus, Tolerances};
use std::time::Instant;

mod encoding;
mod search;

pub(crate) fn solve(
    problem: &Problem,
    tol: Tolerances,
    deadline: Option<Instant>,
    node_limit: u64,
    accept: &mut dyn FnMut(&[f64]),
) -> Result<Option<(SolveStatus, Option<f64>, BbStats)>> {
    let Some(net) = encoding::Encoding::detect(problem, tol, deadline) else {
        return Ok(None);
    };
    let n = net.columns.len();
    let mut stats = BbStats::default();
    // 识别器只接受精确整数尺度上的零项或单项二元目标。原变量盒给出中断时可用的界。
    let box_bound = problem.obj_offset
        + problem
            .obj
            .terms
            .first()
            .map_or(0., |&(_, c)| match problem.sense {
                Sense::Minimize => c.min(0.),
                Sense::Maximize => c.max(0.),
            });
    let mut best = None;
    let mut assumption = None;
    let (status, bound) = loop {
        let result = search::solve(
            n,
            &net.clauses,
            assumption,
            node_limit,
            deadline,
            &mut stats.nodes,
            &mut |_| {},
        );
        match result {
            search::Outcome::Limit => break (SolveStatus::LimitReached, Some(box_bound)),
            search::Outcome::Unsat if assumption.is_some() => break (SolveStatus::Optimal, best),
            search::Outcome::Unsat => break (SolveStatus::Infeasible, None),
            search::Outcome::Sat(bits) => {
                if bits.len() != n {
                    return Err(MilpError::Numerical("布尔赋值维数与主列数不一致".into()));
                }
                let preferred = net
                    .preferred
                    .is_none_or(|lit| bits[lit.unsigned_abs() as usize - 1] == (lit > 0));
                if assumption.is_some() && !preferred {
                    return Err(MilpError::Numerical("布尔赋值违反目标假设".into()));
                }
                let mut values = net.restore(&bits);
                if values.len() != problem.vars.len()
                    || values
                        .iter()
                        .zip(&problem.vars)
                        .any(|(&x, v)| !(x == 0. || x == 1.) || x < v.lb || x > v.ub)
                {
                    return Err(MilpError::Numerical("回填的布尔解违反原变量界".into()));
                }
                crate::solver::final_check(problem, &mut values)?;
                let objective = problem.obj.eval(&values) + problem.obj_offset;
                accept(&values);
                best = Some(objective);
                if preferred {
                    break (SolveStatus::Optimal, best);
                }
                assumption = net.preferred;
            }
        }
    };
    Ok(Some((status, bound, stats)))
}
