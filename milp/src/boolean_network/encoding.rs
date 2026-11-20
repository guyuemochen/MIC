//! 从完整的二层布尔网络中消去固定列与合取列，只保留决定网络状态的输出列。
//! 每一对同支撑约束都必须在整个 0/1 真值表上恰好定义一列，否则交还通用搜索。

use crate::model::{Constraint, Problem, RelOp, Sense, Tolerances, VarType};
use std::collections::BTreeMap;
use std::time::Instant;

enum Definition {
    Fixed(bool),
    // 每项记录使该文字为真的原列值。
    All(Vec<(usize, bool)>),
    Any(Vec<(usize, bool)>),
}

pub(super) struct Encoding {
    pub(super) columns: Vec<usize>,
    pub(super) clauses: Vec<Vec<i32>>,
    pub(super) preferred: Option<i32>,
    definitions: Vec<Definition>,
}

// 只接纳能在同一二进制定点格上精确表示的原 f64，乘 2 的幂保持精确。
// 单项至多 2^82，九项之和及与右端之差均严格落在 i128 内。
const UNIT: f64 = (1u64 << 62) as f64;

fn scaled(v: f64) -> Option<i128> {
    let t = v * UNIT;
    (v.abs() <= (1u32 << 20) as f64 && t.fract() == 0.0).then_some(t as i128)
}

struct Row {
    coeff: Vec<i128>,
    rhs: i128,
    op: RelOp,
    margin: i128,
}

impl Row {
    fn new(r: &Constraint, tol: Tolerances) -> Option<Self> {
        let norm: f64 = r.expr.terms.iter().map(|&(_, a)| a.abs()).sum();
        let scale = 1.0 + r.rhs.abs() + norm;
        let allowance = (10.0 * tol.feasibility + tol.integrality * norm)
            .max(crate::solver::FINAL_VIOLATION_REL * scale);
        // 向外留足九项浮点求和、积分容差扰动和最后阈值换算的舍入余量。
        let margin = ((allowance + 64.0 * f64::EPSILON * (scale + allowance)) * UNIT).ceil();
        Some(Self {
            coeff: r
                .expr
                .terms
                .iter()
                .map(|&(_, a)| scaled(a))
                .collect::<Option<_>>()?,
            rhs: scaled(r.rhs)?,
            op: r.op,
            margin: margin as i128,
        })
    }

    fn accepts(&self, mask: usize) -> Option<bool> {
        let delta = self
            .coeff
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1 << i) != 0)
            .map(|(_, a)| a)
            .sum::<i128>()
            - self.rhs;
        let violation = match self.op {
            RelOp::Le => delta,
            RelOp::Ge => -delta,
            RelOp::Eq => delta.abs(),
        };
        if violation <= 0 {
            Some(true)
        } else if violation > self.margin {
            Some(false)
        } else {
            // 精确数学上违反、求解容差却可能接纳，不能借此删去赋值。
            None
        }
    }
}

// 若输出为 1 的整行可行赋值只有一种，该输出是若干带符号输入的合取；
// 若为 0 的可行赋值只有一种，则是析取。还须逐个输入验证输出存在且唯一。
fn gate(cols: &[usize], pair: &[&Constraint], tol: Tolerances) -> Option<(usize, Definition)> {
    if pair.len() != 2 {
        return None;
    }
    let a = Row::new(pair[0], tol)?;
    let b = Row::new(pair[1], tol)?;
    let count = 1 << cols.len();
    let mut table = Vec::with_capacity(count);
    for m in 0..count {
        table.push(a.accepts(m)? && b.accepts(m)?);
    }
    if table.iter().filter(|&&v| v).count() != count / 2 {
        return None;
    }
    for (i, &output) in cols.iter().enumerate() {
        let bit = 1 << i;
        let ones = (0..count).filter(|&m| table[m] && m & bit != 0).count();
        let all = ones == 1;
        if !all && count / 2 - ones != 1 {
            continue;
        }
        if (0..count).any(|m| table[m] == table[m ^ bit]) {
            continue;
        }
        let rare = (0..count).find(|&m| table[m] && (m & bit != 0) == all)?;
        let args = cols
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != i)
            .map(|(k, &j)| (j, (rare & (1 << k) != 0) == all))
            .collect();
        return Some((
            output,
            if all {
                Definition::All(args)
            } else {
                Definition::Any(args)
            },
        ));
    }
    None
}

impl Encoding {
    pub(super) fn detect(p: &Problem, tol: Tolerances, deadline: Option<Instant>) -> Option<Self> {
        let n = p.vars.len();
        if !(3..=16_384).contains(&n)
            || !(n..=2 * n).contains(&p.constraints.len())
            || p.obj.terms.len() > 1
            || !crate::model::objective_is_integral(p)
            || p.obj_offset.abs() + p.obj.terms.first().map_or(0.0, |&(_, c)| c.abs())
                >= (1u64 << 52) as f64
            || !tol.feasibility.is_finite()
            || !tol.integrality.is_finite()
            || tol.feasibility < 0.0
            || tol.integrality < 0.0
            || p.vars
                .iter()
                .any(|v| v.var_type != VarType::Binary || v.lb != 0.0 || v.ub != 1.0)
        {
            return None;
        }
        let expired = || deadline.is_some_and(|t| Instant::now() >= t);
        let mut pairs: BTreeMap<Vec<usize>, Vec<&Constraint>> = BTreeMap::new();
        let mut defs: Vec<Option<Definition>> = (0..n).map(|_| None).collect();
        for row in &p.constraints {
            if expired() {
                return None;
            }
            let terms = &row.expr.terms;
            if terms.len() == 1 {
                let r = Row::new(row, tol)?;
                let (zero, one) = (r.accepts(0)?, r.accepts(1)?);
                if zero == one || defs[terms[0].0].replace(Definition::Fixed(one)).is_some() {
                    return None;
                }
            } else if (3..=9).contains(&terms.len()) {
                let key: Vec<_> = terms.iter().map(|&(j, _)| j).collect();
                let pair = pairs.entry(key).or_default();
                if pair.len() == 2 {
                    return None;
                }
                pair.push(row);
            } else {
                return None;
            }
        }
        for (cols, pair) in pairs {
            if expired() {
                return None;
            }
            let (j, def) = gate(&cols, &pair, tol)?;
            if defs[j].replace(def).is_some() {
                return None;
            }
        }
        let definitions: Vec<_> = defs.into_iter().collect::<Option<_>>()?;
        let columns: Vec<_> = definitions
            .iter()
            .enumerate()
            .filter_map(|(j, d)| matches!(d, Definition::Any(_)).then_some(j))
            .collect();
        if columns.is_empty() || !definitions.iter().any(|d| matches!(d, Definition::All(_))) {
            return None;
        }
        for d in &definitions {
            if let Definition::All(args) = d {
                if args
                    .iter()
                    .any(|&(j, _)| !matches!(definitions[j], Definition::Any(_)))
                {
                    return None;
                }
            }
        }
        let preferred = match p.obj.terms.first() {
            None => None,
            Some(&(j, c)) => match definitions[j] {
                Definition::Fixed(_) => None,
                Definition::Any(_) => {
                    let lit = columns.binary_search(&j).ok()? as i32 + 1;
                    Some(if (c > 0.0) == (p.sense == Sense::Maximize) {
                        lit
                    } else {
                        -lit
                    })
                }
                Definition::All(_) => return None,
            },
        };
        let mut enc = Self {
            columns,
            clauses: Vec::new(),
            preferred,
            definitions,
        };
        for i in 0..enc.columns.len() {
            if expired() || enc.clauses.len() > 200_000 {
                return None;
            }
            enc.encode_equation(i)?;
        }
        enc.clauses.sort_unstable();
        enc.clauses.dedup();
        (!expired()).then_some(enc)
    }

    fn encode_equation(&mut self, i: usize) -> Option<()> {
        let output = self.columns[i];
        let Definition::Any(args) = &self.definitions[output] else {
            return None;
        };
        let mut scope = vec![output];
        for &(j, _) in args {
            match &self.definitions[j] {
                Definition::Fixed(_) => {}
                Definition::All(inputs) => scope.extend(inputs.iter().map(|&(k, _)| k)),
                Definition::Any(_) => return None,
            }
        }
        scope.sort_unstable();
        scope.dedup();
        // 至多六个状态位的整张真值表可放进一个 u64，较大的并集交还通用搜索。
        let width = scope.len();
        if width > 6 {
            return None;
        }
        let size = 1usize << width;
        let universe = u64::MAX >> (64 - size);
        let mut ones = [0u64; 6];
        for (j, t) in ones[..width].iter_mut().enumerate() {
            for state in 0..size {
                if state & (1 << j) != 0 {
                    *t |= 1u64 << state;
                }
            }
        }
        let mut rhs = 0u64;
        for &(j, sign) in args {
            let term = match &self.definitions[j] {
                Definition::Fixed(v) => {
                    if *v {
                        universe
                    } else {
                        0
                    }
                }
                Definition::All(inputs) => {
                    let mut term = universe;
                    for &(k, value) in inputs {
                        let t = ones[scope.binary_search(&k).ok()?];
                        term &= if value { t } else { !t };
                    }
                    term
                }
                Definition::Any(_) => return None,
            };
            rhs |= if sign { term } else { universe ^ term };
        }
        let allowed = universe & !(rhs ^ ones[scope.binary_search(&output).ok()?]);
        let lits: Vec<_> = scope
            .iter()
            .map(|j| self.columns.binary_search(j).ok().map(|k| k as i32 + 1))
            .collect::<Option<_>>()?;
        // 对每个禁配生成一条子句。只有去掉一个条件后仍不含任何允许取值，
        // 才可缩短它。所有原禁配都被覆盖，任何允许取值都不会被排除。
        for state in 0..size {
            if allowed & (1u64 << state) != 0 {
                continue;
            }
            let mut keep = (1usize << width) - 1;
            for k in 0..width {
                let trial = keep ^ (1 << k);
                let mut possible = allowed;
                for (j, &t) in ones[..width].iter().enumerate() {
                    if trial & (1 << j) != 0 {
                        possible &= if state & (1 << j) != 0 { t } else { !t };
                    }
                }
                if possible == 0 {
                    keep = trial;
                }
            }
            self.clauses.push(
                lits.iter()
                    .enumerate()
                    .filter(|&(j, _)| keep & (1 << j) != 0)
                    .map(|(j, &lit)| if state & (1 << j) == 0 { lit } else { -lit })
                    .collect(),
            );
        }
        Some(())
    }

    pub(super) fn restore(&self, bits: &[bool]) -> Vec<f64> {
        let mut x = vec![0.0; self.definitions.len()];
        for (&j, &v) in self.columns.iter().zip(bits) {
            x[j] = f64::from(v);
        }
        for (j, def) in self.definitions.iter().enumerate() {
            match def {
                Definition::Fixed(v) => x[j] = f64::from(*v),
                Definition::All(args) => {
                    x[j] = f64::from(args.iter().all(|&(k, sign)| (x[k] == 1.0) == sign));
                }
                Definition::Any(_) => {}
            }
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // 用整数的两条不等式定义一个门，只靠完整真值识别输出与输入极性。
    fn add_gate(p: &mut Problem, output: usize, args: &[(usize, bool)], all: bool) {
        let n = args.len() as f64;
        let neg = args.iter().filter(|(_, v)| !v).count() as f64;
        let expr = |factor: f64, out: f64| {
            let mut terms: Vec<_> = args
                .iter()
                .map(|&(j, sign)| (j, if sign { factor } else { -factor }))
                .collect();
            terms.push((output, out));
            terms
        };
        let ((a, b), (c, d)) = if all {
            ((expr(-1.0, n), neg), (expr(1.0, -1.0), n - 1.0 - neg))
        } else {
            ((expr(1.0, -n), -neg), (expr(-1.0, 1.0), neg))
        };
        p.add_constraint(None, a, RelOp::Le, b).unwrap();
        p.add_constraint(None, c, RelOp::Le, d).unwrap();
    }

    fn network(positive: bool) -> Problem {
        let mut p = Problem::new("布尔网络", Sense::Minimize);
        for k in 0..5 {
            p.add_var(&format!("x{k}"), VarType::Binary).unwrap();
        }
        p.set_objective(vec![(0, 2.0)], -7.0).unwrap();
        p.add_constraint(None, vec![(4, 0.5)], RelOp::Eq, 0.0)
            .unwrap();
        add_gate(&mut p, 2, &[(0, true), (1, false)], true);
        add_gate(&mut p, 3, &[(0, false), (1, true)], true);
        add_gate(&mut p, 0, &[(2, true), (3, true)], false);
        add_gate(&mut p, 1, &[(3, positive), (4, true)], false);
        p
    }

    fn satisfies(p: &Problem, x: &[f64]) -> bool {
        p.constraints.iter().all(|r| {
            let a = r.expr.eval(x);
            match r.op {
                RelOp::Le => a <= r.rhs,
                RelOp::Ge => a >= r.rhs,
                RelOp::Eq => a == r.rhs,
            }
        })
    }

    #[test]
    fn reconstruction_and_clauses_match_every_assignment() {
        for positive in [true, false] {
            let mut p = network(positive);
            let enc = Encoding::detect(&p, Tolerances::default(), None).unwrap();
            assert_eq!(enc.columns, [0, 1]);
            assert_eq!(enc.preferred, Some(-1));
            let mut count = 0;
            for mask in 0..1u32 << p.vars.len() {
                let point: Vec<_> = (0..p.vars.len())
                    .map(|j| f64::from(mask & (1 << j) != 0))
                    .collect();
                let bits: Vec<_> = enc.columns.iter().map(|&j| point[j] == 1.0).collect();
                let cnf = enc.clauses.iter().all(|c| {
                    c.iter()
                        .any(|&lit| bits[lit.unsigned_abs() as usize - 1] == (lit > 0))
                });
                let restored = enc.restore(&bits);
                assert_eq!(cnf, satisfies(&p, &restored));
                assert_eq!(satisfies(&p, &point), cnf && point == restored);
                count += usize::from(satisfies(&p, &point));
            }
            assert_eq!(count, if positive { 2 } else { 0 });
            p.sense = Sense::Maximize;
            assert_eq!(
                Encoding::detect(&p, Default::default(), None)
                    .unwrap()
                    .preferred,
                Some(1)
            );
            p.set_objective(vec![(1, -3.0)], 4.0).unwrap();
            assert_eq!(
                Encoding::detect(&p, Default::default(), None)
                    .unwrap()
                    .preferred,
                Some(-2)
            );
            p.set_objective(vec![(4, 2.0)], 1.0).unwrap();
            assert_eq!(
                Encoding::detect(&p, Default::default(), None)
                    .unwrap()
                    .preferred,
                None
            );
        }
    }

    #[test]
    fn controlled_search_keeps_the_first_solution_across_objective_searches() {
        use crate::boolean_network::search::{self, Outcome};
        use crate::branch_bound::SearchControl;
        use crate::solver::{solve_controlled_with_incumbent, SolveOptions};
        use std::cell::RefCell;

        let mut p = network(true);
        let encoding = Encoding::detect(&p, Default::default(), None).unwrap();
        let mut first_nodes = 0;
        let Outcome::Sat(first) = search::solve(
            encoding.columns.len(),
            &encoding.clauses,
            None,
            u64::MAX,
            None,
            &mut first_nodes,
            &mut |_| {},
        ) else {
            panic!("测试网络有两个可行解");
        };
        // 让第一个赋值的目标较差。恰好用完首轮额度时要保留该解和原尺度的界。
        p.sense = if first[0] {
            Sense::Minimize
        } else {
            Sense::Maximize
        };
        let first_cost = -7. + 2. * f64::from(first[0]);
        let best_cost = -7. + 2. * f64::from(!first[0]);
        for limit in [first_nodes, 100] {
            let mut opts = SolveOptions::default();
            opts.bb.node_limit = limit;
            let seen = RefCell::new(Vec::new());
            let r = solve_controlled_with_incumbent(
                &p,
                &opts,
                &SearchControl::default(),
                None,
                &|_, z, x| {
                    assert!(satisfies(&p, x));
                    assert_eq!(z, p.obj_offset + p.obj.eval(x));
                    seen.borrow_mut().push(z);
                },
            )
            .unwrap();
            assert!(r.report.stats.nodes <= limit);
            assert_eq!(r.bound, Some(best_cost));
            if limit == first_nodes {
                assert_eq!(r.report.status, crate::model::SolveStatus::LimitReached);
                assert_eq!(r.report.objective, Some(first_cost));
                assert_eq!(seen.into_inner(), [first_cost]);
            } else {
                assert_eq!(r.report.status, crate::model::SolveStatus::Optimal);
                assert_eq!(r.report.objective, Some(best_cost));
                assert_eq!(seen.into_inner(), [first_cost, best_cost]);
            }
        }
        // 第二个主列在两个可行解里都为零。较优假设不可行只证明原目标最优。
        p.sense = Sense::Maximize;
        p.set_objective(vec![(1, 3.)], -5.).unwrap();
        let r = crate::solver::solve(&p, &SolveOptions::default()).unwrap();
        assert_eq!(r.status, crate::model::SolveStatus::Optimal);
        assert_eq!(r.objective, Some(-5.));
    }

    #[test]
    fn ambiguous_grid_extra_rows_and_incomplete_gates_decline() {
        let valid = network(true);
        let detect = |p| Encoding::detect(p, Tolerances::default(), None).is_some();
        let mut p = valid.clone();
        p.constraints[0].rhs = 2.0_f64.powi(-32);
        assert!(!detect(&p));
        let mut p = valid.clone();
        p.constraints[0].expr.terms[0].1 = 2.0_f64.powi(-80);
        assert!(!detect(&p));
        let mut p = valid.clone();
        p.constraints[1].rhs += 1.0;
        assert!(!detect(&p));
        let mut p = valid.clone();
        p.constraints.push(p.constraints[1].clone());
        assert!(!detect(&p));
        let mut p = valid.clone();
        p.constraints.pop();
        assert!(!detect(&p));
        let mut p = valid.clone();
        p.set_objective(vec![(3, 2.0)], 0.0).unwrap();
        assert!(!detect(&p));
        assert!(Encoding::detect(&valid, Default::default(), Some(Instant::now())).is_none());
        assert!(Encoding::detect(
            &valid,
            Tolerances {
                integrality: 1.0,
                ..Default::default()
            },
            None
        )
        .is_none());
    }
}
