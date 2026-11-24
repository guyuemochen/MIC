//! 四个独立业务模型，供 Rust 示例和发布版验收共用。

use planner::business::{Expr, ListSpec, Model, Rule};

fn variable(m: &mut Model, domain: &[i64]) -> usize {
    let index = m.domains.len();
    m.domains.push(domain.to_vec());
    m.expr(Expr::Var(index))
}

fn constraint(m: &mut Model, condition: usize, name: String, hard: bool, weight: i64) {
    let zero = m.expr(Expr::Const(0));
    let one = m.expr(Expr::Const(1));
    let violation = m.expr(Expr::If(condition, zero, one));
    m.rules.push(Rule {
        name,
        hard,
        weight,
        expr: violation,
    });
}

fn roster() -> Model {
    let mut m = Model::default();
    let shifts = [
        variable(&mut m, &[0, 1]),
        variable(&mut m, &[0, 1]),
        variable(&mut m, &[1]),
    ];
    for employee in 0..2 {
        let id = m.expr(Expr::Const(employee));
        let matches = shifts.iter().map(|&s| m.expr(Expr::Eq(s, id))).collect();
        let count = m.expr(Expr::Sum(matches));
        let cap = m.expr(Expr::Const(2));
        let bounded = m.expr(Expr::Le(count, cap));
        constraint(&mut m, bounded, format!("capacity{employee}"), true, 1);
    }
    for (i, shift) in shifts.into_iter().enumerate() {
        let preferred = m.expr(Expr::Const(1));
        let matched = m.expr(Expr::Eq(shift, preferred));
        constraint(&mut m, matched, format!("pref{i}"), false, 2);
    }
    m
}

fn timetable() -> Model {
    let mut m = Model::default();
    let lessons: Vec<_> = (0..3)
        .map(|_| (variable(&mut m, &[0, 1]), variable(&mut m, &[0, 1])))
        .collect();
    let separate = m.expr(Expr::Ne(lessons[0].0, lessons[1].0));
    constraint(&mut m, separate, "teacher".into(), true, 1);
    for i in 0..3 {
        for j in i + 1..3 {
            let same_time = m.expr(Expr::Eq(lessons[i].0, lessons[j].0));
            let same_room = m.expr(Expr::Eq(lessons[i].1, lessons[j].1));
            let violation = m.expr(Expr::Mul(same_time, same_room));
            m.rules.push(Rule {
                name: format!("room{i}{j}"),
                hard: true,
                weight: 1,
                expr: violation,
            });
        }
    }
    m.objective = Some(m.expr(Expr::Sum(lessons.iter().map(|&(time, _)| time).collect())));
    m
}

fn delivery() -> Model {
    let mut m = Model::default();
    m.lists.push(ListSpec {
        owners: 2,
        eligible: vec![vec![0, 1]; 3],
        durations: vec![0; 3],
        precedences: vec![],
    });
    let positions: [i64; 4] = [1, 3, 8, 0]; // 仓库在最后一行和最后一列。
    let matrix: Vec<i64> = positions
        .iter()
        .flat_map(|x| positions.iter().map(move |y| (x - y).abs()))
        .collect();
    let mut distances = Vec::new();
    for owner in 0..2 {
        let load = m.expr(Expr::Load {
            group: 0,
            owner,
            demands: vec![1; 3],
        });
        let cap = m.expr(Expr::Const(2));
        let bounded = m.expr(Expr::Le(load, cap));
        constraint(&mut m, bounded, format!("capacity{owner}"), true, 1);
        distances.push(m.expr(Expr::Distance {
            group: 0,
            owner,
            matrix: matrix.clone(),
        }));
    }
    m.objective = Some(m.expr(Expr::Sum(distances)));
    m
}

fn factory() -> Model {
    let mut m = Model::default();
    m.lists.push(ListSpec {
        owners: 2,
        eligible: vec![vec![0], vec![1], vec![1], vec![0]],
        durations: vec![2, 1, 2, 1],
        precedences: vec![(0, 1), (2, 3)],
    });
    let ends = (0..4)
        .map(|item| m.expr(Expr::Finish { group: 0, item }))
        .collect();
    m.objective = Some(m.expr(Expr::Max(ends)));
    m
}

pub fn models() -> Vec<(&'static str, Model)> {
    vec![
        ("roster", roster()),
        ("timetable", timetable()),
        ("delivery", delivery()),
        ("factory", factory()),
    ]
}
