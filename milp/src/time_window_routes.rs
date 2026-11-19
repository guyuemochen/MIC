//! 从原模型提取带时间窗的路径覆盖。行、列和变量名均不承担识别语义。
//!
//! 每个二元列连接两条恰选一行。时间差分行里的负系数列是选中的
//! 行程，正系数列表示跳过同一出发任务。两列共享的恰选一行确定
//! 出发端，跳过列的另一端则把同一任务的入边与出边配为一个节点。

use crate::model::{Constraint, Problem, RelOp, Sense, VarType};
use std::collections::VecDeque;
use std::time::Instant;

/// 上限约束识别器的二次空间；更大的模型留给通用求解器。
const MAX_NODES: usize = 256;
const MISSING: usize = usize::MAX;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Window {
    pub(crate) lower: f64,
    pub(crate) upper: f64,
}

#[derive(Debug)]
pub(crate) struct RouteModel {
    pub(crate) arcs: Vec<Vec<Option<usize>>>,
    pub(crate) timing: Vec<Option<Window>>,
    /// 保存原行右端和激活列系数，实际时延是 rhs - coefficient。
    pub(crate) duration: Vec<Vec<Option<(f64, f64)>>>,
    /// 最小化口径的原列目标系数。
    pub(crate) costs: Vec<f64>,
    clock_columns: Vec<Option<usize>>,
    arc_ends: Vec<Option<(usize, usize)>>,
}

/// 浮点和与舍入误差（TwoSum）。有限和的上、下端都保留精确相等。
pub(crate) fn sum_bounds(a: f64, b: f64) -> Option<(f64, f64)> {
    let sum = a + b;
    let portion_b = sum - a;
    let correction = (a - (sum - portion_b)) + (b - portion_b);
    if !sum.is_finite() || !correction.is_finite() {
        return None;
    }
    let lower = if correction < 0.0 {
        sum.next_down()
    } else {
        sum
    };
    let upper = if correction > 0.0 { sum.next_up() } else { sum };
    (lower.is_finite() && upper.is_finite()).then_some((lower, upper))
}

fn set_once<T: Copy + PartialEq>(slot: &mut Option<T>, value: T) -> Option<()> {
    if slot.is_some_and(|existing| existing != value) {
        return None;
    }
    *slot = Some(value);
    Some(())
}

/// 每个方向的时间行统一为 >=，以单位系数直接读取时间差。
fn ge_sign(row: &Constraint) -> Option<f64> {
    match row.op {
        RelOp::Ge => Some(1.0),
        RelOp::Le => Some(-1.0),
        RelOp::Eq => None,
    }
}

struct DurationRow {
    earlier_clock: usize,
    later_clock: usize,
    selected: usize,
    skipped: usize,
    rhs: f64,
    selected_coefficient: f64,
}

/// 所有形似时间差的行必须获得完整解释。
fn duration_row(problem: &Problem, row: &Constraint) -> Option<DurationRow> {
    if row.expr.terms.len() != 4 {
        return None;
    }
    let sign = ge_sign(row)?;
    let mut earlier = None;
    let mut later = None;
    let mut selected = None;
    let mut skipped = None;
    for &(j, a) in &row.expr.terms {
        let a = sign * a;
        let cell = match problem.vars[j].var_type {
            VarType::Continuous if a == -1.0 => &mut earlier,
            VarType::Continuous if a == 1.0 => &mut later,
            VarType::Binary if a < 0.0 => &mut selected,
            VarType::Binary if a > 0.0 => &mut skipped,
            _ => return None,
        };
        if cell.replace((j, a)).is_some() {
            return None;
        }
    }
    let (selected, selected_coefficient) = selected?;
    Some(DurationRow {
        earlier_clock: earlier?.0,
        later_clock: later?.0,
        selected,
        skipped: skipped?.0,
        rhs: sign * row.rhs,
        selected_coefficient,
    })
}

/// 活动行的两个互斥列必须只共享一条恰选一行。
fn shared_endpoint(selected: [usize; 2], skipped: [usize; 2]) -> Option<usize> {
    let common = selected.map(|g| skipped.contains(&g));
    match common {
        [true, false] => Some(selected[0]),
        [false, true] => Some(selected[1]),
        _ => None,
    }
}

impl RouteModel {
    /// 只接受从每条原始行、变量界和目标都能核实的结构。
    pub(crate) fn detect(problem: &Problem) -> Option<Self> {
        let clock_count = problem
            .vars
            .iter()
            .filter(|v| v.var_type == VarType::Continuous)
            .count();
        if clock_count == 0
            || clock_count >= MAX_NODES
            || problem.vars.len() > MAX_NODES * MAX_NODES + MAX_NODES
            || problem
                .obj
                .terms
                .iter()
                .any(|&(j, _)| problem.vars[j].var_type != VarType::Binary)
        {
            return None;
        }

        let mut groups = Vec::new();
        let mut endpoints = vec![[MISSING; 2]; problem.vars.len()];
        for row in problem.constraints.iter().filter(|r| r.op == RelOp::Eq) {
            if row.rhs == 0.0 || row.expr.terms.is_empty() || groups.len() == 2 * MAX_NODES {
                return None;
            }
            for &(j, a) in &row.expr.terms {
                if a != row.rhs || problem.vars[j].var_type != VarType::Binary {
                    return None;
                }
                let ends = &mut endpoints[j];
                if ends[0] == MISSING {
                    ends[0] = groups.len();
                } else if ends[1] == MISSING {
                    ends[1] = groups.len();
                } else {
                    return None;
                }
            }
            groups.push(row);
        }
        if groups.len() <= 2 * clock_count || groups.len() % 2 != 0 {
            return None;
        }
        for (j, var) in problem.vars.iter().enumerate() {
            match var.var_type {
                VarType::Continuous => {}
                VarType::Binary
                    if var.lb == 0.0
                        && (var.ub == 0.0 || var.ub == 1.0)
                        && endpoints[j][1] != MISSING => {}
                _ => return None,
            }
        }

        let mut windows: Vec<Window> = problem
            .vars
            .iter()
            .map(|v| Window {
                lower: v.lb,
                upper: v.ub,
            })
            .collect();
        for row in problem
            .constraints
            .iter()
            .filter(|r| r.expr.terms.len() == 1)
        {
            let (j, a) = row.expr.terms[0];
            if problem.vars[j].var_type == VarType::Continuous {
                let sign = ge_sign(row)?;
                if sign * a == 1.0 {
                    windows[j].lower = windows[j].lower.max(sign * row.rhs);
                } else if sign * a == -1.0 {
                    windows[j].upper = windows[j].upper.min(-sign * row.rhs);
                } else {
                    return None;
                }
            }
        }
        for (j, var) in problem.vars.iter().enumerate() {
            if var.var_type == VarType::Continuous
                && (!windows[j].upper.is_finite() || windows[j].lower > windows[j].upper)
            {
                return None;
            }
        }

        let mut difference_rows = Vec::new();
        for row in problem.constraints.iter().filter(|r| r.op != RelOp::Eq) {
            let clock_terms = row
                .expr
                .terms
                .iter()
                .filter(|&&(j, _)| problem.vars[j].var_type == VarType::Continuous)
                .count();
            if clock_terms == 2 {
                let edge = duration_row(problem, row)?;
                let duration = sum_bounds(edge.rhs, -edge.selected_coefficient)?.0;
                let min_inactive = sum_bounds(
                    windows[edge.later_clock].lower,
                    -windows[edge.earlier_clock].upper,
                )?
                .0;
                if duration <= 0.0 || min_inactive < edge.rhs {
                    return None;
                }
                difference_rows.push(edge);
            } else if clock_terms > 2 || !redundant_on_windows(problem, row, &windows)? {
                return None;
            }
        }

        // 从已知出发端进行二分染色。未连通或方向冲突时不猜测。
        let first = difference_rows.first()?;
        let departure = shared_endpoint(endpoints[first.selected], endpoints[first.skipped])?;
        let mut directions = vec![None; groups.len()];
        directions[departure] = Some(true);
        let mut todo = VecDeque::from([departure]);
        while let Some(g) = todo.pop_front() {
            let opposite = !directions[g]?;
            for &(col, _) in &groups[g].expr.terms {
                let [a, b] = endpoints[col];
                let h = if a == g { b } else { a };
                if directions[h].is_none() {
                    todo.push_back(h);
                }
                set_once(&mut directions[h], opposite)?;
            }
        }
        if directions.contains(&None) {
            return None;
        }
        let count = groups.len() / 2;
        let mut nodes_by_group = vec![None; groups.len()];
        let mut clock_columns = Vec::with_capacity(count);
        for (g, out) in directions.iter().enumerate() {
            if *out == Some(true) {
                nodes_by_group[g] = Some(clock_columns.len());
                clock_columns.push(None);
            }
        }
        if clock_columns.len() != count {
            return None;
        }
        let orient = |j| -> (usize, usize) {
            let [a, b] = endpoints[j];
            if directions[a] == Some(true) {
                (a, b)
            } else {
                (b, a)
            }
        };
        let mut clock_nodes = vec![None; problem.vars.len()];
        let mut duration = vec![vec![None; count]; count];
        for edge in &difference_rows {
            let (depart, _) = orient(edge.selected);
            let (self_depart, self_arrive) = orient(edge.skipped);
            if depart != self_depart {
                return None;
            }
            let from = nodes_by_group[depart]?;
            set_once(&mut clock_columns[from], edge.earlier_clock)?;
            set_once(&mut clock_nodes[edge.earlier_clock], from)?;
            set_once(&mut nodes_by_group[self_arrive], from)?;
            // 等所有出边锚点读完后，再核对入组和到达时刻。
        }
        for (j, v) in problem.vars.iter().enumerate() {
            if v.var_type == VarType::Continuous && clock_nodes[j].is_none() {
                return None;
            }
        }
        for edge in difference_rows {
            let (depart, arrive) = orient(edge.selected);
            let from = nodes_by_group[depart]?;
            let to = clock_nodes[edge.later_clock]?;
            if nodes_by_group[arrive] != Some(to)
                || from == to
                || duration[from][to]
                    .replace((edge.rhs, edge.selected_coefficient))
                    .is_some()
            {
                return None;
            }
        }

        // 无时间节点的两侧可以按各自行号配对。需要实际站点对称的算法仍须
        // 检查这些节点的所有边及成本，它们的配对本身不增加原模型约束。
        let mut stations = clock_columns
            .iter()
            .enumerate()
            .filter_map(|(node, clock)| clock.is_none().then_some(node));
        let mut arrival_seen = vec![false; count];
        for (g, out) in directions.iter().enumerate() {
            if *out == Some(false) {
                let node = match nodes_by_group[g] {
                    Some(node) => node,
                    None => {
                        let node = stations.next()?;
                        nodes_by_group[g] = Some(node);
                        node
                    }
                };
                if std::mem::replace(&mut arrival_seen[node], true) {
                    return None;
                }
            }
        }
        let mut arcs = vec![vec![None; count]; count];
        let mut arc_ends = vec![None; problem.vars.len()];
        for (j, var) in problem
            .vars
            .iter()
            .enumerate()
            .filter(|(_, v)| v.var_type == VarType::Binary)
        {
            let (depart, arrive) = orient(j);
            let (from, to) = (nodes_by_group[depart]?, nodes_by_group[arrive]?);
            if var.ub == 0.0 {
                continue;
            }
            if (from != to
                && clock_columns[from].is_some()
                && clock_columns[to].is_some()
                && duration[from][to].is_none())
                || arcs[from][to].replace(j).is_some()
            {
                return None;
            }
            arc_ends[j] = Some((from, to));
        }
        let mut costs = vec![0.0; problem.vars.len()];
        let sign = if problem.sense == Sense::Minimize {
            1.0
        } else {
            -1.0
        };
        for &(j, a) in &problem.obj.terms {
            costs[j] = sign * a;
        }
        let timing = clock_columns
            .iter()
            .map(|clock| clock.map(|j| windows[j]))
            .collect();
        Some(Self {
            arcs,
            timing,
            duration,
            costs,
            clock_columns,
            arc_ends,
        })
    }

    pub(crate) fn idle_candidate(&self, problem: &Problem) -> Option<Vec<f64>> {
        let choices = self
            .arcs
            .iter()
            .enumerate()
            .map(|(i, row)| row[i])
            .collect::<Option<Vec<_>>>()?;
        self.candidate_from_arcs(problem, &choices)
    }

    /// 按节约值拼接任务链。仅在线路起终点与具体站点无关时采用，
    /// 拼接时沿后半条链传播到达时刻，每个窗口都必须能容纳到达时间。
    pub(crate) fn greedy_candidate(
        &self,
        problem: &Problem,
        deadline: Option<Instant>,
    ) -> Option<Vec<f64>> {
        let (stations, jobs): (Vec<usize>, Vec<usize>) =
            (0..self.arcs.len()).partition(|&node| self.timing[node].is_none());
        let station = *stations.first()?;
        if stations.iter().any(|&a| {
            stations
                .iter()
                .any(|&b| self.arcs[a][b].is_none_or(|j| self.costs[j] != 0.0))
        }) {
            return None;
        }
        let mut depart_cost = vec![0.0; self.arcs.len()];
        let mut return_cost = vec![0.0; self.arcs.len()];
        for &job in &jobs {
            depart_cost[job] = self.costs[self.arcs[station][job]?];
            return_cost[job] = self.costs[self.arcs[job][station]?];
            if stations.iter().any(|&s| {
                self.arcs[s][job].is_none_or(|j| self.costs[j] != depart_cost[job])
                    || self.arcs[job][s].is_none_or(|j| self.costs[j] != return_cost[job])
            }) {
                return None;
            }
        }

        let mut routes: Vec<Vec<usize>> = jobs.iter().map(|&i| vec![i]).collect();
        let mut arrival: Vec<f64> = jobs
            .iter()
            .map(|&i| self.timing[i].map(|w| w.lower))
            .collect::<Option<_>>()?;
        loop {
            if deadline.is_some_and(|limit| Instant::now() >= limit) {
                return None;
            }
            let mut merge = None;
            for (a, left) in routes.iter().enumerate() {
                let Some(&end) = left.last() else { continue };
                for (b, right) in routes.iter().enumerate() {
                    if a == b {
                        continue;
                    }
                    let Some(&begin) = right.first() else {
                        continue;
                    };
                    let Some(arc) = self.arcs[end][begin] else {
                        continue;
                    };
                    let saving = return_cost[end] + depart_cost[begin] - self.costs[arc];
                    if saving <= 0.0
                        || !saving.is_finite()
                        || merge.is_some_and(|(score, col, _, _, _)| {
                            saving < score || (saving == score && arc >= col)
                        })
                    {
                        continue;
                    }
                    if let Some(last_arrival) = self.extend_time(end, arrival[a], right) {
                        merge = Some((saving, arc, a, b, last_arrival));
                    }
                }
            }
            let Some((_, _, a, b, last_arrival)) = merge else {
                break;
            };
            let suffix = std::mem::take(&mut routes[b]);
            routes[a].extend(suffix);
            arrival[a] = last_arrival;
        }

        // 每条链最多占一站。优先选能省下最多闲置成本的链，
        // 余下任务保持闲置，空站使用本身的入、出配对边。
        let mut ranked = Vec::new();
        for route in routes.iter().filter(|route| !route.is_empty()) {
            let start = route[0];
            let end = *route.last()?;
            let mut saving = -depart_cost[start] - return_cost[end];
            for &job in route {
                saving += self.costs[self.arcs[job][job]?];
            }
            for edge in route.windows(2) {
                saving -= self.costs[self.arcs[edge[0]][edge[1]]?];
            }
            if saving.is_finite() && saving > 0.0 {
                ranked.push((saving, self.arcs[start][start]?, route));
            }
        }
        ranked.sort_by(|a, b| b.0.total_cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
        let mut choices: Vec<_> = self
            .arcs
            .iter()
            .enumerate()
            .map(|(i, row)| row[i])
            .collect::<Option<Vec<_>>>()?;
        for (station, (_, _, route)) in stations.into_iter().zip(ranked) {
            choices[station] = self.arcs[station][route[0]]?;
            for edge in route.windows(2) {
                choices[edge[0]] = self.arcs[edge[0]][edge[1]]?;
            }
            let end = *route.last()?;
            choices[end] = self.arcs[end][station]?;
        }
        self.candidate_from_arcs(problem, &choices)
    }

    fn extend_time(&self, mut from: usize, mut earliest: f64, suffix: &[usize]) -> Option<f64> {
        for &to in suffix {
            let window = self.timing[to]?;
            let (rhs, coefficient) = self.duration[from][to]?;
            let duration = sum_bounds(rhs, -coefficient)?.1;
            earliest = sum_bounds(earliest, duration)?.1.max(window.lower);
            if earliest > window.upper {
                return None;
            }
            from = to;
        }
        Some(earliest)
    }

    /// 从候选原列复原最早时刻，走完整条链并逐点核对时间窗。
    /// 最后直接重算原模型的全部行、界与目标。
    pub(crate) fn candidate_from_arcs(
        &self,
        problem: &Problem,
        choices: &[usize],
    ) -> Option<Vec<f64>> {
        let count = self.arcs.len();
        if choices.len() != count {
            return None;
        }
        let mut values = vec![0.0; problem.vars.len()];
        let mut departures = vec![false; count];
        let mut arrivals = vec![false; count];
        let mut next = vec![None; count];
        let mut has_previous = vec![false; count];
        for &j in choices {
            let (from, to) = (*self.arc_ends.get(j)?)?;
            if std::mem::replace(&mut departures[from], true)
                || std::mem::replace(&mut arrivals[to], true)
            {
                return None;
            }
            values[j] = 1.0;
            if from != to && self.timing[from].is_some() && self.timing[to].is_some() {
                let (rhs, coefficient) = self.duration[from][to]?;
                next[from] = Some((to, sum_bounds(rhs, -coefficient)?.1));
                has_previous[to] = true;
            }
        }
        let mut visited = 0usize;
        for (start, &previous) in has_previous.iter().enumerate() {
            if previous || self.timing[start].is_none() {
                continue;
            }
            let mut at = start;
            let mut earliest = self.timing[start]?.lower;
            loop {
                let window = self.timing[at]?;
                earliest = earliest.max(window.lower);
                if earliest > window.upper {
                    return None;
                }
                values[self.clock_columns[at]?] = earliest;
                visited += 1;
                let Some((to, duration)) = next[at] else {
                    break;
                };
                earliest = sum_bounds(earliest, duration)?.1;
                at = to;
            }
        }
        if visited != self.timing.iter().flatten().count() {
            return None; // 正时延环不可能满足所有差分行。
        }
        verify_candidate(problem, &values).then_some(values)
    }
}

/// 对单时刻行和纯二元不等式，以盒下界核实整行已经被蕴含。
fn redundant_on_windows(problem: &Problem, row: &Constraint, windows: &[Window]) -> Option<bool> {
    let sign = ge_sign(row)?;
    let mut lower = 0.0;
    for &(j, coefficient) in &row.expr.terms {
        let a = sign * coefficient;
        let activity = match problem.vars[j].var_type {
            VarType::Binary if a < 0.0 && problem.vars[j].ub == 1.0 => a,
            VarType::Binary => 0.0,
            VarType::Continuous if a == 1.0 => windows[j].lower,
            VarType::Continuous if a == -1.0 => -windows[j].upper,
            _ => return None,
        };
        lower = sum_bounds(lower, activity)?.0;
    }
    Some(lower >= sign * row.rhs)
}

fn verify_candidate(problem: &Problem, values: &[f64]) -> bool {
    const ROW_CHECK_TOL: f64 = 1.0e-7;
    if !(problem.obj.eval(values) + problem.obj_offset).is_finite()
        || problem
            .vars
            .iter()
            .zip(values)
            .any(|(v, &x)| !x.is_finite() || x < v.lb || x > v.ub)
    {
        return false;
    }
    problem.constraints.iter().all(|row| {
        let lhs = row.expr.eval(values);
        lhs.is_finite()
            && match row.op {
                RelOp::Le => lhs - row.rhs <= ROW_CHECK_TOL,
                RelOp::Ge => row.rhs - lhs <= ROW_CHECK_TOL,
                RelOp::Eq => (lhs - row.rhs).abs() <= ROW_CHECK_TOL,
            }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy() -> Problem {
        let mut p = Problem::new("三个作业", Sense::Minimize);
        let arcs: Vec<Vec<_>> = (0..4)
            .map(|from| {
                (0..4)
                    .map(|to| {
                        p.add_var(&format!("y_{from}_{to}"), VarType::Binary)
                            .unwrap()
                    })
                    .collect()
            })
            .collect();
        let times: Vec<_> = (0..3)
            .map(|i| {
                let j = p.add_var(&format!("t_{i}"), VarType::Continuous).unwrap();
                p.set_var_bounds(j, 0.0, 10.0).unwrap();
                j
            })
            .collect();
        let mut costs = Vec::new();
        for i in 0..4 {
            p.add_constraint(
                None,
                (0..4).map(|j| (arcs[i][j], 1.0)).collect(),
                RelOp::Eq,
                1.0,
            )
            .unwrap();
            p.add_constraint(
                None,
                (0..4).map(|j| (arcs[j][i], -2.0)).collect(),
                RelOp::Eq,
                -2.0,
            )
            .unwrap();
            for j in 0..4 {
                if i == 0 && j > 0 || i > 0 && j == 0 {
                    costs.push((arcs[i][j], 1.0));
                } else if i == j && i > 0 {
                    costs.push((arcs[i][j], 10.0));
                } else if i > 0 && j > 0 {
                    p.add_constraint(
                        None,
                        vec![
                            (arcs[i][j], 100.0),
                            (arcs[i][i], -5.0),
                            (times[i - 1], 1.0),
                            (times[j - 1], -1.0),
                        ],
                        RelOp::Le,
                        94.0,
                    )
                    .unwrap();
                }
            }
        }
        p.set_objective(costs, 2.0).unwrap();
        p
    }

    #[test]
    fn cumulative_time_rejects_chain_and_accepts_skipping_one_job() {
        let p = toy();
        let m = RouteModel::detect(&p).unwrap();
        assert!(
            m.candidate_from_arcs(&p, &[1, 6, 11, 12]).is_none(),
            "两条各耗六单位的弧不能塞入十单位窗口"
        );
        assert!(
            m.candidate_from_arcs(&p, &[0, 6, 11, 13]).is_none(),
            "正时延环不可行"
        );
        let feasible = m.candidate_from_arcs(&p, &[1, 6, 8, 15]).unwrap();
        assert_eq!(feasible[17], 6.0);
        assert_eq!(p.obj.eval(&feasible) + p.obj_offset, 14.0);
        let greedy = m.greedy_candidate(&p, None).unwrap();
        assert_eq!(greedy, feasible, "并列时取最小弧列");
        assert!(
            m.greedy_candidate(&p, Some(Instant::now())).is_none(),
            "已用尽时间不能继续构造"
        );
    }

    #[test]
    fn renumbered_rows_columns_and_names_preserve_timing() {
        let mut p = toy();
        p.vars.reverse();
        p.constraints.reverse();
        let n = p.vars.len();
        for (j, variable) in p.vars.iter_mut().enumerate() {
            variable.name = format!("v_{j}");
        }
        for expr in p
            .constraints
            .iter_mut()
            .map(|r| &mut r.expr)
            .chain([&mut p.obj])
        {
            expr.terms.reverse();
            for (j, _) in &mut expr.terms {
                *j = n - 1 - *j;
            }
        }
        p.name = "重排后".into();
        p.validate().unwrap();
        let m = RouteModel::detect(&p).unwrap();
        let greedy = m.greedy_candidate(&p, None).unwrap();
        assert_eq!(p.obj.eval(&greedy) + p.obj_offset, 14.0);
    }

    #[test]
    fn detects_only_verified_equations_and_obeys_fixed_zero_bounds() {
        let mut p = toy();
        p.set_var_bounds(6, 0.0, 0.0).unwrap();
        let m = RouteModel::detect(&p).unwrap();
        assert!(m.candidate_from_arcs(&p, &[1, 6, 8, 15]).is_none());
        assert_eq!(m.greedy_candidate(&p, None).unwrap()[6], 0.0);
        // 将非活动行右侧从 94 收紧到 9，窗口盒无法蕴含该行。
        p.constraints
            .iter_mut()
            .find(|r| r.expr.terms.len() == 4 && r.op == RelOp::Le)
            .unwrap()
            .rhs = 9.0;
        assert!(RouteModel::detect(&p).is_none());
    }
}
