//! 布尔子句搜索。监视两个文字发现单元蕴含；发生矛盾时，沿蕴含原因消去
//! 当前层的赋值，直到只剩一个当前层文字。得到的子句在回跳后立即强制该文字。
use std::time::Instant;

#[derive(Debug, PartialEq, Eq)]
pub(super) enum Outcome {
    Sat(Vec<bool>),
    Unsat,
    Limit,
}

const NONE: usize = usize::MAX;

fn var(lit: i32) -> usize {
    lit.unsigned_abs() as usize - 1
}

fn index(lit: i32) -> usize {
    2 * var(lit) + usize::from(lit < 0)
}

fn value(assign: &[i8], lit: i32) -> i8 {
    let v = assign[var(lit)];
    if lit > 0 {
        v
    } else {
        -v
    }
}

#[derive(Clone, Copy)]
struct Watch {
    clause: usize,
    other: i32,
}

struct Clause {
    lits: Vec<i32>,
    glue: usize,
    used: usize,
}

struct Search {
    clauses: Vec<Clause>,
    watches: Vec<Vec<Watch>>,
    assign: Vec<i8>,
    phase: Vec<i8>,
    levels: Vec<usize>,
    reasons: Vec<usize>,
    trail: Vec<i32>,
    head: usize,
    starts: Vec<usize>,
    seen: Vec<bool>,
    activity: Vec<f64>,
    increment: f64,
    conflicts: usize,
    original: usize,
}

impl Search {
    fn new(n: usize) -> Self {
        Self {
            clauses: Vec::new(),
            watches: vec![Vec::new(); 2 * n],
            assign: vec![0; n],
            phase: vec![-1; n],
            levels: vec![0; n],
            reasons: vec![NONE; n],
            trail: Vec::new(),
            head: 0,
            starts: Vec::new(),
            seen: vec![false; n],
            activity: vec![0.; n],
            increment: 1.,
            conflicts: 0,
            original: 0,
        }
    }

    fn enqueue(&mut self, lit: i32, reason: usize) -> bool {
        let j = var(lit);
        let sign = if lit > 0 { 1 } else { -1 };
        if self.assign[j] != 0 {
            return self.assign[j] == sign;
        }
        self.assign[j] = sign;
        self.phase[j] = sign;
        self.levels[j] = self.starts.len();
        self.reasons[j] = reason;
        self.trail.push(lit);
        true
    }

    // 单元子句直接加入根层或回跳后的轨迹。其他子句总是监视前两个位置。
    fn store(&mut self, lits: Vec<i32>, glue: usize) -> usize {
        if lits.len() == 1 {
            return NONE;
        }
        let id = self.clauses.len();
        self.watches[index(lits[0])].push(Watch {
            clause: id,
            other: lits[1],
        });
        self.watches[index(lits[1])].push(Watch {
            clause: id,
            other: lits[0],
        });
        self.clauses.push(Clause {
            lits,
            glue,
            used: self.conflicts,
        });
        id
    }

    fn propagate(&mut self, end: Option<Instant>) -> Result<Option<usize>, ()> {
        while self.head < self.trail.len() {
            if self.head % 64 == 0 && end.is_some_and(|t| Instant::now() >= t) {
                return Err(());
            }
            let false_lit = -self.trail[self.head];
            self.head += 1;
            let w = index(false_lit);
            let mut i = 0;
            while i < self.watches[w].len() {
                let watch = self.watches[w][i];
                if value(&self.assign, watch.other) > 0 {
                    i += 1;
                    continue;
                }
                let lits = &mut self.clauses[watch.clause].lits;
                if lits[0] == false_lit {
                    lits.swap(0, 1);
                }
                debug_assert_eq!(lits[1], false_lit);
                let other = lits[0];
                let other_value = value(&self.assign, other);
                if other_value != 1 {
                    if let Some(pos) = lits[2..].iter().position(|&l| value(&self.assign, l) >= 0) {
                        lits.swap(1, pos + 2);
                        let new_w = index(lits[1]);
                        self.watches[w].swap_remove(i);
                        self.watches[new_w].push(Watch {
                            clause: watch.clause,
                            other,
                        });
                        continue;
                    }
                    if other_value < 0 {
                        return Ok(Some(watch.clause));
                    }
                    self.enqueue(other, watch.clause);
                }
                self.watches[w][i].other = other;
                i += 1;
            }
        }
        Ok(None)
    }

    fn backtrack(&mut self, level: usize) {
        if self.starts.len() <= level {
            return;
        }
        let start = self.starts[level];
        for lit in self.trail.drain(start..) {
            let j = var(lit);
            self.assign[j] = 0;
            self.reasons[j] = NONE;
        }
        self.head = start;
        self.starts.truncate(level);
    }

    fn analyze(&mut self, mut id: usize) -> (Vec<i32>, usize, usize) {
        let level = self.starts.len();
        let mut kept = vec![0];
        let mut pending = 0;
        let mut cursor = self.trail.len();
        let mut pivot = NONE;
        loop {
            self.clauses[id].used = self.conflicts;
            for &lit in &self.clauses[id].lits {
                let j = var(lit);
                if j == pivot || self.seen[j] || self.levels[j] == 0 {
                    continue;
                }
                self.seen[j] = true;
                self.activity[j] += self.increment;
                if self.levels[j] == level {
                    pending += 1;
                } else {
                    kept.push(lit);
                }
            }
            loop {
                cursor -= 1;
                pivot = var(self.trail[cursor]);
                if self.seen[pivot] {
                    break;
                }
            }
            self.seen[pivot] = false;
            pending -= 1;
            if pending == 0 {
                kept[0] = -self.trail[cursor];
                break;
            }
            id = self.reasons[pivot];
            debug_assert_ne!(id, NONE);
        }
        for &lit in &kept[1..] {
            self.seen[var(lit)] = false;
        }
        self.increment /= 0.95;
        if self.increment > 1e90 {
            for v in &mut self.activity {
                *v *= 1e-90;
            }
            self.increment *= 1e-90;
        }
        // 次高决策层是回跳目标。把该层文字放到第二个监视位置，后续回退也保持监视不变式。
        let mut back = 0;
        let mut position = 0;
        let mut levels = vec![level];
        for (i, &lit) in kept.iter().enumerate().skip(1) {
            let d = self.levels[var(lit)];
            levels.push(d);
            if d > back {
                back = d;
                position = i;
            }
        }
        if kept.len() > 1 {
            kept.swap(1, position);
        }
        levels.sort_unstable();
        levels.dedup();
        (kept, back, levels.len())
    }

    fn decision(&self) -> Option<i32> {
        let mut best = None;
        let mut score = -1.;
        for (j, &a) in self.activity.iter().enumerate() {
            if self.assign[j] == 0 && a > score {
                score = a;
                best = Some((j as i32 + 1) * i32::from(self.phase[j]));
            }
        }
        best
    }

    // 原始子句、短学习子句和轨迹当前引用的原因都保留。删除其余子句时统一重写
    // 监视表与原因索引，避免重启后累计的无用子句长期占用内存。
    fn reduce(&mut self) {
        let n = self.clauses.len();
        let mut locked = vec![false; n];
        for &id in &self.reasons {
            if id != NONE {
                locked[id] = true;
            }
        }
        let mut candidates: Vec<_> = (self.original..n)
            .filter(|&i| !locked[i] && self.clauses[i].glue > 2 && self.clauses[i].lits.len() > 2)
            .collect();
        candidates.sort_unstable_by_key(|&i| {
            let c = &self.clauses[i];
            (std::cmp::Reverse(c.glue), c.used, i)
        });
        let mut map: Vec<_> = (0..n).collect();
        for &i in &candidates[..candidates.len() / 2] {
            map[i] = NONE;
        }
        let mut i = 0;
        let mut next = 0;
        self.clauses.retain(|_| {
            let keep = map[i] != NONE;
            if keep {
                map[i] = next;
                next += 1;
            }
            i += 1;
            keep
        });
        for watches in &mut self.watches {
            watches.retain_mut(|w| {
                w.clause = map[w.clause];
                w.clause != NONE
            });
        }
        for r in &mut self.reasons {
            if *r != NONE {
                *r = map[*r];
            }
        }
    }
}

// Luby 重启间隔依次是 1, 1, 2, 1, 1, 2, 4。重启保留已学子句和变量活跃度。
fn restart_interval(mut n: usize) -> usize {
    let mut k = n.next_power_of_two();
    while n != 2 * k - 1 {
        if n < k {
            k /= 2;
        } else {
            n -= k - 1;
            k = n.next_power_of_two();
        }
    }
    k
}

// 开始搜索只计一个根节点，之后每次变量决策计一个节点。两轮目标假设共享
// nodes 的累计值；重启清除决策但不清零计数。假设在本轮内是一条根单位子句。
pub(super) fn solve(
    n: usize,
    clauses: &[Vec<i32>],
    assumption: Option<i32>,
    node_limit: u64,
    deadline: Option<Instant>,
    nodes: &mut u64,
    learn: &mut dyn FnMut(&[i32]),
) -> Outcome {
    let expired = || deadline.is_some_and(|t| Instant::now() >= t);
    let contradiction = |learn: &mut dyn FnMut(&[i32])| {
        if expired() {
            return Outcome::Limit;
        }
        learn(&[]);
        if expired() {
            Outcome::Limit
        } else {
            Outcome::Unsat
        }
    };
    if *nodes >= node_limit || expired() {
        return Outcome::Limit;
    }
    *nodes += 1;
    let mut s = Search::new(n);
    let mut marks = vec![0i8; n];
    for (idx, c) in clauses.iter().enumerate() {
        if idx % 256 == 0 && expired() {
            return Outcome::Limit;
        }
        let mut lits = Vec::with_capacity(c.len());
        let mut tautology = false;
        for &lit in c {
            let j = var(lit);
            let sign = if lit > 0 { 1 } else { -1 };
            if marks[j] == -sign {
                tautology = true;
                break;
            }
            if marks[j] == 0 {
                marks[j] = sign;
                lits.push(lit);
            }
        }
        for &l in &lits {
            marks[var(l)] = 0;
        }
        if tautology {
            continue;
        }
        if lits.is_empty() || (lits.len() == 1 && !s.enqueue(lits[0], NONE)) {
            return contradiction(learn);
        }
        for &l in &lits {
            s.activity[var(l)] += 1.;
        }
        s.store(lits, 0);
    }
    s.original = s.clauses.len();
    if assumption.is_some_and(|lit| !s.enqueue(lit, NONE)) {
        return contradiction(learn);
    }
    let mut restart = 1;
    let mut next_restart = 64;
    let mut next_reduce = 3200;
    loop {
        if expired() {
            return Outcome::Limit;
        }
        match s.propagate(deadline) {
            Err(()) => return Outcome::Limit,
            Ok(Some(conflict)) => {
                if s.starts.is_empty() {
                    return contradiction(learn);
                }
                s.conflicts += 1;
                let (lits, level, glue) = s.analyze(conflict);
                learn(&lits);
                let asserting = lits[0];
                s.backtrack(level);
                let id = s.store(lits, glue);
                s.enqueue(asserting, id);
            }
            Ok(None) => {
                if s.conflicts >= next_reduce {
                    s.reduce();
                    next_reduce += 3200;
                }
                if s.conflicts >= next_restart {
                    s.backtrack(0);
                    restart += 1;
                    next_restart = s.conflicts + 64 * restart_interval(restart);
                }
                let Some(lit) = s.decision() else {
                    return if expired() {
                        Outcome::Limit
                    } else {
                        Outcome::Sat(s.assign.iter().map(|&a| a == 1).collect())
                    };
                };
                if *nodes >= node_limit {
                    return Outcome::Limit;
                }
                *nodes += 1;
                s.starts.push(s.trail.len());
                s.enqueue(lit, NONE);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_solve(n: usize, rows: &[Vec<i32>], assumption: Option<i32>) -> Outcome {
        solve(n, rows, assumption, u64::MAX, None, &mut 0, &mut |_| {})
    }

    fn satisfies(rows: &[Vec<i32>], x: &[bool]) -> bool {
        rows.iter().all(|c| c.iter().any(|&l| x[var(l)] == (l > 0)))
    }

    // 用独立的逐行扫描核对学习步骤。假设新子句的所有文字为假，旧子句
    // 须通过反复强制其唯一剩余文字导出矛盾；最终的空子句也是同一个检查。
    fn reverse_propagation(n: usize, rows: &[Vec<i32>], learned: &[i32]) -> bool {
        let mut x = vec![None; n];
        for &l in learned {
            x[var(l)] = Some(l < 0);
        }
        loop {
            let mut changed = false;
            for c in rows {
                let mut open = None;
                let mut settled = false;
                for &l in c {
                    match x[var(l)] {
                        Some(b) if b == (l > 0) => {
                            settled = true;
                            break;
                        }
                        None if open.is_some() => {
                            settled = true;
                            break;
                        }
                        None => open = Some(l),
                        _ => {}
                    }
                }
                if !settled {
                    let Some(l) = open else {
                        return true;
                    };
                    x[var(l)] = Some(l > 0);
                    changed = true;
                }
            }
            if !changed {
                return false;
            }
        }
    }

    #[test]
    fn small_formulas_match_enumeration_and_every_learned_clause_has_a_proof() {
        let mut state = 211u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as usize
        };
        let mut proof_count = 0;
        let mut learned_nonempty = 0;
        for trial in 0..160 {
            let n = 2 + next() % 7;
            let mut rows = Vec::new();
            for _ in 0..2 + next() % (n * 6) {
                let len = 1 + next() % n.min(4);
                let mut c: Vec<i32> = Vec::new();
                while c.len() < len {
                    let j = (1 + next() % n) as i32;
                    if !c.iter().any(|&l| l.abs() == j) {
                        c.push(if next() % 2 == 0 { j } else { -j });
                    }
                }
                rows.push(c);
            }
            for assumption in [
                None,
                Some((trial % n + 1) as i32),
                Some(-((trial % n + 1) as i32)),
            ] {
                let mut oracle = rows.clone();
                if let Some(l) = assumption {
                    oracle.push(vec![l]);
                }
                let possible = (0..1usize << n).any(|m| {
                    let x: Vec<_> = (0..n).map(|j| m >> j & 1 == 1).collect();
                    satisfies(&oracle, &x)
                });
                let mut proof = oracle.clone();
                let mut closed = false;
                let r = solve(n, &rows, assumption, u64::MAX, None, &mut 0, &mut |c| {
                    assert!(reverse_propagation(n, &proof, c));
                    closed = c.is_empty();
                    learned_nonempty += usize::from(!closed);
                    proof.push(c.to_vec());
                });
                match r {
                    Outcome::Sat(x) => assert!(possible && !closed && satisfies(&oracle, &x)),
                    Outcome::Unsat => {
                        assert!(!possible && closed);
                        proof_count += 1;
                    }
                    Outcome::Limit => panic!("无预算限制的小模型应当完成搜索"),
                }
            }
        }
        assert!(proof_count > 0 && learned_nonempty > 0);
    }

    #[test]
    fn duplicates_tautologies_and_root_limits() {
        assert_eq!(small_solve(0, &[], None), Outcome::Sat(vec![]));
        assert_eq!(small_solve(0, &[vec![]], None), Outcome::Unsat);
        let rows = vec![vec![1, 1, 1], vec![-2, 2]];
        assert_eq!(small_solve(2, &rows, Some(-1)), Outcome::Unsat);
        assert!(matches!(small_solve(2,&rows,None),Outcome::Sat(x) if x[0]));
        let mut nodes = 0;
        let clauses = [vec![1, 2], vec![-1, 2]];
        assert_eq!(
            solve(2, &clauses, None, 1, None, &mut nodes, &mut |_| {}),
            Outcome::Limit
        );
        assert_eq!(nodes, 1);
        assert_eq!(
            solve(2, &clauses, None, 1, None, &mut nodes, &mut |_| {}),
            Outcome::Limit
        );
        assert_eq!(nodes, 1);
        assert_eq!(
            solve(
                2,
                &clauses,
                None,
                2,
                Some(Instant::now()),
                &mut nodes,
                &mut |_| {}
            ),
            Outcome::Limit
        );
        assert_eq!(nodes, 1);
        assert_eq!(
            solve(2, &clauses, Some(-2), 2, None, &mut nodes, &mut |_| {}),
            Outcome::Unsat
        );
        assert_eq!(nodes, 2);
    }
}
