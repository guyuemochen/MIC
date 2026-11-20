// 独立重放逆向单位传播证明。仅用标准库；每个待验子句都从空赋值重新检验。
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

struct Formula {
    clauses: Vec<Vec<i32>>,
    watches: Vec<Vec<usize>>,
    units: Vec<i32>,
    empty: bool,
    values: Vec<i8>,
    queue: Vec<i32>,
}

fn slot(lit: i32) -> usize {
    2 * (lit.unsigned_abs() as usize - 1) + usize::from(lit < 0)
}

fn value(values: &[i8], lit: i32) -> i8 {
    values[lit.unsigned_abs() as usize] * lit.signum() as i8
}

impl Formula {
    fn new(n: usize) -> Self {
        Self {
            clauses: Vec::new(),
            watches: vec![Vec::new(); 2 * n],
            units: Vec::new(),
            empty: false,
            values: vec![0; n + 1],
            queue: Vec::new(),
        }
    }

    fn add(&mut self, mut clause: Vec<i32>) {
        clause.sort_unstable();
        clause.dedup();
        match clause.as_slice() {
            [] => self.empty = true,
            [unit] => self.units.push(*unit),
            [a, b, ..] => {
                let id = self.clauses.len();
                self.watches[slot(*a)].push(id);
                self.watches[slot(*b)].push(id);
                self.clauses.push(clause);
            }
        }
    }

    // 原公式加上待验子句的逐文字取反，能否只靠单位传播推出矛盾。
    fn rup(&mut self, clause: &[i32]) -> bool {
        if self.empty {
            return true;
        }
        self.values.fill(0);
        self.queue.clear();
        self.queue.extend(self.units.iter().copied());
        self.queue.extend(clause.iter().map(|&lit| -lit));
        let mut head = 0;
        while head < self.queue.len() {
            let lit = self.queue[head];
            head += 1;
            match value(&self.values, lit) {
                -1 => return true,
                1 => continue,
                _ => self.values[lit.unsigned_abs() as usize] = lit.signum() as i8,
            }
            let w = slot(-lit);
            let mut i = 0;
            while i < self.watches[w].len() {
                let id = self.watches[w][i];
                let c = &mut self.clauses[id];
                if c[0] == -lit {
                    c.swap(0, 1);
                }
                debug_assert_eq!(c[1], -lit);
                if value(&self.values, c[0]) != 1 {
                    if let Some(k) = (2..c.len()).find(|&k| value(&self.values, c[k]) >= 0) {
                        c.swap(1, k);
                        self.watches[w].swap_remove(i);
                        self.watches[slot(c[1])].push(id);
                        continue;
                    }
                    if value(&self.values, c[0]) == -1 {
                        return true;
                    }
                    self.queue.push(c[0]);
                }
                i += 1;
            }
        }
        false
    }
}

fn parse_clause(input: &str, n: usize) -> Result<Vec<i32>, String> {
    let mut clause = input
        .split_whitespace()
        .map(|s| s.parse::<i32>().map_err(|_| format!("无效整数 {s}")))
        .collect::<Result<Vec<_>, _>>()?;
    if clause.pop() != Some(0)
        || clause
            .iter()
            .any(|l| *l == 0 || l.unsigned_abs() as usize > n)
    {
        return Err("子句越界或没有零结尾".into());
    }
    Ok(clause)
}

fn lines(path: &Path) -> io::Result<impl Iterator<Item = io::Result<String>>> {
    Ok(BufReader::new(File::open(path)?).lines())
}

fn check(n: usize, input: &Path, proof: &Path, assumption: i32) -> Result<(usize, usize), String> {
    if n == 0 || n > i32::MAX as usize || assumption.unsigned_abs() as usize > n {
        return Err("变量数或假设越界".into());
    }
    let mut formula = Formula::new(n);
    let mut source = 0;
    for line in lines(input).map_err(|e| e.to_string())? {
        let line = line.map_err(|e| e.to_string())?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('c') || line.starts_with('p') {
            continue;
        }
        formula.add(parse_clause(line, n)?);
        source += 1;
    }
    if assumption != 0 {
        formula.add(vec![assumption]);
    }
    let mut steps = 0;
    let mut finished = false;
    for line in lines(proof).map_err(|e| e.to_string())? {
        let line = line.map_err(|e| e.to_string())?;
        let line = line.trim();
        if line.starts_with('c') || line.is_empty() {
            continue;
        }
        if finished {
            return Err("空子句之后还有证明步骤".into());
        }
        let clause = parse_clause(line, n)?;
        if !formula.rup(&clause) {
            return Err(format!("第 {} 个推导不能由单位传播验证", steps + 1));
        }
        steps += 1;
        finished = clause.is_empty();
        formula.add(clause);
    }
    if !finished {
        return Err("证明缺少经核验的最终空子句".into());
    }
    Ok((source, steps))
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("用法：verify_rup 变量数 原子句文件 证明文件 假设（0 表示无）");
        std::process::exit(2);
    }
    match check(
        args[1].parse().expect("变量数须为正整数"),
        Path::new(&args[2]),
        Path::new(&args[3]),
        args[4].parse().expect("假设须为有符号整数"),
    ) {
        Ok((source, steps)) => println!("{source} {steps}"),
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn root_conflict_requires_both_directions() {
        let mut f = Formula::new(3);
        for c in [vec![1, 2, 3], vec![1, -2], vec![-1, 2], vec![-1, -2]] {
            f.add(c);
        }
        assert!(f.rup(&[-1]));
        assert!(!f.rup(&[1]));
        f.add(vec![-1]);
        assert!(!f.rup(&[]));
        assert!(f.rup(&[3]));
        f.add(vec![-3]);
        assert!(f.rup(&[]));
    }

    #[test]
    fn matches_full_scan_on_small_formulas() {
        let mut state = 17011_u64;
        let mut draw = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as usize
        };
        for _ in 0..200 {
            let n = 6;
            let mut f = Formula::new(n);
            let mut clauses = Vec::new();
            for _ in 0..18 {
                let mut c = Vec::new();
                for j in 1..=n {
                    if draw() % 3 == 0 {
                        c.push(if draw() % 2 == 0 {
                            j as i32
                        } else {
                            -(j as i32)
                        });
                    }
                }
                if !c.is_empty() {
                    f.add(c.clone());
                    clauses.push(c);
                }
            }
            for _ in 0..32 {
                let test: Vec<_> = (1..=n as i32)
                    .filter_map(|j| match draw() % 4 {
                        0 => Some(j),
                        1 => Some(-j),
                        _ => None,
                    })
                    .collect();
                let mut values = vec![0; n + 1];
                for &lit in &test {
                    values[lit.unsigned_abs() as usize] = -lit.signum() as i8;
                }
                let expected = 'propagate: loop {
                    let mut changed = false;
                    for c in &clauses {
                        if c.iter().any(|&lit| value(&values, lit) > 0) {
                            continue;
                        }
                        let mut unknown = c.iter().filter(|&&lit| value(&values, lit) == 0);
                        let Some(&lit) = unknown.next() else {
                            break 'propagate true;
                        };
                        if unknown.next().is_none() {
                            values[lit.unsigned_abs() as usize] = lit.signum() as i8;
                            changed = true;
                        }
                    }
                    if !changed {
                        break false;
                    }
                };
                assert_eq!(f.rup(&test), expected);
            }
        }
    }
}
