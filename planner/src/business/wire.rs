//! 本地绑定一次传入一个完整模型。数值用十进制，约束名按 UTF-8 字节十六进制编码。
//! 所有长度由实际行内数据决定；解析错误包含输入行号。

use super::{optimize, Config, ConstructionHeuristic, Expr, ListSpec, Metaheuristic, Model, Rule};
use std::{fmt::Write, str::FromStr, time::Duration};

fn number<T: FromStr>(s: &str) -> Result<T, String> {
    s.parse().map_err(|_| format!("无效数值 {s}"))
}

fn numbers<T: FromStr>(xs: &[&str]) -> Result<Vec<T>, String> {
    xs.iter().map(|x| number(x)).collect()
}

fn name(s: &str) -> Result<String, String> {
    fn digit(b: u8) -> Result<u8, String> {
        (b as char)
            .to_digit(16)
            .map(|v| v as u8)
            .ok_or_else(|| "约束名应为十六进制".into())
    }
    if !s.len().is_multiple_of(2) {
        return Err("约束名长度错误".into());
    }
    let bytes = s
        .as_bytes()
        .chunks_exact(2)
        .map(|c| Ok(digit(c[0])? * 16 + digit(c[1])?))
        .collect::<Result<Vec<_>, String>>()?;
    String::from_utf8(bytes).map_err(|_| "约束名应为 UTF-8".into())
}

fn expr(xs: &[&str]) -> Result<Expr, String> {
    use Expr::*;
    Ok(match xs {
        ["const", v] => Const(number(v)?),
        ["var", v] => Var(number(v)?),
        ["if", a, b, c] => If(number(a)?, number(b)?, number(c)?),
        ["sum", xs @ ..] => Sum(numbers(xs)?),
        ["max", xs @ ..] => Max(numbers(xs)?),
        ["lookup", idx, xs @ ..] => Lookup(numbers(xs)?, number(idx)?),
        ["load", g, o, xs @ ..] => Load {
            group: number(g)?,
            owner: number(o)?,
            demands: numbers(xs)?,
        },
        ["distance", g, o, xs @ ..] => Distance {
            group: number(g)?,
            owner: number(o)?,
            matrix: numbers(xs)?,
        },
        ["adjacent", group, owner, before, after] => Adjacent {
            group: number(group)?,
            owner: number(owner)?,
            before: number(before)?,
            after: number(after)?,
        },
        ["finish", group, item] => Finish {
            group: number(group)?,
            item: number(item)?,
        },
        [op, a, b] => {
            let (a, b) = (number(a)?, number(b)?);
            match *op {
                "add" => Add(a, b),
                "sub" => Sub(a, b),
                "mul" => Mul(a, b),
                "eq" => Eq(a, b),
                "ne" => Ne(a, b),
                "lt" => Lt(a, b),
                "le" => Le(a, b),
                "gt" => Gt(a, b),
                "ge" => Ge(a, b),
                _ => return Err(format!("未知表达式 {op}")),
            }
        }
        _ => return Err("表达式参数错误".into()),
    })
}

pub fn parse(input: &str) -> Result<(Model, Config), String> {
    let mut model = Model::default();
    let mut cfg = Config::default();
    let mut configured = false;
    let mut remaining = 0usize;
    for (line, raw) in input.lines().enumerate() {
        let tokens: Vec<_> = raw.split_ascii_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }
        let mut read = || -> Result<(), String> {
            if remaining > 0 {
                if let ["I", d, owners @ ..] = tokens.as_slice() {
                    let spec = model.lists.last_mut().ok_or("缺少列表声明")?;
                    spec.durations.push(number(d)?);
                    spec.eligible.push(numbers(owners)?);
                    remaining -= 1;
                    return Ok(());
                }
                return Err("列表项目数量不足".into());
            }
            match tokens.as_slice() {
                ["V", xs @ ..] => model.domains.push(numbers(xs)?),
                ["L", owners, n] => {
                    remaining = number(n)?;
                    model.lists.push(ListSpec {
                        owners: number(owners)?,
                        eligible: Vec::new(),
                        durations: Vec::new(),
                        precedences: Vec::new(),
                    });
                }
                ["P", g, a, b] => model
                    .lists
                    .get_mut(number::<usize>(g)?)
                    .ok_or("未知列表组")?
                    .precedences
                    .push((number(a)?, number(b)?)),
                ["E", xs @ ..] => model.nodes.push(expr(xs)?),
                ["R", n, t, w, e] => model.rules.push(Rule {
                    name: name(n)?,
                    hard: match *t {
                        "H" => true,
                        "S" => false,
                        _ => return Err("未知约束类型".into()),
                    },
                    weight: number(w)?,
                    expr: number(e)?,
                }),
                ["O", e] if model.objective.is_none() => model.objective = Some(number(e)?),
                ["C", seed, steps, ms, threads, algo, init] if !configured => {
                    cfg.seed = number(seed)?;
                    cfg.steps = number(steps)?;
                    cfg.time_limit = if *ms == "-1" {
                        None
                    } else {
                        Some(Duration::from_millis(number(ms)?))
                    };
                    cfg.threads = number(threads)?;
                    cfg.algorithm = match *algo {
                        "tabu" => Metaheuristic::TabuSearch { tabu_size: 8 },
                        "la" => Metaheuristic::LateAcceptance { history_len: 30 },
                        "sa" => Metaheuristic::SimulatedAnnealing {
                            start_temp: 1000.0,
                            end_temp: 1.0,
                        },
                        _ => return Err("未知局部搜索算法".into()),
                    };
                    cfg.construction = match *init {
                        "ff" => ConstructionHeuristic::FirstFit,
                        "bf" => ConstructionHeuristic::BestFit,
                        _ => return Err("未知构造策略".into()),
                    };
                    configured = true;
                }
                _ => return Err("指令重复、未知或参数数量不符".into()),
            }
            Ok(())
        };
        read().map_err(|e| format!("第 {} 行：{e}", line + 1))?;
    }
    if remaining > 0 {
        return Err("列表项目数量不足".into());
    }
    Ok((model, cfg))
}

/// 执行本地模型并拷出自包含结果；同一输入无共享可变状态。
pub fn run(input: &str) -> Result<String, String> {
    let (model, cfg) = parse(input)?;
    let r = optimize(&model, &cfg)?;
    let status = r.reason.status();
    let mut out = format!(
        "S {} {} {} {} {}\nT {:?}\nV",
        r.score.hard,
        r.score.soft,
        r.steps,
        status,
        u8::from(r.complete),
        r.reason
    );
    for v in r.state.values {
        match v {
            Some(x) => write!(out, " {x}").unwrap(),
            None => out.push_str(" _"),
        }
    }
    out.push('\n');
    for (g, owners) in r.state.lists.iter().enumerate() {
        for (o, route) in owners.iter().enumerate() {
            write!(out, "L {g} {o}").unwrap();
            for item in route {
                write!(out, " {item}").unwrap();
            }
            out.push('\n');
        }
    }
    for (rule, penalty) in model.rules.iter().zip(&r.score.penalties) {
        out.push_str("R ");
        for byte in rule.name.as_bytes() {
            write!(out, "{byte:02x}").unwrap();
        }
        writeln!(out, " {penalty}").unwrap();
    }
    Ok(out)
}
