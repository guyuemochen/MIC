//! .milp 人类可读算例格式的解析与序列化（WP1 契约）。
//!
//! 文法（行式解析，'#' 起注释至行尾，空白分隔 token，大小写敏感）：
//!
//! ```text
//! [name <算例名>]          # 可选；算例名=单个非空白 token，允许 '-''_''.'，
//! sense max|min            #   仅显示用途，不得含 '#'
//! vars
//!   <变量名> continuous|integer|binary   # 每行一个；变量名 [A-Za-z_][A-Za-z0-9_.]*，
//!                                        # 不重复、不得为保留字
//! obj                      # <系数> <变量名> 对，可跨多行直到下一节关键字；
//!                          #   数值为 Rust f64 语法（允许 -6、2.5、1e-2）
//! cons                     # 每行一条：[约束名:] <系数> <变量名> ... <=|>=|= <右端>
//!                          #   约束名=去尾 ':' 的非空白序列（允许 '-'）；
//!                          #   系数 0.0 的项解析层直接丢弃
//! [bounds                  # 可选；每行：<变量名> <=|>=|= <数值>，
//!   ...]                   #   叠加在类型默认界上取交集；Binary 只允许 [0,1] 内
//! end                      # 必须存在
//! ```
//!
//! 保留字（变量名撞之 → `Err(Parse)`）：name、sense、vars、obj、cons、bounds、
//! end、continuous、integer、binary、max、min。
//!
//! bounds 字面量必为有限数；默认上界 +∞ 合法存在于模型中（validate 允许
//! ub=+∞、要求 lb 有限），文法不提供无穷字面量。
//!
//! 本实现的确定性补充约定（均以 `Err(Parse{line, ..})` 报告，1 起始行号）：
//! - 一切数值字面量（obj/cons 系数、右端、bounds 值）都必须解析为有限 f64，
//!   `inf`/`NaN`/`1e999` 等一律拒绝；
//! - 节头行（vars/obj/cons/bounds/end）必须独占一行（单 token）；
//! - obj 节内、单条约束行内同一变量重复出现即报错（不做静默归并）；
//! - bounds 行与现有界交集为空（lb>ub）或 Binary 界值超出 [0,1] 即报错；
//!   整数变量非空但“取整后为空”的区间（如 [2.5, 2.7])合法通过，由预处理
//!   R1 / Simplex 按空界判定不可行；
//! - `end` 之后不得再出现任何非注释 token。

use std::collections::HashSet;
use std::path::Path;

use crate::error::{MilpError, Result};
use crate::model::{default_bounds, is_reserved_word, Problem, RelOp, Sense, VarType};

/// 节关键字：行首出现即切换节（这些词不得用作变量名，故切分无歧义）。
const SECTION_KEYWORDS: [&str; 7] = ["name", "sense", "vars", "obj", "cons", "bounds", "end"];

/// `token` 是否为节关键字。
fn is_section_keyword(token: &str) -> bool {
    SECTION_KEYWORDS.contains(&token)
}

/// 算例名/约束名字符集：非空的 ASCII 字母/数字/'-'/'_'/'.' 序列。
fn valid_instance_name(name: &str) -> bool {
    !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.')
}

/// 变量名字符集：`[A-Za-z_][A-Za-z0-9_.]*`（不含 '-'）。
fn valid_var_name(name: &str) -> bool {
    let mut chars = name.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() || c == '_' => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.')
}

/// 构造带 1 起始行号的解析错误。
fn parse_err(line: usize, message: impl Into<String>) -> MilpError {
    MilpError::Parse {
        line,
        message: message.into(),
    }
}

/// 解析 Rust f64 语法数值；不可解析或非有限数均为 `Err(Parse)`。
fn parse_number(token: &str, line: usize, what: &str) -> Result<f64> {
    match token.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(v),
        Ok(v) => Err(parse_err(
            line,
            format!("{what} `{token}` = {v} 非有限数（文法不提供无穷字面量）"),
        )),
        Err(_) => Err(parse_err(
            line,
            format!("{what} `{token}` 不是合法的 f64 数值"),
        )),
    }
}

/// 解析关系运算符 token（`<=`/`>=`/`=`）。
fn parse_rel_op(token: &str, line: usize, what: &str) -> Result<RelOp> {
    match token {
        "<=" => Ok(RelOp::Le),
        ">=" => Ok(RelOp::Ge),
        "=" => Ok(RelOp::Eq),
        other => Err(parse_err(
            line,
            format!("{what}运算符非法：`{other}`（期望 <=/>=/=）"),
        )),
    }
}

/// 数值写出：f64 `Display` 给出最短且可精确往返的十进制表示
/// （写出前已 validate，不会出现 inf/NaN）。
fn format_number(v: f64) -> String {
    format!("{v}")
}

/// 从文本解析 .milp 算例；任何 token 数不符、未知节、未知变量名、系数不可
/// 解析、变量名撞保留字 → `Err(Parse{line, ..})`（1 起始行号）。
pub fn parse_instance(text: &str) -> Result<Problem> {
    // 预处理：逐行截断 '#' 注释、按空白切分、丢弃空行；记录 1 起始行号。
    // str::lines 会剥除 \r，故 CRLF 文件同样可解析。
    let lines: Vec<(usize, Vec<&str>)> = text
        .lines()
        .enumerate()
        .filter_map(|(i, raw)| {
            let content = match raw.find('#') {
                Some(pos) => &raw[..pos],
                None => raw,
            };
            let tokens: Vec<&str> = content.split_whitespace().collect();
            if tokens.is_empty() {
                None
            } else {
                Some((i + 1, tokens))
            }
        })
        .collect();
    let last_line = lines.last().map(|&(line, _)| line).unwrap_or(1);

    /// 校验 `lines[pos]` 为独占一行的 `keyword` 节头。
    fn expect_section(
        lines: &[(usize, Vec<&str>)],
        pos: usize,
        keyword: &str,
        last_line: usize,
    ) -> Result<()> {
        match lines.get(pos) {
            Some(&(line, ref toks)) => {
                if toks[0] != keyword {
                    Err(parse_err(
                        line,
                        format!("期望 `{keyword}` 节头，得到 `{}`", toks.join(" ")),
                    ))
                } else if toks.len() != 1 {
                    Err(parse_err(
                        line,
                        format!(
                            "`{keyword}` 节头必须独占一行（多余 token：`{}`）",
                            toks[1..].join(" ")
                        ),
                    ))
                } else {
                    Ok(())
                }
            }
            None => Err(parse_err(
                last_line,
                format!("算例意外结束：缺少 `{keyword}` 节头"),
            )),
        }
    }

    let mut pos = 0usize;

    // —— name 行（可选；若存在必须是首个内容行）——
    let mut problem_name = String::new();
    if let Some(&(line, ref toks)) = lines.get(pos) {
        if toks[0] == "name" {
            if toks.len() != 2 {
                return Err(parse_err(
                    line,
                    format!(
                        "name 行 token 数不符：期望 `name <算例名>` 共 2 个，得到 {}",
                        toks.len()
                    ),
                ));
            }
            if !valid_instance_name(toks[1]) {
                return Err(parse_err(
                    line,
                    format!(
                        "算例名 `{}` 含非法字符（允许字母/数字/'-'/'_'/'.'，不得含 '#'）",
                        toks[1]
                    ),
                ));
            }
            problem_name = toks[1].to_string();
            pos += 1;
        }
    }

    // —— sense 行（必需）——
    let sense = match lines.get(pos) {
        Some(&(line, ref toks)) => {
            if toks.len() != 2 || toks[0] != "sense" {
                return Err(parse_err(
                    line,
                    format!("期望 `sense max|min` 行，得到 `{}`", toks.join(" ")),
                ));
            }
            match toks[1] {
                "max" => Sense::Maximize,
                "min" => Sense::Minimize,
                other => {
                    return Err(parse_err(
                        line,
                        format!("sense 取值非法：`{other}`（期望 max 或 min）"),
                    ))
                }
            }
        }
        None => {
            return Err(parse_err(
                last_line,
                "算例意外结束：缺少 sense 行".to_string(),
            ))
        }
    };
    pos += 1;
    let mut problem = Problem::new(&problem_name, sense);

    // —— vars 节（必需）——
    expect_section(&lines, pos, "vars", last_line)?;
    pos += 1;
    while let Some(&(line, ref toks)) = lines.get(pos) {
        if is_section_keyword(toks[0]) {
            break;
        }
        if toks.len() != 2 {
            return Err(parse_err(
                line,
                format!(
                    "变量声明 token 数不符：期望 `<变量名> continuous|integer|binary` 共 2 个，得到 {}",
                    toks.len()
                ),
            ));
        }
        let (name, type_token) = (toks[0], toks[1]);
        if !valid_var_name(name) {
            return Err(parse_err(
                line,
                format!("变量名 `{name}` 非法（须匹配 [A-Za-z_][A-Za-z0-9_.]*，不含 '-'）"),
            ));
        }
        if is_reserved_word(name) {
            return Err(parse_err(line, format!("变量名 `{name}` 撞 .milp 保留字")));
        }
        let var_type = match type_token {
            "continuous" => VarType::Continuous,
            "integer" => VarType::Integer,
            "binary" => VarType::Binary,
            other => {
                return Err(parse_err(
                    line,
                    format!("未知变量类型 `{other}`（期望 continuous/integer/binary）"),
                ))
            }
        };
        if problem.var_index(name).is_some() {
            return Err(parse_err(line, format!("变量名 `{name}` 重复声明")));
        }
        problem
            .add_var(name, var_type)
            .map_err(|e| parse_err(line, format!("添加变量失败：{e}")))?;
        pos += 1;
    }

    // —— obj 节（必需；token 流可跨行，直到下一节关键字行）——
    expect_section(&lines, pos, "obj", last_line)?;
    pos += 1;
    let mut obj_tokens: Vec<(usize, &str)> = Vec::new();
    while let Some(&(line, ref toks)) = lines.get(pos) {
        if is_section_keyword(toks[0]) {
            break;
        }
        for &token in toks {
            obj_tokens.push((line, token));
        }
        pos += 1;
    }
    if obj_tokens.len() % 2 != 0 {
        let (line, token) = obj_tokens[obj_tokens.len() - 1];
        return Err(parse_err(
            line,
            format!("obj 节 token 数为奇数（末尾 `{token}` 无法与变量名配对）"),
        ));
    }
    let mut obj_terms: Vec<(usize, f64)> = Vec::new();
    let mut obj_seen: HashSet<&str> = HashSet::new();
    for pair in obj_tokens.chunks_exact(2) {
        let (coeff_line, coeff_token) = pair[0];
        let (var_line, var_token) = pair[1];
        let coeff = parse_number(coeff_token, coeff_line, "obj 系数")?;
        let var = problem
            .var_index(var_token)
            .ok_or_else(|| parse_err(var_line, format!("obj 中未知变量名 `{var_token}`")))?;
        if !obj_seen.insert(var_token) {
            return Err(parse_err(
                var_line,
                format!("obj 中变量 `{var_token}` 重复出现"),
            ));
        }
        obj_terms.push((var, coeff));
    }
    problem
        .set_objective(obj_terms, 0.0)
        .map_err(|e| MilpError::Internal(format!("内部一致性问题：set_objective 失败：{e}")))?;

    // —— cons 节（必需；每行一条约束）——
    expect_section(&lines, pos, "cons", last_line)?;
    pos += 1;
    while let Some(&(line, ref toks)) = lines.get(pos) {
        if is_section_keyword(toks[0]) {
            break;
        }
        // 可选约束名：首 token 以 ':' 结尾即视为约束名（不可能与节关键字混淆）。
        let mut rest: &[&str] = toks;
        let mut cons_name: Option<&str> = None;
        if let Some(first) = rest.first() {
            if let Some(stripped) = first.strip_suffix(':') {
                if !valid_instance_name(stripped) {
                    return Err(parse_err(
                        line,
                        format!("约束名 `{stripped}` 非法（允许字母/数字/'-'/'_'/'.'，非空）"),
                    ));
                }
                cons_name = Some(stripped);
                rest = &rest[1..];
            }
        }
        if rest.len() < 2 {
            return Err(parse_err(
                line,
                format!(
                    "约束行 token 数不符：缺少运算符/右端（得到 {} 个 token）",
                    rest.len()
                ),
            ));
        }
        let op_token = rest[rest.len() - 2];
        let rhs_token = rest[rest.len() - 1];
        let op = parse_rel_op(op_token, line, "约束")?;
        let rhs = parse_number(rhs_token, line, "约束右端")?;
        let pair_tokens = &rest[..rest.len() - 2];
        if pair_tokens.len() % 2 != 0 {
            return Err(parse_err(line, "约束行系数/变量名未成对".to_string()));
        }
        let mut terms: Vec<(usize, f64)> = Vec::new();
        let mut row_seen: HashSet<&str> = HashSet::new();
        for pair in pair_tokens.chunks_exact(2) {
            let coeff = parse_number(pair[0], line, "约束系数")?;
            let var = problem
                .var_index(pair[1])
                .ok_or_else(|| parse_err(line, format!("约束中未知变量名 `{}`", pair[1])))?;
            if !row_seen.insert(pair[1]) {
                return Err(parse_err(
                    line,
                    format!("约束中变量 `{}` 重复出现", pair[1]),
                ));
            }
            terms.push((var, coeff));
        }
        problem
            .add_constraint(cons_name, terms, op, rhs)
            .map_err(|e| parse_err(line, format!("添加约束失败：{e}")))?;
        pos += 1;
    }

    // —— bounds 节（可选）——
    if let Some(&(line, ref toks)) = lines.get(pos) {
        if toks[0] == "bounds" {
            if toks.len() != 1 {
                return Err(parse_err(line, "`bounds` 节头必须独占一行".to_string()));
            }
            pos += 1;
            while let Some(&(line, ref toks)) = lines.get(pos) {
                if is_section_keyword(toks[0]) {
                    break;
                }
                if toks.len() != 3 {
                    return Err(parse_err(
                        line,
                        format!(
                            "bounds 行 token 数不符：期望 `<变量名> <=|>=|= <数值>` 共 3 个，得到 {}",
                            toks.len()
                        ),
                    ));
                }
                let (name, op_token, value_token) = (toks[0], toks[1], toks[2]);
                let op = parse_rel_op(op_token, line, "bounds ")?;
                let value = parse_number(value_token, line, "bounds 数值")?;
                let var = problem
                    .var_index(name)
                    .ok_or_else(|| parse_err(line, format!("bounds 中未知变量名 `{name}`")))?;
                if problem.vars[var].var_type == VarType::Binary && !(0.0 <= value && value <= 1.0)
                {
                    return Err(parse_err(
                        line,
                        format!("Binary 变量 `{name}` 的界 {value} 超出 [0,1]"),
                    ));
                }
                // 叠加在现有界上取交集。
                let (cur_lb, cur_ub) = (problem.vars[var].lb, problem.vars[var].ub);
                let (new_lb, new_ub) = match op {
                    RelOp::Le => (cur_lb, cur_ub.min(value)),
                    RelOp::Ge => (cur_lb.max(value), cur_ub),
                    RelOp::Eq => (cur_lb.max(value), cur_ub.min(value)),
                };
                problem
                    .set_var_bounds(var, new_lb, new_ub)
                    .map_err(|e| parse_err(line, format!("bounds 行施加失败：{e}")))?;
                pos += 1;
            }
        }
    }

    // —— end 行（必需）——
    expect_section(&lines, pos, "end", last_line)?;
    pos += 1;
    if let Some(&(line, ref toks)) = lines.get(pos) {
        return Err(parse_err(
            line,
            format!("end 之后仍存在内容：`{}`", toks.join(" ")),
        ));
    }

    // 解析路径已保证全部数值有限、变量名合法、terms 规范化；此处为安全网。
    problem.validate()?;
    Ok(problem)
}

/// 读入并解析指定路径的 .milp 文件；IO 错误包装为 `MilpError::Io`。
pub fn load_instance(path: &Path) -> Result<Problem> {
    let text = std::fs::read_to_string(path)?;
    parse_instance(&text)
}

/// 规范化序列化（round-trip：`parse_instance(&write_instance(p)?)` 与 `p`
/// 语义一致；ub=+∞ 的变量以省略上界行的方式写出）。
///
/// 文法无目标常数 token，故 `p.obj_offset != 0.0` 时返回
/// `Err(InvalidModel)`（明确限制，而非静默丢失）。
///
/// 规范化细节：
/// - 空算例名不写 name 行；全默认界的变量不写 bounds 行（默认上界 +∞ 即以
///   省略上界行表达），单点界（lb==ub 且非默认）写为 `x = v`；
/// - 数值用 f64 `Display`（最短可精确往返表示）；
/// - 系数为 0 的项已由 [`crate::model::LinearExpr`] 不变量在建模时丢弃，
///   故不会写出；空表达式约束写出为 `名: <= 右端` 形式，可无损读回；
/// - 变量名/约束名/算例名不符合文法字符集、或变量名重复/撞保留字时返回
///   `Err(InvalidModel)`（无法写出可解析文件，宁可显式失败）；
/// - bounds 节语义为与类型默认界（lb=0）取交集，故低于默认下界的下界
///   （如负下界）同样返回 `Err(InvalidModel)`。
pub fn write_instance(p: &Problem) -> Result<String> {
    if p.obj_offset != 0.0 {
        return Err(MilpError::InvalidModel(format!(
            "obj_offset = {} 非零：.milp 文法无目标常数 token，无法无损写出",
            p.obj_offset
        )));
    }
    if !p.name.is_empty() && !valid_instance_name(&p.name) {
        return Err(MilpError::InvalidModel(format!(
            "算例名 `{}` 含无法按文法写出的字符（允许字母/数字/'-'/'_'/'.'）",
            p.name
        )));
    }
    let mut seen_vars: HashSet<&str> = HashSet::new();
    for v in &p.vars {
        if !valid_var_name(&v.name) {
            return Err(MilpError::InvalidModel(format!(
                "变量名 `{}` 不匹配 [A-Za-z_][A-Za-z0-9_.]*，无法写出",
                v.name
            )));
        }
        if is_reserved_word(&v.name) {
            return Err(MilpError::InvalidModel(format!(
                "变量名 `{}` 撞 .milp 保留字，无法写出",
                v.name
            )));
        }
        if !seen_vars.insert(v.name.as_str()) {
            return Err(MilpError::InvalidModel(format!(
                "变量名 `{}` 重复，无法写出",
                v.name
            )));
        }
    }
    for c in &p.constraints {
        if let Some(name) = &c.name {
            if !valid_instance_name(name) {
                return Err(MilpError::InvalidModel(format!(
                    "约束名 `{name}` 含无法按文法写出的字符（允许字母/数字/'-'/'_'/'.'）"
                )));
            }
        }
    }
    // validate 保证全部系数/右端/offset 有限、下标有效、terms 规范化、界合法，
    // 从而写出的数值 token 均可被 parse_instance 读回。
    p.validate()?;
    // bounds 节语义是“与类型默认界取交集”，默认 lb=0，故低于默认下界的下界
    // （如负下界）无法按文法无损写出——显式拒绝而非静默丢失。
    for v in &p.vars {
        let (def_lb, _) = default_bounds(v.var_type);
        if v.lb < def_lb {
            return Err(MilpError::InvalidModel(format!(
                "变量 `{}` 下界 {} 低于类型默认下界 {}：bounds 节只能与默认界取交集，该下界无法按文法写出",
                v.name, v.lb, def_lb
            )));
        }
    }

    let mut out = String::new();
    if !p.name.is_empty() {
        out.push_str("name ");
        out.push_str(&p.name);
        out.push('\n');
    }
    out.push_str(match p.sense {
        Sense::Maximize => "sense max\n",
        Sense::Minimize => "sense min\n",
    });
    out.push_str("vars\n");
    for v in &p.vars {
        out.push_str("  ");
        out.push_str(&v.name);
        out.push(' ');
        out.push_str(match v.var_type {
            VarType::Continuous => "continuous",
            VarType::Integer => "integer",
            VarType::Binary => "binary",
        });
        out.push('\n');
    }
    out.push_str("obj\n");
    if !p.obj.terms.is_empty() {
        out.push_str("  ");
        for &(var, coeff) in &p.obj.terms {
            let Some(v) = p.vars.get(var) else {
                return Err(MilpError::Internal(format!(
                    "obj 变量下标 {var} 越界（validate 已通过，不应发生）"
                )));
            };
            out.push_str(&format_number(coeff));
            out.push(' ');
            out.push_str(&v.name);
            out.push(' ');
        }
        out.pop(); // 去掉行尾多余空格
        out.push('\n');
    }
    out.push_str("cons\n");
    for c in &p.constraints {
        out.push_str("  ");
        if let Some(name) = &c.name {
            out.push_str(name);
            out.push_str(": ");
        }
        for &(var, coeff) in &c.expr.terms {
            let Some(v) = p.vars.get(var) else {
                return Err(MilpError::Internal(format!(
                    "约束变量下标 {var} 越界（validate 已通过，不应发生）"
                )));
            };
            out.push_str(&format_number(coeff));
            out.push(' ');
            out.push_str(&v.name);
            out.push(' ');
        }
        out.push_str(match c.op {
            RelOp::Le => "<=",
            RelOp::Ge => ">=",
            RelOp::Eq => "=",
        });
        out.push(' ');
        out.push_str(&format_number(c.rhs));
        out.push('\n');
    }
    // bounds：只写与类型默认界不同的变量。
    let mut bounds_lines: Vec<String> = Vec::new();
    for v in &p.vars {
        let (def_lb, def_ub) = default_bounds(v.var_type);
        if v.lb == v.ub {
            if !(v.lb == def_lb && v.ub == def_ub) {
                bounds_lines.push(format!("  {} = {}", v.name, format_number(v.lb)));
            }
        } else {
            if v.lb != def_lb {
                bounds_lines.push(format!("  {} >= {}", v.name, format_number(v.lb)));
            }
            if v.ub != def_ub {
                bounds_lines.push(format!("  {} <= {}", v.name, format_number(v.ub)));
            }
        }
    }
    if !bounds_lines.is_empty() {
        out.push_str("bounds\n");
        for line in bounds_lines {
            out.push_str(&line);
            out.push('\n');
        }
    }
    out.push_str("end\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 契约示例算例（逐字取自 WP1 brief）。
    const SAMPLE: &str = "\
sense max
vars
  x1 integer
  y1 binary
obj
  5 x1 10 y1
cons
  c1: 3 x1 1 y1 <= 12
  c2: 1 x1 -2 y1 >= -4
bounds
  x1 <= 9
end
";

    /// 逐字段比较两个 Problem（f64 精确比较：写出用最短可往返表示）。
    fn assert_same_problem(a: &Problem, b: &Problem) {
        assert_eq!(a.name, b.name, "name 不一致");
        assert_eq!(a.sense, b.sense, "sense 不一致");
        assert_eq!(a.obj_offset, b.obj_offset, "obj_offset 不一致");
        assert_eq!(a.obj.terms, b.obj.terms, "obj 不一致");
        assert_eq!(a.vars.len(), b.vars.len(), "变量数不一致");
        for (va, vb) in a.vars.iter().zip(b.vars.iter()) {
            assert_eq!(va.name, vb.name, "变量名不一致");
            assert_eq!(va.var_type, vb.var_type, "变量类型不一致");
            assert_eq!(va.lb, vb.lb, "变量 {} 的 lb 不一致", va.name);
            assert_eq!(va.ub, vb.ub, "变量 {} 的 ub 不一致", va.name);
        }
        assert_eq!(a.constraints.len(), b.constraints.len(), "约束数不一致");
        for (ca, cb) in a.constraints.iter().zip(b.constraints.iter()) {
            assert_eq!(ca.name, cb.name, "约束名不一致");
            assert_eq!(ca.expr.terms, cb.expr.terms, "约束项不一致");
            assert_eq!(ca.op, cb.op, "约束运算符不一致");
            assert_eq!(ca.rhs, cb.rhs, "约束右端不一致");
        }
    }

    /// 解析契约示例：结构、bounds 交集、整数变量集合全部正确。
    #[test]
    fn parse_contract_example() {
        let p = parse_instance(SAMPLE).unwrap();
        assert_eq!(p.name, "");
        assert_eq!(p.sense, Sense::Maximize);
        assert_eq!(p.num_vars(), 2);
        assert_eq!(p.vars[0].name, "x1");
        assert_eq!(p.vars[0].var_type, VarType::Integer);
        assert_eq!(
            (p.vars[0].lb, p.vars[0].ub),
            (0.0, 9.0),
            "bounds 行 x1 <= 9"
        );
        assert_eq!(p.vars[1].var_type, VarType::Binary);
        assert_eq!((p.vars[1].lb, p.vars[1].ub), (0.0, 1.0));
        assert_eq!(p.obj.terms, vec![(0, 5.0), (1, 10.0)]);
        assert_eq!(p.obj_offset, 0.0);
        assert_eq!(p.constraints.len(), 2);
        assert_eq!(p.constraints[0].name.as_deref(), Some("c1"));
        assert_eq!(p.constraints[0].expr.terms, vec![(0, 3.0), (1, 1.0)]);
        assert_eq!(p.constraints[0].op, RelOp::Le);
        assert_eq!(p.constraints[0].rhs, 12.0);
        assert_eq!(p.constraints[1].name.as_deref(), Some("c2"));
        assert_eq!(p.constraints[1].expr.terms, vec![(0, 1.0), (1, -2.0)]);
        assert_eq!(p.constraints[1].op, RelOp::Ge);
        assert_eq!(p.constraints[1].rhs, -4.0);
        assert_eq!(p.integer_vars(), vec![0, 1]);
        assert!(p.validate().is_ok());
    }

    /// 完整 round-trip，含默认 +∞ 上界变量（以省略上界行写出）。
    #[test]
    fn round_trip_with_default_infinite_ub() {
        // 示例算例的 round-trip。
        let p = parse_instance(SAMPLE).unwrap();
        let text = write_instance(&p).unwrap();
        let q = parse_instance(&text).unwrap();
        assert_same_problem(&p, &q);
        assert_eq!(q.vars[0].ub, 9.0);

        // API 构造的模型：含 '-' 的算例名/约束名 + 默认 +∞ 上界变量。
        let mut r = Problem::new("lp-prod-3", Sense::Maximize);
        r.add_var("x", VarType::Continuous).unwrap();
        r.add_var("k", VarType::Integer).unwrap();
        r.set_var_bounds(0, 0.0, 2.5).unwrap(); // x ∈ [0, 2.5]
                                                // k 保持默认 [0, +∞)：不写出上界行。
        r.set_objective(vec![(0, 3.0), (1, 2.0)], 0.0).unwrap();
        r.add_constraint(Some("cap-1"), vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        let text = write_instance(&r).unwrap();
        assert!(text.contains("k integer"), "写出 vars 行：\n{text}");
        // 默认 +∞ 上界不写出上界行（约束行中的 `k <= 4` 不算 bounds 行）。
        assert!(
            !text.lines().any(|l| l.starts_with("  k <=")),
            "默认 +∞ 上界不写出上界行：\n{text}"
        );
        let s = parse_instance(&text).unwrap();
        assert_same_problem(&r, &s);
        assert_eq!(s.vars[1].ub, f64::INFINITY, "读回后仍为 +∞");

        // 小数/科学记数系数的最短往返表示。
        let mut t = Problem::new("t", Sense::Minimize);
        t.add_var("x", VarType::Continuous).unwrap();
        t.set_objective(vec![(0, 0.1), (0, 0.2)], 0.0).unwrap(); // 归并为 0.30000000000000004
        t.add_constraint(None, vec![(0, 1e-2)], RelOp::Eq, 2.5)
            .unwrap();
        let u = parse_instance(&write_instance(&t).unwrap()).unwrap();
        assert_same_problem(&t, &u);
    }

    /// 解析错误必须携带正确的 1 起始行号。
    #[test]
    fn error_line_numbers() {
        // 重复变量名：第 4 行。
        let text = "sense max\nvars\n  x1 integer\n  x1 integer\nobj\ncons\nend\n";
        match parse_instance(text) {
            Err(MilpError::Parse { line, .. }) => assert_eq!(line, 4),
            other => panic!("期望 Parse 错误，得到 {other:?}"),
        }
        // obj 中未知变量名：第 5 行。
        let text = "sense max\nvars\n  x1 integer\nobj\n  5 z9\ncons\nend\n";
        match parse_instance(text) {
            Err(MilpError::Parse { line, .. }) => assert_eq!(line, 5),
            other => panic!("期望 Parse 错误，得到 {other:?}"),
        }
        // 缺 sense 行：第 1 行。
        let text = "vars\n  x1 integer\nobj\ncons\nend\n";
        match parse_instance(text) {
            Err(MilpError::Parse { line, .. }) => assert_eq!(line, 1),
            other => panic!("期望 Parse 错误，得到 {other:?}"),
        }
        // 缺 end：报在最后一个内容行。
        let text = "sense min\nvars\n  x1 continuous\nobj\ncons\n";
        match parse_instance(text) {
            Err(MilpError::Parse { line, .. }) => assert_eq!(line, 5),
            other => panic!("期望 Parse 错误，得到 {other:?}"),
        }
    }

    /// 保留字作变量名一律 Parse 错误（第 3 行）。
    ///
    /// 节关键字类保留字（name/sense/vars/obj/cons/bounds/end）在行首即触发
    /// 节切换，由节序校验拦截；非节关键字保留字（max/min/continuous 等）由
    /// 变量声明路径的保留字检查拦截——两者均报 `Parse{line:3}`。
    #[test]
    fn reserved_word_var_names_rejected() {
        for bad in ["end", "obj", "name", "sense", "vars", "cons", "bounds"] {
            let text = format!("sense min\nvars\n  {bad} integer\nobj\ncons\nend\n");
            match parse_instance(&text) {
                Err(MilpError::Parse { line, .. }) => {
                    assert_eq!(line, 3, "保留字 {bad} 应在第 3 行报错")
                }
                other => panic!("保留字 {bad} 应报 Parse，得到 {other:?}"),
            }
        }
        for bad in ["max", "min", "continuous", "integer", "binary"] {
            let text = format!("sense min\nvars\n  {bad} integer\nobj\ncons\nend\n");
            match parse_instance(&text) {
                Err(MilpError::Parse { line, message }) => {
                    assert_eq!(line, 3, "保留字 {bad} 应在第 3 行报错");
                    assert!(message.contains("保留字"), "信息应说明撞保留字：{message}");
                }
                other => panic!("保留字 {bad} 应报 Parse，得到 {other:?}"),
            }
        }
    }

    /// 非法 token（系数、运算符、类型、sense、非有限数值字面量）报 Parse。
    #[test]
    fn invalid_tokens_rejected() {
        // obj 系数不可解析。
        let t = "sense min\nvars\n  x1 integer\nobj\n  1.2.3 x1\ncons\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 5, .. })
        ));
        // 约束运算符非法。
        let t = "sense min\nvars\n  x1 integer\nobj\ncons\n  1 x1 < 5\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 6, .. })
        ));
        // 未知变量类型。
        let t = "sense min\nvars\n  x1 int\nobj\ncons\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 3, .. })
        ));
        // 非有限数值字面量。
        for lit in ["inf", "NaN", "1e999", "-inf"] {
            let t = format!("sense min\nvars\n  x1 continuous\nobj\n  {lit} x1\ncons\nend\n");
            assert!(
                matches!(parse_instance(&t), Err(MilpError::Parse { .. })),
                "字面量 {lit} 应被拒绝"
            );
        }
        // sense 取值非法。
        let t = "sense maxi\nvars\n  x1 integer\nobj\ncons\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 1, .. })
        ));
        // 变量名含 '-'（不匹配变量名正则）。
        let t = "sense min\nvars\n  foo-bar integer\nobj\ncons\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 3, .. })
        ));
    }

    /// bounds 叠加在类型默认界上取交集。
    #[test]
    fn bounds_intersection_semantics() {
        let text = "\
sense min
vars
  x1 integer
  y1 binary
  z1 continuous
obj
  1 x1
cons
  c: 1 x1 >= 1
bounds
  x1 >= 2
  x1 <= 5
  y1 = 0
  z1 <= 2.5
end
";
        let p = parse_instance(text).unwrap();
        assert_eq!(
            (p.vars[0].lb, p.vars[0].ub),
            (2.0, 5.0),
            "整数变量双边界交集"
        );
        assert_eq!(
            (p.vars[1].lb, p.vars[1].ub),
            (0.0, 0.0),
            "Binary 固定为 0 合法"
        );
        assert_eq!((p.vars[2].lb, p.vars[2].ub), (0.0, 2.5));
    }

    /// bounds 交集为空与 Binary 越界均报 Parse（附行号）。
    #[test]
    fn bounds_empty_intersection_and_binary_range_rejected() {
        let t1 = "sense min\nvars\n  x1 continuous\nobj\ncons\nbounds\n  x1 <= 3\n  x1 >= 5\nend\n";
        match parse_instance(t1) {
            Err(MilpError::Parse { line, .. }) => assert_eq!(line, 8, "使交集为空的行"),
            other => panic!("期望 Parse，得到 {other:?}"),
        }
        let t2 = "sense min\nvars\n  y1 binary\nobj\ncons\nbounds\n  y1 <= 2\nend\n";
        assert!(matches!(
            parse_instance(t2),
            Err(MilpError::Parse { line: 7, .. })
        ));
        let t3 = "sense min\nvars\n  y1 binary\nobj\ncons\nbounds\n  y1 >= -1\nend\n";
        assert!(matches!(
            parse_instance(t3),
            Err(MilpError::Parse { line: 7, .. })
        ));
    }

    /// 系数为 0 的项在解析层直接丢弃；整行皆 0 成为空行照常入模型。
    #[test]
    fn zero_coefficient_terms_dropped() {
        let text = "\
sense max
vars
  x1 integer
  y1 integer
obj
  0 x1 5 y1
cons
  c1: 0 x1 3 y1 <= 7
  c2: 0 x1 <= 2
  c3: 0 x1 0 y1 >= -3
end
";
        let p = parse_instance(text).unwrap();
        assert_eq!(p.obj.terms, vec![(1, 5.0)], "obj 中 `0 x1` 不得出现");
        assert_eq!(p.constraints[0].expr.terms, vec![(1, 3.0)]);
        assert!(
            p.constraints[1].expr.terms.is_empty(),
            "单变量零系数行成空行"
        );
        assert_eq!(p.constraints[1].rhs, 2.0);
        assert!(p.constraints[2].expr.terms.is_empty(), "整行零系数成空行");
        assert_eq!(p.num_constraints(), 3, "空行照常进入模型");
        assert!(p.validate().is_ok());
        // 空行 round-trip。
        let q = parse_instance(&write_instance(&p).unwrap()).unwrap();
        assert_same_problem(&p, &q);
    }

    /// obj_offset != 0 时 write_instance 报 InvalidModel。
    #[test]
    fn write_instance_rejects_nonzero_obj_offset() {
        let mut p = Problem::new("t", Sense::Minimize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.set_objective(vec![(0, 1.0)], 2.5).unwrap();
        match write_instance(&p) {
            Err(MilpError::InvalidModel(message)) => {
                assert!(
                    message.contains("obj_offset"),
                    "信息应指明 obj_offset：{message}"
                )
            }
            other => panic!("期望 InvalidModel，得到 {other:?}"),
        }
    }

    /// 含 '-' 的算例名与约束名可解析，且含 '.' 的变量名合法。
    #[test]
    fn names_with_hyphen_parse() {
        let text = "\
name knapsack-10
sense max
vars
  x_1.5 binary
obj
  3 x_1.5
cons
  cap-1: 2 x_1.5 <= 5
end
";
        let p = parse_instance(text).unwrap();
        assert_eq!(p.name, "knapsack-10");
        assert_eq!(p.vars[0].name, "x_1.5");
        assert_eq!(p.constraints[0].name.as_deref(), Some("cap-1"));
        assert_eq!(p.constraints[0].expr.terms, vec![(0, 2.0)]);
        // round-trip 保持。
        let q = parse_instance(&write_instance(&p).unwrap()).unwrap();
        assert_same_problem(&p, &q);
    }

    /// obj 的 <系数 变量> 对可跨行；obj 节可为空。
    #[test]
    fn obj_pairs_span_lines_and_can_be_empty() {
        let text =
            "sense min\nvars\n  x1 continuous\n  x2 continuous\nobj\n  1\n  x1\n  -2.5e0 x2\ncons\nend\n";
        let p = parse_instance(text).unwrap();
        assert_eq!(p.obj.terms, vec![(0, 1.0), (1, -2.5)]);
        let empty = "sense min\nvars\n  x1 continuous\nobj\ncons\nend\n";
        let p2 = parse_instance(empty).unwrap();
        assert!(p2.obj.terms.is_empty(), "obj 节可为空");
        assert_eq!(p2.obj.eval(&[3.0]), 0.0, "空目标常数 0");
        assert_same_problem(&p2, &parse_instance(&write_instance(&p2).unwrap()).unwrap());
    }

    /// obj 节与单条约束内同一变量重复出现必须报 Parse。
    #[test]
    fn duplicate_variable_in_obj_and_constraint_rejected() {
        let t1 = "sense min\nvars\n  x1 continuous\nobj\n  1 x1 2 x1\ncons\nend\n";
        assert!(matches!(
            parse_instance(t1),
            Err(MilpError::Parse { line: 5, .. })
        ));
        let t2 = "sense min\nvars\n  x1 continuous\nobj\ncons\n  c: 1 x1 2 x1 <= 3\nend\n";
        assert!(matches!(
            parse_instance(t2),
            Err(MilpError::Parse { line: 6, .. })
        ));
    }

    /// 结构错误（缺节、未知内容、节序错乱、token 数不符）报 Parse。
    #[test]
    fn structural_errors_reported() {
        // 未知节内容出现在 cons 处。
        let t = "sense min\nvars\n  x1 continuous\nobj\ncons\nfoo\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 6, .. })
        ));
        // end 之后仍有内容。
        let t = "sense min\nvars\n  x1 continuous\nobj\ncons\nend\nextra\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 7, .. })
        ));
        // name 行不在首位。
        let t = "sense min\nname late\nvars\n  x1 continuous\nobj\ncons\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 2, .. })
        ));
        // name 行 token 数不符。
        let t = "name a b\nsense min\nvars\n  x1 continuous\nobj\ncons\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 1, .. })
        ));
        // 变量声明 token 数不符。
        let t = "sense min\nvars\n  x1 integer extra\nobj\ncons\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 3, .. })
        ));
        // bounds 行 token 数不符。
        let t = "sense min\nvars\n  x1 continuous\nobj\ncons\nbounds\n  x1 <=\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 7, .. })
        ));
        // 节序错乱：obj 后直接 bounds。
        let t = "sense min\nvars\n  x1 continuous\nobj\nbounds\n  x1 <= 1\nend\n";
        assert!(matches!(
            parse_instance(t),
            Err(MilpError::Parse { line: 5, .. })
        ));
        // 空文本与纯注释。
        assert!(matches!(
            parse_instance(""),
            Err(MilpError::Parse { line: 1, .. })
        ));
        assert!(matches!(
            parse_instance("# 仅注释\n\n"),
            Err(MilpError::Parse { line: 1, .. })
        ));
    }

    /// 注释（整行/行尾/节内）与空行被忽略。
    #[test]
    fn comments_and_blank_lines_tolerated() {
        let text = "\
# 全行注释
name t # 行尾注释
sense max   # 注释
vars
  x1 binary # 注释

obj
  4 x1 # 注释
cons
  c1: 1 x1 <= 1 # 注释
bounds
  # 空注释行
end
";
        let p = parse_instance(text).unwrap();
        assert_eq!(p.name, "t");
        assert_eq!(p.num_vars(), 1);
        assert_eq!(p.obj.terms, vec![(0, 4.0)]);
        assert_eq!(p.constraints.len(), 1);
    }

    /// CRLF 行尾可解析（str::lines 剥除 \r）。
    #[test]
    fn crlf_line_endings_tolerated() {
        let text = "sense max\r\nvars\r\n  x1 binary\r\nobj\r\n  1 x1\r\ncons\r\nend\r\n";
        let p = parse_instance(text).unwrap();
        assert_eq!(p.vars[0].name, "x1", "\\r 不得混入 token");
    }

    /// 全默认界不输出 bounds 节；单点界写为 `x = v` 并可读回。
    #[test]
    fn write_instance_omits_default_bounds_section() {
        let mut p = Problem::new("t", Sense::Maximize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.add_var("y", VarType::Binary).unwrap();
        p.set_objective(vec![(0, 1.0), (1, 2.0)], 0.0).unwrap();
        let text = write_instance(&p).unwrap();
        assert!(
            !text.contains("bounds"),
            "全默认界不应输出 bounds 节：\n{text}"
        );
        assert_same_problem(&p, &parse_instance(&text).unwrap());

        // Binary 固定为 1：写出 `y = 1`（默认为 [0,1]，单点界非默认）。
        p.set_var_bounds(1, 1.0, 1.0).unwrap();
        let text2 = write_instance(&p).unwrap();
        assert!(text2.contains("y = 1"), "单点界应写为等式行：\n{text2}");
        assert_same_problem(&p, &parse_instance(&text2).unwrap());

        // 有限双边界：分别写 `>=` 与 `<=` 行。
        p.set_var_bounds(0, 0.5, 3.0).unwrap();
        let text3 = write_instance(&p).unwrap();
        assert!(
            text3.contains("x >= 0.5") && text3.contains("x <= 3"),
            "\n{text3}"
        );
        assert_same_problem(&p, &parse_instance(&text3).unwrap());
    }

    /// 低于类型默认下界（默认 lb=0）的下界无法按文法写出，write_instance 显式拒绝。
    #[test]
    fn write_instance_rejects_lb_below_default() {
        let mut p = Problem::new("t", Sense::Minimize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        p.set_var_bounds(0, -1.5, 3.0).unwrap();
        match write_instance(&p) {
            Err(MilpError::InvalidModel(message)) => {
                assert!(
                    message.contains("默认下界"),
                    "信息应说明默认下界限制：{message}"
                )
            }
            other => panic!("期望 InvalidModel，得到 {other:?}"),
        }
    }

    /// 无 name 行算例 round-trip（写出同样无 name 行）。
    #[test]
    fn no_name_line_round_trip() {
        let text = "sense min\nvars\n  x1 continuous\nobj\n  2 x1\ncons\n  1 x1 >= 1\nend\n";
        let p = parse_instance(text).unwrap();
        assert_eq!(p.name, "");
        assert!(p.constraints[0].name.is_none(), "无约束名");
        let written = write_instance(&p).unwrap();
        assert!(!written.starts_with("name"));
        assert_same_problem(&p, &parse_instance(&written).unwrap());
    }

    /// write_instance 拒绝无法按文法写出的名字。
    #[test]
    fn write_instance_rejects_unserializable_names() {
        let mut p = Problem::new("t", Sense::Minimize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
        let mut q = p.clone();
        q.vars[0].name = "foo-bar".to_string();
        assert!(matches!(
            write_instance(&q),
            Err(MilpError::InvalidModel(_))
        ));
        let mut r = p.clone();
        r.add_constraint(Some("c:1"), vec![(0, 1.0)], RelOp::Le, 1.0)
            .unwrap();
        assert!(matches!(
            write_instance(&r),
            Err(MilpError::InvalidModel(_))
        ));
    }

    /// load_instance 从文件读入；文件不存在报 Io。
    #[test]
    fn load_instance_file_io() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("wp1_load_test.milp");
        std::fs::write(&path, SAMPLE).unwrap();
        let p = load_instance(&path).unwrap();
        assert_eq!(p.sense, Sense::Maximize);
        assert_eq!(p.num_vars(), 2);
        let _ = std::fs::remove_file(&path);

        let missing = dir.join("definitely_missing_9x7.milp");
        assert!(matches!(load_instance(&missing), Err(MilpError::Io(_))));
    }
}
