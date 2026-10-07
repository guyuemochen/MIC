//! 算例文本格式解析：行导向纯文本（'#' 整行注释、空白分词、大小写敏感；拒绝 JSON，
//! 避免任何解析依赖——自写 JSON 子集解析器代码量也不低于本格式解析器）。
//!
//! ## 指令总序（normative，唯一合法排布，无其他合法位置）
//! 1. `NAME <单词>`：可选，至多一次，**只能出现在 VALUES 之前**（重复或出现在 VALUES
//!    之后均为 Parse 错误带行号）。
//! 2. `VALUES <m>`：恰好一次；除 NAME 外的首条指令；`1 ≤ m ≤ 2^31`。
//! 3. `ENTITIES <n>`：恰好一次，**紧随 VALUES，VALUES 与 ENTITIES 之间不允许任何其他
//!    指令**（含 NAME、RES、CAP 及实体级指令）；`1 ≤ n ≤ 2^20`。
//! 4. 其后区域（七类任意交错均合法，唯一强约束：**所有 RES/CAP 行必须先于所有 DEMAND 行**，
//!    违者 Parse 错误带行号）：
//!     - `COSTROW <e> <c0> .. <c_{m-1}>`：整行软代价，恰好 m 个非负整数，单值 ≤ 2^40
//!     - `COST <e> <v> <c>`：单点软代价，非负 ≤ 2^40；与 COSTROW 按行序后写覆盖先写
//!     - `FORBID <e> <v>`：从实体 e 的值域删除 v
//!     - `CONF <e1> <e2>`：两实体不得同值（硬约束）；允许单向书写，模型侧对称化
//!     - `RES <name> <c0> .. <c_{m-1}>`：新增资源，恰好 m 个非负容量 ≤ 2^40；
//!       名字不可重复且不可为 "cap"
//!     - `CAP <v> <c>`：糖——首次出现创建资源 "cap"（每值容量 i64::MAX/2）并设
//!       cap[v] = c ≤ 2^40；后续 CAP 行复用同一资源
//!     - `DEMAND <e> <d0> .. <d_{R-1}>`：实体对全部 R 个资源（含 CAP 产生的、按声明顺序）
//!       的需求，元数必须等于 R，单值非负 ≤ 2^31
//!
//! ## 规模硬界（解析期以 i128 求和判定，超界信息注明是哪项超界）
//! 资源总数 R ≤ 64；`Σ_{e,v} cost ≤ 2^40`；每资源 `Σ_e demand ≤ 2^40`。
//!
//! ## 缺省与语义
//! 未声明的 cost 为 0；值域默认全值 [0, m) 减 FORBID；未写 DEMAND 的实体对 "cap" 资源
//! 需求 1（若存在）、其余资源 0；写了 DEMAND 的实体以该行为准（含 cap 维度）。
//! 解析完成后构建 [`ProblemModel`] 并调用 `validate()`（双保险）；
//! 任何实体值域为空 → `PlannerError::Model`（信息含实体 id）。
//! 实现要点：先逐行解析进中间结构 + 已见状态（NAME/VALUES/ENTITIES 是否出现、DEMAND 是否
//! 已出现），再组装模型；不得用 HashMap 迭代序参与决策。

use std::path::Path;

use crate::model::{
    PlanningEntity, PlannerError, ProblemModel, Resource, ValueId, CAP_SENTINEL,
};

/// 规模硬界（与 model::validate 同源；解析期先拦，行级错误带行号）。
const MAX_VALUES: i128 = 1 << 31;
const MAX_ENTITIES: i128 = 1 << 20;
const MAX_RESOURCES: usize = 64;
const MAX_COST_OR_CAP: i128 = 1 << 40;
const MAX_DEMAND: i128 = 1 << 31;
const MAX_TOTAL: i128 = 1 << 40;

/// 解析算例文本。
///
/// 语法 / 指令顺序 / id 越界 / 数值超界错误 → `PlannerError::Parse { line, .. }`
/// （line 为 1 起始行号）；组装期校验（如值域为空）→ `PlannerError::Model(..)`。
/// 用户输入不引发 panic。
pub fn parse(text: &str) -> Result<ProblemModel, PlannerError> {
    // 容忍文件头单个 U+FEFF BOM（剥去再解析）；行尾由 lines() 统一处理 LF / CRLF
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);

    // ---- 已见状态 + 中间结构（全部有序容器，无 HashMap 迭代序参与决策） ----
    let mut name: Option<String> = None;
    let mut values_m: Option<usize> = None;
    let mut entities_n: Option<usize> = None;
    let mut cost: Vec<Vec<i64>> = Vec::new();
    let mut forbids: Vec<Vec<ValueId>> = Vec::new();
    let mut confs: Vec<Vec<usize>> = Vec::new();
    let mut resources: Vec<Resource> = Vec::new();
    let mut cap_res_index: Option<usize> = None; // CAP 糖资源 "cap" 在 resources 中的下标
    let mut demand_rows: Vec<Option<Vec<i64>>> = Vec::new();
    let mut demand_seen = false;
    let mut last_cost_line = 0usize; // Σcost 超界时归因的最后一条 cost 指令行
    let mut last_demand_line = 0usize; // Σdemand 超界时归因的最后一条 DEMAND 行

    let total_lines = text.lines().count().max(1);

    for (idx, raw) in text.lines().enumerate() {
        let line_no = idx + 1;
        let tokens: Vec<&str> = raw.split_whitespace().collect();
        if tokens.is_empty() || tokens[0].starts_with('#') {
            continue; // 空行或整行注释
        }
        let directive = tokens[0];
        match directive {
            "NAME" => {
                if values_m.is_some() {
                    return Err(parse_err(line_no, "NAME must appear before VALUES"));
                }
                if name.is_some() {
                    return Err(parse_err(line_no, "duplicate NAME"));
                }
                if tokens.len() != 2 {
                    return Err(parse_err(line_no, "NAME expects exactly 1 argument"));
                }
                name = Some(tokens[1].to_string());
            }
            "VALUES" => {
                if values_m.is_some() {
                    return Err(parse_err(line_no, "duplicate VALUES"));
                }
                if tokens.len() != 2 {
                    return Err(parse_err(line_no, "VALUES expects exactly 1 argument"));
                }
                let m = parse_count(tokens[1], line_no, "VALUES count")?;
                if m < 1 {
                    return Err(parse_err(line_no, "VALUES count must be >= 1"));
                }
                if m > MAX_VALUES {
                    return Err(parse_err(line_no, "VALUES count exceeds 2^31"));
                }
                values_m = Some(m as usize);
            }
            "ENTITIES" => {
                if values_m.is_none() {
                    return Err(parse_err(line_no, "ENTITIES before VALUES"));
                }
                if entities_n.is_some() {
                    return Err(parse_err(line_no, "duplicate ENTITIES"));
                }
                if tokens.len() != 2 {
                    return Err(parse_err(line_no, "ENTITIES expects exactly 1 argument"));
                }
                let n = parse_count(tokens[1], line_no, "ENTITIES count")?;
                if n < 1 {
                    return Err(parse_err(line_no, "ENTITIES count must be >= 1"));
                }
                if n > MAX_ENTITIES {
                    return Err(parse_err(line_no, "ENTITIES count exceeds 2^20"));
                }
                entities_n = Some(n as usize);
                // 初始化按实体索引的中间结构（此后才允许实体级 / RES / CAP / DEMAND 指令）
                let (m, n) = (values_m.unwrap(), n as usize);
                cost = vec![vec![0i64; m]; n];
                forbids = vec![Vec::new(); n];
                confs = vec![Vec::new(); n];
                demand_rows = vec![None; n];
            }
            "COSTROW" => {
                let (m, n) = expect_entities(&values_m, &entities_n, directive, line_no)?;
                if tokens.len() != 2 + m {
                    return Err(parse_err(
                        line_no,
                        format!("COSTROW expects exactly {m} cost values, got {}", tokens.len() - 2),
                    ));
                }
                let e = parse_index(tokens[1], line_no, n, "COSTROW entity id")?;
                let mut row = Vec::with_capacity(m);
                for tok in &tokens[2..] {
                    let c = parse_count(tok, line_no, "COSTROW cost")?;
                    if c < 0 {
                        return Err(parse_err(line_no, "COSTROW cost must be non-negative"));
                    }
                    if c > MAX_COST_OR_CAP {
                        return Err(parse_err(line_no, "COSTROW cost exceeds 2^40"));
                    }
                    row.push(c as i64);
                }
                cost[e] = row;
                last_cost_line = line_no;
            }
            "COST" => {
                let (m, n) = expect_entities(&values_m, &entities_n, directive, line_no)?;
                if tokens.len() != 4 {
                    return Err(parse_err(line_no, "COST expects exactly 3 arguments"));
                }
                let e = parse_index(tokens[1], line_no, n, "COST entity id")?;
                let v = parse_index(tokens[2], line_no, m, "COST value id")?;
                let c = parse_count(tokens[3], line_no, "COST value")?;
                if c < 0 {
                    return Err(parse_err(line_no, "COST value must be non-negative"));
                }
                if c > MAX_COST_OR_CAP {
                    return Err(parse_err(line_no, "COST value exceeds 2^40"));
                }
                cost[e][v] = c as i64;
                last_cost_line = line_no;
            }
            "FORBID" => {
                let (m, n) = expect_entities(&values_m, &entities_n, directive, line_no)?;
                if tokens.len() != 3 {
                    return Err(parse_err(line_no, "FORBID expects exactly 2 arguments"));
                }
                let e = parse_index(tokens[1], line_no, n, "FORBID entity id")?;
                let v = parse_index(tokens[2], line_no, m, "FORBID value id")?;
                forbids[e].push(v as u32);
            }
            "CONF" => {
                let (_, n) = expect_entities(&values_m, &entities_n, directive, line_no)?;
                if tokens.len() != 3 {
                    return Err(parse_err(line_no, "CONF expects exactly 2 arguments"));
                }
                let e1 = parse_index(tokens[1], line_no, n, "CONF entity id")?;
                let e2 = parse_index(tokens[2], line_no, n, "CONF entity id")?;
                if e1 == e2 {
                    return Err(parse_err(line_no, "CONF self-loop is not allowed"));
                }
                confs[e1].push(e2);
                confs[e2].push(e1); // 单向书写，模型侧对称化
            }
            "RES" => {
                let (m, _) = expect_entities(&values_m, &entities_n, directive, line_no)?;
                if demand_seen {
                    return Err(parse_err(line_no, "RES after DEMAND is not allowed"));
                }
                if tokens.len() != 2 + m {
                    return Err(parse_err(
                        line_no,
                        format!("RES expects exactly {m} capacity values, got {}", tokens.len() - 2),
                    ));
                }
                let res_name = tokens[1];
                if res_name == "cap" {
                    return Err(parse_err(line_no, "RES name 'cap' is reserved for CAP sugar"));
                }
                if resources.iter().any(|r| r.name == res_name) {
                    return Err(parse_err(line_no, format!("duplicate RES name '{res_name}'")));
                }
                if resources.len() >= MAX_RESOURCES {
                    return Err(parse_err(line_no, "resource count exceeds 64"));
                }
                let mut capacity = Vec::with_capacity(m);
                for tok in &tokens[2..] {
                    let c = parse_count(tok, line_no, "RES capacity")?;
                    if c < 0 {
                        return Err(parse_err(line_no, "RES capacity must be non-negative"));
                    }
                    if c > MAX_COST_OR_CAP {
                        return Err(parse_err(line_no, "RES capacity exceeds 2^40"));
                    }
                    capacity.push(c as i64);
                }
                resources.push(Resource { name: res_name.to_string(), capacity });
            }
            "CAP" => {
                let (m, _) = expect_entities(&values_m, &entities_n, directive, line_no)?;
                if demand_seen {
                    return Err(parse_err(line_no, "CAP after DEMAND is not allowed"));
                }
                if tokens.len() != 3 {
                    return Err(parse_err(line_no, "CAP expects exactly 2 arguments"));
                }
                let v = parse_index(tokens[1], line_no, m, "CAP value id")?;
                let c = parse_count(tokens[2], line_no, "CAP capacity")?;
                if c < 0 {
                    return Err(parse_err(line_no, "CAP capacity must be non-negative"));
                }
                if c > MAX_COST_OR_CAP {
                    return Err(parse_err(line_no, "CAP capacity exceeds 2^40"));
                }
                // 记录 "cap" 资源下标而非 last_mut()：首个 CAP 之后、DEMAND 之前仍可出现
                // RES 行并追加到 resources 尾部，此时 last_mut() 会误写最后那个 RES 资源。
                let cap_at = match cap_res_index {
                    Some(i) => i,
                    None => {
                        if resources.len() >= MAX_RESOURCES {
                            return Err(parse_err(line_no, "resource count exceeds 64"));
                        }
                        resources.push(Resource {
                            name: "cap".to_string(),
                            capacity: vec![CAP_SENTINEL; m],
                        });
                        let i = resources.len() - 1;
                        cap_res_index = Some(i);
                        i
                    }
                };
                resources[cap_at].capacity[v] = c as i64;
            }
            "DEMAND" => {
                let (_, n) = expect_entities(&values_m, &entities_n, directive, line_no)?;
                let r = resources.len();
                if tokens.len() != 2 + r {
                    return Err(parse_err(
                        line_no,
                        format!("DEMAND expects exactly {r} demand values, got {}", tokens.len() - 2),
                    ));
                }
                let e = parse_index(tokens[1], line_no, n, "DEMAND entity id")?;
                let mut row = Vec::with_capacity(r);
                for tok in &tokens[2..] {
                    let d = parse_count(tok, line_no, "DEMAND value")?;
                    if d < 0 {
                        return Err(parse_err(line_no, "DEMAND value must be non-negative"));
                    }
                    if d > MAX_DEMAND {
                        return Err(parse_err(line_no, "DEMAND value exceeds 2^31"));
                    }
                    row.push(d as i64);
                }
                demand_rows[e] = Some(row);
                demand_seen = true;
                last_demand_line = line_no;
            }
            other => {
                return Err(parse_err(line_no, format!("unknown directive '{other}'")));
            }
        }
    }

    // ---- 必需指令存在性 ----
    let m = match values_m {
        Some(m) => m,
        None => {
            return Err(parse_err(total_lines, "missing VALUES instruction"));
        }
    };
    let n = match entities_n {
        Some(n) => n,
        None => {
            return Err(parse_err(total_lines, "missing ENTITIES instruction"));
        }
    };

    // ---- 组装模型 ----
    let mut entities = Vec::with_capacity(n);
    for e in 0..n {
        // 值域 = 全值 [0, m) 减 FORBID（升序天然保持；重复 FORBID 天然幂等）
        let domain: Vec<ValueId> = (0..m as u32)
            .filter(|&v| !forbids[e].contains(&v))
            .collect();
        if domain.is_empty() {
            return Err(PlannerError::Model(format!("entity {e} has empty domain")));
        }
        // 需求缺省：未写 DEMAND 的实体对 "cap" 资源需求 1（若存在）、其余 0
        let demands: Vec<i64> = match &demand_rows[e] {
            Some(row) => row.clone(),
            None => resources
                .iter()
                .map(|res| if res.name == "cap" { 1 } else { 0 })
                .collect(),
        };
        // 冲突对称化后去重（排序 = 确定性，无 HashMap）
        let mut conflicts = confs[e].clone();
        conflicts.sort_unstable();
        conflicts.dedup();
        entities.push(PlanningEntity { domain, demands, conflicts });
    }

    // ---- 解析期总和硬界（i128 求和；错误归因到最后一条相关指令行） ----
    let cost_total: i128 = cost
        .iter()
        .flat_map(|row| row.iter().map(|&c| c as i128))
        .sum();
    if cost_total > MAX_TOTAL {
        let at = if last_cost_line > 0 { last_cost_line } else { total_lines };
        return Err(parse_err(
            at,
            format!("total cost {cost_total} exceeds 2^40"),
        ));
    }
    for (r, res) in resources.iter().enumerate() {
        let demand_total: i128 = entities.iter().map(|ent| ent.demands[r] as i128).sum();
        if demand_total > MAX_TOTAL {
            let at = if last_demand_line > 0 { last_demand_line } else { total_lines };
            return Err(parse_err(
                at,
                format!("resource '{}': total demand {demand_total} exceeds 2^40", res.name),
            ));
        }
    }

    let model = ProblemModel {
        name: name.unwrap_or_default(),
        value_count: m,
        resources,
        entities,
        cost,
    };
    // 双保险：程序化组装结果同样必须过结构 + 硬界校验
    model.validate()?;
    Ok(model)
}

/// 读取文件后调用 [`parse`]。
///
/// 文件须为 UTF-8；容忍文件头单个 U+FEFF BOM（剥去再解析）与 CRLF/LF 行尾；
/// 非 UTF-8 内容 → `Err(PlannerError::Io(..))`。
pub fn parse_file(path: &Path) -> Result<ProblemModel, PlannerError> {
    let bytes = std::fs::read(path).map_err(PlannerError::Io)?;
    let text = String::from_utf8(bytes).map_err(|_| {
        PlannerError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("{}: file is not valid UTF-8", path.display()),
        ))
    })?;
    parse(&text)
}

fn parse_err(line: usize, message: impl Into<String>) -> PlannerError {
    PlannerError::Parse { line, message: message.into() }
}

/// 解析整数（含负号以给出"负值"错误而非"非法整数"）；失败为该行 Parse 错误。
fn parse_count(tok: &str, line: usize, what: &str) -> Result<i128, PlannerError> {
    match tok.parse::<i128>() {
        Ok(v) => Ok(v),
        Err(_) => Err(parse_err(line, format!("invalid integer for {what}: '{tok}'"))),
    }
}

/// 解析非负 id（实体/值下标），越界为该行 Parse 错误。
fn parse_index(tok: &str, line: usize, bound: usize, what: &str) -> Result<usize, PlannerError> {
    let v = parse_count(tok, line, what)?;
    if v < 0 {
        return Err(parse_err(line, format!("{what} must be non-negative")));
    }
    if v as i128 >= bound as i128 {
        return Err(parse_err(line, format!("{what} {v} out of range (bound {bound})")));
    }
    Ok(v as usize)
}

/// 实体级 / 资源级指令的前置：VALUES 与 ENTITIES 均已出现。
/// （由此自然排除"插入 VALUES 与 ENTITIES 之间"与"先于 ENTITIES"两类非法排布。）
fn expect_entities(
    values_m: &Option<usize>,
    entities_n: &Option<usize>,
    directive: &str,
    line: usize,
) -> Result<(usize, usize), PlannerError> {
    if entities_n.is_none() {
        if values_m.is_some() {
            return Err(parse_err(
                line,
                format!(
                    "'{directive}' is not allowed between VALUES and ENTITIES \
                     (ENTITIES must immediately follow VALUES)"
                ),
            ));
        }
        return Err(parse_err(line, format!("'{directive}' before ENTITIES")));
    }
    Ok((values_m.unwrap(), entities_n.unwrap()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{UNASSIGNED};

    fn parse_ok(text: &str) -> ProblemModel {
        match parse(text) {
            Ok(m) => m,
            Err(e) => panic!("expected parse success, got: {e}"),
        }
    }

    fn parse_err_at(text: &str) -> (usize, String) {
        match parse(text) {
            Err(PlannerError::Parse { line, message }) => (line, message),
            other => panic!("expected Parse error, got: {other:?}"),
        }
    }

    fn parse_model_err(text: &str) -> String {
        match parse(text) {
            Err(PlannerError::Model(msg)) => msg,
            other => panic!("expected Model error, got: {other:?}"),
        }
    }

    #[test]
    fn full_instance_with_every_directive() {
        let text = "\
# 排版注释与空行均被忽略

NAME demo
VALUES 3
ENTITIES 3
COSTROW 0 5 6 7
COST 1 2 1
FORBID 2 0
CONF 0 1
RES gpu 4 4 4
CAP 0 2
DEMAND 0 1 0
DEMAND 2 3 1
";
        let m = parse_ok(text);
        assert_eq!(m.name, "demo");
        assert_eq!(m.value_count, 3);
        assert_eq!(m.entity_count(), 3);
        // COSTROW + COST：后写覆盖先写（行内其余值保留）
        assert_eq!(m.cost[0], vec![5, 6, 7]);
        assert_eq!(m.cost[1], vec![0, 0, 1]);
        assert_eq!(m.cost[2], vec![0, 0, 0]);
        // FORBID：值域 = 全值减禁值
        assert_eq!(m.domain(2), &[1u32, 2]);
        // CONF 单向书写被对称化
        assert_eq!(m.entities[0].conflicts, vec![1]);
        assert_eq!(m.entities[1].conflicts, vec![0]);
        // 资源顺序：RES 先于 CAP 创建的 "cap"
        assert_eq!(m.resources.len(), 2);
        assert_eq!(m.resources[0].name, "gpu");
        assert_eq!(m.resources[1].name, "cap");
        // CAP 糖：容量哨兵 + 指明值覆盖
        assert_eq!(m.resources[1].capacity[0], 2);
        assert_eq!(m.resources[1].capacity[1], i64::MAX / 2);
        assert_eq!(m.resources[0].capacity, vec![4, 4, 4]);
        // DEMAND：写者以行为准，缺省 = cap 维度 1、其余 0
        assert_eq!(m.entities[0].demands, vec![1, 0]);
        assert_eq!(m.entities[1].demands, vec![0, 1]);
        assert_eq!(m.entities[2].demands, vec![3, 1]);
        assert!(m.validate().is_ok());
    }

    #[test]
    fn seven_directives_interleave_freely() {
        // 七类指令在 ENTITIES 之后任意交错（唯一强约束 RES/CAP 全先于 DEMAND 在此满足）
        let text = "\
VALUES 2
ENTITIES 3
COST 0 1 4
CONF 1 2
RES r0 2 2
CAP 1 1
DEMAND 0 1 1
COSTROW 2 9 9
FORBID 0 1
DEMAND 1 0 0
";
        let m = parse_ok(text);
        assert_eq!(m.entities[2].conflicts, vec![1]);
        assert_eq!(m.domain(0), &[0u32]);
        assert_eq!(m.cost[2], vec![9, 9]);
        assert_eq!(m.entities[1].demands, vec![0, 0]);
    }

    #[test]
    fn cap_sugar_default_capacity_and_demand() {
        let text = "\
VALUES 3
ENTITIES 2
CAP 2 5
";
        let m = parse_ok(text);
        assert_eq!(m.resources.len(), 1);
        assert_eq!(m.resources[0].name, "cap");
        assert_eq!(m.resources[0].capacity, vec![i64::MAX / 2, i64::MAX / 2, 5]);
        // 缺省需求：cap 维度 1
        assert_eq!(m.entities[0].demands, vec![1]);
        assert_eq!(m.entities[1].demands, vec![1]);
        // 后续 CAP 复用同一资源
        let m = parse_ok("VALUES 2\nENTITIES 1\nCAP 0 3\nCAP 1 4\n");
        assert_eq!(m.resources.len(), 1);
        assert_eq!(m.resources[0].capacity, vec![3, 4]);
    }

    #[test]
    fn cap_after_intervening_res_targets_cap_resource() {
        // 回归：首个 CAP 之后插入 RES，再写 CAP 时容量必须进 "cap" 资源，
        // 而非 resources 尾部（当时恰为该 RES 新建的资源）。
        let text = "\
VALUES 3
ENTITIES 2
CAP 0 2
RES gpu 7 7 7
CAP 2 5
";
        let m = parse_ok(text);
        assert_eq!(m.resources.len(), 2);
        assert_eq!(m.resources[0].name, "cap");
        assert_eq!(m.resources[1].name, "gpu");
        assert_eq!(m.resources[0].capacity, vec![2, i64::MAX / 2, 5]);
        assert_eq!(m.resources[1].capacity, vec![7, 7, 7]);
        assert!(m.validate().is_ok());
    }

    #[test]
    fn cost_overrides_costrow() {
        let text = "VALUES 2\nENTITIES 1\nCOST 0 1 9\nCOSTROW 0 1 2\nCOST 0 0 7\n";
        let m = parse_ok(text);
        assert_eq!(m.cost[0], vec![7, 2]);
    }

    #[test]
    fn bom_and_crlf_are_tolerated() {
        let text = "\u{feff}NAME bom\r\nVALUES 2\r\nENTITIES 1\r\nCOST 0 0 3\r\n";
        let m = parse_ok(text);
        assert_eq!(m.name, "bom");
        assert_eq!(m.cost[0][0], 3);
    }

    #[test]
    fn parse_file_roundtrip() {
        let dir = std::env::temp_dir().join("planner_fmt_tests");
        std::fs::create_dir_all(&dir).expect("create temp dir");
        let path = dir.join("ok-instance.txt");
        std::fs::write(&path, "VALUES 2\nENTITIES 1\nCAP 0 1\n").expect("write");
        let m = parse_file(&path).expect("parse file");
        assert_eq!(m.entity_count(), 1);
        // 非 UTF-8 → Io 错误
        let bad = dir.join("bad-utf8.txt");
        std::fs::write(&bad, [0xFF, 0xFE, 0x00, 0x31]).expect("write");
        match parse_file(&bad) {
            Err(PlannerError::Io(e)) => assert_eq!(e.kind(), std::io::ErrorKind::InvalidData),
            other => panic!("expected Io error, got {other:?}"),
        }
        // 不存在的文件 → Io
        match parse_file(&dir.join("no-such-file.txt")) {
            Err(PlannerError::Io(e)) => assert_eq!(e.kind(), std::io::ErrorKind::NotFound),
            other => panic!("expected Io error, got {other:?}"),
        }
    }

    // ---- 错误例：逐条断言行号与类别 ----

    #[test]
    fn err_unknown_directive() {
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nFROB 0 0\n");
        assert_eq!((line, msg.as_str()), (3, "unknown directive 'FROB'"));
    }

    #[test]
    fn err_name_duplicate_and_after_values() {
        let (line, msg) = parse_err_at("NAME a\nNAME b\nVALUES 2\nENTITIES 1\n");
        assert_eq!((line, msg.as_str()), (2, "duplicate NAME"));
        let (line, msg) = parse_err_at("VALUES 2\nNAME b\nENTITIES 1\n");
        assert_eq!((line, msg.as_str()), (2, "NAME must appear before VALUES"));
        let (line, _) = parse_err_at("VALUES 2\nENTITIES 1\nNAME b\n");
        assert_eq!(line, 3);
    }

    #[test]
    fn err_values_missing_duplicate_zero_over() {
        let (line, msg) = parse_err_at("# only a comment\nNAME x\n");
        assert_eq!(line, 2);
        assert!(msg.contains("missing VALUES"), "got: {msg}");
        let (line, _) = parse_err_at("VALUES 2\nENTITIES 1\nVALUES 3\n");
        assert_eq!(line, 3);
        let (line, msg) = parse_err_at("VALUES 0\nENTITIES 1\n");
        assert_eq!((line, msg.as_str()), (1, "VALUES count must be >= 1"));
        let (line, msg) = parse_err_at(&format!("VALUES {}\nENTITIES 1\n", (1i128 << 31) + 1));
        assert_eq!((line, msg.as_str()), (1, "VALUES count exceeds 2^31"));
        // 边界值本身合法：m = 2^31 通过该检查（后续内存压力由调用方算例规模负责）
        // —— 不在此实例化，因会真的分配 2^31 宽的行。
        let (line, msg) = parse_err_at("VALUES x\nENTITIES 1\n");
        assert_eq!(line, 1);
        assert!(msg.starts_with("invalid integer"), "got: {msg}");
    }

    #[test]
    fn err_entities_missing_duplicate_zero_over() {
        let (line, msg) = parse_err_at("VALUES 2\n");
        assert_eq!(line, 1);
        assert!(msg.contains("missing ENTITIES"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nENTITIES 2\n");
        assert_eq!((line, msg.as_str()), (3, "duplicate ENTITIES"));
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 0\n");
        assert_eq!((line, msg.as_str()), (2, "ENTITIES count must be >= 1"));
        let (line, msg) = parse_err_at(&format!("VALUES 2\nENTITIES {}\n", (1i128 << 20) + 1));
        assert_eq!((line, msg.as_str()), (2, "ENTITIES count exceeds 2^20"));
        let (line, msg) = parse_err_at("ENTITIES 1\n");
        assert_eq!((line, msg.as_str()), (1, "ENTITIES before VALUES"));
    }

    #[test]
    fn err_instruction_between_values_and_entities() {
        // 实体级 / 资源级 / NAME 插在 VALUES 与 ENTITIES 之间均错
        let (line, msg) = parse_err_at("VALUES 2\nCOSTROW 0 1 1\nENTITIES 1\n");
        assert_eq!(line, 2);
        assert!(msg.contains("between VALUES and ENTITIES"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nRES r 1 1\nENTITIES 1\n");
        assert_eq!(line, 2);
        assert!(msg.contains("between VALUES and ENTITIES"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nNAME x\nENTITIES 1\n");
        assert_eq!((line, msg.as_str()), (2, "NAME must appear before VALUES"));
    }

    #[test]
    fn err_entity_level_before_entities() {
        let (line, msg) = parse_err_at("COST 0 0 1\nVALUES 2\nENTITIES 1\n");
        assert_eq!(line, 1);
        assert!(msg.contains("before ENTITIES"), "got: {msg}");
        let (line, msg) =
            parse_err_at("VALUES 2\nENTITIES 1\nCAP 0 5\nDEMAND 0 1\nRES r 1 1\n");
        assert_eq!((line, msg.as_str()), (5, "RES after DEMAND is not allowed"));
        let (line, msg) =
            parse_err_at("VALUES 2\nENTITIES 1\nCAP 0 5\nDEMAND 0 1\nCAP 1 1\n");
        assert_eq!((line, msg.as_str()), (5, "CAP after DEMAND is not allowed"));
    }

    #[test]
    fn err_costrow_arity_and_negative() {
        let (line, msg) = parse_err_at("VALUES 3\nENTITIES 1\nCOSTROW 0 1 2\n");
        assert_eq!(line, 3);
        assert!(msg.contains("expects exactly 3 cost values"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCOSTROW 0 1 -2\n");
        assert_eq!((line, msg.as_str()), (3, "COSTROW cost must be non-negative"));
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCOST 0 1 -2\n");
        assert_eq!((line, msg.as_str()), (3, "COST value must be non-negative"));
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nRES r -1 1\n");
        assert_eq!((line, msg.as_str()), (3, "RES capacity must be non-negative"));
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCAP 0 -1\n");
        assert_eq!((line, msg.as_str()), (3, "CAP capacity must be non-negative"));
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCAP 0 5\nDEMAND 0 -1\n");
        assert_eq!((line, msg.as_str()), (4, "DEMAND value must be non-negative"));
    }

    #[test]
    fn err_single_value_over_bound() {
        let over_cost = (1i128 << 40) + 1;
        let (line, msg) =
            parse_err_at(&format!("VALUES 2\nENTITIES 1\nCOST 0 0 {over_cost}\n"));
        assert_eq!((line, msg.as_str()), (3, "COST value exceeds 2^40"));
        let (line, msg) = parse_err_at(&format!(
            "VALUES 2\nENTITIES 1\nCOSTROW 0 {over_cost} 0\n"
        ));
        assert_eq!((line, msg.as_str()), (3, "COSTROW cost exceeds 2^40"));
        let (line, msg) =
            parse_err_at(&format!("VALUES 2\nENTITIES 1\nRES r {over_cost} 0\n"));
        assert_eq!((line, msg.as_str()), (3, "RES capacity exceeds 2^40"));
        let (line, msg) = parse_err_at(&format!("VALUES 2\nENTITIES 1\nCAP 0 {over_cost}\n"));
        assert_eq!((line, msg.as_str()), (3, "CAP capacity exceeds 2^40"));
        let over_demand = (1i128 << 31) + 1;
        let (line, msg) = parse_err_at(&format!(
            "VALUES 2\nENTITIES 1\nCAP 0 5\nDEMAND 0 {over_demand}\n"
        ));
        assert_eq!((line, msg.as_str()), (4, "DEMAND value exceeds 2^31"));
        // 恰在界上合法
        parse_ok(&format!("VALUES 2\nENTITIES 1\nCOST 0 0 {}\n", 1i128 << 40));
    }

    #[test]
    fn err_total_cost_sum_over_bound() {
        // 单值均 ≤ 2^40，总和超界：归因到最后一条 cost 指令行
        let half = 1i128 << 40;
        let (line, msg) = parse_err_at(&format!(
            "VALUES 2\nENTITIES 2\nCOST 0 0 {half}\nCOST 1 1 1\n"
        ));
        assert_eq!(line, 4);
        assert!(msg.contains("total cost") && msg.contains("exceeds 2^40"), "got: {msg}");
    }

    #[test]
    fn err_total_demand_sum_over_bound() {
        // 每资源 Σdemand > 2^40（单值 ≤ 2^31）：513 个实体各 2^31（513·2^31 = 2^40 + 2^31），
        // 归因到最后一条 DEMAND 行（第 3+513 = 516 行）
        let mut text = String::from("VALUES 1\nENTITIES 513\nRES r 9\n");
        for e in 0..513 {
            text.push_str(&format!("DEMAND {e} {}\n", 1i128 << 31));
        }
        let (line, msg) = parse_err_at(&text);
        assert_eq!(line, 3 + 513);
        assert!(msg.contains("total demand") && msg.contains("exceeds 2^40"), "got: {msg}");
    }

    #[test]
    fn err_demand_arity_mismatch() {
        // R=1（cap 糖）但给了 2 个需求
        let (line, msg) =
            parse_err_at("VALUES 2\nENTITIES 1\nCAP 0 5\nDEMAND 0 1 2\n");
        assert_eq!(line, 4);
        assert!(msg.contains("expects exactly 1 demand values"), "got: {msg}");
        // R=0 时 DEMAND 必须无需求值
        let (line, _) = parse_err_at("VALUES 2\nENTITIES 1\nDEMAND 0 1\n");
        assert_eq!(line, 3);
        parse_ok("VALUES 2\nENTITIES 1\nDEMAND 0\n"); // R=0：零元合法
    }

    #[test]
    fn err_id_out_of_range() {
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCOST 1 0 3\n");
        assert_eq!(line, 3);
        assert!(msg.contains("out of range"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCOST 0 2 3\n");
        assert_eq!(line, 3);
        assert!(msg.contains("out of range"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nFORBID 0 5\n");
        assert_eq!(line, 3);
        assert!(msg.contains("out of range"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 2\nCONF 0 2\n");
        assert_eq!(line, 3);
        assert!(msg.contains("out of range"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCOSTROW 3 1 1\n");
        assert_eq!(line, 3);
        assert!(msg.contains("out of range"), "got: {msg}");
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nCAP 9 1\n");
        assert_eq!(line, 3);
        assert!(msg.contains("out of range"), "got: {msg}");
    }

    #[test]
    fn err_conf_self_loop() {
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 2\nCONF 1 1\n");
        assert_eq!((line, msg.as_str()), (3, "CONF self-loop is not allowed"));
    }

    #[test]
    fn err_forbid_empties_domain_is_model_error() {
        let msg = parse_model_err("VALUES 2\nENTITIES 1\nFORBID 0 0\nFORBID 0 1\n");
        assert!(msg.contains("empty domain") && msg.contains("entity 0"), "got: {msg}");
    }

    #[test]
    fn err_res_name_rules() {
        let (line, msg) = parse_err_at("VALUES 2\nENTITIES 1\nRES cap 1 1\n");
        assert_eq!(line, 3);
        assert!(msg.contains("reserved"), "got: {msg}");
        let (line, msg) =
            parse_err_at("VALUES 2\nENTITIES 1\nRES r 1 1\nRES r 2 2\n");
        assert_eq!(line, 4);
        assert!(msg.contains("duplicate RES name"), "got: {msg}");
    }

    #[test]
    fn err_resource_count_exceeds_64() {
        // m=1：65 行 RES 各 1 个容量；第 65 个 RES（文件第 67 行）触发超界
        let mut text = String::from("VALUES 1\nENTITIES 1\n");
        for i in 0..65 {
            text.push_str(&format!("RES r{i} 9\n"));
        }
        let (line, msg) = parse_err_at(&text);
        assert_eq!(line, 2 + 65);
        assert!(msg.contains("resource count exceeds 64"), "got: {msg}");
    }

    #[test]
    fn unassigned_is_never_a_value() {
        // m ≤ 2^31 保证值 id 与 UNASSIGNED 不冲突：解析层不接受 u32::MAX 级 id
        let (line, _) = parse_err_at(&format!(
            "VALUES 2\nENTITIES 1\nCOST 0 {} 3\n",
            u32::MAX as i128
        ));
        assert_eq!(line, 3);
        let _ = UNASSIGNED;
    }
}
