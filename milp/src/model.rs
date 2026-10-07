//! 问题建模与集中容差定义（WP1 契约）。
//!
//! `Problem` 为不可变数据核心：分支定界只改 [`crate::simplex::Bounds`]，
//! 不改本结构。零系数项由所有构造入口强制丢弃（见 [`LinearExpr`]）。

use crate::error::{MilpError, Result};

/// 变量类型。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarType {
    /// 连续变量。
    Continuous,
    /// 0-1 变量（界固定 [0,1]）。
    Binary,
    /// 一般整数变量。
    Integer,
}

/// 目标优化方向。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sense {
    /// 最小化。
    Minimize,
    /// 最大化。
    Maximize,
}

/// 约束关系运算符。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RelOp {
    /// 小于等于（≤）。
    Le,
    /// 大于等于（≥）。
    Ge,
    /// 等于（=）。
    Eq,
}

/// 求解终态（诚实性铁律：只有满足最优性终止条件才报 Optimal；
/// 触发迭代/节点/时间上限一律报 LimitReached）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveStatus {
    /// 已证明最优。
    Optimal,
    /// 不可行。
    Infeasible,
    /// 无界。
    Unbounded,
    /// 触发上限（附已有 incumbent，绝不冒充最优）。
    LimitReached,
}

impl SolveStatus {
    /// 是否为 Optimal。
    pub fn is_optimal(&self) -> bool {
        matches!(self, SolveStatus::Optimal)
    }
}

/// .milp 文法保留字（不得用作变量名；[`Problem::add_var`] 与
/// [`crate::instance`] 解析层共同强制）。
pub(crate) const RESERVED_WORDS: [&str; 12] = [
    "name",
    "sense",
    "vars",
    "obj",
    "cons",
    "bounds",
    "end",
    "continuous",
    "integer",
    "binary",
    "max",
    "min",
];

/// 判断 `word` 是否为 .milp 文法保留字。
pub(crate) fn is_reserved_word(word: &str) -> bool {
    RESERVED_WORDS.contains(&word)
}

/// 变量类型的默认界 `(lb, ub)`：Continuous/Integer 为 `[0, +∞)`，
/// Binary 为 `[0, 1]`。
pub(crate) fn default_bounds(var_type: VarType) -> (f64, f64) {
    match var_type {
        VarType::Continuous | VarType::Integer => (0.0, f64::INFINITY),
        VarType::Binary => (0.0, 1.0),
    }
}

/// 模型变量。
#[derive(Debug, Clone)]
pub struct Variable {
    /// 变量名（匹配 [A-Za-z_][A-Za-z0-9_.]*，不含 '-'，不得为保留字）。
    pub name: String,
    /// 变量类型。
    pub var_type: VarType,
    /// 下界（必须有限；本 crate 求解路径不支持 -∞）。
    pub lb: f64,
    /// 上界（允许 f64::INFINITY 或有限数）。
    pub ub: f64,
}

/// 稀疏线性表达式。
///
/// 不变量：`terms` 按变量下标严格升序、无重复下标、无零系数；由所有构造入口
/// 强制（`add_term`/`add_constraint`/`set_objective` 丢弃 `coeff == 0.0` 的项
/// 并对同下标归并求和，归并结果为 0.0 也删除）。
#[derive(Debug, Clone, Default)]
pub struct LinearExpr {
    /// (变量下标, 系数) 对，按下标严格升序。
    pub terms: Vec<(usize, f64)>,
}

impl LinearExpr {
    /// 空表达式。
    pub fn new() -> Self {
        LinearExpr { terms: Vec::new() }
    }

    /// 追加一项：`coeff == 0.0` 直接丢弃；同下标归并求和（归并为 0.0 则删除该项）；
    /// 保持按下标严格升序。
    pub fn add_term(&mut self, var: usize, coeff: f64) {
        if coeff == 0.0 {
            return; // +0.0 与 -0.0 一并丢弃（== 比较对二者同为真）
        }
        match self.terms.binary_search_by_key(&var, |&(v, _)| v) {
            Ok(pos) => {
                let merged = self.terms[pos].1 + coeff;
                if merged == 0.0 {
                    self.terms.remove(pos);
                } else {
                    self.terms[pos].1 = merged;
                }
            }
            Err(pos) => self.terms.insert(pos, (var, coeff)),
        }
    }

    /// 在取值 `x` 上求值：`Σ coeff_j * x[j]`。
    ///
    /// `x` 长度不足时返回 NaN（以值示警，不 panic）；NaN 在浮点加法中粘滞。
    pub fn eval(&self, x: &[f64]) -> f64 {
        let mut sum = 0.0;
        for &(var, coeff) in &self.terms {
            let value = match x.get(var) {
                Some(v) => *v,
                None => return f64::NAN,
            };
            sum += coeff * value;
        }
        sum
    }
}

/// 线性约束：`expr (op) rhs`。
#[derive(Debug, Clone)]
pub struct Constraint {
    /// 约束名（可空）。
    pub name: Option<String>,
    /// 约束左端稀疏表达式。
    pub expr: LinearExpr,
    /// 关系运算符。
    pub op: RelOp,
    /// 右端（必须有限）。
    pub rhs: f64,
}

/// MILP 问题模型（不可变数据核心）。
#[derive(Debug, Clone)]
pub struct Problem {
    /// 算例名（仅显示用途）。
    pub name: String,
    /// 目标方向。
    pub sense: Sense,
    /// 变量表（下标即列号）。
    pub vars: Vec<Variable>,
    /// 目标函数（稀疏）。
    pub obj: LinearExpr,
    /// 目标常数偏移（.milp 文法无目标常数 token，解析结果恒为 0.0）。
    pub obj_offset: f64,
    /// 约束表（下标即行号）。
    pub constraints: Vec<Constraint>,
}

impl Problem {
    /// 创建空问题。
    pub fn new(name: &str, sense: Sense) -> Self {
        Problem {
            name: name.to_string(),
            sense,
            vars: Vec::new(),
            obj: LinearExpr::new(),
            obj_offset: 0.0,
            constraints: Vec::new(),
        }
    }

    /// 添加变量，返回其下标。
    ///
    /// 默认界：Continuous/Integer 为 [0,+∞)，Binary 为 [0,1]；
    /// 变量重名或撞保留字（见 .milp 文法）返回 `Err(InvalidModel)`。
    pub fn add_var(&mut self, name: &str, var_type: VarType) -> Result<usize> {
        if is_reserved_word(name) {
            return Err(MilpError::InvalidModel(format!(
                "变量名 `{name}` 撞 .milp 文法保留字"
            )));
        }
        if self.var_index(name).is_some() {
            return Err(MilpError::InvalidModel(format!("变量名 `{name}` 重复")));
        }
        let (lb, ub) = default_bounds(var_type);
        let index = self.vars.len();
        self.vars.push(Variable {
            name: name.to_string(),
            var_type,
            lb,
            ub,
        });
        Ok(index)
    }

    /// 设置变量上下界。
    ///
    /// 模型层宽松：允许 `ub = f64::INFINITY`（lb≤ub 校验仍执行）；
    /// Binary 的界只允许 [0,1] 内；lb>ub 返回 `Err(InvalidModel)`。
    pub fn set_var_bounds(&mut self, var: usize, lb: f64, ub: f64) -> Result<()> {
        let variable = self.vars.get(var).ok_or_else(|| {
            MilpError::InvalidModel(format!("变量下标 {var} 越界（变量数 {}）", self.vars.len()))
        })?;
        if !(lb <= ub) {
            // 同时覆盖 lb>ub 与 lb/ub 为 NaN 的情形（NaN 的任何比较均为假）
            return Err(MilpError::InvalidModel(format!(
                "变量 `{}` 的界非法：lb = {lb}，ub = {ub}（要求 lb ≤ ub）",
                variable.name
            )));
        }
        if variable.var_type == VarType::Binary && !(0.0 <= lb && ub <= 1.0) {
            return Err(MilpError::InvalidModel(format!(
                "Binary 变量 `{}` 的界 [{lb}, {ub}] 超出 [0,1]",
                variable.name
            )));
        }
        let variable = &mut self.vars[var];
        variable.lb = lb;
        variable.ub = ub;
        Ok(())
    }

    /// 设置目标（稀疏项 + 常数偏移）；下标越界返回 `Err(InvalidModel)`；
    /// 零系数项按 [`LinearExpr`] 不变量丢弃。
    pub fn set_objective(&mut self, terms: Vec<(usize, f64)>, offset: f64) -> Result<()> {
        let num_vars = self.vars.len();
        for &(var, _) in &terms {
            if var >= num_vars {
                return Err(MilpError::InvalidModel(format!(
                    "目标项变量下标 {var} 越界（变量数 {num_vars}）"
                )));
            }
        }
        let mut expr = LinearExpr::new();
        for (var, coeff) in terms {
            expr.add_term(var, coeff);
        }
        self.obj = expr;
        self.obj_offset = offset;
        Ok(())
    }

    /// 添加约束，返回其下标。
    ///
    /// `terms` 未排序/含重复下标时归并排序；下标越界返回 `Err(InvalidModel)`；
    /// 零系数项按 [`LinearExpr`] 不变量丢弃。
    pub fn add_constraint(
        &mut self,
        name: Option<&str>,
        terms: Vec<(usize, f64)>,
        op: RelOp,
        rhs: f64,
    ) -> Result<usize> {
        let num_vars = self.vars.len();
        for &(var, _) in &terms {
            if var >= num_vars {
                return Err(MilpError::InvalidModel(format!(
                    "约束项变量下标 {var} 越界（变量数 {num_vars}）"
                )));
            }
        }
        let mut expr = LinearExpr::new();
        for (var, coeff) in terms {
            expr.add_term(var, coeff);
        }
        let index = self.constraints.len();
        self.constraints.push(Constraint {
            name: name.map(str::to_string),
            expr,
            op,
            rhs,
        });
        Ok(index)
    }

    /// 按名查变量下标；不存在返回 `None`。
    pub fn var_index(&self, name: &str) -> Option<usize> {
        self.vars.iter().position(|v| v.name == name)
    }

    /// 变量数。
    pub fn num_vars(&self) -> usize {
        self.vars.len()
    }

    /// 约束数。
    pub fn num_constraints(&self) -> usize {
        self.constraints.len()
    }

    /// Binary + Integer 的全部下标（升序）。
    pub fn integer_vars(&self) -> Vec<usize> {
        self.vars
            .iter()
            .enumerate()
            .filter(|(_, v)| matches!(v.var_type, VarType::Binary | VarType::Integer))
            .map(|(i, _)| i)
            .collect()
    }

    /// 模型校验（与默认 +∞ 上界自洽）：
    ///
    /// - 下标有效；
    /// - 全部系数/右端/obj_offset 有限（`is_finite`）；
    /// - 每个变量 lb 有限（`is_finite`）、ub 为有限数或 `f64::INFINITY`
    ///   （其余 ±∞ 均非法）、lb ≤ ub；
    /// - terms 规范化（升序、无重复、无零系数）。
    ///
    /// 注意：不因默认 ub=+∞ 报错——无 bounds 节的 Continuous/Integer 变量
    /// （如 lp-prod-3、milp-wolsey-2）必须通过校验。
    pub fn validate(&self) -> Result<()> {
        if !self.obj_offset.is_finite() {
            return Err(MilpError::InvalidModel(format!(
                "目标常数 obj_offset = {} 非有限数",
                self.obj_offset
            )));
        }
        validate_terms(&self.obj, self.vars.len(), "目标")?;
        for (i, v) in self.vars.iter().enumerate() {
            if !v.lb.is_finite() {
                return Err(MilpError::InvalidModel(format!(
                    "变量 `{}`（下标 {i}）下界 {} 非有限数（本 crate 不支持 -∞ 下界）",
                    v.name, v.lb
                )));
            }
            let ub_legal = v.ub.is_finite() || v.ub == f64::INFINITY;
            if !ub_legal {
                return Err(MilpError::InvalidModel(format!(
                    "变量 `{}`（下标 {i}）上界 {} 非法（须为有限数或 +∞）",
                    v.name, v.ub
                )));
            }
            if !(v.lb <= v.ub) {
                return Err(MilpError::InvalidModel(format!(
                    "变量 `{}`（下标 {i}）lb = {} > ub = {}",
                    v.name, v.lb, v.ub
                )));
            }
        }
        for (row, c) in self.constraints.iter().enumerate() {
            if !c.rhs.is_finite() {
                return Err(MilpError::InvalidModel(format!(
                    "约束 {row}（{}）右端 {} 非有限数",
                    c.name.as_deref().unwrap_or("<未命名>"),
                    c.rhs
                )));
            }
            validate_terms(&c.expr, self.vars.len(), &format!("约束 {row}"))?;
        }
        Ok(())
    }
}

/// 校验一条稀疏表达式：下标有效、系数有限且非零、下标严格升序（蕴含无重复）。
fn validate_terms(expr: &LinearExpr, num_vars: usize, what: &str) -> Result<()> {
    for window in expr.terms.windows(2) {
        if window[0].0 >= window[1].0 {
            return Err(MilpError::InvalidModel(format!(
                "{what} 的 terms 未按下标严格升序（{} 后跟 {}）",
                window[0].0, window[1].0
            )));
        }
    }
    for &(var, coeff) in &expr.terms {
        if var >= num_vars {
            return Err(MilpError::InvalidModel(format!(
                "{what} 的变量下标 {var} 越界（变量数 {num_vars}）"
            )));
        }
        if !coeff.is_finite() {
            return Err(MilpError::InvalidModel(format!(
                "{what} 中变量 {var} 的系数 {coeff} 非有限数"
            )));
        }
        if coeff == 0.0 {
            return Err(MilpError::InvalidModel(format!(
                "{what} 中变量 {var} 的系数为 0（违反 LinearExpr 不变量）"
            )));
        }
    }
    Ok(())
}

/// 集中数值容差（库内禁止散落 magic number）。
#[derive(Debug, Clone, Copy)]
pub struct Tolerances {
    /// 原始/对偶可行性容差（缩放后）。默认 1e-9。
    pub feasibility: f64,
    /// 检验数最优性容差。默认 1e-9。
    pub optimality: f64,
    /// 主元最小绝对值。默认 1e-9。
    pub pivot: f64,
    /// 整数变量判定容差。默认 1e-6。
    pub integrality: f64,
    /// 系数清零阈值。默认 1e-12。
    pub zero: f64,
}

impl Default for Tolerances {
    /// 默认值：feasibility=optimality=pivot=1e-9、integrality=1e-6、zero=1e-12。
    fn default() -> Self {
        Tolerances {
            feasibility: 1e-9,
            optimality: 1e-9,
            pivot: 1e-9,
            integrality: 1e-6,
            zero: 1e-12,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 默认容差逐项精确符合契约。
    #[test]
    fn tolerances_default_values() {
        let t = Tolerances::default();
        assert_eq!(t.feasibility, 1e-9);
        assert_eq!(t.optimality, 1e-9);
        assert_eq!(t.pivot, 1e-9);
        assert_eq!(t.integrality, 1e-6);
        assert_eq!(t.zero, 1e-12);
    }

    /// is_optimal 仅对 Optimal 为真。
    #[test]
    fn solve_status_is_optimal() {
        assert!(SolveStatus::Optimal.is_optimal());
        assert!(!SolveStatus::Infeasible.is_optimal());
        assert!(!SolveStatus::Unbounded.is_optimal());
        assert!(!SolveStatus::LimitReached.is_optimal());
    }

    /// 乱序插入仍保持下标严格升序。
    #[test]
    fn linear_expr_insert_keeps_ascending_order() {
        let mut e = LinearExpr::new();
        e.add_term(2, 1.0);
        e.add_term(0, 3.0);
        e.add_term(1, -1.5);
        assert_eq!(e.terms, vec![(0, 3.0), (1, -1.5), (2, 1.0)]);
        assert!(LinearExpr::default().terms.is_empty());
    }

    /// ±0 系数直接丢弃；同下标归并；归并结果为 0 也删除。
    #[test]
    fn linear_expr_drops_zero_and_merges_to_zero() {
        let mut e = LinearExpr::new();
        e.add_term(0, 0.0);
        e.add_term(1, -0.0);
        assert!(e.terms.is_empty(), "±0 系数直接丢弃");
        e.add_term(0, 2.5);
        e.add_term(0, 1.5);
        assert_eq!(e.terms, vec![(0, 4.0)], "同下标归并求和");
        e.add_term(0, -4.0);
        assert!(e.terms.is_empty(), "归并结果为 0 也删除");
    }

    /// eval 正确求值；空表达式为 0；x 长度不足返回 NaN 而非 panic。
    #[test]
    fn linear_expr_eval() {
        let mut e = LinearExpr::new();
        e.add_term(0, 2.0);
        e.add_term(2, -1.0);
        assert_eq!(e.eval(&[1.0, 5.0, 3.0]), -1.0);
        assert_eq!(LinearExpr::new().eval(&[]), 0.0);
        assert!(e.eval(&[1.0, 1.0]).is_nan(), "x 长度不足返回 NaN");
    }

    /// 构造 x/y/z（Integer/Binary/Continuous）的测试问题。
    fn small_problem() -> Problem {
        let mut p = Problem::new("t", Sense::Maximize);
        assert_eq!(p.num_vars(), 0);
        assert!(p.add_var("x", VarType::Integer).is_ok());
        assert!(p.add_var("y", VarType::Binary).is_ok());
        assert!(p.add_var("z", VarType::Continuous).is_ok());
        p
    }

    /// 默认界、返回下标、var_index 查询。
    #[test]
    fn add_var_defaults_and_lookup() {
        let p = small_problem();
        assert_eq!(p.num_vars(), 3);
        assert_eq!(p.vars[0].name, "x");
        assert_eq!(p.vars[0].var_type, VarType::Integer);
        assert_eq!((p.vars[0].lb, p.vars[0].ub), (0.0, f64::INFINITY));
        assert_eq!((p.vars[1].lb, p.vars[1].ub), (0.0, 1.0));
        assert_eq!((p.vars[2].lb, p.vars[2].ub), (0.0, f64::INFINITY));
        assert_eq!(p.var_index("x"), Some(0));
        assert_eq!(p.var_index("z"), Some(2));
        assert_eq!(p.var_index("nope"), None);
    }

    /// 重名与保留字返回 InvalidModel。
    #[test]
    fn add_var_rejects_duplicate_and_reserved() {
        let mut p = small_problem();
        assert!(matches!(
            p.add_var("x", VarType::Continuous),
            Err(MilpError::InvalidModel(_))
        ));
        for word in ["end", "obj", "max", "min", "bounds", "continuous"] {
            assert!(
                matches!(
                    p.add_var(word, VarType::Binary),
                    Err(MilpError::InvalidModel(_))
                ),
                "保留字 {word} 应被拒绝"
            );
        }
    }

    /// 界设置规则：+∞ 上界合法、lb>ub/NaN 拒绝、Binary 仅 [0,1]。
    #[test]
    fn set_var_bounds_rules() {
        let mut p = small_problem();
        assert!(
            p.set_var_bounds(0, 0.0, f64::INFINITY).is_ok(),
            "ub=+∞ 合法"
        );
        assert!(p.set_var_bounds(0, -3.5, 7.0).is_ok(), "有限界合法");
        assert!(matches!(
            p.set_var_bounds(0, 5.0, 4.0),
            Err(MilpError::InvalidModel(_))
        ));
        assert!(matches!(
            p.set_var_bounds(0, f64::NAN, 1.0),
            Err(MilpError::InvalidModel(_))
        ));
        // Binary 固定为 0 或 1 合法；越出 [0,1] 拒绝。
        assert!(p.set_var_bounds(1, 0.0, 0.0).is_ok());
        assert!(p.set_var_bounds(1, 1.0, 1.0).is_ok());
        assert!(matches!(
            p.set_var_bounds(1, 0.0, 2.0),
            Err(MilpError::InvalidModel(_))
        ));
        assert!(matches!(
            p.set_var_bounds(1, -0.5, 1.0),
            Err(MilpError::InvalidModel(_))
        ));
        // 下标越界拒绝。
        assert!(matches!(
            p.set_var_bounds(9, 0.0, 1.0),
            Err(MilpError::InvalidModel(_))
        ));
    }

    /// 目标项归并排序、零系数丢弃、下标越界拒绝。
    #[test]
    fn set_objective_merges_and_bounds_check() {
        let mut p = small_problem();
        p.set_objective(vec![(2, 1.0), (0, 5.0), (0, -2.0), (1, 0.0)], 0.0)
            .unwrap();
        assert_eq!(p.obj.terms, vec![(0, 3.0), (2, 1.0)]);
        assert_eq!(p.obj_offset, 0.0);
        assert!(matches!(
            p.set_objective(vec![(7, 1.0)], 0.0),
            Err(MilpError::InvalidModel(_))
        ));
    }

    /// 约束项归并排序；全零系数行为空表达式约束照常入模型；越界拒绝。
    #[test]
    fn add_constraint_merges_and_allows_empty_row() {
        let mut p = small_problem();
        let i = p
            .add_constraint(
                Some("c1"),
                vec![(1, 2.0), (0, 1.0), (0, 1.0)],
                RelOp::Le,
                12.0,
            )
            .unwrap();
        assert_eq!(i, 0);
        assert_eq!(p.constraints[0].expr.terms, vec![(0, 2.0), (1, 2.0)]);
        assert_eq!(p.constraints[0].op, RelOp::Le);
        assert_eq!(p.constraints[0].rhs, 12.0);
        assert_eq!(p.num_constraints(), 1);
        p.add_constraint(Some("c2"), vec![(0, 0.0), (1, 0.0)], RelOp::Ge, -1.0)
            .unwrap();
        assert!(
            p.constraints[1].expr.terms.is_empty(),
            "整行系数皆 0 时空行照常进入模型"
        );
        assert!(matches!(
            p.add_constraint(None, vec![(99, 1.0)], RelOp::Le, 0.0),
            Err(MilpError::InvalidModel(_))
        ));
    }

    /// integer_vars 返回 Binary+Integer 的升序下标。
    #[test]
    fn integer_vars_ascending() {
        let p = small_problem();
        assert_eq!(p.integer_vars(), vec![0, 1]);
        let q = Problem::new("q", Sense::Minimize);
        assert!(q.integer_vars().is_empty());
    }

    /// 无 bounds 节的连续/整数变量（ub=+∞）必须通过校验。
    #[test]
    fn validate_accepts_default_infinite_ub() {
        let mut p = Problem::new("lp-prod-3", Sense::Maximize);
        p.add_var("x", VarType::Continuous).unwrap();
        p.add_var("k", VarType::Integer).unwrap();
        p.set_objective(vec![(0, 1.0), (1, 2.0)], 0.0).unwrap();
        p.add_constraint(None, vec![(0, 1.0), (1, 1.0)], RelOp::Le, 4.0)
            .unwrap();
        assert!(p.validate().is_ok(), "默认 ub=+∞ 不得报错");
    }

    /// validate 拒绝各类非法状态（直接改 pub 字段构造）。
    #[test]
    fn validate_rejects_invalid_states() {
        let base = || {
            let mut p = Problem::new("t", Sense::Minimize);
            p.add_var("x", VarType::Continuous).unwrap();
            p.set_objective(vec![(0, 1.0)], 0.0).unwrap();
            p
        };
        assert!(base().validate().is_ok(), "基线必须通过");

        let mut q = base();
        q.obj_offset = f64::NAN;
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "obj_offset NaN"
        );

        let mut q = base();
        q.obj.terms = vec![(0, f64::INFINITY)];
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "系数 inf"
        );

        let mut q = base();
        q.add_constraint(None, vec![(0, 1.0)], RelOp::Le, 1.0)
            .unwrap();
        q.constraints[0].rhs = f64::INFINITY;
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "右端 inf"
        );

        let mut q = base();
        q.vars[0].lb = f64::NEG_INFINITY;
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "lb=-∞"
        );

        let mut q = base();
        q.vars[0].ub = f64::NEG_INFINITY;
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "ub=-∞"
        );

        let mut q = base();
        q.vars[0].lb = 5.0;
        q.vars[0].ub = 1.0;
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "lb>ub"
        );

        let mut q = base();
        q.obj.terms = vec![(0, 0.0)];
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "零系数项"
        );

        let mut q = base();
        q.obj.terms = vec![(0, 1.0), (0, 2.0)];
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "重复下标"
        );

        let mut q = base();
        q.obj.terms = vec![(3, 1.0)];
        assert!(
            matches!(q.validate(), Err(MilpError::InvalidModel(_))),
            "下标越界"
        );
    }
}
