//! 业务规划模型。变量和值表均按注册顺序编号，同一实体可注册多个变量。
//! 模型经校验后只读，各搜索线程独立维护工作解和评分缓存。

pub mod score;
mod search;
mod sequence;
pub mod wire;

use crate::construct::ConstructionHeuristic;
use crate::localsearch::Metaheuristic;
pub use search::{optimize, optimize_timed, Config, FirstFeasible, Outcome, StopReason};

/// 每种列表变量把元素分给若干所有者，同时规定各自的访问或加工顺序。
#[derive(Debug, Clone)]
pub struct ListSpec {
    pub owners: usize,
    pub eligible: Vec<Vec<usize>>,
    pub durations: Vec<i64>,
    pub precedences: Vec<(usize, usize)>,
}

/// 节点只能引用在它之前声明的节点。比较的结果为 0 或 1。
#[derive(Debug, Clone)]
pub enum Expr {
    Const(i64),
    Var(usize),
    Add(usize, usize),
    Sub(usize, usize),
    Mul(usize, usize),
    Eq(usize, usize),
    Ne(usize, usize),
    Lt(usize, usize),
    Le(usize, usize),
    Gt(usize, usize),
    Ge(usize, usize),
    If(usize, usize, usize),
    Sum(Vec<usize>),
    Max(Vec<usize>),
    Lookup(Vec<i64>, usize),
    Load {
        group: usize,
        owner: usize,
        demands: Vec<i64>,
    },
    Distance {
        group: usize,
        owner: usize,
        matrix: Vec<i64>,
    },
    /// 同一列表中，after 是否紧跟在 before 后面。
    Adjacent {
        group: usize,
        owner: usize,
        before: usize,
        after: usize,
    },
    Finish {
        group: usize,
        item: usize,
    },
}

#[derive(Debug, Clone)]
pub struct Rule {
    pub name: String,
    pub hard: bool,
    pub weight: i64,
    pub expr: usize,
}

#[derive(Debug, Default, Clone)]
pub struct Model {
    pub domains: Vec<Vec<i64>>,
    pub lists: Vec<ListSpec>,
    pub nodes: Vec<Expr>,
    pub rules: Vec<Rule>,
    pub objective: Option<usize>,
}

impl Model {
    /// 返回节点编号，供后续约束或目标引用。
    pub fn expr(&mut self, expr: Expr) -> usize {
        self.nodes.push(expr);
        self.nodes.len() - 1
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct State {
    pub values: Vec<Option<i64>>,
    pub lists: Vec<Vec<Vec<usize>>>,
}

impl State {
    fn empty(model: &Model) -> Self {
        Self {
            values: vec![None; model.domains.len()],
            lists: model
                .lists
                .iter()
                .map(|g| vec![Vec::new(); g.owners])
                .collect(),
        }
    }

    fn complete(&self, model: &Model) -> bool {
        self.values.iter().all(Option::is_some)
            && self.lists.iter().zip(&model.lists).all(|(owners, spec)| {
                owners.iter().map(Vec::len).sum::<usize>() == spec.eligible.len()
            })
    }
}
