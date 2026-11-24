//! 与 Python 和 Java 的示例使用相同模型、配置与输出顺序。
#[path = "business_models/cases.rs"]
mod cases;
use planner::{
    business::{optimize, Config},
    ConstructionHeuristic,
};

fn joined<T: ToString>(values: impl IntoIterator<Item = T>) -> String {
    values
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn main() -> Result<(), String> {
    let config = Config {
        seed: 7,
        steps: 120,
        construction: ConstructionHeuristic::BestFit,
        ..Config::default()
    };
    for (name, model) in cases::models() {
        let r = optimize(&model, &config)?;
        let values = joined(
            r.state
                .values
                .iter()
                .map(|v| v.map_or_else(|| "_".into(), |v| v.to_string())),
        );
        let lists = r
            .state
            .lists
            .iter()
            .flatten()
            .map(joined)
            .collect::<Vec<_>>()
            .join("/");
        let status = r.reason.status();
        println!(
            "{name}\t{}\t{}\t{}\t{status}\t{:?}\t{}\t{values}\t{lists}\t{}",
            r.score.hard,
            r.score.soft,
            r.steps,
            r.reason,
            u8::from(r.complete),
            joined(r.score.penalties)
        );
    }
    Ok(())
}
