//! 给跨引擎基准提供可预热的单线程运行入口。读入业务模型和写出最终解均在计时区间外。

use planner::business::{optimize_timed, wire, Config, State};
use planner::{ConstructionHeuristic, Metaheuristic};
use std::io::{self, BufRead, Write};
use std::time::{Duration, Instant};

fn solution_json(state: &State) -> (String, String) {
    let values = state
        .values
        .iter()
        .map(|v| v.map_or_else(|| "null".to_owned(), |n| n.to_string()))
        .collect::<Vec<_>>()
        .join(",");
    let lists: Vec<_> = state.lists.iter().flatten().collect();
    (format!("[{values}]"), format!("{lists:?}"))
}

fn run(fields: &[&str]) -> Result<(), String> {
    let [file, seed, budget] = fields else {
        return Err("参数应为 模型文件 随机种子 毫秒预算".into());
    };
    let input = std::fs::read_to_string(file).map_err(|e| e.to_string())?;
    let (model, _) = wire::parse(&input)?;
    let cfg = Config {
        construction: ConstructionHeuristic::BestFit,
        algorithm: Metaheuristic::TabuSearch { tabu_size: 8 },
        seed: seed.parse().map_err(|_| "随机种子必须是非负整数")?,
        steps: usize::MAX,
        threads: 1,
        time_limit: Some(Duration::from_millis(
            budget.parse().map_err(|_| "预算必须是非负毫秒数")?,
        )),
    };
    let start = Instant::now();
    let (result, first) = optimize_timed(&model, &cfg)?;
    let elapsed = start.elapsed();
    let (first_ms, first_values, first_lists) = match first {
        Some(event) => {
            let (v, l) = solution_json(&event.state);
            (format!("{:.6}", event.at.as_secs_f64() * 1000.0), v, l)
        }
        None => ("null".to_owned(), "null".to_owned(), "null".to_owned()),
    };
    let (values, lists) = solution_json(&result.state);
    println!(
        r#"RESULT {{"hard":{},"soft":{},"complete":{},"first_ms":{first_ms},"elapsed_ms":{:.6},"values":{values},"lists":{lists},"first_values":{first_values},"first_lists":{first_lists}}}"#,
        result.score.hard,
        result.score.soft,
        result.complete,
        elapsed.as_secs_f64() * 1000.0,
    );
    io::stdout().flush().map_err(|e| e.to_string())
}

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !args.is_empty() {
        return run(&args.iter().map(String::as_str).collect::<Vec<_>>());
    }
    for line in io::stdin().lock().lines() {
        let line = line.map_err(|e| e.to_string())?;
        if !line.trim().is_empty() {
            run(&line.split_whitespace().collect::<Vec<_>>())?;
        }
    }
    Ok(())
}
