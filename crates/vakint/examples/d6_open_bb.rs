//! D6 (helper VH, lane d6-open-indices): the background-field `B -> B` DOT
//! graphs through the OPEN-INDEX route (`refpath::run_graph_open`: ports left
//! open, tensor reduction with the external momentum declared, the open
//! `Π^{μν}` returned) against the PROJECTED scalar route (`refpath::run_graph`,
//! metric and transverse + longitudinal), per graph: exact equality of the open
//! result closed afterwards with each projection, and the stage seconds of
//! both routes (medians over `--reps` paired repetitions on fresh engines, the
//! routes alternating).
//!
//! ```text
//! d6_open_bb <prepared dir> --out <dir> [--graphs 1,2,7-9 | --graph-list <file>] [--reps <n>] [--p-library <dir>]
//! ```
//! Writes `<out>/results.jsonl` (one object per graph).
use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::PathBuf;
use std::time::Instant;

use serde_json::json;
use symbolica::atom::Atom;
use vakint::refpath::{self, Projection, RunOptions, StageSeconds};
use vakint::{Vakint, vakint_parse, vakint_symbol};

#[path = "../tests/support/open_indices.rs"]
mod open_indices;

fn parse_graphs(text: &str) -> Vec<usize> {
    let mut out = Vec::new();
    for part in text.split(',').filter(|part| !part.is_empty()) {
        match part.split_once('-') {
            Some((a, b)) => out.extend(a.parse::<usize>().unwrap()..=b.parse::<usize>().unwrap()),
            None => out.push(part.parse().unwrap()),
        }
    }
    out
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n == 0 {
        0.0
    } else if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}

fn stages(reps: &[StageSeconds]) -> serde_json::Value {
    let m = |get: fn(&StageSeconds) -> f64| median(reps.iter().map(get).collect());
    json!({
        "parse": m(|s| s.parse), "colour": m(|s| s.colour), "contract": m(|s| s.contract),
        "ingress": m(|s| s.ingress), "lower": m(|s| s.lower), "reduce": m(|s| s.reduce),
        "finalise": m(|s| s.finalise), "split": m(|s| s.split), "lowering": m(|s| s.lowering),
    })
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let (mut directory, mut out, mut graphs, mut reps) = (None, None, None, 1usize);
    let mut library = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").map(PathBuf::from);
    let mut index = 0;
    while index < args.len() {
        let value = || args.get(index + 1).cloned().expect("missing value");
        match args[index].as_str() {
            "--out" => {
                out = Some(PathBuf::from(value()));
                index += 1;
            }
            "--graphs" => {
                graphs = Some(parse_graphs(&value()));
                index += 1;
            }
            "--graph-list" => {
                let text = std::fs::read_to_string(value())?;
                graphs = Some(
                    text.lines()
                        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
                        .map(|line| line.trim().parse().unwrap())
                        .collect(),
                );
                index += 1;
            }
            "--reps" => {
                reps = value().parse()?;
                index += 1;
            }
            "--p-library" => {
                library = Some(PathBuf::from(value()));
                index += 1;
            }
            other if directory.is_none() && !other.starts_with("--") => {
                directory = Some(PathBuf::from(other))
            }
            other => eyre::bail!("unknown argument {other}"),
        }
        index += 1;
    }
    let directory = directory.ok_or_else(|| eyre::eyre!("no prepared directory"))?;
    let out = out.ok_or_else(|| eyre::eyre!("no --out"))?;
    std::fs::create_dir_all(&out)?;
    let model = refpath::load_model(&directory)?;
    let mut files = std::fs::read_dir(&directory)?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .is_some_and(|name| name.starts_with("diagram-") && name.ends_with(".dot"))
        })
        .collect::<Vec<_>>();
    files.sort();
    let number = |path: &PathBuf| -> usize {
        path.file_stem().unwrap().to_str().unwrap().trim_start_matches("diagram-").parse().unwrap()
    };
    let selected = match &graphs {
        Some(list) => {
            let by_number = files.iter().map(|p| (number(p), p.clone())).collect::<BTreeMap<_, _>>();
            list.iter()
                .map(|n| by_number.get(n).cloned().ok_or_else(|| eyre::eyre!("no diagram {n}")))
                .collect::<eyre::Result<Vec<_>>>()?
        }
        None => files.clone(),
    };
    let mut template = Vakint::new()?;
    if let Some(library) = &library {
        template.load_p_library_directory(library, false)?;
    }
    let mut results = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(out.join("results.jsonl"))?;
    eprintln!("[d6_open_bb] {} graphs of {}, reps {reps}, library {library:?}", selected.len(), directory.display());
    let q2 = vakint_parse!("dot(p(1),p(1))").unwrap();
    let g = vakint_parse!("g(1,2)").unwrap();
    let pp = vakint_parse!("p(1,1)*p(1,2)").unwrap();
    let mut mismatches = 0usize;
    for path in &selected {
        let graph = number(path);
        let dot = std::fs::read_to_string(path)?;
        let options = |projection| RunOptions { projection, extra_terms: 1, count_terms: false, contract_only: false };
        let (mut metric_s, mut both_s, mut open_s) = (Vec::new(), Vec::new(), Vec::new());
        let (mut canon_probe, mut tensor_probe) = (Vec::new(), Vec::new());
        let mut references = BTreeMap::<&'static str, Atom>::new();
        let mut open_value = None;
        let mut loops = 0;
        let mut zero = None;
        let mut counts = json!({});
        for rep in 0..reps {
            let mut run_projected = |projection| -> eyre::Result<()> {
                let engine = refpath::fresh_engine(&template)?;
                let run = refpath::run_graph(&engine, &model, &dot, &options(projection))?;
                loops = run.loops;
                for result in &run.results {
                    references.entry(result.projection).or_insert_with(|| result.value.clone());
                }
                match projection {
                    Projection::Metric => metric_s.push(run.seconds.clone()),
                    _ => both_s.push(run.seconds.clone()),
                }
                Ok(())
            };
            let mut run_open = || -> eyre::Result<()> {
                let engine = refpath::fresh_engine(&template)?;
                let run = refpath::run_graph_open(&engine, &model, &dot, 1)?;
                zero = run.zero_reason.clone();
                canon_probe.push(run.canonicalize_probe);
                tensor_probe.push(run.tensor_probe);
                open_s.push(run.seconds.clone());
                if let Some(observation) = &run.observation {
                    counts = json!({"reductions": observation.reductions, "channels": observation.channels,
                        "pairs": observation.pairs, "folded": observation.folded, "direct_lowered": observation.direct_lowered});
                }
                if open_value.is_none() {
                    open_value = Some(run.value);
                }
                Ok(())
            };
            // Paired repetitions, the routes alternating.
            if rep % 2 == 0 {
                run_projected(Projection::Metric)?;
                run_projected(Projection::Both)?;
                run_open()?;
            } else {
                run_open()?;
                run_projected(Projection::Metric)?;
                run_projected(Projection::Both)?;
            }
        }
        // Exactness: the open value closed with each projection.
        let open = open_value.unwrap_or(Atom::Zero);
        let settings = refpath::vakint_settings(loops.max(1) as i64 + 1);
        let epsilon = vakint_symbol!(settings.epsilon_symbol.as_str());
        let d = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
        let i = Atom::i();
        let mut check = BTreeMap::new();
        let started = Instant::now();
        for (name, closing, factor) in [
            ("metric", g.clone(), (&d - Atom::num(1)) * &i * &q2),
            ("transverse", &g - &pp / &q2, (&d - Atom::num(1)) * &i * &q2),
            ("longitudinal", pp.clone(), &i * &q2 * &q2),
        ] {
            let reference = references.get(name).cloned().unwrap_or(Atom::Zero);
            let closed = open_indices::close_indices(&(&open * &closing), epsilon);
            let last = open_indices::last_order(&reference, &closed, epsilon);
            let ok = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                open_indices::assert_series_through(&closed, &(&factor * &reference), last, name, epsilon)
            }))
            .is_ok();
            if !ok {
                mismatches += 1;
            }
            check.insert(name, ok);
        }
        let check_seconds = started.elapsed().as_secs_f64();
        let record = json!({
            "graph": graph, "file": path.display().to_string(), "loops": loops, "zero": zero, "reps": reps,
            "metric": stages(&metric_s), "both": stages(&both_s), "open": stages(&open_s),
            "open_canonicalize_probe": median(canon_probe), "open_tensor_probe": median(tensor_probe),
            "open_counts": counts, "check": check, "check_seconds": check_seconds,
        });
        writeln!(results, "{record}")?;
        eprintln!("[d6_open_bb] diagram {graph}: {check:?}");
    }
    eprintln!("[d6_open_bb] done: {} graphs, {mismatches} projection mismatches", selected.len());
    if mismatches > 0 {
        eyre::bail!("{mismatches} projection mismatches");
    }
    Ok(())
}
