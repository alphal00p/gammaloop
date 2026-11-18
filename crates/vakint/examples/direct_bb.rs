//! The reference path R2 on background-field `B -> B` self-energy DOT graphs
//! (helper VH, numerator assignment, lane `refpath`; `vakint::refpath`).
//!
//! ```text
//! direct_bb <prepared dir with model.json + diagram-NNNNNN.dot> --out <dir>
//!     [--graphs 1,2,7-9 | --graph-list <file with one index per line>]
//!     [--projection metric|transverse|longitudinal|both]   (default metric)
//!     [--reps <n>]            repetitions per graph, each on a fresh engine (default 1)
//!     [--extra-terms <k>]     Laurent terms = loops + k (default 1: poles and finite part)
//!     [--p-library <dir>]     RustRed p package (default $VAKINT_TEST_P_LIBRARY_1_4, else built-in 1L-3L)
//!     [--no-count]            skip the term count of the expanded numerator
//!     [--contract-only]       stop after the contraction
//!     [--warmup]              one untimed run of the first graph before the timed ones
//!     [--no-gate]             do not write gate dumps
//! ```
//!
//! Writes `<out>/results.jsonl` (one JSON object per graph and repetition
//! summary), `<out>/gate/<graph>.<projection>.gate` (the gate view without
//! zero-sector keys, `vakint.p-gate-view.v1`) and `<out>/values/<graph>.txt`
//! (the Laurent expansion per projection, canonical strings).
use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::PathBuf;
use std::time::Instant;

use serde_json::json;
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;
use vakint::refpath::{self, RunOptions};
use vakint::rustred_evaluation::p_gate;

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

fn median(values: &mut [f64]) -> f64 {
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

fn cpu_seconds() -> f64 {
    // utime + stime of this process (clock ticks of 100 Hz).
    let stat = std::fs::read_to_string("/proc/self/stat").unwrap_or_default();
    let fields = stat
        .rsplit_once(')')
        .map(|(_, rest)| rest.split_whitespace().collect::<Vec<_>>())
        .unwrap_or_default();
    let tick = |i: usize| fields.get(i).and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    (tick(11) + tick(12)) / 100.0
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let mut directory = None;
    let mut out = None;
    let mut graphs = None;
    let mut options = RunOptions::default();
    let mut reps = 1usize;
    let mut library = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").map(PathBuf::from);
    let mut gate = true;
    let mut warmup = false;
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
            "--projection" => {
                options.projection = value().parse().map_err(|e: String| eyre::eyre!(e))?;
                index += 1;
            }
            "--reps" => {
                reps = value().parse()?;
                index += 1;
            }
            "--extra-terms" => {
                options.extra_terms = value().parse()?;
                index += 1;
            }
            "--p-library" => {
                library = Some(PathBuf::from(value()));
                index += 1;
            }
            "--no-count" => options.count_terms = false,
            "--contract-only" => options.contract_only = true,
            "--no-gate" => gate = false,
            "--warmup" => warmup = true,
            other if directory.is_none() && !other.starts_with("--") => {
                directory = Some(PathBuf::from(other))
            }
            other => eyre::bail!("unknown argument {other}"),
        }
        index += 1;
    }
    let directory = directory.ok_or_else(|| eyre::eyre!("no prepared directory"))?;
    let out = out.ok_or_else(|| eyre::eyre!("no --out"))?;
    std::fs::create_dir_all(out.join("gate"))?;
    std::fs::create_dir_all(out.join("values"))?;
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
        let name = path.file_stem().unwrap().to_str().unwrap();
        name.trim_start_matches("diagram-").parse().unwrap()
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
    let mut results = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(out.join("results.jsonl"))?;
    let host = std::fs::read_to_string("/proc/sys/kernel/hostname").unwrap_or_default();
    eprintln!(
        "[direct_bb] {} graphs of {} projection {} reps {reps} library {:?}",
        selected.len(),
        directory.display(),
        options.projection.name(),
        library
    );
    if warmup && let Some(path) = selected.first() {
        // One untimed run of the first graph: process-wide lazy statics and
        // symbol registrations are not charged to the first timed graph.
        let mut engine = Vakint::new()?;
        if let Some(library) = &library {
            engine.load_p_library_directory(library, false)?;
        }
        let engine = refpath::fresh_engine(&engine)?;
        let started = Instant::now();
        let _ = refpath::run_graph(&engine, &model, &std::fs::read_to_string(path)?, &options);
        eprintln!("[direct_bb] warm-up on {} in {:.3} s", path.display(), started.elapsed().as_secs_f64());
    }
    let setup = Instant::now();
    let mut template = Vakint::new()?;
    if let Some(library) = &library {
        template.load_p_library_directory(library, false)?;
    }
    eprintln!("[direct_bb] engine and p library loaded in {:.3} s", setup.elapsed().as_secs_f64());
    for path in &selected {
        let graph_started = Instant::now();
        let dot = std::fs::read_to_string(path)?;
        let graph = number(path);
        let mut rep_seconds = Vec::new();
        let mut first = None;
        let mut rss = Vec::new();
        let mut cpu = Vec::new();
        let mut wall = Vec::new();
        let mut setup_seconds = Vec::new();
        for rep in 0..reps {
            let setup = Instant::now();
            let engine = refpath::fresh_engine(&template)?;
            setup_seconds.push(setup.elapsed().as_secs_f64());
            refpath::start_peak_rss();
            let cpu_start = cpu_seconds();
            let started = Instant::now();
            let run = refpath::run_graph(&engine, &model, &dot, &options);
            wall.push(started.elapsed().as_secs_f64());
            cpu.push(cpu_seconds() - cpu_start);
            rss.push(refpath::peak_rss_bytes().unwrap_or(0));
            let run = match run {
                Ok(run) => run,
                Err(error) => {
                    let record = json!({"graph": graph, "file": path.display().to_string(), "rep": rep, "error": format!("{error:?}")});
                    writeln!(results, "{record}")?;
                    eprintln!("[direct_bb] diagram {graph}: ERROR {error:?}");
                    break;
                }
            };
            rep_seconds.push(run.seconds.clone());
            if first.is_none() {
                first = Some(run);
            } else if let Some(reference) = &first {
                for (a, b) in reference.results.iter().zip(&run.results) {
                    if !(&a.value - &b.value).expand().is_zero() {
                        eprintln!("[direct_bb] diagram {graph}: repetition {rep} differs");
                    }
                }
            }
        }
        let Some(run) = first else { continue };
        let stage = |get: fn(&refpath::StageSeconds) -> f64| {
            let mut values = rep_seconds.iter().map(get).collect::<Vec<_>>();
            median(&mut values)
        };
        let seconds = json!({
            "parse": stage(|s| s.parse), "colour": stage(|s| s.colour), "contract": stage(|s| s.contract),
            "lower": stage(|s| s.lower), "reduce": stage(|s| s.reduce), "finalise": stage(|s| s.finalise),
            "ingress": stage(|s| s.ingress), "split": stage(|s| s.split), "lowering": stage(|s| s.lowering),
            "prepare": stage(|s| s.prepare), "library": stage(|s| s.library), "render": stage(|s| s.render),
            "finalize": stage(|s| s.finalize), "count": stage(|s| s.count),
        });
        let all_reps = rep_seconds
            .iter()
            .map(|s| json!([s.parse, s.colour, s.contract, s.lower, s.reduce, s.finalise]))
            .collect::<Vec<_>>();
        let mut projections = Vec::new();
        let mut values_text = String::new();
        let post = Instant::now();
        let mut view_seconds = 0.0;
        for result in &run.results {
            let started = Instant::now();
            let (view, count) = refpath::gate_view(result);
            view_seconds += started.elapsed().as_secs_f64();
            let all_keys = count.total;
            if gate {
                let header = BTreeMap::from([
                    ("graph".to_string(), run.name.clone()),
                    ("file".to_string(), path.display().to_string()),
                    ("projection".to_string(), result.projection.to_string()),
                    ("loops".to_string(), run.loops.to_string()),
                    ("zero_sectors_dropped".to_string(), count.zero.to_string()),
                    ("keys_before_drop".to_string(), count.total.to_string()),
                    ("census_unknown_kept".to_string(), count.unknown.to_string()),
                    ("producer".to_string(), "vakint refpath R2 (helper VH lane refpath)".to_string()),
                ]);
                std::fs::write(
                    out.join("gate").join(format!("{graph:06}.{}.gate", result.projection)),
                    p_gate::write_gate_view(&view, &header),
                )?;
            }
            values_text.push_str(&format!(
                "{}\t{}\n",
                result.projection,
                result.value.to_canonical_string()
            ));
            projections.push(json!({
                "projection": result.projection,
                "value": result.value.to_canonical_string(),
                "keys": all_keys,
                "keys_zero_sector": count.zero,
                "keys_census_unknown": count.unknown,
                "keys_kept": view.len(),
                "reductions": result.observation.reductions,
                "pairs": result.observation.pairs,
                "folded_channels": result.observation.folded,
                "split_channels": result.observation.channels,
            }));
        }
        for (name, numerator) in &run.numerators {
            values_text.push_str(&format!("numerator.{name}.bytes\t{}\n", numerator.to_canonical_string().len()));
        }
        std::fs::write(out.join("values").join(format!("{graph:06}.txt")), values_text)?;
        let post_seconds = post.elapsed().as_secs_f64();
        let record = json!({
            "post_seconds": post_seconds,
            "graph_total_seconds": graph_started.elapsed().as_secs_f64(),
            "gate_view_seconds": view_seconds,
            "graph": graph,
            "file": path.display().to_string(),
            "name": run.name,
            "loops": run.loops,
            "quark_loops": run.quark_loops,
            "overall_factor": run.overall_factor.to_canonical_string(),
            "branches_total": run.branches_total,
            "branches": run.branches,
            "zero_reason": run.zero_reason,
            "contracted_terms": run.contracted_terms,
            "reps": reps,
            "seconds_median": seconds,
            "seconds_reps": all_reps,
            "wall_reps": wall,
            "engine_setup_reps": setup_seconds,
            "cpu_reps": cpu,
            "peak_rss_bytes": rss.iter().max(),
            "projections": projections,
            "host": host.trim(),
        });
        writeln!(results, "{record}")?;
        results.flush()?;
        eprintln!(
            "[direct_bb] diagram {graph} L{} q{} branches {}/{} parse {:.3} colour {:.3} contract {:.3} lower {:.3} reduce {:.3} finalise {:.3} terms {} rss {:.2} GB {}",
            run.loops,
            run.quark_loops,
            run.branches,
            run.branches_total,
            stage(|s| s.parse),
            stage(|s| s.colour),
            stage(|s| s.contract),
            stage(|s| s.lower),
            stage(|s| s.reduce),
            stage(|s| s.finalise),
            run.contracted_terms,
            rss.iter().max().copied().unwrap_or(0) as f64 / (1u64 << 30) as f64,
            run.zero_reason.as_deref().unwrap_or("")
        );
        let _ = Atom::Zero;
    }
    Ok(())
}
