//! D6 through the hand-over (helper VH, lane idenso-audit): the background-field
//! `B -> B` DOT graphs with the two external Lorentz ports left OPEN, through
//! vakint's open-index route (`refpath::run_graph_open`: the open numerator
//! through the FeynKit ingress, `Vakint::evaluate`) and through R3-open
//! (`refpath::direct::lower_graph_open`: idenso's contraction with the ports
//! open, ONE denominator match, the open rename, the tensor reduction of the
//! factorised numerator with `p(1)` external, the direct reader), per graph:
//! the open values exactly equal, the folded inputs key by key (views without
//! zero-sector keys), and the stage seconds of both routes (medians over
//! `--reps` paired repetitions on fresh engines, the routes alternating).
//!
//! ```text
//! r3_open_bb <prepared dir> --out <dir> [--graphs 1,2,7-9 | --graph-list <file>] [--reps <n>]
//!     [--p-library <dir>] [--merge-colour] [--no-compact] [--no-dots] [--no-prune]
//! ```
//! Writes `<out>/results.jsonl` (one object per graph).
use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use serde_json::json;
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;
use vakint::refpath::{self, direct};
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

fn median(mut values: Vec<f64>) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        0.5 * (values[n / 2 - 1] + values[n / 2])
    }
}

/// One R3-open evaluation of a graph: value, gate view and stage seconds.
struct R3Open {
    zero: Option<String>,
    value: Atom,
    view: Option<p_gate::GateView>,
    seconds: BTreeMap<&'static str, f64>,
}

fn r3_open(
    engine: &Vakint,
    model: &Arc<feynkit_model::Model>,
    dot: &str,
    merge: bool,
    options: direct::DirectOptions,
) -> eyre::Result<R3Open> {
    let s = &*refpath::SYMBOLS;
    let mut seconds = BTreeMap::new();
    let started = Instant::now();
    let diagram = feynkit_graph::FeynmanDiagram::from_dot(Arc::clone(model), dot)?;
    seconds.insert("parse", started.elapsed().as_secs_f64());
    let zero = |reason: &str| R3Open {
        zero: Some(reason.to_string()),
        value: Atom::Zero,
        view: None,
        seconds: BTreeMap::new(),
    };
    if refpath::scaleless_contact(&diagram) {
        return Ok(zero(
            "no internal edge carries the external momentum: scaleless",
        ));
    }
    let started = Instant::now();
    let (terms, _) = direct::colour_terms(&diagram, merge)?;
    seconds.insert("colour", started.elapsed().as_secs_f64());
    if terms.is_empty() {
        return Ok(zero("every colour branch vanishes"));
    }
    let dimension = Atom::var(s.d);
    let kinematics = feynkit_kinematics::Kinematics::in_dimension(&dimension)?;
    let family = diagram.propagator_family(&kinematics)?;
    let [momentum] = family.external_momenta() else {
        eyre::bail!("a two-point graph needs one independent external momentum");
    };
    let momentum = momentum.clone();
    let weight =
        feynkit_graph::expressions::evaluate_overall_factor(diagram.overall_factor().as_view())
            * Atom::var(s.nf).pow(Atom::num(refpath::closed_quark_loops(&diagram)? as i64));
    // run_graph_open's settings: the loop count + 1 extra term
    let settings = refpath::vakint_settings(diagram.loop_count() as i64 + 1);
    let mut stages = direct::DirectSeconds::default();
    let outcome = direct::lower_graph_open(
        engine,
        &settings,
        &diagram,
        &terms,
        &momentum,
        &weight,
        options,
        &mut stages,
    )?;
    for (name, value) in [
        ("tree", stages.routing),
        ("contract", stages.contract),
        ("match", stages.matching),
        ("rename", stages.handover),
        ("tensor", stages.tensor),
        ("lower", stages.lower),
    ] {
        seconds.insert(name, value);
    }
    let projection = match outcome {
        direct::DirectOutcome::Declined(reason) => eyre::bail!("R3-open declined: {reason}"),
        direct::DirectOutcome::Zero(reason) => {
            let mut z = zero(reason);
            z.seconds = seconds;
            return Ok(z);
        }
        direct::DirectOutcome::Lowered(mut projections) => projections.remove(0),
    };
    let (view, _) = direct::gate_view(engine, &settings, &projection)?;
    let (value, timed) = direct::reduce_projection_timed(engine, &settings, &projection)?;
    seconds.insert("reduce", timed.reduce);
    seconds.insert("finalise", timed.finalise);
    Ok(R3Open {
        zero: None,
        value,
        view: Some(view),
        seconds,
    })
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let (mut directory, mut out, mut graphs, mut reps) = (None, None, None, 1usize);
    let mut library = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").map(PathBuf::from);
    let mut merge = false;
    let mut options = direct::DirectOptions::default();
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
            "--merge-colour" => merge = true,
            "--no-compact" => options.compact = false,
            "--no-dots" => options.dots = false,
            "--no-prune" => options.prune_early = false,
            other if directory.is_none() && !other.starts_with("--") => {
                directory = Some(PathBuf::from(other))
            }
            other => eyre::bail!("unknown argument {other}"),
        }
        index += 1;
    }
    let directory = directory.ok_or_else(|| eyre::eyre!("no prepared directory"))?;
    let out = out.ok_or_else(|| eyre::eyre!("no --out"))?;
    let graphs = graphs.ok_or_else(|| eyre::eyre!("no --graphs / --graph-list"))?;
    std::fs::create_dir_all(&out)?;
    let model = refpath::load_model(&directory)?;
    let mut template = Vakint::new()?;
    if let Some(library) = &library {
        template.load_p_library_directory(library, false)?;
    }
    let mut results = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(out.join("results.jsonl"))?;
    eprintln!(
        "[r3_open_bb] {} graphs of {}, reps {reps}, merge {merge}, {options:?}",
        graphs.len(),
        directory.display()
    );
    let (mut value_mismatch, mut key_mismatch) = (0usize, 0usize);
    for graph in &graphs {
        let path = directory.join(format!("diagram-{graph:06}.dot"));
        let dot = std::fs::read_to_string(&path)?;
        let mut open_s: Vec<BTreeMap<&'static str, f64>> = Vec::new();
        let mut r3_s: Vec<BTreeMap<&'static str, f64>> = Vec::new();
        let (mut open_value, mut open_view, mut open_zero) = (None, None, None);
        let (mut r3_value, mut r3_view, mut r3_zero) = (None, None, None);
        let mut open_counts = json!({});
        for rep in 0..reps {
            let mut run_open = || -> eyre::Result<()> {
                let engine = refpath::fresh_engine(&template)?;
                let run = refpath::run_graph_open(&engine, &model, &dot, 1)?;
                open_zero = run.zero_reason.clone();
                let s = &run.seconds;
                open_s.push(BTreeMap::from([
                    ("parse", s.parse),
                    ("colour", s.colour),
                    ("contract", s.contract),
                    ("ingress", s.ingress),
                    ("lower", s.lower),
                    ("reduce", s.reduce),
                    ("finalise", s.finalise),
                    ("probe_match", run.canonicalize_probe),
                    ("probe_tensor", run.tensor_probe),
                ]));
                if let Some(observation) = &run.observation {
                    if open_view.is_none() {
                        open_view = Some(observation.view_without_zero_sectors().0);
                    }
                    open_counts = json!({"reductions": observation.reductions,
                        "channels": observation.channels, "folded": observation.folded});
                }
                if open_value.is_none() {
                    open_value = Some(run.value);
                }
                Ok(())
            };
            let mut run_r3 = || -> eyre::Result<()> {
                let engine = refpath::fresh_engine(&template)?;
                let run = r3_open(&engine, &model, &dot, merge, options)?;
                r3_zero = run.zero.clone();
                r3_s.push(run.seconds);
                if r3_value.is_none() {
                    r3_value = Some(run.value);
                    r3_view = run.view;
                }
                Ok(())
            };
            if rep % 2 == 0 {
                run_open()?;
                run_r3()?;
            } else {
                run_r3()?;
                run_open()?;
            }
        }
        let open_value = open_value.unwrap_or(Atom::Zero);
        let r3_value = r3_value.unwrap_or(Atom::Zero);
        let values_equal = direct::values_equal(&open_value, &r3_value);
        let keys = match (&open_view, &r3_view) {
            (Some(a), Some(b)) => {
                let difference = p_gate::compare_gate_views(a, b);
                json!({"identical": difference.is_empty(), "open_keys": a.len(), "r3_keys": b.len(),
                    "difference": format!("{difference:?}").chars().take(400).collect::<String>()})
            }
            // no view or an empty one: a zero on that side (the open route reduces
            // nothing; R3-open stops at "scaleless after the match")
            (a, b) => json!({
                "identical": a.as_ref().is_none_or(|v| v.is_empty())
                    && b.as_ref().is_none_or(|v| v.is_empty()),
                "open_keys": a.as_ref().map(|v| v.len()),
                "r3_keys": b.as_ref().map(|v| v.len()),
            }),
        };
        if !values_equal {
            value_mismatch += 1;
        }
        if keys["identical"] != true {
            key_mismatch += 1;
        }
        let stages = |runs: &[BTreeMap<&'static str, f64>]| {
            let names = runs
                .iter()
                .flat_map(|r| r.keys().copied())
                .collect::<std::collections::BTreeSet<_>>();
            names
                .into_iter()
                .map(|n| {
                    (
                        n,
                        median(
                            runs.iter()
                                .map(|r| r.get(n).copied().unwrap_or(0.0))
                                .collect(),
                        ),
                    )
                })
                .collect::<BTreeMap<_, _>>()
        };
        let record = json!({
            "graph": graph, "reps": reps, "open_zero": open_zero, "r3_zero": r3_zero,
            "values_equal": values_equal, "keys": keys, "open_counts": open_counts,
            "open": stages(&open_s), "r3": stages(&r3_s),
        });
        writeln!(results, "{record}")?;
        eprintln!(
            "[r3_open_bb] diagram {graph}: values_equal {values_equal}, keys {}",
            keys["identical"]
        );
    }
    eprintln!(
        "[r3_open_bb] done: {} graphs, {value_mismatch} value mismatches, {key_mismatch} key mismatches",
        graphs.len()
    );
    if value_mismatch > 0 {
        eyre::bail!("{value_mismatch} value mismatches");
    }
    Ok(())
}
