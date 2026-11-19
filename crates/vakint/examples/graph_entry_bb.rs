//! The GRAPH ENTRY through the hand-over (helper VH, lane idenso-audit): the
//! background-field `B -> B` DOT graphs through
//! `p_integral::factorized::evaluate_graph` (one hand-over call per graph on
//! the input type's colour branches) and through the ADAPTER route it took
//! before (`evaluate_graph_adapter`: per colour branch, idenso then
//! `Vakint::evaluate`), per graph and projector choice: the values exactly
//! equal in the caller's notation, the route taken, the seconds of both
//! (medians over `--reps` paired repetitions on fresh engines, the routes
//! alternating). With `--ref-dir`, the metric projection is also gated
//! against lane refpath's STORED dumps (`gate/<stem>.metric.gate`): the
//! entry's own contraction (couplings as given, unnormalised projector) and
//! lowering, fed refpath's branch coefficients times R2's weight and
//! normalisation, must hand RustRed the same folded inputs key by key.
//!
//! ```text
//! graph_entry_bb <prepared dir> --out <dir> [--graphs 1,2,7-9 | --graph-list <file>] [--reps <n>]
//!     [--p-library <dir>] [--projectors metric,transverse,none] [--merge-colour] [--ref-dir <dir>]
//! ```
//! Writes `<out>/results.jsonl` (one object per graph).
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use serde_json::json;
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;
use vakint::p_integral::FactorizedNumerator;
use vakint::p_integral::factorized::{
    FromFeynkitOptions, GraphRoute, ProjectorChoice, evaluate_graph, evaluate_graph_adapter,
};
use vakint::refpath::{self, ColourBranch, direct};
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

/// The metric projection through the entry's hand-over with refpath's
/// coefficients, gated against the stored dump of the graph.
fn dump_gate(
    engine: &Vakint,
    settings: &vakint::VakintSettings,
    diagram: &feynkit_graph::FeynmanDiagram,
    model: &feynkit_model::Model,
    reference: &Path,
    stem: &str,
) -> eyre::Result<serde_json::Value> {
    let s = &*refpath::SYMBOLS;
    let stored =
        std::fs::read_to_string(reference.join("gate").join(format!("{stem}.metric.gate"))).ok();
    let options = FromFeynkitOptions {
        dimension: s.d,
        projector: ProjectorChoice::None,
        ..FromFeynkitOptions::default()
    };
    let values = FactorizedNumerator::from_feynkit(diagram, model, &options)?;
    if values.is_empty() || refpath::scaleless_contact(diagram) {
        return Ok(
            json!({"status": "zero before the hand-over", "stored": stored.is_some(),
            "pass": stored.is_none()}),
        );
    }
    let weight = refpath::factorized_weight(&values)?;
    let d = Atom::var(s.d);
    let kinematics = feynkit_kinematics::Kinematics::in_dimension(&d)?;
    let family = diagram.propagator_family(&kinematics)?;
    let [momentum] = family.external_momenta() else {
        eyre::bail!("a two-point graph needs one independent external momentum");
    };
    let momentum = momentum.clone();
    let q_squared = feynkit_kinematics::Kinematics::in_dimension(&d)?
        .with_momenta([momentum.clone()])?
        .scalar_product(&momentum, &momentum)?;
    // R2's metric projection: g . Pi / ((D - 1) i q^2) times the weight
    let normalisation = weight / ((&d - Atom::num(1)) * Atom::i() * &q_squared);
    let couplings = |value: &Atom| {
        model
            .expand_couplings(value)
            .replace(Atom::var(symbolica::symbol!("UFO::gs")))
            .with(Atom::num(1))
    };
    let terms = values
        .iter()
        .map(|value| ColourBranch {
            coefficient: couplings(&(refpath::factorized_colour(value) * &value.prefactor.scalars))
                * &normalisation,
            residuals: vec![value.factor_product()],
        })
        .collect::<Vec<_>>();
    let mut seconds = direct::DirectSeconds::default();
    let outcome = direct::lower_entry(
        engine,
        settings,
        diagram,
        &terms,
        &momentum,
        direct::EntryProjector::Metric,
        direct::DirectOptions::default(),
        &mut seconds,
    )?;
    Ok(match outcome {
        direct::DirectOutcome::Declined(reason) => {
            json!({"status": format!("declined: {reason}"), "pass": false})
        }
        direct::DirectOutcome::Zero(reason) => {
            let empty = match &stored {
                None => true,
                Some(text) => p_gate::read_gate_view(text)
                    .map_err(|e| eyre::eyre!(e))?
                    .0
                    .is_empty(),
            };
            json!({"status": format!("zero: {reason}"), "stored": stored.is_some(), "pass": empty})
        }
        direct::DirectOutcome::Lowered(projections) => {
            let (view, _) = direct::gate_view(engine, settings, &projections[0])?;
            match &stored {
                None => json!({"status": "lowered", "stored": false, "keys": view.len(),
                    "pass": view.is_empty()}),
                Some(text) => {
                    let (stored_view, _) =
                        p_gate::read_gate_view(text).map_err(|e| eyre::eyre!(e))?;
                    let difference = p_gate::compare_gate_views(&stored_view, &view);
                    json!({"status": "lowered", "stored": true, "keys": view.len(),
                        "identical": difference.is_empty(), "pass": difference.is_empty(),
                        "difference": format!("{difference:?}").chars().take(300).collect::<String>()})
                }
            }
        }
    })
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let (mut directory, mut out, mut graphs, mut reps) = (None, None, None, 1usize);
    let mut library = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").map(PathBuf::from);
    let mut merge_colour = false;
    let mut reference = None;
    let mut projectors = vec![
        ProjectorChoice::Metric,
        ProjectorChoice::Transverse,
        ProjectorChoice::None,
    ];
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
            "--ref-dir" => {
                reference = Some(PathBuf::from(value()));
                index += 1;
            }
            "--projectors" => {
                projectors = value()
                    .split(',')
                    .map(|name| match name {
                        "metric" => Ok(ProjectorChoice::Metric),
                        "transverse" => Ok(ProjectorChoice::Transverse),
                        "none" | "open" => Ok(ProjectorChoice::None),
                        other => Err(eyre::eyre!("unknown projector {other}")),
                    })
                    .collect::<eyre::Result<_>>()?;
                index += 1;
            }
            "--merge-colour" => merge_colour = true,
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
        "[graph_entry_bb] {} graphs of {}, reps {reps}, projectors {projectors:?}, merge {merge_colour}",
        graphs.len(),
        directory.display()
    );
    let (mut value_mismatch, mut not_handed, mut dump_fail) = (0usize, 0usize, 0usize);
    for graph in &graphs {
        let stem = format!("{graph:06}");
        let path = directory.join(format!("diagram-{stem}.dot"));
        let dot = std::fs::read_to_string(&path)?;
        let diagram = feynkit_graph::FeynmanDiagram::from_dot(Arc::clone(&model), &dot)?;
        let settings = refpath::vakint_settings(diagram.loop_count() as i64 + 1);
        let mut record = json!({"graph": graph, "reps": reps});
        for &projector in &projectors {
            let options = FromFeynkitOptions {
                projector,
                merge_colour,
                ..FromFeynkitOptions::default()
            };
            let mut entry_s: Vec<[f64; 3]> = Vec::new();
            let mut adapter_s: Vec<[f64; 3]> = Vec::new();
            let (mut entry, mut adapter) = (None, None);
            for rep in 0..reps {
                let mut run_entry = || -> eyre::Result<()> {
                    let engine = refpath::fresh_engine(&template)?;
                    let started = Instant::now();
                    let got = evaluate_graph(&engine, &settings, &diagram, &model, &options)?;
                    entry_s.push([
                        started.elapsed().as_secs_f64(),
                        got.contract_seconds,
                        got.evaluate_seconds,
                    ]);
                    entry.get_or_insert(got);
                    Ok(())
                };
                let mut run_adapter = || -> eyre::Result<()> {
                    let engine = refpath::fresh_engine(&template)?;
                    let started = Instant::now();
                    let got =
                        evaluate_graph_adapter(&engine, &settings, &diagram, &model, &options)?;
                    adapter_s.push([
                        started.elapsed().as_secs_f64(),
                        got.contract_seconds,
                        got.evaluate_seconds,
                    ]);
                    adapter.get_or_insert(got);
                    Ok(())
                };
                if rep % 2 == 0 {
                    run_entry()?;
                    run_adapter()?;
                } else {
                    run_adapter()?;
                    run_entry()?;
                }
            }
            let (entry, adapter) = (entry.unwrap(), adapter.unwrap());
            let equal = direct::values_equal(&entry.value, &adapter.value);
            if !equal {
                value_mismatch += 1;
            }
            let handed = entry.route == GraphRoute::HandOver;
            if !handed && entry.branches > 0 {
                not_handed += 1;
            }
            let medians = |runs: &[[f64; 3]]| {
                json!({
                    "wall": median(runs.iter().map(|r| r[0]).collect()),
                    "contract": median(runs.iter().map(|r| r[1]).collect()),
                    "evaluate": median(runs.iter().map(|r| r[2]).collect()),
                })
            };
            let name = match projector {
                ProjectorChoice::Metric => "metric",
                ProjectorChoice::Transverse => "transverse",
                ProjectorChoice::None => "open",
            };
            record[name] = json!({
                "values_equal": equal, "route": format!("{:?}", entry.route),
                "branches": entry.branches, "zero": adapter.value.is_zero(),
                "entry": medians(&entry_s), "adapter": medians(&adapter_s),
            });
            if projector == ProjectorChoice::Metric
                && let Some(reference) = &reference
            {
                let engine = refpath::fresh_engine(&template)?;
                let gate = dump_gate(&engine, &settings, &diagram, &model, reference, &stem)?;
                if gate["pass"] != true {
                    dump_fail += 1;
                }
                record["metric"]["dump_gate"] = gate;
            }
        }
        writeln!(results, "{record}")?;
        eprintln!("[graph_entry_bb] diagram {graph}: {record}");
    }
    eprintln!(
        "[graph_entry_bb] done: {} graphs, {value_mismatch} value mismatches, {not_handed} not handed over, {dump_fail} stored-dump failures",
        graphs.len()
    );
    if value_mismatch > 0 || dump_fail > 0 {
        eyre::bail!("{value_mismatch} value mismatches, {dump_fail} stored-dump failures");
    }
    Ok(())
}
