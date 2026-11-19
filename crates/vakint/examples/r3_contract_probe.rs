//! R3's contraction alone (helper VH, lane idenso-audit): per graph the colour
//! terms, the compact spanning tree and `refpath::direct::contract_edges` with
//! the hand-over's default options, timed, with a marker line on stderr before
//! and after each graph so that idenso's `IDENSO_FRONTIER_TRACE=1` lines can be
//! attributed to it. No match, no lowering, no reduction.
//!
//! ```text
//! r3_contract_probe <prepared dir> --graph-list <file> [--no-merge] [--no-compact] [--no-dots]
//! ```
//! Prints one JSON line per graph on stdout.
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use serde_json::json;
use symbolica::atom::{Atom, AtomCore};
use vakint::refpath::{self, Projection, direct};

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let (mut directory, mut graphs) = (None, Vec::<usize>::new());
    let mut merge = true;
    let mut options = direct::DirectOptions::default();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--graph-list" => {
                let text = std::fs::read_to_string(&args[index + 1])?;
                graphs = text
                    .split_whitespace()
                    .map(|t| t.parse().unwrap())
                    .collect();
                index += 1;
            }
            "--no-merge" => merge = false,
            "--no-compact" => options.compact = false,
            "--no-dots" => options.dots = false,
            other if directory.is_none() && !other.starts_with("--") => {
                directory = Some(PathBuf::from(other))
            }
            other => eyre::bail!("unknown argument {other}"),
        }
        index += 1;
    }
    let directory = directory.ok_or_else(|| eyre::eyre!("no prepared directory"))?;
    let model = refpath::load_model(&directory)?;
    for graph in graphs {
        if let Err(error) = probe(&directory, &model, graph, merge, options) {
            eprintln!("PROBE skip graph={graph} {error}");
        }
    }
    Ok(())
}

fn probe(
    directory: &std::path::Path,
    model: &Arc<feynkit_model::Model>,
    graph: usize,
    merge: bool,
    options: direct::DirectOptions,
) -> eyre::Result<()> {
    let s = &*refpath::SYMBOLS;
    {
        let dot = std::fs::read_to_string(directory.join(format!("diagram-{graph:06}.dot")))?;
        let diagram = feynkit_graph::FeynmanDiagram::from_dot(Arc::clone(model), &dot)?;
        if refpath::scaleless_contact(&diagram) {
            eyre::bail!("scaleless contact");
        }
        let (terms, branches) = direct::colour_terms(&diagram, merge)?;
        if terms.is_empty() {
            eyre::bail!("colour zero");
        }
        let dimension = Atom::var(s.d);
        let kinematics = feynkit_kinematics::Kinematics::in_dimension(&dimension)?;
        let family = diagram.propagator_family(&kinematics)?;
        let momentum = family.external_momenta()[0].clone();
        let q_squared = feynkit_kinematics::Kinematics::in_dimension(&dimension)?
            .with_momenta([momentum.clone()])?
            .scalar_product(&momentum, &momentum)?;
        let lines = direct::lines_of(&diagram, &momentum).map_err(|e| eyre::eyre!(e))?;
        let chosen = if options.compact {
            Some(lines.compact_basis().map_err(|e| eyre::eyre!(e))?)
        } else {
            None
        };
        eprintln!(
            "PROBE begin graph={graph} terms={} branches={branches}",
            terms.len()
        );
        let started = Instant::now();
        let numerators = direct::contract_edges(
            &diagram,
            &terms,
            &momentum,
            &q_squared,
            Projection::Metric,
            options,
            chosen.as_ref().map(|chosen| (&lines, chosen)),
        )?;
        let seconds = started.elapsed().as_secs_f64();
        let scalar_terms: usize = numerators
            .iter()
            .map(|(_, v)| match v.expand().as_view() {
                symbolica::atom::AtomView::Add(a) => a.get_nargs(),
                _ => 1,
            })
            .sum();
        eprintln!("PROBE end graph={graph} contract_s={seconds:.3}");
        println!(
            "{}",
            json!({"graph": graph, "terms": terms.len(), "branches": branches,
                "contract_s": seconds, "scalar_terms": scalar_terms})
        );
    }
    Ok(())
}
