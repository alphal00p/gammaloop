//! R3's contraction alone (helper VH, lane idenso-audit): per graph the colour
//! terms, the compact spanning tree and `refpath::direct::contract_edges` with
//! the hand-over's default options, timed, with a marker line on stderr before
//! and after each graph so that idenso's `IDENSO_FRONTIER_TRACE=1` lines can be
//! attributed to it. No match, no lowering, no reduction.
//!
//! ```text
//! r3_contract_probe <prepared dir> --graph-list <file> [--no-merge] [--no-compact] [--no-dots] [--reps <n>]
//!     [--order boundary|close|merge|merge:<cap>] [--orders <order>,<order>,..]
//!     [--closings off,on]  (idenso's frontier-level closing, DirectOptions::closing)
//! (`--orders` x `--closings`: the variants interleaved per repetition, one line per graph and variant)
//! ```
//! Prints one JSON line per graph on stdout: the median contraction seconds
//! over `--reps`, the expanded scalar's term count and a hash of its canonical
//! form (to compare idenso configurations, e.g. the frontier order
//! `DirectOptions::frontier_order`).
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use serde_json::json;
use symbolica::atom::{Atom, AtomCore};
use vakint::refpath::{self, Projection, direct};

fn parse_order(text: &str) -> eyre::Result<idenso::tensor::FrontierOrder> {
    use idenso::tensor::FrontierOrder;
    Ok(match text {
        "boundary" | "off" => FrontierOrder::Boundary,
        "close" => FrontierOrder::Close,
        "merge" => FrontierOrder::Merge { cap: None },
        other => match other.strip_prefix("merge:") {
            Some(cap) => FrontierOrder::Merge {
                cap: Some(cap.parse()?),
            },
            None => eyre::bail!("unknown order {other}"),
        },
    })
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let (mut directory, mut graphs) = (None, Vec::<usize>::new());
    let mut merge = true;
    let mut reps = 1usize;
    let mut orders = Vec::new();
    let mut closings: Vec<bool> = Vec::new();
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
            "--reps" => {
                reps = args[index + 1].parse()?;
                index += 1;
            }
            "--order" => {
                options.frontier_order = parse_order(&args[index + 1])?;
                index += 1;
            }
            "--closings" => {
                closings = args[index + 1]
                    .split(',')
                    .map(|c| match c {
                        "off" => Ok(false),
                        "on" => Ok(true),
                        other => Err(eyre::eyre!("unknown closing {other}")),
                    })
                    .collect::<eyre::Result<_>>()?;
                index += 1;
            }
            "--orders" => {
                orders = args[index + 1]
                    .split(',')
                    .map(parse_order)
                    .collect::<eyre::Result<_>>()?;
                index += 1;
            }
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
        let orders = if orders.is_empty() {
            vec![options.frontier_order]
        } else {
            orders.clone()
        };
        let closings = if closings.is_empty() {
            vec![options.closing]
        } else {
            closings.clone()
        };
        let variants = closings
            .iter()
            .flat_map(|&closing| orders.iter().map(move |&order| (order, closing)))
            .collect::<Vec<_>>();
        if let Err(error) = probe(&directory, &model, graph, merge, options, reps, &variants) {
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
    reps: usize,
    orders: &[(idenso::tensor::FrontierOrder, bool)],
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
        // The orders interleaved per repetition (rotated, so that none always
        // runs first): paired timings robust to the host's load drifting.
        let mut seconds = vec![Vec::with_capacity(reps.max(1)); orders.len()];
        let mut numerators = vec![Vec::new(); orders.len()];
        for rep in 0..reps.max(1) {
            for step in 0..orders.len() {
                let which = (step + rep) % orders.len();
                let mut options = options;
                (options.frontier_order, options.closing) = orders[which];
                let started = Instant::now();
                numerators[which] = direct::contract_edges(
                    &diagram,
                    &terms,
                    &momentum,
                    &q_squared,
                    Projection::Metric,
                    options,
                    chosen.as_ref().map(|chosen| (&lines, chosen)),
                )?;
                seconds[which].push(started.elapsed().as_secs_f64());
            }
        }
        for (which, (order, closing)) in orders.iter().enumerate() {
            let seconds = &mut seconds[which];
            seconds.sort_by(f64::total_cmp);
            let median = seconds[seconds.len() / 2];
            let mut scalar_terms = 0usize;
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            for (name, value) in &numerators[which] {
                let expanded = value.expand();
                scalar_terms += match expanded.as_view() {
                    symbolica::atom::AtomView::Add(a) => a.get_nargs(),
                    _ => 1,
                };
                std::hash::Hash::hash(name, &mut hasher);
                std::hash::Hash::hash(&expanded.to_canonical_string(), &mut hasher);
            }
            let hash = format!("{:016x}", std::hash::Hasher::finish(&hasher));
            eprintln!(
                "PROBE end graph={graph} order={order:?} closing={closing} contract_s={median:.3}"
            );
            println!(
                "{}",
                json!({"graph": graph, "terms": terms.len(), "branches": branches,
                    "contract_s": median, "contract_all_s": seconds, "scalar_terms": scalar_terms,
                    "hash": hash, "order": format!("{order:?}"), "closing": closing})
            );
        }
    }
    Ok(())
}
