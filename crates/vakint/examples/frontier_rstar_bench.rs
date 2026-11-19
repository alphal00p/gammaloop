//! idenso's contraction on the beta thread's readable R* dumps (helper VH,
//! lane idenso-audit): to compare idenso's frontier orders
//! (`AlgebraSettings::frontier_order`) on the R* inputs, per item the median
//! seconds over `--reps`, a hash of the expanded result's canonical form and
//! idenso's certificates of the result (closed scalar, terminal polynomial).
//! The contraction is the R* production chain of
//! `p_integral::factorized::contract_with_idenso` (cook with reversible
//! encoding and payloads, `infer_indexed`, `simplify_algebra` with gamma
//! `canonical()` and epsilon, `to_dots`, uncook, expand) with the chosen order.
//! Read-only on its inputs.
//!
//! ```text
//! frontier_rstar_bench records <record file>...   pint6.factorized-numerator.v1 records: projector . factors per record
//! frontier_rstar_bench whole <dump>...            the whole expression
//! frontier_rstar_bench summands <dump> [from to]  per summand (open slots kept)
//!     [--reps <n>] [--order boundary|close|merge|merge:<cap>]
//! ```
//! A dump is a canonical text (`.txt`) or a Symbolica binary export. One JSON
//! line per item on stdout.
use std::time::Instant;

use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, FrontierOrder, SymbolicTensor};
use idenso::{CookMode, CookSettings};
use serde_json::json;
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::p_integral::factorized::read_factorized_records;

fn import(path: &str) -> eyre::Result<Atom> {
    Ok(if path.ends_with(".txt") {
        let text = std::fs::read_to_string(path)?;
        symbolica::try_parse!(text.trim()).map_err(|e| eyre::eyre!("parse {path}: {e}"))?
    } else {
        let mut reader = std::io::BufReader::new(std::fs::File::open(path)?);
        Atom::import(&mut reader, None).map_err(|e| eyre::eyre!("import {path}: {e}"))?
    })
}

fn parse_order(text: &str) -> eyre::Result<FrontierOrder> {
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

/// The contraction and idenso's certificates of its result.
struct Contracted {
    value: Atom,
    closed_scalar: bool,
    terminal_polynomial: bool,
}

fn contract(expression: &Atom, order: FrontierOrder) -> eyre::Result<Contracted> {
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking.try_cook(expression.as_view())?;
    let simplified = SymbolicTensor::infer_indexed(cooked)?.simplify_algebra(&AlgebraSettings {
        gamma: Some(GammaSimplifySettings::canonical()),
        epsilon: true,
        frontier_order: order,
        ..Default::default()
    })?;
    let closed_scalar = simplified.is_certified_closed_scalar();
    let terminal_polynomial = simplified.is_certified_terminal_polynomial();
    let reduced = simplified.to_dots()?;
    Ok(Contracted {
        value: cooking.uncook(reduced.expression().as_view()).expand(),
        closed_scalar,
        terminal_polynomial,
    })
}

fn hash_terms(value: &Atom) -> (String, usize) {
    let terms = match value.as_view() {
        AtomView::Add(sum) => sum.get_nargs(),
        _ => 1,
    };
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    std::hash::Hash::hash(&value.to_canonical_string(), &mut hasher);
    (
        format!("{:016x}", std::hash::Hasher::finish(&hasher)),
        terms,
    )
}

fn run(item: &str, expression: &Atom, order: FrontierOrder, reps: usize) -> eyre::Result<()> {
    let mut seconds = Vec::new();
    let mut last = None;
    for _ in 0..reps.max(1) {
        let started = Instant::now();
        last = Some(contract(expression, order)?);
        seconds.push(started.elapsed().as_secs_f64());
    }
    seconds.sort_by(f64::total_cmp);
    let contracted = last.unwrap();
    let (hash, terms) = hash_terms(&contracted.value);
    println!(
        "{}",
        json!({"item": item, "contract_s": seconds[seconds.len() / 2], "terms": terms, "hash": hash,
            "closed_scalar": contracted.closed_scalar,
            "terminal_polynomial": contracted.terminal_polynomial,
            "order": format!("{order:?}")})
    );
    Ok(())
}

fn main() -> eyre::Result<()> {
    let mut args = std::env::args().skip(1).collect::<Vec<_>>();
    let mut reps = 1usize;
    let mut order = FrontierOrder::Boundary;
    if let Some(i) = args.iter().position(|a| a == "--reps") {
        reps = args[i + 1].parse()?;
        args.drain(i..=i + 1);
    }
    if let Some(i) = args.iter().position(|a| a == "--order") {
        order = parse_order(&args[i + 1])?;
        args.drain(i..=i + 1);
    }
    let mode = args.first().cloned().unwrap_or_default();
    match mode.as_str() {
        "records" => {
            for path in &args[1..] {
                let records = read_factorized_records(std::path::Path::new(path))?;
                for (i, record) in records.iter().enumerate() {
                    let mut product = record.factor_product();
                    if let Some(projector) = &record.projector {
                        product *= &projector.tensor;
                    }
                    run(&format!("{path}#{i}"), &product, order, reps)?;
                }
            }
        }
        "whole" => {
            for path in &args[1..] {
                run(path, &import(path)?, order, reps)?;
            }
        }
        "summands" => {
            let path = &args[1];
            let atom = import(path)?;
            let summands = match atom.as_view() {
                AtomView::Add(sum) => sum.iter().map(|s| s.to_owned()).collect::<Vec<_>>(),
                single => vec![single.to_owned()],
            };
            let from = args.get(2).map_or(Ok(0), |a| a.parse())?;
            let to = args.get(3).map_or(Ok(summands.len()), |a| a.parse())?;
            for (i, summand) in summands
                .iter()
                .enumerate()
                .take(to.min(summands.len()))
                .skip(from)
            {
                run(&format!("{path}#{i}"), summand, order, reps)?;
            }
        }
        other => eyre::bail!("unknown mode {other:?} (records | whole | summands)"),
    }
    Ok(())
}
