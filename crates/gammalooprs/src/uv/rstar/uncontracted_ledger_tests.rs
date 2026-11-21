//! U7a delta, lane A (accounting): the per-diagram ledger of
//! `GAMMALOOP_UNCONTRACTED_NUMERATOR` ([`super::uncontracted_ledger`]):
//!
//! * U7a.4 (Astra's U7a verdict finding 4, R2-M2, R2-M3): one disposition line per diagram and
//!   kind (also when nothing is expected), the summary per diagram, the
//!   child records' evaluation path;

use std::{collections::BTreeSet, sync::Arc};

use eyre::{Result, ensure, eyre};
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
    symbol,
};
use vakint::{Vakint, p_integral::factorized::FactorizedNumerator};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    utils::{GS, W_},
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    CountertermCache, RStarPlan,
    boundary::{BoundaryMode, ClosedRootNormalization, with_boundary_mode},
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{self, Captured, Kind, Mode, Switch, with_capture, with_switch},
    uncontracted_ledger::{self as ledger, DiagramLedger},
};

const L2_22: &str = include_str!("testdata/u7a/l2-d22.dot");
const L3_36: &str = include_str!("testdata/u7a/l3-d036.dot");
const ROOT_CHILD: [Kind; 2] = [Kind::Root, Kind::Child];
const ALL_KINDS: [Kind; 3] = [Kind::Root, Kind::Child, Kind::IrHard];

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// The beta's closed metric-root normalisation (as `uncontracted_tests`).
fn normalization(graph: &Graph) -> ClosedRootNormalization {
    let external = graph
        .iter_edges()
        .filter_map(|(pair, edge, _)| match pair {
            HedgePair::Unpaired { flow, .. } => Some((edge, flow)),
            _ => None,
        })
        .collect::<Vec<_>>();
    let (p, p_flow) = external[0];
    let external_routing = external
        .iter()
        .map(|(edge, flow)| {
            let sign = if *edge == p || *flow != p_flow { 1 } else { -1 };
            Replacement::new(
                function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                (Atom::num(sign) * function!(GS.emr_mom, usize::from(p), W_.x___)).to_pattern(),
            )
        })
        .collect::<Vec<_>>();
    let x0: EdgeIndex = p;
    let mink = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
    let p2 = function!(
        SPENSO_TAG.dot,
        function!(GS.emr_mom, usize::from(x0), mink.clone()),
        function!(GS.emr_mom, usize::from(x0), mink)
    );
    let factor = Atom::var(symbol!("beta_projection::transverse"))
        * graph
            .global_prefactor
            .num
            .replace_multiple(&external_routing)
        / (Atom::i() * p2)
        / (Atom::var(GS.dim) - 1);
    ClosedRootNormalization {
        external: x0,
        external_routing,
        factor,
    }
}

/// One run: the root's value, its attempts and the ledger's lines printed
/// inside it.
struct Run {
    value: Result<Atom>,
    captured: Vec<Captured>,
    lines: Vec<String>,
}

/// `evaluations` diagram roots (the beta's closed metric root) of `dot`
/// under `switch`, the boundary path on (production), one fresh plan;
/// `cache`: a Δ cache shared with other runs; `diagram`: a diagram ledger
/// opened around them as the beta opens one per diagram job (else each root
/// evaluation opens its own).
fn run_n(
    dot: &str,
    switch: Switch,
    cache: Option<Arc<CountertermCache>>,
    diagram: Option<&str>,
    evaluations: usize,
) -> Result<Run> {
    let uv = settings();
    let loaded = colour_factored(&load(dot)?)?;
    let graph = loaded.to_gamma_loop_graph(None, true)?;
    let model = loaded.model();
    let mut plan = RStarPlan::with_automatic_edge(&graph, model, &uv)?;
    if let Some(cache) = cache {
        plan = plan.with_counterterm_cache(cache);
    }
    let closed_root = normalization(&graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    let vakint = Vakint::new()?;
    let _ = ledger::take_lines();
    let _ = uncontracted::take_summaries();
    let (value, captured) = with_switch(switch, || {
        with_boundary_mode(BoundaryMode::On, false, || {
            with_capture(|| {
                let opened = diagram.and_then(|key| DiagramLedger::open(&graph.name, key));
                let mut value = Ok(Atom::Zero);
                for _ in 0..evaluations {
                    value = plan.projected_counterterm_closing_root_scalar(
                        &graph,
                        model,
                        &vakint,
                        &uv,
                        &callback,
                        &closed_root,
                    );
                    if value.is_err() {
                        break;
                    }
                }
                match opened {
                    None => value,
                    Some(opened) => {
                        let closed = opened.close(value.is_ok());
                        match (value, closed) {
                            (Ok(value), Ok(())) => Ok(value),
                            (Err(error), _) | (Ok(_), Err(error)) => Err(error),
                        }
                    }
                }
            })
        })
    });
    let _ = uncontracted::take_summaries();
    Ok(Run {
        value,
        captured,
        lines: ledger::take_lines(),
    })
}

fn run(
    dot: &str,
    switch: Switch,
    cache: Option<Arc<CountertermCache>>,
    diagram: Option<&str>,
) -> Result<Run> {
    run_n(dot, switch, cache, diagram, 1)
}

/// The `<key>=` value of a log line.
fn field<'l>(line: &'l str, key: &str) -> Option<&'l str> {
    line.split_whitespace()
        .find_map(|token| token.strip_prefix(key)?.strip_prefix('='))
}

fn count(line: &str, key: &str) -> usize {
    field(line, key)
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no count {key}= in {line}"))
}

fn lines_of<'l>(lines: &'l [String], head: &str) -> Vec<&'l str> {
    lines
        .iter()
        .filter(|line| line.starts_with(head))
        .map(String::as_str)
        .collect()
}

/// The disposition line of `kind` (exactly one).
fn disposition<'l>(lines: &'l [String], kind: Kind) -> Result<&'l str> {
    let found = lines_of(lines, "UNCONTRACTED_DISPOSITION ")
        .into_iter()
        .filter(|line| field(line, "kind") == Some(kind.as_str()))
        .collect::<Vec<_>>();
    ensure!(
        found.len() == 1,
        "{} disposition lines of kind {}: {lines:#?}",
        found.len(),
        kind.as_str()
    );
    Ok(found[0])
}

// ---------------------------------------------------------------------------
// U7a.4: one line per diagram and kind
// ---------------------------------------------------------------------------

/// One disposition line per enabled kind and one summary per diagram, also
/// for a kind with nothing expected (at the base a diagram without
/// attempts printed nothing and every root evaluation its own summary);
/// two root evaluations in one diagram ledger give one line per kind.
#[test]
fn one_line_per_diagram_and_kind() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let check = Switch::new(Mode::Check, &ALL_KINDS);
        let single = run(L2_22, check, None, None)?;
        let _ = single.value?;
        for kind in ALL_KINDS {
            let line = disposition(&single.lines, kind)?;
            let accounted = count(line, "compared");
            ensure!(
                field(line, "outcome") == Some("ok")
                    && field(line, "declined") == Some("none")
                    && accounted == count(line, "expected"),
                "{line}"
            );
            ensure!(
                field(line, "diagram").is_some_and(|key| key.starts_with("root#")),
                "{line}"
            );
        }
        let ir_hard = disposition(&single.lines, Kind::IrHard)?;
        ensure!(count(ir_hard, "expected") == 0, "{ir_hard}");
        let summaries = lines_of(&single.lines, "UNCONTRACTED_SUMMARY ");
        ensure!(
            summaries.len() == 1
                && count(summaries[0], "attempted") == single.captured.len()
                && field(summaries[0], "kinds") == Some("root,child,ir-hard")
                && field(summaries[0], "outcome") == Some("ok"),
            "{summaries:#?}"
        );
        // a kind without records in the diagram: its lines still exist
        let none = run(L2_22, Switch::new(Mode::Check, &[Kind::IrHard]), None, None)?;
        let _ = none.value?;
        let only = lines_of(&none.lines, "UNCONTRACTED_DISPOSITION ");
        ensure!(
            only.len() == 1 && count(only[0], "expected") == 0,
            "{only:#?}"
        );
        ensure!(lines_of(&none.lines, "UNCONTRACTED_SUMMARY ").len() == 1);
        // two root evaluations of one diagram (as an insertion surrogate's or
        // a colour branch's next to the root): one ledger, one line per kind
        let twice = run_n(L2_22, check, None, Some("000022"), 2)?;
        let _ = twice.value?;
        ensure!(
            lines_of(&twice.lines, "UNCONTRACTED_SUMMARY ").len() == 1,
            "{:#?}",
            twice.lines
        );
        for kind in ROOT_CHILD {
            let line = disposition(&twice.lines, kind)?;
            ensure!(field(line, "diagram") == Some("000022"), "{line}");
            ensure!(
                count(line, "expected") == 2 * count(disposition(&single.lines, kind)?, "expected"),
                "{line}"
            );
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// U7a.4 (R2-M3)
// ---------------------------------------------------------------------------

/// A child record's label names its evaluation path, the forest of its
/// parents (helper VH R2-M3); a root record's path is `root`.
#[test]
fn child_records_name_their_parent_forest() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let checked = run(L3_36, Switch::new(Mode::Check, &ROOT_CHILD), None, None)?;
        let _ = checked.value?;
        let mut kinds = BTreeSet::new();
        for attempt in &checked.captured {
            let Some(text) = &attempt.text else { continue };
            let record = FactorizedNumerator::from_text(text).map_err(|e| eyre!("{e}"))?;
            let path =
                field(&record.label, "path").ok_or_else(|| eyre!("no path in {}", record.label))?;
            match attempt.kind {
                Kind::Root => ensure!(path == "root", "{}", record.label),
                Kind::Child => ensure!(
                    path.starts_with("root/child{") && path.ends_with('}'),
                    "{}",
                    record.label
                ),
                _ => {}
            }
            kinds.insert(attempt.kind);
        }
        ensure!(kinds.contains(&Kind::Root) && kinds.contains(&Kind::Child));
        Ok(())
    })
}

/// A K3 term that builds no record by design (a zero-loop term of an IR hard
/// part: Δ(γ) inserted, no vakint call today) is listed as `no_record=` and
/// is not expected: it neither counts as a record nor as a decline.
#[test]
fn k3_terms_without_a_record_are_not_expected() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let checked = run(L3_36, Switch::new(Mode::Check, &ALL_KINDS), None, None)?;
        let _ = checked.value?;
        let line = disposition(&checked.lines, Kind::IrHard)?;
        let records = checked
            .captured
            .iter()
            .filter(|attempt| attempt.kind == Kind::IrHard)
            .count();
        ensure!(
            field(line, "no_record").is_some_and(|n| n.starts_with("zero_loops:"))
                && count(line, "expected") == records
                && count(line, "compared") == records
                && records > 0
                && field(line, "declined") == Some("none")
                && field(line, "outcome") == Some("ok"),
            "{records} K3 records: {line}"
        );
        Ok(())
    })
}
