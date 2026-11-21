//! U7a delta, lane A (accounting): the per-diagram ledger of
//! `GAMMALOOP_UNCONTRACTED_NUMERATOR` ([`super::uncontracted_ledger`]):
//!
//! * R2-I1 (helper VH round 2): strict is the only behaviour of on and
//!   check; `GAMMALOOP_UNCONTRACTED_STRICT=0` is refused;
//! * U7a.1 (Astra's U7a verdict finding 1, (iii) and (iv)): check mode's
//!   coverage rule (an
//!   all-declined diagram-kind fails unless named), Δ-cache hits as their own
//!   disposition, check mode with the PINT6 dump refused;
//! * U7a.4 (finding 4, R2-M2, R2-M3): one disposition line per diagram and
//!   kind (also when nothing is expected), the summary per diagram, the
//!   child records' evaluation path;
//!
//! The environment tests set the switch's variables in their own process
//! (cargo nextest runs every test in one; they refuse to run elsewhere).

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
    uncontracted_ledger::{self as ledger, DiagramLedger, Exception, Exceptions},
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

/// The graph name of a DOT's diagram.
fn graph_name(dot: &str) -> Result<String> {
    Ok(colour_factored(&load(dot)?)?
        .to_gamma_loop_graph(None, true)?
        .name)
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

fn agree(a: &Atom, b: &Atom) -> bool {
    a == b || (a - b).expand().is_zero() || uncontracted::agree(a, b).unwrap_or(false)
}

/// Tests that set the switch's environment: only in a process of their own.
fn own_process() -> Result<()> {
    ensure!(
        std::env::var_os("NEXTEST").is_some(),
        "this test sets process environment variables: run it under cargo nextest (one process per test)"
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// R2-I1: strict is the only behaviour
// ---------------------------------------------------------------------------

#[test]
fn strict_is_the_only_behaviour() {
    for accepted in [None, Some(""), Some("1"), Some(" on ")] {
        assert!(ledger::parse_strict(accepted).is_ok(), "{accepted:?}");
    }
    for refused in ["0", "off", " 0 ", "2", "yes"] {
        let error = ledger::parse_strict(Some(refused)).unwrap_err();
        assert!(error.contains(uncontracted::STRICT_VARIABLE), "{error}");
    }
    // refused in every mode, off included (the variable has no meaning left)
    for switch in [
        Switch::OFF,
        Switch::new(Mode::On, &ROOT_CHILD),
        Switch::new(Mode::Check, &ROOT_CHILD),
    ] {
        assert!(ledger::initialise(switch, Some("0"), || false, None).is_err());
        assert!(ledger::initialise(switch, Some("1"), || false, None).is_ok());
    }
}

/// `GAMMALOOP_UNCONTRACTED_STRICT=0` refuses the switch where it
/// initialises: the beta's start-up check and the first read (at the base,
/// on with STRICT=0 was the non-strict mode, whose errors fell back).
#[test]
fn strict_zero_in_the_environment_is_refused() -> Result<()> {
    own_process()?;
    unsafe {
        std::env::set_var(uncontracted::VARIABLE, "on:root,child");
        std::env::set_var(uncontracted::STRICT_VARIABLE, "0");
    }
    let error = super::validate_uncontracted()
        .err()
        .ok_or_else(|| eyre!("the switch initialised with STRICT=0"))?;
    ensure!(
        format!("{error}").contains("R2-I1"),
        "the refusal does not say why: {error}"
    );
    let read = std::panic::catch_unwind(uncontracted::switch);
    ensure!(read.is_err(), "the switch was read with STRICT=0: {read:?}");
    Ok(())
}

/// An error outcome in mode on fails the diagram (at the base it fell back
/// to today's path unless STRICT=1), and so in check mode; a designed
/// decline still falls back to today's path.
#[test]
fn on_mode_error_fails_the_diagram() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let on = Switch::new(Mode::On, &ROOT_CHILD);
        let off = run(L2_22, Switch::OFF, None, None)?.value?;
        let clean = run(L2_22, on, None, None)?;
        ensure!(agree(&clean.value?, &off), "on differs from off");
        for mode in [Mode::On, Mode::Check] {
            let failed = ledger::with_fault(Kind::Root, "routing", false, || {
                run(L2_22, Switch::new(mode, &ROOT_CHILD), None, None)
            })?;
            let error = match failed.value {
                Ok(value) => {
                    return Err(eyre!(
                        "{}: an error outcome took today's path: {value}",
                        mode.as_str()
                    ));
                }
                Err(error) => format!("{error:#}"),
            };
            ensure!(
                error.contains("outcome=error:routing"),
                "{}: {error}",
                mode.as_str()
            );
            let summaries = lines_of(&failed.lines, "UNCONTRACTED_SUMMARY ");
            ensure!(
                summaries.len() == 1
                    && field(summaries[0], "outcome") == Some("error")
                    && field(summaries[0], "error") == Some("routing:1"),
                "{}: the summary does not report the error: {summaries:#?}",
                mode.as_str()
            );
            let root = disposition(&failed.lines, Kind::Root)?;
            ensure!(
                field(root, "outcome") == Some("error"),
                "{}: {root}",
                mode.as_str()
            );
        }
        let declined = ledger::with_fault(Kind::Root, "zero_loops", true, || {
            run(L2_22, on, None, None)
        })?;
        ensure!(
            agree(&declined.value?, &off),
            "a designed decline must take today's path"
        );
        let root = disposition(&declined.lines, Kind::Root)?;
        ensure!(
            field(root, "outcome") == Some("ok")
                && field(root, "declined").is_some_and(|d| d.starts_with("zero_loops:")),
            "{root}"
        );
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// U7a.1 (iv): refusals at the switch's initialisation
// ---------------------------------------------------------------------------

#[test]
fn check_with_the_pint6_dump_and_inert_lists_are_refused() {
    let check = Switch::new(Mode::Check, &ROOT_CHILD);
    let error = ledger::initialise(check, None, || true, None).unwrap_err();
    assert!(error.contains("PINT6_CLOSURE_DUMP"), "{error}");
    // mode on keeps the PINT6 dump as a designed decline (today's records)
    assert!(ledger::initialise(Switch::new(Mode::On, &ROOT_CHILD), None, || true, None).is_ok());
    // a list without a kind this build emits would build no record
    assert!(
        ledger::initialise(
            Switch::new(Mode::Check, &[Kind::Quotient, Kind::Surrogate]),
            None,
            || false,
            None
        )
        .is_err()
    );
    // an exception list is read in check mode only
    assert!(
        ledger::initialise(
            Switch::new(Mode::On, &ROOT_CHILD),
            None,
            || false,
            Some("/nonexistent")
        )
        .is_err()
    );
    assert!(ledger::initialise(check, None, || false, Some("/nonexistent")).is_err());
    assert!(ledger::initialise(Switch::OFF, None, || true, Some("/nonexistent")).is_ok());
}

/// Check mode with `PINT6_CLOSURE_DUMP` from the environment: refused where
/// the switch initialises (at the base every record declined `pint6_dump`
/// and the check passed with zero comparisons).
#[test]
fn check_with_the_pint6_dump_in_the_environment_is_refused() -> Result<()> {
    own_process()?;
    let dump = std::env::temp_dir().join(format!("u7a-ledger-pint6-{}", std::process::id()));
    unsafe {
        std::env::set_var(uncontracted::VARIABLE, "check:root,child");
        std::env::set_var("PINT6_CLOSURE_DUMP", &dump);
    }
    let refused = super::validate_uncontracted();
    let read = std::panic::catch_unwind(uncontracted::switch);
    let _ = std::fs::remove_dir_all(&dump);
    let error = refused
        .err()
        .ok_or_else(|| eyre!("check mode initialised with the PINT6 dump"))?;
    ensure!(format!("{error}").contains("pint6_dump"), "{error}");
    ensure!(
        read.is_err(),
        "the switch was read with the PINT6 dump: {read:?}"
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// U7a.4: one line per diagram and kind; U7a.1 (iii) coverage
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
            let accounted = count(line, "compared") + count(line, "cache_hit");
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

/// A diagram-kind whose every record declines fails check mode (at the base
/// it passed with zero comparisons, Astra's pint6 receipt), unless the
/// run's exception list names each of its codes for that graph and kind.
#[test]
fn all_declined_kind_fails_check_unless_named() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let check = Switch::new(Mode::Check, &ROOT_CHILD);
        let graph = graph_name(L2_22)?;
        let declined = || {
            ledger::with_fault(Kind::Root, "zero_loops", true, || {
                run(L2_22, check, None, None)
            })
        };
        let failed = declined()?;
        let error = failed
            .value
            .err()
            .ok_or_else(|| eyre!("an all-declined root passed check mode"))?;
        ensure!(format!("{error:#}").contains("none compared"), "{error:#}");
        let root = disposition(&failed.lines, Kind::Root)?;
        ensure!(
            field(root, "outcome") == Some("uncovered") && count(root, "compared") == 0,
            "{root}"
        );
        let named = |graph: &str, code: &str| Exceptions {
            entries: vec![Exception {
                graph: graph.to_owned(),
                kind: Kind::Root,
                code: code.to_owned(),
                reason: "the test's designed decline".to_owned(),
            }],
            digest: Some("test".to_owned()),
        };
        // another code or another graph does not name it
        for exceptions in [
            named(&graph, "pint6_dump"),
            named("another_graph", "zero_loops"),
        ] {
            let other = ledger::with_exceptions(exceptions, declined)?;
            ensure!(
                other.value.is_err(),
                "an exception for another code or graph"
            );
        }
        let excepted = ledger::with_exceptions(named(&graph, "zero_loops"), declined)?;
        let _ = excepted.value?;
        let root = disposition(&excepted.lines, Kind::Root)?;
        ensure!(
            field(root, "outcome") == Some("excepted")
                && field(root, "excepted").is_some_and(|e| e.starts_with("zero_loops:")),
            "{root}"
        );
        // the child kind is compared as usual
        let child = disposition(&excepted.lines, Kind::Child)?;
        ensure!(
            field(child, "outcome") == Some("ok") && count(child, "compared") > 0,
            "{child}"
        );
        Ok(())
    })
}

/// Astra's all-declined CHECK receipt at the full hook: with the PINT6 dump
/// every record declines `pint6_dump`; check mode now fails the diagram (at
/// the base it returned today's value with `compared=0 ... outcome=ok`;
/// the environment's check switch refuses the combination outright, above).
#[test]
fn pint6_dump_check_with_every_record_declined_fails() -> Result<()> {
    own_process()?;
    let dump = std::env::temp_dir().join(format!("u7a-ledger-declined-{}", std::process::id()));
    unsafe {
        std::env::set_var("PINT6_CLOSURE_DUMP", &dump);
    }
    let outcome = with_production_piece_limit(|| -> Result<Run> {
        test_initialise()?;
        run(L2_22, Switch::new(Mode::Check, &ROOT_CHILD), None, None)
    });
    let _ = std::fs::remove_dir_all(&dump);
    let declined = outcome?;
    ensure!(
        declined
            .captured
            .iter()
            .all(|attempt| attempt.outcome == "declined:pint6_dump")
            && !declined.captured.is_empty(),
        "every record declines pint6_dump"
    );
    ensure!(
        declined.value.is_err(),
        "a check with zero comparisons passed"
    );
    for kind in ROOT_CHILD {
        let line = disposition(&declined.lines, kind)?;
        ensure!(
            field(line, "outcome") == Some("uncovered")
                && count(line, "compared") == 0
                && field(line, "declined").is_some_and(|d| d.starts_with("pint6_dump:")),
            "{line}"
        );
    }
    Ok(())
}

/// A Δ-cache hit is its own disposition: the hit counts the fill's records
/// with the fill's dispositions (expected does not depend on the cache
/// order; at the base a hit skipped them silently); a hit whose fill left no
/// note (an off fill) is unchecked and fails check mode.
#[test]
fn delta_cache_hits_count_the_fill() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let check = Switch::new(Mode::Check, &ROOT_CHILD);
        let cache = Arc::new(CountertermCache::default());
        let fill = run(L2_22, check, Some(cache.clone()), None)?;
        let _ = fill.value?;
        let hit = run(L2_22, check, Some(cache), None)?;
        let _ = hit.value?;
        let filled = disposition(&fill.lines, Kind::Child)?;
        let served = disposition(&hit.lines, Kind::Child)?;
        ensure!(count(filled, "compared") > 0, "{filled}");
        ensure!(
            count(served, "expected") == count(filled, "expected")
                && count(served, "attempted") == 0
                && count(served, "cache_hit") == count(filled, "expected")
                && field(served, "cache_hit_fill")
                    == Some(format!("compared:{}", count(filled, "expected")).as_str())
                && field(served, "outcome") == Some("ok"),
            "{filled}\n{served}"
        );
        let hits = lines_of(&hit.lines, "UNCONTRACTED_CACHE_HIT ");
        ensure!(
            !hits.is_empty()
                && hits
                    .iter()
                    .all(|line| field(line, "fill_diagram").is_some_and(|d| d.starts_with("root#"))),
            "{:#?}",
            hit.lines
        );
        // the root kind never hits (a root value never reaches the cache)
        ensure!(count(disposition(&hit.lines, Kind::Root)?, "cache_hit") == 0);
        // a cache an off evaluation filled (helper VH R2-M1): loud, not silent
        let mixed = Arc::new(CountertermCache::default());
        let _ = run(L2_22, Switch::OFF, Some(mixed.clone()), None)?.value?;
        let unchecked = run(L2_22, check, Some(mixed), None)?;
        ensure!(
            unchecked.value.is_err(),
            "a check served by an off fill passed"
        );
        let served = disposition(&unchecked.lines, Kind::Child)?;
        ensure!(
            field(served, "cache_hit_fill").is_some_and(|f| f.contains("unchecked"))
                && field(served, "outcome") == Some("uncovered"),
            "{served}"
        );
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
