//! U7a delta, lane A (accounting): the per-diagram ledger of
//! `GAMMALOOP_UNCONTRACTED_NUMERATOR` ([`super::uncontracted_ledger`]):
//!
//! * R2-I1 (helper VH round 2): strict is the only behaviour of on and
//!   check; `GAMMALOOP_UNCONTRACTED_STRICT=0` is refused;
//! * U7a.1 (Astra's U7a verdict finding 1): the switch's checkpoint
//!   identity, the resume evidence, check mode's coverage rule (an
//!   all-declined diagram-kind fails unless named), Δ-cache hits as their own
//!   disposition, check mode with the PINT6 dump refused;
//! * U7a.4 (finding 4, R2-M2, R2-M3): one disposition line per diagram and
//!   kind (also when nothing is expected), the summary per diagram, the
//!   child records' evaluation path;
//! * the timing note: unique occurrence keys that pair an off and an on run.
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
    uncontracted_ledger::{self as ledger, DiagramLedger, Exception, Exceptions, Resumed},
};

pub(super) const L2_22: &str = include_str!("testdata/u7a/l2-d22.dot");
pub(super) const L3_36: &str = include_str!("testdata/u7a/l3-d036.dot");
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

/// One run: the root's value, its attempts, the ledger's lines and the
/// timing lines printed inside it, and the dispositions of the diagram
/// ledger it was given (`diagram`).
pub(super) struct Run {
    pub(super) value: Result<Atom>,
    pub(super) captured: Vec<Captured>,
    pub(super) lines: Vec<String>,
    pub(super) timing: Vec<String>,
    pub(super) evidence: Option<serde_json::Value>,
}

/// The graph name of a DOT's diagram.
pub(super) fn graph_name(dot: &str) -> Result<String> {
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
    let _ = (ledger::take_lines(), ledger::take_timing_lines());
    let _ = uncontracted::take_summaries();
    let ((value, evidence), captured) = with_switch(switch, || {
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
                    None => (value, None),
                    Some(opened) => {
                        let closed = opened.close(value.is_ok());
                        match (value, closed) {
                            (Ok(value), Ok(evidence)) => (Ok(value), evidence),
                            (Err(error), _) | (Ok(_), Err(error)) => (Err(error), None),
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
        timing: ledger::take_timing_lines(),
        evidence,
    })
}

pub(super) fn run(
    dot: &str,
    switch: Switch,
    cache: Option<Arc<CountertermCache>>,
    diagram: Option<&str>,
) -> Result<Run> {
    run_n(dot, switch, cache, diagram, 1)
}

/// The `<key>=` value of a log line.
pub(super) fn field<'l>(line: &'l str, key: &str) -> Option<&'l str> {
    line.split_whitespace()
        .find_map(|token| token.strip_prefix(key)?.strip_prefix('='))
}

pub(super) fn count(line: &str, key: &str) -> usize {
    field(line, key)
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no count {key}= in {line}"))
}

pub(super) fn lines_of<'l>(lines: &'l [String], head: &str) -> Vec<&'l str> {
    lines
        .iter()
        .filter(|line| line.starts_with(head))
        .map(String::as_str)
        .collect()
}

/// The disposition line of `kind` (exactly one).
pub(super) fn disposition<'l>(lines: &'l [String], kind: Kind) -> Result<&'l str> {
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

pub(super) fn agree(a: &Atom, b: &Atom) -> bool {
    a == b || (a - b).expand().is_zero() || uncontracted::agree(a, b).unwrap_or(false)
}

/// [`ledger::initialise`] with the production routes.
fn init(
    switch: Switch,
    strict: Option<&str>,
    timing: Option<&str>,
    pint6: impl FnOnce() -> bool,
    exceptions: Option<&str>,
) -> std::result::Result<Exceptions, String> {
    ledger::initialise(
        switch,
        strict,
        timing,
        || ledger::Routes::PRODUCTION,
        pint6,
        exceptions,
    )
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
    // refused with the switch on or check; off reads nothing but the timing
    // switch, as before the switch existed (review of the U7a delta,
    // values-and-off finding 2)
    for switch in [
        Switch::new(Mode::On, &ROOT_CHILD),
        Switch::new(Mode::Check, &ROOT_CHILD),
    ] {
        assert!(init(switch, Some("0"), None, || false, None).is_err());
        assert!(init(switch, Some("1"), None, || false, None).is_ok());
    }
    assert!(init(Switch::OFF, Some("0"), None, || false, None).is_ok());
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
        // every root record declined: mode on applies the all-declined rule
        // (pre-review of the U7a delta, finding 7), so the decline is named
        let named = Exceptions {
            entries: vec![Exception {
                graph: graph_name(L2_22)?,
                kind: Kind::Root,
                code: "zero_loops".into(),
                reason: "the test's designed decline".into(),
            }],
            digest: Some("test".into()),
        };
        let declined = ledger::with_exceptions(named, || {
            ledger::with_fault(Kind::Root, "zero_loops", true, || {
                run(L2_22, on, None, None)
            })
        })?;
        ensure!(
            agree(&declined.value?, &off),
            "a designed decline must take today's path"
        );
        let root = disposition(&declined.lines, Kind::Root)?;
        ensure!(
            field(root, "outcome") == Some("excepted")
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
    let error = init(check, None, None, || true, None).unwrap_err();
    assert!(error.contains("PINT6_CLOSURE_DUMP"), "{error}");
    // mode on too (pre-review of the U7a delta, finding 7: every record
    // would decline, a setting-wide fallback)
    let error = init(
        Switch::new(Mode::On, &ROOT_CHILD),
        None,
        None,
        || true,
        None,
    )
    .unwrap_err();
    assert!(error.contains("PINT6_CLOSURE_DUMP"), "{error}");
    // a list without a kind this build emits would build no record
    assert!(
        init(
            Switch::new(Mode::Check, &[Kind::Quotient, Kind::Insertion]),
            None,
            None,
            || false,
            None
        )
        .is_err()
    );
    // an exception list is read in on and check (a missing file refused)
    assert!(
        init(
            Switch::new(Mode::On, &ROOT_CHILD),
            None,
            None,
            || false,
            Some("/nonexistent")
        )
        .is_err()
    );
    assert!(init(check, None, None, || false, Some("/nonexistent")).is_err());
    assert!(init(Switch::OFF, None, None, || true, Some("/nonexistent")).is_ok());
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
// U7a.1 (i) and (ii): the checkpoint identity and the resume evidence
// ---------------------------------------------------------------------------

#[test]
fn identity_names_mode_kinds_and_strictness() {
    let none = Exceptions::default();
    assert_eq!(ledger::identity_of(Switch::OFF, &none), None);
    assert_eq!(
        ledger::identity_of(Switch::new(Mode::Check, &ROOT_CHILD), &none).as_deref(),
        Some("check:root,child:strict")
    );
    assert_eq!(
        ledger::identity_of(Switch::new(Mode::On, &[Kind::Root]), &none).as_deref(),
        Some("on:root:strict")
    );
    // inert kinds do not change what is built
    assert_eq!(
        ledger::identity_of(Switch::new(Mode::On, &[Kind::Root, Kind::Quotient]), &none),
        ledger::identity_of(Switch::new(Mode::On, &[Kind::Root]), &none)
    );
    // L4 49 id 8, the single-propagator shape: the emitter's own decline
    // `one_loop_tree_cut` (lane B; at the base vakint's `all_massive`)
    let named = ledger::parse_exceptions(
        "# graph\tkind\tcode\treason\nB_self_energy_L4_0049\troot\tone_loop_tree_cut\ta massive tadpole\n",
    )
    .unwrap();
    assert_eq!(
        named.entries,
        vec![Exception {
            graph: "B_self_energy_L4_0049".into(),
            kind: Kind::Root,
            code: "one_loop_tree_cut".into(),
            reason: "a massive tadpole".into(),
        }]
    );
    for designed in [
        "vakint:taylor_degree",
        "pint6_dump",
        "insertion_mass_symbol",
    ] {
        let text = format!("g\tchild\t{designed}\twhy\n");
        assert!(ledger::parse_exceptions(&text).is_ok(), "{designed}");
    }
    // an error is never excepted: a code that is not a designed refusal
    // (requirement 13's `vakint:all_massive`, lane B's `routing`, `crown`,
    // `sigma`, a vakint variant) is refused, never silently inert
    for error in [
        "vakint:all_massive",
        "vakint:Closure",
        "routing",
        "crown",
        "sigma",
        "unknown_masters",
    ] {
        let text = format!("g\troot\t{error}\twhy\n");
        let refused = ledger::parse_exceptions(&text).unwrap_err();
        assert!(refused.contains("not a designed refusal"), "{refused}");
    }
    let check = ledger::identity_of(Switch::new(Mode::Check, &ROOT_CHILD), &named).unwrap();
    assert!(
        check.starts_with("check:root,child:strict:exceptions=") && check.len() > 40,
        "{check}"
    );
    // on reads the list too (its all-declined rule, pre-review finding 7)
    let on = ledger::identity_of(Switch::new(Mode::On, &ROOT_CHILD), &named).unwrap();
    assert!(on.starts_with("on:root,child:strict:exceptions="), "{on}");
    for bad in [
        "g\troot\tcode\n",
        "g\tnot-a-kind\tcode\twhy\n",
        "g\troot\t\twhy\n",
        "g\troot\tcode\t \n",
    ] {
        assert!(ledger::parse_exceptions(bad).is_err(), "{bad:?}");
    }
}

/// A reused record is reported checked only with matching dispositions;
/// without them (an off record, another mode's, a failed one) it is
/// `outcome=unchecked`, never ok.
#[test]
fn resumed_records_need_matching_evidence() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let check = Switch::new(Mode::Check, &ROOT_CHILD);
        let on = Switch::new(Mode::On, &ROOT_CHILD);
        let fresh = run(L2_22, check, None, Some("000022"))?;
        let _ = fresh.value?;
        let evidence = fresh
            .evidence
            .ok_or_else(|| eyre!("a closed check ledger returns its dispositions"))?;
        ensure!(
            evidence["schema"] == ledger::EVIDENCE_SCHEMA
                && evidence["identity"] == "check:root,child:strict"
                && evidence["outcome"] == "ok",
            "{evidence:#}"
        );
        let on_evidence = run(L2_22, on, None, Some("000022"))?
            .evidence
            .ok_or_else(|| eyre!("no on-mode dispositions"))?;
        let mut failed = evidence.clone();
        failed["outcome"] = serde_json::json!("error");
        let graph = graph_name(L2_22)?;
        with_switch(check, || -> Result<()> {
            let _ = ledger::take_lines();
            ensure!(ledger::resumed(&graph, "000022", Some(&evidence)) == Resumed::Checked);
            let lines = ledger::take_lines();
            ensure!(
                lines
                    .iter()
                    .all(|line| line.ends_with(" source=checkpoint"))
            );
            for kind in ROOT_CHILD {
                let saved = disposition(&lines, kind)?;
                let original = disposition(&fresh.lines, kind)?;
                ensure!(
                    count(saved, "expected") == count(original, "expected")
                        && count(saved, "compared") == count(original, "compared")
                        && field(saved, "outcome") == Some("ok"),
                    "{saved}\n{original}"
                );
            }
            for (saved, reason) in [
                (None, "no_evidence"),
                (Some(&on_evidence), "identity:on:root,child:strict"),
                (Some(&failed), "outcome:error"),
            ] {
                ensure!(
                    ledger::resumed(&graph, "000022", saved)
                        == Resumed::Unchecked(reason.to_owned()),
                    "{reason}"
                );
                // one line: the beta recomputes the diagram (pre-review of
                // the U7a delta, finding 2)
                let lines = ledger::take_lines();
                ensure!(
                    lines.len() == 1
                        && lines[0].starts_with("UNCONTRACTED_RESUME ")
                        && field(&lines[0], "outcome") == Some("recompute")
                        && field(&lines[0], "reason") == Some(reason),
                    "{reason}: {lines:#?}"
                );
            }
            Ok(())
        })?;
        // off: nothing to say
        ensure!(ledger::resumed(&graph, "000022", None) == Resumed::Off);
        Ok(())
    })
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
// U7a.4 (R2-M3) and the timing note
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

/// The occurrence keys of `GAMMALOOP_UNCONTRACTED_TIMING`: unique within a
/// run and the same in an off and an on run of one diagram, so that every
/// occurrence pairs one to one, while graph and integration id alone repeat
/// (Astra: "Those identifiers repeat").
#[test]
fn occurrences_pair_off_and_on_uniquely() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let keys = |timing: &[String]| -> Result<Vec<(String, String, String)>> {
            timing
                .iter()
                .map(|line| {
                    let get = |key| {
                        field(line, key)
                            .map(str::to_owned)
                            .ok_or_else(|| eyre!("no {key} in {line}"))
                    };
                    Ok((get("kind")?, get("occurrence")?, get("route")?))
                })
                .collect()
        };
        let (off, on) = ledger::with_timing(true, || -> Result<(Run, Run)> {
            Ok((
                run(L3_36, Switch::OFF, None, Some("000036"))?,
                run(
                    L3_36,
                    Switch::new(Mode::On, &ROOT_CHILD),
                    None,
                    Some("000036"),
                )?,
            ))
        })?;
        let _ = off.value?;
        let _ = on.value?;
        ensure!(off.captured.is_empty(), "off builds no record");
        let (off_keys, on_keys) = (keys(&off.timing)?, keys(&on.timing)?);
        let unique = |keys: &[(String, String, String)]| {
            keys.iter()
                .map(|(kind, key, _)| (kind.clone(), key.clone()))
                .collect::<BTreeSet<_>>()
        };
        ensure!(
            unique(&off_keys).len() == off_keys.len() && unique(&on_keys).len() == on_keys.len(),
            "repeated occurrence keys: off {off_keys:#?} on {on_keys:#?}"
        );
        ensure!(
            unique(&off_keys) == unique(&on_keys),
            "the off and on occurrences differ:\noff {off_keys:#?}\non {on_keys:#?}"
        );
        ensure!(off_keys.iter().all(|(_, _, route)| route == "today"));
        let records = on.captured.len();
        let taken = on_keys
            .iter()
            .filter(|(_, _, route)| route == "uncontracted" || route == "zero")
            .count();
        ensure!(
            records > 0
                && taken
                    + on_keys
                        .iter()
                        .filter(|(_, _, r)| r.starts_with("declined:"))
                        .count()
                    == records,
            "{records} records, routes {on_keys:#?}"
        );
        // graph and integration id alone repeat among the child records
        let ids = on
            .timing
            .iter()
            .filter(|line| field(line, "kind") == Some("child"))
            .map(|line| field(line, "id").unwrap_or_default())
            .collect::<Vec<_>>();
        ensure!(
            ids.iter().collect::<BTreeSet<_>>().len() < ids.len(),
            "no repeated child id to disambiguate: {ids:?}"
        );
        // without the timing lines mode off computes no key and prints nothing
        let quiet = run(L3_36, Switch::OFF, None, None)?;
        ensure!(quiet.timing.is_empty() && quiet.lines.is_empty());
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

// ---------------------------------------------------------------------------
// Review of the U7a delta at 6c492c9ddc: the fixes that need new interfaces
// (the base-interface reproducers are `uncontracted_review_tests`)
// ---------------------------------------------------------------------------

/// zero-comparison 2: an enabled kind whose route is off is refused at the
/// switch's initialisation (at 6c492c9ddc `check:child` with
/// `GAMMALOOP_UV_CHILDREN=massive`, and `check:ir-hard` with
/// `GAMMALOOP_IR_HARD=off`, passed with `expected=0 outcome=ok`).
#[test]
fn a_kind_whose_route_is_off_is_refused() {
    let routes = |uv_children_rstar, ir_hard_child| ledger::Routes {
        uv_children_rstar,
        ir_hard_child,
        ..ledger::Routes::PRODUCTION
    };
    for mode in [Mode::Check, Mode::On] {
        let child = Switch::new(mode, &[Kind::Child]);
        let error = ledger::initialise(child, None, None, || routes(false, true), || false, None)
            .unwrap_err();
        assert!(error.contains("GAMMALOOP_UV_CHILDREN=massive"), "{error}");
        let ir_hard = Switch::new(mode, &[Kind::IrHard]);
        let error = ledger::initialise(ir_hard, None, None, || routes(true, false), || false, None)
            .unwrap_err();
        assert!(error.contains("GAMMALOOP_IR_HARD"), "{error}");
        // the kinds whose routes are on, and root (always produced)
        assert!(
            ledger::initialise(child, None, None, || routes(true, false), || false, None).is_ok()
        );
        assert!(
            ledger::initialise(ir_hard, None, None, || routes(false, true), || false, None).is_ok()
        );
        assert!(
            ledger::initialise(
                Switch::new(mode, &[Kind::Root]),
                None,
                None,
                || routes(false, false),
                || false,
                None
            )
            .is_ok()
        );
    }
    // off reads no route
    assert!(
        ledger::initialise(
            Switch::OFF,
            None,
            None,
            || panic!("off reads the routes"),
            || false,
            None
        )
        .is_ok()
    );
}

/// Helper VH's C7, settings exception E2 (c7-normalisation.md sections
/// 7-8): with the switch on or check, every quotient route but `rstar` is
/// refused at the initialisation, for every kind list (each kind's records
/// carry child Δs whose R̄ sees the quotient), and so is a value
/// `GAMMALOOP_IR_QUOTIENT` itself refuses; off reads no route.
/// (`uncontracted_init_refusal_tests` has the environment's half.)
#[test]
fn a_quotient_route_other_than_rstar_is_refused() {
    use super::top_mass::IrQuotient;
    let routes = |ir_quotient| ledger::Routes {
        ir_quotient,
        ..ledger::Routes::PRODUCTION
    };
    for mode in [Mode::Check, Mode::On] {
        for kinds in [
            &[Kind::Root][..],
            &[Kind::Child],
            &[Kind::IrHard],
            &ALL_KINDS,
        ] {
            let switch = Switch::new(mode, kinds);
            for (route, name) in [
                (IrQuotient::Massive, "massive"),
                (IrQuotient::Off, "off"),
                (IrQuotient::Check, "check"),
                (IrQuotient::Auto, "auto"),
            ] {
                let error =
                    ledger::initialise(switch, None, None, || routes(Ok(route)), || false, None)
                        .unwrap_err();
                assert!(
                    error.contains(&format!("GAMMALOOP_IR_QUOTIENT={name} is refused"))
                        && error.contains("E2"),
                    "{error}"
                );
            }
            assert!(
                ledger::initialise(
                    switch,
                    None,
                    None,
                    || routes(Ok(IrQuotient::Rstar)),
                    || false,
                    None
                )
                .is_ok()
            );
            let error = ledger::initialise(
                switch,
                None,
                None,
                || routes(Err("GAMMALOOP_IR_QUOTIENT must be one of ...".into())),
                || false,
                None,
            )
            .unwrap_err();
            assert!(error.contains("GAMMALOOP_IR_QUOTIENT must be"), "{error}");
        }
    }
}

/// Helper VH's C7, settings exception E1: every normalisation but exactly
/// `MSbar` is refused as an error (never designed), and by the start-up
/// check [`super::validate_uncontracted_settings`] in on and check, not off.
/// `additional_normalization` is not a settings refusal: it stays the
/// designed decline of the record (`uncontracted_init_refusal_tests`).
#[test]
fn the_normalisation_must_be_msbar() {
    let mut uv = settings();
    assert!(uncontracted::settings_refusal(&uv.vakint).is_ok());
    for mode in [Mode::Check, Mode::On, Mode::Off] {
        let switch = if mode == Mode::Off {
            Switch::OFF
        } else {
            Switch::new(mode, &ROOT_CHILD)
        };
        assert!(with_switch(switch, || super::validate_uncontracted_settings(&uv)).is_ok());
    }
    for name in [
        "GScheme",
        "FMFTandMATAD",
        "pySecDec",
        "MSbar ",
        "msbar",
        "MS",
    ] {
        uv.vakint.normalization = name.into();
        let refusal = uncontracted::settings_refusal(&uv.vakint).unwrap_err();
        assert!(
            refusal.code == "vakint_normalization" && !refusal.designed,
            "{name}: {refusal:?}"
        );
        for mode in [Mode::Check, Mode::On] {
            let error = with_switch(Switch::new(mode, &ROOT_CHILD), || {
                super::validate_uncontracted_settings(&uv)
            })
            .unwrap_err();
            let text = format!("{error}");
            assert!(
                text.contains("vakint_normalization") && text.contains("E1"),
                "{name}: {text}"
            );
        }
        assert!(with_switch(Switch::OFF, || super::validate_uncontracted_settings(&uv)).is_ok());
    }
    let mut uv = settings();
    uv.vakint.additional_normalization = "-1".into();
    assert!(uncontracted::settings_refusal(&uv.vakint).is_ok());
}

/// Helper VH's C7, settings exception E3: the evaluation methods must lie
/// in {rustred, forcer} (either order, either alone); any other entry is an
/// error (never designed) and refused by the start-up check in on and
/// check, not off.
#[test]
fn the_evaluation_methods_must_be_rustred_or_forcer() {
    let methods = |list: &[&str]| list.iter().map(|m| m.to_string()).collect::<Vec<_>>();
    let mut uv = settings();
    for accepted in [
        &["rustred", "forcer"][..],
        &["forcer", "rustred"],
        &["rustred"],
        &["forcer"],
    ] {
        uv.vakint.evaluation_methods = methods(accepted);
        assert!(
            uncontracted::settings_refusal(&uv.vakint).is_ok(),
            "{accepted:?}"
        );
        assert!(
            with_switch(Switch::new(Mode::Check, &ROOT_CHILD), || {
                super::validate_uncontracted_settings(&uv)
            })
            .is_ok()
        );
    }
    for refused in [
        &["matad"][..],
        &["rustred", "forcer", "fmft"],
        &["alphaloop", "matad", "fmft"],
        &["rustred", "pysecdec"],
        &["rustred", "RustRed"],
    ] {
        uv.vakint.evaluation_methods = methods(refused);
        let refusal = uncontracted::settings_refusal(&uv.vakint).unwrap_err();
        assert!(
            refusal.code == "vakint_evaluation_methods" && !refusal.designed,
            "{refused:?}: {refusal:?}"
        );
        for mode in [Mode::Check, Mode::On] {
            let error = with_switch(Switch::new(mode, &ROOT_CHILD), || {
                super::validate_uncontracted_settings(&uv)
            })
            .unwrap_err();
            let text = format!("{error}");
            assert!(
                text.contains("vakint_evaluation_methods") && text.contains("E3"),
                "{refused:?}: {text}"
            );
        }
        assert!(with_switch(Switch::OFF, || super::validate_uncontracted_settings(&uv)).is_ok());
    }
}

/// values-and-off 2: the timing switch is validated at the initialisation in
/// every mode (at 6c492c9ddc a bad value started and then panicked in every
/// diagram); STRICT is read only with the switch on or check.
#[test]
fn the_timing_switch_is_validated_at_start_up() {
    for switch in [
        Switch::OFF,
        Switch::new(Mode::On, &ROOT_CHILD),
        Switch::new(Mode::Check, &ROOT_CHILD),
    ] {
        let error = init(switch, None, Some("yes"), || false, None).unwrap_err();
        assert!(error.contains(ledger::TIMING_VARIABLE), "{error}");
        for accepted in [None, Some(""), Some("0"), Some("1"), Some(" on ")] {
            assert!(
                init(switch, None, accepted, || false, None).is_ok(),
                "{accepted:?}"
            );
        }
    }
    assert_eq!(ledger::parse_timing(Some("1")), Ok(true));
    assert_eq!(ledger::parse_timing(Some("off")), Ok(false));
}

/// The dispositions of a section's K3 records by the correction's end
/// (refusal-integrity 2 and 4, values-and-off 1; pre-review findings 3 and
/// 4): in mode on a value the route then discarded is the error
/// `ir_hard_route:<head>` and fails the diagram, never `taken`; check mode
/// keeps its comparisons there; an exhausted guard budget declines the
/// correction's records (`ir_hard_guard:budget`, designed); a panic or an
/// error the guard caught makes every record of the correction the error
/// `ir_hard_guard:<panic|error>` and fails the diagram.
#[test]
fn a_declined_correction_reclassifies_its_records() -> Result<()> {
    use super::uncontracted_ledger::{Section, SectionEnd};
    test_initialise()?;
    let ir_hard = [Kind::Root, Kind::Child, Kind::IrHard];
    // (mode, end, completed records, interrupted, expected line fields, close ok)
    let cases: [(Mode, SectionEnd, usize, usize, &[(&str, &str)], bool); 8] = [
        (
            Mode::On,
            SectionEnd::Kept,
            2,
            0,
            &[("declined", "none"), ("outcome", "ok")],
            true,
        ),
        // pre-review finding 4: a route decline after taken values is an
        // error in mode on
        (
            Mode::On,
            SectionEnd::RouteDeclined("moments:not_a_polynomial".into()),
            2,
            0,
            &[
                ("declined", "none"),
                ("error", "ir_hard_route:moments:2"),
                ("outcome", "error"),
            ],
            false,
        ),
        (
            Mode::Check,
            SectionEnd::RouteDeclined("depth:3".into()),
            2,
            0,
            &[("declined", "none"), ("compared", "2"), ("outcome", "ok")],
            true,
        ),
        (
            Mode::On,
            SectionEnd::GuardBudget,
            1,
            1,
            // every record declined: mode on's all-declined rule (pre-review
            // finding 7) fails it unless named
            &[
                ("declined", "ir_hard_guard:budget:2"),
                ("outcome", "uncovered"),
            ],
            false,
        ),
        (
            Mode::Check,
            SectionEnd::GuardBudget,
            1,
            1,
            &[
                ("declined", "ir_hard_guard:budget:1"),
                ("compared", "1"),
                ("outcome", "ok"),
            ],
            true,
        ),
        // check mode: the record's today's path failed on the budget inside
        // the guard (counted unchecked by the hook): declined with the budget
        (
            Mode::Check,
            SectionEnd::GuardBudget,
            0,
            0,
            &[
                ("declined", "ir_hard_guard:budget:1"),
                ("unchecked", "0"),
                ("outcome", "uncovered"),
            ],
            false,
        ),
        // every record of a section the guard failed is the error, inside
        // a record or after them (pre-review finding 3)
        (
            Mode::On,
            SectionEnd::GuardFailure("panic:index_out_of_bounds".into()),
            1,
            1,
            &[("error", "ir_hard_guard:panic:2"), ("outcome", "error")],
            false,
        ),
        (
            Mode::Check,
            SectionEnd::GuardFailure("error:taylor".into()),
            2,
            0,
            &[
                ("error", "ir_hard_guard:error:2"),
                ("compared", "0"),
                ("outcome", "error"),
            ],
            false,
        ),
    ];
    for (mode, end, completed, interrupted, fields, close_ok) in cases {
        with_switch(Switch::new(mode, &ir_hard), || -> Result<()> {
            let _ = ledger::take_lines();
            let diagram = DiagramLedger::open("review_section", "000009")
                .ok_or_else(|| eyre!("no ledger"))?;
            let section = Section::open().ok_or_else(|| eyre!("no section"))?;
            for _ in 0..completed + interrupted {
                ledger::expect(Kind::IrHard);
            }
            for _ in 0..completed {
                ledger::record_value(Kind::IrHard, mode);
            }
            if completed + interrupted == 0 {
                // the hook's report of an attempt whose today's path failed
                ledger::expect(Kind::IrHard);
                ledger::record_unchecked(Kind::IrHard);
            }
            let _ = section.close(end.clone());
            let closed = diagram.close(true);
            let lines = ledger::take_lines();
            let line = disposition(&lines, Kind::IrHard)?;
            for (key, value) in fields {
                ensure!(
                    field(line, key) == Some(*value),
                    "{} {end:?}: {key}={value} expected: {line}",
                    mode.as_str()
                );
            }
            ensure!(
                count(line, "expected") == (completed + interrupted).max(1),
                "{} {end:?}: {line}",
                mode.as_str()
            );
            ensure!(
                closed.is_ok() == close_ok,
                "{} {end:?}: close {closed:?}",
                mode.as_str()
            );
            Ok(())
        })?;
    }
    // the guard's budget is the one designed section code (an exception
    // may name it); the route's declines and the guard's failures are errors
    // (pre-review findings 3 and 4)
    assert!(uncontracted::is_designed_code("ir_hard_guard:budget"));
    for code in [
        "ir_hard_guard:panic",
        "ir_hard_route:",
        "ir_hard_route:a:b",
        "ir_hard_route:moments",
        "ir_hard_route:depth",
        "ir_hard_route:panic",
    ] {
        assert!(!uncontracted::is_designed_code(code), "{code}");
    }
    // vakint's text of an exhausted scoped budget (the guard's budget)
    assert!(uncontracted::is_budget_text(
        "evaluation: vacuum: process RSS 243523584 bytes plus requested allocation 0 bytes exceeds RAM budget 243513052 bytes"
    ));
    assert!(!uncontracted::is_budget_text(
        "every line is massive: not a closure"
    ));
    Ok(())
}

/// refusal-integrity 2: vakint's unknown-master refusal is a diagnostic
/// failure, so the IR hard child guard (`GAMMALOOP_IR_HARD_CHILD_GUARD=1`)
/// propagates it instead of declining the correction to today's terms (at
/// 6c492c9ddc it was a plain error: `Ok(Err("error:..."))`).
#[test]
fn unknown_masters_survive_the_child_guard() {
    let failure = || {
        uncontracted::unknown_masters_failure(
            Kind::IrHard,
            "g",
            1,
            "expected master coefficient(s) that are not known yet",
        )
    };
    assert!(super::ir_dirac::is_diagnostic_failure(&failure()));
    let guarded = super::ir_dirac::guarded_unless_diagnostic(None, None, || {
        Err::<std::result::Result<(), String>, _>(failure())
    });
    assert!(
        guarded.is_err(),
        "the guard declined the unknown masters: {guarded:?}"
    );
}

/// zero-comparison 3: a remainder of a kind the switch does not enable is
/// counted inert and listed on the summary line (`inert=<kind>:<n>`), so a
/// remainder classified into an inert kind is visible.
#[test]
fn inert_remainders_are_listed() -> Result<()> {
    test_initialise()?;
    with_switch(Switch::new(Mode::Check, &ROOT_CHILD), || -> Result<()> {
        let _ = ledger::take_lines();
        let diagram =
            DiagramLedger::open("review_inert", "000004").ok_or_else(|| eyre!("no ledger"))?;
        ledger::inert(Kind::Quotient);
        ledger::inert(Kind::Quotient);
        ledger::inert(Kind::Surrogate);
        let evidence = diagram.close(true)?.ok_or_else(|| eyre!("no evidence"))?;
        let lines = ledger::take_lines();
        let summary = lines_of(&lines, "UNCONTRACTED_SUMMARY ");
        ensure!(
            summary.len() == 1 && field(summary[0], "inert") == Some("quotient:2,surrogate:1"),
            "{summary:#?}"
        );
        ensure!(evidence["inert"]["quotient"] == 2, "{evidence:#}");
        Ok(())
    })
}

/// refusal-integrity 1, the key: an R* child Δ looked up inside a scope of a
/// kind the switch does not enable carries `#uncontracted-inert=<kind>`
/// (on and check); an enabled kind's key is unchanged, and off has no suffix.
#[test]
fn inert_scope_child_keys_are_their_own() {
    let check = Switch::new(Mode::Check, &ROOT_CHILD);
    let on = Switch::new(Mode::On, &ROOT_CHILD);
    assert_eq!(with_switch(check, uncontracted::cache_key_suffix), None);
    assert_eq!(
        with_switch(check, || super::expression::in_dressed_scope(
            uncontracted::cache_key_suffix
        ))
        .as_deref(),
        Some("#uncontracted-inert=surrogate")
    );
    assert_eq!(
        with_switch(on, || super::expression::in_dressed_scope(
            uncontracted::cache_key_suffix
        ))
        .as_deref(),
        Some("#uncontracted=child#uncontracted-inert=surrogate")
    );
    assert_eq!(
        with_switch(Switch::OFF, || super::expression::in_dressed_scope(
            uncontracted::cache_key_suffix
        )),
        None
    );
}

// ---------------------------------------------------------------------------
// Pre-review of the U7a delta at 86de2ec899: the fixes that need new
// interfaces (the 86de2ec899-interface reproducers are
// `uncontracted_prereview_tests`)
// ---------------------------------------------------------------------------

/// Finding 5: a comparison counted `compared` that did not agree is an
/// error outcome and fails the diagram, whether or not the hook raised it.
#[test]
fn a_compared_record_that_did_not_agree_is_an_error() -> Result<()> {
    test_initialise()?;
    with_switch(Switch::new(Mode::Check, &ROOT_CHILD), || -> Result<()> {
        for (equal, outcome) in [(true, "ok"), (false, "error")] {
            let _ = ledger::take_lines();
            let diagram = DiagramLedger::open("prereview_equal", "000005")
                .ok_or_else(|| eyre!("no ledger"))?;
            ledger::expect(Kind::Root);
            ledger::record_compared(Kind::Root, equal);
            ledger::expect(Kind::Child);
            ledger::record_compared(Kind::Child, true);
            let closed = diagram.close(true);
            let lines = ledger::take_lines();
            let root = disposition(&lines, Kind::Root)?;
            ensure!(
                field(root, "outcome") == Some(outcome) && closed.is_ok() == equal,
                "equal={equal}: {root} ({closed:?})"
            );
            if !equal {
                let error = format!("{:#}", closed.unwrap_err());
                ensure!(error.contains("did not agree"), "{error}");
            }
        }
        Ok(())
    })
}

/// Finding 6: the run-level exception row is `*<TAB><kind><TAB>none_compared`
/// only (the graph `*` takes that code, the code takes graph `*`).
#[test]
fn run_level_exception_rows() {
    let row =
        ledger::parse_exceptions("*\tir-hard\tnone_compared\tno IR correction of order >= 1\n")
            .unwrap();
    assert_eq!(row.entries[0].graph, "*");
    for bad in [
        "*\troot\tone_loop_tree_cut\twhy\n",
        "g\troot\tnone_compared\twhy\n",
    ] {
        assert!(ledger::parse_exceptions(bad).is_err(), "{bad:?}");
    }
}

/// Finding 6: the run tally. An enabled kind with nothing compared over the
/// run fails it unless the run-level row names it; a nonzero diagram with
/// nothing compared in any kind fails it; a record without the run's
/// evidence fails it; the totals are printed per kind.
#[test]
fn the_run_fails_a_kind_with_nothing_compared() -> Result<()> {
    test_initialise()?;
    let check = Switch::new(Mode::Check, &ALL_KINDS);
    let evidence = |identity: Option<String>, ir_hard: serde_json::Value| {
        serde_json::json!({
            "schema": ledger::EVIDENCE_SCHEMA, "identity": identity, "graph": "g",
            "diagram": "000001", "outcome": "ok", "inert": {},
            "kinds": {
                "root": {"expected": 2, "attempted": 2, "uncontracted": 2, "zero": 0, "equal": 2,
                         "classes": {"compared": 2}, "outcome": "ok"},
                "child": {"expected": 3, "attempted": 1, "uncontracted": 1, "zero": 0, "equal": 1,
                          "hit_equal": 2, "classes": {"compared": 1, "hit:compared": 2}, "outcome": "ok"},
                "ir-hard": ir_hard,
            },
        })
    };
    let none = serde_json::json!({"expected": 0, "attempted": 0, "uncontracted": 0, "zero": 0,
                                  "equal": 0, "classes": {}, "outcome": "ok"});
    let compared = serde_json::json!({"expected": 1, "attempted": 1, "uncontracted": 1, "zero": 0,
                                      "equal": 1, "classes": {"compared": 1}, "outcome": "ok"});
    with_switch(check, || -> Result<()> {
        let identity = ledger::current_identity();
        // nothing ir-hard over the run: a failure naming the kind
        let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add(
            "g/000001",
            Some(&evidence(identity.clone(), none.clone())),
            true,
        );
        let failures = tally.failures();
        ensure!(
            failures.len() == 1 && failures[0].contains("kind ir-hard: no record compared"),
            "{failures:#?}"
        );
        ensure!(
            tally
                .fields()
                .contains("child:expected=3,compared=1,hit_compared=2,"),
            "{}",
            tally.fields()
        );
        // one compared ir-hard record anywhere in the run passes it
        tally.add(
            "g/000002",
            Some(&evidence(identity.clone(), compared.clone())),
            true,
        );
        ensure!(tally.failures().is_empty(), "{:#?}", tally.failures());
        // a run-level row accepts the kind
        let row = ledger::parse_exceptions("*\tir-hard\tnone_compared\tthe test\n")
            .map_err(|error| eyre!(error))?;
        ledger::with_exceptions(row, || -> Result<()> {
            let identity = ledger::current_identity();
            let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
            tally.add("g/000001", Some(&evidence(identity, none.clone())), true);
            ensure!(tally.failures().is_empty(), "{:#?}", tally.failures());
            Ok(())
        })?;
        // a nonzero diagram with nothing compared in any kind; a record
        // without evidence; evidence of another identity
        let empty = serde_json::json!({
            "schema": ledger::EVIDENCE_SCHEMA, "identity": identity, "graph": "z",
            "diagram": "000009", "outcome": "ok", "inert": {},
            "kinds": {"root": none, "child": none, "ir-hard": none},
        });
        let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add(
            "g/000002",
            Some(&evidence(identity.clone(), compared.clone())),
            true,
        );
        tally.add("z/000009", Some(&empty), false);
        ensure!(
            tally.failures().is_empty(),
            "a zero diagram: {:#?}",
            tally.failures()
        );
        tally.add("z/000009", Some(&empty), true);
        tally.add("y/000010", None, true);
        tally.add(
            "x/000011",
            Some(&evidence(
                Some("on:root,child,ir-hard:strict".into()),
                compared,
            )),
            true,
        );
        let failures = tally.failures();
        // with a partial kind list the diagram rule does not apply (the
        // kind-level rule does)
        let partial = with_switch(Switch::new(Mode::Check, &[Kind::IrHard]), || {
            let identity = ledger::current_identity();
            let mut tally = ledger::RunTally::new().expect("a tally");
            let empty = serde_json::json!({
                "schema": ledger::EVIDENCE_SCHEMA, "identity": identity, "graph": "z",
                "diagram": "000009", "outcome": "ok", "inert": {},
                "kinds": {"ir-hard": none},
            });
            tally.add("z/000009", Some(&empty), true);
            tally.failures()
        });
        ensure!(
            partial.len() == 1 && partial[0].contains("kind ir-hard: no record compared"),
            "{partial:#?}"
        );
        ensure!(
            failures.len() == 3
                && failures[0].contains("z/000009: the counterterm is nonzero")
                && failures[1].contains("y/000010: a record without the run's evidence")
                && failures[2].contains("x/000011: a record without the run's evidence"),
            "{failures:#?}"
        );
        Ok(())
    })
}

/// Item B (helper VH 2's `k5_beta_ledger_misses_settings_wide_root_declines`,
/// without the K5 code): one kind with records of two classes (frame shapes,
/// as a K5 surrogate's K1 roots and K2 children). Every K1 record declines
/// by a setting's code (`open_root`) while the K2 records compare: the
/// kind-wide all-declined rule passes it (compared > 0), the class rule
/// must fail the diagram, in check and in on. Controls: a class that
/// declines by a graph code (`one_loop_tree_cut`) or that has one record
/// compared passes; an exception naming the code accepts it (`excepted`).
#[test]
fn a_setting_wide_decline_of_one_record_class_fails() -> Result<()> {
    test_initialise()?;
    let named = |code: &str| Exceptions {
        entries: vec![Exception {
            graph: "setting_wide".to_owned(),
            kind: Kind::Root,
            code: code.to_owned(),
            reason: "the test's setting".to_owned(),
        }],
        digest: Some("test".to_owned()),
    };
    for mode in [Mode::Check, Mode::On] {
        with_switch(Switch::new(mode, &ROOT_CHILD), || -> Result<()> {
            // (K1 declines, K1 evaluated, exceptions) -> the diagram passes
            let cases: [(&[&str], usize, Option<Exceptions>, Option<&str>); 5] = [
                (&["open_root"; 3], 0, None, None),
                (
                    &["open_root", "open_root", "one_loop_tree_cut"],
                    0,
                    None,
                    None,
                ),
                (&["one_loop_tree_cut"; 3], 0, None, Some("ok")),
                (&["open_root"; 2], 1, None, Some("ok")),
                (
                    &["open_root"; 3],
                    0,
                    Some(named("open_root")),
                    Some("excepted"),
                ),
            ];
            for (declines, evaluated, exceptions, passes) in cases {
                let what = format!("{} {declines:?} + {evaluated} evaluated", mode.as_str());
                let run = || -> Result<(Result<Option<serde_json::Value>>, Vec<String>)> {
                    let _ = ledger::take_lines();
                    let diagram = DiagramLedger::open("setting_wide", "000007")
                        .ok_or_else(|| eyre!("no ledger"))?;
                    for code in declines {
                        ledger::record_shaped(Kind::Root, "K1", Some(code), mode);
                    }
                    for _ in 0..evaluated {
                        ledger::record_shaped(Kind::Root, "K1", None, mode);
                    }
                    for _ in 0..2 {
                        ledger::record_shaped(Kind::Root, "K2", None, mode);
                    }
                    ledger::record_shaped(Kind::Child, "K2", None, mode);
                    let closed = diagram.close(true);
                    Ok((closed, ledger::take_lines()))
                };
                let (closed, lines) = match exceptions {
                    Some(exceptions) => ledger::with_exceptions(exceptions, run)?,
                    None => run()?,
                };
                let root = disposition(&lines, Kind::Root)?;
                match passes {
                    Some(outcome) => ensure!(
                        closed.is_ok() && field(root, "outcome") == Some(outcome),
                        "{what}: {root} ({closed:?})"
                    ),
                    None => {
                        let error = format!(
                            "{:#}",
                            closed
                                .err()
                                .ok_or_else(|| eyre!("{what}: the diagram passed: {root}"))?
                        );
                        ensure!(
                            error.contains("K1-shaped records declined")
                                && error.contains("setting-caused")
                                && field(root, "outcome") == Some("uncovered"),
                            "{what}: {error} / {root}"
                        );
                    }
                }
            }
            Ok(())
        })?;
    }
    Ok(())
}

/// Item B, the run tally and the resume: saved evidence whose K1 class all
/// declined by a setting's code (K2 compared) fails the run tally and is
/// not reused (recomputed); one K1 record compared, or the kind excepted,
/// passes.
#[test]
fn the_run_and_the_resume_fail_a_setting_wide_class() -> Result<()> {
    test_initialise()?;
    with_switch(Switch::new(Mode::Check, &ROOT_CHILD), || -> Result<()> {
        let identity = ledger::current_identity();
        let evidence = |k1: serde_json::Value, outcome: &str| {
            let k1_records = k1.as_object().map_or(0, |classes| {
                classes.values().filter_map(|n| n.as_u64()).sum::<u64>()
            });
            let mut classes = serde_json::Map::new();
            for (class, n) in k1.as_object().into_iter().flatten() {
                classes.insert(class.clone(), n.clone());
            }
            let k1_compared = k1["compared"].as_u64().unwrap_or(0);
            classes.insert("compared".to_owned(), serde_json::json!(2 + k1_compared));
            serde_json::json!({
                "schema": ledger::EVIDENCE_SCHEMA, "identity": identity, "graph": "g",
                "diagram": "000001", "outcome": "ok", "inert": {},
                "kinds": {
                    "root": {"expected": 2 + k1_records, "attempted": 2 + k1_records,
                             "uncontracted": 2 + k1_compared, "zero": 0, "equal": 2 + k1_compared,
                             "classes": classes, "shapes": {"K1": k1, "K2": {"compared": 2}},
                             "outcome": outcome},
                    "child": {"expected": 1, "attempted": 1, "uncontracted": 1, "zero": 0,
                              "equal": 1, "classes": {"compared": 1},
                              "shapes": {"K2": {"compared": 1}}, "outcome": "ok"},
                },
            })
        };
        let setting = evidence(serde_json::json!({"declined:open_root": 3}), "ok");
        let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add("g/000001", Some(&setting), true);
        let failures = tally.failures();
        ensure!(
            failures.len() == 1
                && failures[0].contains("kind root: all 3 K1-shaped records over the run declined")
                && failures[0].contains("open_root:3"),
            "{failures:#?}"
        );
        match ledger::resumed("g", "000001", Some(&setting)) {
            Resumed::Unchecked(reason) if reason == "verdict:root:uncovered" => {}
            other => return Err(eyre!("the setting-wide evidence was reused: {other:?}")),
        }
        // one K1 record compared, anywhere in the run, passes; so does an
        // excepted kind
        let mixed = evidence(
            serde_json::json!({"declined:open_root": 2, "compared": 1}),
            "ok",
        );
        ensure!(
            ledger::resumed("g", "000001", Some(&mixed)) == Resumed::Checked,
            "the mixed evidence"
        );
        let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add("g/000001", Some(&mixed), true);
        ensure!(tally.failures().is_empty(), "{:#?}", tally.failures());
        let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add(
            "g/000001",
            Some(&evidence(
                serde_json::json!({"declined:open_root": 3}),
                "excepted",
            )),
            true,
        );
        ensure!(tally.failures().is_empty(), "{:#?}", tally.failures());
        Ok(())
    })
}

/// Review of the follow-ups (ledger lens, finding 1): evidence written
/// before the record classes existed (`u7a-disposition.v2`, no `shapes`:
/// e.g. helper VH 2's K5 binaries, 16 K1 roots `declined:open_root` and 35
/// K2 children compared in one kind) deserialised with an empty shape map,
/// so the class rule saw nothing and the resume and the run tally accepted
/// it (an import across pipeline fingerprints reaches both). The schema is
/// bumped: such evidence is recomputed (`reason=schema:`) and fails the
/// tally.
#[test]
fn evidence_without_record_classes_is_not_reused() -> Result<()> {
    test_initialise()?;
    ensure!(
        ledger::EVIDENCE_SCHEMA != "u7a-disposition.v2",
        "the schema that wrote no record classes is still current"
    );
    with_switch(Switch::new(Mode::Check, &ROOT_CHILD), || -> Result<()> {
        let identity = ledger::current_identity();
        let old = serde_json::json!({
            "schema": "u7a-disposition.v2", "identity": identity, "graph": "g",
            "diagram": "000001", "outcome": "ok", "inert": {},
            "kinds": {
                "root": {"expected": 51, "attempted": 51, "uncontracted": 35, "zero": 0,
                         "equal": 35, "classes": {"compared": 35, "declined:open_root": 16},
                         "outcome": "ok"},
                "child": {"expected": 1, "attempted": 1, "uncontracted": 1, "zero": 0,
                          "equal": 1, "classes": {"compared": 1}, "outcome": "ok"},
            },
        });
        match ledger::resumed("g", "000001", Some(&old)) {
            Resumed::Unchecked(reason) if reason.starts_with("schema:") => {}
            other => {
                return Err(eyre!(
                    "evidence without record classes was reused: {other:?}"
                ));
            }
        }
        let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add("g/000001", Some(&old), true);
        let failures = tally.failures();
        ensure!(
            failures
                .iter()
                .any(|failure| failure.contains("g/000001: a record without the run's evidence")),
            "{failures:#?}"
        );
        // the same counts under the current schema, with their classes, fail
        // by the class rule (not reused, and the tally fails the class)
        let mut current = old.clone();
        current["schema"] = serde_json::json!(ledger::EVIDENCE_SCHEMA);
        current["kinds"]["root"]["shapes"] =
            serde_json::json!({"K1": {"declined:open_root": 16}, "K2": {"compared": 35}});
        current["kinds"]["child"]["shapes"] = serde_json::json!({"K2": {"compared": 1}});
        ensure!(
            matches!(
                ledger::resumed("g", "000001", Some(&current)),
                Resumed::Unchecked(_)
            ),
            "the current evidence with a setting-wide class was reused"
        );
        Ok(())
    })
}

/// Review of the follow-ups (ledger lens, finding 2): a [`Section`]'s
/// reclassification moves the records in the record classes too. Mode on,
/// one kind with two classes: K1 taken outside the section (the kind-wide
/// rule sees a taken record), and inside an IR correction whose guard
/// budget ran out one K2 record declined by a setting's code (`open_root`)
/// and one K2 record taken. The budget makes the taken K2 record
/// `declined:ir_hard_guard:budget`; with the class map left at `taken`
/// the K2 class looked evidenced and the diagram passed. Now every K2
/// record declined, one by a setting's code: the class rule fails it. The
/// saved evidence carries the moved classes. Control: the same with the
/// K2 decline a graph code passes.
#[test]
fn a_section_reclassifies_the_record_classes_too() -> Result<()> {
    use super::uncontracted_ledger::{Section, SectionEnd};
    test_initialise()?;
    with_switch(Switch::new(Mode::On, &ROOT_CHILD), || -> Result<()> {
        for (code, fails) in [("open_root", true), ("one_loop_tree_cut", false)] {
            let _ = ledger::take_lines();
            let diagram = DiagramLedger::open("section_shapes", "000005")
                .ok_or_else(|| eyre!("no ledger"))?;
            ledger::record_shaped(Kind::Root, "K1", None, Mode::On);
            ledger::record_shaped(Kind::Child, "K2", None, Mode::On);
            let section = Section::open().ok_or_else(|| eyre!("no section"))?;
            ledger::record_shaped(Kind::Root, "K2", Some(code), Mode::On);
            ledger::record_shaped(Kind::Root, "K2", None, Mode::On);
            section.close(SectionEnd::GuardBudget)?;
            let closed = diagram.close(true);
            let lines = ledger::take_lines();
            let root = disposition(&lines, Kind::Root)?;
            if fails {
                let error = format!(
                    "{:#}",
                    closed
                        .err()
                        .ok_or_else(|| eyre!("{code}: the diagram passed: {root}"))?
                );
                ensure!(
                    error.contains("K2-shaped records declined")
                        && field(root, "outcome") == Some("uncovered"),
                    "{code}: {error} / {root}"
                );
            } else {
                let evidence = closed?.ok_or_else(|| eyre!("no evidence"))?;
                let k2 = &evidence["kinds"]["root"]["shapes"]["K2"];
                ensure!(
                    field(root, "outcome") == Some("ok")
                        && k2["declined:ir_hard_guard:budget"] == 1
                        && k2["taken"].is_null(),
                    "{code}: {root} / {k2}"
                );
            }
        }
        Ok(())
    })
}

/// Review of the follow-ups (ledger lens, finding 3): with
/// `GAMMALOOP_IR_DIRAC_CLOSURE=off` the hard part's open quark line stays
/// open, so every quark-line K3 record declines the designed
/// `vakint:open_spinor` because of a setting while gluon-line K3 records
/// compare (one K3 class: neither the kind-wide nor the class rule sees
/// it). With `ir-hard` enabled in on or check it is refused at the
/// initialisation; `on` and `check` (production: on) start, a kind list
/// without `ir-hard` does not read it, and an invalid value is refused.
#[test]
fn the_ir_dirac_closure_off_is_refused_with_ir_hard() {
    use super::ir_dirac::IrDiracMode;
    let routes = |ir_dirac| ledger::Routes {
        ir_dirac,
        ..ledger::Routes::PRODUCTION
    };
    for mode in [Mode::Check, Mode::On] {
        for kinds in [&[Kind::IrHard][..], &ALL_KINDS] {
            let switch = Switch::new(mode, kinds);
            let error = ledger::initialise(
                switch,
                None,
                None,
                || routes(Ok(IrDiracMode::Off)),
                || false,
                None,
            )
            .unwrap_err();
            assert!(
                error.contains("GAMMALOOP_IR_DIRAC_CLOSURE=off is refused")
                    && error.contains("open_spinor"),
                "{error}"
            );
            for accepted in [IrDiracMode::On, IrDiracMode::Check] {
                assert!(
                    ledger::initialise(switch, None, None, || routes(Ok(accepted)), || false, None)
                        .is_ok(),
                    "{accepted:?}"
                );
            }
            let error = ledger::initialise(
                switch,
                None,
                None,
                || routes(Err("GAMMALOOP_IR_DIRAC_CLOSURE must be unset, ...".into())),
                || false,
                None,
            )
            .unwrap_err();
            assert!(
                error.contains("GAMMALOOP_IR_DIRAC_CLOSURE must be"),
                "{error}"
            );
        }
        assert!(
            ledger::initialise(
                Switch::new(mode, &ROOT_CHILD),
                None,
                None,
                || routes(Ok(IrDiracMode::Off)),
                || false,
                None
            )
            .is_ok()
        );
    }
}
