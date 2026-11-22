//! Pick 8 (`crown_legs_at_node`) at a POSITIVE Taylor degree: helper VH's
//! caution on the `100000 + e` leg of a paired crown edge.
//!
//! The record is the K2 child record of L5 006961 id 5 that the gate binary
//! of 0484ddd9f6 wrote (`testdata/u7a/l5-006961-id5-paired-leg.record.txt`):
//! crown edge 13 has both ends in the child, its legs are 13 (out of vertex
//! 4) and 100013 (into the `ct:` vertex 3), and the insertion factor at
//! vertex 3 names `Q(100013, ..)` (the rename of pick 8). Its Taylor degree
//! is 0, as for every record with a paired crown edge (a child or IR hard
//! part with one has at least four legs, `local_dod <= 0`); the variants
//! below raise it to 1 and 2 to test what a positive degree would do:
//!
//! * the leg is scaled exactly as its partner: the record whose
//!   `scaled_edges` also lists 100013 or lists 100013 instead of 13 has the
//!   same value (vakint scales the momentum `p(j)` that a leg carries,
//!   `scaled_edges` only has to cover every `p(j)`); the record with the two
//!   legs' roles swapped (13 into vertex 3, 100013 out of vertex 4, every
//!   `Q(13)` <-> `Q(100013)`) is REFUSED by vakint's paired-line rule (edge
//!   `e` leaves the record, edge `100000 + e` enters it: the emitter's
//!   form), at every degree;
//! * an undeclared `Q(100000 + e)` (the leg's line dropped, or the factor
//!   naming an id without a line) is REFUSED, at degree 0 and above: every
//!   record the emitter writes has a `taylor` line, and vakint's closure
//!   validation refuses a constant external vector with a Taylor step;
//! * the pre-pick-8 form `Q(13)` in the insertion at vertex 3 is refused
//!   (`not a boundary momentum of the insertion`).
//! * a leg that does not mirror its partner's momentum is REFUSED by vakint's
//!   paired-line rule (the two legs of a paired line carry one momentum); the
//!   emitter's guard (`declared_momenta`) refuses it on the record value too.

use eyre::{Result, ensure, eyre};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
use vakint::{
    Vakint,
    p_integral::factorized::{FactorizedKind, FactorizedNumerator, evaluate_factorized_closure},
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    initialisation::test_initialise,
    utils::{GS, W_},
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    RStarPlan,
    boundary::{BoundaryMode, with_boundary_mode},
    plan::{uv_child_edge_ids, with_child_taylor_order},
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{self, Captured, Kind, Mode, Switch, with_capture, with_strict, with_switch},
    uncontracted_oracle_tests::atoms_equal,
};

const RECORD: &str = include_str!("testdata/u7a/l5-006961-id5-paired-leg.record.txt");

fn lines_map(text: &str, f: impl Fn(&str) -> Option<String>) -> String {
    let mut out = String::new();
    for line in text.lines() {
        if let Some(line) = f(line) {
            out.push_str(&line);
            out.push('\n');
        }
    }
    out
}

fn with_degree(text: &str, degree: u32) -> String {
    let out = text.replace("taylor\tdegree=0\t", &format!("taylor\tdegree={degree}\t"));
    assert!(out.contains(&format!("taylor\tdegree={degree}\t")));
    out
}

fn with_scaled(text: &str, scaled: &str) -> String {
    let out = text.replace(
        "\tscaled_edges=0,1,13\t",
        &format!("\tscaled_edges={scaled}\t"),
    );
    assert!(out.contains(&format!("\tscaled_edges={scaled}\t")));
    out
}

/// The legs' roles swapped: 13 enters vertex 3, 100013 leaves vertex 4, every
/// momentum renamed with them; `external 3 edge=13` now names the incoming leg.
fn swapped(text: &str) -> String {
    let out = text
        .replace("::Q(13,", "::Q(@@,")
        .replace("::Q(100013,", "::Q(13,")
        .replace("::Q(@@,", "::Q(100013,")
        .replace(
            "edge\t13\tsource=4\ttarget=ext\thedges=24,-\t",
            "edge\t13\tsource=ext\ttarget=3\thedges=-,25\t",
        )
        .replace(
            "edge\t100013\tsource=ext\ttarget=3\thedges=-,25\t",
            "edge\t100013\tsource=4\ttarget=ext\thedges=24,-\t",
        );
    assert!(out.contains("edge\t13\tsource=ext\ttarget=3\thedges=-,25\t"));
    assert!(out.contains("edge\t100013\tsource=4\ttarget=ext\thedges=24,-\t"));
    out
}

/// The insertion factor's `Q(100013, ..)` renamed to `Q(id, ..)`.
fn insertion_names(text: &str, id: usize) -> String {
    lines_map(text, |line| {
        Some(if line.starts_with("insertion_factor\t") {
            line.replace("::Q(100013,", &format!("::Q({id},"))
        } else {
            line.to_owned()
        })
    })
}

/// The leg's line dropped (the record keeps `Q(100013)` in its insertion).
fn leg_dropped(text: &str) -> String {
    lines_map(text, |line| {
        if line.starts_with("edge\t100013\t") || line.starts_with("edge_factor\t100013\t") {
            None
        } else if line == "edges\t10" {
            Some("edges\t9".to_owned())
        } else if line == "lines\t61" {
            // the record's own line count, so that the refusal is not the count's
            Some("lines\t59".to_owned())
        } else {
            Some(line.to_owned())
        }
    })
}

fn read(text: &str) -> std::result::Result<FactorizedNumerator, String> {
    let value = FactorizedNumerator::from_text(text).map_err(|e| e.to_string())?;
    value.validate().map_err(|e| e.to_string())?;
    Ok(value)
}

/// vakint's raw closure value of a record text (`Err`: the reader's, the
/// validation's or the evaluation's error text).
fn evaluate(vakint: &Vakint, text: &str) -> std::result::Result<Atom, String> {
    let value = read(text)?;
    let FactorizedKind::Closure(spec) = &value.kind else {
        return Err("not a closure record".into());
    };
    let mut vakint_settings = settings().vakint.true_settings();
    vakint_settings.number_of_terms_in_epsilon_expansion =
        i64::try_from(spec.depth.expect("depth")).expect("depth");
    let prefactor = value.prefactor.total();
    let evaluated = evaluate_factorized_closure(
        vakint,
        &vakint_settings,
        std::slice::from_ref(&value),
        std::slice::from_ref(&prefactor),
    )
    .map_err(|e| e.to_string())?;
    evaluated
        .records
        .into_iter()
        .next()
        .map(|record| record.value)
        .ok_or_else(|| "no record value".into())
}

fn refused(what: &str, result: std::result::Result<Atom, String>, needle: &str) -> Result<()> {
    match result {
        Ok(value) => Err(eyre!("{what}: evaluated, expected a refusal: {value:.300}")),
        Err(error) => {
            eprintln!("PICK8 {what}: refused: {error:.400}");
            ensure!(
                error.contains(needle),
                "{what}: refused, but not by `{needle}`: {error}"
            );
            Ok(())
        }
    }
}

#[test]
fn pick8_paired_leg_scales_like_its_partner_at_a_positive_degree() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    ensure!(RECORD.contains("taylor\tdegree=0\tcompensate_mass=1\tscaled_edges=0,1,13\t"));
    ensure!(RECORD.contains("edge\t100013\tsource=ext\ttarget=3\thedges=-,25\t"));
    ensure!(RECORD.contains("insertion\t3\thedges=5,18,25\t"));
    // degree 0 (the record as written); the legs swapped are refused by
    // vakint's paired-line rule (not the emitter's form)
    let base = evaluate(&vakint, RECORD).map_err(|e| eyre!("degree 0: {e}"))?;
    eprintln!(
        "PICK8 degree 0: zero={} len={}",
        base.is_zero(),
        base.to_string().len()
    );
    refused(
        "degree 0 legs swapped",
        evaluate(&vakint, &swapped(RECORD)),
        "the two legs of a paired line carry one momentum",
    )?;
    let mut nonzero = false;
    for degree in [1, 2] {
        let text = with_degree(RECORD, degree);
        let value = evaluate(&vakint, &text).map_err(|e| eyre!("degree {degree}: {e}"))?;
        let zero = value.is_zero();
        nonzero |= !zero;
        let differs_from_0 = !atoms_equal(&value, &base);
        eprintln!(
            "PICK8 degree {degree}: zero={zero} len={} differs_from_degree_0={differs_from_0}",
            value.to_string().len()
        );
        refused(
            &format!("degree {degree} legs swapped"),
            evaluate(&vakint, &swapped(&text)),
            "the two legs of a paired line carry one momentum",
        )?;
        for (what, variant) in [
            (
                "scaled_edges=0,1,13,100013",
                with_scaled(&text, "0,1,13,100013"),
            ),
            ("scaled_edges=0,1,100013", with_scaled(&text, "0,1,100013")),
        ] {
            let other =
                evaluate(&vakint, &variant).map_err(|e| eyre!("degree {degree} {what}: {e}"))?;
            let equal = atoms_equal(&value, &other);
            eprintln!("PICK8 degree {degree} {what}: equal={equal}");
            ensure!(equal, "degree {degree}: {what} differs");
        }
    }
    ensure!(
        nonzero,
        "every positive-degree value is zero: the test shows nothing"
    );
    // an undeclared Q(100000 + e): refused at every degree
    for degree in [0, 1, 2] {
        let text = with_degree(RECORD, degree);
        refused(
            &format!("degree {degree} factor names Q(200013) (no line)"),
            evaluate(&vakint, &insertion_names(&text, 200_013)),
            "constant external vectors",
        )?;
        refused(
            &format!("degree {degree} leg line 100013 dropped"),
            evaluate(&vakint, &leg_dropped(&text)),
            "",
        )?;
        refused(
            &format!("degree {degree} pre-pick-8 Q(13) at vertex 3"),
            evaluate(&vakint, &insertion_names(&text, 13)),
            "not a boundary momentum of the insertion",
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// The emitter's guard (`uncontracted::declared_momenta`), refutation round of
// pick 8: every momentum a record names is declared, an insertion names only
// its boundary momenta, a `100000 + e` leg mirrors its partner.
// ---------------------------------------------------------------------------

fn guard(text: &str) -> std::result::Result<(), String> {
    let value = FactorizedNumerator::from_text(text).map_err(|e| format!("reader: {e}"))?;
    uncontracted::declared_momenta(&value).map_err(|d| format!("{}: {}", d.code, d.detail))
}

/// The leg 100013 given the momentum p(2) instead of its partner's p(3).
fn leg_momentum_changed(text: &str) -> String {
    let from = "edge\t100013\tsource=ext\ttarget=3\thedges=-,25\tparticle=Q\tline=vector\tpower=1\tmomentum=0,0,0,0,1\t";
    let to = "edge\t100013\tsource=ext\ttarget=3\thedges=-,25\tparticle=Q\tline=vector\tpower=1\tmomentum=0,0,0,1,0\t";
    assert!(text.contains(from));
    text.replace(from, to)
}

fn guard_value(value: &FactorizedNumerator) -> std::result::Result<(), String> {
    uncontracted::declared_momenta(value).map_err(|d| format!("{}: {}", d.code, d.detail))
}

fn value_refused(what: &str, value: &FactorizedNumerator, code: &str) -> Result<()> {
    match guard_value(value) {
        Ok(()) => Err(eyre!("{what}: the guard passed it")),
        Err(error) => {
            eprintln!("PICK8 GUARD {what}: {error:.300}");
            ensure!(error.starts_with(&format!("{code}: ")), "{what}: {error}");
            Ok(())
        }
    }
}

/// What vakint's READER says of a text form (`from_text` validates).
fn reader_says(what: &str, text: &str) {
    match FactorizedNumerator::from_text(text) {
        Ok(_) => eprintln!("PICK8 READER {what}: read"),
        Err(e) => eprintln!("PICK8 READER {what}: refused: {e:.200}"),
    }
}

fn renamed_in_insertions(
    value: &FactorizedNumerator,
    from: usize,
    to: usize,
) -> FactorizedNumerator {
    let mut value = value.clone();
    for insertion in &mut value.insertions {
        insertion.factor = insertion
            .factor
            .replace(function!(GS.emr_mom, from, W_.x___))
            .with(function!(GS.emr_mom, to, W_.x___));
    }
    value
}

/// The guard on record VALUES mutated after reading (vakint's reader refuses
/// most of these forms as text already: printed as `PICK8 READER`).
#[test]
fn pick8_guard_refuses_undeclared_and_unmirrored_legs() -> Result<()> {
    test_initialise()?;
    for degree in [0, 2] {
        let text = with_degree(RECORD, degree);
        for (what, form) in [
            ("legs swapped", swapped(&text)),
            ("Q(200013)", insertion_names(&text, 200_013)),
            ("leg line dropped", leg_dropped(&text)),
            ("pre-pick-8 Q(13) at vertex 3", insertion_names(&text, 13)),
            ("leg momentum p(2)", leg_momentum_changed(&text)),
            ("partner not scaled", with_scaled(&text, "0,1,100013")),
        ] {
            reader_says(&format!("degree {degree} {what}"), &form);
        }
        let good =
            FactorizedNumerator::from_text(&text).map_err(|e| eyre!("degree {degree}: {e}"))?;
        guard_value(&good).map_err(|e| eyre!("degree {degree}: {e}"))?;
        // vakint's reader refuses the swapped form (its paired-line rule); the
        // guard is checked on the unvalidated value: it is not the emitter's
        // form either (100000 + e is always the INCOMING leg)
        let read = FactorizedNumerator::from_text(&swapped(&text));
        ensure!(
            read.as_ref().is_err_and(|e| e
                .to_string()
                .contains("the two legs of a paired line carry one momentum")),
            "degree {degree} swapped: vakint's reader did not refuse it: {:?}",
            read.map(|_| ())
        );
        let swapped_value = FactorizedNumerator::parse_text(&swapped(&text))
            .map_err(|e| eyre!("degree {degree} swapped: {e}"))?;
        value_refused(
            &format!("degree {degree} legs swapped"),
            &swapped_value,
            "paired_leg",
        )?;
        value_refused(
            &format!("degree {degree} Q(200013)"),
            &renamed_in_insertions(&good, 100_013, 200_013),
            "undeclared_momentum",
        )?;
        let mut dropped = good.clone();
        dropped.edges.retain(|edge| edge.id != 100_013);
        value_refused(
            &format!("degree {degree} leg line dropped"),
            &dropped,
            "undeclared_momentum",
        )?;
        value_refused(
            &format!("degree {degree} pre-pick-8 Q(13) at vertex 3"),
            &renamed_in_insertions(&good, 100_013, 13),
            "insertion_momentum",
        )?;
        let mut moved = good.clone();
        for edge in moved.edges.iter_mut().filter(|edge| edge.id == 100_013) {
            edge.momentum = vec![0, 0, 0, 1, 0];
        }
        value_refused(
            &format!("degree {degree} leg momentum p(2)"),
            &moved,
            "paired_leg",
        )?;
        let mut unscaled = good.clone();
        if let FactorizedKind::Closure(spec) = &mut unscaled.kind
            && let Some(taylor) = &mut spec.taylor
        {
            taylor.scaled_edges = vec![0, 1, 100_013];
        }
        value_refused(
            &format!("degree {degree} partner not scaled"),
            &unscaled,
            "paired_leg",
        )?;
    }
    Ok(())
}

/// The negative control the guard was written for: a leg that does NOT
/// mirror its partner (100013 carrying p(2), not p(3)). vakint's paired-line
/// rule refuses it at every degree (the two legs of a paired line carry one
/// momentum), so it can no longer change a value silently; the emitter's
/// guard refuses the same value (`pick8_guard_refuses_undeclared_and_unmirrored_legs`).
#[test]
fn pick8_an_unmirrored_leg_is_refused_by_vakint() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    for degree in [0, 2] {
        let text = with_degree(RECORD, degree);
        evaluate(&vakint, &text).map_err(|e| eyre!("degree {degree}: {e}"))?;
        refused(
            &format!("degree {degree} leg momentum p(2)"),
            evaluate(&vakint, &leg_momentum_changed(&text)),
            "the two legs of a paired line carry one momentum",
        )?;
    }
    Ok(())
}

/// The guard over dumped records: `PICK8_GUARD_LIST` names a file of record
/// paths (one per line); every v3 record must parse and pass.
#[test]
#[ignore = "reads record dumps outside the tree: PICK8_GUARD_LIST"]
fn pick8_guard_passes_every_dumped_record() -> Result<()> {
    test_initialise()?;
    let list = std::env::var("PICK8_GUARD_LIST").map_err(|_| eyre!("PICK8_GUARD_LIST unset"))?;
    let (mut v3, mut other, mut legs, mut failed) = (0usize, 0usize, 0usize, Vec::new());
    for path in std::fs::read_to_string(&list)?
        .lines()
        .filter(|l| !l.trim().is_empty())
    {
        let text = std::fs::read_to_string(path.trim())?;
        if !text.starts_with("schema\tpint6.factorized-numerator.v3") {
            other += 1;
            continue;
        }
        v3 += 1;
        if text.contains("\nedge\t1000") {
            legs += 1;
        }
        if let Err(error) = guard(&text) {
            eprintln!("PICK8 GUARD FAIL {path}: {error:.300}");
            failed.push(path.to_owned());
        }
    }
    eprintln!(
        "PICK8 GUARD list={list} v3={v3} other={other} with_100000_leg={legs} failed={}",
        failed.len()
    );
    ensure!(v3 > 0, "no v3 record");
    ensure!(failed.is_empty(), "{} records fail the guard", failed.len());
    Ok(())
}

// ---------------------------------------------------------------------------
// A REAL positive-degree case: a physical child with a paired crown edge has
// `local_dod <= 0`, so its Taylor order is forced (`with_child_taylor_order`)
// on both paths: check compares today's contracted path with the record, on
// takes the record.
// ---------------------------------------------------------------------------

const L2_22_DOT: &str = include_str!("testdata/u7a/l2-d22.dot");
const L4_3966_DOT: &str = include_str!("testdata/u7a/l4-d3966.dot");

fn forced_child_run(
    diagram: &feynkit_graph::FeynmanDiagram,
    support: &[usize],
    order: Option<i32>,
    switch: Switch,
) -> Result<(Result<Atom>, Vec<Captured>)> {
    let settings = settings();
    let vakint = Vakint::new()?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let child = plan
        .counterterms
        .iter()
        .map(|call| &call.argument.support)
        .find(|child| uv_child_edge_ids(&graph, child) == support)
        .ok_or_else(|| eyre!("{}: no UV child {support:?}", graph.name))?;
    let cache = std::sync::Arc::new(super::CountertermCache::default());
    let run = || {
        with_switch(switch, || {
            with_strict(true, || {
                with_boundary_mode(BoundaryMode::On, false, || {
                    with_capture(|| {
                        RStarPlan::rstar_child_counterterm(
                            &cache,
                            &graph,
                            child,
                            model,
                            &vakint,
                            &settings,
                            plan.uv_children,
                            plan.uv_child_closure,
                        )
                    })
                })
            })
        })
    };
    Ok(match order {
        Some(order) => with_child_taylor_order(support, order, run),
        None => run(),
    })
}

fn check_switch() -> Switch {
    Switch::new(Mode::Check, &[Kind::Root, Kind::Child])
}

fn on_switch() -> Switch {
    Switch::new(Mode::On, &[Kind::Root, Kind::Child])
}

/// The child of `diagram` with edge ids `support` at the forced Taylor order
/// `order`: off, check (== off, every record compared equal) and on (equal
/// to off); at least one evaluated record has `taylor degree=order`, a
/// `100000 + e` leg, and (with `renamed`) an insertion factor naming it.
/// Returns off.
fn positive_degree(
    name: &str,
    dot: &str,
    support: &[usize],
    order: i32,
    renamed: bool,
) -> Result<Atom> {
    let diagram = colour_factored(&load(dot)?)?;
    let (off, nothing) = forced_child_run(&diagram, support, Some(order), Switch::OFF)?;
    let off = off.map_err(|e| eyre!("{name} order {order} off: {e:#}"))?;
    ensure!(nothing.is_empty());
    for (mode, switch) in [("check", check_switch()), ("on", on_switch())] {
        let (value, captured) = forced_child_run(&diagram, support, Some(order), switch)?;
        let value = value.map_err(|e| eyre!("{name} order {order} {mode}: {e:#}"))?;
        let mut paired = 0;
        let mut named = 0;
        for attempt in &captured {
            eprintln!(
                "PICK8 REAL {name} order={order} {mode} kind={} id={} outcome={} equal={:?} first={:?}",
                attempt.kind.as_str(),
                attempt.id,
                attempt.outcome,
                attempt.equal,
                attempt.first_difference
            );
            ensure!(
                attempt.outcome == "uncontracted" || attempt.outcome == "zero",
                "{name} order {order} {mode}: id {} outcome {}",
                attempt.id,
                attempt.outcome
            );
            if mode == "check" && attempt.outcome == "uncontracted" {
                ensure!(
                    attempt.equal == Some(true),
                    "{name} order {order} check: id {} not compared equal",
                    attempt.id
                );
            }
            let Some(text) = attempt.text.as_deref() else {
                continue;
            };
            if attempt.kind == Kind::Child
                && text.contains(&format!("\ntaylor\tdegree={order}\t"))
                && text.contains("\nedge\t1000")
            {
                paired += 1;
                guard(text)
                    .map_err(|e| eyre!("{name} order {order} {mode} id {}: {e}", attempt.id))?;
                if text
                    .lines()
                    .any(|l| l.starts_with("insertion_factor\t") && l.contains("::Q(1000"))
                {
                    named += 1;
                    eprintln!(
                        "PICK8 REAL {name} order={order} {mode} id={} RENAMED leg at degree {order}",
                        attempt.id
                    );
                }
            }
        }
        eprintln!(
            "PICK8 REAL {name} order={order} {mode}: paired_records={paired} renamed={named} zero={}",
            value.is_zero()
        );
        ensure!(
            paired > 0,
            "{name} order {order} {mode}: no evaluated record with a paired leg at that degree"
        );
        ensure!(
            !renamed || named > 0,
            "{name} order {order} {mode}: no renamed leg at that degree"
        );
        if mode == "check" {
            ensure!(
                value == off,
                "{name} order {order}: check does not return today's value"
            );
        } else {
            ensure!(
                atoms_equal(&value, &off),
                "{name} order {order}: on differs from off"
            );
        }
    }
    Ok(off)
}

/// L2 22, child {2,3,5} (crown 0, 1 and the paired gluon 4: legs 4 and
/// 100004), at its own degree 0 and forced to 1 and 2.
#[test]
fn pick8_real_positive_degree_l2_22_child() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let support = [2, 3, 5];
        let zero = positive_degree("L2 22", L2_22_DOT, &support, 0, false)?;
        let one = positive_degree("L2 22", L2_22_DOT, &support, 1, false)?;
        let two = positive_degree("L2 22", L2_22_DOT, &support, 2, false)?;
        eprintln!(
            "PICK8 REAL L2 22: order0 zero={} order1 zero={} order2 zero={} two==zero {}",
            zero.is_zero(),
            one.is_zero(),
            two.is_zero(),
            atoms_equal(&two, &zero)
        );
        ensure!(!two.is_zero(), "order 2 is zero: the test shows nothing");
        ensure!(
            !atoms_equal(&two, &zero),
            "order 2 equals order 0: the test shows nothing"
        );
        Ok(())
    })
}

/// L4 3966, child {2,3,4,5,6,8,9,10} with {9,10} contracted onto the `ct:`
/// vertex 5 next to a `100000 + e` leg, whose insertion factor names it (pick
/// 8's rename), forced to degree 2. About 6 s.
#[test]
fn pick8_real_positive_degree_l4_3966_renamed_child() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let support = [2, 3, 4, 5, 6, 8, 9, 10];
        let two = positive_degree("L4 3966", L4_3966_DOT, &support, 2, true)?;
        ensure!(!two.is_zero(), "order 2 is zero: the test shows nothing");
        Ok(())
    })
}
