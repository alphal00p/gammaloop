//! Remainder-boundary unit tests (DESIGN.md v2, section 11c). The mode is set
//! per test with `with_boundary_mode`, never through the environment.

use eyre::{Result, ensure, eyre};
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Replacement,
    symbol,
};
use vakint::{Vakint, symbols::S};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    utils::{GS, W_},
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    RStarPlan, TaylorExpansion,
    boundary::{
        BoundaryMode, BoundaryStats, ClosedRootNormalization, Grammar, Role,
        denominator_orders, fallback_max_bytes, fallback_max_summands, large_reason, split_into_blocks,
        take_logged_summaries, validate_remainder_boundary, verify_vakint_scalar, with_boundary_mode,
        with_census_dump,
    },
    expression::RStarIntegrand,
    recursive_tests::{colour_factored, load, settings},
};

fn minkowski() -> Atom {
    Minkowski {}.new_rep(GS.dim).to_symbolic([])
}

/// `Q(e, mink(dim))`.
fn q(edge: usize) -> Atom {
    function!(GS.emr_mom, edge, minkowski())
}

/// The spenso dot of two edge momenta.
fn dot(a: usize, b: usize) -> Atom {
    function!(SPENSO_TAG.dot, q(a), q(b))
}

fn reference() -> Atom {
    Atom::var(*super::REFERENCE_MASS_SQUARED)
}

fn auxiliary() -> Atom {
    Atom::var(symbol!("rstar_mass_squared"; Scalar))
}

const BASIS: [usize; 2] = [2, 3];
const X0: usize = 0;

/// Vakint spelling back to spenso: `k(j+1) -> Q(b_j, mink(dim))`,
/// `p(x0) -> Q(x0, mink(dim))`, vakint `dot` -> spenso `dot`.
fn to_spenso(atom: &Atom) -> Atom {
    let mut rules = BASIS
        .iter()
        .enumerate()
        .map(|(j, edge)| Replacement::new(function!(S.k, j as i64 + 1).to_pattern(), q(*edge)))
        .collect::<Vec<_>>();
    rules.push(Replacement::new(function!(S.p, X0 as i64).to_pattern(), q(X0)));
    atom.replace_multiple(&rules)
        .replace(function!(S.dot, W_.a_, W_.b_))
        .with(function!(SPENSO_TAG.dot, W_.a_, W_.b_))
}

fn regulated(atom: &Atom) -> Atom {
    crate::uv::approx::integrated::Integrated::dimensionally_regularized(atom)
}

/// One denominator as `RStarIntegrand::new_with_closure` builds it.
fn den(edge: usize, hard: Atom, momentum: Atom, mass_squared: Atom) -> Atom {
    function!(
        GS.den,
        edge,
        hard,
        mass_squared.sqrt(),
        function!(SPENSO_TAG.dot, &momentum, &momentum) - &mass_squared
    )
    .pow(-1)
}

fn synthetic_denominators() -> Vec<Atom> {
    let bare = |edge: usize| function!(GS.emr_mom, edge);
    vec![
        den(2, bare(2), q(2), Atom::Zero),
        den(3, bare(3), q(3) + q(X0), reference()),
        den(4, bare(2) - bare(3), q(2) - q(3) - q(X0), Atom::Zero),
    ]
}

/// Σ_{i+j≤deg} P_i D_j reassembled from the blocks and the denominator
/// orders: equal to the outer Taylor polynomial of P·D (DESIGN 2c).
#[test]
fn graded_split_matches_the_outer_taylor_polynomial() -> Result<()> {
    test_initialise()?;
    let degree = 2;
    let ms = reference();
    let nf = Atom::var(symbol!("nf"));
    let dim = Atom::var(GS.dim);
    // A remainder numerator after σ: basis Q(2), Q(3), external Q(0); the
    // reference mass to the powers 0, 1, 2; the scalar dimension; an opaque
    // sum holding Q(0) inside a product, which is graded and not expanded.
    let opaque = (dot(2, X0) + Atom::num(2) * dot(3, 3)) * dot(2, 3);
    let p = Atom::num(3) * dot(2, 2) * dot(3, X0) * &ms
        + &dim * dot(2, X0) * dot(3, X0)
        + (dot(2, 3) + dot(2, X0)) * dot(X0, X0) * ms.pow(2) * &nf
        + Atom::num(5) * dot(3, 3) * dot(2, X0).pow(3)
        + opaque * (&dim - 2) * &ms
        + Atom::num(7);
    let denominators = synthetic_denominators();
    let root_edges = [EdgeIndex::from(X0), EdgeIndex::from(1)];
    let reference_symbol = *super::REFERENCE_MASS_SQUARED;
    let auxiliary = auxiliary();

    let mut grammar = Grammar::new(&BASIS, X0, degree, reference_symbol);
    let mut stats = BoundaryStats::default();
    let blocks = split_into_blocks(p.clone(), &mut grammar, &mut stats)?
        .map_err(|decline| eyre!("{}:{}", decline.stage, decline.reason))?;
    ensure!(stats.p_truncated > 0, "the rank-3 term must be truncated");
    let (by_monomial, _) = denominator_orders(
        &denominators,
        &root_edges,
        degree,
        reference_symbol,
        &auxiliary,
    )?
    .map_err(|decline| eyre!("{}:{}", decline.stage, decline.reason))?;
    // D_j summed over the denominator monomials.
    let orders = (0..=degree)
        .map(|j| {
            Atom::add_many(
                by_monomial
                    .iter()
                    .map(|(monomial, orders)| monomial * &orders[j])
                    .collect::<Vec<_>>()
                    .iter(),
            )
        })
        .collect::<Vec<_>>();
    let mut reassembled = Atom::Zero;
    for (grade, block) in &blocks {
        let a = grade.degree();
        let b = grade.mass_power;
        let mut sum = Atom::Zero;
        for c in 0..=b {
            let Some(limit) = degree.checked_sub(a + 2 * c) else {
                break;
            };
            let binomial = (0..c).fold(1i64, |v, i| v * (b - i) as i64 / (i as i64 + 1));
            let sign = if c % 2 == 0 { 1 } else { -1 };
            for order in orders.iter().take(limit + 1) {
                sum = sum + Atom::num(binomial * sign) * order;
            }
        }
        let structure = to_spenso(&grade.structure_atom(grammar.external()));
        reassembled = reassembled + to_spenso(block) * structure * auxiliary.pow(b as i64) * sum;
    }
    let parameter = symbol!("rstar_outer_taylor"; Scalar);
    let taylor = TaylorExpansion::new(parameter, degree);
    let product = &p * Atom::mul_many(denominators.iter());
    let scaled = RStarIntegrand::scale_edges(&product, &root_edges, parameter);
    let expected = taylor.evaluate(&taylor.compensate_mass(&scaled, reference_symbol, &auxiliary))?;
    let difference = (regulated(&reassembled) - regulated(&expected)).expand();
    ensure!(
        difference.is_zero(),
        "the graded split differs from the outer Taylor polynomial: {difference}"
    );
    eprintln!(
        "graded split: {} blocks, {} entries, {} truncated",
        blocks.len(),
        stats.graded_entries,
        stats.p_truncated
    );
    Ok(())
}

fn decline_reason(p: Atom) -> Result<String> {
    let mut grammar = Grammar::new(&BASIS, X0, 2, *super::REFERENCE_MASS_SQUARED);
    let mut stats = BoundaryStats::default();
    match split_into_blocks(p, &mut grammar, &mut stats)? {
        Ok(_) => Ok("accepted".into()),
        Err(decline) => Ok(decline.reason),
    }
}

#[test]
fn graded_split_declines() -> Result<()> {
    test_initialise()?;
    use crate::numerator::aind::Aind;
    use spenso::structure::slot::ParseableAind;
    let indexed = |edge: usize, index: usize| {
        function!(
            GS.emr_mom,
            edge,
            Minkowski {}
                .new_rep(GS.dim)
                .to_symbolic([Aind::Hedge(index as u16, 1).to_atom()])
        )
    };
    // (a) a kept sum holding a non-dot Q(b_j, mink(dim, μ)) with an open index.
    let reason = decline_reason((indexed(2, 1) + indexed(3, 1)) * dot(2, X0) + dot(3, 3))?;
    ensure!(reason.starts_with("indexed"), "(a) {reason}");
    // (b) a trace power.
    let trace = function!(SPENSO_TAG.trace, Atom::var(symbol!("boundary_test_x"))).pow(2);
    let reason = decline_reason(trace * dot(2, 2))?;
    ensure!(reason.starts_with("indexed"), "(b) {reason}");
    // (c) a closed indexed pair: the retry maps it to a dot; an open index
    // declines.
    let pair = indexed(2, 5) * indexed(3, 5);
    let notated = super::boundary::retry_contraction(pair.as_view())?;
    let reason = decline_reason(pair * dot(2, X0) + dot(3, 3))?;
    ensure!(reason == "accepted", "(c) closed pair: {reason}; the notation pass gives {notated}");
    let reason = decline_reason(indexed(2, 5) * indexed(3, 6) * dot(2, X0) + dot(3, 3))?;
    ensure!(reason.starts_with("indexed"), "(c) open pair: {reason}");
    // (d) dot(Q(x0),Q(x0))^-1.
    let reason = decline_reason(dot(X0, X0).pow(-1) * dot(2, 2) + dot(3, 3))?;
    ensure!(reason == "graded_power", "(d) {reason}");
    // (e) (M_ref^2)^(1/2).
    let reason = decline_reason(reference().sqrt() * dot(2, 2) + dot(3, 3))?;
    ensure!(reason == "fractional_reference_mass_power", "(e) {reason}");
    // (f) Q(x1) (σ not applied).
    let reason = decline_reason(dot(2, 1) + dot(3, 3))?;
    ensure!(reason == "dot_argument:Q(1)", "(f) {reason}");
    // (g) an unrouted Q(7).
    let reason = decline_reason(dot(2, 7) + dot(3, 3))?;
    ensure!(reason == "dot_argument:Q(7)", "(g) {reason}");
    // (h) a den in P.
    let reason = decline_reason(synthetic_denominators()[0].clone() * dot(2, 2) + dot(3, 3))?;
    // `den` is feynkit's `denom` head.
    ensure!(reason == "function:denom", "(h) {reason}");
    Ok(())
}

/// DESIGN 9: the split checkpoints inside a factorized P (one top-level
/// product whose sum has more than 4,096 terms, the shape idenso returns), and
/// an exhausted resource budget is a hard error, never a decline.
#[test]
fn split_checkpoints_inside_a_factorized_p() -> Result<()> {
    use std::time::Duration;
    use symbolica::utils::resources::{ResourceBudget, ResourceLimits};
    test_initialise()?;
    let terms = (1..=5000).map(|i| dot(2, 2).pow(i as i64)).collect::<Vec<_>>();
    let p = Atom::add_many(terms.iter()) * dot(3, 3);
    ensure!(matches!(p.as_view(), AtomView::Mul(_)), "P must be one top-level product");
    let split = |p: Atom| -> Result<std::result::Result<usize, String>> {
        let mut grammar = Grammar::new(&BASIS, X0, 2, *super::REFERENCE_MASS_SQUARED);
        let mut stats = BoundaryStats::default();
        Ok(split_into_blocks(p, &mut grammar, &mut stats)?
            .map(|blocks| blocks.len())
            .map_err(|decline| decline.reason))
    };
    let unbounded = split(p.clone())?;
    ensure!(unbounded == Ok(1), "without a budget: {unbounded:?}");
    let budget = ResourceBudget::new(ResourceLimits {
        max_elapsed: Some(Duration::ZERO),
        ..Default::default()
    });
    let exhausted = budget.run(|| split(p));
    ensure!(exhausted.is_err(), "an exhausted budget must be a hard error: {exhausted:?}");
    Ok(())
}

/// The thread-local scopes of check mode (closure capture, callback check)
/// and of the dump channel are restored when the integration inside them
/// panics: beta catches a panicking diagram and the worker thread goes on
/// with the next one, which must not inherit them.
#[test]
fn scopes_are_restored_when_an_integration_panics() -> Result<()> {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    let normalization = ClosedRootNormalization {
        external: EdgeIndex::from(X0),
        external_routing: Vec::new(),
        factor: Atom::num(1),
    };
    let channel = || super::pint6_dump::ClosureChannel {
        channel: "projection metric",
        taylor_order: 2,
        normalisation: "1/(i p^2 d_A (d-1))",
    };
    let caught = catch_unwind(AssertUnwindSafe(|| {
        super::integrate::capturing_closures(|| {
            super::boundary::with_callback_check(&normalization, || {
                super::pint6_dump::with_channel(channel(), || -> usize {
                    panic!("test: an exhausted budget inside a contraction")
                })
            })
        })
    }));
    ensure!(caught.is_err(), "the panic must reach the caller");
    ensure!(!super::integrate::capture_active(), "the closure capture is left set");
    ensure!(!super::boundary::callback_check_active(), "the callback check is left set");
    ensure!(!super::pint6_dump::channel_active(), "the channel is left set");
    // A normal end returns the scope's values and restores the slots too.
    let ((value, outcome), capture) = super::integrate::capturing_closures(|| {
        super::boundary::with_callback_check(&normalization, || {
            super::pint6_dump::with_channel(channel(), || {
                (super::integrate::capture_active(), super::pint6_dump::channel_active())
            })
        })
    });
    ensure!(value == (true, true), "the scopes are not active inside: {value:?}");
    ensure!(outcome == "none" && capture.closures.is_empty(), "{outcome}");
    ensure!(!super::integrate::capture_active() && !super::boundary::callback_check_active());
    ensure!(!super::pint6_dump::channel_active());
    Ok(())
}

#[test]
fn postconditions_reject_unmapped_loop_momenta() -> Result<()> {
    test_initialise()?;
    let reference = *super::REFERENCE_MASS_SQUARED;
    let k = |n: i64| function!(S.k, n);
    let p = |n: i64| function!(S.p, n);
    let vdot = |a: &Atom, b: &Atom| function!(S.dot, a, b);
    let check = |atom: Atom, role: Role| {
        verify_vakint_scalar(atom.as_view(), role, 2, X0, 1, reference)
            .err()
            .map(|decline| decline.reason)
            .unwrap_or_else(|| "accepted".into())
    };
    // Accepted shapes.
    ensure!(check(vdot(&k(1), &k(2)) * Atom::var(symbol!("nf")), Role::Block) == "accepted");
    ensure!(check(vdot(&k(1), &p(0)) * vdot(&p(0), &p(0)).pow(-1), Role::Coefficient) == "accepted");
    ensure!(check(k(1) - p(1), Role::Routing) == "accepted");
    // Rejected ones.
    let reason = check(function!(GS.emr_mom, 2) * vdot(&k(1), &k(1)), Role::Block);
    ensure!(reason.starts_with("postcondition:block:"), "emr_mom: {reason}");
    let reason = check(k(1) * vdot(&k(1), &k(2)), Role::Block);
    ensure!(reason == "postcondition:block:k(1)", "bare k: {reason}");
    let reason = check(vdot(&p(1), &k(1)), Role::Block);
    ensure!(reason == "postcondition:block:p(1)", "p in a block: {reason}");
    let reason = check(vdot(&p(2), &k(1)), Role::Coefficient);
    ensure!(reason == "postcondition:coefficient:p(2)", "p(b_j) in R: {reason}");
    let reason = check(vdot(&k(3), &k(1)), Role::Coefficient);
    ensure!(reason == "postcondition:coefficient:k(3)", "k beyond L: {reason}");
    let reason = check(k(2) + p(1), Role::Routing);
    ensure!(reason == "postcondition:routing:k(2)", "vacuum label in a routing image: {reason}");
    let reason = check(Atom::var(GS.dim) * vdot(&k(1), &k(1)), Role::Block);
    ensure!(reason.starts_with("postcondition:block:"), "dim: {reason}");
    Ok(())
}

/// DESIGN 10c with Addendum C and the step 1-3 review: one predicate, before
/// and after the contraction. Six loops or a census dump are always large;
/// after the contraction, P above the byte or the summand limit. A factorized
/// P (one top-level product, summand count 1) is large by its loops, the
/// dump or its bytes, never by its summand count alone.
#[test]
fn large_rule_holds_before_and_after_the_contraction() {
    let mut stats = BoundaryStats::default();
    // Before the contraction.
    assert!(large_reason(&stats, 5, false).is_none());
    assert!(large_reason(&stats, 6, false).is_some());
    assert!(large_reason(&stats, 2, true).is_some());
    // After it, with a factorized P: one summand.
    stats.contracted = true;
    stats.p_summands = 1;
    stats.p_bytes = 3_341_388;
    assert!(large_reason(&stats, 5, false).is_none());
    assert!(large_reason(&stats, 6, false).is_some(), "six loops after the contraction");
    assert!(large_reason(&stats, 5, true).is_some(), "a census dump after the contraction");
    stats.p_bytes = fallback_max_bytes().saturating_add(1);
    assert!(fallback_max_bytes() == usize::MAX || large_reason(&stats, 5, false).is_some());
    stats.p_bytes = 10;
    stats.p_summands = fallback_max_summands().saturating_add(1);
    assert!(fallback_max_summands() == usize::MAX || large_reason(&stats, 2, false).is_some());
}

/// The physical external rules of a two-point graph, as beta builds them:
/// `(x0, Q(e, x___) -> ±Q(x0, x___))`.
fn external_rules(graph: &Graph) -> (EdgeIndex, Vec<Replacement>) {
    let external = graph
        .iter_edges()
        .filter_map(|(pair, edge, _)| match pair {
            HedgePair::Unpaired { flow, .. } => Some((edge, flow)),
            _ => None,
        })
        .collect::<Vec<_>>();
    let (p, p_flow) = external[0];
    let rules = external
        .iter()
        .map(|(edge, flow)| {
            let sign = if *edge == p || *flow != p_flow { 1 } else { -1 };
            Replacement::new(
                function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                (Atom::num(sign) * function!(GS.emr_mom, usize::from(p), W_.x___)).to_pattern(),
            )
        })
        .collect();
    (p, rules)
}

/// A closed-root normalisation that is `factor · σ` by construction (a
/// stand-in for beta's callback on closed scalars).
fn test_normalization(graph: &Graph) -> ClosedRootNormalization {
    let (external, external_routing) = external_rules(graph);
    let x0 = usize::from(external);
    let q_squared = function!(SPENSO_TAG.dot, q(x0), q(x0));
    let factor = Atom::var(symbol!("boundary_test_channel"))
        * graph.global_prefactor.num.replace_multiple(&external_routing)
        / (Atom::i() * q_squared)
        / (Atom::var(GS.dim) - 1);
    ClosedRootNormalization {
        external,
        external_routing,
        factor,
    }
}

/// Off == On == Check on production background diagrams (colour factored),
/// through the closing-root `_scalar` entry point, under STRICT: every
/// admitted remainder is built and evaluated by the boundary (On), or built
/// and compared pair by pair with today's closures, value and callback
/// (Check, which errors on any MISMATCH and returns today's value);
/// children and one-loop remainders are silently not admitted.
#[test]
fn boundary_matches_old_path_two_and_three_loop() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let settings = settings();
        let mut total_boundaries = 0;
        for (name, dot) in [
            ("L2/0", include_str!("testdata/l2-d0.dot")),
            ("L2/8", include_str!("testdata/l2-d8.dot")),
            ("L2/12", include_str!("testdata/l2-d12.dot")),
            ("L3/185", include_str!("testdata/l3-d185.dot")),
            ("L3/241", include_str!("testdata/l3-d241.dot")),
        ] {
            let diagram = colour_factored(&load(dot)?)?;
            let graph = diagram.to_gamma_loop_graph(None, true)?;
            let model = diagram.model();
            let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
            let closed_root = test_normalization(&graph);
            let callback = |tensor: &Atom| -> Result<Atom> {
                Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
            };
            let run = |mode: BoundaryMode| {
                with_boundary_mode(mode, true, || {
                    plan.projected_counterterm_closing_root_scalar(
                        &graph,
                        model,
                        &vakint,
                        &settings,
                        &callback,
                        &closed_root,
                    )
                })
            };
            // The summaries of a run, summed: (boundary, zero, checks,
            // fallback, errors, uncertified); every one must end `ok`.
            let totals = |mode: BoundaryMode| -> Result<[usize; 6]> {
                let mut totals = [0; 6];
                for (root, outcome, summary) in take_logged_summaries() {
                    ensure!(outcome == "ok", "{name} {mode:?}: root {root} ended {outcome}");
                    for (total, count) in totals.iter_mut().zip([
                        summary.boundary,
                        summary.zero,
                        summary.checks,
                        summary.fallback,
                        summary.errors,
                        summary.uncertified,
                    ]) {
                        *total += count;
                    }
                }
                Ok(totals)
            };
            let _ = take_logged_summaries();
            let off = run(BoundaryMode::Off)?;
            ensure!(take_logged_summaries().is_empty(), "{name}: mode off logs no summary");
            let on = run(BoundaryMode::On)?;
            let [boundaries, _, _, fallback, errors, uncertified] = totals(BoundaryMode::On)?;
            // Without this, On == Off would hold trivially if nothing were
            // admitted (A7-A9 are silent under STRICT).
            ensure!(boundaries >= 1, "{name}: no remainder took the boundary path in mode on");
            ensure!(
                fallback == 0 && errors == 0 && uncertified == 0,
                "{name} on: fallback={fallback} errors={errors} uncertified={uncertified}"
            );
            let check = run(BoundaryMode::Check)?;
            let [checked, zero, checks, fallback, errors, _] = totals(BoundaryMode::Check)?;
            ensure!(
                checked == boundaries && checks == checked + zero,
                "{name} check: {checked} built and {checks} checks against {boundaries} in mode on"
            );
            ensure!(fallback == 0 && errors == 0, "{name} check: fallback={fallback} errors={errors}");
            total_boundaries += boundaries;
            ensure!(!off.expand().is_zero(), "{name}: zero counterterm");
            let identical = off == on;
            let difference = (&off - &on).expand();
            ensure!(
                identical || difference.is_zero(),
                "{name}: the boundary path differs from today's: {difference}"
            );
            ensure!(check == off, "{name}: check mode must return today's value");
            eprintln!(
                "BOUNDARY_TEST {name} on identical={identical} boundaries={boundaries} check ok"
            );
        }
        eprintln!("BOUNDARY_TEST boundaries={total_boundaries}");
        Ok(())
    })
}

/// Mode on under a census dump is accepted at start-up since the boundary
/// dump exists (DESIGN 6b; the refusal of 6e49bae5e5 is lifted), in every mode.
#[test]
fn census_dump_is_accepted_in_every_mode() -> Result<()> {
    test_initialise()?;
    for mode in [BoundaryMode::On, BoundaryMode::Check, BoundaryMode::Off] {
        for dump in [true, false] {
            let startup = with_census_dump(dump, || {
                with_boundary_mode(mode, false, validate_remainder_boundary)
            });
            ensure!(startup.is_ok(), "{mode:?} dump={dump}: start-up {startup:?}");
        }
    }
    Ok(())
}

/// The boundary record (DESIGN 6b, Addendum B): written with a streamed block
/// file and an inline block, read back by the streaming reader into the same
/// routing (plus the closing route, A1), blocks and patterns; the block file
/// has one line per summand (A3); a newline inside an atom is refused (A2);
/// the boundary and v1 counters are independent.
#[test]
fn boundary_record_round_trips() -> Result<()> {
    test_initialise()?;
    let parse = |text: &str| vakint::vakint_parse!(text).map_err(|error| eyre!("{error}"));
    let directory = std::env::temp_dir().join(format!("pint6-boundary-test-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&directory);
    let dump = super::pint6_dump::test_settings(&directory)?;
    // A block of 100 summands (file) and a small one (inline).
    let large = Atom::add_many(
        (1..=100)
            .map(|n| parse(&format!("{n}*dot(k(1),k(2))^{n}*dot(k(2),k(2))")))
            .collect::<Result<Vec<_>>>()?
            .iter(),
    );
    let small = parse("dot(k(1),k(1))+3/7*dot(k(1),k(2))")?;
    let topology = parse("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,2))")?;
    let routing = vec![(parse("k(1)")?, parse("k(1)")?), (parse("k(2)")?, parse("k(1)-p(1)")?)];
    let boundary = vakint::p_integral::RemainderBoundary::new(
        "graph=test marked_edge=2 integration_id=1 loops=2",
        parse("p(1)")?,
        Atom::var(symbol!("rstar_mass_squared"; Scalar)),
        routing.clone(),
        vec![large, small],
        vec![vakint::p_integral::BoundaryPattern {
            topology: topology.clone(),
            closure_power: 1,
            terms: vec![(0, parse("dot(p(0),p(0))")?), (1, parse("2*ε")?)],
        }],
    )
    .map_err(|error| eyre!("{error}"))?;
    let meta = super::boundary::PatternMeta {
        order: 0,
        topology: topology.clone(),
        key: topology.clone(),
        integral: topology,
        closure_power: 1,
        unprojected_nonzero: true,
    };
    let built = super::boundary::BuiltBoundary {
        boundary,
        patterns: vec![meta],
        zero_patterns: Vec::new(),
        grades: vec![super::boundary::Grade::default(), super::boundary::Grade::default()],
        loops: 2,
        basis: vec![2, 3],
        normalisation: Atom::num(1),
        external: parse("p(0)")?,
        degree: 2,
        integration_id: 1,
    };
    let vakint_settings = settings().vakint.true_settings();
    let references = [(0usize, 7usize)].into_iter().collect();
    let first = super::pint6_dump::write_boundary(
        &dump,
        &vakint_settings,
        &built,
        &serde_json::json!({"test": true}),
        false,
        &references,
    )?;
    let second = super::pint6_dump::write_boundary(
        &dump,
        &vakint_settings,
        &built,
        &serde_json::json!({}),
        false,
        &Default::default(),
    )?;
    ensure!(second.sequence == first.sequence + 1, "boundary sequence {first:?}");
    ensure!(first.block_files == 1, "one streamed block file expected");
    super::pint6_dump::verify_boundary(&built.boundary, &first.header)?;
    let (read, keys) = super::pint6_dump::read_boundary(&first.header)?;
    ensure!(read.routing.len() == routing.len() + 1, "the closing route (A1)");
    ensure!(keys.get("routing").map(String::as_str) == Some("3"), "routing count {keys:?}");
    let header = std::fs::read_to_string(&first.header)?;
    ensure!(header.contains("\tref_closure=7\t"), "ref_closure missing:\n{header}");
    ensure!(header.contains("\nroute "), "route lines missing");
    // A3: the block file has exactly `summands` lines after its three head lines.
    let block = std::fs::read_dir(&directory)?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .find(|path| path.to_string_lossy().ends_with(".block-000.txt.zst"))
        .ok_or_else(|| eyre!("no block file"))?;
    let text = std::process::Command::new("zstd").arg("-dcq").arg(&block).output()?.stdout;
    let text = String::from_utf8(text)?;
    ensure!(text.lines().count() == 3 + 100, "block lines");
    ensure!(text.lines().nth(2) == Some("summands 100"), "summands line");
    let _ = std::fs::remove_dir_all(&directory);
    Ok(())
}

/// The vacuum projection of a boundary coefficient does not depend on the
/// integral it is attached to (one image per remainder), and odd ranks in
/// `k·p` vanish (DESIGN 4a, OPEN-7).
#[test]
fn feynkit_images_are_pattern_independent_and_odd_images_vanish() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    let mut tensor_settings = settings().vakint.true_settings();
    tensor_settings.tensor_reduction_method = vakint::TensorReductionMethod::FeynKit;
    tensor_settings.use_dot_product_notation = true;
    let parse = |text: &str| vakint::vakint_parse!(text).map_err(|error| eyre!("{error}"));
    let first = parse("topo(prop(1,edge(1,2),k(0),M2,1)*prop(2,edge(1,2),k(1),0,1)*prop(3,edge(2,1),k(0)+k(1),0,1))")?;
    let second = parse("topo(prop(1,edge(1,2),k(0),M2,2)*prop(2,edge(1,2),k(1),0,3)*prop(3,edge(2,1),k(0)+k(1),0,1))")?;
    let project = |integral: &Atom, numerator: &str| -> Result<Atom> {
        let numerator = parse(numerator)?;
        let vectors = vakint::VakintTerm::identify_vectors_in_numerator(numerator.as_view())
            .map_err(|error| eyre!("{error}"))?;
        let term = vakint::VakintTerm {
            integral: integral.clone(),
            numerator,
            vectors,
        };
        let mut images = std::collections::HashMap::new();
        Ok(super::integrate::tensor_reduce_by_factors(&vakint, &tensor_settings, &term, &mut images)?
            .unwrap_or_default()
            .expand())
    };
    let rank2 = "dot(k(1),p(0))*dot(k(2),p(0))*dot(k(1),k(2))";
    let a = project(&first, rank2)?;
    let b = project(&second, rank2)?;
    ensure!(!a.is_zero(), "a zero rank-2 image");
    ensure!(a == b, "the image depends on the integral: {a} against {b}");
    for odd in [
        "dot(k(1),p(0))*dot(k(2),k(2))",
        "dot(k(1),p(0))*dot(k(2),p(0))^2",
        "dot(k(1),p(0))^3*dot(k(1),k(2))",
    ] {
        let image = project(&first, odd)?;
        ensure!(image.is_zero(), "odd image {odd} -> {image}");
    }
    Ok(())
}

/// The `AtomView` kind, for messages.
#[allow(dead_code)]
fn kind(view: AtomView<'_>) -> &'static str {
    match view {
        AtomView::Num(_) => "number",
        AtomView::Var(_) => "symbol",
        AtomView::Add(_) => "sum",
        AtomView::Mul(_) => "product",
        AtomView::Pow(_) => "power",
        AtomView::Fun(_) => "function",
    }
}
