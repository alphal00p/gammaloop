//! `GAMMALOOP_UV_CHILD_TRUNCATE`: the root-grade truncation of the ordinary
//! remainders of a UV child evaluated as its own R* root
//! (results/uv-degree-truncation/DESIGN.md, "Tests and gate").
//!
//! The mode is read once per process; these tests call the builder with an
//! explicit mode, so they hold whatever the environment says. Running the
//! whole `uv::rstar` suite with the switch set checks every child against the
//! massive R operation in that mode.

use std::{collections::BTreeSet, sync::Arc};

use eyre::{Result, ensure, eyre};
use linnet::half_edge::involution::EdgeIndex;
use spenso::network::tags::SPENSO_TAG;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function, symbol,
};
use vakint::Vakint;

use super::{
    ChildTruncate, IntegrationTimings, RStarIntegrand, TaylorExpansion, UvStream,
    truncate_routed, truncation_degree, truncation_dump_label, truncation_key_suffix,
    uv_child_truncation_key,
};
use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    utils::{GS, W_},
    uv::{UltravioletGraph, approx::streamed::with_production_piece_limit},
};

use super::super::{
    ChildMass, CountertermCache, REFERENCE_MASS_SQUARED, RStarPlan, UvChildClosure,
    UvChildrenRoute,
    recursive_tests::{colour_factored, load, settings},
    tests::canonical_tensor,
};

/// `Q(edge, slot)`.
fn q(edge: i64, slot: &str) -> Atom {
    function!(GS.emr_mom, Atom::num(edge), Atom::var(symbol!(slot)))
}

/// `Q(edge)`, a vector without a slot (as inside `dot`).
fn p(edge: i64) -> Atom {
    function!(GS.emr_mom, Atom::num(edge))
}

fn dot(a: &Atom, b: &Atom) -> Atom {
    function!(SPENSO_TAG.dot, a, b)
}

const MU: &str = "uv_truncate_test::mu";
const NU: &str = "uv_truncate_test::nu";

/// `den(edge, Q(hard), 0, dot(k, k))^-1` with `k` mixing loop and root momenta.
fn den(edge: i64, hard: i64, momentum: &Atom) -> Atom {
    function!(GS.den, Atom::num(edge), p(hard), Atom::num(0), dot(momentum, momentum)).pow(-1)
}

fn edges(ids: &[usize]) -> Vec<EdgeIndex> {
    ids.iter().map(|id| EdgeIndex::from(*id)).collect()
}

fn loops(ids: &[usize]) -> BTreeSet<usize> {
    ids.iter().copied().collect()
}

/// The denominators of the toy remainders: their fourth argument carries the
/// root momenta, as after routing.
fn denominators() -> Atom {
    den(6, 6, &(p(6) + p(2))) * den(7, 7, &(p(7) - p(3)))
}

fn truncate(
    numerator: &Atom,
    roots: &[usize],
    degree: usize,
    max_parts: usize,
) -> std::result::Result<Vec<Atom>, &'static str> {
    let routed = numerator * denominators();
    truncate_routed(&routed, &edges(roots), &loops(&[6, 7]), degree, max_parts)
        .expect("no expansion error")
        .map(|truncated| truncated.parts)
}

fn sum(parts: &[Atom]) -> Atom {
    Atom::add_many(parts.iter())
}

fn assert_equal(actual: &Atom, expected: &Atom) {
    assert!(
        (actual - expected).expand().is_zero(),
        "actual: {actual}\nexpected: {expected}"
    );
}

/// Every `Q(root, ..)` set to zero, the reference for degree 0.
fn zero_roots(expression: &Atom, roots: &[usize]) -> Atom {
    let mut zeroed = expression.clone();
    for root in roots {
        zeroed = zeroed
            .replace(function!(GS.emr_mom, Atom::num(*root as i64), W_.x___))
            .with(Atom::Zero);
    }
    zeroed
}

/// The outer Taylor polynomial of the numerator in the root momenta, at one.
fn graded_reference(numerator: &Atom, roots: &[usize], degree: usize) -> Result<Atom> {
    let parameter = symbol!("uv_truncate_test::t"; Scalar);
    let scaled = RStarIntegrand::scale_edges(numerator, &edges(roots), parameter);
    TaylorExpansion::new(parameter, degree).evaluate(&scaled)
}

#[test]
fn truncation_degree_follows_the_mode_the_scope_and_the_order() {
    use ChildTruncate::*;
    for order in [-1, 0, 1, 2] {
        for mode in [Off, DegreeZero, Graded] {
            assert_eq!(truncation_degree(mode, false, order), None, "{mode:?} {order}");
        }
        assert_eq!(truncation_degree(Off, true, order), None);
        assert_eq!(
            truncation_degree(DegreeZero, true, order),
            (order == 0).then_some(0)
        );
        assert_eq!(
            truncation_degree(Graded, true, order),
            usize::try_from(order).ok()
        );
    }
}

#[test]
fn the_switch_accepts_only_its_modes() {
    use ChildTruncate::*;
    assert_eq!(ChildTruncate::parse(None), Ok(Off));
    for (value, mode) in [
        ("", Off),
        ("  ", Off),
        ("0", Off),
        (" 0\n", Off),
        ("1", DegreeZero),
        (" 1 ", DegreeZero),
        ("graded", Graded),
        ("\tgraded ", Graded),
    ] {
        assert_eq!(ChildTruncate::parse(Some(value)), Ok(mode), "{value:?}");
    }
    for value in ["off", "on", "2", "01", "true", "Graded", "GRADED", "1 graded", "graded=1"] {
        let error = ChildTruncate::parse(Some(value)).expect_err(value);
        assert!(error.contains("GAMMALOOP_UV_CHILD_TRUNCATE"), "{error}");
        assert!(error.contains(value.trim()), "{error}");
    }
    // One line per process, only when on: the off-mode logs are unchanged.
    assert_eq!(Off.mode_line(32), None);
    assert_eq!(
        DegreeZero.mode_line(32).as_deref(),
        Some("UV_TRUNCATE_MODE mode=1 max_parts=32")
    );
    assert_eq!(
        Graded.mode_line(7).as_deref(),
        Some("UV_TRUNCATE_MODE mode=graded max_parts=7")
    );
}

#[test]
fn closure_dump_label_is_present_only_when_on() {
    use ChildTruncate::*;
    assert_eq!(truncation_dump_label(Off, None), None);
    assert_eq!(truncation_dump_label(Off, Some(0)), None);
    assert_eq!(truncation_dump_label(DegreeZero, Some(0)).as_deref(), Some("1@0"));
    assert_eq!(truncation_dump_label(DegreeZero, None).as_deref(), Some("1@-"));
    assert_eq!(truncation_dump_label(Graded, Some(2)).as_deref(), Some("graded@2"));
    assert_eq!(truncation_dump_label(Graded, None).as_deref(), Some("graded@-"));
}

#[test]
fn cache_key_suffix_is_present_only_when_on() {
    use ChildTruncate::*;
    assert_eq!(truncation_key_suffix(Off, 0), None);
    assert_eq!(truncation_key_suffix(Off, 2), None);
    assert_eq!(
        truncation_key_suffix(DegreeZero, 0).as_deref(),
        Some("#uv-child-truncate=1@0")
    );
    assert_eq!(
        truncation_key_suffix(DegreeZero, 1).as_deref(),
        Some("#uv-child-truncate=1@1")
    );
    assert_eq!(
        truncation_key_suffix(Graded, 2).as_deref(),
        Some("#uv-child-truncate=graded@2")
    );
    // The off keys are byte-identical to the keys before the switch.
    let route = UvChildrenRoute::Rstar { min_loops: 1 };
    let closure = UvChildClosure::Isotropic { min_loops: 1 };
    let settings = r#"{"a":1}"#.to_string();
    assert_eq!(
        RStarPlan::rstar_child_key(settings.clone(), route, UvChildClosure::Off, None),
        format!("{settings}#uv-children={route}")
    );
    assert_eq!(
        RStarPlan::rstar_child_key(settings.clone(), route, closure, None),
        format!("{settings}#uv-children={route}#child-closure={closure}")
    );
    // On: the suffix follows, and keys of the two modes or degrees differ.
    let on = |mode, order| {
        RStarPlan::rstar_child_key(
            settings.clone(),
            route,
            closure,
            truncation_key_suffix(mode, order),
        )
    };
    assert_eq!(
        on(DegreeZero, 0),
        format!("{settings}#uv-children={route}#child-closure={closure}#uv-child-truncate=1@0")
    );
    let keys = [on(Off, 0), on(DegreeZero, 0), on(Graded, 0), on(Graded, 1)];
    assert_eq!(keys.iter().collect::<BTreeSet<_>>().len(), keys.len());
    // The process mode: unset means no suffix and no degree computed.
    let mut called = false;
    let suffix = uv_child_truncation_key(|| {
        called = true;
        0
    });
    match ChildTruncate::parse(std::env::var("GAMMALOOP_UV_CHILD_TRUNCATE").ok().as_deref()) {
        Ok(Off) => assert!(suffix.is_none() && !called),
        _ => assert!(suffix.is_some() && called),
    }
}

#[test]
fn slotted_root_vectors_are_graded() {
    // Q(2,μ) Q(6,μ): grade 1. Zero at degree 0, unchanged at degree 1.
    let numerator = q(2, MU) * q(6, MU);
    assert_eq!(truncate(&numerator, &[2], 0, 32), Ok(vec![]));
    let parts = truncate(&numerator, &[2], 1, 32).unwrap();
    assert_eq!(parts.len(), 1);
    assert_equal(&sum(&parts), &(&numerator * denominators()));
    // A sum mixing root and loop vectors in one slot.
    let numerator = (q(2, MU) + q(6, MU)) * (q(3, NU) - q(7, NU)) * q(7, MU) * q(6, NU);
    let parts = truncate(&numerator, &[2, 3], 0, 32).unwrap();
    assert_eq!(parts.len(), 1);
    assert_equal(
        &parts[0],
        &(Atom::num(-1) * q(6, MU) * q(7, NU) * q(7, MU) * q(6, NU) * denominators()),
    );
}

#[test]
fn root_vectors_inside_dot_are_graded() {
    // dot(Q2 + Q6, Q7) at degree 0 is dot(Q6, Q7): dot is linear under the
    // scaling, as in the outer Taylor step.
    let numerator = dot(&(p(2) + p(6)), &p(7));
    let parts = truncate(&numerator, &[2], 0, 32).unwrap();
    assert_eq!(parts.len(), 1);
    assert_equal(&parts[0], &(dot(&p(6), &p(7)) * denominators()));
    // dot(Q2, Q3) has grade 2, dot(Q2, Q6) grade 1.
    let numerator = dot(&p(2), &p(3)) + dot(&p(6), &p(7));
    let at = |numerator: &Atom, degree| sum(&truncate(numerator, &[2, 3], degree, 32).unwrap());
    assert_equal(&at(&numerator, 1), &(dot(&p(6), &p(7)) * denominators()));
    assert_equal(&at(&numerator, 2), &(&numerator * denominators()));
    let numerator = dot(&p(2), &p(6)) + dot(&p(6), &p(7));
    assert_equal(&at(&numerator, 0), &(dot(&p(6), &p(7)) * denominators()));
    assert_equal(&at(&numerator, 1), &(&numerator * denominators()));
    // A squared root vector and a mixed product.
    let numerator = dot(&(p(2) - p(6)), &(p(2) + p(7))) * q(3, MU) * q(7, MU)
        + Atom::num(3) * dot(&p(6), &p(6)) * q(6, MU) * q(7, MU);
    for degree in 0..=3 {
        let truncated = sum(&truncate(&numerator, &[2, 3], degree, 32).unwrap());
        assert_equal(
            &truncated,
            &(graded_reference(&numerator, &[2, 3], degree).unwrap() * denominators()),
        );
    }
}

#[test]
fn nested_combinations_cancel_before_the_rename() {
    // Q(2,μ) - (Q(2,μ) + Q(3,μ) - Q(6,μ)) with R = {2, 3}: Q(6,μ), one vector.
    let nested = q(2, MU) - (q(2, MU) + q(3, MU) - q(6, MU));
    let numerator = &nested * q(7, MU);
    let parts = truncate(&numerator, &[2, 3], 0, 32).unwrap();
    assert_eq!(parts.len(), 1);
    assert_equal(&parts[0], &(q(6, MU) * q(7, MU) * denominators()));
    let mut combinations = Default::default();
    let renamed = super::rename_edge_basis(&parts[0], true, &mut combinations);
    // Only the denominators' momenta remain as combinations.
    let vectors = renamed
        .as_view()
        .to_owned()
        .replace(function!(GS.den, W_.x___))
        .with(Atom::one());
    assert_eq!(vectors, q(6, MU) * q(7, MU), "{renamed}");
}

#[test]
fn denominators_are_kept_byte_identical() {
    let numerator = (q(2, MU) + q(6, MU)) * q(7, MU) * dot(&p(3), &p(6));
    let routed = &numerator * denominators();
    let AtomView::Mul(product) = routed.as_view() else {
        panic!("a product");
    };
    let originals = product
        .iter()
        .filter(|factor| factor.contains_symbol(GS.den))
        .map(|factor| factor.to_owned())
        .collect::<Vec<_>>();
    assert_eq!(originals.len(), 2);
    for degree in [1, 2] {
        for part in truncate(&numerator, &[2, 3], degree, 32).unwrap() {
            let AtomView::Mul(product) = part.as_view() else {
                panic!("a product part: {part}");
            };
            let kept = product
                .iter()
                .filter(|factor| factor.contains_symbol(GS.den))
                .map(|factor| factor.to_owned())
                .collect::<Vec<_>>();
            assert_eq!(kept, originals);
        }
    }
}

#[test]
fn masses_are_grade_zero() {
    let reference = Atom::var(*REFERENCE_MASS_SQUARED);
    let m_uv = Atom::var(GS.m_uv_vacuum);
    let numerator = (&reference + dot(&p(2), &p(6))) * (&m_uv * q(6, MU) + q(3, MU)) * q(7, MU);
    let parts = truncate(&numerator, &[2, 3], 0, 32).unwrap();
    assert_eq!(parts.len(), 1);
    assert_equal(
        &parts[0],
        &(&reference * &m_uv * q(6, MU) * q(7, MU) * denominators()),
    );
}

#[test]
fn degree_zero_equals_the_zeroth_order() {
    let a = Atom::var(symbol!("uv_truncate_test::a"));
    let numerators = [
        (q(2, MU) + q(6, MU)) * (q(6, NU) - q(3, NU) + q(7, NU)) * q(7, MU) * q(2, NU),
        (dot(&p(2), &p(7)) + &a * dot(&p(6), &p(7))).pow(2) * (q(3, MU) - q(6, MU)) * q(7, MU),
        (q(6, MU) - (q(6, MU) + q(2, MU) - q(7, MU))) * q(6, MU) + &a,
        dot(&(p(2) + p(3) + p(6)), &(p(7) - p(2))) * (Atom::num(2) * q(2, MU) + q(7, MU)) * q(6, MU),
    ];
    for numerator in &numerators {
        let truncated = sum(&truncate(numerator, &[2, 3], 0, 32).unwrap());
        let order = graded_reference(numerator, &[2, 3], 0).unwrap();
        assert_equal(&truncated, &(&order * denominators()));
        // The slotted vectors zeroed by hand, after expanding the dots.
        let zeroed = zero_roots(&numerator.expand(), &[2, 3]);
        assert_equal(&order, &zeroed);
    }
}

#[test]
fn graded_parts_add_up_to_the_truncated_series() {
    // k dependent factors, each with nonzero grades 0, 1 and 2.
    let factor = |i: i64| {
        Atom::num(i + 1) * dot(&p(6), &p(7))
            + dot(&p(2), &p(6 + i % 2))
            + Atom::num(i) * dot(&p(2), &p(3))
            + dot(&p(3), &p(3))
    };
    let binomial = |n: usize, k: usize| (1..=k).fold(1usize, |acc, i| acc * (n + 1 - i) / i);
    for k in 1..=4usize {
        let numerator = (0..k as i64).fold(q(6, MU) * q(7, MU), |product, i| product * factor(i));
        for degree in 1..=2 {
            let parts = truncate(&numerator, &[2, 3], degree, 64).unwrap();
            assert_eq!(parts.len(), binomial(k + degree, degree), "k={k} d={degree}");
            assert_equal(
                &sum(&parts),
                &(graded_reference(&numerator, &[2, 3], degree).unwrap() * denominators()),
            );
            // The cap declines above it and admits exactly at it.
            let count = parts.len();
            assert_eq!(truncate(&numerator, &[2, 3], degree, count - 1), Err("parts"));
            assert_eq!(truncate(&numerator, &[2, 3], degree, count).map(|p| p.len()), Ok(count));
        }
    }
    // Vanishing graded factors are dropped from the count.
    let numerator = q(2, MU) * q(3, NU) * q(6, MU) * q(7, NU);
    assert_eq!(truncate(&numerator, &[2, 3], 1, 32), Ok(vec![]));
    let parts = truncate(&numerator, &[2, 3], 2, 32).unwrap();
    assert_eq!(parts.len(), 1);
}

#[test]
fn declines_keep_the_routed_form() {
    // 1/dot(Q2, Q2) inside N: not polynomial in the roots.
    let numerator = dot(&p(2), &p(2)).pow(-1) * q(6, MU) * q(7, MU);
    assert_eq!(truncate(&numerator, &[2], 0, 32), Err("non_polynomial"));
    // A root momentum under another function.
    let f = symbol!("uv_truncate_test::f");
    let numerator = function!(f, q(2, MU)) * q(6, MU);
    assert_eq!(truncate(&numerator, &[2], 0, 32), Err("non_polynomial"));
    // A root-free function is fine.
    let numerator = function!(f, q(6, MU)) * q(2, MU) * q(7, NU);
    assert_eq!(truncate(&numerator, &[2], 0, 32), Ok(vec![]));
    // A label outside the loop edges and the roots.
    let numerator = q(9, MU) * q(2, MU);
    assert_eq!(truncate(&numerator, &[2], 0, 32), Err("labels"));
    // A denominator inside a numerator factor.
    let inside = Atom::one() + den(6, 6, &p(2));
    let routed = &inside * q(2, MU) * q(6, MU) * denominators();
    let truncated = truncate_routed(&routed, &edges(&[2]), &loops(&[6, 7]), 0, 32).unwrap();
    assert_eq!(truncated.map(|t| t.parts), Err("nonsplit"));
    // Not a product.
    let routed = q(2, MU) * q(6, MU) * denominators() + q(6, MU) * q(6, MU) * denominators();
    let truncated = truncate_routed(&routed, &edges(&[2]), &loops(&[6, 7]), 0, 32).unwrap();
    assert_eq!(truncated.map(|t| t.parts), Err("not_product"));
    // No factor depends on a root.
    let numerator = q(6, MU) * q(7, MU);
    assert_eq!(truncate(&numerator, &[2], 0, 32), Err("no_root"));
}

/// The routed remainders of the support roots of the children of `dot`'s
/// automatic plan: (degree, root edges, untagged, tagged by `mode`).
type Remainders = Vec<(i32, Vec<EdgeIndex>, RStarIntegrand, RStarIntegrand)>;

fn child_remainders(dot: &str, vakint: &Vakint, mode: ChildTruncate) -> Result<(Graph, Remainders)> {
    let diagram = colour_factored(&load(dot)?)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let settings = settings();
    let mass = Atom::var(*REFERENCE_MASS_SQUARED).sqrt();
    let parent = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let cache = Arc::new(CountertermCache::default());
    let mut remainders = Vec::new();
    for call in &parent.counterterms {
        let support = &call.argument.support;
        if graph.compute_dod(support) < 0 {
            continue;
        }
        let Ok(plan) = RStarPlan::for_support(&graph, model, support, &settings) else {
            continue;
        };
        let plan = plan
            .with_uv_children(UvChildrenRoute::Rstar { min_loops: 1 })
            .with_counterterm_cache(Arc::clone(&cache));
        let children = plan.evaluate_uv_children_for(
            &graph,
            model,
            &mass,
            vakint,
            &settings,
            ChildMass::Reference,
        )?;
        let boundary = plan.root_boundary.as_ref().ok_or_else(|| eyre!("a support root"))?;
        let root_edges = graph
            .iter_edges_of(boundary)
            .map(|(_, edge, _)| edge)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        for term in &plan.ordinary_uv {
            if term.remainder.is_scaleless(&graph, plan.marked_edge) {
                continue;
            }
            let insertions = term
                .counterterms
                .iter()
                .map(|child| children[&child.argument.support].clone())
                .collect::<Vec<_>>();
            let build = || {
                RStarIntegrand::new_with_closure(
                    &graph,
                    model,
                    &term.remainder,
                    &insertions,
                    plan.marked_edge,
                    Atom::var(*REFERENCE_MASS_SQUARED),
                    None,
                )
            };
            let tagged = build()?.with_root_truncation_in(
                mode,
                usize::MAX,
                true,
                &root_edges,
                plan.taylor_order,
            )?;
            remainders.push((plan.taylor_order, root_edges.clone(), build()?, tagged));
        }
    }
    Ok((graph, remainders))
}

fn outer_polynomial(
    remainder: &RStarIntegrand,
    root_edges: &[EdgeIndex],
    degree: usize,
) -> Result<Atom> {
    let polynomial = remainder.outer_taylor(
        remainder.expression()?,
        root_edges,
        degree,
        *REFERENCE_MASS_SQUARED,
        &Atom::var(symbol!("rstar_mass_squared"; Scalar)),
        &mut IntegrationTimings::new(),
    )?;
    RStarIntegrand::simplify_summands(polynomial.as_view())
}

#[test]
fn truncated_child_remainders_keep_the_outer_taylor_polynomial() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let mut degrees = BTreeSet::new();
        let (mut tagged, mut untagged, mut vanishing, mut several) = (0, 0, 0, 0);
        for dot in [
            include_str!("testdata/l2-d0.dot"),
            include_str!("testdata/l2-d8.dot"),
            include_str!("testdata/l3-d166.dot"),
            include_str!("testdata/l3-d185.dot"),
        ] {
            let (_, remainders) = child_remainders(dot, &vakint, ChildTruncate::Graded)?;
            for (order, root_edges, full, truncated) in remainders {
                let Some(truncation) = &truncated.truncation else {
                    untagged += 1;
                    continue;
                };
                let degree = usize::try_from(order)?;
                assert_eq!(truncation.degree, degree);
                tagged += 1;
                degrees.insert(degree);
                vanishing += usize::from(truncation.parts.is_empty());
                several += usize::from(truncation.parts.len() > 1);
                let expected = outer_polynomial(&full, &root_edges, degree)?;
                let actual = outer_polynomial(&truncated, &root_edges, degree)?;
                let difference = (&expected - &actual).expand();
                if !difference.is_zero() {
                    let reduced = canonical_tensor(&difference)?;
                    ensure!(
                        reduced.together().is_zero(),
                        "degree {degree}: truncated {actual}\nfull {expected}\ndifference {reduced}"
                    );
                }
            }
        }
        eprintln!(
            "UV_TRUNCATE_TEST tagged={tagged} untagged={untagged} vanishing={vanishing} several_parts={several} degrees={degrees:?}"
        );
        ensure!(degrees.is_superset(&BTreeSet::from([0, 1, 2])), "{degrees:?}");
        ensure!(several > 0, "no graded remainder in several parts");
        Ok(())
    })
}

#[test]
fn truncated_remainders_guard_their_consumers_and_stream_their_parts() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let (_, remainders) =
            child_remainders(include_str!("testdata/l3-d166.dot"), &vakint, ChildTruncate::Graded)?;
        let mut checked = false;
        let mut streamed = false;
        for (order, root_edges, full, truncated) in remainders {
            let Some(truncation) = &truncated.truncation else {
                continue;
            };
            let degree = usize::try_from(order)?;
            // Another outer Taylor: fewer roots or a higher degree is an error.
            let expression = truncated.expression()?;
            let mut timing = IntegrationTimings::new();
            let mass = Atom::var(symbol!("rstar_mass_squared"; Scalar));
            let reference = *REFERENCE_MASS_SQUARED;
            assert!(
                truncated
                    .outer_taylor(expression, &root_edges[1..], degree, reference, &mass, &mut timing)
                    .is_err()
            );
            assert!(
                truncated
                    .outer_taylor(expression, &root_edges, degree + 1, reference, &mass, &mut timing)
                    .is_err()
            );
            let _ = truncated.outer_taylor(expression, &root_edges, degree, reference, &mass, &mut timing)?;
            // The consumers that need the whole remainder.
            assert!(truncated.taylor(0).is_err());
            assert!(truncated.taylor_contracted(0).is_err());
            assert!(
                truncated
                    .closed_moments(0, &mut Default::default())
                    .is_err()
            );
            assert!(truncated.physical_numerator().is_err());
            assert!(truncated.project_hard(expression).is_err());
            assert!(truncated.split_routed().is_err());
            // The truncation after a contraction is an error; off is inert.
            full.expression()?;
            assert!(
                full.with_root_truncation_in(ChildTruncate::Graded, 32, true, &root_edges, order)
                    .is_err()
            );
            checked = true;
            // The stream of several graded parts adds up to their contraction.
            if truncation.parts.len() > 1 && !streamed {
                let mut emitted = Vec::new();
                let applied = truncated.stream_parts(&truncation.parts, UvStream::On, &mut |batch| {
                    emitted.push(batch);
                    Ok(())
                })?;
                if applied {
                    let mut whole = Vec::new();
                    for part in &truncation.parts {
                        whole.push(RStarIntegrand::simplify_edge_basis(part)?);
                    }
                    assert_equal(&sum(&emitted), &sum(&whole));
                    streamed = true;
                }
            }
        }
        ensure!(checked && streamed, "checked={checked} streamed={streamed}");
        Ok(())
    })
}

#[test]
fn off_mode_is_inert() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let (_, remainders) =
            child_remainders(include_str!("testdata/l2-d8.dot"), &vakint, ChildTruncate::Off)?;
        ensure!(!remainders.is_empty());
        for (order, root_edges, full, off) in remainders {
            assert!(off.truncation.is_none());
            assert_eq!(off.contraction_inputs(), std::slice::from_ref(&full.routed));
            // Not a support root: untouched in every mode.
            let root = off.with_root_truncation_in(
                ChildTruncate::Graded,
                32,
                false,
                &root_edges,
                order,
            )?;
            assert!(root.truncation.is_none());
            assert_eq!(root.contracted_certified()?, full.contracted_certified()?);
            if let Ok(degree) = usize::try_from(order) {
                assert_eq!(
                    outer_polynomial(&root, &root_edges, degree)?,
                    outer_polynomial(&full, &root_edges, degree)?
                );
            }
        }
        Ok(())
    })
}
