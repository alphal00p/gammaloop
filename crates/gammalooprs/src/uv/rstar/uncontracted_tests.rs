//! U7a (`GAMMALOOP_UNCONTRACTED_NUMERATOR`, [`super::uncontracted`]): the
//! switch, the Δ-cache suffix, the record pieces (insertion split, the Δ
//! momentum-conservation step), the comparator, and the K1 records of real
//! beta diagrams against the VERIFIED fixtures of the beta project
//! (results/consolidation/wh-fixtures): the emitted record equals the
//! fixture's `input.txt` up to canonical order (atoms by expanded
//! difference), its value equals `value.txt` through `valid_through`, and
//! check mode (both paths, pole parts compared) passes and returns today's
//! value.
//!
//! The diagrams' colour branch is the fixture's: every vertex and propagator
//! factor that a fixture record states replaces the colour-factored rule
//! (`recursive_tests::colour_factored` keeps the first colour row, the beta
//! weights rows), so the in-process remainders are the beta's.
//!
//! Data: `testdata/u7a/` (k1-root, k1-root-rc and its r_c = 0 control, and the
//! insertion-split records); k1-quick is read in place from
//! `GAMMALOOP_U7A_FIXTURES_DIR` (default the beta's wh-fixtures directory),
//! skipped when absent.

use std::collections::{BTreeMap, BTreeSet};

use eyre::{Result, WrapErr, ensure, eyre};
use idenso::tensor::SymbolicTensor;
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
use vakint::{
    Vakint,
    p_integral::factorized::{
        Endpoint, FactorizedEdge, FactorizedKind, FactorizedNumerator, LineKind,
    },
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{GS, W_},
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    RStarPlan, UvChildClosure, UvChildrenRoute,
    boundary::{BoundaryMode, ClosedRootNormalization, with_boundary_mode},
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{
        self, Captured, Kind, Mode, Switch, agree, conserve, is_unknown_master_text,
        split_insertion, split_root_factor, with_capture, with_switch,
    },
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse(text: &str) -> Atom {
    symbolica::try_parse!(text).unwrap_or_else(|error| panic!("{error}: {text:.200}"))
}

fn minkowski() -> Atom {
    Minkowski {}.new_rep(GS.dim).to_symbolic([])
}

fn epsilon() -> symbolica::atom::Symbol {
    GS.dim_epsilon
}

/// The `value` line of a fixture value file.
fn fixture_value(text: &str) -> Atom {
    parse(
        text.lines()
            .find_map(|line| line.strip_prefix("value\t"))
            .expect("a value line"),
    )
}

/// The order `n` of a fixture's `valid_through eps^n` line (`-1` by default).
fn valid_through(text: &str) -> i64 {
    text.lines()
        .find_map(|line| line.strip_prefix("valid_through\t"))
        .and_then(|line| line.trim_start_matches("eps^").split_whitespace().next())
        .and_then(|order| order.parse().ok())
        .unwrap_or(-1)
}

/// `got == expected` through ε^order (the fixtures' identity rule).
fn equal_through(got: &Atom, expected: &Atom, order: i64) -> Result<bool> {
    let difference = (got - expected).expand();
    if difference.is_zero() {
        return Ok(true);
    }
    let truncated = difference.series(epsilon(), 0, order)?.to_atom().expand();
    Ok(truncated.is_zero() || truncated.together().is_zero())
}

/// The physical external rules of a two-point graph, as the beta builds
/// them: `(x0, Q(e, x___) -> ±Q(x0, x___))`.
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

/// The beta's closed metric-root normalisation (beta src/integrate.rs
/// `closed_root_scalar_normalization`): `c = transverse σ(prefactor) / (𝑖 p²
/// (d - 1))`.
fn beta_normalization(graph: &Graph) -> ClosedRootNormalization {
    let (external, external_routing) = external_rules(graph);
    let x0 = usize::from(external);
    let p2 = function!(
        SPENSO_TAG.dot,
        function!(GS.emr_mom, x0, minkowski()),
        function!(GS.emr_mom, x0, minkowski())
    );
    let factor = Atom::var(symbol!("beta_projection::transverse"))
        * graph
            .global_prefactor
            .num
            .replace_multiple(&external_routing)
        / (Atom::i() * p2)
        / (Atom::var(GS.dim) - 1);
    ClosedRootNormalization {
        external,
        external_routing,
        factor,
    }
}

/// The colour branch of a fixture: the colour-factored diagram with every
/// vertex and internal-edge factor that `records` state (mink(dim) back to
/// the DOT's mink(4)).
fn fixture_branch(
    dot: &str,
    records: &[&FactorizedNumerator],
) -> Result<feynkit_graph::FeynmanDiagram> {
    let diagram = colour_factored(&load(dot)?)?;
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for record in records {
        for vertex in &record.vertices {
            if !vertex.rule.starts_with("ct:") {
                vertices.insert(vertex.id, vertex.factor.map_mink_dim(Atom::num(4)));
            }
        }
        for edge in &record.edges {
            if !edge.is_external() {
                edges.insert(edge.id, edge.factor.map_mink_dim(Atom::num(4)));
            }
        }
    }
    let vertex_fragments = diagram
        .vertices()
        .filter_map(|(id, _)| vertices.get(&id.0).map(|factor| (id, factor.clone())))
        .map(|(id, factor)| Ok((id, SymbolicTensor::infer(factor)?)))
        .collect::<Result<Vec<_>>>()?;
    let edge_fragments = diagram
        .edges()
        .filter_map(|(id, _, _)| edges.get(&id.0).map(|factor| (id, factor.clone())))
        .map(|(id, factor)| Ok((id, SymbolicTensor::infer(factor)?)))
        .collect::<Result<Vec<_>>>()?;
    diagram
        .with_numerator_fragments(vertex_fragments, edge_fragments)
        .map_err(|error| eyre!("{error}"))
}

/// One diagram root (`projected_counterterm_closing_root_scalar`, the beta's
/// full-graph closing call) under `switch`, with the remainder-boundary path
/// on (production), and every attempt reported inside it.
fn root_run(
    diagram: &feynkit_graph::FeynmanDiagram,
    switch: Switch,
) -> Result<(Atom, Vec<Captured>, usize)> {
    let vakint = Vakint::new()?;
    let settings = settings();
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let closed_root = beta_normalization(&graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    let (value, captured) = with_switch(switch, || {
        with_boundary_mode(BoundaryMode::On, false, || {
            with_capture(|| {
                plan.projected_counterterm_closing_root_scalar(
                    &graph,
                    model,
                    &vakint,
                    &settings,
                    &callback,
                    &closed_root,
                )
            })
        })
    });
    Ok((value?, captured, usize::from(plan.marked_edge)))
}

/// The internal edge ids and the `ct:` rules of a record (its remainder).
fn record_key(record: &FactorizedNumerator) -> (Vec<usize>, Vec<String>) {
    let edges = record
        .edges
        .iter()
        .filter(|edge| !edge.is_external())
        .map(|edge| edge.id)
        .collect();
    let mut rules = record
        .vertices
        .iter()
        .filter(|vertex| vertex.rule.starts_with("ct:"))
        .map(|vertex| vertex.rule.clone())
        .collect::<Vec<_>>();
    rules.sort();
    (edges, rules)
}

fn atoms_agree(a: &Atom, b: &Atom) -> bool {
    a == b || (a - b).expand().is_zero() || agree(a, b).unwrap_or(false)
}

/// The comparator of U7a-DESIGN 11.2: the emitted record against a fixture
/// "up to canonical order" (the label ignored; structural lines equal; atoms
/// by expanded difference, then the tensor canonical form; an insertion by
/// its product coefficient × factor). The differences found, empty when they
/// agree.
fn record_differences(ours: &FactorizedNumerator, theirs: &FactorizedNumerator) -> Vec<String> {
    let mut out = Vec::new();
    let mut check = |ok: bool, what: String| {
        if !ok {
            out.push(what);
        }
    };
    match (&ours.kind, &theirs.kind) {
        (FactorizedKind::Closure(a), FactorizedKind::Closure(b)) => {
            check(
                a.edge == b.edge,
                format!("closure edge {} vs {}", a.edge, b.edge),
            );
            check(a.mass_squared == b.mass_squared, "closure mass".into());
            check(
                a.auxiliary_mass_squared == b.auxiliary_mass_squared,
                "auxiliary mass".into(),
            );
            check(a.power == b.power, "closure power".into());
            check(
                a.independent_externals == b.independent_externals,
                "externals".into(),
            );
            check(
                a.taylor == b.taylor,
                format!("taylor {:?} vs {:?}", a.taylor, b.taylor),
            );
            check(a.softs == b.softs, "softs".into());
            check(
                a.original_cycles == b.original_cycles,
                format!(
                    "original_cycles {:?} vs {:?}",
                    a.original_cycles, b.original_cycles
                ),
            );
            check(
                a.depth == b.depth,
                format!("depth {:?} vs {:?}", a.depth, b.depth),
            );
        }
        _ => check(false, "kind".into()),
    }
    check(
        ours.loop_edges == theirs.loop_edges,
        format!(
            "loop edges {:?} vs {:?}",
            ours.loop_edges, theirs.loop_edges
        ),
    );
    check(
        ours.external_edges == theirs.external_edges,
        format!(
            "external edges {:?} vs {:?}",
            ours.external_edges, theirs.external_edges
        ),
    );
    fn edges(r: &FactorizedNumerator) -> BTreeMap<usize, &FactorizedEdge> {
        r.edges.iter().map(|e| (e.id, e)).collect()
    }
    let (a, b) = (edges(ours), edges(theirs));
    check(
        a.keys().eq(b.keys()),
        format!(
            "edge ids {:?} vs {:?}",
            a.keys().collect::<Vec<_>>(),
            b.keys().collect::<Vec<_>>()
        ),
    );
    for (id, x) in &a {
        let Some(y) = b.get(id) else { continue };
        check(
            (x.source, x.target, x.source_hedge, x.target_hedge)
                == (y.source, y.target, y.source_hedge, y.target_hedge),
            format!("edge {id} ends"),
        );
        check(
            x.particle == y.particle,
            format!("edge {id} particle {} vs {}", x.particle, y.particle),
        );
        check(x.line == y.line, format!("edge {id} line"));
        check(x.power == y.power, format!("edge {id} power"));
        check(
            x.momentum == y.momentum,
            format!("edge {id} momentum {:?} vs {:?}", x.momentum, y.momentum),
        );
        check(x.mass_squared == y.mass_squared, format!("edge {id} mass"));
        check(
            atoms_agree(&x.factor, &y.factor),
            format!("edge {id} factor"),
        );
    }
    fn vertices(
        r: &FactorizedNumerator,
    ) -> BTreeMap<usize, &vakint::p_integral::factorized::FactorizedVertex> {
        r.vertices.iter().map(|v| (v.id, v)).collect()
    }
    let (a, b) = (vertices(ours), vertices(theirs));
    check(
        a.keys().eq(b.keys()),
        format!(
            "vertex ids {:?} vs {:?}",
            a.keys().collect::<Vec<_>>(),
            b.keys().collect::<Vec<_>>()
        ),
    );
    for (id, x) in &a {
        let Some(y) = b.get(id) else { continue };
        check(
            x.rule == y.rule,
            format!("vertex {id} rule {} vs {}", x.rule, y.rule),
        );
        check(
            atoms_agree(&x.factor, &y.factor),
            format!("vertex {id} factor"),
        );
    }
    let ports = |r: &FactorizedNumerator| {
        r.ports
            .iter()
            .map(|p| p.to_canonical_string())
            .collect::<BTreeSet<_>>()
    };
    check(ports(ours) == ports(theirs), "ports".into());
    match (&ours.projector, &theirs.projector) {
        (None, None) => {}
        (Some(x), Some(y)) => {
            check(
                x.p2_power == y.p2_power,
                format!("p2_power {} vs {}", x.p2_power, y.p2_power),
            );
            check(atoms_agree(&x.tensor, &y.tensor), "projector tensor".into());
        }
        _ => check(false, "projector presence".into()),
    }
    let (p, q) = (&ours.prefactor, &theirs.prefactor);
    check(
        (p.loop_sign, p.nf_power) == (q.loop_sign, q.nf_power),
        "prefactor integers".into(),
    );
    check(
        atoms_agree(&p.overall, &q.overall),
        "prefactor_overall".into(),
    );
    check(
        atoms_agree(&p.scalars, &q.scalars),
        format!("prefactor_scalars {} vs {}", p.scalars, q.scalars),
    );
    check(atoms_agree(&p.colour, &q.colour), "prefactor_colour".into());
    fn insertions(
        r: &FactorizedNumerator,
    ) -> BTreeMap<usize, &vakint::p_integral::factorized::Insertion> {
        r.insertions.iter().map(|i| (i.vertex, i)).collect()
    }
    let (a, b) = (insertions(ours), insertions(theirs));
    check(
        a.keys().eq(b.keys()),
        format!(
            "insertion vertices {:?} vs {:?}",
            a.keys().collect::<Vec<_>>(),
            b.keys().collect::<Vec<_>>()
        ),
    );
    for (vertex, x) in &a {
        let Some(y) = b.get(vertex) else { continue };
        check(
            x.hedges == y.hedges,
            format!("insertion {vertex} hedges {:?} vs {:?}", x.hedges, y.hedges),
        );
        check(x.cycles == y.cycles, format!("insertion {vertex} cycles"));
        check(
            atoms_agree(&(&x.coefficient * &x.factor), &(&y.coefficient * &y.factor)),
            format!("insertion {vertex} coefficient x factor"),
        );
    }
    out
}

/// One fixture record against the attempt captured for the same remainder:
/// the record (comparator), the value through `valid_through`, the check.
fn compare_fixture(what: &str, captured: &[Captured], input: &str, value_text: &str) -> Result<()> {
    let theirs = FactorizedNumerator::from_text(input).map_err(|e| eyre!("{what}: {e}"))?;
    let key = record_key(&theirs);
    let mut matches = Vec::new();
    for attempt in captured.iter().filter(|a| a.kind == Kind::Root) {
        let Some(text) = &attempt.text else { continue };
        let ours = FactorizedNumerator::from_text(text).map_err(|e| eyre!("{what}: {e}"))?;
        if record_key(&ours) == key {
            matches.push((attempt, ours));
        }
    }
    ensure!(
        matches.len() == 1,
        "{what}: {} captured records for the remainder {key:?} (outcomes {:?})",
        matches.len(),
        captured
            .iter()
            .map(|a| format!("{}:{}", a.id, a.outcome))
            .collect::<Vec<_>>()
    );
    let (attempt, ours) = &matches[0];
    // vakint's frozen writer and reader: a byte-identical round trip
    ensure!(
        ours.to_text().map_err(|e| eyre!("{e}"))? == *attempt.text.as_ref().unwrap(),
        "{what}: the record's text is not a fixed point of vakint's writer"
    );
    let differences = record_differences(ours, &theirs);
    ensure!(
        differences.is_empty(),
        "{what}: the emitted record differs from the fixture: {differences:#?}\n{}",
        attempt.text.as_ref().unwrap()
    );
    ensure!(
        attempt.outcome == "uncontracted",
        "{what}: outcome {}",
        attempt.outcome
    );
    let value = attempt
        .value
        .as_ref()
        .ok_or_else(|| eyre!("{what}: no value"))?;
    let expected = fixture_value(value_text);
    let order = valid_through(value_text);
    if let Some(contracted) = &attempt.contracted {
        ensure!(
            equal_through(contracted, &expected, order)?,
            "{what}: today's value differs from value.txt through eps^{order} (the fixture is not this branch)\ntoday  {contracted}\ntheirs {expected}"
        );
    }
    ensure!(
        equal_through(value, &expected, order)?,
        "{what}: the value differs from value.txt through eps^{order}\nours  {value}\ntheirs {expected}"
    );
    if let Some(equal) = attempt.equal {
        ensure!(equal, "{what}: check mode equal=false");
        eprintln!(
            "U7A_FIXTURE {what} id={} identical={} equal={equal}",
            attempt.id,
            attempt.identical.unwrap_or(false)
        );
    }
    Ok(())
}

fn fixtures_dir() -> std::path::PathBuf {
    std::env::var("GAMMALOOP_U7A_FIXTURES_DIR")
        .unwrap_or_else(|_| "/common/ben/beta/results/consolidation/wh-fixtures".into())
        .into()
}

fn check_switch() -> Switch {
    Switch::new(Mode::Check, &[Kind::Root, Kind::Child])
}

// ---------------------------------------------------------------------------
// The switch and the cache key
// ---------------------------------------------------------------------------

#[test]
fn switch_parses_every_form() {
    for off in [None, Some(""), Some("0"), Some("off"), Some(" off ")] {
        assert_eq!(Switch::parse(off), Ok(Switch::OFF), "{off:?}");
    }
    let on = Switch::parse(Some("on")).unwrap();
    assert_eq!(on.mode, Mode::On);
    assert!(on.enabled(Kind::Root) && on.enabled(Kind::Child));
    assert!(!on.enabled(Kind::IrHard) && !on.enabled(Kind::Quotient));
    assert_eq!(Switch::parse(Some("1")), Ok(on));
    let check = Switch::parse(Some("check")).unwrap();
    assert_eq!(check.mode, Mode::Check);
    assert_eq!(check.enabled_kinds(), vec![Kind::Root, Kind::Child]);
    assert_eq!(
        Switch::parse(Some("check:a,b")),
        Ok(check),
        "letter aliases"
    );
    let root = Switch::parse(Some("on:root, root")).unwrap();
    assert_eq!(root.enabled_kinds(), vec![Kind::Root], "duplicates ignored");
    let inert = Switch::parse(Some("on:ir-hard,quotient,surrogate,insertion")).unwrap();
    assert!(inert.enabled_kinds().is_empty());
    assert_eq!(
        inert.inert_kinds(),
        vec![
            Kind::IrHard,
            Kind::Quotient,
            Kind::Surrogate,
            Kind::Insertion
        ]
    );
    assert_eq!(
        inert.cache_key_suffix(),
        None,
        "nothing enabled, today's key"
    );
    let mixed = Switch::parse(Some("check:c,root")).unwrap();
    assert_eq!(mixed.enabled_kinds(), vec![Kind::Root]);
    assert_eq!(
        mixed.line(),
        "UNCONTRACTED_SWITCH mode=check kinds=root inert=ir-hard"
    );
    assert_eq!(
        on.line(),
        "UNCONTRACTED_SWITCH mode=on kinds=root,child inert=none"
    );
    for bad in [
        "yes",
        "on:",
        "on:bogus",
        "on:root,",
        "off:root",
        "check:root,,child",
        "checks",
    ] {
        let error = Switch::parse(Some(bad)).expect_err(bad);
        assert!(error.contains(uncontracted::VARIABLE), "{error}");
    }
}

#[test]
fn cache_suffix_only_when_enabled() {
    let key = || {
        RStarPlan::rstar_child_key(
            "settings".into(),
            UvChildrenRoute::Rstar { min_loops: 1 },
            UvChildClosure::Off,
            Some("#uv-child-truncate=graded@1".into()),
        )
    };
    let off = with_switch(Switch::OFF, key);
    assert!(!off.contains("#uncontracted"), "{off}");
    assert!(off.ends_with("#uv-child-truncate=graded@1"), "{off}");
    let on = with_switch(Switch::new(Mode::On, &[Kind::Root, Kind::Child]), key);
    assert_eq!(on, format!("{off}#uncontracted=root,child"));
    let check = with_switch(Switch::new(Mode::Check, &[Kind::Child, Kind::Root]), key);
    assert_eq!(check, on, "the mode is not part of the key");
    let inert = with_switch(Switch::new(Mode::On, &[Kind::Quotient]), key);
    assert_eq!(inert, off);
    assert_eq!(
        with_switch(Switch::OFF, uncontracted::cache_key_suffix),
        None
    );
}

#[test]
fn unknown_master_text_is_recognised() {
    let unknown = vakint::UnknownMasterCoefficients {
        backend: "rustred".into(),
        requested_order: Some(0),
        coefficients: vec![vakint::UnknownMasterCoefficient {
            master: "p5_m1".into(),
            order: 1,
            first_affected_result_order: 0,
        }],
    };
    assert!(is_unknown_master_text(&unknown.to_string()), "{unknown}");
    // as vakint's factorised entry wraps it
    let wrapped = vakint::p_integral::factorized::FactorizedError::Evaluation {
        detail: unknown.to_string(),
    };
    assert!(is_unknown_master_text(&wrapped.to_string()));
    let missing = vakint::MissingMasterOrders {
        backend: "rustred".into(),
        requested_order: Some(0),
        masters: Vec::new(),
    };
    assert!(is_unknown_master_text(&missing.to_string()), "{missing}");
    assert!(!is_unknown_master_text(
        "closure: degree 3 exceeds the maximum Taylor degree"
    ));
}

// ---------------------------------------------------------------------------
// Record pieces
// ---------------------------------------------------------------------------

#[test]
fn split_matches_the_converter() {
    test_initialise().unwrap();
    for (what, text) in [
        (
            "k1-quick L4-42 rr2",
            include_str!("testdata/u7a/k1-quick-L4-42-rr2.input.txt"),
        ),
        (
            "k2-child-ins",
            include_str!("testdata/u7a/k2-child-ins.input.txt"),
        ),
    ] {
        let record = FactorizedNumerator::from_text(text).unwrap_or_else(|e| panic!("{what}: {e}"));
        assert!(!record.insertions.is_empty(), "{what}");
        for insertion in &record.insertions {
            let delta = (&insertion.coefficient * &insertion.factor).expand();
            let (coefficient, factor) = split_insertion(&delta);
            assert!(
                (&coefficient * &factor - &delta).expand().is_zero(),
                "{what}: coefficient x factor is not the counterterm"
            );
            let ratio = &coefficient / &insertion.coefficient;
            let ratio = ratio.to_canonical_string();
            assert!(
                ratio
                    .trim_start_matches('-')
                    .bytes()
                    .all(|b| b.is_ascii_digit() || b == b'/'),
                "{what}: coefficient {coefficient} is not a rational multiple of the converter's {}",
                insertion.coefficient
            );
            assert_ne!(
                coefficient,
                Atom::num(1),
                "{what}: the converter's split is non-trivial"
            );
        }
    }
    // a counterterm whose terms' scalars are not rational multiples: 1 x Δ
    let g = parse(
        "spenso::g(spenso::mink(gammalooprs::dim,gammalooprs::hedge(1,1)),spenso::mink(gammalooprs::dim,gammalooprs::hedge(2,1)))",
    );
    let delta = (Atom::var(symbol!("u7a_test_a")) + Atom::var(symbol!("u7a_test_b")) * &g).expand();
    let (coefficient, factor) = split_insertion(&delta);
    assert_eq!(coefficient, Atom::num(1));
    assert_eq!(factor, delta);
}

fn edge(id: usize, source: Endpoint, target: Endpoint, momentum: &[i64]) -> FactorizedEdge {
    FactorizedEdge {
        id,
        source,
        target,
        source_hedge: None,
        target_hedge: None,
        particle: "Q".into(),
        line: LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum: momentum.to_vec(),
        factor: Atom::num(1),
    }
}

#[test]
fn conservation_removes_zero_mass_logs() {
    test_initialise().unwrap();
    let q = |e: usize, h: usize| {
        function!(
            GS.emr_mom,
            e,
            function!(
                symbol!("spenso::mink"),
                GS.dim,
                function!(GS.hedgeaind, h, 1)
            )
        )
    };
    let mass = Atom::var(symbol!("rstar_mass_squared"; Scalar));
    let log = function!(symbolica::atom::Symbol::LOG, &mass);
    let g = function!(
        symbol!("spenso::g"),
        function!(
            symbol!("spenso::mink"),
            GS.dim,
            function!(GS.hedgeaind, 7, 1)
        ),
        function!(
            symbol!("spenso::mink"),
            GS.dim,
            function!(GS.hedgeaind, 8, 1)
        )
    );
    // node 5: edge 1 in (k1), edge 2 out (k1 - p1), edge 3 out (p1)
    let v = Endpoint::Vertex(5);
    let edges = [
        edge(1, Endpoint::Vertex(0), v, &[1, 0]),
        edge(2, v, Endpoint::Vertex(1), &[1, -1]),
        edge(3, v, Endpoint::External, &[0, 1]),
    ];
    let delta = (Atom::num(21) / 2048
        * &log
        * (q(1, 9) - q(2, 9) - q(3, 9))
        * Atom::var(symbol!("u7a_test_c"))
        + &g * Atom::var(symbol!("u7a_test_d")))
    .expand();
    assert!(delta.contains_symbol(mass.get_symbol().unwrap()));
    let rewritten = conserve(&delta, 5, &edges).expect("a conserving node");
    assert!(
        !rewritten.contains_symbol(mass.get_symbol().unwrap()),
        "the zero-mass log survives: {rewritten}"
    );
    assert!(
        (&rewritten - &g * Atom::var(symbol!("u7a_test_d")))
            .expand()
            .is_zero()
    );
    // a K2 crown vertex where independent legs meet: unchanged (None)
    let independent = [
        edge(1, Endpoint::External, v, &[0, 1, 0]),
        edge(2, Endpoint::External, v, &[0, 0, 1]),
        edge(3, v, Endpoint::Vertex(1), &[1, 0, 0]),
    ];
    assert!(conserve(&delta, 5, &independent).is_none());
}

#[test]
fn root_factor_split() {
    test_initialise().unwrap();
    let p2 = function!(
        SPENSO_TAG.dot,
        function!(GS.emr_mom, 0, minkowski()),
        function!(GS.emr_mom, 0, minkowski())
    );
    let t = Atom::var(symbol!("beta_projection::transverse"));
    let c = &t / (Atom::i() * &p2) / (Atom::var(GS.dim) - 1);
    let (k, s) = split_root_factor(&c, 0).unwrap();
    assert_eq!(k, -1);
    let theirs = parse("(-1+gammalooprs::{}::dim)^(-1)*(-𝑖)*beta_projection::{}::transverse");
    assert!((&s - &theirs).expand().is_zero(), "{s}");
    let (k, s) = split_root_factor(&Atom::num(1), 0).unwrap();
    assert_eq!((k, s), (0, Atom::num(1)));
    let momentum = function!(GS.emr_mom, 0, minkowski());
    assert!(split_root_factor(&(&t * momentum), 0).is_err());
}

#[test]
fn comparator_catches_mutations() {
    test_initialise().unwrap();
    for (what, text) in [
        ("k1-root", include_str!("testdata/u7a/k1-root.input.txt")),
        (
            "k1-quick L4-42 rr2",
            include_str!("testdata/u7a/k1-quick-L4-42-rr2.input.txt"),
        ),
    ] {
        let record = FactorizedNumerator::from_text(text).unwrap();
        assert!(record_differences(&record, &record).is_empty(), "{what}");
        // a re-labelled copy still agrees (the label is ignored)
        let mut relabelled = record.clone();
        relabelled.label = "other".into();
        assert!(record_differences(&relabelled, &record).is_empty());
        let mutations: Vec<(&str, Box<dyn Fn(&mut FactorizedNumerator)>)> = vec![
            (
                "a wrong vector",
                Box::new(|r: &mut FactorizedNumerator| {
                    let e = r.edges.iter_mut().find(|e| !e.is_external()).unwrap();
                    e.momentum[0] += 1;
                }),
            ),
            (
                "a dropped (-i) in the scalars",
                Box::new(|r: &mut FactorizedNumerator| {
                    r.prefactor.scalars = &r.prefactor.scalars * Atom::i();
                }),
            ),
            (
                "a wrong p2_power",
                Box::new(|r: &mut FactorizedNumerator| {
                    r.projector.as_mut().unwrap().p2_power += 1;
                }),
            ),
            (
                "a wrong vertex factor",
                Box::new(|r: &mut FactorizedNumerator| {
                    let v = r
                        .vertices
                        .iter_mut()
                        .find(|v| !v.rule.starts_with("ct:"))
                        .unwrap();
                    v.factor = &v.factor * Atom::num(2);
                }),
            ),
        ];
        for (mutation, apply) in &mutations {
            let mut mutated = record.clone();
            apply(&mut mutated);
            assert!(
                !record_differences(&mutated, &record).is_empty(),
                "{what}: {mutation} is not caught"
            );
        }
        if !record.insertions.is_empty() {
            let mut removed = record.clone();
            removed.insertions.pop();
            assert!(
                !record_differences(&removed, &record).is_empty(),
                "{what}: a removed insertion"
            );
            let mut hedges = record.clone();
            hedges.insertions[0].hedges.pop();
            assert!(
                !record_differences(&hedges, &record).is_empty(),
                "{what}: wrong hedges"
            );
            let mut factor = record.clone();
            factor.insertions[0].coefficient = &factor.insertions[0].coefficient * Atom::num(3);
            assert!(
                !record_differences(&factor, &record).is_empty(),
                "{what}: wrong coefficient"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// K1 records of real diagrams
// ---------------------------------------------------------------------------

/// One closed root with a single-record fixture: emitted record == the
/// fixture's input.txt, value == value.txt, check equal, and the check
/// mode's diagram result == off's.
fn single_root_fixture(what: &str, dot: &str, input: &str, value_text: &str) -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let theirs = FactorizedNumerator::from_text(input).map_err(|e| eyre!("{e}"))?;
        let diagram = fixture_branch(dot, &[&theirs])?;
        let (off, nothing, marked) = root_run(&diagram, Switch::OFF)?;
        ensure!(nothing.is_empty(), "{what}: mode off reports nothing");
        let FactorizedKind::Closure(spec) = &theirs.kind else {
            return Err(eyre!("{what}: not a closure record"));
        };
        ensure!(
            marked == spec.edge,
            "{what}: marked edge {marked}, the fixture's closure edge {}",
            spec.edge
        );
        let (checked, captured, _) = root_run(&diagram, check_switch())?;
        ensure!(
            checked == off,
            "{what}: check mode returns today's value\noff   {off}\ncheck {checked}"
        );
        for attempt in captured.iter().filter(|a| a.kind == Kind::Root) {
            eprintln!(
                "U7A_ROOT {what} graph={} id={} outcome={} identical={:?} equal={:?}",
                attempt.graph, attempt.id, attempt.outcome, attempt.identical, attempt.equal
            );
            ensure!(
                attempt.equal != Some(false),
                "{what}: id {} check equal=false",
                attempt.id
            );
        }
        compare_fixture(what, &captured, input, value_text)?;
        // mode on: the diagram's poles equal off's
        let (on, _, _) = root_run(&diagram, Switch::new(Mode::On, &[Kind::Root]))?;
        ensure!(
            atoms_agree(&on, &off),
            "{what}: mode on differs from off\noff {off}\non  {on}"
        );
        Ok(())
    })
}

#[test]
fn k1_root_fixture() -> Result<()> {
    single_root_fixture(
        "k1-root (L2 22)",
        include_str!("testdata/u7a/l2-d22.dot"),
        include_str!("testdata/u7a/k1-root.input.txt"),
        include_str!("testdata/u7a/k1-root.value.txt"),
    )
}

#[test]
fn k1_root_rc_fixture() -> Result<()> {
    if std::env::var("GAMMALOOP_ROOT_ROUTING").is_ok_and(|v| v.trim() == "default") {
        eprintln!("SKIP: k1-root-rc needs the num-est root routing");
        return Ok(());
    }
    single_root_fixture(
        "k1-root-rc (L2 21, r_c = +p(1))",
        include_str!("testdata/u7a/l2-d21.dot"),
        include_str!("testdata/u7a/k1-root-rc.input.txt"),
        include_str!("testdata/u7a/k1-root-rc.value.txt"),
    )
}

/// The r_c = 0 control of k1-root-rc: the default root routing, read once
/// per process, so run alone with `GAMMALOOP_ROOT_ROUTING=default`.
#[test]
#[ignore = "needs GAMMALOOP_ROOT_ROUTING=default in its own process"]
fn k1_root_rc0_fixture() -> Result<()> {
    if !std::env::var("GAMMALOOP_ROOT_ROUTING").is_ok_and(|v| v.trim() == "default") {
        eprintln!("SKIP: k1-root-rc0 needs GAMMALOOP_ROOT_ROUTING=default");
        return Ok(());
    }
    single_root_fixture(
        "k1-root-rc0 (L2 21, r_c = 0)",
        include_str!("testdata/u7a/l2-d21.dot"),
        include_str!("testdata/u7a/k1-root-rc0.input.txt"),
        include_str!("testdata/u7a/k1-root-rc0.value.txt"),
    )
}

/// k1-quick L4-100: both root evaluations (colour branches b0 / b1), seven
/// boundary remainders each (`rr1`, bare, states every vertex factor of its
/// branch): every fixture record emitted and equal, every value equal, every
/// check equal; the one-loop remainders (no fixture) checked equal too.
#[test]
fn k1_quick_l4_100_fixtures() -> Result<()> {
    let root = fixtures_dir().join("k1-quick/L4-100");
    if !root.is_dir() {
        eprintln!("SKIP: no k1-quick data under {}", root.display());
        return Ok(());
    }
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let read = |path: std::path::PathBuf| {
            std::fs::read_to_string(&path).wrap_err_with(|| format!("{}", path.display()))
        };
        let dot = read(root.join("b0-rr1/raw/diagram-000100.dot"))?;
        let mut records = 0;
        for branch in ["b0", "b1"] {
            let bare = read(root.join(format!("{branch}-rr1/input.txt")))?;
            let bare = FactorizedNumerator::from_text(&bare).map_err(|e| eyre!("{e}"))?;
            let diagram = fixture_branch(&dot, &[&bare])?;
            let started = std::time::Instant::now();
            let (_, captured, _) = root_run(&diagram, check_switch())?;
            eprintln!(
                "U7A_L4_100 {branch} seconds={:.1} attempts={}",
                started.elapsed().as_secs_f64(),
                captured.len()
            );
            for attempt in captured.iter().filter(|a| a.kind == Kind::Root) {
                eprintln!(
                    "U7A_L4_100 {branch} id={} outcome={} identical={:?} equal={:?}",
                    attempt.id, attempt.outcome, attempt.identical, attempt.equal
                );
                ensure!(attempt.equal != Some(false), "{branch} id {}", attempt.id);
            }
            for rr in 1..=7 {
                let dir = root.join(format!("{branch}-rr{rr}"));
                compare_fixture(
                    &format!("L4-100 {branch}-rr{rr}"),
                    &captured,
                    &read(dir.join("input.txt"))?,
                    &read(dir.join("value.txt"))?,
                )?;
                records += 1;
            }
        }
        ensure!(records == 14);
        Ok(())
    })
}

/// Check mode == off on the in-repo production diagrams (colour factored,
/// beta-shaped closed root): every root remainder is attempted, every
/// check is equal, and the diagram's result is today's.
#[test]
fn k1_check_equals_off() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        for (name, dot) in [
            ("L2/0", include_str!("testdata/l2-d0.dot")),
            ("L2/8", include_str!("testdata/l2-d8.dot")),
            ("L2/12", include_str!("testdata/l2-d12.dot")),
            ("L3/185", include_str!("testdata/l3-d185.dot")),
        ] {
            let diagram = colour_factored(&load(dot)?)?;
            let (off, _, _) = root_run(&diagram, Switch::OFF)?;
            let (checked, captured, _) = root_run(&diagram, check_switch())?;
            ensure!(checked == off, "{name}: check mode returns today's value");
            let roots = captured
                .iter()
                .filter(|a| a.kind == Kind::Root)
                .collect::<Vec<_>>();
            ensure!(!roots.is_empty(), "{name}: no root record attempted");
            for attempt in &roots {
                eprintln!(
                    "U7A_CHECK {name} id={} outcome={} identical={:?} equal={:?}",
                    attempt.id, attempt.outcome, attempt.identical, attempt.equal
                );
                ensure!(
                    attempt.equal != Some(false),
                    "{name}: id {} equal=false",
                    attempt.id
                );
            }
            ensure!(
                roots.iter().any(|a| a.outcome == "uncontracted"),
                "{name}: no record was evaluated by vakint"
            );
            for child in captured.iter().filter(|a| a.kind == Kind::Child) {
                ensure!(
                    child.outcome == "declined:child_not_built",
                    "{name}: child outcome {}",
                    child.outcome
                );
            }
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// K2 (after the skeleton)
// ---------------------------------------------------------------------------

/// TODO(U7a K2 emitter): the isotropic closed child {2,3,5} of L2 22 against
/// wh-fixtures/k2-child (record and value).
#[test]
#[ignore = "TODO: the K2 emitter is not in the skeleton"]
fn k2_child_fixture() {
    eprintln!("TODO: k2-child (L2 22, isotropic child {{2,3,5}})");
}

/// TODO(U7a K2 emitter): the child {5,6,7,8} of L3 36 with {7,8} contracted
/// against wh-fixtures/k2-child-ins.
#[test]
#[ignore = "TODO: the K2 emitter is not in the skeleton"]
fn k2_child_ins_fixture() {
    eprintln!("TODO: k2-child-ins (L3 36, child {{5,6,7,8}}, insertion {{7,8}})");
}

/// TODO(U7a K2 emitter): an open child (L2 22, degree 1), a Dirac vertex
/// (L2 28), a Dirac self-energy with its open nu (L2 26) and a closed child
/// through the child boundary (L3 352): values against the raw b_uv_child
/// values through eps^-1, and check equal.
#[test]
#[ignore = "TODO: the K2 emitter is not in the skeleton"]
fn k2_shapes_values() {
    eprintln!("TODO: K2 open / Dirac vertex / Dirac self-energy / child boundary");
}
