#![cfg(feature = "feynkit-ingress")]
//! The physics row of benchmark v1 (`docs/benchmark/benchmark.typ`): the
//! massless-quark (L-1)-rung gluon self-energy at L = 2, 3, 4, the p analogue
//! of `tests/fixtures/qcd_ladder`. The numerator comes from the fixture's
//! native model rules (Feynman gauge), the external momentum is kept, and the
//! transverse projector closes the two gluon legs. The inputs are frozen in
//! `tests/fixtures/benchmark_p_physics.json`.
//!
//! Regenerate with
//! `VAKINT_P_PHYSICS=write cargo test --features feynkit-ingress --test
//! benchmark_p_tests -- --ignored`.

use ::vakint;
#[path = "fixtures/qcd_ladder/mod.rs"]
#[allow(dead_code)]
mod fixture;
#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "../examples/benchmark_p_integrals/results.rs"]
#[allow(dead_code)]
mod results;

use std::collections::BTreeMap;

use feynkit_graph::{
    DiagramEndpoint, DiagramHalfEdge, FeynmanDiagram, NumeratorInterface, symbols,
};
use feynkit_kinematics::Kinematics;
use idenso::{
    CookSettings, IndexTooling,
    color::ColorSimplifySettings,
    tensor::{AlgebraContraction, AlgebraSettings, ReductionStatus, SymbolicTensor},
};
use rustred::algebra::CoefficientContext;
use serde_json::json;
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AbstractIndex,
        representation::{LibrarySlot, Minkowski},
        slot::{IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    printer::PrintOptions,
    symbol,
};
use vakint::{Vakint, VakintExpression, VakintSettings, diagram_integral::DiagramIntegralOptions};

use fixture::Result;

const FROZEN: &str = "tests/fixtures/benchmark_p_physics.json";
const GENERATOR: &str = "crates/vakint/tests/benchmark_p_tests.rs";
const CONVENTIONS: &str = "One massless down-quark loop with L-1 gluon rungs, native qcd.json \
    rules, Feynman gauge, tr_Dirac(1) = 4, SU(3), G = 1. Each input is the graph numerator N with \
    Pi = delta_ab (g_munu q^2 - q_mu q_nu) N^{ab,munu} / [8 (d-1)], d = 4 - 2 eps, over the L-loop \
    graph with all propagators at power 1, so that Pi/(q^2)^2 is the transverse form factor. \
    `stripped_phase` is the constant phase of the numerator left out of the input.";

/// Every gluon-leg index of a Lorentz representation in `d` dimensions.
fn in_dimension(expression: &Atom, d: &Atom) -> Atom {
    expression.replace_map(|view, _, output| {
        if let AtomView::Fun(function) = view
            && function.get_symbol() == symbol!("spenso::mink")
        {
            let mut arguments = function.iter();
            if arguments.next().is_some() {
                **output = FunctionBuilder::new(function.get_symbol())
                    .add_arg(d)
                    .add_args(arguments)
                    .finish();
            }
        }
    })
}

/// The numerator of the graph with its open gluon legs, one (Lorentz, colour)
/// slot pair per leg.
fn amputated(
    diagram: &FeynmanDiagram,
    d: &Atom,
) -> Result<(Atom, Vec<Vec<LibrarySlot<AbstractIndex>>>)> {
    let raw = in_dimension(
        &(NumeratorInterface::export(diagram.numerator()).0
            * NumeratorInterface::export(diagram.numerator_prefactor()).0),
        d,
    );
    let cooking = CookSettings::reversible();
    let numerator = cooking
        .try_cook_indices(raw.as_view())
        .map_err(|error| format!("failed to encode the generated indices: {error:?}"))?;
    let ports = numerator
        .list_dangling::<AbstractIndex>()?
        .iter()
        .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let mut legs = BTreeMap::new();
    for (id, ends, edge) in diagram.edges() {
        let Some(external) = &edge.external else {
            continue;
        };
        let endpoint = if ends.source.is_some() {
            DiagramEndpoint::Source
        } else {
            DiagramEndpoint::Target
        };
        let hedge = diagram
            .half_edge_id(DiagramHalfEdge { edge: id, endpoint })
            .ok_or("an external edge has no half-edge")?;
        let index = symbols::hedge_index().call((hedge.0, 1));
        let mut slots = ports
            .iter()
            .filter(|slot| cooking.uncook(slot.aind.to_atom().as_view()) == index)
            .copied()
            .collect::<Vec<_>>();
        slots.sort();
        if slots.len() != 2 || legs.insert(external.index, slots).is_some() {
            return Err("each gluon leg needs one Lorentz and one colour port".into());
        }
    }
    if legs.len() != 2 || ports.len() != 4 {
        return Err("the amputated graph has unassigned external ports".into());
    }
    Ok((numerator, legs.into_values().collect()))
}

/// The Dirac and colour algebra of a fully contracted expression.
fn contract(expression: Atom) -> Result<Atom> {
    let reduced = SymbolicTensor::infer(expression)?.simplify_algebra(&AlgebraSettings {
        gamma: Some(Default::default()),
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::Dots,
        ..Default::default()
    })?;
    if reduced.reduction_status() != ReductionStatus::Complete {
        return Err(format!("incomplete contraction: {:?}", reduced.reduction_status()).into());
    }
    let expression = reduced.expression().clone();
    if !expression.list_dangling::<AbstractIndex>()?.is_empty() {
        return Err("the contraction left open indices".into());
    }
    Ok(expression)
}

/// The generation factors of the graph, without their provenance tags.
fn overall_factor(diagram: &FeynmanDiagram) -> Atom {
    let value = Atom::var(symbol!("benchmark_p_physics::factor_"));
    let mut result = diagram.overall_factor().to_owned();
    for head in [
        "AutG",
        "CouplingsMultiplicity",
        "InternalFermionLoopSign",
        "ExternalFermionOrderingSign",
        "AntiFermionSpinSumSign",
        "NumeratorIndependentSymmetryGrouping",
    ] {
        result = result
            .replace(symbol!(&format!("feynkit_generator_factor::{head}")).call(&value))
            .with(value.clone());
    }
    result.expand()
}

/// The product of dot products in `numerator` has exactly `degree` factors in
/// every term.
fn homogeneous(numerator: &Atom, degree: i64) -> bool {
    let scale = Atom::var(symbol!("benchmark_p_physics::scale"));
    let scaled = numerator.replace_map(|view, _, output| {
        if let AtomView::Fun(function) = view
            && function.get_symbol() == SPENSO_TAG.dot
        {
            **output = view.to_owned() * scale.clone();
        }
    });
    (scaled - numerator * scale.pow(degree)).expand().is_zero()
}

/// The ladder of `loop_count` loops as p-integral inputs.
#[derive(Debug)]
struct Lowered {
    /// The transverse projection, without its constant phase.
    transverse: Atom,
    /// `q_mu q_nu N^{ab,munu}` contracted with `delta_ab`, same phase removed.
    longitudinal: Atom,
    stripped_phase: &'static str,
}

fn lower(loop_count: usize) -> Result<Lowered> {
    Vakint::initialize_vakint_symbols();
    let evidence = tempfile::tempdir()?;
    let generated = fixture::dress_at_loop_count(evidence.path(), loop_count)?;
    let diagram = &generated.diagrams[0];
    let d = CoefficientContext::try_new(["d"])?
        .parameter("d")
        .ok_or("missing dimension")?
        .to_expression();
    let (numerator, legs) = amputated(diagram, &d)?;
    let routed = diagram.loop_momentum_basis().route_expression(&numerator);
    let family = diagram.propagator_family(&Kinematics::in_dimension(&d)?)?;
    let [q] = family.external_momenta() else {
        return Err("a self-energy has one independent external momentum".into());
    };
    let basis = diagram.loop_momentum_basis();
    let axis = basis
        .external_edges
        .iter()
        .position(|edge| !basis.dependent_externals.contains(edge))
        .ok_or("no independent external leg")?;

    let mut colour = Atom::one();
    let mut metric = Atom::one();
    let mut longitudinal = Atom::one();
    for slot in &legs[0] {
        let other = legs[1]
            .iter()
            .find(|other| other.rep_name() == slot.rep_name())
            .ok_or("external representations differ")?;
        if slot.rep_name() == (Minkowski {}).into() {
            let left = in_dimension(&slot.to_atom(), &d);
            let right = in_dimension(&other.to_atom(), &d);
            metric *= ETS.metric(&left, &right);
            longitudinal *= symbols::external_momentum().call((axis, left))
                * symbols::external_momentum().call((axis, right));
        } else {
            colour *= ETS.metric(slot.to_atom(), other.to_atom());
        }
    }
    let kinematics = Kinematics::in_dimension(&d)?
        .with_momenta(family.loop_momenta().iter().chain([q]).cloned())?;
    let square = kinematics.scalar_product(q, q)?;
    // G = 1, with the model's coupling GC_11 = i G.
    let couplings = symbol!("UFO::GC_11");
    let factor = overall_factor(diagram);
    let unit = symbolica::parse!("𝑖");
    let imaginary = |atom: &Atom| atom.to_canonical_string().contains('𝑖');
    let closed = |tensor: Atom| -> Result<Atom> {
        Ok((contract(tensor)? * &factor)
            .replace(Atom::var(couplings))
            .with(unit.clone())
            .expand())
    };
    let metric = closed(&routed * &colour * metric)?;
    let longitudinal = closed(routed * colour * longitudinal)?;

    // The numerator is a constant phase times a real polynomial.
    let (phase, stripped_phase) = if imaginary(&metric) {
        (-&unit, "i")
    } else {
        (Atom::one(), "1")
    };
    let real = |tensor: Atom| -> Result<Atom> {
        let tensor = (tensor * &phase).expand();
        if imaginary(&tensor) {
            return Err("the numerator is not a constant phase times a real polynomial".into());
        }
        Ok(tensor)
    };
    let metric = real(metric)?;
    let longitudinal = real(longitudinal)?;
    let degree = loop_count as i64;
    if !homogeneous(&metric, degree) || !homogeneous(&longitudinal, degree + 1) {
        return Err("the projected numerator fails momentum homogeneity".into());
    }

    let epsilon = Atom::var(vakint::vakint_symbol!(
        &VakintSettings::default().epsilon_symbol
    ));
    let dimension = Atom::num(4) - Atom::num(2) * epsilon;
    let ingress = |numerator: &Atom| -> Result<Atom> {
        let input: Atom = VakintExpression::from_diagram(
            diagram,
            &family,
            numerator,
            &DiagramIntegralOptions::default(),
        )?
        .into();
        Ok(input.replace(d.clone()).with(dimension.clone()))
    };
    let transverse = (metric * square - &longitudinal) / (Atom::num(8) * (&d - Atom::one()));
    Ok(Lowered {
        transverse: ingress(&transverse.expand())?,
        longitudinal: ingress(&longitudinal)?,
        stripped_phase,
    })
}

fn bare(input: &Atom) -> String {
    input.printer(PrintOptions::file_no_namespace()).to_string()
}

#[test]
#[ignore = "generator: writes the frozen fixture when VAKINT_P_PHYSICS=write"]
fn write_frozen_physics_inputs() -> Result<()> {
    let mut cells = Vec::new();
    for loops in [2, 3, 4] {
        let lowered = lower(loops)?;
        let input = bare(&lowered.transverse);
        eprintln!(
            "physics L{loops}: {} characters, {} terms",
            input.len(),
            lowered.transverse.expand().nterms()
        );
        cells.push(json!({
            "name": format!("physics_L{loops}"),
            "loops": loops,
            "rungs": loops - 1,
            "stripped_phase": lowered.stripped_phase,
            "input": input,
        }));
    }
    let document = json!({
        "schema": "vakint.p-physics-inputs.v1",
        "generator": GENERATOR,
        "conventions": CONVENTIONS,
        "cells": cells,
    });
    if std::env::var("VAKINT_P_PHYSICS").as_deref() == Ok("write") {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(FROZEN);
        std::fs::write(path, serde_json::to_string_pretty(&document)? + "\n")?;
    }
    Ok(())
}

mod physics {
    use super::*;
    use p_integrals::{forcer_evaluate, forcer_form, forcer_settings, parse};
    use serde_json::Value;
    use vakint::{LoopNormalizationFactor, RustRedMasterMode};

    fn frozen() -> Result<Vec<Value>> {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(FROZEN);
        let document: Value = serde_json::from_str(&std::fs::read_to_string(path)?)?;
        assert_eq!(document["schema"], "vakint.p-physics-inputs.v1");
        Ok(document["cells"].as_array().ok_or("no cells")?.clone())
    }

    fn retained(input: &Atom) -> Option<Atom> {
        let settings = forcer_settings(
            &forcer_form()?,
            RustRedMasterMode::Retain,
            LoopNormalizationFactor::GScheme,
            4,
        );
        Some(forcer_evaluate(input, &settings))
    }

    #[test]
    fn the_frozen_inputs_are_what_the_generator_writes() -> Result<()> {
        for cell in frozen()? {
            let loops = cell["loops"].as_u64().ok_or("no loops")? as usize;
            let lowered = lower(loops)?;
            let input = cell["input"].as_str().ok_or("no input")?;
            assert_eq!(input, bare(&lowered.transverse), "L = {loops} drifted");
            assert_eq!(parse(input), lowered.transverse, "L = {loops}");
            assert_eq!(cell["stripped_phase"], lowered.stripped_phase);
        }
        Ok(())
    }

    /// Without a rung the graph alone is gauge invariant: `q_mu q_nu` gives
    /// zero and the transverse part is the closed form of the quark bubble.
    #[test]
    fn the_one_loop_quark_bubble_is_transverse() -> Result<()> {
        let lowered = lower(1)?;
        let (Some(longitudinal), Some(transverse), Some(bubble)) = (
            retained(&lowered.longitudinal),
            retained(&lowered.transverse),
            retained(&parse(
                "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            )),
        ) else {
            return Ok(());
        };
        assert!(longitudinal.expand().is_zero(), "{longitudinal}");
        // N^{ab,munu} g_munu delta_ab = -4 tr[gamma^mu k gamma_mu (k - q)]
        // = 8 (2-d) q^2 up to scaleless terms, with colour tr(T^a T^a) = 4.
        let square = parse("dot(p(1),p(1))");
        let epsilon = parse("ε");
        let expected = (Atom::num(2) * &epsilon - Atom::num(2)) * square.pow(2)
            / (Atom::num(3) - Atom::num(2) * epsilon);
        assert!(
            (transverse - expected * bubble)
                .expand()
                .together()
                .is_zero(),
            "the quark bubble differs from its closed form"
        );
        Ok(())
    }

    #[test]
    fn the_two_loop_ladder_reduces_with_forcer() -> Result<()> {
        let cells = frozen()?;
        let Some(value) = retained(&parse(cells[0]["input"].as_str().ok_or("no input")?)) else {
            return Ok(());
        };
        assert!(!value.expand().is_zero());
        Ok(())
    }
}

mod results_file {
    use super::*;
    use serde_json::Value;

    use super::results::{self, Key, Source};

    fn measurement(cell: &str, repetition: usize) -> Value {
        json!({"cell": cell, "backend": "forcer", "arithmetic": "none", "threads": 1,
               "state": "cold", "repetition": repetition, "status": "ok", "evaluate_seconds": 0.5})
    }

    fn protocol(gammaloop: &str, rustred: &str) -> Value {
        json!({"gammaloop_sha": gammaloop, "rustred_sha": rustred, "build_profile": "ci-optim",
               "source_dirty": false, "p_library": {"identity": "p-library:scope=1-3"}})
    }

    fn source<'a>(gammaloop: &'a str, rustred: &'a str, build_profile: &'a str) -> Source<'a> {
        Source {
            gammaloop_sha: gammaloop,
            rustred_sha: rustred,
            build_profile,
            source_dirty: false,
            p_library: "p-library:scope=1-3",
        }
    }

    #[test]
    fn measurements_append_and_a_rerun_skips_what_this_source_recorded() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("p_results.json");
        results::begin_run(&path, "first", &protocol("aaa", "bbb"), true)?;
        results::append(&path, "first", measurement("L1_d0_r0", 0))?;
        results::begin_run(&path, "second", &protocol("aaa", "bbb"), false)?;
        results::append(&path, "second", measurement("L1_d0_r0", 1))?;
        results::begin_run(&path, "other", &protocol("ccc", "bbb"), false)?;
        results::append(&path, "other", measurement("L1_d0_r2", 0))?;

        let document = results::read(&path)?;
        assert_eq!(document["selected_run"], "first");
        let key = |cell: &str, repetition| Key {
            cell: cell.into(),
            backend: "forcer".into(),
            arithmetic: "none".into(),
            threads: 1,
            state: "cold".into(),
            repetition,
        };
        let recorded = results::recorded(&document, source("aaa", "bbb", "ci-optim"), None)?;
        assert_eq!(
            recorded,
            [key("L1_d0_r0", 0), key("L1_d0_r0", 1)]
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>()
        );
        let recorded_for = |source| results::recorded(&document, source, None);
        assert!(recorded_for(source("aaa", "bbb", "release"))?.is_empty());
        assert!(recorded_for(source("aaa", "zzz", "ci-optim"))?.is_empty());
        let dirty = Source {
            source_dirty: true,
            ..source("aaa", "bbb", "ci-optim")
        };
        assert!(recorded_for(dirty)?.is_empty());

        let before = std::fs::read(&path)?;
        assert!(results::append(&path, "first", measurement("L1_d0_r0", 0)).is_err());
        assert!(results::append(&path, "missing", measurement("L1_d0_r0", 2)).is_err());
        assert!(results::begin_run(&path, "first", &protocol("aaa", "bbb"), false).is_err());
        assert_eq!(std::fs::read(&path)?, before, "a refusal changed the file");
        Ok(())
    }

    #[test]
    fn the_p_library_is_part_of_the_source_and_a_child_must_have_used_it() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("p_results.json");
        results::begin_run(&path, "m2", &protocol("aaa", "bbb"), false)?;
        results::append(&path, "m2", measurement("L4_d0_r0", 0))?;
        let mut listed = protocol("aaa", "bbb");
        listed["p_library"]["identity"] = json!("p-library:scope=1-4");
        results::begin_run(&path, "m3", &listed, true)?;
        assert!(results::cite_run(&path, "missing").is_err());
        results::cite_run(&path, "m2")?;

        let document = results::read(&path)?;
        assert_eq!(document["cited_run"], "m2");
        let one_to_four = Source {
            p_library: "p-library:scope=1-4",
            ..source("aaa", "bbb", "ci-optim")
        };
        assert!(results::recorded(&document, one_to_four, None)?.is_empty());
        let built_in = results::recorded(&document, source("aaa", "bbb", "ci-optim"), None)?;
        assert_eq!(built_in.len(), 1);

        let library = json!({"origin": "external-data", "directory": "p-library-1-4"});
        let child = |program: &str, origin: &str| {
            json!({"cache": {"programs": [
                {"program": "1L-unit-mass-vacuum", "data_version": "v0.2.2"},
                {"program": program, "data_version": origin}]}})
        };
        assert!(
            results::check_p_program(&child("p-library-1-4", "external-data"), &library).is_ok()
        );
        assert!(results::check_p_program(&child("p-library-1-3", "built-in"), &library).is_err());
        assert!(results::check_p_program(&child("p-library-1-4", "built-in"), &library).is_err());
        assert!(results::check_p_program(&json!({"cache": {"programs": []}}), &library).is_err());
        Ok(())
    }

    #[test]
    fn the_arithmetic_is_part_of_what_a_rerun_skips() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("p_results.json");
        results::begin_run(&path, "run", &protocol("aaa", "bbb"), true)?;
        let rustred = |arithmetic: &str| {
            json!({"cell": "L3_d1_r4", "backend": "rustred", "arithmetic": arithmetic,
                   "threads": 16, "state": "cold", "repetition": 0, "status": "ok"})
        };
        results::append(&path, "run", rustred("exact"))?;
        results::append(&path, "run", rustred("modular_polynomial"))?;
        let mut without = rustred("exact");
        without.as_object_mut().unwrap().remove("arithmetic");
        assert!(results::append(&path, "run", without).is_err());

        let document = results::read(&path)?;
        let recorded = results::recorded(&document, source("aaa", "bbb", "ci-optim"), None)?;
        let arithmetics = recorded
            .iter()
            .map(|key| key.arithmetic.as_str())
            .collect::<Vec<_>>();
        assert_eq!(arithmetics, ["exact", "modular_polynomial"]);
        Ok(())
    }

    #[test]
    fn a_capped_record_is_final_only_up_to_its_cap() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("p_results.json");
        results::begin_run(&path, "run", &protocol("aaa", "bbb"), true)?;
        let capped = |cell, status| {
            json!({"cell": cell, "backend": "forcer", "arithmetic": "none", "threads": 16, "state": "warm",
                   "repetition": 0, "status": status, "cap_seconds": 100.0})
        };
        results::append(&path, "run", capped("L4_d0_r0", "capped"))?;
        results::append(&path, "run", capped("L4_d0_r2", "skipped"))?;

        let document = results::read(&path)?;
        let count = |cap| {
            results::recorded(&document, source("aaa", "bbb", "ci-optim"), cap)
                .map(|keys| keys.len())
        };
        assert_eq!(count(Some(50.0))?, 2);
        assert_eq!(count(Some(100.0))?, 2);
        assert_eq!(count(Some(1800.0))?, 0);
        assert_eq!(count(None)?, 0);
        Ok(())
    }

    #[test]
    fn exact_marks_must_come_from_the_measured_source() {
        let measured = source("aaa", "bbb", "release");
        let exact = |gammaloop, rustred| {
            json!({"schema": results::EXACT_SCHEMA, "gammaloop_sha": gammaloop,
                   "rustred_sha": rustred, "cells": []})
        };
        assert!(results::check_exact(&exact("aaa", "bbb"), measured).is_ok());
        assert!(results::check_exact(&exact("ccc", "bbb"), measured).is_err());
        assert!(results::check_exact(&exact("aaa", "ddd"), measured).is_err());
        let no_shas = json!({"schema": results::EXACT_SCHEMA, "cells": []});
        assert!(results::check_exact(&no_shas, measured).is_err());
        let stub = json!({"schema": results::EXACT_SCHEMA, "stub": true, "gammaloop_sha": null,
                          "rustred_sha": null, "cells": []});
        assert!(results::check_exact(&stub, measured).is_ok());
    }

    #[test]
    fn a_release_measurement_needs_a_clean_stamped_build_of_the_real_child() {
        let release = source("aaa", "bbb", "release");
        let dirty = Source {
            source_dirty: true,
            ..release
        };
        assert!(results::check_build(release, Some("aaa"), false).is_ok());
        assert!(results::check_build(release, None, false).is_err());
        assert!(results::check_build(release, Some("older"), false).is_err());
        assert!(results::check_build(release, Some("aaa"), true).is_err());
        assert!(results::check_build(dirty, Some("aaa-dirty"), false).is_err());

        let draft = source("aaa", "bbb", "ci-optim");
        let dirty_draft = Source {
            source_dirty: true,
            ..draft
        };
        assert!(results::check_build(draft, None, true).is_ok());
        assert!(results::check_build(draft, Some("aaa"), false).is_ok());
        assert!(results::check_build(draft, Some("older"), false).is_err());
        assert!(results::check_build(draft, Some("aaa-dirty"), false).is_err());
        assert!(results::check_build(dirty_draft, Some("aaa-dirty"), false).is_ok());
        assert!(results::check_build(dirty_draft, Some("aaa"), false).is_err());
    }

    #[test]
    fn a_reduction_without_a_math_span_has_no_seconds() {
        assert_eq!(results::reduction_seconds(true, 3, Some(1.5)), Some(1.5));
        assert_eq!(results::reduction_seconds(true, 0, Some(0.0)), None);
        assert_eq!(results::reduction_seconds(false, 3, Some(1.5)), None);
        assert_eq!(results::reduction_seconds(false, 3, None), None);
    }

    #[test]
    fn the_integral_count_is_one_per_top_on_the_grid_and_ignores_projector_powers() -> Result<()> {
        Vakint::initialize_vakint_symbols();
        let count = |fixture: &str, name: &str| -> Result<usize> {
            let cells: Value = serde_json::from_str(fixture)?;
            let cell = cells["cells"]
                .as_array()
                .and_then(|cells| cells.iter().find(|cell| cell["name"] == name))
                .ok_or("no such cell")?;
            let input = vakint::vakint_parse!(cell["input"].as_str().ok_or("no input")?)?;
            Ok(results::integral_count(&input))
        };
        let grid = include_str!("fixtures/benchmark_p_integrals.json");
        for (name, tops) in [
            ("L1_d0_r0", 1),
            ("L2_d1_r4", 1),
            ("L3_d2_r8", 3),
            ("L4_d0_r2", 11),
        ] {
            assert_eq!(count(grid, name)?, tops, "{name}");
        }
        let physics = include_str!("fixtures/benchmark_p_physics.json");
        for (name, integrals) in [("physics_L2", 11), ("physics_L3", 53), ("physics_L4", 293)] {
            assert_eq!(count(physics, name)?, integrals, "{name}");
        }
        Ok(())
    }

    #[test]
    fn a_file_of_another_schema_is_refused_and_kept() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("p_results.json");
        std::fs::write(&path, br#"{"schema": "something else", "runs": []}"#)?;
        assert!(results::begin_run(&path, "run", &protocol("a", "b"), false).is_err());
        std::fs::write(&path, b"{invalid")?;
        assert!(results::read(&path).is_err());
        assert_eq!(std::fs::read(&path)?, b"{invalid");
        Ok(())
    }
}
