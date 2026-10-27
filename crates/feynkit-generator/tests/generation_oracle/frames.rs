//! Frame-dependent numerator factors through grouping.
//!
//! The prefactor pairs internal and external momenta in provisional edge and vertex indices,
//! and every diagram requests its own loop-momentum basis by name, so grouping must move both
//! into each re-prepared comparison and representative frame. Grouping up to a sign reads the
//! comparison frames, whose fermion flows are not yet normalized; exact particles in the
//! canonical frame move the representatives away from their ungrouped frames. The cases cover
//! an amplitude with fermion-flow reversals and the two-pass left/right canonicalization of a
//! cross section. Their rows were recorded on the API line (f9dc5531), which tracked these
//! frames explicitly instead of deriving them from the composed permutations. The FeynKit port
//! re-recorded only the numerator digests (tensor numerators) and the cross section's content
//! (its loop momenta are routed from the outgoing sewn attachment); the diagram sets, frames,
//! names and ids are still the API line's.
//!
//! A requested projector replaces the external-state projector and moves through the same
//! relabelling as the prefactor; a bare dot product would leave the external spinor indices
//! open, so it is not sampled here.

use feynkit_generator::{Classes, GraphGroupingOptions, Grouping, Naming, NumeratorOptions, Run};
use feynkit_graph::momentum_symbol;
use feynkit_graph::symbols::vertex_index;
use feynkit_model::Model;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, FunctionBuilder};

use super::cases::Case;
use super::census::{Census, Scope, Summary};
use super::{Row, threads};

/// `Q(a, mu) Q(b, mu)`, with the dummy index `mu` owned by `vertex`.
fn dot(a: i64, b: i64, vertex: i64) -> Atom {
    let mu = Minkowski {}
        .new_rep(4)
        .to_symbolic([vertex_index().call((vertex, 9))]);
    let q = |edge: i64| {
        FunctionBuilder::new(momentum_symbol())
            .add_arg(edge)
            .add_arg(mu.clone())
            .finish()
    };
    q(a) * q(b)
}

/// Grouping that compares exact particles in the canonical frame.
fn exact_particles() -> Grouping {
    Grouping::Identical(GraphGroupingOptions {
        differentiate_particle_masses_only: false,
        ..GraphGroupingOptions::default()
    })
}

/// The case's QED model with massless u, d and s quarks, as in the Standard Model data the rows
/// were recorded with. Their equal loops are the identifications that grouping must find
/// through the comparison and representative frames; physical light-quark masses leave none.
fn model(case: &Case) -> Model {
    const LIGHT: [&str; 6] = ["u", "u~", "d", "d~", "s", "s~"];
    let light = |value: &serde_json::Value| LIGHT.contains(&value.as_str().unwrap());
    let mut definition = serde_json::to_value(case.model()).unwrap();
    for particle in definition["particles"].as_array_mut().unwrap() {
        if light(&particle["name"]) {
            particle["mass"] = serde_json::json!("ZERO");
        }
    }
    for propagator in definition["propagators"].as_array_mut().unwrap() {
        if light(&propagator["particle"]) {
            // No other symbol of these propagators contains the mass names.
            for field in ["numerator", "denominator"] {
                let massless = ["MU", "MD", "MS"].iter().fold(
                    propagator[field].as_str().unwrap().to_owned(),
                    |expression, mass| expression.replace(mass, "0"),
                );
                propagator[field] = serde_json::json!(massless);
            }
        }
    }
    Model::from_json(&definition.to_string()).unwrap()
}

fn census(case: &Case, grouping: &Grouping, threads: usize) -> Census {
    let model = model(case);
    let run = Run::default().threads(threads);
    let prepared = case
        .request(&model)
        .with_classes(Classes::Auto)
        .generate(model.into(), &run)
        .unwrap();
    // Edges 0 and 1 are internal, 4 and 6 are external legs.
    let named = NumeratorOptions::default().with_naming(Naming::prefix("D"));
    let mut options = named.clone().with_prefactor(dot(0, 4, 0) + dot(1, 6, 1));
    for diagram in prepared.build(&named, &run).unwrap().iter() {
        options = options.with_named_loop_momentum_basis(
            diagram.name().unwrap(),
            diagram.loop_momentum_basis().loop_edges.clone(),
        );
    }
    let grouped = prepared
        .build(&options, &run)
        .and_then(|built| built.group(grouping, &run))
        .unwrap_or_else(|error| panic!("{}: {error}", case.name));
    assert!(grouped.report().completed, "{} did not complete", case.name);
    Census::of(&grouped)
}

/// Require identical censuses at every thread count and compare their summary with `row`.
fn check(name: &str, grouping: Grouping, row: Row) {
    let case = Case::named(name);
    let mut counts = threads().into_iter();
    let reference = census(case, &grouping, counts.next().unwrap());
    for count in counts {
        reference.assert_same(
            &census(case, &grouping, count),
            Scope::Exact,
            &format!("{name} frames at {count} threads"),
        );
    }
    let summary = reference.summary();
    if std::env::var_os("FEYNKIT_ORACLE_RECORD").is_some() {
        println!("    frames {name}: {summary};");
    } else {
        let expected = Summary::from(&row);
        if summary != expected {
            reference.capture_failure(&reference, &format!("{name} recorded frame row"));
        }
        assert_eq!(
            summary, expected,
            "{name}: frame-dependent factors differ from their oracle row"
        );
    }
}

#[test]
fn comparison_frames_with_fermion_flow_reversals() {
    check(
        "qed_eeaa_sign_L1",
        Grouping::UpToSign(GraphGroupingOptions::default()),
        Row {
            diagrams: 60,
            topologies: 20,
            assignments: 66,
            inverse_symmetry: (60, 1),
            weights: (42, 1),
            symmetry_factors: &[(1, 60)],
            content: 0x9052838b19af8d40,
            numerators: 0xbbb6f06b55b1dd03,
            named: 0x2cfd98a195f641a8,
            ids: 0xaa53fdeb847997d1,
        },
    );
}

#[test]
fn representative_frames_with_fermion_flow_reversals() {
    check(
        "qed_eeaa_sign_L1",
        exact_particles(),
        Row {
            diagrams: 66,
            topologies: 20,
            assignments: 66,
            inverse_symmetry: (66, 1),
            weights: (42, 1),
            symmetry_factors: &[(1, 66)],
            content: 0x3e1fda8c899ed426,
            numerators: 0xbd1d483fc596a496,
            named: 0x0ca9d9be0107a4f0,
            ids: 0x79a2cf2343cb733c,
        },
    );
}

#[test]
fn representative_frames_with_left_right_symmetrization() {
    check(
        "xs_qed_aaee_lr_L2",
        exact_particles(),
        Row {
            diagrams: 130,
            topologies: 75,
            assignments: 3456,
            inverse_symmetry: (130, 1),
            weights: (176, 1),
            symmetry_factors: &[(1, 130)],
            content: 0x14725973c2aa6413,
            numerators: 0x87988a044c993bcf,
            named: 0x188d06a04da1a1e8,
            ids: 0x9fdbcd539a998b4a,
        },
    );
}
