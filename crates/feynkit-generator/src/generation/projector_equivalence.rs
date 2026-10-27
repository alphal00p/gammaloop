//! The automatic external-state projector is built in the finalized frame of every diagram and
//! memoised per external-state position. These tests rebuild it the way a requested projector is
//! written: in process external-leg labels, which generation resolves in each finalized frame.
//! Both constructions share the per-state factors, so the tests check the frame. They require
//! byte-identical diagrams, numerators and projectors for the comparison and the representative
//! frame of every diagram of representative amplitudes and cross sections.
//!
//! One generator serves all diagrams of a case, so a memo key that misses an input of the
//! projector hands a later diagram a stale projector. All diagrams of an amplitude share one key;
//! only cross sections, whose keys vary with the internal edge count, guard the key.

use super::cases::Case;
use super::*;

/// The projector of `graph` in process labels: every external state's carrier and half-edge
/// index name its external leg.
fn supplied_projector(generator: &Generator, graph: &GeneratedGraph) -> Atom {
    let mut projector = Atom::one();
    for (_, external, _) in graph.externals() {
        if let Some(factor) = generator
            .external_state_factor(
                external.particle,
                external.state,
                external.index,
                external.index,
            )
            .unwrap()
        {
            projector *= factor;
        }
    }
    projector
}

/// The diagram inputs of `case` as [`Request::generate`] prepares them, from a sequential
/// search, with the resolved options they were prepared with.
fn diagram_inputs(generator: &Generator, case: &Case) -> (NumeratorOptions, Vec<DiagramInput>) {
    let PreparedSearch {
        request,
        resolved,
        signatures,
        externals,
        plan,
        ..
    } = generator
        .search(&case.request(&generator.model), |_, _| {}, |_, _| {})
        .unwrap();
    let run = EngineRun::default();
    let inputs = plan
        .as_ref()
        .unwrap()
        .diagrams(&run)
        .unwrap()
        .filter_map(|diagram| {
            let multiplicity = externals.sewn_pair_orbit(&diagram)?;
            let graph = externals.graph(&diagram, &signatures);
            let symmetry = diagram.symmetry_factor().unwrap();
            generator
                .diagram_input(
                    &resolved,
                    graph,
                    symmetry,
                    multiplicity,
                    &request,
                    &case.grouping(),
                )
                .unwrap()
        })
        .collect();
    (NumeratorOptions::default(), inputs)
}

/// The atoms of a diagram: numerator, prefactor, projector, fragments and overall factor.
fn atoms(diagram: &FeynmanDiagram) -> Vec<Atom> {
    [
        diagram.numerator().unwrap(),
        diagram.numerator_prefactor(),
        diagram.projector(),
    ]
    .into_iter()
    .chain(diagram.vertices().map(|(_, vertex)| &vertex.numerator))
    .chain(diagram.edges().map(|(_, _, edge)| &edge.numerator))
    .map(|tensor| tensor.expression().clone())
    .chain([diagram.overall_factor().clone()])
    .collect()
}

/// Finalize every diagram input of `case` with the automatic projector and with its projector in
/// process labels, and require identical diagrams, down to the bytes of every atom. Returns the number
/// of memoised automatic projectors.
fn check(name: &str) -> usize {
    let case = Case::named(name);
    let generator = Generator::new(case.model());
    let (options, inputs) = diagram_inputs(&generator, case);
    assert!(!inputs.is_empty(), "{name} has no diagram");
    for (index, input) in inputs.into_iter().enumerate() {
        let context = format!("{name}_{index}");
        for topology in [input.comparison, input.representative] {
            let diagram = |options: &NumeratorOptions| {
                generator
                    .to_diagram(
                        Some(context.clone()),
                        case.kind,
                        topology.clone(),
                        input.signs,
                        input.multiplicity,
                        options,
                    )
                    .unwrap()
            };
            let projector = supplied_projector(&generator, &topology.graph);
            let automatic = diagram(&options);
            let supplied = diagram(&options.clone().with_projector(projector));
            assert_ne!(
                automatic.projector().expression(),
                &Atom::one(),
                "{context}"
            );
            assert_eq!(automatic.projector(), supplied.projector(), "{context}");
            assert_eq!(atoms(&automatic), atoms(&supplied), "{context}");
            let json = |diagram: &FeynmanDiagram| diagram.to_json().unwrap();
            assert_eq!(json(&automatic), json(&supplied), "{context}");
        }
    }
    generator.cache.projectors.read().unwrap().len()
}

#[test]
fn amplitude_projectors_match_supplied_process_label_projectors() {
    for name in [
        // Quarks and gluons, with interchangeable and with symmetrized legs.
        "qcd_gggg_L1",
        "qcd_uugg_L1",
        "qcd_gggg_sym_L1",
        // Fermions, with outgoing antiparticles.
        "qed_eemm_L2",
        "qed_ee4f_L0",
        // Massive vectors, and tops with gluons.
        "sm_eeww_L0",
        "sm_ggttg_L0",
        // Numerator grouping finalizes a comparison frame besides the representative.
        "qed_eeaa_identical_L1",
    ] {
        check(name);
    }
}

#[test]
fn cross_section_projectors_match_supplied_process_label_projectors() {
    let entries = [
        // Sewn fermion and gluon initial states.
        "xs_qed_eemm_L2",
        "xs_qcd_gguu_L1",
        "xs_qcd_gggg_L1",
        "xs_sm_uudd_L1",
        "xs_qcd_gguu_sym_L1",
        "xs_qed_eemm_identical_L2",
    ]
    .map(check);
    // Only a case with several memo keys can catch a key that misses an input.
    assert!(entries.iter().any(|&entries| entries > 1), "{entries:?}");
}
