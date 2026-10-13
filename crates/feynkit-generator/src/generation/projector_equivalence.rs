//! The automatic external-state projector is built in the finalized frame of every diagram and
//! memoised per external-state position. These tests rebuild it the way a requested projector is
//! built: in the provisional frame of the generated graph, then finalized with the diagram. Both
//! constructions share the per-state factors, so the tests check the frame. They require
//! byte-identical diagrams, numerators and projectors for the comparison and the representative
//! frame of every diagram of representative amplitudes and cross sections.
//!
//! One generator serves all diagrams of a case, so a memo key that misses an input of the
//! projector hands a later diagram a stale projector. All diagrams of an amplitude share one key;
//! only cross sections, whose keys vary with the internal edge count, guard the key.

use super::{cases::Case, *};

/// The projector in the provisional frame of `graph`: every external state on its provisional
/// edge, indexed at the half-edge where it attaches (`2 * edge`, plus one at the sink).
fn provisional_projector(generator: &Generator, graph: &GeneratedGraph) -> Atom {
    let positions = graph
        .externals()
        .map(|(edge, external, _)| ExternalStatePosition {
            external,
            edge,
            hedge: 2 * edge + usize::from(graph.edge(edge).source.external().is_some()),
        })
        .collect::<Vec<_>>();
    generator.external_state_projector(&positions).unwrap()
}

/// The diagram inputs of `case` as [`Generator::generate`] prepares them, from a sequential
/// search, with the resolved options they were prepared with.
fn diagram_inputs(generator: &Generator, case: &Case) -> (GenerationOptions, Vec<DiagramInput>) {
    let PreparedSearch {
        options,
        resolved,
        signatures,
        externals,
        generator: search,
    } = generator
        .search(
            &case.process(),
            &case.options(&generator.model, 1),
            case.kind,
            |_, _| {},
        )
        .unwrap();
    let inputs = search
        .iter()
        .unwrap()
        .filter_map(|diagram| {
            let multiplicity = externals.sewn_pair_orbit(&diagram)?;
            let graph = externals.graph(&diagram, &signatures);
            let symmetry = u64::try_from(diagram.automorphism_count.unwrap()).unwrap();
            generator
                .diagram_input(&resolved, graph, symmetry, multiplicity, &options)
                .unwrap()
        })
        .collect();
    (options, inputs)
}

/// The atoms of a diagram: numerator, prefactor, projector, overall factor and fragments.
fn atoms(diagram: &FeynmanDiagram) -> Vec<Atom> {
    [
        diagram.numerator(),
        diagram.numerator_prefactor(),
        diagram.projector(),
        diagram.overall_factor(),
    ]
    .into_iter()
    .chain(diagram.vertices().map(|(_, vertex)| &vertex.numerator))
    .chain(diagram.edges().map(|(_, _, edge)| &edge.numerator))
    .cloned()
    .collect()
}

/// Finalize every diagram input of `case` with the automatic projector and with its provisional
/// projector, and require identical diagrams, down to the bytes of every atom. Returns the number
/// of memoised automatic projectors.
fn check(name: &str) -> usize {
    let case = Case::named(name);
    let generator = Generator::new(case.model());
    let (options, inputs) = diagram_inputs(&generator, case);
    assert!(!inputs.is_empty(), "{name} has no diagram");
    for (index, input) in inputs.into_iter().enumerate() {
        let context = format!("{name}_{index}");
        for topology in [input.comparison, input.representative] {
            let diagram = |options: &GenerationOptions| {
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
            let projector = provisional_projector(&generator, &topology.graph);
            let automatic = diagram(&options);
            let provisional = diagram(&options.clone().projector(projector));
            assert_ne!(automatic.projector(), &Atom::one(), "{context}");
            assert_eq!(automatic.projector(), provisional.projector(), "{context}");
            assert_eq!(atoms(&automatic), atoms(&provisional), "{context}");
            let json = |diagram: &FeynmanDiagram| diagram.to_json().unwrap();
            assert_eq!(json(&automatic), json(&provisional), "{context}");
        }
    }
    generator.cache.projectors.read().unwrap().len()
}

#[test]
fn amplitude_projectors_match_finalized_provisional_projectors() {
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
fn cross_section_projectors_match_finalized_provisional_projectors() {
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
