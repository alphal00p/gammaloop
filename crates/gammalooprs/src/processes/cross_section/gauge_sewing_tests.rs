//! Complete model-generated covariant cuts, tested against three physical
//! polarizations. Component contraction never expands a physical numerator.

use std::collections::BTreeMap;

use itertools::Itertools;
use linnet::half_edge::{
    involution::{Flow, HedgePair},
    subgraph::Inclusion,
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::{
            LibraryTensor,
            symbolic::{ExplicitKey, TensorLibrary},
        },
        parsing::ParseSettings,
    },
    structure::{
        PermutedStructure,
        representation::{Minkowski, RepName},
    },
    tensors::parametric::{ParamOrConcrete, ParamTensor},
};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    feyngen::{FeynGenFilter, FeynGenFilters, GenerationType},
    initialisation::test_initialise,
    model::UFOSymbol,
    numerator::{ParsingNet, aind::Aind},
    processes::{Process, ProcessCollection, ProcessDefinition},
    settings::GlobalSettings,
    utils::{F, FUN_LIB, GS, load_generic_model},
};

#[test]
fn declared_covariant_cut_states_preserve_multisets_and_physical_labels() -> eyre::Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    for (state, count) in [(vec![24, -24], 16), (vec![23, 23], 10)] {
        let process = ProcessDefinition {
            generation_type: GenerationType::CrossSection,
            final_pdgs_lists: vec![state.clone()],
            ..Default::default()
        };
        let closure = process.covariant_cut_states(&model)?;
        assert_eq!(closure.len(), count);
        assert_eq!(model.covariant_cut_states(&closure)?, closure);
        let labels = process.covariant_cut_representatives(&model);
        for members in closure {
            assert_eq!(
                members
                    .iter()
                    .map(|pdg| labels[&(*pdg as isize)] as i64)
                    .sorted()
                    .collect_vec(),
                state.iter().copied().sorted().collect_vec()
            );
        }
    }
    let diagnostic = ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        final_pdgs_lists: vec![vec![250, 250]],
        ..Default::default()
    };
    assert!(diagnostic.covariant_cut_representatives(&model).is_empty());
    assert_eq!(
        diagnostic.covariant_cut_states(&model)?,
        diagnostic.final_pdgs_lists
    );
    let mixed = ProcessDefinition {
        final_pdgs_lists: vec![vec![23, 23], vec![250, 25]],
        ..diagnostic.clone()
    };
    assert!(
        mixed
            .covariant_cut_states(&model)
            .unwrap_err()
            .to_string()
            .contains("observable labels remain unambiguous")
    );
    let unresolved = ProcessDefinition {
        cross_section_filters: FeynGenFilters(vec![FeynGenFilter::PerturbativeOrders(
            [("QED".into(), 1)].into_iter().collect(),
        )]),
        ..diagnostic
    };
    let mut unresolved_model = model.clone();
    // The default unresolved sets are massless. Exercise the same closure
    // boundary for a declared class that also contains a massive vector.
    unresolved_model
        .unresolved_particles
        .entry("QED".into())
        .or_default()
        .insert(model.get_particle_from_pdg(23));
    assert!(unresolved.covariant_cut_states(&unresolved_model).is_err());

    let incomplete = ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        final_pdgs_lists: vec![vec![24, -24]],
        cross_section_filters: FeynGenFilters(vec![FeynGenFilter::ParticleVeto(vec![251])]),
        ..Default::default()
    };
    assert!(
        incomplete
            .covariant_cut_states(&model)
            .unwrap_err()
            .to_string()
            .contains("complete vector/Goldstone/ghost sector")
    );
    let mut invalid = model.clone();
    invalid.covariant_cut_multiplets.get_mut(&24).unwrap()[3] = -9000003;
    assert!(invalid.validate_covariant_cut_multiplets().is_err());
    let mut invalid = model.clone();
    let position = invalid
        .propagators
        .iter()
        .position(|p| p.particle.name.as_str() == "G0")
        .unwrap();
    std::sync::Arc::make_mut(&mut invalid.propagators[position]).denominator += Atom::one();
    assert!(invalid.validate_covariant_cut_multiplets().is_err());
    Ok(())
}

#[test]
fn generated_higgs_covariant_cuts_equal_three_physical_vector_polarizations() -> eyre::Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let mut settings = GlobalSettings::default();
    settings.n_cores.feyngen = 1;

    let dot = |a: &[Atom; 4], b: &[Atom; 4]| {
        &a[0] * &b[0] - &a[1] * &b[1] - &a[2] * &b[2] - &a[3] * &b[3]
    };
    for (physical, vector_mass, graph_count, symmetry) in [
        (vec![24, -24], Atom::num(3), 6, Atom::one()),
        (vec![23, 23], Atom::num((15, 4)), 4, Atom::num((1, 2))),
    ] {
        let process = ProcessDefinition {
            generation_type: GenerationType::CrossSection,
            initial_pdgs: vec![25],
            final_pdgs_lists: vec![physical],
            loop_count_range: (1, 1),
            cross_section_filters: FeynGenFilters(vec![FeynGenFilter::CouplingOrders(
                [("QED".into(), (2, Some(2)))].into_iter().collect(),
            )]),
            ..Default::default()
        };
        let graphs = process.generate(&model, &settings)?;
        assert_eq!(
            graphs.len(),
            graph_count,
            "the complete covariant Higgs cut inventory"
        );
        // Raw imported states cannot establish physical-vector intent. The
        // same collection with its physical process declaration retains every
        // partner through the import owner and the cut matcher below.
        let inferred =
            ProcessDefinition::from_graph_list(&graphs, GenerationType::CrossSection, &model)?;
        assert!(
            inferred
                .covariant_cut_states(&model)
                .unwrap_err()
                .to_string()
                .contains("--process-spec")
        );
        let imported = Process::from_graph_list(
            "higgs_gauge".into(),
            "born".into(),
            graphs,
            GenerationType::CrossSection,
            Some(process),
            None,
            &model,
        )?;
        let ProcessCollection::CrossSections(cross_sections) = &imported.collection else {
            panic!("forward graphs must import as a cross section");
        };
        let process = &imported.definition;
        for (energy_ratio, momentum_ratio) in [((5, 3), (4, 3)), ((13, 5), (12, 5))] {
            let energy = &vector_mass * Atom::num(energy_ratio);
            let momentum = &vector_mass * Atom::num(momentum_ratio);
            let higgs_mass = Atom::num(2) * &energy;
            // ee=sw=3/5, cw=4/5, MW=3 and vev=6 obey the same action
            // relations as the model. Set MH^2=2 lambda vev^2 on shell.
            let coupling = vector_mass.pow(2) / Atom::num(3);
            for boosted in [false, true] {
                let boost = |p: [Atom; 4]| -> [Atom; 4] {
                    if !boosted {
                        return p;
                    }
                    [
                        Atom::num((5, 4)) * &p[0] + Atom::num((3, 4)) * &p[1],
                        Atom::num((3, 4)) * &p[0] + Atom::num((5, 4)) * &p[1],
                        p[2].clone(),
                        p[3].clone(),
                    ]
                };
                let q = [
                    boost([energy.clone(), Atom::Zero, Atom::Zero, momentum.clone()]),
                    boost([energy.clone(), Atom::Zero, Atom::Zero, -&momentum]),
                ];
                let total: [Atom; 4] = std::array::from_fn(|mu| &q[0][mu] + &q[1][mu]);
                let polarizations = [1, -1].map(|direction| {
                    [
                        boost([Atom::Zero, Atom::one(), Atom::Zero, Atom::Zero]),
                        boost([Atom::Zero, Atom::Zero, Atom::one(), Atom::Zero]),
                        boost([
                            &momentum / &vector_mass,
                            Atom::Zero,
                            Atom::Zero,
                            Atom::num(direction) * &energy / &vector_mass,
                        ]),
                    ]
                });
                let physical_norm = polarizations[0]
                    .iter()
                    .flat_map(|left| polarizations[1].iter().map(|right| dot(left, right).pow(2)))
                    .fold(Atom::Zero, |sum, norm| sum + norm)
                    * coupling.pow(2)
                    * &symmetry;
                assert_eq!(
                    physical_norm,
                    coupling.pow(2)
                        * &symmetry
                        * (Atom::num(2) + dot(&q[0], &q[1]).pow(2) / vector_mass.pow(4))
                );

                let mut channels = BTreeMap::<Vec<i64>, Atom>::new();
                for forward in &cross_sections["born"].supergraphs {
                    let graph = &forward.graph;
                    let cuts = forward.process_valid_cuts(&model, process, &settings.generation)?;
                    assert_eq!(cuts.len(), 1, "each Higgs Born graph has one cut");
                    let cut = &cuts[super::CutId(0)];
                    // Build the generic metric with Atom entries so its
                    // spatial signs remain exact before component contraction.
                    let mut library = spenso_hep_lib::hep_lib_atom::<Aind, F<f64>>();
                    library.insert_generic(
                        TensorLibrary::<ParamTensor<ExplicitKey<Aind>>, Aind>::id(
                            Minkowski {}.into(),
                        ),
                        |key| {
                            ParamOrConcrete::Param(TensorLibrary::<
                                ParamTensor<ExplicitKey<Aind>>,
                                Aind,
                            >::diag_unimodular_metric(
                                key
                            ))
                        },
                    );
                    let cut_edges = graph.iter_edges_of(&cut.cut).collect_vec();
                    let mut particles = Vec::new();
                    for (position, (pair, eid, edge)) in cut_edges.iter().enumerate() {
                        let flow = match pair {
                            HedgePair::Split { split, .. }
                            | HedgePair::Unpaired { flow: split, .. } => *split,
                            HedgePair::Paired { .. } => panic!("cut edge must be split"),
                        };
                        let particle = edge.data.particle().unwrap();
                        particles.push(if flow == Flow::Source {
                            particle.pdg_code
                        } else {
                            particle.get_anti_particle(&model).pdg_code
                        } as i64);
                        let sign = if flow == Flow::Source { 1 } else { -1 };
                        let key = ExplicitKey::from_iter(
                            [Minkowski {}.new_rep(4)],
                            GS.emr_mom,
                            Some(vec![Atom::num(eid.0)]),
                        );
                        let tensor = ParamTensor::from_dense(
                            key.structure,
                            q[position].iter().map(|p| Atom::num(sign) * p).collect(),
                        )?;
                        library.insert_explicit(PermutedStructure::identity(
                            ParamOrConcrete::Param(tensor),
                        ));
                    }
                    for (pair, eid, _) in graph.iter_edges_of(&graph.initial_state_cut) {
                        let source = match pair {
                            HedgePair::Paired { source, .. } | HedgePair::Split { source, .. } => {
                                source
                            }
                            HedgePair::Unpaired { hedge, .. } => hedge,
                        };
                        let sign = if cut.left.includes(&source) { -1 } else { 1 };
                        let key = ExplicitKey::from_iter(
                            [Minkowski {}.new_rep(4)],
                            GS.emr_mom,
                            Some(vec![Atom::num(eid.0)]),
                        );
                        let tensor = ParamTensor::from_dense(
                            key.structure,
                            total.iter().map(|p| Atom::num(sign) * p).collect(),
                        )?;
                        library.insert_explicit(PermutedStructure::identity(
                            ParamOrConcrete::Param(tensor),
                        ));
                    }
                    let mut expression = model.apply_coupling_replacement_rules(
                        &graph.production_numerator_atom_for_full_3d_expression(),
                    );
                    for (name, value) in [
                        ("ee", Atom::num((3, 5))),
                        ("sw", Atom::num((3, 5))),
                        ("cw", Atom::num((4, 5))),
                        ("MW", Atom::num(3)),
                        ("MZ", Atom::num((15, 4))),
                        ("vev", Atom::num(6)),
                        ("lam", higgs_mass.pow(2) / Atom::num(72)),
                        ("MH", higgs_mass.clone()),
                    ] {
                        expression = expression
                            .replace(Atom::from(UFOSymbol::from(name)).to_pattern())
                            .with(value);
                    }
                    let mut network = ParsingNet::try_from_view(
                        expression.as_view(),
                        &library,
                        &ParseSettings::default(),
                    )?;
                    network.execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB)?;
                    // Read the scalar Atom directly instead of materializing
                    // a rank-zero tensor and extracting its component again.
                    let value = match network.result_scalar()? {
                        ExecutionResult::Zero => Atom::Zero,
                        ExecutionResult::One => Atom::one(),
                        ExecutionResult::Val(value) => value.into_owned(),
                    };
                    *channels
                        .entry(particles.into_iter().sorted().collect())
                        .or_default() += value;
                }
                assert_eq!(channels.len(), graph_count);
                for (particles, value) in &channels {
                    let spins = particles
                        .iter()
                        .map(|pdg| model.get_particle_from_pdg(*pdg as isize).spin)
                        .sorted()
                        .collect_vec();
                    let expected = match spins.as_slice() {
                        [3, 3] => Atom::num(4) * coupling.pow(2) * &symmetry,
                        [1, 1] => higgs_mass.pow(4) / Atom::num(36) * &symmetry,
                        [-1, -1] => -coupling.pow(2) / Atom::num(4),
                        [1, 3] => {
                            -coupling.pow(2)
                                * (Atom::num(2) * higgs_mass.pow(2) + vector_mass.pow(2))
                                / (Atom::num(4) * vector_mass.pow(2))
                        }
                        _ => panic!("unexpected Higgs gauge channel {particles:?}"),
                    };
                    assert_eq!(*value, expected, "channel {particles:?}, boosted={boosted}");
                }
                let covariant_norm = channels.values().fold(Atom::Zero, |sum, value| sum + value);
                assert_eq!(
                    covariant_norm, physical_norm,
                    "{channels:?}, boosted={boosted}"
                );
                // The VV channel alone is the incomplete covariant metric sum.
                let vector_channel = process.final_pdgs_lists[0]
                    .iter()
                    .copied()
                    .sorted()
                    .collect_vec();
                assert_eq!(
                    channels[&vector_channel],
                    Atom::num(4) * coupling.pow(2) * &symmetry
                );
                assert_ne!(channels[&vector_channel], physical_norm);
            }
        }
    }
    Ok(())
}
