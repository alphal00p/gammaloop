//! The amplitude boundary consumes finalized diagrams, not generator state.

use feynkit_amplitude::{Amplitude, AmplitudeError};
use feynkit_generator::{
    GenerationFilter, GenerationOptions, Generator, NumeratorGrouping, Process,
};
use feynkit_model::Model;
use idenso::{IndexTooling, shorthands::metric::MetricSimplifier};
use spenso::structure::abstract_index::AbstractIndex;
use std::{collections::BTreeMap, sync::Arc};
use symbolica::{
    atom::{Atom, AtomCore},
    parse,
};

fn photons() -> Vec<Arc<feynkit_graph::FeynmanDiagram>> {
    let model =
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
    let process = Process::amplitude(["e-", "e+"], ["a", "a"])
        .with_loop_count(0, 0)
        .unwrap();
    let options = GenerationOptions::default()
        .threads(1)
        .max_vertices(2)
        .numerator_grouping(NumeratorGrouping::None)
        .with_graph_filter(GenerationFilter::VertexAllow(vec!["V_98".into()]));
    let generated = Generator::new(model).generate(&process, &options).unwrap();
    assert_eq!(generated.diagrams.len(), 2);
    generated.diagrams.into_iter().map(Arc::new).collect()
}

#[test]
fn generated_amplitude_aligns_ports_and_conjugates_involutively() {
    let amplitude = Amplitude::from_diagrams(photons()).unwrap();
    assert_eq!(amplitude.legs().len(), 4);
    assert_eq!(
        amplitude
            .expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .len(),
        4
    );
    let adjoint = amplitude.conjugate().unwrap();
    assert!(adjoint.is_conjugated());
    assert!(!adjoint.expression().to_plain_string().contains("conj"));
    let roundtrip = adjoint.conjugate().unwrap();
    let canonical = |a: Atom| {
        a.expand()
            .simplify_metrics()
            .canonize::<AbstractIndex>(AbstractIndex::from)
            .unwrap()
    };
    assert_eq!(
        canonical(roundtrip.expression()),
        canonical(amplitude.expression())
    );
}

#[test]
fn squared_ports_remain_open_until_explicit_state_sums() {
    let amplitude = Amplitude::from_diagrams(photons()).unwrap();
    let squared = amplitude.squared().unwrap();
    assert_eq!(
        squared
            .expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .len(),
        8
    );
    let fermions = squared
        .sum_spins(&[0, 1], true, &BTreeMap::new(), &BTreeMap::new())
        .unwrap();
    assert_eq!(
        fermions
            .expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .len(),
        4
    );
    let all = fermions
        .sum_spins(&[2, 3], false, &BTreeMap::new(), &BTreeMap::new())
        .unwrap();
    assert!(
        all.expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .is_empty()
    );
    assert!(matches!(
        all.sum_spins(&[0], true, &BTreeMap::new(), &BTreeMap::new()),
        Err(AmplitudeError::AlreadySummed { .. })
    ));
    assert!(squared.spin_summed().is_empty());
}

#[test]
fn coherent_square_includes_cross_diagram_interference() {
    let diagrams = photons();
    let one = Amplitude::from_diagram(diagrams[0].clone()).unwrap();
    let doubled = Amplitude::from_diagrams([diagrams[0].clone(), diagrams[0].clone()]).unwrap();
    let normalize = |a: &Amplitude| {
        a.squared()
            .unwrap()
            .sum_spins(&[0, 1, 2, 3], true, &BTreeMap::new(), &BTreeMap::new())
            .unwrap()
            .expression()
            .expand()
            .simplify_metrics()
            .canonize::<AbstractIndex>(AbstractIndex::from)
            .unwrap()
    };
    // Two identical diagrams give |2 A|² = 4 |A|², not the diagonal-only 2 |A|².
    let difference = normalize(&doubled) - normalize(&one) * Atom::num(4);
    assert!(difference.expand().is_zero());
}

#[test]
fn rejects_mixed_external_states_and_unknown_sum_labels() {
    assert!(matches!(
        Amplitude::from_diagrams([]),
        Err(AmplitudeError::Empty)
    ));
    let diagrams = photons();
    let changed = diagrams[0]
        .map_data(
            |_, v| v.clone(),
            |_, _, e| {
                let mut e = e.clone();
                if let Some(leg) = &mut e.external {
                    leg.index += 10;
                }
                e
            },
        )
        .unwrap();
    assert!(matches!(
        Amplitude::from_diagrams([diagrams[0].clone(), Arc::new(changed)]),
        Err(AmplitudeError::DifferentExternals(_))
    ));
    let squared = Amplitude::from_diagrams(diagrams)
        .unwrap()
        .squared()
        .unwrap();
    assert!(matches!(
        squared.sum_colors(&[99], false),
        Err(AmplitudeError::UnknownLeg(99))
    ));
}

#[test]
fn scalar_couplings_are_conjugated_without_assuming_they_are_real() {
    let model = Arc::new(
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
    );
    let diagram = feynkit_graph::FeynmanDiagram::from_dot(
        model,
        r#"digraph {
        a [num="amplitude_test::z"];
        ext [style=invis];
        ext -> a [particle="H"];
        a -> ext [particle="H"];
        a -> ext [particle="H"];
    }"#,
    )
    .unwrap();
    let amplitude = Amplitude::from_diagram(diagram).unwrap();
    assert_eq!(
        amplitude.conjugate().unwrap().expression(),
        parse!("spenso::conj(amplitude_test::z)")
    );
    assert_eq!(
        amplitude.squared().unwrap().expression(),
        &parse!("amplitude_test::z*spenso::conj(amplitude_test::z)")
    );
}

#[test]
fn colored_fermion_interference_preserves_physical_leg_pairings() {
    let model =
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
    let mut species = [
        model.particle_id("b").unwrap(),
        model.particle_id("b~").unwrap(),
        model.particle_id("g").unwrap(),
    ];
    species.sort();
    let vertices = model
        .vertex_rules()
        .iter()
        .filter(|rule| {
            let mut particles = rule.particles.clone();
            particles.sort();
            particles == species
        })
        .map(|rule| rule.name.clone().into())
        .collect();
    let process = Process::amplitude(["b", "b"], ["b", "b"])
        .with_loop_count(0, 0)
        .unwrap();
    let options = GenerationOptions::default()
        .threads(1)
        .max_vertices(2)
        .with_graph_filter(GenerationFilter::VertexAllow(vertices));
    let generated = Generator::new(model).generate(&process, &options).unwrap();
    assert_eq!(generated.diagrams.len(), 2);
    let amplitude = Amplitude::from_diagrams(generated.diagrams.into_iter().map(Arc::new)).unwrap();
    assert!(amplitude.legs().iter().all(|leg| leg.slots.len() == 2));
    let squared = amplitude.squared().unwrap();
    assert_eq!(
        squared
            .expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .len(),
        16
    );
    let colors = squared.sum_colors(&[0, 1, 2, 3], true).unwrap();
    assert_eq!(
        colors
            .expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .len(),
        8
    );
    let summed = colors
        .sum_spins(&[0, 1, 2, 3], true, &BTreeMap::new(), &BTreeMap::new())
        .unwrap();
    assert!(
        summed
            .expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .is_empty()
    );
}

#[test]
fn loop_square_keeps_independent_integration_momenta() {
    let model = Arc::new(
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
    );
    let diagram = feynkit_graph::FeynmanDiagram::from_dot(
        model,
        r#"digraph triangle {
        ext [style=invis];
        ext -> a [id=0, particle="H"];
        b -> ext [id=1, particle="H"];
        c -> ext [id=2, particle="H"];
        a -> b [id=3, particle="H", lmb_id=0];
        b -> c [id=4, particle="H"];
        c -> a [id=5, particle="H"];
    }"#,
    )
    .unwrap();
    let squared = Amplitude::from_diagram(diagram).unwrap().squared().unwrap();
    assert!(
        squared
            .expression()
            .contains(&parse!("gammalooprs::K(0,spenso::mink(4))"))
    );
    assert!(squared.expression().contains(&parse!(
        "gammalooprs::K(feynkit_amplitude::bra(0),spenso::mink(4))"
    )));
}
