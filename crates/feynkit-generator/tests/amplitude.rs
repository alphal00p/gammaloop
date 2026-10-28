//! The amplitude boundary consumes finalized diagrams, not generator state.

use feynkit_amplitude::{Amplitude, AmplitudeError, AmplitudeOptions};
use feynkit_generator::{GenerationFilter, GenerationOptions, NumeratorGrouping, Process};
use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use idenso::{
    IndexTooling,
    tensor::{ContractSettings, SymbolicTensor},
};
use serde_json::{Value, json};
use spenso::{
    shadowing::TensorCollectFilter,
    structure::{abstract_index::AbstractIndex, partial::PartialStructure},
};
use std::{collections::BTreeMap, sync::Arc};
use symbolica::atom::{Atom, AtomCore};
use symbolica::id::{ConditionResult, Pattern, Replacement};
use symbolica::{parse, symbol};

fn canonical_tensor(source: SymbolicTensor<PartialStructure>) -> SymbolicTensor<PartialStructure> {
    let contracted = source
        .contract(ContractSettings::default().without_rank_one_tensors())
        .unwrap()
        .simplify_algebra(&idenso::tensor::AlgebraSettings {
            gamma: Some(idenso::dirac::GammaSimplifySettings::default()),
            epsilon: true,
            ..Default::default()
        })
        .unwrap()
        .contract(idenso::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap();
    let canonical = contracted
        .expression()
        .canonize::<AbstractIndex>(AbstractIndex::from)
        .unwrap();
    SymbolicTensor::infer(canonical).unwrap()
}

fn photons() -> Vec<Arc<feynkit_graph::FeynmanDiagram>> {
    let model =
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap();
    let process = Process::new(["e-", "e+"], ["a", "a"]);
    let options = GenerationOptions::default()
        .with_loop_count(0, 0)
        .unwrap()
        .threads(1)
        .max_vertices(2)
        .numerator_grouping(NumeratorGrouping::None)
        .with_graph_filter(GenerationFilter::VertexAllow(vec!["V_98".into()]));
    let generated = process.generate_diagrams(model, &options).unwrap();
    assert_eq!(generated.diagrams.len(), 2);
    generated.diagrams.into_iter().map(Arc::new).collect()
}

#[test]
fn generated_amplitude_aligns_ports_and_conjugates_involutively() {
    let amplitude = Amplitude::from_diagrams(photons()).unwrap();
    assert_eq!(amplitude.legs().len(), 4);
    for term in amplitude.terms() {
        assert_eq!(term.structure(), &amplitude.structure());
    }
    assert_eq!(
        amplitude
            .expression()
            .unwrap()
            .expression()
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .len(),
        4
    );
    let adjoint = amplitude.conjugate().unwrap();
    assert!(adjoint.is_conjugated());
    assert!(
        !adjoint
            .expression()
            .unwrap()
            .expression()
            .to_plain_string()
            .contains("conj")
    );
    let roundtrip = adjoint.conjugate().unwrap();
    assert_eq!(
        canonical_tensor(roundtrip.expression().unwrap())
            .coefficients_equal(
                &canonical_tensor(amplitude.expression().unwrap()),
                TensorCollectFilter::<0>::Tensors,
            )
            .unwrap(),
        ConditionResult::True,
        "conjugation must be involutive by exact coefficient proof"
    );
}

#[test]
fn zero_amplitude_keeps_physical_leg_order_through_adjoint_and_state_sums() {
    let diagrams = photons();
    let original = Amplitude::from_diagram(diagrams[0].clone()).unwrap();
    let zero = diagrams[0]
        .numerator()
        .with_checked_expression(Atom::Zero)
        .unwrap();
    let diagram = diagrams[0].as_ref().clone().with_numerator(zero).unwrap();
    let amplitude = Amplitude::from_diagram(diagram).unwrap();
    assert_eq!(amplitude.legs(), original.legs());
    let expression = amplitude.expression().unwrap();
    assert!(expression.is_zero());
    assert_eq!(expression.structure(), &original.structure());
    let adjoint = amplitude.conjugate().unwrap();
    assert!(adjoint.expression().unwrap().is_zero());
    assert_eq!(
        adjoint.expression().unwrap().structure(),
        &adjoint.structure()
    );
    let squared = amplitude.squared().unwrap();
    assert!(squared.expression().is_zero());
    assert_eq!(squared.expression().rank(), 8);
    let summed = squared
        .sum_spins(&[0, 1, 2, 3], true, &BTreeMap::new(), &BTreeMap::new())
        .unwrap();
    assert!(summed.expression().is_zero());
    assert_eq!(summed.expression().rank(), 0);
}

#[test]
fn squared_ports_remain_open_until_explicit_state_sums() {
    let amplitude = Amplitude::from_diagrams(photons()).unwrap();
    let squared = amplitude.squared().unwrap();
    assert_eq!(
        squared
            .expression()
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
            .expression()
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
            .clone()
    };
    let single = normalize(&one);
    let expected = single
        .with_checked_expression(single.expression() * Atom::num(4))
        .unwrap();
    // Two identical diagrams give |2 A|² = 4 |A|², not the diagonal-only 2 |A|².
    assert_eq!(
        canonical_tensor(normalize(&doubled))
            .coefficients_equal(
                &canonical_tensor(expected),
                TensorCollectFilter::<0>::Tensors,
            )
            .unwrap(),
        ConditionResult::True,
        "the coherent square must include every interference term"
    );
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
        amplitude
            .conjugate()
            .unwrap()
            .expression()
            .unwrap()
            .expression()
            .clone(),
        parse!("spenso::conj(amplitude_test::z)")
    );
    assert_eq!(
        amplitude.squared().unwrap().expression().expression(),
        &parse!("amplitude_test::z*spenso::conj(amplitude_test::z)")
    );
}

#[test]
fn real_parameter_wildcard_suffixes_preserve_complex_scalar_norms() {
    for name in ["g_", "g__", "g___"] {
        let mut definition: Value =
            serde_json::from_str(&Model::phi4().to_json().unwrap()).unwrap();
        definition["parameters"].as_array_mut().unwrap().extend([
            json!({
                "name": name, "nature": "external", "parameter_type": "real",
                "value": [1.0, 0.0], "expression": null
            }),
            json!({
                "name": "h", "nature": "external", "parameter_type": "complex",
                "value": [1.0, 2.0], "expression": null
            }),
        ]);
        definition["couplings"][0]["expression"] = format!("{name}*h").into();
        let model = Model::from_json(&definition.to_string()).unwrap();
        let amplitude = Process::new(["phi", "phi"], ["phi", "phi"])
            .generate_amplitude(
                model,
                &GenerationOptions::default().threads(1).max_vertices(1),
                AmplitudeOptions::default(),
            )
            .unwrap();
        assert_eq!(amplitude.diagrams().len(), 1);
        let squared = amplitude
            .squared()
            .unwrap()
            .sum_spins(&[0, 1, 2, 3], true, &BTreeMap::new(), &BTreeMap::new())
            .unwrap()
            .sum_colors(&[0, 1, 2, 3], true)
            .unwrap();
        let g = Atom::var(symbol!(&format!("UFO::{name}")));
        let h = Atom::var(symbol!("UFO::h"));
        let conjugate_h = symbol!("spenso::conj").call(&h);
        let expected = g.pow(2) * &h * &conjugate_h;
        assert!(
            (squared.expression().expression() - expected)
                .expand()
                .is_zero(),
            "{name}: {}",
            squared.expression().expression()
        );
        let norm = squared.expression().expression().replace_multiple([
            Replacement::new(
                Pattern::Literal(conjugate_h),
                Pattern::Literal(Atom::one() - 2 * Atom::i()),
            ),
            Replacement::new(
                Pattern::Literal(h),
                Pattern::Literal(Atom::one() + 2 * Atom::i()),
            ),
            Replacement::new(Pattern::Literal(g), Pattern::Literal(Atom::one())),
        ]);
        assert_eq!(norm.expand(), Atom::num(5), "{name}");
    }
}

#[test]
fn explicit_reality_assumptions_are_literal_expressions() {
    let real = parse!("amplitude_reality_test::r_");
    let diagram = FeynmanDiagram::from_dot(
        Arc::new(Model::phi3()),
        r#"digraph {
        a [num="amplitude_reality_test::D_*amplitude_reality_test::r_*amplitude_reality_test::h"];
        ext [style=invis];
        ext -> a [particle="phi"];
        a -> ext [particle="phi"];
        a -> ext [particle="phi"];
    }"#,
    )
    .unwrap();
    let diagram = Arc::new(diagram);
    for (dimension, expected) in [
        (
            AmplitudeOptions::default().dimension,
            parse!(
                "amplitude_reality_test::r_*spenso::conj(amplitude_reality_test::D_)*spenso::conj(amplitude_reality_test::h)"
            ),
        ),
        (
            symbol!("amplitude_reality_test::D_").into(),
            parse!(
                "amplitude_reality_test::D_*amplitude_reality_test::r_*spenso::conj(amplitude_reality_test::h)"
            ),
        ),
    ] {
        let amplitude = Amplitude::new(
            [diagram.clone()],
            AmplitudeOptions {
                dimension,
                real: vec![real.clone()],
            },
        )
        .unwrap();
        assert_eq!(
            amplitude
                .conjugate()
                .unwrap()
                .expression()
                .unwrap()
                .expression(),
            &expected
        );
    }
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
    let process = Process::new(["b", "b"], ["b", "b"]);
    let options = GenerationOptions::default()
        .with_loop_count(0, 0)
        .unwrap()
        .threads(1)
        .max_vertices(2)
        .with_graph_filter(GenerationFilter::VertexAllow(vertices));
    let generated = process.generate_diagrams(model, &options).unwrap();
    assert_eq!(generated.diagrams.len(), 2);
    let amplitude = Amplitude::from_diagrams(generated.diagrams.into_iter().map(Arc::new)).unwrap();
    assert!(amplitude.legs().iter().all(|leg| leg.slots.len() == 2));
    let squared = amplitude.squared().unwrap();
    assert_eq!(
        squared
            .expression()
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
            .expression()
            .contains(parse!("gammalooprs::K(0,spenso::mink(4))"))
    );
    assert!(squared.expression().expression().contains(parse!(
        "gammalooprs::K(feynkit_amplitude::bra(0),spenso::mink(4))"
    )));
}

#[test]
fn one_process_generates_diagrams_amplitudes_and_cross_sections() {
    let model = Arc::new(
        Model::from_json(include_str!("../../feynkit-model/tests/fixtures/sm.json")).unwrap(),
    );
    let process = Process::new(["e-", "e+"], ["a", "a"]);
    let options = GenerationOptions::default().threads(1).max_vertices(4);
    let diagrams = process.generate_diagrams(model.clone(), &options).unwrap();
    assert_eq!(diagrams.diagrams.len(), 2);
    let amplitude = process
        .generate_amplitude(model.clone(), &options, AmplitudeOptions::default())
        .unwrap();
    assert_eq!(amplitude.diagrams().len(), diagrams.diagrams.len());
    assert_eq!(
        amplitude.expression().unwrap(),
        Amplitude::from_diagrams(diagrams.diagrams.into_iter().map(Arc::new))
            .unwrap()
            .expression()
            .unwrap()
    );
    let forward = process
        .generate_cross_section(
            model.clone(),
            &options.clone().with_loop_count(1, 1).unwrap(),
        )
        .unwrap();
    assert!(!forward.diagrams.is_empty());
    assert!(
        forward
            .diagrams
            .iter()
            .all(|diagram| diagram.loop_count() == 1 && !diagram.cuts().is_empty())
    );
    assert_eq!(options.loop_count(), 0..=0);
    let token = feynkit_generator::CancellationToken::new();
    token.cancel();
    let cancelled = options.cancellation_token(token);
    assert!(matches!(
        process.generate_amplitude(model, &cancelled, AmplitudeOptions::default()),
        Err(feynkit_generator::GenerationError::IncompleteAmplitude)
    ));
}
