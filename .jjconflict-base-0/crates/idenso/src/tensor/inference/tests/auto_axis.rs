use super::*;
use spenso::{
    network::{
        ExecutionResult, Network, Sequential, SmallestDegree,
        library::{panicing::ErroringLibrary, symbolic::TensorLibrary},
        parsing::{ParseSettings, ShadowedStructure},
        store::NetworkStore,
    },
    structure::{representation::Euclidean, slot::ParseableAind},
    tensors::parametric::ParamTensor,
};
use symbolica::function;

#[test]
fn changed_auto_axis_order_cannot_relabel_asymmetric_components() {
    crate::test_support::test_initialize();
    type Tensor = ParamTensor<ShadowedStructure<AbstractIndex>>;
    type Lib = TensorLibrary<ParamTensor<ExplicitKey<AbstractIndex>>, AbstractIndex>;
    type Net = Network<NetworkStore<Tensor, Atom>, ExplicitKey<AbstractIndex>, Symbol>;

    let rep = LibraryRep::from(Euclidean {}).new_rep(2);
    let slot = |index| rep.slot::<AbstractIndex, _>(index).to_atom();
    let mut heads = [
        spenso::vector_symbol!("auto_axis_components::first"),
        spenso::vector_symbol!("auto_axis_components::second"),
        spenso::vector_symbol!("auto_axis_components::third"),
    ];
    heads.sort_by_key(|&head| function!(head, slot(AbstractIndex::Normal(0))));
    let [low, middle, high] = heads;
    let selectors = [
        spenso::vector_symbol!("auto_axis_components::basis_zero"),
        spenso::vector_symbol!("auto_axis_components::basis_one"),
    ];
    let mut library = Lib::new();
    for (head, components) in [
        (middle, [2, 3]),
        (low, [5, 7]),
        (high, [5, 7]),
        (selectors[0], [1, 0]),
        (selectors[1], [0, 1]),
    ] {
        library
            .insert_explicit_sparse(
                ExplicitKey::from_iter([rep], head, None),
                components
                    .into_iter()
                    .enumerate()
                    .map(|(i, v)| (vec![i], Atom::num(v))),
                Atom::Zero,
            )
            .unwrap();
    }
    let indices = [74821, 74823].map(AbstractIndex::Normal);
    let component = |expression: Atom| {
        let projection = expression
            * function!(selectors[0], slot(indices[0]))
            * function!(selectors[1], slot(indices[1]));
        let mut network = Net::try_from_view::<ShadowedStructure<AbstractIndex>, _>(
            projection.as_view(),
            &library,
            &ParseSettings::default(),
        )
        .unwrap();
        network
            .execute::<Sequential, SmallestDegree, _, _, _>(&library, &ErroringLibrary::new())
            .unwrap();
        match network.result_scalar().unwrap() {
            ExecutionResult::Val(value) => value.into_owned(),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::One => Atom::one(),
        }
    };
    let owner = AbstractIndex::fresh_open_owner();
    let ends = [0, 1].map(|axis| AbstractIndex::Open { owner, axis });
    let original = function!(middle, slot(ends[0])) * function!(high, slot(ends[1]));
    // The replacement vector has the same components, but its normalized factor
    // moves before the other vector. The encoded axes have not changed identity.
    let changed = function!(middle, slot(ends[0])) * function!(low, slot(ends[1]));
    let (before, _) =
        SymbolicTensor::observed_interface(&original, LeafInference::ObserveEncoded).unwrap();
    let (after, _) =
        SymbolicTensor::observed_interface(&changed, LeafInference::ObserveEncoded).unwrap();
    assert_eq!(
        before
            .logical_slots()
            .iter()
            .map(|s| s.aind)
            .collect::<Vec<_>>(),
        ends.map(PartialIndex::Explicit)
    );
    assert_eq!(
        after
            .logical_slots()
            .iter()
            .map(|s| s.aind)
            .collect::<Vec<_>>(),
        [ends[1], ends[0]].map(PartialIndex::Explicit)
    );

    let source = SymbolicTensor::infer(original).unwrap();
    let bindings = HashMap::from([(0, indices[0]), (1, indices[1])]);
    assert_eq!(
        component(source.materialize_interface_ports(&bindings).unwrap()),
        Atom::num(14)
    );
    let physical = changed
        .replace(ends[0].to_atom())
        .with(indices[0].to_atom())
        .replace(ends[1].to_atom())
        .with(indices[1].to_atom());
    assert_eq!(component(physical), Atom::num(14));
    // This is the old unsafe result: keeping only the anonymous layout causes
    // positional materialization to swap the two component coordinates.
    let positional = SymbolicTensor::new(changed.clone(), source.structure.clone());
    assert_eq!(
        component(positional.materialize_interface_ports(&bindings).unwrap()),
        Atom::num(15)
    );
    assert!(source.with_rewritten_expression(changed.clone()).is_err());
    assert!(source.with_algebra_result(changed.clone()).is_err());
    let mixed = source.expression() + &changed;
    assert!(source.with_rewritten_expression(mixed.clone()).is_err());
    assert!(source.with_algebra_result(mixed).is_err());

    for finish in [
        SymbolicTensor::with_rewritten_expression,
        SymbolicTensor::with_algebra_result,
    ] {
        let doubled = finish(&source, Atom::num(2) * source.expression()).unwrap();
        assert_eq!(
            component(doubled.materialize_interface_ports(&bindings).unwrap()),
            Atom::num(28)
        );
        assert_eq!(
            finish(&source, Atom::Zero).unwrap().structure,
            source.structure
        );
        assert_eq!(finish(&source, source.expression.clone()).unwrap(), source);
    }
}

#[test]
fn alias_admission_preserves_scalar_metadata_and_mixed_auto_representations() {
    crate::test_support::test_initialize();
    let euc = LibraryRep::from(Euclidean {}).new_rep(2);
    let mink = LibraryRep::from(Minkowski {}).new_rep(4);
    let source = function!(
        spenso::tensor_symbol!("auto_axis_components::T"),
        Atom::num(7),
        euc.to_symbolic([]),
        mink.to_symbolic([])
    );
    let value = SymbolicTensor::infer(source).unwrap();
    assert_eq!(value.rank(), 2);
    assert_eq!(
        value
            .structure
            .logical_slots()
            .iter()
            .map(|s| s.rep())
            .collect::<Vec<_>>(),
        [euc, mink]
    );
    let handle = value.alias_handle().unwrap();
    let registered = handle
        .clone()
        .with_aliases([(handle, value.clone())])
        .unwrap();
    let restored = registered.resolved().unwrap();
    assert_eq!(restored.structure, value.structure);
    assert_eq!(restored.expression, value.expression);
}

#[test]
fn one_anonymous_axis_stays_unambiguous_across_explicit_port_movement() {
    crate::test_support::test_initialize();
    let rep = LibraryRep::from(Euclidean {}).new_rep(2);
    let explicit = rep
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(74831))
        .to_atom();
    let head = spenso::tensor_symbol!("auto_axis_components::mixed_anonymous");
    let atom = function!(head, explicit.clone(), rep.to_symbolic([]));
    let source = SymbolicTensor::infer(atom).unwrap();
    let changed = function!(head, rep.to_symbolic([]), explicit);
    for finish in [
        SymbolicTensor::with_rewritten_expression,
        SymbolicTensor::with_algebra_result,
    ] {
        let result = finish(&source, changed.clone()).unwrap();
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.expression, changed);
    }
    let ambiguous =
        SymbolicTensor::infer(function!(head, rep.to_symbolic([]), rep.to_symbolic([]))).unwrap();
    assert!(
        ambiguous
            .with_rewritten_expression(Atom::num(2) * ambiguous.expression())
            .is_err()
    );
}

#[test]
fn unchanged_ordered_products_keep_raw_and_encoded_auto_ports() {
    use crate::shorthands::UndoShorthands;

    crate::test_support::test_initialize();
    let rep = LibraryRep::from(Euclidean {}).new_rep(2);
    let mut heads = [
        spenso::vector_symbol!("ordered_noop_components::first"),
        spenso::vector_symbol!("ordered_noop_components::second"),
    ];
    heads.sort_by_key(|&head| {
        function!(
            head,
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(0))
                .to_atom()
        )
    });
    let [low, high] = heads;
    let indices = [74891, 74892].map(AbstractIndex::Normal);
    let indexed = |head, index| function!(head, rep.slot::<AbstractIndex, _>(index).to_atom());
    let components = [
        (indexed(high, indices[0]), 2),
        (indexed(high, indices[1]), 3),
        (indexed(low, indices[0]), 5),
        (indexed(low, indices[1]), 7),
    ];
    for encoded in [false, true] {
        let owner = AbstractIndex::fresh_open_owner();
        let axes = [0, 1].map(|axis| AbstractIndex::Open { owner, axis });
        let factors = if encoded {
            [indexed(high, axes[0]), indexed(low, axes[1])]
        } else {
            [
                function!(high, rep.to_symbolic([])),
                function!(low, rep.to_symbolic([])),
            ]
        };
        let original = spenso::bracket!(&factors[0], &factors[1]);
        let source = SymbolicTensor::infer(original.clone()).unwrap();
        assert_eq!(source.expression(), &original);
        assert_eq!(
            crate::shorthands::bracket::BracketNormalizer::normalize(source.expression().as_view()),
            original,
            "raw and encoded AUTO ports both retain their ordered bracket"
        );
        let result = source
            .contract(Default::default())
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(result.expression(), source.expression());
        assert_eq!(result.structure(), source.structure());
        if encoded {
            for value in [&source, &result] {
                let observed =
                    InterfaceInference::replacement_interface(value.expression().as_view())
                        .unwrap();
                assert_eq!(
                    observed
                        .logical_slots()
                        .iter()
                        .map(|s| s.aind)
                        .collect::<Vec<_>>(),
                    axes.map(PartialIndex::Explicit)
                );
            }
            let reordered = &factors[0] * &factors[1];
            let observed = InterfaceInference::replacement_interface(reordered.as_view()).unwrap();
            assert_eq!(
                observed
                    .logical_slots()
                    .iter()
                    .map(|s| s.aind)
                    .collect::<Vec<_>>(),
                [axes[1], axes[0]].map(PartialIndex::Explicit)
            );
            assert!(source.with_rewritten_expression(reordered).is_err());
        }
        for (positions, expected_component) in [(indices, 14), ([indices[1], indices[0]], 15)] {
            let bindings = HashMap::from([(0, positions[0]), (1, positions[1])]);
            for value in [&source, &result] {
                let materialized = value.materialize_interface_ports(&bindings).unwrap();
                assert_eq!(
                    materialized,
                    spenso::bracket!(indexed(high, positions[0]), indexed(low, positions[1]))
                );
                // Port materialization retains ordered syntax. Explicit shorthand
                // materialization now has named ports, so it can emit an ordinary product.
                let materialized = materialized.undo_all::<AbstractIndex>().unwrap();
                assert_eq!(
                    materialized,
                    indexed(high, positions[0]) * indexed(low, positions[1])
                );
                let component = materialized.replace_map_bottom_up(|node, _, out| {
                    if let Some((_, number)) =
                        components.iter().find(|(atom, _)| atom.as_view() == node)
                    {
                        **out = Atom::num(*number);
                    }
                });
                assert_eq!(component, Atom::num(expected_component));
            }
        }
    }
}
