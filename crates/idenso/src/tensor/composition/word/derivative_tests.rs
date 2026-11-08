use super::*;
use crate::tensor::AlgebraSettings;
use crate::test_support;
use crate::{bis, gamma};
use spenso::mink;
use spenso::shadowing::FactorProjector;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::representation::Minkowski;
use spenso::structure::representation::{LibraryRep, RepName, Representation};
use spenso::structure::slot::IsAbstractSlot;
use std::slice;
use symbolica::{function, symbol};

type Tensor = SymbolicTensor<PartialStructure>;

fn head() -> Symbol {
    test_support::test_initialize();
    symbol!("word_derivative::exterior"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word)
}

fn leaf(name: &str, slots: &[Atom]) -> Tensor {
    Tensor::infer(
        FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
            .add_args(slots)
            .finish(),
    )
    .unwrap()
}

fn target(slots: &[Atom]) -> PartialStructure {
    leaf("word_derivative::target", slots).structure
}

/// A test domain's field-family policy; all tensor algebra remains native.
fn tangent(
    generator: &Tensor,
    selected: Symbol,
    target: &PartialStructure,
) -> Result<Option<Tensor>, TensorInferenceError> {
    if generator.expression.get_symbol() != Some(selected) {
        return Ok(None);
    }
    let occurrence = generator.structure.logical_slots();
    let variable = target.logical_slots();
    assert_eq!(occurrence.len(), variable.len());
    let mut expression = Atom::one();
    for (source, target) in occurrence.iter().zip(&variable) {
        assert_eq!(source.rep(), target.rep());
        let (PartialIndex::Explicit(source_index), PartialIndex::Explicit(target_index)) =
            (source.aind, target.aind)
        else {
            unreachable!()
        };
        expression *= source
            .rep()
            .id(target_index.to_atom(), source_index.to_atom());
    }
    Ok(Some(Tensor::from_declared_parts(
        expression,
        PartialStructure::from_logical_slots(
            (occurrence)
                .into_iter()
                .chain(variable.into_iter().map(|slot| slot.dual())),
        ),
    )?))
}

fn linear_coefficient(value: &Atom, parameter: Symbol) -> Atom {
    value
        .coefficient_list::<i32>(&[Atom::var(parameter)])
        .into_iter()
        .find(|(power, _)| power == &Atom::var(parameter))
        .map(|(_, coefficient)| coefficient)
        .unwrap_or(Atom::Zero)
}

fn canonical(value: Tensor) -> Atom {
    value
        .simplify_algebra(&AlgebraSettings::default())
        .unwrap()
        .canonize()
        .unwrap()
        .expression
}

#[test]
fn word_derivative_left_and_right_reconstruct_native_one_to_four_and_ghost_shifts() {
    let head = head();
    let t = symbol!("word_derivative::t"; Scalar);
    for count in 1..=4 {
        let factors = (0..count)
            .map(|i| leaf(&format!("word_derivative::f{i}"), &[mink!(4, 10 + i)]))
            .collect::<Vec<_>>();
        let source = Tensor::word(head, &factors).unwrap();
        for (position, generator) in factors.iter().enumerate() {
            let variable = target(&[mink!(4, 100 + position)]);
            let response = leaf("word_derivative::response", &[mink!(4, 100 + position)]);
            let local_response = leaf("word_derivative::response", &[mink!(4, 10 + position)]);
            let shifted = generator
                .try_add(
                    &local_response
                        .with_algebra_result(Atom::var(t) * local_response.expression())
                        .unwrap(),
                )
                .unwrap();
            let shifted = source
                .expression
                .replace(generator.expression.to_pattern())
                .with(shifted.expression);
            let expected = canonical(
                source
                    .with_algebra_result(linear_coefficient(&shifted, t))
                    .unwrap(),
            );
            for side in [WordDerivativeSide::Left, WordDerivativeSide::Right] {
                let cut = source
                    .differentiate_word_with(head, &variable, side, |g| {
                        tangent(g, generator.expression.get_symbol().unwrap(), &variable)
                    })
                    .unwrap();
                assert_eq!(
                    cut.structure.logical_slots(),
                    PartialStructure::from_logical_slots(
                        (source.structure.logical_slots())
                            .into_iter()
                            .chain(variable.logical_slots().into_iter().map(|slot| slot.dual()))
                    )
                    .logical_slots()
                );
                let response = Tensor::word(head, slice::from_ref(&response)).unwrap();
                let reconstructed = match side {
                    WordDerivativeSide::Left => response.multiply_word(head, &cut),
                    WordDerivativeSide::Right => cut.multiply_word(head, &response),
                }
                .unwrap();
                assert_eq!(
                    canonical(reconstructed),
                    expected,
                    "count={count}, position={position}, side={side:?}"
                );
            }
        }
    }
    let ghost = leaf("word_derivative::ghost", &[]);
    let anti = leaf("word_derivative::antighost", &[]);
    let source = Tensor::word(head, &[anti.clone(), ghost.clone()]).unwrap();
    let variable = target(&[]);
    let eta = leaf("word_derivative::eta", &[]);
    let shifted = ghost
        .try_add(
            &eta.with_algebra_result(Atom::var(t) * eta.expression())
                .unwrap(),
        )
        .unwrap();
    let expected = linear_coefficient(
        &source
            .expression
            .replace(ghost.expression.to_pattern())
            .with(shifted.expression),
        t,
    );
    for side in [WordDerivativeSide::Left, WordDerivativeSide::Right] {
        let cut = source
            .differentiate_word_with(head, &variable, side, |g| {
                tangent(g, ghost.expression.get_symbol().unwrap(), &variable)
            })
            .unwrap();
        let response = Tensor::word(head, slice::from_ref(&eta)).unwrap();
        let reconstructed = match side {
            WordDerivativeSide::Left => response.multiply_word(head, &cut),
            WordDerivativeSide::Right => cut.multiply_word(head, &response),
        }
        .unwrap();
        assert_eq!(reconstructed.expression.expand(), expected.expand());
    }
}

#[test]
fn word_derivative_keeps_dual_target_variance_and_spectator_ports() {
    let head = head();
    let charge = LibraryRep::new_dual("word_derivative::charge")
        .unwrap()
        .new_rep(3);
    let flavor = LibraryRep::new_dual("word_derivative::flavor")
        .unwrap()
        .new_rep(5);
    let slot = |rep: Representation<LibraryRep>, index| {
        rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    };
    let field = leaf(
        "word_derivative::charged",
        &[slot(charge, 701), slot(flavor.dual(), 702)],
    );
    let partner = leaf(
        "word_derivative::partner",
        &[slot(charge.dual(), 706), slot(flavor, 707)],
    );
    let coefficient = leaf("word_derivative::spectator_coefficient", &[mink!(4, 703)]);
    let source = coefficient
        .multiply_word(
            head,
            &Tensor::word(head, &[field.clone(), partner.clone()]).unwrap(),
        )
        .unwrap();
    let variable = target(&[slot(charge, 704), slot(flavor.dual(), 705)]);
    let response = leaf(
        "word_derivative::charged_response",
        &[slot(charge, 704), slot(flavor.dual(), 705)],
    );
    let expected = coefficient
        .multiply_word(
            head,
            &Tensor::word(
                head,
                &[
                    leaf(
                        "word_derivative::charged_response",
                        &[slot(charge, 701), slot(flavor.dual(), 702)],
                    ),
                    partner,
                ],
            )
            .unwrap(),
        )
        .unwrap();
    let selected = field.expression.get_symbol().unwrap();
    for side in [WordDerivativeSide::Left, WordDerivativeSide::Right] {
        let cut = source
            .differentiate_word_with(head, &variable, side, |g| tangent(g, selected, &variable))
            .unwrap();
        assert_eq!(
            cut.structure.logical_slots(),
            PartialStructure::from_logical_slots(
                (source.structure.logical_slots())
                    .into_iter()
                    .chain(variable.logical_slots().into_iter().map(|s| s.dual()))
            )
            .logical_slots()
        );
        assert_eq!(
            cut.structure.logical_slots()[source.rank()].rep(),
            charge.dual()
        );
        assert_eq!(
            cut.structure.logical_slots()[source.rank() + 1].rep(),
            flavor
        );
        let response = Tensor::word(head, slice::from_ref(&response)).unwrap();
        let reconstructed = match side {
            WordDerivativeSide::Left => response.multiply_word(head, &cut),
            WordDerivativeSide::Right => cut.multiply_word(head, &response),
        }
        .unwrap();
        assert_eq!(canonical(reconstructed), canonical(expected.clone()));
    }
    let zero = source
        .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |g| {
            tangent(g, selected, &variable)?
                .map(|map| map.with_algebra_result(Atom::Zero))
                .transpose()
        })
        .unwrap();
    assert!(zero.is_zero());
    assert_eq!(zero.rank(), source.rank() + variable.logical_slots().len());
    let wrong_variance = Tensor::from_declared_parts(
        Atom::Zero,
        PartialStructure::from_logical_slots(
            (field.structure.logical_slots())
                .into_iter()
                .chain(variable.logical_slots()),
        ),
    )
    .unwrap();
    assert!(
        source
            .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |_| Ok(Some(
                wrong_variance.clone()
            )))
            .is_err()
    );
}

#[test]
fn word_derivative_retains_gamma_coefficients_and_explicit_spinor_target() {
    let head = head();
    let factors = (0..4)
        .map(|i| leaf(&format!("word_derivative::spinor{i}"), &[bis!(4, 200 + i)]))
        .collect::<Vec<_>>();
    let coefficient = Tensor::infer(
        (Atom::num((3, 5)) + Atom::i()) * gamma!(200, 201, 210) * gamma!(202, 203, 210),
    )
    .unwrap();
    let source = coefficient
        .multiply_word(head, &Tensor::word(head, &factors).unwrap())
        .unwrap();
    assert_eq!(source.rank(), 0);
    let variable = target(&[bis!(4, 220)]);
    let selected = factors[1].expression.get_symbol().unwrap();
    for side in [WordDerivativeSide::Left, WordDerivativeSide::Right] {
        let cut = source
            .differentiate_word_with(head, &variable, side, |g| tangent(g, selected, &variable))
            .unwrap();
        assert_eq!(
            cut.structure.logical_slots(),
            variable
                .logical_slots()
                .into_iter()
                .map(|slot| slot.dual())
                .collect::<Vec<_>>()
        );
        let response = Tensor::word(
            head,
            &[leaf("word_derivative::spinor_response", &[bis!(4, 220)])],
        )
        .unwrap();
        let reconstructed = match side {
            WordDerivativeSide::Left => response.multiply_word(head, &cut),
            WordDerivativeSide::Right => cut.multiply_word(head, &response),
        }
        .unwrap();
        let mut expected = factors.clone();
        expected[1] = leaf("word_derivative::spinor_response", &[bis!(4, 201)]);
        let expected = coefficient
            .multiply_word(head, &Tensor::word(head, &expected).unwrap())
            .unwrap();
        assert_eq!(canonical(reconstructed), canonical(expected));
    }
}

#[test]
fn word_derivative_sums_repeated_components_and_preserves_exact_coefficients() {
    let head = head();
    let a = leaf("word_derivative::repeated", &[mink!(4, 301)]);
    let b = leaf("word_derivative::repeated", &[mink!(4, 302)]);
    let word = Tensor::word(head, &[a.clone(), b.clone()]).unwrap();
    let scalar = Tensor::infer(Atom::num((2, 7)) + Atom::i()).unwrap();
    let source = scalar
        .multiply_word(head, &word)
        .unwrap()
        .try_add(&word)
        .unwrap();
    let variable = target(&[mink!(4, 303)]);
    let cut = source
        .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |g| {
            tangent(g, a.expression.get_symbol().unwrap(), &variable)
        })
        .unwrap();
    assert!(!cut.is_zero());
    let response = Tensor::word(
        head,
        &[leaf("word_derivative::repeat_response", &[mink!(4, 303)])],
    )
    .unwrap();
    let reconstructed = response.multiply_word(head, &cut).unwrap();
    let t = symbol!("word_derivative::t"; Scalar);
    let shifted = [a, b].map(|g| {
        let xi = leaf(
            "word_derivative::repeat_response",
            &[composition::port_atom(g.structure.logical_slots()[0])],
        );
        g.try_add(
            &xi.with_algebra_result(Atom::var(t) * xi.expression())
                .unwrap(),
        )
        .unwrap()
    });
    let shifted = Tensor::word(head, &shifted).unwrap();
    let expected = scalar
        .multiply_word(head, &shifted)
        .unwrap()
        .try_add(&shifted)
        .unwrap();
    let expected = source
        .with_algebra_result(linear_coefficient(&expected.expression, t))
        .unwrap();
    assert_eq!(canonical(reconstructed), canonical(expected));
}

#[test]
fn word_derivative_materializes_open_coefficients_and_retains_no_match_zero_shape() {
    let head = head();
    let rep = Minkowski {}.new_rep(4).cast::<LibraryRep>();
    let coefficient = Tensor::from_declared_parts(
        function!(
            SPENSO_TAG.tensor_symbol("word_derivative::even"),
            rep.to_symbolic([])
        ),
        PartialStructure::from_logical_slots([rep.slot(PartialIndex::open(0))]),
    )
    .unwrap();
    let a = leaf("word_derivative::a", &[mink!(4, 401)]);
    let b = leaf("word_derivative::b", &[mink!(4, 401)]);
    let source = coefficient
        .multiply_word(head, &Tensor::word(head, &[a.clone(), b]).unwrap())
        .unwrap();
    let variable = target(&[mink!(4, 402)]);
    let actual = source
        .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |g| {
            tangent(g, a.expression.get_symbol().unwrap(), &variable)
        })
        .unwrap();
    assert_eq!(actual.rank(), 2);
    assert!(actual.structure.open_positions().is_empty());
    assert_eq!(
        actual.structure.logical_slots()[1],
        variable.logical_slots()[0].dual()
    );
    for input in [
        &source,
        &coefficient,
        &source.with_algebra_result(Atom::Zero).unwrap(),
    ] {
        let actual = input
            .differentiate_word_with(head, &variable, WordDerivativeSide::Right, |_| Ok(None))
            .unwrap();
        assert!(actual.is_zero());
        assert_eq!(actual.rank(), 2);
        assert!(actual.structure.open_positions().is_empty());
        assert_eq!(
            actual.structure.logical_slots()[1],
            variable.logical_slots()[0].dual()
        );
    }
    let nilpotent = Tensor::word(
        head,
        &[
            a.clone(),
            a,
            leaf("word_derivative::free", &[mink!(4, 403)]),
        ],
    )
    .unwrap();
    assert!(nilpotent.is_zero());
    let actual = nilpotent
        .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |_| {
            panic!("ranked zero has no generators")
        })
        .unwrap();
    assert!(actual.is_zero());
    assert_eq!(actual.rank(), 2);
}

#[test]
fn word_derivative_accepts_explicit_native_symmetric_and_antisymmetric_tangent_maps() {
    let head = head();
    for (field_head, response_head, kind) in [
        (
            symbol!("word_derivative::symmetric"; Symmetric; tag=&SPENSO_TAG.tensor),
            symbol!("word_derivative::symmetric_response"; Symmetric; tag=&SPENSO_TAG.tensor),
            FactorProjector::Symmetric,
        ),
        (
            symbol!("word_derivative::antisymmetric"; Antisymmetric; tag=&SPENSO_TAG.tensor),
            symbol!("word_derivative::antisymmetric_response"; Antisymmetric; tag=&SPENSO_TAG.tensor),
            FactorProjector::Antisymmetric,
        ),
    ] {
        let field = Tensor::infer(function!(field_head, mink!(4, 501), mink!(4, 502))).unwrap();
        let spectator = leaf("word_derivative::spectator", &[]);
        let source = Tensor::word(head, &[field.clone(), spectator.clone()]).unwrap();
        let variable = target(&[mink!(4, 503), mink!(4, 504)]);
        let rep = Minkowski {}.new_rep(4);
        let expression: Atom = kind
            .permutations(2)
            .unwrap()
            .into_iter()
            .map(|(c, p)| {
                c * rep.id(Atom::num(503), Atom::num(501 + p[0]))
                    * rep.id(Atom::num(504), Atom::num(501 + p[1]))
            })
            .sum();
        let map = Tensor::from_declared_parts(
            expression,
            PartialStructure::from_logical_slots(
                (field.structure.logical_slots())
                    .into_iter()
                    .chain(variable.logical_slots().into_iter().map(|slot| slot.dual())),
            ),
        )
        .unwrap();
        let cut = source
            .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |g| {
                Ok((g.expression.get_symbol() == Some(field_head)).then(|| map.clone()))
            })
            .unwrap();
        let response = Tensor::word(
            head,
            &[Tensor::infer(function!(response_head, mink!(4, 503), mink!(4, 504))).unwrap()],
        )
        .unwrap();
        let actual = response.multiply_word(head, &cut).unwrap();
        let expected = Tensor::word(
            head,
            &[
                Tensor::infer(function!(response_head, mink!(4, 501), mink!(4, 502))).unwrap(),
                spectator,
            ],
        )
        .unwrap();
        assert_eq!(canonical(actual), canonical(expected));
    }
}

#[test]
fn word_derivative_rejects_collisions_wrong_policies_and_unsupported_grammar() {
    let head = head();
    let a = leaf("word_derivative::reject_a", &[mink!(4, 601)]);
    let b = leaf("word_derivative::reject_b", &[mink!(4, 602)]);
    let source = Tensor::word(head, &[a.clone(), b.clone()]).unwrap();
    let variable = target(&[mink!(4, 603)]);
    for invalid in [
        target(&[mink!(4, 601)]),
        PartialStructure::from_logical_slots([Minkowski {}
            .new_rep(4)
            .cast::<LibraryRep>()
            .slot(PartialIndex::open(0))]),
        PartialStructure::from_logical_slots(
            [Minkowski {}
                .new_rep(4)
                .cast::<LibraryRep>()
                .slot(PartialIndex::Explicit(AbstractIndex::Normal(604))); 2],
        ),
    ] {
        assert!(
            source
                .differentiate_word_with(head, &invalid, WordDerivativeSide::Left, |_| panic!(
                    "invalid target must fail before callback"
                ))
                .is_err()
        );
    }
    let closed = Tensor::word(
        head,
        &[
            a.clone(),
            leaf("word_derivative::reject_b", &[mink!(4, 601)]),
        ],
    )
    .unwrap();
    assert!(
        closed
            .differentiate_word_with(
                head,
                &target(&[mink!(4, 601)]),
                WordDerivativeSide::Left,
                |_| panic!("bound collision must fail before callback")
            )
            .is_err()
    );
    let map = tangent(&a, a.expression.get_symbol().unwrap(), &variable)
        .unwrap()
        .unwrap();
    for bad in [
        Tensor::infer(Atom::one()).unwrap(),
        leaf("word_derivative::bad", &[mink!(4, 601), mink!(4, 604)]),
        map.permuted(&[1, 0]).unwrap(),
        Tensor::word(
            head,
            &[leaf(
                "word_derivative::bad",
                &[mink!(4, 601), mink!(4, 603)],
            )],
        )
        .unwrap(),
    ] {
        assert!(
            source
                .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |_| Ok(Some(
                    bad.clone()
                )))
                .is_err()
        );
    }
    assert!(
        matches!(source.differentiate_word_with(head,&variable,WordDerivativeSide::Left,|_|Err(TensorInferenceError::invalid("test policy error"))),Err(TensorInferenceError::Invalid(message)) if message=="test policy error")
    );
    let ordinary = symbol!("word_derivative::ordinary"; Flat,Linear; tag=&SPENSO_TAG.tensor_word);
    let foreign =
        symbol!("word_derivative::foreign"; Flat,Linear,Antisymmetric; tag=&SPENSO_TAG.tensor_word);
    assert!(
        source
            .differentiate_word_with(ordinary, &variable, WordDerivativeSide::Left, |_| Ok(None))
            .is_err()
    );
    let foreign = Tensor::word(foreign, &[a.clone(), b.clone()]).unwrap();
    assert!(
        foreign
            .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |_| Ok(None))
            .is_err()
    );
    let multiple = Tensor::word(head, &[a])
        .unwrap()
        .multiply(&Tensor::word(head, &[b]).unwrap())
        .unwrap();
    assert!(
        multiple
            .differentiate_word_with(head, &variable, WordDerivativeSide::Left, |_| Ok(None))
            .is_err()
    );
}

#[test]
fn word_derivative_rejects_hidden_target_collisions_and_retains_no_match_history() {
    let head = head();
    let a = leaf("word_derivative_history::a", &[mink!(4, 991)]);
    let zero = Tensor::word(head, &[a.clone(), a]).unwrap();
    let nonzero = zero.with_algebra_result(Atom::one()).unwrap();
    let collision = target(&[mink!(4, 991)]);
    let fresh = target(&[mink!(4, 992)]);
    let probe = leaf("word_derivative_history::probe", &[mink!(4, 991)]);
    for source in [zero, nonzero] {
        assert!(
            source
                .differentiate_word_with(head, &collision, WordDerivativeSide::Left, |_| panic!(
                    "target freshness must precede callbacks"
                ))
                .is_err()
        );
        let derivative = source
            .differentiate_word_with(head, &fresh, WordDerivativeSide::Left, |_| Ok(None))
            .unwrap();
        assert!(derivative.is_zero());
        assert!(matches!(
            derivative.multiply_word(head, &probe),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
    }
}
