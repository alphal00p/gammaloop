use crate::{
    Cookable, IndexTooling,
    color::CS,
    representations::{ColorAdjoint, ColorFundamental},
    test_support::test_initialize,
};
use spenso::{antisym, bracket, broadcast_symbol, cyclic, dind, euc, lor, mink, sym, tensor};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function, parse_lit, symbol,
};

use spenso::structure::{
    abstract_index::AbstractIndex, representation::RepName, slot::IsAbstractSlot,
};

fn canonicalize(expression: &Atom) -> Atom {
    expression
        .canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .expect("test expression should canonicalize")
}

fn canonicalize_formal(expression: &Atom) -> Atom {
    let cooking = crate::CookSettings::reversible().with_representation_payloads(true, true);
    let admitted = cooking.try_cook(expression.as_view()).unwrap();
    let canonical = canonicalize(&admitted);
    cooking.uncook(canonical.as_view())
}

#[test]
fn formal_color_dimensions_canonicalize_closed_quark_trace_partners() {
    test_initialize();
    let nc = Atom::var(CS.nc);
    let adjoint_dimension = nc.clone().pow(2) - Atom::one();
    let fundamental = ColorFundamental {}.to_symbolic([nc]);
    let trace = |indices: [usize; 6]| {
        spenso::shadowing::trace(
            fundamental.clone(),
            indices.map(|index| {
                CS.chain_t(ColorAdjoint {}.to_symbolic([
                    adjoint_dimension.clone(),
                    function!(symbol!("formal_trace_hedge"), index),
                ]))
            }),
        )
        .cook_indices()
    };

    // Actual sampled color traces of XS GL204/205 (also GL206/207) and GL282/283.
    // The two contracted adjoint names may be exchanged without changing T slots.
    for (left, right) in [
        ([6, 10, 12, 6, 10, 12], [6, 12, 10, 6, 12, 10]),
        ([6, 14, 10, 6, 14, 10], [6, 10, 14, 6, 10, 14]),
    ] {
        let canonical = canonicalize_formal(&trace(left));
        assert_eq!(canonical, canonicalize_formal(&trace(right)));
        assert_eq!(canonical, canonicalize_formal(&canonical));
        assert!(canonical.contains(&adjoint_dimension));
    }
}

#[test]
fn formal_color_dimensions_preserve_free_adjoint_and_fundamental_slots() {
    test_initialize();
    let open_adjoint = spenso::shadowing::trace(
        parse_lit!(cof(Nc), default_namespace = "spenso"),
        [
            parse_lit!(
                coad(Nc ^ 2 - 1, formal_free_a),
                default_namespace = "spenso"
            ),
            parse_lit!(
                coad(Nc ^ 2 - 1, formal_free_b),
                default_namespace = "spenso"
            ),
            parse_lit!(
                coad(Nc ^ 2 - 1, formal_free_c),
                default_namespace = "spenso"
            ),
        ]
        .map(|adjoint| CS.chain_t(adjoint)),
    );
    let renamed_adjoint = open_adjoint
        .replace(parse_lit!(formal_free_b, default_namespace = "spenso"))
        .with(parse_lit!(formal_free_d, default_namespace = "spenso"));
    let canonical = canonicalize_formal(&open_adjoint);
    assert_ne!(canonical, canonicalize_formal(&renamed_adjoint));
    assert!(canonical.contains(&parse_lit!(
        coad(Nc ^ 2 - 1, formal_free_b),
        default_namespace = "spenso"
    )));

    let open_fundamental = parse_lit!(
        t(
            coad(Nc ^ 2 - 1, formal_paired_a),
            cof(Nc, formal_free_i),
            dind(cof(Nc, formal_paired_j))
        ) * t(
            coad(Nc ^ 2 - 1, formal_paired_a),
            cof(Nc, formal_paired_j),
            dind(cof(Nc, formal_free_k))
        ),
        default_namespace = "spenso"
    );
    let renamed_fundamental = open_fundamental
        .replace(parse_lit!(formal_free_i, default_namespace = "spenso"))
        .with(parse_lit!(formal_free_l, default_namespace = "spenso"));
    let canonical = canonicalize_formal(&open_fundamental);
    assert_ne!(canonical, canonicalize_formal(&renamed_fundamental));
    for external in [
        parse_lit!(cof(Nc, formal_free_i), default_namespace = "spenso"),
        parse_lit!(dind(cof(Nc, formal_free_k)), default_namespace = "spenso"),
    ] {
        assert!(canonical.contains(&external));
    }
    assert_eq!(canonical, canonicalize_formal(&canonical));
}

fn odd_cycle() -> Atom {
    let a = mink!(4, signed_cycle_a);
    let b = mink!(4, signed_cycle_b);
    let c = mink!(4, signed_cycle_c);
    antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a)
}

#[test]
fn odd_antisymmetric_cycle_canonicalizes_to_zero() {
    test_initialize();
    assert!(canonicalize(&odd_cycle()).is_zero());
}

#[test]
fn canonicalization_preserves_spectators_after_pruning() {
    test_initialize();
    assert!(canonicalize(&odd_cycle()).is_zero());
    let spectator = tensor!(signed_metric_spectator_a) + tensor!(signed_metric_spectator_b);
    let expression = &spectator * (Atom::one() + odd_cycle());
    let canonical = canonicalize(&expression);
    assert_eq!(canonical, spectator);
    assert_eq!(canonicalize(&canonical), canonical);
}

#[test]
fn nested_odd_summand_is_removed() {
    test_initialize();
    let scalar = tensor!(signed_nested_scalar);
    let expression = &scalar * (Atom::num(1) + odd_cycle());

    assert_eq!(canonicalize(&expression), scalar);
}

#[test]
fn pruning_nested_zero_preserves_unrelated_sum_factors() {
    test_initialize();
    let spectator = tensor!(signed_spectator_a) + tensor!(signed_spectator_b);
    let expected = tensor!(signed_nested_scalar) * spectator;
    let expression = &expected * (Atom::one() + odd_cycle());

    assert_eq!(canonicalize(&expression), expected);
}

#[test]
fn products_of_tensor_sums_keep_closed_indices_and_factorization() {
    test_initialize();
    let expression = |i: Atom, j: Atom| {
        (tensor!(factorized_a, i.clone()) + tensor!(factorized_b, i.clone()))
            * (tensor!(factorized_c, i.clone()) + tensor!(factorized_d, i))
            + tensor!(factorized_e, j.clone()) * tensor!(factorized_f, j)
    };
    let original = expression(mink!(4, factorized_i), mink!(4, factorized_j));
    let renamed = expression(mink!(4, factorized_l), mink!(4, factorized_k));
    let canonical = canonicalize(&original);

    assert_eq!(canonical, canonicalize(&renamed));
    assert_eq!(canonical, canonicalize(&canonical));
    assert!(
        canonical
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .is_empty()
    );
    let AtomView::Add(sum) = canonical.as_view() else {
        panic!("the two original scalar summands must remain present");
    };
    assert!(sum.iter().any(|term| matches!(term, AtomView::Mul(product)
        if product.iter().filter(|factor| matches!(factor, AtomView::Add(_))).count() == 2)));
    assert_ne!(canonical, canonicalize(&(original + Atom::one())));
}

#[test]
fn products_of_tensor_sums_preserve_genuine_external_indices() {
    test_initialize();
    // Native allocation must reserve true external names even when they sort
    // before all dummy candidates, as well as when they sort after them.
    for external in [
        mink!(4, a_factorized_external),
        mink!(4, z_factorized_external),
    ] {
        let expression = |i: Atom| {
            (tensor!(factorized_open_a, i.clone()) + tensor!(factorized_open_b, i.clone()))
                * (tensor!(factorized_open_c, i.clone(), external.clone())
                    + tensor!(factorized_open_d, i, external.clone()))
        };
        let original = expression(mink!(4, factorized_open_i));
        let renamed = expression(mink!(4, factorized_open_j));
        let canonical = canonicalize(&original);

        assert_eq!(canonical, canonicalize(&renamed));
        assert_eq!(canonical, canonicalize(&canonical));
        assert_eq!(
            canonical.list_dangling::<AbstractIndex>().unwrap(),
            vec![external]
        );
        assert!(matches!(canonical.as_view(), AtomView::Mul(product)
            if product.iter().filter(|factor| matches!(factor, AtomView::Add(_))).count() == 2));
    }
}

#[test]
fn factorized_canonicalization_rejects_mismatched_external_indices() {
    test_initialize();
    let expression = tensor!(factorized_bad_a, mink!(4, factorized_bad_i))
        + tensor!(factorized_bad_b, mink!(4, factorized_bad_j));

    assert!(std::panic::catch_unwind(|| canonicalize(&expression)).is_err());
}

#[test]
fn factorized_canonicalization_rejects_overcontracted_indices() {
    test_initialize();
    let index = mink!(4, factorized_overcontracted);
    let expression = tensor!(factorized_over_a, index.clone())
        * tensor!(factorized_over_b, index.clone())
        * tensor!(factorized_over_c, index);

    assert!(std::panic::catch_unwind(|| canonicalize(&expression)).is_err());
}

#[test]
fn factored_expression_is_preserved_when_no_term_is_removed() {
    test_initialize();
    let a = mink!(4, signed_factored_a);
    let b = mink!(4, signed_factored_b);
    let scalar = tensor!(signed_factored_scalar);
    let expression = scalar * (Atom::num(1) + antisym!(a, b).pow(2));
    let canonical = canonicalize(&expression);

    assert!(!canonical.is_zero());
    assert!(matches!(canonical.as_view(), AtomView::Mul(product)
        if product.iter().any(|factor| matches!(factor, AtomView::Add(_)))));
}

#[test]
fn sum_inside_function_remains_opaque() {
    test_initialize();
    let expression = tensor!(signed_opaque, Atom::num(1) + odd_cycle());
    let canonical = canonicalize(&expression);

    assert!(!canonical.is_zero());
    assert!(matches!(canonical.as_view(), AtomView::Fun(function)
        if function.iter().any(|argument| matches!(argument, AtomView::Add(_)))));
}

#[test]
fn broadcast_boundary_does_not_hide_independent_zero() {
    test_initialize();
    let a = mink!(4, signed_opaque_a);
    let b = mink!(4, signed_opaque_b);
    let c = mink!(4, signed_opaque_c);
    let inner = antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a);
    let opaque = function!(
        broadcast_symbol!(signed_opaque_broadcast),
        Atom::num(1) + inner
    );

    // A nonlinear function must retain the vanishing summand rather than
    // applying its zero before evaluating the enclosing function.
    let canonical = canonicalize(&opaque);
    assert!(!canonical.is_zero());
    assert!(matches!(canonical.as_view(), AtomView::Fun(function)
        if function.iter().any(|argument| matches!(argument, AtomView::Add(_)))));

    // The opaque function boundary must not suppress an independent zero.
    assert!(canonicalize(&(odd_cycle() * opaque)).is_zero());
}

#[test]
fn even_antisymmetric_automorphism_survives_canonicalization() {
    test_initialize();
    let a = mink!(4, signed_square_a);
    let b = mink!(4, signed_square_b);

    assert!(!canonicalize(&antisym!(a, b).pow(2)).is_zero());
}

#[test]
fn power_copies_closed_components_as_a_whole() {
    test_initialize();

    // Keep the closed contraction as one Power child until Spenso parses it.
    assert!(canonicalize(&bracket!(odd_cycle()).pow(2)).is_zero());
}

#[test]
fn external_indices_prevent_slot_automorphism() {
    test_initialize();
    let a = mink!(4, signed_external_a);
    let b = mink!(4, signed_external_b);

    assert!(!canonicalize(&antisym!(a, b)).is_zero());
}

#[test]
fn ordered_slots_prevent_antisymmetric_automorphism() {
    test_initialize();
    let a = mink!(4, signed_ordered_a);
    let b = mink!(4, signed_ordered_b);
    let contraction = tensor!(signed_ordered, a.clone(), b.clone()) * antisym!(a, b);

    assert!(!canonicalize(&contraction).is_zero());
}

#[test]
fn symmetric_slots_allow_antisymmetric_automorphism() {
    test_initialize();
    let a = mink!(4, signed_symmetric_a);
    let b = mink!(4, signed_symmetric_b);
    let contraction = sym!(a.clone(), b.clone()) * antisym!(a, b);

    assert!(canonicalize(&contraction).is_zero());
}

#[test]
fn different_representation_groups_break_odd_automorphism() {
    test_initialize();
    let a = mink!(4, signed_group_a);
    let b = euc!(4, signed_group_b);
    let c = lor!(4, signed_group_c);
    let contraction =
        antisym!(a.clone(), b.clone()) * antisym!(b, c.clone()) * antisym!(dind!(c), a);

    assert!(!canonicalize(&contraction).is_zero());
}

#[test]
fn cyclic_symmetry_preserves_only_rotations() {
    test_initialize();
    let a = mink!(4, signed_cyclic_a);
    let b = mink!(4, signed_cyclic_b);
    let c = mink!(4, signed_cyclic_c);
    let rank_three =
        cyclic!(a.clone(), b.clone(), c.clone()) * antisym!(a.clone(), b.clone(), c.clone());

    // A three-cycle is even, so the permitted cyclic rotation has no sign.
    assert!(!canonicalize(&rank_three).is_zero());

    let d = mink!(4, signed_cyclic_d);
    let rank_four = cyclic!(a.clone(), b.clone(), c.clone(), d.clone()) * antisym!(a, b, c, d);

    // A four-cycle is odd. Reflections are not part of cyclic symmetry.
    assert!(canonicalize(&rank_four).is_zero());
}

#[test]
fn inverse_scalar_contractions_preserve_their_boundary_during_canonicalization() {
    test_initialize();
    let p = spenso::vector_symbol!("inverse_scalar_canon::p");
    let q = spenso::vector_symbol!("inverse_scalar_canon::q");
    let minkowski = spenso::structure::representation::Minkowski {};
    use spenso::structure::representation::RepName;
    let denominator = minkowski.new_rep(4).inner_product(p, q) - Atom::one();
    let numerator = tensor!(inverse_scalar_a, mink!(4, inverse_scalar_i))
        * tensor!(inverse_scalar_b, mink!(4, inverse_scalar_i));
    let expression = numerator / denominator;
    let canonical = canonicalize(&expression);
    assert_eq!(canonical, canonicalize(&canonical));
    assert!(
        canonical
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .is_empty()
    );
    assert!(!canonical.to_plain_string().contains("scalar_store"));
}

#[test]
fn opaque_scalar_payloads_do_not_share_outer_dummy_indices() {
    test_initialize();
    let metadata = symbol!("opaque_scalar_canon::metadata"; Scalar);
    let a = spenso::vector_symbol!("opaque_scalar_canon::a");
    let b = spenso::vector_symbol!("opaque_scalar_canon::b");
    let c = spenso::vector_symbol!("opaque_scalar_canon::c");
    let d = spenso::vector_symbol!("opaque_scalar_canon::d");
    let i = mink!(4, 701);
    let payload = function!(a, &i) * function!(b, &i);
    let scalar = function!(metadata, &payload);
    let outer = function!(c, &i) * function!(d, &i);
    let canonical = canonicalize(&(&scalar * &outer));
    assert_eq!(canonical, scalar * canonicalize(&outer));
    assert_eq!(canonical, canonicalize(&canonical));
    assert!(!canonical.contains_symbol(spenso::network::tags::SPENSO_TAG.scalar));
}

#[test]
fn nonlinear_scalar_payloads_preserve_their_literal_binding_on_repeated_canon() {
    test_initialize();
    let a = spenso::vector_symbol!("nonlinear_scalar_canon::a");
    let b = spenso::vector_symbol!("nonlinear_scalar_canon::b");
    let c = spenso::vector_symbol!("nonlinear_scalar_canon::c");
    let d = spenso::vector_symbol!("nonlinear_scalar_canon::d");
    let i = mink!(4, 702);
    let payload = function!(a, &i) * function!(b, &i);
    let outer = function!(c, &i) * function!(d, &i);
    for payload in [payload.clone(), canonicalize(&payload)] {
        let scalar = function!(symbolica::atom::Symbol::COS, payload);
        let canonical = canonicalize(&(&scalar * &outer));
        assert_eq!(canonical, scalar * canonicalize(&outer));
        assert_eq!(canonical, canonicalize(&canonical));
        assert!(!canonical.contains_symbol(spenso::network::tags::SPENSO_TAG.scalar));
    }
}

#[test]
fn opaque_scalar_hidden_powers_are_unchanged_by_outer_power_preparation() {
    use crate::tensor::SymbolicTensor;
    use spenso::structure::partial::PartialStructure;
    test_initialize();
    let metadata = symbol!("hidden_power_canon::metadata"; Scalar);
    let i = mink!(4, 709);
    let inner = SymbolicTensor::<PartialStructure>::infer(
        tensor!(hidden_power_a, i.clone()) * tensor!(hidden_power_b, &i),
    )
    .unwrap();
    let outer = tensor!(hidden_power_c, &i) * tensor!(hidden_power_d, i);
    for exponent in [Atom::num(2), Atom::num(-1), Atom::num((1, 2))] {
        let exponent = SymbolicTensor::infer(exponent).unwrap();
        let power = inner.pow(&exponent).unwrap().into_expression();
        for head in [metadata, symbolica::atom::Symbol::COS] {
            let payload = function!(head, &power);
            let canonical = canonicalize(&(&payload * &outer));
            assert_eq!(canonical, payload * canonicalize(&outer));
            assert_eq!(canonical, canonicalize(&canonical));
        }
    }
}

#[test]
fn scalar_leaf_aliases_do_not_expose_temporary_metadata_to_tensor_hooks() {
    use crate::tensor::SymbolicTensor;
    use spenso::structure::partial::PartialStructure;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    test_initialize();
    let metadata = symbol!("scalar_roles_canon::metadata"; Scalar);
    let i = mink!(4, 710);
    let inner = SymbolicTensor::<PartialStructure>::infer(
        tensor!(scalar_roles_inner_a, &i) * tensor!(scalar_roles_inner_b, &i),
    )
    .unwrap();
    for (case, exponent) in [None, Some(2), Some(-1)].into_iter().enumerate() {
        let payload = if let Some(exponent) = exponent {
            let exponent = SymbolicTensor::infer(Atom::num(exponent)).unwrap();
            function!(metadata, inner.pow(&exponent).unwrap().expression())
        } else {
            function!(metadata, symbol!("scalar_roles_canon::x"))
        };
        let altered = Arc::new(AtomicUsize::new(0));
        let observed = altered.clone();
        let expected_metadata = payload.clone();
        let head = spenso::tensor_symbol!(
            format!("scalar_roles_canon::hooked_{case}"),
            norm = move |value, output| {
                if !matches!(value, symbolica::atom::AtomView::Fun(function)
                    if function.iter().next() == Some(expected_metadata.as_view()))
                {
                    observed.fetch_add(1, Ordering::SeqCst);
                    output.to_num(0);
                }
            }
        );
        let outer = function!(head, &payload, &i) * tensor!(scalar_roles_b, &i);
        let source = &payload * &outer;
        let canonical = canonicalize(&source);
        assert_eq!(altered.load(Ordering::SeqCst), 0, "case {case}");
        assert!(!canonical.is_zero(), "case {case}");
        assert_eq!(canonical, payload * canonicalize(&outer), "case {case}");
        assert_eq!(canonical, canonicalize(&canonical), "case {case}");
        assert_eq!(altered.load(Ordering::SeqCst), 0, "case {case}");
    }
}

#[test]
fn opaque_scalar_aliases_do_not_capture_literal_native_handles() {
    test_initialize();
    let metadata = symbol!("literal_alias_canon::metadata"; Scalar);
    let i = mink!(4, 703);
    let outer = tensor!(literal_alias_a, &i) * tensor!(literal_alias_b, &i);
    let native_handle = function!(spenso::network::tags::SPENSO_TAG.scalar, 0);
    let canonical_name = Atom::var(symbol!("idenso::canonical_scalar_0"));
    let payload = function!(metadata, native_handle, canonical_name);
    let canonical = canonicalize(&(&payload * &outer));
    assert_eq!(canonical, payload * canonicalize(&outer));
    assert_eq!(canonical, canonicalize(&canonical));
}

#[cfg(not(target_arch = "wasm32"))]
#[test]
fn canonical_scalar_aliases_skip_preexisting_global_heads() {
    const CHILD: &str = "IDENSO_CANONICAL_ALIAS_FRESHNESS_CHILD";
    if std::env::var_os(CHILD).is_none() {
        for mode in ["opaque", "inverse"] {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "tensor::tests::canonicalization::canonical_scalar_aliases_skip_preexisting_global_heads",
                    "--nocapture",
                    "--test-threads=1",
                ])
                .env(CHILD, mode)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{mode}: {}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
        }
        return;
    }
    use symbolica::{
        atom::{Symbol, SymbolAttribute, SymbolBuilder},
        wrap_symbol,
    };
    test_initialize();
    // The isolated process starts the allocator at zero, so this exercises an
    // occupied candidate without relying on parallel tests' global counters.
    let occupied = SymbolBuilder::build_group(vec![
        SymbolBuilder::new(wrap_symbol!("idenso::canonical_scalar_0"))
            .with_attributes(vec![SymbolAttribute::Real, SymbolAttribute::Positive])
            .with_normalization_function(|_, output| {
                output.to_num(0);
            }),
    ])
    .unwrap()[0];
    let metadata = symbol!("fresh_alias_canon::metadata"; Scalar);
    let payload = function!(metadata, symbol!("fresh_alias_canon::x"));
    let i = mink!(4, 708);
    let outer = tensor!(fresh_alias_a, i.clone()) * tensor!(fresh_alias_b, i);
    let (source, expected) = if std::env::var(CHILD).unwrap() == "inverse" {
        use crate::tensor::SymbolicTensor;
        use spenso::structure::partial::PartialStructure;
        let exponent = SymbolicTensor::<PartialStructure>::infer(Atom::num(-1)).unwrap();
        let source = SymbolicTensor::<PartialStructure>::infer(outer.clone())
            .unwrap()
            .pow(&exponent)
            .unwrap()
            .into_expression();
        let expected = SymbolicTensor::<PartialStructure>::infer(canonicalize(&outer))
            .unwrap()
            .pow(&exponent)
            .unwrap()
            .into_expression();
        (source, expected)
    } else {
        (&payload * &outer, payload * canonicalize(&outer))
    };
    assert!(!source.contains_symbol(occupied));
    let canonical = canonicalize(&source);
    assert_eq!(canonical, expected);
    assert_eq!(canonical, canonicalize(&canonical));
    let fresh = Symbol::get_symbol(wrap_symbol!("idenso::canonical_scalar_1")).unwrap();
    assert!(!fresh.is_real());
    assert!(!fresh.is_positive());
    assert!(fresh.get_normalization_function().is_none());
    assert!(fresh.get_tags().is_empty());
}

#[test]
fn equal_composite_scalar_coefficients_keep_signed_tensor_cancellation() {
    test_initialize();
    let metadata = symbol!("signed_alias_canon::metadata"; Scalar);
    let coefficient = function!(metadata, symbol!("signed_alias_canon::x"));
    let i = mink!(4, 704);
    let j = mink!(4, 705);
    let first = tensor!(signed_alias_a, i.clone()) * tensor!(signed_alias_b, i);
    let second = tensor!(signed_alias_a, j.clone()) * tensor!(signed_alias_b, j);
    assert!(canonicalize(&(coefficient * (first - second))).is_zero());
}

#[test]
fn opaque_scalar_aliasing_retains_independent_inverse_and_fractional_base_canon() {
    use crate::tensor::SymbolicTensor;
    use spenso::structure::partial::PartialStructure;
    test_initialize();
    let i = mink!(4, 706);
    let j = mink!(4, 707);
    let first = SymbolicTensor::<PartialStructure>::infer(
        tensor!(fractional_alias_a, i.clone()) * tensor!(fractional_alias_b, i),
    )
    .unwrap();
    let second = SymbolicTensor::<PartialStructure>::infer(
        tensor!(fractional_alias_a, j.clone()) * tensor!(fractional_alias_b, j),
    )
    .unwrap();
    for exponent in [Atom::num(-1), Atom::num((1, 2)), Atom::num((-1, 2))] {
        let exponent = SymbolicTensor::infer(exponent).unwrap();
        let canonical = canonicalize(first.pow(&exponent).unwrap().expression());
        assert_eq!(
            canonical,
            canonicalize(second.pow(&exponent).unwrap().expression())
        );
        assert_eq!(canonical, canonicalize(&canonical));
    }
}

#[test]
fn pruning_zero_arm_exposes_outer_signed_contraction() {
    test_initialize();
    let a = mink!(4, signed_exposed_a);
    let b = mink!(4, signed_exposed_b);
    let sum = sym!(a.clone(), b.clone())
        + odd_cycle() * tensor!(signed_exposed_other, a.clone(), b.clone());
    let expression = antisym!(a, b) * sum;
    assert!(canonicalize(&expression).is_zero());
}

#[test]
fn pruning_open_zero_arm_preserves_free_and_dual_ports() {
    test_initialize();
    let a = lor!(4, signed_pruned_dual_a);
    let b = mink!(4, signed_pruned_free_b);
    let survivor = tensor!(signed_pruned_survivor, a.clone(), b.clone());
    let sum = &survivor + odd_cycle() * tensor!(signed_pruned_zero, a.clone(), b.clone());
    let closing = tensor!(signed_pruned_close, dind!(a));
    let expected = canonicalize(&(&closing * survivor));
    let result = canonicalize(&(closing * sum));
    assert_eq!(result, expected);
    assert_eq!(canonicalize(&result), result);
    assert_eq!(result.list_dangling::<AbstractIndex>().unwrap(), vec![b]);
}

#[test]
fn pruning_signed_term_does_not_mutate_a_shared_tensor_occurrence() {
    test_initialize();
    let a = mink!(4, signed_cycle_a);
    let b = mink!(4, signed_cycle_b);
    let even = antisym!(a, b).pow(2);
    let expected = canonicalize(&even);
    assert!(!expected.is_zero());
    assert_eq!(canonicalize(&(odd_cycle() + even)), expected);
}

#[test]
fn unparseable_indices_are_rejected_before_signed_zero_pruning() {
    test_initialize();
    use crate::{CanonicalizationError, Cookable, NetworkToolingError};
    use symbolica::symbol;
    let slot = |n| mink!(4, function!(symbol!("signed_raw_index"), Atom::num(n)));
    let expression =
        antisym!(slot(1), slot(2)) * antisym!(slot(2), slot(3)) * antisym!(slot(3), slot(1));
    assert!(!expression.is_zero());
    let raw = expression.canonize::<AbstractIndex>(AbstractIndex::Dummy);
    assert!(
        matches!(
            raw,
            Err(CanonicalizationError::Network(
                NetworkToolingError::Parse { .. }
            ))
        ),
        "structured index admission result: {raw:?}"
    );
    // A caller with structured labels explicitly admits them once. Zero and
    // nonzero terms obey the same abstract-index boundary.
    assert!(canonicalize(&expression.cook_indices()).is_zero());
}

#[test]
fn real_outer_power_preserves_scalar_metadata_without_hook_replay() {
    use crate::tensor::SymbolicTensor;
    use spenso::structure::{partial::PartialStructure, representation::Minkowski};
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    test_initialize();
    let scope = symbol!("real_power_metadata_canon::scope");
    let private_index = Minkowski {}
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(1701).scoped(scope))
        .to_atom();
    let outer_index = mink!(4, 1702);
    let metadata = function!(
        symbol!("real_power_metadata_canon::metadata"; Scalar),
        tensor!(real_power_metadata_a, &private_index)
            * tensor!(real_power_metadata_b, &private_index)
    );
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = calls.clone();
    let expected = metadata.clone();
    let head = spenso::tensor_symbol!(
        "real_power_metadata_canon::T",
        norm = move |value, output| {
            if !matches!(value, AtomView::Fun(function) if function.iter().next() == Some(expected.as_view()))
            {
                observed.fetch_add(1, Ordering::SeqCst);
                output.to_num(0);
            }
        }
    );
    let source =
        SymbolicTensor::<PartialStructure>::infer(function!(head, &metadata, &outer_index))
            .unwrap();
    let power = source
        .pow(&SymbolicTensor::infer(Atom::num(2)).unwrap())
        .unwrap();
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    let canonical = power.canonize().unwrap();
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    assert!(!canonical.expression().is_zero());
    assert!(canonical.expression().contains(&metadata));
    assert_eq!(
        canonical.expression(),
        canonical.canonize().unwrap().expression()
    );
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

#[test]
fn exterior_word_scalar_coefficients_keep_independent_inner_bindings() {
    use crate::tensor::SymbolicTensor;
    use spenso::{network::tags::SPENSO_TAG, structure::partial::PartialStructure};
    test_initialize();
    type Tensor = SymbolicTensor<PartialStructure>;
    let i = mink!(4, 1801);
    let j = mink!(4, 1802);
    let a = Tensor::infer(tensor!(word_scalar_canon_a, &i)).unwrap();
    let b = Tensor::infer(tensor!(word_scalar_canon_b, &i)).unwrap();
    let head = symbol!("word_scalar_canon::exterior"; Flat, Linear, Antisymmetric; tag=&SPENSO_TAG.tensor_word);
    let word = Tensor::word(head, &[a, b]).unwrap();
    let payload = function!(
        symbolica::atom::Symbol::COS,
        tensor!(word_scalar_canon_x, &i) * tensor!(word_scalar_canon_y, &i)
    );
    let weighted = Tensor::infer(&payload * word.expression()).unwrap();
    assert_eq!(weighted.expression(), &(&payload * word.expression()));
    let outer =
        Tensor::infer(tensor!(word_scalar_canon_c, &j) * tensor!(word_scalar_canon_d, &j)).unwrap();
    let source = weighted.multiply(&outer).unwrap();
    let canonical = source.canonize().unwrap();
    assert!(canonical.expression().contains(&payload));
    assert_eq!(
        canonical.expression(),
        canonical.canonize().unwrap().expression()
    );
    let expected = &payload
        * word
            .multiply(&outer)
            .unwrap()
            .canonize()
            .unwrap()
            .expression();
    assert_eq!(canonical.expression(), &expected);
}

#[test]
fn raw_and_typed_non_scalar_inverse_powers_remain_rejected() {
    use crate::tensor::SymbolicTensor;
    use spenso::structure::partial::PartialStructure;
    test_initialize();
    let source = tensor!(open_inverse_canon_a, mink!(4, 1901));
    let typed = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
    for exponent in [Atom::num(-1), Atom::num(-2), Atom::num((1, 2))] {
        assert!(
            typed
                .pow(&SymbolicTensor::infer(exponent.clone()).unwrap())
                .is_err()
        );
        let error = source
            .clone()
            .pow(exponent)
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap_err();
        assert!(error.to_string().contains("non-scalar"), "{error}");
    }
    let square = typed
        .pow(&SymbolicTensor::infer(Atom::num(2)).unwrap())
        .unwrap();
    assert!(square.canonize().is_ok());
}

#[test]
fn operation_power_preparation_keeps_generic_index_storage() {
    use spenso::structure::{
        abstract_index::AbstractIndexError,
        representation::Minkowski,
        slot::{AbsInd, DummyAind, ParseableAind},
    };
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    struct WrappedIndex(AbstractIndex);
    impl std::fmt::Display for WrappedIndex {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            std::fmt::Display::fmt(&self.0, f)
        }
    }
    impl AbsInd for WrappedIndex {}
    impl ParseableAind for WrappedIndex {
        type Error = AbstractIndexError;
        fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
            AbstractIndex::from_view(view).map(Self)
        }
        fn to_atom(&self) -> Atom {
            self.0.to_atom()
        }
    }
    impl DummyAind for WrappedIndex {
        fn new_dummy() -> Self {
            Self(AbstractIndex::new_dummy())
        }
        fn new_dummy_at(i: usize) -> Self {
            Self(AbstractIndex::new_dummy_at(i))
        }
        fn is_dummy(&self) -> bool {
            self.0.is_dummy()
        }
    }
    test_initialize();
    let scope = symbol!("generic_power_canon::scope");
    let index = Minkowski {}
        .new_rep(4)
        .slot::<WrappedIndex, _>(WrappedIndex(AbstractIndex::Normal(2001).scoped(scope)))
        .to_atom();
    let metadata = function!(
        symbol!("generic_power_canon::metadata"; Scalar),
        tensor!(generic_power_canon_x, &index) * tensor!(generic_power_canon_y, &index)
    );
    let source = tensor!(generic_power_canon_a, &metadata, &index).pow(2);
    let canonical = source
        .canonize::<WrappedIndex>(|i| WrappedIndex(AbstractIndex::Dummy(i)))
        .unwrap();
    assert_eq!(canonical, canonicalize(&source));
    assert_eq!(
        canonical,
        canonical
            .canonize::<WrappedIndex>(|i| WrappedIndex(AbstractIndex::Dummy(i)))
            .unwrap()
    );
}
