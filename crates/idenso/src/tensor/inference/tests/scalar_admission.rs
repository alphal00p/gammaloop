use super::*;
use crate::color::CS;
use crate::test_support;
use spenso::network::tags;
use spenso::{bracket, dot, tensor_symbol, vector_symbol};
use symbolica::atom::{EvaluationInfo, SymbolBuilder};
use symbolica::state::State;
use symbolica::{function, parse, symbol, wrap_symbol};

/// General admission, bypassing the plain-scalar shortcut.
fn general(atom: &Atom) -> InferenceResult<(SymbolicTensor<PartialStructure>, bool)> {
    SymbolicTensor::<PartialStructure>::admit_with_contractions(atom.clone())
}

/// `infer` must publish exactly what general admission publishes: payload
/// bytes, logical interface, carrier flags, the contraction report, every
/// proof cell, and the same error for invalid syntax.
#[track_caller]
fn assert_admits_like_general(atom: Atom, plain: bool) {
    assert_eq!(
        InterfaceInference::is_plain_scalar(atom.as_view()),
        plain,
        "{atom}"
    );
    let fast = SymbolicTensor::<PartialStructure>::infer_with_contractions(atom.clone());
    match (fast, general(&atom)) {
        (Ok((mut fast, fast_contracted)), Ok((general, general_contracted))) => {
            assert_eq!(
                fast.expression.as_view().get_data(),
                general.expression.as_view().get_data(),
                "{atom}"
            );
            assert_eq!(fast.structure, general.structure, "{atom}");
            assert_eq!(
                format!("{:?}", fast.structure),
                format!("{:?}", general.structure),
                "{atom}"
            );
            assert_eq!(
                (fast.is_metric, fast.is_composite, fast_contracted),
                (general.is_metric, general.is_composite, general_contracted),
                "{atom}"
            );
            assert_eq!(
                fast.explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>(),
                general
                    .explicit_index_occurrences()
                    .unwrap()
                    .oriented_entries()
                    .collect::<Vec<_>>(),
                "{atom}"
            );
            // Independent admissions allocate independent history Arcs. After
            // comparing contents, share that identity to compare every cache.
            fast.proofs.occurrences = general.proofs.occurrences.clone();
            assert!(fast.same_reduction_state(&general), "{atom}");
            if plain {
                // The shortcut only claims inputs that admission fully certifies.
                assert!(general.is_scalar(), "{atom}");
                for proof in [
                    &general.proofs.validated,
                    &general.proofs.intrinsic,
                    &general.proofs.rewrite,
                    &general.proofs.algebra,
                ] {
                    assert_eq!(proof.get(), Some(&true), "{atom}");
                }
            }
        }
        (Err(fast), Err(general)) => {
            assert_eq!(format!("{fast:?}"), format!("{general:?}"), "{atom}");
        }
        (fast, general) => panic!("{atom}: {fast:?} != {general:?}"),
    }
}

#[test]
fn plain_scalars_skip_admission_with_identical_results() {
    test_support::test_initialize();
    let namespace = |source: &str| parse!(source, default_namespace = "scalar_admission");
    let x = namespace("x");
    let user = symbol!("scalar_admission::f");
    let parameters = namespace("-𝑖*GC_11*ee^2*cw/sw");
    for atom in [
        Atom::Zero,
        Atom::num(1),
        Atom::num(-7),
        parse!("2^100"),
        Atom::num(Rational::new(3, 7)),
        parse!("1.25"),
        parse!("3.0e-2"),
        Atom::i(),
        parse!("2+3𝑖"),
        x.clone(),
        Atom::var(Symbol::PI),
        Atom::var(Symbol::E),
        parameters.clone(),
        parameters.pow(2) / (&x + 3),
        namespace("(x+y)^2*z^(-1/2)+3"),
        namespace("(a+b)^(c+d)"),
        namespace("sqrt(x)*cos(y)+exp(log(z))*abs(x)*conj(y)*sin(z)"),
        function!(user, function!(user, &x), namespace("cos(y)")),
        function!(symbol!("scalar_admission::s"; Scalar), &x, 2),
        function!(symbol!("scalar_admission::sym"; Symmetric), &parameters, &x),
        function!(
            symbol!("scalar_admission::anti"; Antisymmetric),
            &x,
            &parameters
        ),
        function!(
            symbol!("scalar_admission::cyc"; Cyclesymmetric),
            &x,
            2,
            &parameters
        ),
        function!(
            symbol!("scalar_admission::lin"; Linear),
            &x + 2 * &parameters
        ),
        FunctionBuilder::new(user).finish(),
        function!(Symbol::COEFF, &x + 1),
        &parameters * function!(user, &x),
        // Spenso-owned variables other than chain placeholders are opaque.
        parse!("z", default_namespace = "spenso"),
        Atom::var(CS.nc).pow(2) - 1,
        Atom::var(SPENSO_TAG.bracket) * &x,
    ] {
        assert_admits_like_general(atom, true);
    }
}

#[test]
fn spenso_syntax_and_callbacks_keep_general_admission() {
    test_support::test_initialize();
    let x = parse!("x", default_namespace = "scalar_admission");
    let user = symbol!("scalar_admission::f");
    let callback = symbol!("scalar_admission::callback", norm = |_, _| {});
    let evaluated = SymbolBuilder::new(wrap_symbol!("scalar_admission::evaluated"))
        .with_evaluation_info(EvaluationInfo::new().register(|arguments: &[f64]| arguments[0]))
        .build()
        .unwrap();
    let tensor = tensor_symbol!("scalar_admission::t");
    let vector = vector_symbol!("scalar_admission::p");
    let mink = Minkowski {}.new_rep(4);
    let slot = |index| {
        mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    };
    let p = |index| function!(vector, slot(index));
    let compact = function!(vector, mink.to_symbolic([]));

    // Scalars outside the shortcut's claim still admit as scalars.
    for atom in [
        function!(symbol!("spenso::scalar_admission_f"), &x),
        Atom::var(symbol!("scalar_admission::mu", tag = SPENSO_TAG.index)),
        Atom::from(AbstractIndex::Dummy(3)),
        function!(callback, &x),
        function!(evaluated, &x),
        function!(symbol!("scalar_admission::s"; Scalar), p(1)),
        tags::scalar_store_alias(0),
        p(1) * p(1) * &x,
        dot!(&compact, &compact),
    ] {
        assert!(
            general(&atom).is_ok_and(|(value, _)| value.is_scalar()),
            "{atom}"
        );
        assert_admits_like_general(atom, false);
    }

    // Tensors are admitted by the general pipeline alone.
    for atom in [
        p(1),
        p(1) * p(2) * &x,
        ETS.metric(slot(1), slot(2)),
        function!(*SYM, p(1), p(2)),
        bracket!(&compact, &compact),
    ] {
        assert!(
            general(&atom).is_ok_and(|(value, _)| !value.is_scalar()),
            "{atom}"
        );
        assert_admits_like_general(atom, false);
    }

    // Invalid syntax, including bare tensor names, is rejected by it alone.
    for atom in [
        Atom::var(tensor),
        &x * Atom::var(tensor),
        function!(user, Atom::var(tensor)),
        Atom::var(SPENSO_TAG.chain_in),
        &x + Atom::var(SPENSO_TAG.chain_out),
        function!(user, slot(1)),
        p(1) * p(1) * p(1),
        FunctionBuilder::new(SPENSO_TAG.bracket).finish(),
    ] {
        assert!(general(&atom).is_err(), "{atom}");
        assert_admits_like_general(atom, false);
    }
}

/// The shortcut trusts that admission only reacts to Spenso tags, Spenso-owned
/// heads, chain placeholders and callbacks. Every registered symbol, as a
/// variable and as a head, must therefore be declined or admitted exactly like
/// general admission.
#[test]
fn registered_symbols_admit_like_general_when_plain() {
    test_support::test_initialize();
    let x = parse!("x", default_namespace = "scalar_admission");
    let y = parse!("y", default_namespace = "scalar_admission");
    let user = symbol!("scalar_admission::f");
    let symbols: Vec<_> = State::symbol_iter().map(|(symbol, _)| symbol).collect();
    for symbol in symbols {
        let variable = Atom::var(symbol);
        let mut forms = vec![
            variable.clone(),
            &variable * &x,
            x.pow(&variable),
            function!(user, &variable),
        ];
        // Building a call runs the head's normalization callback, which may
        // belong to a concurrently running test; the shortcut declines such
        // heads anyway.
        if symbol.get_normalization_function().is_none() {
            let call = function!(symbol, &x);
            forms.extend([
                FunctionBuilder::new(symbol).finish(),
                function!(symbol, &x, &y),
                function!(symbol, 4, 1),
                call.pow(2) * &y,
                function!(user, &call),
                call,
            ]);
        }
        for atom in forms {
            if InterfaceInference::is_plain_scalar(atom.as_view()) {
                assert_admits_like_general(atom, true);
            }
        }
    }
}

#[test]
fn plain_scalar_skips_interface_discovery() {
    test_support::test_initialize();
    let atom = parse!(
        "-𝑖*GC_11*ee^2*cw/sw",
        default_namespace = "scalar_admission"
    );
    INFERENCE_CALLS.with(|count| count.set(0));
    let fast = SymbolicTensor::<PartialStructure>::infer(atom.clone()).unwrap();
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    let (full, _) = general(&atom).unwrap();
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 1);
    assert_eq!(fast.expression, full.expression);
    assert_eq!(fast.structure, full.structure);
}

#[test]
fn eligible_zero_payload_does_not_erase_retained_tensor_metadata() {
    test_support::test_initialize();
    let mink = Minkowski {}.new_rep(4);
    let slot = mink
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(93271))
        .to_atom();
    let vector = function!(vector_symbol!("scalar_admission::zero_vector"), &slot);
    // Distinct ordinary heads retain a shared explicit pair. A squared vector
    // instead creates a closed power scope whose internal pair is local.
    let closed = SymbolicTensor::<PartialStructure>::infer(
        function!(tensor_symbol!("scalar_admission::history_left"), &slot)
            * function!(tensor_symbol!("scalar_admission::history_right"), &slot),
    )
    .unwrap();
    // Certified rewriting retains construction history. Addition/subtraction
    // deliberately bind internal pairs in their separate summand scopes.
    let closed_zero = closed
        .with_algebra_result(Atom::Zero)
        .unwrap()
        .finish_parts()
        .unwrap();
    assert!(closed_zero.is_scalar());
    assert!(closed_zero.expression.is_zero());
    assert!(!closed_zero.explicit_index_occurrences().unwrap().is_empty());
    assert_eq!(
        closed_zero
            .explicit_index_occurrences()
            .unwrap()
            .oriented_entries()
            .collect::<Vec<_>>(),
        closed
            .explicit_index_occurrences()
            .unwrap()
            .oriented_entries()
            .collect::<Vec<_>>()
    );
    let tensor = SymbolicTensor::<PartialStructure>::infer(vector).unwrap();
    let zero = tensor.try_sub(&tensor).unwrap();
    assert!(InterfaceInference::is_plain_scalar(
        zero.expression.as_view()
    ));
    assert!(zero.expression.is_zero());
    assert!(!zero.is_scalar());
    assert_eq!(zero.structure, tensor.structure);
    assert_eq!(
        zero.explicit_indices().unwrap().collect::<Vec<_>>(),
        tensor.explicit_indices().unwrap().collect::<Vec<_>>()
    );
    // The syntax shortcut applies only to fresh Atom admission. Publishing
    // established algebra must retain nominal shape and cancellation history.
    let finished = zero.finish_parts().unwrap();
    assert!(!finished.is_scalar());
    assert!(!finished.explicit_index_occurrences().unwrap().is_empty());
}
