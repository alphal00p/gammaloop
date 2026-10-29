use super::*;
use crate::color::CS;
use crate::tensor::AlgebraSettings;
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
        (Ok((fast, fast_contracted)), Ok((general, general_contracted))) => {
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
            for (fast, general) in [
                (&fast.proofs.validated, &general.proofs.validated),
                (&fast.proofs.intrinsic, &general.proofs.intrinsic),
                (
                    &fast.proofs.normalized_product,
                    &general.proofs.normalized_product,
                ),
                (&fast.proofs.algebra, &general.proofs.algebra),
                (&fast.proofs.rewrite, &general.proofs.rewrite),
            ] {
                assert_eq!(fast.get(), general.get(), "{atom}");
            }
            assert_eq!(fast.proofs.contracted, general.proofs.contracted, "{atom}");
            assert_eq!(fast.proofs.frontier, general.proofs.frontier, "{atom}");
            assert_eq!(fast.proofs.reduction, general.proofs.reduction, "{atom}");
            let fast_occurrences = fast.proofs.occurrences.get().unwrap();
            let general_occurrences = general.proofs.occurrences.get().unwrap();
            fast_occurrences.validate().unwrap();
            general_occurrences.validate().unwrap();
            assert_eq!(
                fast_occurrences.oriented_entries().collect::<Vec<_>>(),
                general_occurrences.oriented_entries().collect::<Vec<_>>(),
                "{atom}"
            );
            assert_eq!(
                fast_occurrences.selected_pair, general_occurrences.selected_pair,
                "{atom}"
            );
            for value in [&fast, &general] {
                assert!(value.proofs.shallow_graph.get().is_none(), "{atom}");
                assert!(value.proofs.leaf_interfaces.get().is_none(), "{atom}");
                assert!(value.proofs.observations.get().is_none(), "{atom}");
                assert!(value.same_reduction_state(&value.clone()), "{atom}");
            }
            if plain {
                // The shortcut only claims inputs that admission fully certifies.
                assert!(general.is_scalar(), "{atom}");
                for proof in [
                    &general.proofs.validated,
                    &general.proofs.intrinsic,
                    &general.proofs.normalized_product,
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
fn admitted_scalars_retain_metadata_through_typed_algebra() {
    test_support::test_initialize();
    let expression = parse!(
        "(x+y)^2/z",
        default_namespace = "scalar_admission_operations"
    );
    let fast = SymbolicTensor::infer(expression.clone()).unwrap();
    let (general, _) = general(&expression).unwrap();
    let fast_product = fast.multiply(&fast).unwrap();
    let general_product = general.multiply(&general).unwrap();
    assert!(fast_product.same_value(&general_product));
    let exponent = SymbolicTensor::infer(Atom::num(2)).unwrap();
    let fast_power = fast.pow(&exponent).unwrap();
    let general_power = general.pow(&exponent).unwrap();
    assert!(fast_power.same_value(&general_power));
    for value in [fast_product, general_product, fast_power, general_power] {
        assert!(value.explicit_index_occurrences().unwrap().is_empty());
        let simplified = value.simplify_algebra(&AlgebraSettings::default()).unwrap();
        assert!(simplified.same_value(&value));
        assert!(simplified.explicit_index_occurrences().unwrap().is_empty());
    }
}
