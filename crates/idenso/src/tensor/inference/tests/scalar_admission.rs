use super::*;
use crate::tensor::AlgebraSettings;
use spenso::{tensor_symbol, vector_symbol};
use std::hint::black_box;
use std::time::Instant;
use symbolica::atom::{EvaluationInfo, SymbolBuilder};
use symbolica::{function, parse, symbol, wrap_symbol};

fn general(atom: &Atom) -> InferenceResult<(SymbolicTensor<PartialStructure>, bool)> {
    SymbolicTensor::<PartialStructure>::admit_with_contractions(atom.clone())
}

#[track_caller]
fn compare(atom: Atom) {
    let (candidate, contracted) =
        SymbolicTensor::<PartialStructure>::infer_with_contractions(atom.clone()).unwrap();
    let (baseline, baseline_contracted) = general(&atom).unwrap();
    assert_eq!(candidate.expression, baseline.expression);
    assert_eq!(candidate.structure, baseline.structure);
    assert_eq!(
        (candidate.is_metric, candidate.is_composite, contracted),
        (
            baseline.is_metric,
            baseline.is_composite,
            baseline_contracted
        )
    );
    assert_eq!(
        candidate
            .explicit_index_occurrences()
            .unwrap()
            .oriented_entries()
            .collect::<Vec<_>>(),
        baseline
            .explicit_index_occurrences()
            .unwrap()
            .oriented_entries()
            .collect::<Vec<_>>()
    );
    // Admission may populate different optimization cells. Compare their lazy outcomes.
    assert_eq!(
        candidate.normalization_is_intrinsic(),
        baseline.normalization_is_intrinsic()
    );
    assert_eq!(
        candidate.borrowed_parts().product_is_normalized(),
        baseline.borrowed_parts().product_is_normalized()
    );
    assert_eq!(
        candidate.borrowed_parts().algebra_preserves_interface(),
        baseline.borrowed_parts().algebra_preserves_interface()
    );
    assert_eq!(
        candidate.rewrites_preserve_interface(),
        baseline.rewrites_preserve_interface()
    );
}

#[test]
fn tensor_free_admission_is_opaque_and_proofs_remain_lazy() {
    test_support::test_initialize();
    let callback = symbol!("scalar_admission::callback", norm = |_, _| {});
    let evaluated = SymbolBuilder::new(wrap_symbol!("scalar_admission::evaluated"))
        .with_evaluation_info(EvaluationInfo::new().register(|arguments: &[f64]| arguments[0]))
        .build()
        .unwrap();
    let mut raw = Atom::new();
    raw.to_fun(callback);
    assert!(raw.as_view().needs_normalization());
    for atom in [
        raw,
        Atom::Zero,
        parse!("(x+y)^2/z", default_namespace = "scalar_admission"),
        function!(symbol!("spenso::scalar_admission_f"), 1),
        Atom::var(symbol!("scalar_admission::index", tag = SPENSO_TAG.index)),
        Atom::var(SPENSO_TAG.chain_in),
        FunctionBuilder::new(SPENSO_TAG.bracket).finish(),
        function!(callback, 1),
        function!(evaluated, 1),
    ] {
        assert!(!InterfaceInference::contains_tensor_symbol(atom.as_view()));
        let value = SymbolicTensor::infer(atom.clone()).unwrap();
        assert_eq!(value.expression, atom);
        assert!(value.is_scalar());
        assert!(value.explicit_index_occurrences().unwrap().is_empty());
        for proof in [
            &value.proofs.intrinsic,
            &value.proofs.normalized_product,
            &value.proofs.algebra,
            &value.proofs.rewrite,
        ] {
            assert!(proof.get().is_none());
        }
    }
    let value = SymbolicTensor::infer(function!(callback, 1)).unwrap();
    assert!(!value.normalization_is_intrinsic());
    assert!(!value.rewrites_preserve_interface());
    compare(function!(callback, 1));
    compare(function!(evaluated, 1));
    compare(parse!("(x+y)^2/z", default_namespace = "scalar_admission"));
}

#[test]
fn opaque_scalar_scoping_preserves_callback_contract() {
    test_support::test_initialize();
    let callback = symbol!("scalar_admission::opaque_scoping", norm = |_, _| {});
    let scope = symbol!("scalar_admission::scope");
    let slot = Minkowski {}
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(98767))
        .to_atom();
    for atom in [
        function!(callback, Atom::var(SPENSO_TAG.chain_in)),
        function!(callback, slot),
    ] {
        let value = SymbolicTensor::infer(atom).unwrap();
        for dummies_only in [false, true] {
            let wrapped = value.wrap_indices(scope, dummies_only).unwrap();
            assert_eq!(wrapped.structure, value.structure);
            assert!(!InterfaceInference::contains_tensor_symbol(
                wrapped.expression.as_view()
            ));
            assert!(wrapped.explicit_index_occurrences().unwrap().is_empty());
        }
        // Arbitrary checked results still use strict observation.
        assert!(
            value
                .with_checked_expression(function!(callback, Atom::var(SPENSO_TAG.chain_out)))
                .is_err()
        );
    }
}

#[test]
fn tensor_presence_scans_variables_heads_metadata_and_exponents() {
    test_support::test_initialize();
    let tensor = tensor_symbol!("scalar_admission::t");
    let bare = Atom::var(tensor);
    let scalar = symbol!("scalar_admission::metadata"; Scalar);
    for atom in [
        bare.clone(),
        function!(tensor, 1),
        function!(scalar, &bare),
        Atom::num(2).pow(&bare),
        function!(symbol!("scalar_admission::f"), &bare),
    ] {
        assert!(
            InterfaceInference::contains_tensor_symbol(atom.as_view()),
            "{atom}"
        );
        let candidate = SymbolicTensor::infer(atom.clone());
        let baseline = general(&atom);
        assert_eq!(candidate.is_ok(), baseline.is_ok(), "{atom}");
        if let (Err(candidate), Err(baseline)) = (candidate, baseline) {
            assert_eq!(candidate.to_string(), baseline.to_string());
        }
    }
}

#[test]
fn tensor_syntax_still_uses_general_admission() {
    test_support::test_initialize();
    let vector = vector_symbol!("scalar_admission::p");
    let slot = |index| {
        Minkowski {}
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    };
    let p = |index| function!(vector, slot(index));
    compare(p(1));
    compare(p(1) * p(1));
    for atom in [
        Atom::var(vector),
        p(1) * p(1) * p(1),
        p(1) * FunctionBuilder::new(SPENSO_TAG.bracket).finish(),
    ] {
        assert!(InterfaceInference::contains_tensor_symbol(atom.as_view()));
        assert!(general(&atom).is_err(), "{atom}");
        assert!(SymbolicTensor::infer(atom).is_err());
    }
}

#[test]
fn scalar_typed_algebra_retains_history_without_inference() {
    test_support::test_initialize();
    let expression = parse!("(x+y)^2/z", default_namespace = "scalar_admission");
    let value = SymbolicTensor::infer(expression).unwrap();
    INFERENCE_CALLS.with(|count| count.set(0));
    let product = value.multiply(&value).unwrap();
    let simplified = product
        .simplify_algebra(&AlgebraSettings::default())
        .unwrap();
    assert!(simplified.is_scalar());
    assert!(simplified.explicit_index_occurrences().unwrap().is_empty());
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    let vector = vector_symbol!("scalar_admission::replacement");
    let slot = Minkowski {}
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(98765))
        .to_atom();
    let expression = function!(vector, slot);
    assert!(value.with_checked_expression(expression.clone()).is_err());
    let (transformed, _) = value.with_transformed_expression(expression).unwrap();
    assert_eq!(transformed.rank(), 1);
    assert_eq!(
        transformed
            .explicit_index_occurrences()
            .unwrap()
            .entries()
            .count(),
        1
    );
    assert!(
        value
            .with_transformed_expression(Atom::var(vector))
            .is_err()
    );
}

/// Same executable, same general-admission implementation as the parent source;
/// excludes parsing and symbol registration. Each iteration owns and drops a new carrier.
#[test]
#[ignore = "targeted admission timing; run explicitly with --nocapture"]
fn scalar_admission_timing() {
    test_support::test_initialize();
    let vector = vector_symbol!("scalar_admission::benchmark_p");
    let tensor = function!(
        vector,
        Minkowski {}
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(98766))
            .to_atom()
    );
    for (name, atom) in [
        (
            "parameter",
            parse!("ee^2*cw/sw", default_namespace = "scalar_admission"),
        ),
        (
            "functions",
            parse!(
                "cos(x)*sin(y)+(a+b)^2/z",
                default_namespace = "scalar_admission"
            ),
        ),
        ("tensor", tensor.clone()),
        (
            "weighted_tensor",
            tensor
                * parse!(
                    "ee^2*cw*sw^(-1)*(cos(x)+sin(y))",
                    default_namespace = "scalar_admission"
                ),
        ),
    ] {
        for proofs in [false, true] {
            for candidate in [false, true] {
                let start = Instant::now();
                for _ in 0..20_000 {
                    let (value, _) = if candidate {
                        SymbolicTensor::<PartialStructure>::infer_with_contractions(black_box(
                            atom.clone(),
                        ))
                    } else {
                        general(black_box(&atom))
                    }
                    .unwrap();
                    if proofs {
                        black_box(value.normalization_is_intrinsic());
                        black_box(value.borrowed_parts().product_is_normalized());
                        black_box(value.borrowed_parts().algebra_preserves_interface());
                        black_box(value.rewrites_preserve_interface());
                    }
                    black_box(value);
                }
                println!(
                    "{name} proofs={proofs} candidate={candidate}: {:.1} ns/admission",
                    start.elapsed().as_nanos() as f64 / 20_000.0
                );
            }
        }
    }
}
