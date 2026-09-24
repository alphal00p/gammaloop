use pyo3::{
    Bound, PyResult,
    types::{PyModule, PyModuleMethods},
};

mod algebra;
pub(crate) mod expansion;
mod pipeline;
mod tooling;

pub(crate) use pipeline::PySimplifySettings;

pub(crate) use algebra::{
    GammaConjugationError, PyColorCasimirSettings, PyColorSimplifySettings, PyGammaSimplifySettings,
};
pub(crate) use tooling::{
    CanonicalizationError, CookingError, DiracAdjointError, DotExpansionError, NetworkToolingError,
    PyCookSettings, PySchoonschipSettings,
};

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySimplifySettings>()?;
    algebra::register(module)?;
    tooling::register(module)
}

#[cfg(test)]
mod tests {
    use crate::expression::TensorExpression;
    use pyo3::{
        Python,
        exceptions::PyValueError,
        types::{PyAnyMethods, PyDictMethods, PyListMethods, PyModuleMethods},
    };
    use spenso::network::tags::SPENSO_TAG;
    use spenso::structure::partial::{PartialStructure, PartialStructureExt};
    use symbolica::{
        atom::{Atom, FunctionBuilder},
        symbol,
    };

    use super::*;

    const PUBLIC_API: &[&str] = &[
        "CanonicalizationError",
        "ColorCasimirSettings",
        "ColorSimplifySettings",
        "CookMode",
        "CookSettings",
        "CookSourceFilter",
        "CookTagFilter",
        "CookingError",
        "DiracAdjointError",
        "DotExpansionError",
        "GammaChainOrdering",
        "GammaConjugationError",
        "GammaSimplifySettings",
        "NetworkToolingError",
        "SchoonschipContractionOrder",
        "SchoonschipMode",
        "SchoonschipSettings",
        "SchoonschipTraversal",
        "SimplifySettings",
        "alias_subtensors",
        "canonize",
        "chainify",
        "collect_chains",
        "collect_color",
        "collect_color_constants",
        "collect_gamma_chains",
        "conjugate_transpose",
        "cook",
        "cook_function",
        "cook_indices",
        "dirac_adjoint",
        "expand_bis",
        "expand_color",
        "expand_dots",
        "expand_in_patterns",
        "expand_metrics",
        "expand_mink",
        "expand_mink_bis",
        "list_dangling",
        "metric_shorthand_to_dot",
        "normalize_chains",
        "normalize_dots",
        "schoonschip",
        "schoonschip_net",
        "simplify_color",
        "simplify",
        "simplify_epsilon",
        "simplify_gamma",
        "simplify_gamma0",
        "simplify_gamma_conjugate",
        "simplify_metrics",
        "spenso_conjugate",
        "to_cof_dimension_invariants",
        "to_color_casimir",
        "to_dots",
        "uncook",
        "undo_all",
        "undo_chain",
        "undo_dots",
        "undo_schoonschip",
        "undo_single_length",
        "undo_trace",
        "wrap_color",
        "wrap_dummies",
        "wrap_indices",
    ];

    #[test]
    fn registers_exact_public_python_surface() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "spenso").unwrap();
            register(&module).unwrap();

            let mut actual = module
                .dict()
                .keys()
                .iter()
                .filter_map(|key| key.extract::<String>().ok())
                .filter(|name| !name.starts_with('_'))
                .collect::<Vec<_>>();
            let mut expected = PUBLIC_API
                .iter()
                .filter(|name| name.starts_with(char::is_uppercase))
                .map(|name| (*name).to_string())
                .collect::<Vec<_>>();
            actual.sort();
            expected.sort();
            assert_eq!(actual, expected);
            let tensor_expression = py.get_type::<TensorExpression>();
            for name in PUBLIC_API {
                if name.starts_with(char::is_lowercase) {
                    assert!(
                        tensor_expression.hasattr(*name).unwrap(),
                        "missing TensorExpression.{name}"
                    );
                    assert!(
                        !module.hasattr(*name).unwrap(),
                        "unexpected free function {name}"
                    );
                } else {
                    assert_eq!(
                        module
                            .getattr(*name)
                            .unwrap()
                            .getattr("__module__")
                            .unwrap()
                            .extract::<String>()
                            .unwrap(),
                        "symbolica.community.spenso",
                        "wrong public module for {name}"
                    );
                }
            }
            assert!(module.getattr("initialize").is_err());
            assert!(module.getattr("initialize_module").is_err());
        });
    }

    #[test]
    fn network_tooling_failures_are_python_value_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "spenso")?;
            register(&module)?;
            let malformed = FunctionBuilder::new(SPENSO_TAG.dot)
                .add_arg(Atom::var(symbol!("malformed_dot_operand")))
                .finish();
            let expression = TensorExpression::from_atom_interface(
                py,
                malformed,
                PartialStructure::from_logical_slots([]),
            )?;

            for name in ["undo_dots", "schoonschip_net"] {
                let error = expression
                    .bind(py)
                    .call_method0(name)
                    .expect_err("malformed dot notation should return an error");
                assert!(error.is_instance_of::<NetworkToolingError>(py));
                assert!(error.is_instance_of::<PyValueError>(py));
                assert!(error.to_string().contains("cannot parse tensor network"));
                assert!(error.to_string().contains("Invalid dot function"));
            }

            let error = expression
                .bind(py)
                .call_method0("canonize")
                .expect_err("malformed dot notation should not canonicalize");
            assert!(error.is_instance_of::<CanonicalizationError>(py));
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("cannot parse tensor network"));
            assert!(error.to_string().contains("Invalid dot function"));

            let error = expression
                .bind(py)
                .call_method0("dirac_adjoint")
                .expect_err("malformed dot notation should not have a Dirac adjoint");
            assert!(error.is_instance_of::<DiracAdjointError>(py));
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("cannot parse tensor network"));
            assert!(error.to_string().contains("Invalid dot function"));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn public_python_surface_has_runtime_docstrings() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "spenso").unwrap();
            register(&module).unwrap();

            let tensor_expression = py.get_type::<TensorExpression>();
            for name in PUBLIC_API {
                let item = if name.starts_with(char::is_uppercase) {
                    module.getattr(*name).unwrap()
                } else {
                    tensor_expression.getattr(*name).unwrap()
                };
                let documentation = item
                    .getattr("__doc__")
                    .unwrap()
                    .extract::<Option<String>>()
                    .unwrap();
                assert!(
                    documentation.is_some_and(|documentation| !documentation.trim().is_empty()),
                    "{name} is missing its Python docstring"
                );
            }
        });
    }

    #[test]
    fn python_signatures_use_public_names_and_concrete_defaults() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "spenso").unwrap();
            register(&module).unwrap();
            let inspect_signature = PyModule::import(py, "inspect")
                .unwrap()
                .getattr("signature")
                .unwrap();

            for (name, expected) in [
                (
                    "GammaSimplifySettings",
                    "(*, chain_ordering=None, evaluate_traces=True, expand_three_gamma_epsilon=False)",
                ),
                (
                    "CookSettings",
                    "(*, mode=None, source=None, output_tags=None, preserve_tags=False)",
                ),
                (
                    "SchoonschipSettings",
                    "(*, depth_limit=1, mode=None, traversal=None, expand_contracted_sums=False, simplify_chain_like_functions=False, schoonschip_rank1_tensors=True, contraction_order=None)",
                ),
            ] {
                let class = module.getattr(name).unwrap();
                let signature = class
                    .getattr("__text_signature__")
                    .unwrap()
                    .extract::<String>()
                    .unwrap();
                assert_eq!(signature, expected, "unexpected signature for {name}");
                let inspected = inspect_signature
                    .call1((&class,))
                    .unwrap()
                    .str()
                    .unwrap()
                    .to_string();
                assert!(
                    !inspected.contains("..."),
                    "inspect.signature retained an ellipsis for {name}: {inspected}"
                );
            }

            for name in [
                "cook_function",
                "cook_indices",
                "dirac_adjoint",
                "expand_bis",
                "expand_color",
                "expand_metrics",
                "expand_mink",
                "expand_mink_bis",
                "list_dangling",
                "simplify_color",
                "simplify_gamma",
                "simplify_metrics",
                "to_dots",
                "wrap_dummies",
                "wrap_indices",
            ] {
                let signature = py
                    .get_type::<TensorExpression>()
                    .getattr(name)
                    .unwrap()
                    .getattr("__text_signature__")
                    .unwrap()
                    .extract::<String>()
                    .unwrap();
                assert!(
                    !signature.contains("expression") && signature.contains("self"),
                    "unexpected signature for {name}: {signature}"
                );
            }
        });
    }
}
