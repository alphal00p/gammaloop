use pyo3::{
    Bound, PyResult,
    types::{PyModule, PyModuleMethods},
};

mod pipeline;
mod tooling;

pub(crate) use pipeline::{PyReductionStatus, algebra_settings, contraction_representations};

pub(crate) use tooling::{
    CanonicalizationError, CookingError, DiracAdjointError, Intern, NetworkToolingError,
};

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyReductionStatus>()?;
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
        "CookingError",
        "DiracAdjointError",
        "NetworkToolingError",
        "ReductionStatus",
        "canonize",
        "dirac_adjoint",
        "list_dangling",
        "simplify_algebra",
        "to_dots",
        "undo_dots",
        "undo_chain",
        "undo_trace",
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
            for name in [
                "CookSettings",
                "CookMode",
                "CookSourceFilter",
                "CookTagFilter",
            ] {
                assert!(!module.hasattr(name).unwrap(), "retired class {name}");
            }
            let tensor_expression = py.get_type::<TensorExpression>();
            for name in [
                "alias_subtensors",
                "simplify",
                "simplify_gamma",
                "simplify_color",
                "simplify_epsilon",
                "undo_all",
                "chainify",
                "undo_schoonschip",
                "undo_single_length",
                "normalize_chains",
                "reinfer",
                "unsafe_from_expression",
                "to_color_casimir",
                "conjugate_transpose",
                "spenso_conjugate",
                "to_cof_dimension_invariants",
                "wrap_dummies",
                "cook",
                "uncook",
                "cook_function",
                "cook_indices",
            ] {
                assert!(
                    !tensor_expression.hasattr(name).unwrap(),
                    "retired method {name}"
                );
            }
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
                        "symbolica.community.tensor",
                        "wrong public module for {name}"
                    );
                }
            }
            assert!(module.getattr("initialize").is_err());
            assert!(module.getattr("initialize_module").is_err());
        });
    }

    #[test]
    fn independent_reduction_settings_and_explicit_notation_runtime() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "spenso")?;
            crate::initialize_spenso(&module)?;
            let globals = pyo3::types::PyDict::new(py);
            globals.set_item("sp", module)?;
            py.run(
                c"
spin = sp.Representation.bis(4)
lorentz = sp.Representation.mink(4)
gamma = sp.TensorExpression.dirac_gamma(4)
source = gamma(spin('i'), spin('j'), lorentz('mu')) * gamma(spin('j'), spin('i'), lorentz('nu'))
assert source.contract(representations=[]) == source
assert source.simplify_algebra(gamma=False, color=False, contract='none') == source
assert source.contract(representations=[lorentz]).to_expression() == source.contract(representations=[lorentz.name]).to_expression()
collected = source.contract()
assert source.contract(expand=False) == collected
assert 'trace' in collected.to_expression().format_plain()
unfolded = collected.undo_trace()
assert 'trace' not in unfolded.to_expression().format_plain()
assert 'chain' in unfolded.to_expression().format_plain()
explicit = unfolded.undo_chain()
assert 'chain' not in explicit.to_expression().format_plain()
assert explicit.components() == collected.components()
reduced = source.simplify_algebra(gamma=True)
assert reduced.reduction_status == sp.ReductionStatus.Complete, (reduced.reduction_status, str(reduced.to_expression()))
assert reduced.expand().components() == collected.components()
assert source.simplify_algebra(gamma=True, max_steps_per_domain=None) == reduced
assert source.simplify_algebra(gamma=True, contract='minimal') == reduced
assert source.contract(max_steps_per_domain=None) == collected
capped = source.simplify_algebra(gamma=True, max_steps_per_domain=0)
assert capped.reduction_status == sp.ReductionStatus.Capped
assert capped == source
generator = sp.TensorExpression.color_t(8, 3)
word = generator('a', 'i', 'j') * generator('a', 'k', 'l')
for mode in ('fully', 'minimal', 'dots', 'selected'):
    options = dict(color=True, color_expand_fierz=False, contract=mode)
    if mode == 'selected':
        options['representations'] = [sp.Representation.cof(3)]
    result = word.simplify_algebra(**options)
    assert result.reduction_status == sp.ReductionStatus.Complete, mode
    assert result.components() == word.components()
    assert result.simplify_algebra(**options).to_expression() == result.to_expression()
assert not hasattr(sp.Tensor, 'contract')
assert hasattr(sp.Tensor, 'contract_ports')
",
                Some(&globals),
                None,
            )
        })
        .unwrap();
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
            // Deliberately inject malformed internals to test tooling error translation.
            let expression = TensorExpression::from_parts_unchecked(
                py,
                malformed,
                PartialStructure::from_logical_slots([]),
                None,
                Vec::new(),
            )?;

            let error = expression
                .bind(py)
                .call_method0("undo_dots")
                .expect_err("malformed dot notation should return an error");
            assert!(error.is_instance_of::<NetworkToolingError>(py));
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("cannot parse tensor network"));
            assert!(error.to_string().contains("Invalid dot function"));

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

            for (class, expected) in [
                (
                    py.get_type::<TensorExpression>().into_any(),
                    "(expression, *, structure=None, intern=None)",
                ),
                (
                    py.get_type::<TensorExpression>()
                        .getattr("contract")
                        .unwrap(),
                    "($self, *, representations=None, metrics=True, rank_one=True, collect_chains=True, collect_traces=True, expand=True, order=None, max_steps_per_domain=None)",
                ),
                (
                    py.get_type::<TensorExpression>()
                        .getattr("simplify_algebra")
                        .unwrap(),
                    "($self, *, gamma=True, color=True, epsilon=False, contract=\"fully\", representations=None, collect_coefficients=True, gamma_output=None, gamma_ordering=None, gamma_evaluate_traces=None, gamma0=None, gamma_conjugate=None, gamma_expand_three_gamma_epsilon=None, color_evaluate_traces=None, color_expand_fierz=None, color_substitute_cof_dimension_invariants=None, max_steps_per_domain=None)",
                ),
            ] {
                let name = class.getattr("__name__").unwrap();
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
                "dirac_adjoint",
                "list_dangling",
                "simplify_algebra",
                "to_dots",
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
