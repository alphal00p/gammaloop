//! Python keyword conversion for the shared tensor-algebra scheduler.
use crate::{metadata::SpensoRepresentationName, structure::SpensoRepresentation};
use idenso::{
    color::ColorSimplifySettings,
    dirac::{GammaChainOrdering, GammaOutput, GammaSimplifySettings},
    tensor::simplification::{AlgebraContraction, AlgebraSettings},
};
use pyo3::{exceptions::PyValueError, prelude::*};
use spenso::structure::representation::LibraryRep;

pub(crate) fn contraction_representations(
    values: Option<Vec<Bound<'_, PyAny>>>,
) -> PyResult<Option<Vec<LibraryRep>>> {
    values
        .map(|values| {
            values
                .into_iter()
                .map(|value| {
                    if let Ok(value) = value.extract::<PyRef<'_, SpensoRepresentation>>() {
                        Ok(value.representation.rep)
                    } else {
                        Ok(value.extract::<PyRef<'_, SpensoRepresentationName>>()?.rep)
                    }
                })
                .collect()
        })
        .transpose()
}

// Optional family keywords distinguish omission from an explicitly supplied default.
// Only this Python boundary translates names; native settings own all semantics.
#[allow(clippy::too_many_arguments)]
pub(crate) fn algebra_settings(
    gamma: bool,
    color: bool,
    epsilon: bool,
    contract: &str,
    representations: Option<Vec<Bound<'_, PyAny>>>,
    collect_coefficients: bool,
    gamma_output: Option<&str>,
    gamma_ordering: Option<&str>,
    gamma_evaluate_traces: Option<bool>,
    gamma0: Option<bool>,
    gamma_conjugate: Option<bool>,
    gamma_expand_three_gamma_epsilon: Option<bool>,
    color_evaluate_traces: Option<bool>,
    color_expand_fierz: Option<bool>,
    color_substitute_cof_dimension_invariants: Option<bool>,
    max_steps_per_domain: Option<usize>,
) -> PyResult<AlgebraSettings> {
    let representations = contraction_representations(representations)?;
    let contract = match (contract, representations) {
        ("fully", None) => AlgebraContraction::Fully,
        ("minimal", None) => AlgebraContraction::Minimal,
        ("dots", None) => AlgebraContraction::Dots,
        ("none", None) => AlgebraContraction::None,
        ("selected", Some(representations)) => AlgebraContraction::Selected(representations),
        ("selected", None) => {
            return Err(PyValueError::new_err(
                "contract='selected' requires representations (an empty list permits no additional contractions)",
            ));
        }
        ("fully" | "minimal" | "dots" | "none", Some(_)) => {
            return Err(PyValueError::new_err(
                "representations requires contract='selected'",
            ));
        }
        _ => {
            return Err(PyValueError::new_err(
                "contract must be 'fully', 'minimal', 'dots', 'selected', or 'none'",
            ));
        }
    };
    for (enabled, options) in [
        (
            gamma,
            vec![
                ("gamma_output", gamma_output.is_some()),
                ("gamma_ordering", gamma_ordering.is_some()),
                ("gamma_evaluate_traces", gamma_evaluate_traces.is_some()),
                ("gamma0", gamma0.is_some()),
                ("gamma_conjugate", gamma_conjugate.is_some()),
                (
                    "gamma_expand_three_gamma_epsilon",
                    gamma_expand_three_gamma_epsilon.is_some(),
                ),
            ],
        ),
        (
            color,
            vec![
                ("color_evaluate_traces", color_evaluate_traces.is_some()),
                ("color_expand_fierz", color_expand_fierz.is_some()),
                (
                    "color_substitute_cof_dimension_invariants",
                    color_substitute_cof_dimension_invariants.is_some(),
                ),
            ],
        ),
    ] {
        if !enabled {
            for (name, supplied) in options {
                if supplied {
                    let family = if name.starts_with("color_") {
                        "color"
                    } else {
                        "gamma"
                    };
                    return Err(PyValueError::new_err(format!(
                        "{name} requires {family}=True"
                    )));
                }
            }
        }
    }
    let mut gamma_settings = GammaSimplifySettings::default();
    if let Some(output) = gamma_output {
        gamma_settings.output = match output {
            "reduced" => GammaOutput::Reduced,
            "chains" => GammaOutput::Chains,
            _ => {
                return Err(PyValueError::new_err(
                    "gamma_output must be 'reduced' or 'chains'",
                ));
            }
        };
    }
    if let Some(ordering) = gamma_ordering {
        gamma_settings.chain_ordering = match ordering {
            "repeated_pairs" => GammaChainOrdering::RepeatedPairs,
            "canonical" => GammaChainOrdering::Canonical,
            _ => {
                return Err(PyValueError::new_err(
                    "gamma_ordering must be 'repeated_pairs' or 'canonical'",
                ));
            }
        };
    }
    if let Some(value) = gamma_evaluate_traces {
        gamma_settings.evaluate_traces = value;
    }
    if let Some(value) = gamma0 {
        gamma_settings.gamma0 = value;
    }
    if let Some(value) = gamma_conjugate {
        gamma_settings.conjugate = value;
    }
    if let Some(value) = gamma_expand_three_gamma_epsilon {
        gamma_settings.expand_three_gamma_epsilon = value;
    }
    let mut color_settings = ColorSimplifySettings::default();
    if let Some(value) = color_evaluate_traces {
        color_settings.evaluate_traces = value;
    }
    if let Some(value) = color_expand_fierz {
        color_settings.expand_cross_chain_fierz = value;
    }
    if let Some(value) = color_substitute_cof_dimension_invariants {
        color_settings.substitute_cof_dimension_invariants = value;
    }
    Ok(AlgebraSettings {
        gamma: gamma.then_some(gamma_settings),
        color: color.then_some(color_settings),
        epsilon,
        contract,
        collect_coefficients,
        max_passes: max_steps_per_domain,
    })
}

/// Completion of a tensor reduction under its selected settings.
///
/// Complete means no eligible work remains. Deferred means exact eligible work
/// remains for another operation. Capped means an explicitly requested work limit
/// stopped processing. A partial result still represents the exact original tensor;
/// completion is relative to the requested operations, not every possible identity.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import TensorExpression, ReductionStatus
/// >>> result = TensorExpression(3).contract()
/// >>> result.reduction_status == ReductionStatus.Complete
/// True
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyclass_enum
)]
#[pyclass(
    frozen,
    from_py_object,
    eq,
    eq_int,
    name = "ReductionStatus",
    module = "symbolica.community.tensor"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyReductionStatus {
    /// No eligible work remains under the requested settings.
    Complete,
    /// Exact eligible work remains for a subsequent operation.
    Deferred,
    /// An explicit work limit stopped processing; the remaining expression is exact.
    Capped,
}

impl From<idenso::tensor::simplification::ReductionStatus> for PyReductionStatus {
    fn from(status: idenso::tensor::simplification::ReductionStatus) -> Self {
        use idenso::tensor::simplification::ReductionStatus;
        match status {
            ReductionStatus::Complete => Self::Complete,
            ReductionStatus::Deferred => Self::Deferred,
            ReductionStatus::Capped => Self::Capped,
        }
    }
}
