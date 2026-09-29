//! Python conversion for the shared tensor-algebra scheduler.
use super::{PyColorSimplifySettings, PyGammaSimplifySettings};
use crate::{metadata::SpensoRepresentationName, structure::SpensoRepresentation};
use idenso::tensor::{ContractSettings, simplification::AlgebraSettings};
use pyo3::prelude::*;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use spenso::structure::representation::{LibraryRep, RepName};

/// Select optional tensor-algebra identities for reduce_algebra().
///
/// Instances are immutable. Dirac, color, and Levi-Civita identities are disabled
/// by default. Enabling a family includes the contractions needed by that family;
/// use ContractSettings for general structural contraction.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import AlgebraSettings, GammaSimplifySettings
/// >>> settings = AlgebraSettings(gamma=GammaSimplifySettings())
/// >>> settings.color is None
/// True
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "AlgebraSettings",
    module = "symbolica.community.tensor"
)]
#[derive(Clone, Copy, Default)]
pub(crate) struct PyAlgebraSettings {
    inner: AlgebraSettings,
}

impl PyAlgebraSettings {
    pub(crate) fn rust(&self) -> &AlgebraSettings {
        &self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyAlgebraSettings {
    /// Choose identity families and a transformation budget.
    ///
    /// Parameters
    /// ----------
    /// gamma : GammaSimplifySettings, optional
    ///     Dirac-algebra options. None disables this family.
    /// color : ColorSimplifySettings, optional
    ///     Color-algebra options. None disables this family.
    /// epsilon : bool, default False
    ///     Enable supported Levi-Civita identities and their required contractions.
    /// max_passes : int, default 16
    ///     Maximum successful transformations in each root or shared-definition domain.
    ///     Zero is allowed. Inspect the result's reduction_status to detect unfinished work.
    ///
    /// Returns
    /// -------
    /// AlgebraSettings
    ///     Immutable options for reduce_algebra().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import AlgebraSettings, GammaSimplifySettings
    /// >>> settings = AlgebraSettings(gamma=GammaSimplifySettings.canonical())
    /// >>> settings.max_passes
    /// 16
    #[new]
    #[pyo3(signature=(*, gamma=None, color=None, epsilon=false, max_passes=16))]
    pub(crate) fn new(
        gamma: Option<PyGammaSimplifySettings>,
        color: Option<PyColorSimplifySettings>,
        epsilon: bool,
        max_passes: usize,
    ) -> Self {
        let inner = AlgebraSettings {
            gamma: gamma.map(|settings| settings.rust()),
            color: color.map(|settings| settings.rust()),
            epsilon,
            max_passes,
        };
        Self { inner }
    }

    /// Enable Dirac, color, and Levi-Civita identities.
    ///
    /// Returns
    /// -------
    /// AlgebraSettings
    ///     Default options for all three families, including their necessary contractions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import AlgebraSettings
    /// >>> AlgebraSettings.hep().epsilon
    /// True
    #[staticmethod]
    fn hep() -> Self {
        Self {
            inner: AlgebraSettings::hep(),
        }
    }

    /// Dirac-algebra options; None disables this family.
    ///
    /// Returns
    /// -------
    /// GammaSimplifySettings or None
    ///     Dirac-algebra options; None disables this family.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import AlgebraSettings
    /// >>> value = AlgebraSettings().gamma
    /// >>> value is None
    /// True
    #[getter]
    fn gamma(&self) -> Option<PyGammaSimplifySettings> {
        self.inner.gamma.map(PyGammaSimplifySettings::from_rust)
    }
    /// Color-algebra options; None disables this family.
    ///
    /// Returns
    /// -------
    /// ColorSimplifySettings or None
    ///     Color-algebra options; None disables this family.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import AlgebraSettings
    /// >>> value = AlgebraSettings().color
    /// >>> value is None
    /// True
    #[getter]
    fn color(&self) -> Option<PyColorSimplifySettings> {
        self.inner.color.map(PyColorSimplifySettings::from_rust)
    }
    /// Whether supported Levi-Civita identities are enabled.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Whether supported Levi-Civita identities are enabled.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import AlgebraSettings
    /// >>> AlgebraSettings().epsilon
    /// False
    #[getter]
    fn epsilon(&self) -> bool {
        self.inner.epsilon
    }
    /// Maximum successful transformations in each root or shared-definition domain. Zero is allowed.
    ///
    /// Returns
    /// -------
    /// int
    ///     Maximum successful transformations in each root or shared-definition domain. Zero is allowed.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import AlgebraSettings
    /// >>> AlgebraSettings().max_passes
    /// 16
    #[getter]
    fn max_passes(&self) -> usize {
        self.inner.max_passes
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.AlgebraSettings()
    /// >>> text = repr(settings)
    fn __repr__(self_: PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = self_.into_pyobject(py)?;
        crate::display::constructor_repr(
            object.as_any(),
            &[
                ("gamma", "gamma"),
                ("color", "color"),
                ("epsilon", "epsilon"),
                ("max_passes", "max_passes"),
            ],
        )
    }
}

/// Select structural contractions and compact tensor forms.
///
/// The default enables metrics, rank-one contractions, matrix-chain collection,
/// and trace collection in every representation. Filters select representation
/// names, independently of dimensions or other axes of a tensor. An empty filter
/// permits no representations. Intrinsic constructor normalization is unaffected.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import ContractSettings, Representation
/// >>> settings = ContractSettings(representations=[Representation.euc(3)], collect_traces=False)
/// >>> settings.collect_traces
/// False
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    from_py_object,
    name = "ContractSettings",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub(crate) struct PyContractSettings {
    representations: Option<Vec<LibraryRep>>,
    metrics: bool,
    rank_one: bool,
    collect_chains: bool,
    collect_traces: bool,
    order: Option<Vec<usize>>,
    max_passes: usize,
}

impl Default for PyContractSettings {
    fn default() -> Self {
        let settings = ContractSettings::default();
        Self {
            representations: None,
            metrics: settings.metrics,
            rank_one: settings.rank_one,
            collect_chains: settings.collect_chains,
            collect_traces: settings.collect_traces,
            order: None,
            max_passes: settings.max_passes,
        }
    }
}

impl PyContractSettings {
    pub(crate) fn rust(&self) -> ContractSettings<'_> {
        ContractSettings {
            representations: self.representations.as_deref(),
            metrics: self.metrics,
            rank_one: self.rank_one,
            collect_chains: self.collect_chains,
            collect_traces: self.collect_traces,
            max_passes: self.max_passes,
            order: self.order.as_deref(),
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyContractSettings {
    /// Choose which compatible indices and matrix channels to contract.
    ///
    /// Parameters
    /// ----------
    /// representations : sequence of Representation or RepresentationName, optional
    ///     Restrict eligible connecting axes by representation name. None permits all
    ///     representations; an empty sequence permits none. Dimensions are not filtered.
    /// metrics : bool, default True
    ///     Contract compatible metric tensors.
    /// rank_one : bool, default True
    ///     Contract compatible rank-one tensors into compact scalar products.
    /// collect_chains : bool, default True
    ///     Collect connected matrix factors into compact chains.
    /// collect_traces : bool, default True
    ///     Collect closed matrix products into compact traces, without evaluating
    ///     Dirac or color trace identities.
    /// order : sequence of int, optional
    ///     Requested order of tensor factors for contraction. None selects the default.
    ///     This does not reorder the tensor's external axes.
    /// max_passes : int, default 16
    ///     Maximum successful transformations in each root or shared-definition domain.
    ///     Zero is allowed. Inspect reduction_status for unfinished work.
    ///
    /// Returns
    /// -------
    /// ContractSettings
    ///     Immutable settings for contract().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings, Representation, TensorExpression
    /// >>> space = Representation.euc(3)
    /// >>> settings = ContractSettings(collect_chains=False, collect_traces=False)
    /// >>> result = TensorExpression.g(space)("i", "i").contract(settings)
    /// >>> result.to_expression().is_scalar
    /// True
    #[new]
    #[pyo3(signature=(*, representations=None, metrics=true, rank_one=true, collect_chains=true, collect_traces=true, order=None, max_passes=16))]
    fn new(
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Sequence[Representation | RepresentationName]]", imports=("typing")))]
        representations: Option<Vec<Bound<'_, PyAny>>>,
        metrics: bool,
        rank_one: bool,
        collect_chains: bool,
        collect_traces: bool,
        order: Option<Vec<usize>>,
        max_passes: usize,
    ) -> PyResult<Self> {
        let representations = representations
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| -> PyResult<LibraryRep> {
                        if let Ok(value) = value.extract::<PyRef<'_, SpensoRepresentation>>() {
                            Ok(value.representation.rep)
                        } else {
                            value
                                .extract::<PyRef<'_, SpensoRepresentationName>>()
                                .map(|value| value.rep)
                                .map_err(Into::into)
                        }
                    })
                    .collect::<PyResult<Vec<_>>>()
            })
            .transpose()?;
        Ok(Self {
            representations,
            metrics,
            rank_one,
            collect_chains,
            collect_traces,
            order,
            max_passes,
        })
    }

    /// Allowed representation names; None permits all and an empty list permits none.
    ///
    /// Returns
    /// -------
    /// list of RepresentationName or None
    ///     Allowed representation names; None permits all and an empty list permits none.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings
    /// >>> ContractSettings().representations is None
    /// True
    #[getter]
    fn representations(&self) -> Option<Vec<SpensoRepresentationName>> {
        self.representations.as_ref().map(|values| {
            values
                .iter()
                .map(|&rep| SpensoRepresentationName { rep })
                .collect()
        })
    }
    /// Whether compatible metric tensors are contracted.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Whether compatible metric tensors are contracted.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings
    /// >>> ContractSettings().metrics
    /// True
    #[getter]
    fn metrics(&self) -> bool {
        self.metrics
    }
    /// Whether compatible rank-one factors form compact scalar products.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Whether compatible rank-one factors form compact scalar products.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings
    /// >>> ContractSettings().rank_one
    /// True
    #[getter]
    fn rank_one(&self) -> bool {
        self.rank_one
    }
    /// Whether compatible matrix products are collected as chains.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Whether compatible matrix products are collected as chains.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings
    /// >>> ContractSettings().collect_chains
    /// True
    #[getter]
    fn collect_chains(&self) -> bool {
        self.collect_chains
    }
    /// Whether closed matrix products are collected as traces. This does not evaluate their algebraic identities.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Whether closed matrix products are collected as traces. This does not evaluate their algebraic identities.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings
    /// >>> ContractSettings().collect_traces
    /// True
    #[getter]
    fn collect_traces(&self) -> bool {
        self.collect_traces
    }
    /// Requested tensor-factor contraction order. None uses the default order.
    ///
    /// Returns
    /// -------
    /// list of int or None
    ///     Requested tensor-factor contraction order. None uses the default order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings
    /// >>> ContractSettings().order is None
    /// True
    #[getter]
    fn order(&self) -> Option<Vec<usize>> {
        self.order.clone()
    }
    /// Maximum successful transformations in each root or shared-definition domain. Zero is allowed.
    ///
    /// Returns
    /// -------
    /// int
    ///     Maximum successful transformations in each root or shared-definition domain. Zero is allowed.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import ContractSettings
    /// >>> ContractSettings().max_passes
    /// 16
    #[getter]
    fn max_passes(&self) -> usize {
        self.max_passes
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> text = repr(settings)
    fn __repr__(self_: PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = self_.into_pyobject(py)?;
        // Filters retain representation identities, so the placeholder dimension
        // only reconstructs the already registered name; no dimension is filtered.
        let representations = object
            .borrow()
            .representations
            .as_ref()
            .map(|values| {
                let names = values
                    .iter()
                    .map(|rep| {
                        format!(
                            "Representation({:?}, 1, is_self_dual={}){}.name",
                            rep.base().symbol().get_name(),
                            if rep.is_self_dual() { "True" } else { "False" },
                            if !rep.is_self_dual() && rep.is_dual() {
                                ".dual()"
                            } else {
                                ""
                            },
                        )
                    })
                    .collect::<Vec<_>>();
                format!("[{}]", names.join(", "))
            })
            .unwrap_or_else(|| "None".to_owned());
        let mut arguments = vec![format!("representations={representations}")];
        for field in [
            "metrics",
            "rank_one",
            "collect_chains",
            "collect_traces",
            "order",
            "max_passes",
        ] {
            arguments.push(format!("{field}={}", object.getattr(field)?.repr()?));
        }
        Ok(format!("ContractSettings({})", arguments.join(", ")))
    }
}

/// Completion of a tensor reduction under its selected settings.
///
/// Complete means no eligible work remains. Deferred means exact eligible work
/// remains for another operation. Capped means the transformation budget stopped
/// processing. A partial result still represents the exact original tensor;
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
    /// The work budget stopped processing; the remaining expression is exact.
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
