//! Python conversion for the shared tensor-algebra scheduler.
use super::{PyColorSimplifySettings, PyGammaSimplifySettings};
use crate::{metadata::SpensoRepresentationName, structure::SpensoRepresentation};
use idenso::tensor::{ContractSettings, simplification::AlgebraSettings};
use pyo3::prelude::*;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use spenso::structure::representation::{LibraryRep, RepName};

/// Immutable selection of tensor identities. Dimensions come from tensor slots.
/// Materialize scalar polynomials explicitly with the result's `expand()` method.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> settings = sp.AlgebraSettings(gamma=sp.GammaSimplifySettings.canonical())
/// >>> gamma = sp.TensorExpression.gamma(4)
/// >>> reduced = (gamma("a", "b", "mu") * gamma("b", "a", "nu")).reduce_algebra(settings)
///
/// Parameters
/// ----------
/// gamma : GammaSimplifySettings, optional
///     Gamma-algebra settings, or None to leave gamma identities disabled.
/// color : ColorSimplifySettings, optional
///     Color-algebra settings, or None to leave color identities disabled.
/// epsilon : bool
///     Enable Levi-Civita identities and the contractions needed to apply them.
/// max_passes : int
///     Maximum successful transformations per root or alias domain. Inspect .reduction_status when using a small limit; a capped result retains the exact unfinished expression.
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
    /// Select identities explicitly. An enabled family includes its necessary contractions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.AlgebraSettings(gamma=sp.GammaSimplifySettings.canonical())
    /// >>> gamma = sp.TensorExpression.gamma(4)
    /// >>> reduced = (gamma("a", "b", "mu") * gamma("b", "a", "nu")).reduce_algebra(settings)
    ///
    /// Parameters
    /// ----------
    /// gamma : GammaSimplifySettings, optional
    ///     Gamma-algebra settings, or None to leave gamma identities disabled.
    /// color : ColorSimplifySettings, optional
    ///     Color-algebra settings, or None to leave color identities disabled.
    /// epsilon : bool
    ///     Enable Levi-Civita identities and the contractions needed to apply them.
    /// max_passes : int
    ///     Maximum successful transformations per root or alias domain. Inspect .reduction_status when using a small limit; a capped result retains the exact unfinished expression.
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

    /// Enable gamma, color, and epsilon identities with their native defaults.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.AlgebraSettings.hep()
    /// >>> assert settings.epsilon
    #[staticmethod]
    fn hep() -> Self {
        Self {
            inner: AlgebraSettings::hep(),
        }
    }

    /// Gamma-algebra settings, or None to leave gamma identities disabled.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.AlgebraSettings()
    /// >>> result = settings.gamma
    #[getter]
    fn gamma(&self) -> Option<PyGammaSimplifySettings> {
        self.inner.gamma.map(PyGammaSimplifySettings::from_rust)
    }
    /// Color-algebra settings, or None to leave color identities disabled.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.AlgebraSettings()
    /// >>> result = settings.color
    #[getter]
    fn color(&self) -> Option<PyColorSimplifySettings> {
        self.inner.color.map(PyColorSimplifySettings::from_rust)
    }
    /// Enable Levi-Civita identities and the contractions needed to apply them.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.AlgebraSettings()
    /// >>> result = settings.epsilon
    #[getter]
    fn epsilon(&self) -> bool {
        self.inner.epsilon
    }
    /// Maximum successful transformations per root or alias domain. Inspect .reduction_status when using a small limit; a capped result retains the exact unfinished expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.AlgebraSettings()
    /// >>> result = settings.max_passes
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

/// Structural contraction policy. Filters select connecting representation names,
/// independently of other ports on the same tensor. None permits every representation;
/// an empty sequence permits none. Intrinsic constructor normalization is unaffected.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> settings = sp.ContractSettings(collect_chains=False, collect_traces=False)
/// >>> contracted = sp.TensorExpression.g(r)("i", "i").contract(settings)
///
/// Parameters
/// ----------
/// representations : Sequence[Representation | RepresentationName], optional
///     Restrict contraction and collection to these Representation or RepresentationName values. None permits all representations.
/// metrics : bool
///     Contract compatible metric tensors. Disable to retain explicit metric factors.
/// rank_one : bool
///     Contract compatible rank-one tensors and recognize their scalar products.
/// collect_chains : bool
///     Collect compatible connected matrix products into compact chains. Disable to retain their indexed factors.
/// collect_traces : bool
///     Collect compatible closed matrix products into compact traces. This does not request gamma or color trace identities.
/// order : Sequence[int], optional
///     Optional tensor-factor contraction order. None lets Spenso choose an order; this does not permute the public tensor axes.
/// max_passes : int
///     Maximum successful transformations per root or alias domain. Inspect .reduction_status when using a small limit; a capped result retains the exact unfinished expression.
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
    /// Choose which index contractions and compact tensor forms to produce.
    /// By default all representations are eligible, and metrics, rank-one contractions,
    /// chains, and traces are enabled. Set a flag to False to leave that operation untouched.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> settings = sp.ContractSettings(collect_chains=False, collect_traces=False)
    /// >>> contracted = sp.TensorExpression.g(r)("i", "i").contract(settings)
    ///
    /// Parameters
    /// ----------
    /// representations : Sequence[Representation | RepresentationName], optional
    ///     Restrict contraction and collection to these Representation or RepresentationName values. None permits all representations.
    /// metrics : bool
    ///     Contract compatible metric tensors. Disable to retain explicit metric factors.
    /// rank_one : bool
    ///     Contract compatible rank-one tensors and recognize their scalar products.
    /// collect_chains : bool
    ///     Collect compatible connected matrix products into compact chains. Disable to retain their indexed factors.
    /// collect_traces : bool
    ///     Collect compatible closed matrix products into compact traces. This does not request gamma or color trace identities.
    /// order : Sequence[int], optional
    ///     Optional tensor-factor contraction order. None lets Spenso choose an order; this does not permute the public tensor axes.
    /// max_passes : int
    ///     Maximum successful transformations per root or alias domain. Inspect .reduction_status when using a small limit; a capped result retains the exact unfinished expression.
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

    /// Restrict contraction and collection to these Representation or RepresentationName values. None permits all representations.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> result = settings.representations
    #[getter]
    fn representations(&self) -> Option<Vec<SpensoRepresentationName>> {
        self.representations.as_ref().map(|values| {
            values
                .iter()
                .map(|&rep| SpensoRepresentationName { rep })
                .collect()
        })
    }
    /// Contract compatible metric tensors. Disable to retain explicit metric factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> result = settings.metrics
    #[getter]
    fn metrics(&self) -> bool {
        self.metrics
    }
    /// Contract compatible rank-one tensors and recognize their scalar products.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> result = settings.rank_one
    #[getter]
    fn rank_one(&self) -> bool {
        self.rank_one
    }
    /// Collect compatible connected matrix products into compact chains. Disable to retain their indexed factors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> result = settings.collect_chains
    #[getter]
    fn collect_chains(&self) -> bool {
        self.collect_chains
    }
    /// Collect compatible closed matrix products into compact traces. This does not request gamma or color trace identities.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> result = settings.collect_traces
    #[getter]
    fn collect_traces(&self) -> bool {
        self.collect_traces
    }
    /// Optional tensor-factor contraction order. None lets Spenso choose an order; this does not permute the public tensor axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> result = settings.order
    #[getter]
    fn order(&self) -> Option<Vec<usize>> {
        self.order.clone()
    }
    /// Maximum successful transformations per root or alias domain. Inspect .reduction_status when using a small limit; a capped result retains the exact unfinished expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> settings = sp.ContractSettings()
    /// >>> result = settings.max_passes
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

/// Completion under the selected capabilities and work budget.
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
