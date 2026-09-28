use idenso::{
    CookMode as RustCookMode, CookSettings as RustCookSettings,
    CookSourceFilter as RustCookSourceFilter, CookTagFilter as RustCookTagFilter,
};
#[cfg(not(feature = "python_stubgen"))]
use pyo3::create_exception;
use pyo3::{
    Bound, IntoPyObject, PyResult, Python,
    exceptions::{PyTypeError, PyValueError},
    pyclass, pymethods,
    types::{PyAnyMethods, PyModule, PyModuleMethods},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::create_exception;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

create_exception!(
    symbolica.community.spenso,
    CanonicalizationError,
    PyValueError,
    "Raised when dummy-index canonicalization fails."
);
create_exception!(
    symbolica.community.spenso,
    CookingError,
    PyTypeError,
    "Raised when a symbolic function or representation index cannot be cooked."
);
create_exception!(
    symbolica.community.spenso,
    DiracAdjointError,
    PyValueError,
    "Raised when a Dirac adjoint cannot be constructed consistently."
);
create_exception!(
    symbolica.community.spenso,
    NetworkToolingError,
    PyValueError,
    "Raised when a symbolic tensor network cannot be parsed or evaluated."
);

/// Selects how cooked function payloads are represented as symbols.
///
/// Available values are `FlattenedSymbol` for readable names and `ReversibleEncoding` for a
/// stable encoding that can later be restored by `uncook`.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "CookMode",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCookMode {
    /// Build a readable symbol name from the function name and its arguments.
    FlattenedSymbol,
    /// Store a stable encoding that can later be restored by `uncook` with matching settings.
    ReversibleEncoding,
}

impl From<PyCookMode> for RustCookMode {
    fn from(value: PyCookMode) -> Self {
        match value {
            PyCookMode::FlattenedSymbol => Self::FlattenedSymbol,
            PyCookMode::ReversibleEncoding => Self::ReversibleEncoding,
        }
    }
}

impl From<RustCookMode> for PyCookMode {
    fn from(value: RustCookMode) -> Self {
        match value {
            RustCookMode::FlattenedSymbol => Self::FlattenedSymbol,
            RustCookMode::ReversibleEncoding => Self::ReversibleEncoding,
        }
    }
}

/// A tag predicate used to select which function heads are cooked.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookTagFilter",
    from_py_object,
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub(crate) struct PyCookTagFilter {
    inner: RustCookTagFilter,
}

impl PyCookTagFilter {
    pub(crate) fn rust(&self) -> RustCookTagFilter {
        self.inner.clone()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCookTagFilter {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let (factory, tags) = match &self.inner {
            RustCookTagFilter::Any(tags) => ("any", tags),
            RustCookTagFilter::All(tags) => ("all", tags),
            RustCookTagFilter::MatchedOutputTags => {
                return Ok("CookTagFilter.matched_output_tags()".into());
            }
        };
        Ok(format!(
            "CookTagFilter.{factory}({})",
            tags.into_pyobject(py)?.repr()?
        ))
    }

    /// Match a function head when it carries at least one listed tag.
    #[staticmethod]
    pub(crate) fn any(tags: Vec<String>) -> Self {
        Self {
            inner: RustCookTagFilter::any(tags),
        }
    }

    /// Match a function head only when it carries every listed tag.
    #[staticmethod]
    pub(crate) fn all(tags: Vec<String>) -> Self {
        Self {
            inner: RustCookTagFilter::all(tags),
        }
    }

    /// Match the explicit output tags configured on the associated `CookSettings`.
    #[staticmethod]
    pub(crate) fn matched_output_tags() -> Self {
        Self {
            inner: RustCookTagFilter::MatchedOutputTags,
        }
    }
}

/// Selects the function occurrences or representation-index payloads to cook.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookSourceFilter",
    from_py_object,
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub(crate) struct PyCookSourceFilter {
    inner: RustCookSourceFilter,
}

impl PyCookSourceFilter {
    pub(crate) fn rust(&self) -> RustCookSourceFilter {
        self.inner.clone()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCookSourceFilter {
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(match &self.inner {
            RustCookSourceFilter::AnyFunction => "CookSourceFilter.any_function()".into(),
            RustCookSourceFilter::FunctionTags(filter) => format!(
                "CookSourceFilter.function_tags({})",
                PyCookTagFilter {
                    inner: filter.clone()
                }
                .__repr__(py)?
            ),
            RustCookSourceFilter::RepresentationPayload { filter, .. } => format!(
                "CookSourceFilter.representation_index_payload({})",
                filter
                    .as_ref()
                    .map(|filter| PyCookTagFilter {
                        inner: filter.clone()
                    }
                    .__repr__(py))
                    .transpose()?
                    .unwrap_or_else(|| "None".into())
            ),
        })
    }

    /// Select every function-like subexpression.
    #[staticmethod]
    pub(crate) fn any_function() -> Self {
        Self {
            inner: RustCookSourceFilter::AnyFunction,
        }
    }

    /// Select function heads accepted by `filter`.
    #[staticmethod]
    pub(crate) fn function_tags(filter: &PyCookTagFilter) -> Self {
        Self {
            inner: RustCookSourceFilter::FunctionTags(filter.rust()),
        }
    }

    /// Select only function payloads inside representation indices.
    ///
    /// When `filter` is supplied, the payload's function head must also match it.
    #[staticmethod]
    #[pyo3(signature = (filter = None))]
    pub(crate) fn representation_index_payload(filter: Option<&PyCookTagFilter>) -> Self {
        Self {
            inner: RustCookSourceFilter::RepresentationPayload {
                indices: true,
                dimensions: false,
                filter: filter.map(PyCookTagFilter::rust),
            },
        }
    }
}

/// Immutable configuration for cooking symbolic functions and index payloads.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookSettings",
    from_py_object,
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub(crate) struct PyCookSettings {
    inner: RustCookSettings,
}

impl PyCookSettings {
    /// Construct settings from explicit Rust-side values.
    pub(crate) fn new(
        mode: PyCookMode,
        source: Option<&PyCookSourceFilter>,
        output_tags: Option<Vec<String>>,
        preserve_tags: bool,
    ) -> Self {
        let mut inner = RustCookSettings::flattened().with_mode(mode.into());
        if let Some(source) = source {
            inner = inner.with_source_filter(source.rust());
        }
        if let Some(output_tags) = output_tags {
            inner = inner.with_output_tags(output_tags);
        }
        if preserve_tags {
            inner = inner.preserve_tags();
        }
        Self { inner }
    }

    pub(crate) fn rust(&self) -> RustCookSettings {
        self.inner.clone()
    }

    pub(crate) fn indices_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(Self::rust)
            .unwrap_or_else(RustCookSettings::indices)
    }

    pub(crate) fn flattened_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(Self::rust)
            .unwrap_or_else(RustCookSettings::flattened)
    }

    pub(crate) fn reversible_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(Self::rust)
            .unwrap_or_else(RustCookSettings::reversible)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCookSettings {
    fn __repr__(self_: pyo3::PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = pyo3::IntoPyObject::into_pyobject(self_, py)?;
        crate::display::constructor_repr(
            object.as_any(),
            &[
                ("mode", "mode"),
                ("source", "source_filter"),
                ("output_tags", "output_tags"),
                ("preserve_tags", "preserve_tags"),
            ],
        )
    }

    /// Configure how functions are selected, encoded, and tagged when cooked.
    ///
    /// Tags must be fully namespaced Symbolica tags, for example `idenso::cooked`.
    /// Selecting `ReversibleEncoding` changes only the encoding; use `reversible()` to also
    /// select the conventional `idenso::cooked` output tag.
    /// `mode=None` selects `CookMode.FlattenedSymbol`.
    #[new]
    #[pyo3(
        signature = (
            *,
            mode = None,
            source = None,
            output_tags = None,
            preserve_tags = false
        ),
        text_signature = "(*, mode=None, source=None, output_tags=None, preserve_tags=False)"
    )]
    pub(crate) fn py_new(
        mode: Option<PyCookMode>,
        source: Option<&PyCookSourceFilter>,
        output_tags: Option<Vec<String>>,
        preserve_tags: bool,
    ) -> Self {
        Self::new(
            mode.unwrap_or(PyCookMode::FlattenedSymbol),
            source,
            output_tags,
            preserve_tags,
        )
    }

    /// Cook all functions into readable flattened names.
    #[staticmethod]
    pub(crate) fn flattened() -> Self {
        Self {
            inner: RustCookSettings::flattened(),
        }
    }

    /// Cook only nested function payloads inside representation indices and preserve tags.
    #[staticmethod]
    pub(crate) fn indices() -> Self {
        Self {
            inner: RustCookSettings::indices(),
        }
    }

    /// Cook all functions into stable symbols that can be restored by `uncook`.
    #[staticmethod]
    pub(crate) fn reversible() -> Self {
        Self {
            inner: RustCookSettings::reversible(),
        }
    }

    /// Symbol-encoding mode.
    #[getter]
    pub(crate) fn mode(&self) -> PyCookMode {
        self.inner.mode().into()
    }

    /// Function occurrences selected for cooking.
    #[getter]
    pub(crate) fn source_filter(&self) -> PyCookSourceFilter {
        PyCookSourceFilter {
            inner: self.inner.source_filter().clone(),
        }
    }

    /// Explicit tags attached to newly created cooked symbols.
    #[getter]
    pub(crate) fn output_tags(&self) -> Vec<String> {
        self.inner.output_tags().to_vec()
    }

    /// Whether matching input tags are preserved on cooked symbols.
    #[getter]
    pub(crate) fn preserve_tags(&self) -> bool {
        self.inner.preserves_tags()
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "CanonicalizationError",
        module.py().get_type::<CanonicalizationError>(),
    )?;
    module.add("CookingError", module.py().get_type::<CookingError>())?;
    module.add(
        "DiracAdjointError",
        module.py().get_type::<DiracAdjointError>(),
    )?;
    module.add(
        "NetworkToolingError",
        module.py().get_type::<NetworkToolingError>(),
    )?;

    module.add_class::<PyCookMode>()?;
    module.add_class::<PyCookTagFilter>()?;
    module.add_class::<PyCookSourceFilter>()?;
    module.add_class::<PyCookSettings>()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_cook_settings_match_rust_defaults() {
        assert_eq!(
            PyCookSettings::new(PyCookMode::FlattenedSymbol, None, None, false).rust(),
            RustCookSettings::default()
        );
    }

    #[test]
    fn network_tooling_exception_is_a_value_error() {
        Python::initialize();
        Python::attach(|py| {
            assert!(
                CanonicalizationError::new_err("canonicalization failure")
                    .is_instance_of::<PyValueError>(py)
            );
            assert!(
                DiracAdjointError::new_err("adjoint failure").is_instance_of::<PyValueError>(py)
            );
            assert!(
                NetworkToolingError::new_err("network failure").is_instance_of::<PyValueError>(py)
            );
        });
    }
}
