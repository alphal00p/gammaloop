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
    symbolica.community.tensor,
    CanonicalizationError,
    PyValueError,
    python_doc!("CanonicalizationError")
);
create_exception!(
    symbolica.community.tensor,
    CookingError,
    PyTypeError,
    python_doc!("CookingError")
);
create_exception!(
    symbolica.community.tensor,
    DiracAdjointError,
    PyValueError,
    python_doc!("DiracAdjointError")
);
create_exception!(
    symbolica.community.tensor,
    NetworkToolingError,
    PyValueError,
    python_doc!("NetworkToolingError")
);

/// Encoding used when compound index labels become symbols.
///
/// ``FlattenedSymbol`` creates readable names from the function and arguments.
/// ``ReversibleEncoding`` retains the payload in a stable encoding. These are
/// encoding choices; the Python API does not currently expose an uncook method.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import CookMode, CookSettings
/// >>> settings = CookSettings(mode=CookMode.ReversibleEncoding)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "CookMode",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.tensor"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCookMode {
    /// Build a readable symbol name from the function name and its arguments.
    FlattenedSymbol,
    /// Retain the original payload in a stable symbol encoding.
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

/// Select function heads for index cooking by their Symbolica tags.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import CookTagFilter
/// >>> selection = CookTagFilter.any(["example::index_label"])
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookTagFilter",
    from_py_object,
    module = "symbolica.community.tensor"
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
#[spenso_macros::track_usage(crate::record_usage)]
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

    /// Select function heads carrying at least one listed tag.
    ///
    /// Parameters
    /// ----------
    /// tags : sequence of str
    ///     Fully namespaced Symbolica tags to test.
    ///
    /// Returns
    /// -------
    /// CookTagFilter
    ///     Tag predicate for a CookSourceFilter.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookTagFilter
    /// >>> selection = CookTagFilter.any(["example::index_label"])
    #[staticmethod]
    pub(crate) fn any(tags: Vec<String>) -> Self {
        Self {
            inner: RustCookTagFilter::any(tags),
        }
    }

    /// Select function heads carrying every listed tag.
    ///
    /// Parameters
    /// ----------
    /// tags : sequence of str
    ///     Fully namespaced Symbolica tags to test.
    ///
    /// Returns
    /// -------
    /// CookTagFilter
    ///     Tag predicate for a CookSourceFilter.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookTagFilter
    /// >>> selection = CookTagFilter.all(["example::index_label"])
    #[staticmethod]
    pub(crate) fn all(tags: Vec<String>) -> Self {
        Self {
            inner: RustCookTagFilter::all(tags),
        }
    }

    /// Select function heads using the associated settings' output tags.
    ///
    /// Returns
    /// -------
    /// CookTagFilter
    ///     Predicate resolved against CookSettings.output_tags.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookTagFilter
    /// >>> selection = CookTagFilter.matched_output_tags()
    #[staticmethod]
    pub(crate) fn matched_output_tags() -> Self {
        Self {
            inner: RustCookTagFilter::MatchedOutputTags,
        }
    }
}

/// Choose which function payloads are eligible for index cooking.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import CookSourceFilter
/// >>> selection = CookSourceFilter.representation_index_payload()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookSourceFilter",
    from_py_object,
    module = "symbolica.community.tensor"
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
#[spenso_macros::track_usage(crate::record_usage)]
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

    /// Allow every function-like payload encountered by the cooking operation.
    ///
    /// Returns
    /// -------
    /// CookSourceFilter
    ///     Unrestricted function selection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSourceFilter
    /// >>> selection = CookSourceFilter.any_function()
    #[staticmethod]
    pub(crate) fn any_function() -> Self {
        Self {
            inner: RustCookSourceFilter::AnyFunction,
        }
    }

    /// Select function payloads whose heads satisfy a tag filter.
    ///
    /// Parameters
    /// ----------
    /// filter : CookTagFilter
    ///     Required tag predicate.
    ///
    /// Returns
    /// -------
    /// CookSourceFilter
    ///     Function-head selection for CookSettings.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSourceFilter, CookTagFilter
    /// >>> selection = CookSourceFilter.function_tags(CookTagFilter.any(["example::index_label"]))
    #[staticmethod]
    pub(crate) fn function_tags(filter: &PyCookTagFilter) -> Self {
        Self {
            inner: RustCookSourceFilter::FunctionTags(filter.rust()),
        }
    }

    /// Select compound index labels inside representation syntax.
    ///
    /// Parameters
    /// ----------
    /// filter : CookTagFilter, optional
    ///     Additional predicate on the index payload's function head.
    ///
    /// Returns
    /// -------
    /// CookSourceFilter
    ///     Selection restricted to index payloads, excluding dimensions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSourceFilter
    /// >>> selection = CookSourceFilter.representation_index_payload()
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

/// Encoding options for compound tensor index labels.
///
/// Instances are immutable. Pass a configured instance to the corresponding
/// tensor operation; its properties expose the selected values.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import CookSettings
/// >>> settings = CookSettings()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookSettings",
    from_py_object,
    module = "symbolica.community.tensor"
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
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
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

    /// Configure encoding options for compound tensor index labels.
    ///
    /// Parameters
    /// ----------
    /// mode : CookMode or None, default None
    ///     FlattenedSymbol makes readable names. ReversibleEncoding stores the
    ///     original payload in a stable encoding. None selects FlattenedSymbol.
    /// source : CookSourceFilter or None, default None
    ///     Predicate selecting function payloads to encode. None selects any function.
    /// output_tags : sequence of str or None, default None
    ///     Fully namespaced tags for encoded symbols. None uses the mode factory
    ///     defaults; the generic constructor starts with no explicit output tags.
    /// preserve_tags : bool, default False
    ///     Preserve matching input tags on newly encoded symbols.
    ///
    /// Returns
    /// -------
    /// CookSettings
    ///     Immutable settings with the requested options.
    ///
    /// Notes
    /// -----
    /// Pass these settings through cook_indices in tensor constructors or indexing.
    /// Cooking replaces compound labels by symbols; it is not tensor simplification.
    /// The Python API currently exposes no uncook method.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings()
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

    /// Use readable names for selected compound labels.
    ///
    /// Returns
    /// -------
    /// CookSettings
    ///     Immutable configuration for index cooking.
    ///
    /// Notes
    /// -----
    /// Selects any function payload, without explicit output tags or tag preservation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings.flattened()
    #[staticmethod]
    pub(crate) fn flattened() -> Self {
        Self {
            inner: RustCookSettings::flattened(),
        }
    }

    /// Encode nested function payloads in indices while preserving tags.
    ///
    /// Returns
    /// -------
    /// CookSettings
    ///     Immutable configuration for index cooking.
    ///
    /// Notes
    /// -----
    /// Uses readable flattened names and selects representation index payloads,
    /// not representation dimensions. This is the usual cook_indices setting.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings.indices()
    #[staticmethod]
    pub(crate) fn indices() -> Self {
        Self {
            inner: RustCookSettings::indices(),
        }
    }

    /// Use a stable encoding that retains the original function payload.
    ///
    /// Returns
    /// -------
    /// CookSettings
    ///     Immutable configuration for index cooking.
    ///
    /// Notes
    /// -----
    /// Selects any function and applies the idenso::cooked output tag. The name
    /// refers to the encoding format; there is no public Python uncook method.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings.reversible()
    #[staticmethod]
    pub(crate) fn reversible() -> Self {
        Self {
            inner: RustCookSettings::reversible(),
        }
    }

    /// FlattenedSymbol makes readable names. ReversibleEncoding stores the
    ///
    /// Returns
    /// -------
    /// CookMode
    ///     FlattenedSymbol makes readable names. ReversibleEncoding stores the
    ///     original payload in a stable encoding. None selects FlattenedSymbol.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings()
    /// >>> value = settings.mode
    #[getter]
    pub(crate) fn mode(&self) -> PyCookMode {
        self.inner.mode().into()
    }

    /// Predicate selecting function payloads to encode. None selects any function.
    ///
    /// Returns
    /// -------
    /// CookSourceFilter
    ///     Predicate selecting function payloads to encode. None selects any function.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings()
    /// >>> value = settings.source_filter
    #[getter]
    pub(crate) fn source_filter(&self) -> PyCookSourceFilter {
        PyCookSourceFilter {
            inner: self.inner.source_filter().clone(),
        }
    }

    /// Fully namespaced tags for encoded symbols. None uses the mode factory
    ///
    /// Returns
    /// -------
    /// list of str
    ///     Fully namespaced tags for encoded symbols. The generic constructor starts with no explicit output tags.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings()
    /// >>> value = settings.output_tags
    #[getter]
    pub(crate) fn output_tags(&self) -> Vec<String> {
        self.inner.output_tags().to_vec()
    }

    /// Preserve matching input tags on newly encoded symbols.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Preserve matching input tags on newly encoded symbols.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import CookSettings
    /// >>> settings = CookSettings()
    /// >>> value = settings.preserve_tags
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
