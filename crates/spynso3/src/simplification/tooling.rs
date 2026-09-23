use idenso::{
    CookMode as RustCookMode, CookSettings as RustCookSettings,
    CookSourceFilter as RustCookSourceFilter, CookTagFilter as RustCookTagFilter,
    shorthands::schoonschip::{
        SchoonschipContractionOrder as RustSchoonschipContractionOrder,
        SchoonschipSettings as RustSchoonschipSettings,
    },
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
    DotExpansionError,
    PyValueError,
    "Raised when compact dot notation cannot be expanded into a tensor expression."
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
            RustCookSourceFilter::RepresentationIndexPayload { filter } => format!(
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
            inner: RustCookSourceFilter::RepresentationIndexPayload {
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

/// Selects whether a Schoonschip pass runs once or recursively.
///
/// Available values are `SinglePass` and `Recursive`.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "SchoonschipMode",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySchoonschipMode {
    /// Visit each eligible expression at most once.
    SinglePass,
    /// Repeat traversal until the configured depth or a fixed point is reached.
    Recursive,
}

/// Selects the recursive traversal order for a Schoonschip pass.
///
/// Available values are `DepthFirst` and `BreadthFirst`.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "SchoonschipTraversal",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySchoonschipTraversal {
    /// Fully simplify each branch before advancing to its siblings.
    DepthFirst,
    /// Advance all branches one level before descending further.
    BreadthFirst,
}

/// Selects the heuristic used to choose the next tensor-network contraction.
///
/// Available values are `SmallestDegree`, `LargestDegree`, `MinLargestOperandBytes`,
/// `MinProductTerms`, `MinProductBytes`, `SmallestDegreeMinLargestOperandBytes`,
/// `SmallestDegreeMinProductTerms`, and `SmallestDegreeMinProductBytes`.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "SchoonschipContractionOrder",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.spenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySchoonschipContractionOrder {
    /// Contract the pair with the fewest paired tensor slots.
    SmallestDegree,
    /// Contract the pair with the most paired tensor slots.
    LargestDegree,
    /// Minimize the larger operand's estimated memory footprint.
    MinLargestOperandBytes,
    /// Minimize the estimated number of terms in the product.
    MinProductTerms,
    /// Minimize the product's estimated memory footprint.
    MinProductBytes,
    /// Minimize the paired-slot count first, then the larger operand's estimated bytes.
    SmallestDegreeMinLargestOperandBytes,
    /// Minimize the paired-slot count first, then the estimated product term count.
    SmallestDegreeMinProductTerms,
    /// Minimize the paired-slot count first, then the estimated product bytes.
    SmallestDegreeMinProductBytes,
}

impl From<PySchoonschipContractionOrder> for RustSchoonschipContractionOrder {
    fn from(value: PySchoonschipContractionOrder) -> Self {
        match value {
            PySchoonschipContractionOrder::SmallestDegree => Self::SmallestDegree,
            PySchoonschipContractionOrder::LargestDegree => Self::LargestDegree,
            PySchoonschipContractionOrder::MinLargestOperandBytes => Self::MinLargestOperandBytes,
            PySchoonschipContractionOrder::MinProductTerms => Self::MinProductTerms,
            PySchoonschipContractionOrder::MinProductBytes => Self::MinProductBytes,
            PySchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes => {
                Self::SmallestDegreeMinLargestOperandBytes
            }
            PySchoonschipContractionOrder::SmallestDegreeMinProductTerms => {
                Self::SmallestDegreeMinProductTerms
            }
            PySchoonschipContractionOrder::SmallestDegreeMinProductBytes => {
                Self::SmallestDegreeMinProductBytes
            }
        }
    }
}

impl From<RustSchoonschipContractionOrder> for PySchoonschipContractionOrder {
    fn from(value: RustSchoonschipContractionOrder) -> Self {
        match value {
            RustSchoonschipContractionOrder::SmallestDegree => Self::SmallestDegree,
            RustSchoonschipContractionOrder::LargestDegree => Self::LargestDegree,
            RustSchoonschipContractionOrder::MinLargestOperandBytes => Self::MinLargestOperandBytes,
            RustSchoonschipContractionOrder::MinProductTerms => Self::MinProductTerms,
            RustSchoonschipContractionOrder::MinProductBytes => Self::MinProductBytes,
            RustSchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes => {
                Self::SmallestDegreeMinLargestOperandBytes
            }
            RustSchoonschipContractionOrder::SmallestDegreeMinProductTerms => {
                Self::SmallestDegreeMinProductTerms
            }
            RustSchoonschipContractionOrder::SmallestDegreeMinProductBytes => {
                Self::SmallestDegreeMinProductBytes
            }
        }
    }
}

/// Immutable configuration for expression and network Schoonschip passes.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "SchoonschipSettings",
    from_py_object,
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub(crate) struct PySchoonschipSettings {
    depth_limit: Option<usize>,
    mode: PySchoonschipMode,
    traversal: PySchoonschipTraversal,
    expand_contracted_sums: bool,
    simplify_chain_like_functions: bool,
    schoonschip_rank1_tensors: bool,
    contraction_order: PySchoonschipContractionOrder,
}

impl PySchoonschipSettings {
    /// Construct settings from explicit Rust-side values.
    pub(crate) fn new(
        depth_limit: Option<usize>,
        mode: PySchoonschipMode,
        traversal: PySchoonschipTraversal,
        expand_contracted_sums: bool,
        simplify_chain_like_functions: bool,
        schoonschip_rank1_tensors: bool,
        contraction_order: PySchoonschipContractionOrder,
    ) -> Self {
        Self {
            depth_limit,
            mode,
            traversal,
            expand_contracted_sums,
            simplify_chain_like_functions,
            schoonschip_rank1_tensors,
            contraction_order,
        }
    }

    pub(crate) fn rust(&self) -> RustSchoonschipSettings {
        let mut settings = match (self.mode, self.traversal) {
            (PySchoonschipMode::SinglePass, _) => {
                RustSchoonschipSettings::single_pass(self.depth_limit)
            }
            (PySchoonschipMode::Recursive, PySchoonschipTraversal::DepthFirst) => {
                RustSchoonschipSettings::depth_first(self.depth_limit)
            }
            (PySchoonschipMode::Recursive, PySchoonschipTraversal::BreadthFirst) => {
                RustSchoonschipSettings::breadth_first(self.depth_limit)
            }
        };
        if self.expand_contracted_sums {
            settings = settings.with_expanded_contracted_sums();
        }
        if self.simplify_chain_like_functions {
            settings = settings.with_chain_like_functions();
        }
        if self.schoonschip_rank1_tensors {
            settings = settings.with_rank1_tensors();
        } else {
            settings = settings.without_rank1_tensors();
        }
        settings.with_contraction_order(self.contraction_order.into())
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySchoonschipSettings {
    fn __repr__(self_: pyo3::PyRef<'_, Self>) -> PyResult<String> {
        let py = self_.py();
        let object = pyo3::IntoPyObject::into_pyobject(self_, py)?;
        crate::display::constructor_repr(
            object.as_any(),
            &[
                ("depth_limit", "depth_limit"),
                ("mode", "mode"),
                ("traversal", "traversal"),
                ("expand_contracted_sums", "expand_contracted_sums"),
                (
                    "simplify_chain_like_functions",
                    "simplify_chain_like_functions",
                ),
                ("schoonschip_rank1_tensors", "schoonschip_rank1_tensors"),
                ("contraction_order", "contraction_order"),
            ],
        )
    }

    /// Configure traversal, depth, shorthand expansion, and network-contraction policies.
    ///
    /// `depth_limit=None` removes the recursion-depth limit. `traversal` is ignored in
    /// `SinglePass` mode.
    /// `mode=None`, `traversal=None`, and `contraction_order=None` select `Recursive`,
    /// `BreadthFirst`, and `SmallestDegree`, respectively.
    #[new]
    #[pyo3(
        signature = (
            *,
            depth_limit = Some(1),
            mode = None,
            traversal = None,
            expand_contracted_sums = false,
            simplify_chain_like_functions = false,
            schoonschip_rank1_tensors = true,
            contraction_order = None
        ),
        text_signature = "(*, depth_limit=1, mode=None, traversal=None, expand_contracted_sums=False, simplify_chain_like_functions=False, schoonschip_rank1_tensors=True, contraction_order=None)"
    )]
    pub(crate) fn py_new(
        depth_limit: Option<usize>,
        mode: Option<PySchoonschipMode>,
        traversal: Option<PySchoonschipTraversal>,
        expand_contracted_sums: bool,
        simplify_chain_like_functions: bool,
        schoonschip_rank1_tensors: bool,
        contraction_order: Option<PySchoonschipContractionOrder>,
    ) -> Self {
        Self::new(
            depth_limit,
            mode.unwrap_or(PySchoonschipMode::Recursive),
            traversal.unwrap_or(PySchoonschipTraversal::BreadthFirst),
            expand_contracted_sums,
            simplify_chain_like_functions,
            schoonschip_rank1_tensors,
            contraction_order.unwrap_or(PySchoonschipContractionOrder::SmallestDegree),
        )
    }

    /// Apply the default shallow recursive expression pass without rank-one tensors.
    #[staticmethod]
    pub(crate) fn partial() -> Self {
        Self::new(
            Some(1),
            PySchoonschipMode::Recursive,
            PySchoonschipTraversal::BreadthFirst,
            false,
            false,
            false,
            PySchoonschipContractionOrder::SmallestDegree,
        )
    }

    /// Apply one unrestricted-depth pass and include rank-one tensors.
    #[staticmethod]
    pub(crate) fn full() -> Self {
        Self::new(
            None,
            PySchoonschipMode::SinglePass,
            PySchoonschipTraversal::DepthFirst,
            false,
            false,
            true,
            PySchoonschipContractionOrder::SmallestDegree,
        )
    }

    /// Use the settings applied by `schoonschip_net` when no settings are supplied.
    #[staticmethod]
    pub(crate) fn default_network() -> Self {
        Self::partial()
    }

    /// Recursively simplify each branch before visiting its siblings.
    #[staticmethod]
    #[pyo3(signature = (depth_limit = None))]
    pub(crate) fn depth_first(depth_limit: Option<usize>) -> Self {
        Self::new(
            depth_limit,
            PySchoonschipMode::Recursive,
            PySchoonschipTraversal::DepthFirst,
            false,
            false,
            true,
            PySchoonschipContractionOrder::SmallestDegree,
        )
    }

    /// Recursively simplify all branches one level at a time.
    #[staticmethod]
    #[pyo3(signature = (depth_limit = None))]
    pub(crate) fn breadth_first(depth_limit: Option<usize>) -> Self {
        Self::new(
            depth_limit,
            PySchoonschipMode::Recursive,
            PySchoonschipTraversal::BreadthFirst,
            false,
            false,
            false,
            PySchoonschipContractionOrder::SmallestDegree,
        )
    }

    /// Visit each eligible expression once, subject to `depth_limit`.
    #[staticmethod]
    #[pyo3(signature = (depth_limit = None))]
    pub(crate) fn single_pass(depth_limit: Option<usize>) -> Self {
        Self::new(
            depth_limit,
            PySchoonschipMode::SinglePass,
            PySchoonschipTraversal::DepthFirst,
            false,
            false,
            true,
            PySchoonschipContractionOrder::SmallestDegree,
        )
    }

    /// Maximum parsing or recursion depth, or `None` for no limit.
    #[getter]
    pub(crate) fn depth_limit(&self) -> Option<usize> {
        self.depth_limit
    }

    /// Whether the pass is single-pass or recursive.
    #[getter]
    pub(crate) fn mode(&self) -> PySchoonschipMode {
        self.mode
    }

    /// Recursive traversal order, or `None` in single-pass mode.
    #[getter]
    pub(crate) fn traversal(&self) -> Option<PySchoonschipTraversal> {
        match self.mode {
            PySchoonschipMode::SinglePass => None,
            PySchoonschipMode::Recursive => Some(self.traversal),
        }
    }

    /// Whether contracted sums are expanded before network execution.
    #[getter]
    pub(crate) fn expand_contracted_sums(&self) -> bool {
        self.expand_contracted_sums
    }

    /// Whether chain-like function payloads are simplified recursively.
    #[getter]
    pub(crate) fn simplify_chain_like_functions(&self) -> bool {
        self.simplify_chain_like_functions
    }

    /// Whether rank-one tensors participate in the Schoonschip pass.
    #[getter]
    pub(crate) fn schoonschip_rank1_tensors(&self) -> bool {
        self.schoonschip_rank1_tensors
    }

    /// Heuristic used to choose network contractions.
    #[getter]
    pub(crate) fn contraction_order(&self) -> PySchoonschipContractionOrder {
        self.contraction_order
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
        "DotExpansionError",
        module.py().get_type::<DotExpansionError>(),
    )?;
    module.add(
        "NetworkToolingError",
        module.py().get_type::<NetworkToolingError>(),
    )?;

    module.add_class::<PyCookMode>()?;
    module.add_class::<PyCookTagFilter>()?;
    module.add_class::<PyCookSourceFilter>()?;
    module.add_class::<PyCookSettings>()?;
    module.add_class::<PySchoonschipMode>()?;
    module.add_class::<PySchoonschipTraversal>()?;
    module.add_class::<PySchoonschipContractionOrder>()?;
    module.add_class::<PySchoonschipSettings>()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_schoonschip_settings_eq(
        left: RustSchoonschipSettings,
        right: RustSchoonschipSettings,
    ) {
        assert_eq!(left.depth_limit, right.depth_limit);
        assert!(left.mode == right.mode);
        assert_eq!(left.expand_contracted_sums, right.expand_contracted_sums);
        assert_eq!(
            left.simplify_chain_like_functions,
            right.simplify_chain_like_functions
        );
        assert_eq!(
            left.schoonschip_rank1_tensors,
            right.schoonschip_rank1_tensors
        );
        assert_eq!(left.contraction_order, right.contraction_order);
    }

    #[test]
    fn python_settings_match_rust_defaults_and_presets() {
        assert_eq!(
            PyCookSettings::new(PyCookMode::FlattenedSymbol, None, None, false).rust(),
            RustCookSettings::default()
        );
        assert_schoonschip_settings_eq(
            PySchoonschipSettings::new(
                Some(1),
                PySchoonschipMode::Recursive,
                PySchoonschipTraversal::BreadthFirst,
                false,
                false,
                true,
                PySchoonschipContractionOrder::SmallestDegree,
            )
            .rust(),
            RustSchoonschipSettings::default(),
        );
        assert_schoonschip_settings_eq(
            PySchoonschipSettings::default_network().rust(),
            RustSchoonschipSettings::default_network(),
        );
        assert_schoonschip_settings_eq(
            PySchoonschipSettings::full().rust(),
            RustSchoonschipSettings::full(),
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
