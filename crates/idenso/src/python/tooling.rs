use std::panic::{AssertUnwindSafe, catch_unwind};

#[cfg(not(feature = "python_stubgen"))]
use pyo3::create_exception;
use pyo3::{
    Borrowed, Bound, FromPyObject, PyAny, PyErr, PyResult,
    exceptions::{PyTypeError, PyValueError},
    pyclass, pyfunction, pymethods,
    types::{PyAnyMethods, PyModule, PyModuleMethods},
    wrap_pyfunction,
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::create_exception;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};
use spenso::structure::{
    abstract_index::AbstractIndex,
    representation::{LibraryRep, Representation},
};
use symbolica::api::python::PythonExpression;

use crate::{
    CookMode as RustCookMode, CookSettings as RustCookSettings,
    CookSourceFilter as RustCookSourceFilter, CookTagFilter as RustCookTagFilter, IndexTooling,
    shorthands::{
        UndoShorthands,
        chain::Chain,
        metric::MetricSimplifier,
        schoonschip::{
            Schoonschip, SchoonschipContractionOrder as RustSchoonschipContractionOrder,
            SchoonschipSettings as RustSchoonschipSettings,
        },
    },
};

create_exception!(
    symbolica.community.idenso,
    CanonicalizationError,
    PyValueError
);
create_exception!(symbolica.community.idenso, CookingError, PyTypeError);
create_exception!(symbolica.community.idenso, DiracAdjointError, PyValueError);
create_exception!(symbolica.community.idenso, DotExpansionError, PyValueError);

/// A registered representation accepted from either Spynso or its symbolic atom.
pub struct RegisteredRepresentation(pub Representation<LibraryRep>);

impl RegisteredRepresentation {
    pub fn rep(&self) -> LibraryRep {
        self.0.rep
    }

    pub fn to_atom(&self) -> symbolica::atom::Atom {
        self.0.to_symbolic([])
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for RegisteredRepresentation {
    type Error = PyErr;

    fn extract(value: Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let expression = if let Ok(expression) = value.extract::<PythonExpression>() {
            expression
        } else {
            value
                .call_method0("to_expression")?
                .extract::<PythonExpression>()?
        };
        let view = expression.expr.as_view();
        if let Ok(representation) = Representation::<LibraryRep>::try_from(view) {
            return Ok(Self(representation));
        }
        Err(PyTypeError::new_err(
            "representation must be a Spynso Representation or its symbolic expression",
        ))
    }
}

#[cfg(feature = "python_stubgen")]
impl pyo3_stub_gen::PyStubType for RegisteredRepresentation {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        pyo3_stub_gen::TypeInfo::locally_defined(
            "Representation",
            "symbolica.community.spenso".into(),
        ) | PythonExpression::type_output()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "CookMode",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyCookMode {
    FlattenedSymbol,
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

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookTagFilter",
    from_py_object,
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PyCookTagFilter {
    inner: RustCookTagFilter,
}

impl PyCookTagFilter {
    pub fn rust(&self) -> RustCookTagFilter {
        self.inner.clone()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCookTagFilter {
    #[staticmethod]
    pub fn any(tags: Vec<String>) -> Self {
        Self {
            inner: RustCookTagFilter::any(tags),
        }
    }

    #[staticmethod]
    pub fn all(tags: Vec<String>) -> Self {
        Self {
            inner: RustCookTagFilter::all(tags),
        }
    }

    #[staticmethod]
    pub fn matched_output_tags() -> Self {
        Self {
            inner: RustCookTagFilter::MatchedOutputTags,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    frozen,
    name = "CookSourceFilter",
    from_py_object,
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PyCookSourceFilter {
    inner: RustCookSourceFilter,
}

impl PyCookSourceFilter {
    pub fn rust(&self) -> RustCookSourceFilter {
        self.inner.clone()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCookSourceFilter {
    #[staticmethod]
    pub fn any_function() -> Self {
        Self {
            inner: RustCookSourceFilter::AnyFunction,
        }
    }

    #[staticmethod]
    pub fn function_tags(filter: &PyCookTagFilter) -> Self {
        Self {
            inner: RustCookSourceFilter::FunctionTags(filter.rust()),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (filter = None))]
    pub fn representation_index_payload(filter: Option<&PyCookTagFilter>) -> Self {
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
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PyCookSettings {
    inner: RustCookSettings,
}

impl PyCookSettings {
    pub fn rust(&self) -> RustCookSettings {
        self.inner.clone()
    }

    pub fn indices_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(Self::rust)
            .unwrap_or_else(RustCookSettings::indices)
    }

    pub fn flattened_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(Self::rust)
            .unwrap_or_else(RustCookSettings::flattened)
    }

    pub fn reversible_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(Self::rust)
            .unwrap_or_else(RustCookSettings::reversible)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCookSettings {
    #[new]
    #[pyo3(signature = (
        *,
        mode = PyCookMode::FlattenedSymbol,
        source = None,
        output_tags = None,
        preserve_tags = false
    ))]
    pub fn new(
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

    #[staticmethod]
    pub fn flattened() -> Self {
        Self {
            inner: RustCookSettings::flattened(),
        }
    }

    #[staticmethod]
    pub fn indices() -> Self {
        Self {
            inner: RustCookSettings::indices(),
        }
    }

    #[staticmethod]
    pub fn reversible() -> Self {
        Self {
            inner: RustCookSettings::reversible(),
        }
    }

    #[getter]
    pub fn mode(&self) -> PyCookMode {
        self.inner.mode().into()
    }

    #[getter]
    pub fn source_filter(&self) -> PyCookSourceFilter {
        PyCookSourceFilter {
            inner: self.inner.source_filter().clone(),
        }
    }

    #[getter]
    pub fn output_tags(&self) -> Vec<String> {
        self.inner.output_tags().to_vec()
    }

    #[getter]
    pub fn preserve_tags(&self) -> bool {
        self.inner.preserves_tags()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "SchoonschipMode",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PySchoonschipMode {
    SinglePass,
    Recursive,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "SchoonschipTraversal",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PySchoonschipTraversal {
    DepthFirst,
    BreadthFirst,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    frozen,
    name = "SchoonschipContractionOrder",
    from_py_object,
    eq,
    eq_int,
    module = "symbolica.community.idenso"
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PySchoonschipContractionOrder {
    SmallestDegree,
    LargestDegree,
    MinLargestOperandBytes,
    MinProductTerms,
    MinProductBytes,
    SmallestDegreeMinLargestOperandBytes,
    SmallestDegreeMinProductTerms,
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
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PySchoonschipSettings {
    depth_limit: Option<usize>,
    mode: PySchoonschipMode,
    traversal: PySchoonschipTraversal,
    expand_contracted_sums: bool,
    simplify_chain_like_functions: bool,
    schoonschip_rank1_tensors: bool,
    contraction_order: PySchoonschipContractionOrder,
}

impl PySchoonschipSettings {
    pub fn rust(&self) -> RustSchoonschipSettings {
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
    #[new]
    #[pyo3(signature = (
        *,
        depth_limit = Some(1),
        mode = PySchoonschipMode::Recursive,
        traversal = PySchoonschipTraversal::BreadthFirst,
        expand_contracted_sums = false,
        simplify_chain_like_functions = false,
        schoonschip_rank1_tensors = true,
        contraction_order = PySchoonschipContractionOrder::SmallestDegree
    ))]
    pub fn new(
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

    #[staticmethod]
    pub fn partial() -> Self {
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

    #[staticmethod]
    pub fn full() -> Self {
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

    #[staticmethod]
    pub fn default_network() -> Self {
        Self::partial()
    }

    #[staticmethod]
    #[pyo3(signature = (depth_limit = None))]
    pub fn depth_first(depth_limit: Option<usize>) -> Self {
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

    #[staticmethod]
    #[pyo3(signature = (depth_limit = None))]
    pub fn breadth_first(depth_limit: Option<usize>) -> Self {
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

    #[staticmethod]
    #[pyo3(signature = (depth_limit = None))]
    pub fn single_pass(depth_limit: Option<usize>) -> Self {
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

    #[getter]
    pub fn depth_limit(&self) -> Option<usize> {
        self.depth_limit
    }

    #[getter]
    pub fn mode(&self) -> PySchoonschipMode {
        self.mode
    }

    #[getter]
    pub fn traversal(&self) -> Option<PySchoonschipTraversal> {
        match self.mode {
            PySchoonschipMode::SinglePass => None,
            PySchoonschipMode::Recursive => Some(self.traversal),
        }
    }

    #[getter]
    pub fn expand_contracted_sums(&self) -> bool {
        self.expand_contracted_sums
    }

    #[getter]
    pub fn simplify_chain_like_functions(&self) -> bool {
        self.simplify_chain_like_functions
    }

    #[getter]
    pub fn schoonschip_rank1_tensors(&self) -> bool {
        self.schoonschip_rank1_tensors
    }

    #[getter]
    pub fn contraction_order(&self) -> PySchoonschipContractionOrder {
        self.contraction_order
    }
}

fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| {
            payload
                .downcast_ref::<&str>()
                .map(|message| (*message).into())
        })
        .unwrap_or_else(|| "tensor canonicalization failed".into())
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn canonize(expression: &PythonExpression) -> PyResult<PythonExpression> {
    catch_unwind(AssertUnwindSafe(|| {
        expression
            .expr
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
    }))
    .map(Into::into)
    .map_err(|payload| CanonicalizationError::new_err(panic_message(payload)))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn alias_subtensors(
    expression: &PythonExpression,
    tensor_name: &str,
) -> (PythonExpression, Vec<(PythonExpression, PythonExpression)>) {
    let (root, aliases) = expression
        .expr
        .alias_subtensors(tensor_name)
        .into_inner_with_aliases();
    let mut aliases = aliases.into_iter().collect::<Vec<_>>();
    aliases.sort_by_cached_key(|(alias, _)| alias.to_string());
    (
        root.into(),
        aliases
            .into_iter()
            .map(|(alias, original)| (alias.into(), original.into()))
            .collect(),
    )
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn spenso_conjugate(expression: &PythonExpression) -> PythonExpression {
    expression.expr.spenso_conj().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn conjugate_transpose(
    expression: &PythonExpression,
    representation: RegisteredRepresentation,
) -> PythonExpression {
    expression
        .expr
        .conjugate_transpose(representation.rep())
        .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
pub fn cook(
    expression: &PythonExpression,
    settings: Option<&PyCookSettings>,
) -> PyResult<PythonExpression> {
    PyCookSettings::reversible_or(settings)
        .try_cook(expression.expr.as_view())
        .map(Into::into)
        .map_err(|error| CookingError::new_err(format!("cannot cook: {error:?}")))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
pub fn uncook(
    expression: &PythonExpression,
    settings: Option<&PyCookSettings>,
) -> PythonExpression {
    PyCookSettings::reversible_or(settings)
        .uncook(expression.expr.as_view())
        .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None))]
pub fn schoonschip(
    expression: &PythonExpression,
    settings: Option<&PySchoonschipSettings>,
) -> PythonExpression {
    match settings {
        Some(settings) => expression.expr.schoonschip_with_settings(&settings.rust()),
        None => expression.expr.schoonschip(),
    }
    .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction(signature = (expression, settings = None, *, expand_contracted_sums = false))]
pub fn schoonschip_net(
    expression: &PythonExpression,
    settings: Option<&PySchoonschipSettings>,
    expand_contracted_sums: bool,
) -> PythonExpression {
    let mut settings = settings
        .map(PySchoonschipSettings::rust)
        .unwrap_or_else(RustSchoonschipSettings::default_network);
    settings.expand_contracted_sums |= expand_contracted_sums;
    expression
        .expr
        .schoonschip_with_net::<false, AbstractIndex>(&settings)
        .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn normalize_dots(expression: &PythonExpression) -> PythonExpression {
    expression.expr.normalize_dots().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_dots(expression: &PythonExpression) -> PyResult<PythonExpression> {
    expression
        .expr
        .expand_dots()
        .map(Into::into)
        .map_err(|error| DotExpansionError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn metric_shorthand_to_dot(expression: &PythonExpression) -> PythonExpression {
    expression.expr.metric_shorthand_to_dot().into()
}

macro_rules! undo_shorthand_binding {
    ($name:ident, $method:ident) => {
        #[cfg_attr(
            feature = "python_stubgen",
            gen_stub_pyfunction(module = "symbolica.community.idenso")
        )]
        #[pyfunction]
        pub fn $name(expression: &PythonExpression) -> PythonExpression {
            expression.expr.$method::<AbstractIndex>().into()
        }
    };
}

undo_shorthand_binding!(undo_all, undo_all);
undo_shorthand_binding!(undo_schoonschip, undo_schoonschip);
undo_shorthand_binding!(undo_dots, undo_dots);
undo_shorthand_binding!(undo_chain, undo_chain);
undo_shorthand_binding!(undo_trace, undo_trace);

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn collect_chains(
    expression: &PythonExpression,
    representation: RegisteredRepresentation,
) -> PythonExpression {
    expression.expr.collect_chains(representation.rep()).into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn chainify(
    expression: &PythonExpression,
    representation: RegisteredRepresentation,
) -> PythonExpression {
    expression.expr.chainify(representation.rep()).into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn normalize_chains(expression: &PythonExpression) -> PythonExpression {
    expression.expr.normalize_chains().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn undo_single_length(expression: &PythonExpression) -> PythonExpression {
    expression.expr.undo_single_length().into()
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

    module.add_class::<PyCookMode>()?;
    module.add_class::<PyCookTagFilter>()?;
    module.add_class::<PyCookSourceFilter>()?;
    module.add_class::<PyCookSettings>()?;
    module.add_class::<PySchoonschipMode>()?;
    module.add_class::<PySchoonschipTraversal>()?;
    module.add_class::<PySchoonschipContractionOrder>()?;
    module.add_class::<PySchoonschipSettings>()?;

    module.add_function(wrap_pyfunction!(canonize, module)?)?;
    module.add_function(wrap_pyfunction!(alias_subtensors, module)?)?;
    module.add_function(wrap_pyfunction!(spenso_conjugate, module)?)?;
    module.add_function(wrap_pyfunction!(conjugate_transpose, module)?)?;
    module.add_function(wrap_pyfunction!(cook, module)?)?;
    module.add_function(wrap_pyfunction!(uncook, module)?)?;
    module.add_function(wrap_pyfunction!(schoonschip, module)?)?;
    module.add_function(wrap_pyfunction!(schoonschip_net, module)?)?;
    module.add_function(wrap_pyfunction!(normalize_dots, module)?)?;
    module.add_function(wrap_pyfunction!(expand_dots, module)?)?;
    module.add_function(wrap_pyfunction!(metric_shorthand_to_dot, module)?)?;
    module.add_function(wrap_pyfunction!(undo_all, module)?)?;
    module.add_function(wrap_pyfunction!(undo_schoonschip, module)?)?;
    module.add_function(wrap_pyfunction!(undo_dots, module)?)?;
    module.add_function(wrap_pyfunction!(undo_chain, module)?)?;
    module.add_function(wrap_pyfunction!(undo_trace, module)?)?;
    module.add_function(wrap_pyfunction!(collect_chains, module)?)?;
    module.add_function(wrap_pyfunction!(chainify, module)?)?;
    module.add_function(wrap_pyfunction!(normalize_chains, module)?)?;
    module.add_function(wrap_pyfunction!(undo_single_length, module)?)?;
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
}
