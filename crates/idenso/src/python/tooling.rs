use crate::{
    CookMode as RustCookMode, CookSettings as RustCookSettings,
    CookSourceFilter as RustCookSourceFilter, CookTagFilter as RustCookTagFilter, IndexTooling,
    shorthands::{
        UndoShorthands,
        chain::Chain,
        metric::MetricSimplifier,
        schoonschip::{
            Schoonschip, SchoonschipContractionOrder as RustSchoonschipContractionOrder,
            SchoonschipMode as RustSchoonschipMode, SchoonschipSettings as RustSchoonschipSettings,
            SchoonschipTraversal as RustSchoonschipTraversal,
        },
    },
    tensor::CanonicalizationError as RustCanonicalizationError,
};
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
    representation::{LibraryRep, RepName, Representation},
};
use symbolica::{api::python::PythonExpression, atom::AtomView};

create_exception!(
    symbolica.community.idenso,
    CanonicalizationError,
    PyValueError
);
create_exception!(symbolica.community.idenso, CookingError, PyTypeError);
create_exception!(symbolica.community.idenso, DiracAdjointError, PyValueError);

pub(super) struct RegisteredRepresentation(LibraryRep);

impl<'a, 'py> FromPyObject<'a, 'py> for RegisteredRepresentation {
    type Error = PyErr;

    fn extract(value: Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let expression = if let Ok(expression) = value.extract::<PythonExpression>() {
            expression
        } else {
            value
                .call_method0("to_expression")?
                .extract::<PythonExpression>()
                .map_err(PyErr::from)?
        };
        let view = expression.expr.as_view();
        if let Ok(representation) = Representation::<LibraryRep>::try_from(view) {
            return Ok(Self(representation.rep));
        }
        if let AtomView::Var(variable) = view {
            let symbol = variable.get_symbol();
            return LibraryRep::try_from_symbol_coerced(symbol)
                .map(Self)
                .map_err(|error| {
                    PyValueError::new_err(format!(
                        "{symbol} is not a registered Spenso representation: {error}"
                    ))
                });
        }
        Err(PyTypeError::new_err(
            "representation must be a Spenso Representation or its symbolic expression",
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
    name = "CookTagFilter",
    from_py_object,
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PyCookTagFilter {
    inner: RustCookTagFilter,
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
    name = "CookSourceFilter",
    from_py_object,
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PyCookSourceFilter {
    inner: RustCookSourceFilter,
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
            inner: RustCookSourceFilter::FunctionTags(filter.inner.clone()),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (filter=None))]
    pub fn representation_index_payload(filter: Option<&PyCookTagFilter>) -> Self {
        Self {
            inner: RustCookSourceFilter::RepresentationIndexPayload {
                filter: filter.map(|filter| filter.inner.clone()),
            },
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CookSettings",
    from_py_object,
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PyCookSettings {
    inner: RustCookSettings,
}

impl PyCookSettings {
    pub(super) fn indices_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(|settings| settings.inner.clone())
            .unwrap_or_else(RustCookSettings::indices)
    }

    pub(super) fn flattened_or(settings: Option<&Self>) -> RustCookSettings {
        settings
            .map(|settings| settings.inner.clone())
            .unwrap_or_else(RustCookSettings::flattened)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCookSettings {
    #[new]
    pub fn new() -> Self {
        Self::flattened()
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

    pub fn with_mode(&self, mode: PyCookMode) -> Self {
        Self {
            inner: self.inner.clone().with_mode(mode.into()),
        }
    }

    pub fn with_output_tags(&self, tags: Vec<String>) -> Self {
        Self {
            inner: self.inner.clone().with_output_tags(tags),
        }
    }

    pub fn with_source_filter(&self, source: &PyCookSourceFilter) -> Self {
        Self {
            inner: self.inner.clone().with_source_filter(source.inner.clone()),
        }
    }

    pub fn with_input_tags(&self, tags: Vec<String>) -> Self {
        Self {
            inner: self.inner.clone().with_input_tags(tags),
        }
    }

    pub fn with_all_input_tags(&self, tags: Vec<String>) -> Self {
        Self {
            inner: self.inner.clone().with_all_input_tags(tags),
        }
    }

    pub fn matched_filter(&self) -> Self {
        Self {
            inner: self.inner.clone().matched_filter(),
        }
    }

    pub fn preserve_tags(&self) -> Self {
        Self {
            inner: self.inner.clone().preserve_tags(),
        }
    }

    #[pyo3(signature = (filter=None))]
    pub fn with_index_payload_filter(&self, filter: Option<&PyCookTagFilter>) -> Self {
        Self {
            inner: self
                .inner
                .clone()
                .with_index_payload_filter(filter.map(|filter| filter.inner.clone())),
        }
    }

    pub fn mode(&self) -> PyCookMode {
        self.inner.mode().into()
    }

    pub fn source_filter(&self) -> PyCookSourceFilter {
        PyCookSourceFilter {
            inner: self.inner.source_filter().clone(),
        }
    }

    pub fn output_tags(&self) -> Vec<String> {
        self.inner.output_tags().to_vec()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
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

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "SchoonschipSettings",
    from_py_object,
    module = "symbolica.community.idenso"
)]
#[derive(Clone)]
pub struct PySchoonschipSettings {
    inner: RustSchoonschipSettings,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PySchoonschipSettings {
    #[new]
    pub fn new() -> Self {
        Self {
            inner: RustSchoonschipSettings::default(),
        }
    }

    #[staticmethod]
    pub fn partial() -> Self {
        Self {
            inner: RustSchoonschipSettings::partial(),
        }
    }

    #[staticmethod]
    pub fn with_depth(depth_limit: usize) -> Self {
        Self {
            inner: RustSchoonschipSettings::with_depth(depth_limit),
        }
    }

    #[staticmethod]
    pub fn breadth_first_with_depth(depth_limit: usize) -> Self {
        Self {
            inner: RustSchoonschipSettings::breadth_first_with_depth(depth_limit),
        }
    }

    #[staticmethod]
    pub fn full() -> Self {
        Self {
            inner: RustSchoonschipSettings::full(),
        }
    }

    #[staticmethod]
    pub fn default_network() -> Self {
        Self {
            inner: RustSchoonschipSettings::default_network(),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (depth_limit=None))]
    pub fn depth_first(depth_limit: Option<usize>) -> Self {
        Self {
            inner: RustSchoonschipSettings::depth_first(depth_limit),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (depth_limit=None))]
    pub fn breadth_first(depth_limit: Option<usize>) -> Self {
        Self {
            inner: RustSchoonschipSettings::breadth_first(depth_limit),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (depth_limit=None))]
    pub fn single_pass(depth_limit: Option<usize>) -> Self {
        Self {
            inner: RustSchoonschipSettings::single_pass(depth_limit),
        }
    }

    pub fn with_expanded_contracted_sums(&self) -> Self {
        Self {
            inner: self.inner.with_expanded_contracted_sums(),
        }
    }

    pub fn with_chain_like_functions(&self) -> Self {
        Self {
            inner: self.inner.with_chain_like_functions(),
        }
    }

    pub fn without_chain_like_functions(&self) -> Self {
        Self {
            inner: self.inner.without_chain_like_functions(),
        }
    }

    pub fn with_rank1_tensors(&self) -> Self {
        Self {
            inner: self.inner.with_rank1_tensors(),
        }
    }

    pub fn without_rank1_tensors(&self) -> Self {
        Self {
            inner: self.inner.without_rank1_tensors(),
        }
    }

    pub fn with_contraction_order(&self, order: PySchoonschipContractionOrder) -> Self {
        Self {
            inner: self.inner.with_contraction_order(order.into()),
        }
    }

    pub fn with_single_pass(&self) -> Self {
        Self {
            inner: self.inner.into_single_pass(),
        }
    }

    pub fn depth_limit(&self) -> Option<usize> {
        self.inner.depth_limit
    }

    pub fn mode(&self) -> PySchoonschipMode {
        match self.inner.mode {
            RustSchoonschipMode::SinglePass => PySchoonschipMode::SinglePass,
            RustSchoonschipMode::Recursive(_) => PySchoonschipMode::Recursive,
        }
    }

    pub fn traversal(&self) -> Option<PySchoonschipTraversal> {
        match self.inner.mode {
            RustSchoonschipMode::SinglePass => None,
            RustSchoonschipMode::Recursive(RustSchoonschipTraversal::DepthFirst) => {
                Some(PySchoonschipTraversal::DepthFirst)
            }
            RustSchoonschipMode::Recursive(RustSchoonschipTraversal::BreadthFirst) => {
                Some(PySchoonschipTraversal::BreadthFirst)
            }
        }
    }

    pub fn expand_contracted_sums(&self) -> bool {
        self.inner.expand_contracted_sums
    }

    pub fn simplify_chain_like_functions(&self) -> bool {
        self.inner.simplify_chain_like_functions
    }

    pub fn schoonschip_rank1_tensors(&self) -> bool {
        self.inner.schoonschip_rank1_tensors
    }

    pub fn contraction_order(&self) -> PySchoonschipContractionOrder {
        self.inner.contraction_order.into()
    }
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn canonize(self_: &PythonExpression) -> PyResult<PythonExpression> {
    self_
        .expr
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .map(Into::into)
        .map_err(|error: RustCanonicalizationError| {
            CanonicalizationError::new_err(error.to_string())
        })
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
/// Replace repeated non-tensor subexpressions with indexed tensor aliases.
///
/// Returns the aliased root expression and a deterministic list of
/// `(alias, original)` definitions. Definitions may contain other aliases.
pub fn alias_subtensors(
    self_: &PythonExpression,
    tensor_name: &str,
) -> (PythonExpression, Vec<(PythonExpression, PythonExpression)>) {
    let (root, aliases) = self_
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
pub fn spenso_conjugate(self_: &PythonExpression) -> PythonExpression {
    self_.expr.spenso_conj().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn conjugate_transpose(
    self_: &PythonExpression,
    representation: RegisteredRepresentation,
) -> PyResult<PythonExpression> {
    Ok(self_.expr.conjugate_transpose(representation.0).into())
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn dirac_adjoint(self_: &PythonExpression) -> PyResult<PythonExpression> {
    self_
        .expr
        .dirac_adjoint::<AbstractIndex>()
        .map(Into::into)
        .map_err(|error| DiracAdjointError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
#[pyo3(signature = (self_, settings=None))]
pub fn cook(
    self_: &PythonExpression,
    settings: Option<&PyCookSettings>,
) -> PyResult<PythonExpression> {
    let default_settings = RustCookSettings::reversible();
    let settings = settings.map_or(&default_settings, |settings| &settings.inner);
    settings
        .try_cook(self_.expr.as_view())
        .map(Into::into)
        .map_err(|error| CookingError::new_err(format!("cannot cook: {error:?}")))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
#[pyo3(signature = (self_, settings=None))]
pub fn uncook(self_: &PythonExpression, settings: Option<&PyCookSettings>) -> PythonExpression {
    let default_settings = RustCookSettings::reversible();
    let settings = settings.map_or(&default_settings, |settings| &settings.inner);
    settings.uncook(self_.expr.as_view()).into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
#[pyo3(signature = (self_, settings=None))]
pub fn schoonschip(
    self_: &PythonExpression,
    settings: Option<&PySchoonschipSettings>,
) -> PythonExpression {
    match settings {
        Some(settings) => self_.expr.schoonschip_with_settings(&settings.inner),
        None => self_.expr.schoonschip(),
    }
    .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
#[pyo3(signature = (self_, settings=None, expand_contracted_sums=false))]
pub fn schoonschip_net(
    self_: &PythonExpression,
    settings: Option<&PySchoonschipSettings>,
    expand_contracted_sums: bool,
) -> PythonExpression {
    let mut settings = settings
        .map(|settings| settings.inner)
        .unwrap_or_else(RustSchoonschipSettings::default_network);
    settings.expand_contracted_sums |= expand_contracted_sums;
    self_
        .expr
        .schoonschip_with_net::<false, AbstractIndex>(&settings)
        .into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn normalize_dots(self_: &PythonExpression) -> PythonExpression {
    self_.expr.normalize_dots().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn expand_dots(self_: &PythonExpression) -> PyResult<PythonExpression> {
    self_
        .expr
        .expand_dots()
        .map(Into::into)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn metric_shorthand_to_dot(self_: &PythonExpression) -> PythonExpression {
    self_.expr.metric_shorthand_to_dot().into()
}

macro_rules! undo_shorthand_binding {
    ($name:ident, $method:ident) => {
        #[cfg_attr(
            feature = "python_stubgen",
            gen_stub_pyfunction(module = "symbolica.community.idenso")
        )]
        #[pyfunction]
        pub fn $name(self_: &PythonExpression) -> PythonExpression {
            self_.expr.$method::<AbstractIndex>().into()
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
    self_: &PythonExpression,
    representation: RegisteredRepresentation,
) -> PythonExpression {
    self_.expr.collect_chains(representation.0).into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn chainify(
    self_: &PythonExpression,
    representation: RegisteredRepresentation,
) -> PythonExpression {
    self_.expr.chainify(representation.0).into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn normalize_chains(self_: &PythonExpression) -> PythonExpression {
    self_.expr.normalize_chains().into()
}

#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.idenso")
)]
#[pyfunction]
pub fn undo_single_length(self_: &PythonExpression) -> PythonExpression {
    self_.expr.undo_single_length().into()
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "CanonicalizationError",
        m.py().get_type::<CanonicalizationError>(),
    )?;
    m.add("CookingError", m.py().get_type::<CookingError>())?;
    m.add("DiracAdjointError", m.py().get_type::<DiracAdjointError>())?;

    m.add_class::<PyCookMode>()?;
    m.add_class::<PyCookTagFilter>()?;
    m.add_class::<PyCookSourceFilter>()?;
    m.add_class::<PyCookSettings>()?;
    m.add_class::<PySchoonschipMode>()?;
    m.add_class::<PySchoonschipTraversal>()?;
    m.add_class::<PySchoonschipContractionOrder>()?;
    m.add_class::<PySchoonschipSettings>()?;

    m.add_function(wrap_pyfunction!(canonize, m)?)?;
    m.add_function(wrap_pyfunction!(alias_subtensors, m)?)?;
    m.add_function(wrap_pyfunction!(spenso_conjugate, m)?)?;
    m.add_function(wrap_pyfunction!(conjugate_transpose, m)?)?;
    m.add_function(wrap_pyfunction!(dirac_adjoint, m)?)?;
    m.add_function(wrap_pyfunction!(cook, m)?)?;
    m.add_function(wrap_pyfunction!(uncook, m)?)?;
    m.add_function(wrap_pyfunction!(schoonschip, m)?)?;
    m.add_function(wrap_pyfunction!(schoonschip_net, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_dots, m)?)?;
    m.add_function(wrap_pyfunction!(expand_dots, m)?)?;
    m.add_function(wrap_pyfunction!(metric_shorthand_to_dot, m)?)?;
    m.add_function(wrap_pyfunction!(undo_all, m)?)?;
    m.add_function(wrap_pyfunction!(undo_schoonschip, m)?)?;
    m.add_function(wrap_pyfunction!(undo_dots, m)?)?;
    m.add_function(wrap_pyfunction!(undo_chain, m)?)?;
    m.add_function(wrap_pyfunction!(undo_trace, m)?)?;
    m.add_function(wrap_pyfunction!(collect_chains, m)?)?;
    m.add_function(wrap_pyfunction!(chainify, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_chains, m)?)?;
    m.add_function(wrap_pyfunction!(undo_single_length, m)?)?;
    Ok(())
}
