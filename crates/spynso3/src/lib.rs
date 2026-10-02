#[macro_use]
mod documentation;

use std::{cell::Cell, collections::HashMap, ops::Deref};

use eyre::eyre;

use broadcast::SpensoBroadcastFunction;
use library::{SpensorFunctionLibrary, SpensorLibrary};
use network::{ConvertibleToSpensoNet, SpensoNet};

use pyo3::{
    PyClass,
    exceptions::{self, PyIndexError, PyOverflowError, PyRuntimeError, PyTypeError},
    prelude::*,
    types::{PyComplex, PyDict, PyFloat, PyInt, PySlice, PyTuple, PyType},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};

#[cfg(feature = "native")]
use spenso::{
    algebra::complex::symbolica_traits::CompiledComplexEvaluatorSpenso,
    tensors::parametric::EvalTensor,
};

use spenso::{
    algebra::complex::{Complex, RealOrComplex},
    tensors::{
        data::{DenseTensor, GetTensorData, SetTensorData, SparseOrDense, SparseTensor},
        parametric::{ConcreteOrParam, ParamOrConcrete, ParamTensor, atomcore::TensorAtomOps},
    },
};

use spenso::{
    network::parsing::ShadowedStructure,
    structure::{
        Canonicalized, HasStructure, OrderedStructure, ScalarTensor, TensorDataLayout,
        TensorDataLayoutError, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialStructure, PartialStructureExt},
        slot::IsAbstractSlot,
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, StorageTensor},
        parametric::{LinearizedEvalTensor, MixedTensor},
    },
};
use symbolica::{
    api::python::{Citation, SymbolicaCommunityModule},
    domains::float::Complex as SymComplex,
    prelude::*,
};

use symbolica::api::python::{ConvertibleToExpression, PythonExpression, PythonFormattedOutput};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, TypeInfo, derive::*, impl_stub_type};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass_enum, gen_stub_pyfunction};

pub mod broadcast;
mod data;
pub mod display;
pub mod expression;
pub mod library;
pub mod metadata;
pub mod network;
pub mod pattern;
pub mod projectors;
mod simplification;
pub mod structure;
pub mod tensor_rule;

use expression::TensorExpression;
use idenso::tensor::SymbolicTensor;
use structure::ConvertibleToSpensoName;

trait ModuleInit: PyClass {
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<Self>()
    }
}

#[cfg(feature = "python_stubgen")]
/// The Python names that must be present in Spenso's generated stub.
///
/// `registered` is the actual community-module surface. `returned_opaque`
/// contains classes that are not module attributes but are returned by that
/// surface and therefore need public type declarations.
pub struct PythonStubSurface {
    pub registered: &'static [&'static str],
    pub returned_opaque: &'static [&'static str],
}

macro_rules! define_spenso_python_surface {
    (
        registered_classes: [$($(#[$class_attr:meta])* $class:ty),+ $(,)?],
        registered_functions: [$($function:ident => $function_name:literal),+ $(,)?],
        registered_modules: [$($module:ident => [$($export:literal),* $(,)?]),* $(,)?],
        returned_opaque_classes: [$($opaque:ty),* $(,)?],
    ) => {
        pub(crate) fn initialize_spenso(m: &Bound<'_, PyModule>) -> PyResult<()> {
            m.add("__doc__", python_doc!("module"))?;
            $($(#[$class_attr])* m.add_class::<$class>()?;)+
            // m.add_function(?)?;
            $(m.add_function(wrap_pyfunction!($function, m)?)?;)+
            $($module::register(m)?;)*
            let exports = m
                .dict()
                .keys()
                .iter()
                .filter_map(|key| key.extract::<String>().ok())
                .filter(|name| {
                    name != "initialize"
                        && (name == "_" || name == "initialize_module" || !name.starts_with('_'))
                })
                .collect::<Vec<_>>();
            m.add("__all__", exports)?;
            Ok(())
        }

        #[cfg(feature = "python_stubgen")]
        pub const PYTHON_STUB_SURFACE: PythonStubSurface = PythonStubSurface {
            registered: &[
                $($(#[$class_attr])* <$class as PyClass>::NAME,)+
                $($function_name,)+
                $($($export,)*)*
            ],
            returned_opaque: &[$(<$opaque as PyClass>::NAME,)*],
        };
    };
}

pub struct SpensoModule;

#[cfg(feature = "python_stubgen")]
mod stubs;

/// Threading policy for operations on symbolic component expressions.
///
/// ``Auto`` permits parallel work when the Symbolica license supports it,
/// with workload-size heuristics. ``Serial`` keeps symbolic work on the calling
/// thread. ``Parallel`` requests parallel work without Auto's license check.
/// This does not set the number of threads or force every operation to parallelize.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import SymbolicParallelism, set_symbolica_rayon_enabled
/// >>> allowed = set_symbolica_rayon_enabled(SymbolicParallelism.Auto)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(from_py_object, eq, eq_int, module = "symbolica.community.tensor")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymbolicParallelism {
    /// Permit Rayon when licensed and use workload heuristics where available.
    Auto,
    /// Keep symbolic operations on the calling thread.
    Serial,
    /// Force Rayon without `Auto`'s Symbolica license safety check.
    Parallel,
}

impl ModuleInit for SymbolicParallelism {}

impl From<SymbolicParallelism> for spenso::symbolic_parallelism::SymbolicParallelism {
    fn from(value: SymbolicParallelism) -> Self {
        match value {
            SymbolicParallelism::Auto => Self::Auto,
            SymbolicParallelism::Serial => Self::Serial,
            SymbolicParallelism::Parallel => Self::Parallel,
        }
    }
}

/// Set the threading policy for Spenso's Symbolica operations.
///
/// Parameters
/// ----------
/// policy : SymbolicParallelism
///     Auto checks the license at this call; Serial disables symbolic
///     parallelism; Parallel enables it without that automatic check.
///
/// Returns
/// -------
/// bool
///     Whether the resolved policy permits parallel execution. Individual
///     operations may still choose serial execution for small workloads.
///
/// Notes
/// -----
/// This changes the process-wide policy, not just one tensor or network.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import SymbolicParallelism, set_symbolica_rayon_enabled
/// >>> allowed = set_symbolica_rayon_enabled(SymbolicParallelism.Auto)
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.tensor")
)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pyfunction]
fn set_symbolica_rayon_enabled(policy: SymbolicParallelism) -> bool {
    spenso::symbolic_parallelism::set_symbolica_rayon_enabled(policy.into());
    spenso::symbolic_parallelism::symbolica_rayon_enabled()
}

impl SymbolicaCommunityModule for SpensoModule {
    fn get_citations() -> Vec<Citation> {
        if !CITATIONS_USED.load(std::sync::atomic::Ordering::Relaxed) {
            return Vec::new();
        }
        // Authorship and software DOIs: docs/products/registry.toml.
        vec![
            Citation {
                id: "10.5281/zenodo.18248388".into(),
                reference:
                    "Lucien Huber, Valentin Hirschi, Ben Ruijl, Mathijs Fraaije. Spenso (2026)."
                        .into(),
                bibtex: r#"@software{spenso,
  author = {Lucien Huber and Valentin Hirschi and Ben Ruijl and Mathijs Fraaije},
  title = {Spenso},
  year = {2026},
  url = {https://github.com/alphal00p/spenso},
  doi = {10.5281/zenodo.18248388}
}"#
                .into(),
                reasons: vec!["Provides the Spenso functionality in this community module.".into()],
                description: "".into(),
                relevance: None,
            },
            Citation {
                id: "10.5281/zenodo.18248409".into(),
                reference:
                    "Lucien Huber, Valentin Hirschi, Ben Ruijl, Mathijs Fraaije. Idenso (2026)."
                        .into(),
                bibtex: r#"@software{idenso,
  author = {Lucien Huber and Valentin Hirschi and Ben Ruijl and Mathijs Fraaije},
  title = {Idenso},
  year = {2026},
  url = {https://github.com/alphal00p/spenso},
  doi = {10.5281/zenodo.18248409}
}"#
                .into(),
                reasons: vec!["Provides the Idenso functionality in this community module.".into()],
                description: "".into(),
                relevance: None,
            },
        ]
    }

    fn get_name() -> String {
        "tensor".to_string()
    }

    fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
        initialize_spenso(m)
    }

    fn initialize(_py: Python) -> PyResult<()> {
        idenso::representations::initialize();
        spenso::symbolic_parallelism::set_symbolica_rayon_enabled(
            spenso::symbolic_parallelism::SymbolicParallelism::Auto,
        );
        Ok(())
    }
}

define_spenso_python_surface! {
    registered_classes: [
        tensor_rule::PyTensorRule,
        SpensoNet,
        network::ExecutionMode,
        network::execution::ExecutionStatus,
        SymbolicParallelism,
        SpensoExpressionEvaluator,
        #[cfg(feature = "native")]
        SpensoCompiledExpressionEvaluator,
        Spensor,
        structure::SpensoName,
        structure::SpensoSlot,
        structure::SpensoRepresentation,
        metadata::SpensoRepresentationName,
        metadata::SpensoTensorStructure,
        SpensorLibrary,
        SpensorFunctionLibrary,
        SpensoBroadcastFunction,
        projectors::PyFactorProjector,
    ],
    registered_functions: [
        set_symbolica_rayon_enabled => "set_symbolica_rayon_enabled",
    ],
    registered_modules: [
        simplification => [
            "CanonicalizationError", "CookingError",
            "DiracAdjointError", "NetworkToolingError",
            "ReductionStatus",
        ],
        display => [
            "DisplaySettings", "format_tensor", "to_typst", "to_html", "to_svg", "formatted",
            "load_math_font",
        ],
        pattern => ["PortPattern", "TensorPattern"],
        expression => [
            "TensorExpression", "_AutoIndex", "AUTO", "_", "Nc",
            "as_tensor", "dot", "chain", "trace",
        ],
    ],
    returned_opaque_classes: [],
}

/// Tensor components together with their symbolic identity and ordered axes.
///
/// Use dense(), sparse(), or from_numpy() to construct one. Components can be
/// real floats, complex floats, or Symbolica expressions. Square brackets access
/// components; calling the tensor assigns abstract index labels and returns a
/// lazy TensorNetwork. Arithmetic also builds networks, retaining component data.
///
/// The default notebook view is an interactive component explorer with a matrix
/// toggle. ``to_typst()`` provides static mathematical source.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import Tensor
/// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
/// >>> tensor[1, 0]
/// 3.0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object, name = "Tensor", module = "symbolica.community.tensor")]
#[derive(Clone)]
pub struct Spensor {
    pub(crate) tensor: MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
    pub(crate) descriptor: SymbolicTensor<PartialStructure>,
    pub(crate) descriptor_name: Option<Symbol>,
    pub(crate) descriptor_args: Vec<Atom>,
}

pub struct TensorDataDescriptor {
    structure: Canonicalized<ShadowedStructure<AbstractIndex>>,
    descriptor: SymbolicTensor<PartialStructure>,
    name: Symbol,
    args: Vec<Atom>,
}

pub enum AtomsOrFloats {
    Atoms(Vec<Atom>),
    Floats(Vec<f64>),
    Complex(Vec<Complex<f64>>),
}

impl<'a, 'py> FromPyObject<'a, 'py> for AtomsOrFloats {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        // Check typed interfaces before inherited numeric coercion can erase them,
        // particularly for tensor-valued zeros.
        let components = value.extract::<Vec<Bound<'py, PyAny>>>()?;
        for component in &components {
            TensorElements::validate_scalar(component)?;
        }
        // Explicit symbolic data must not go through Expression.__float__.
        // Integer-only input is exact too; explicitly floating-point input
        // continues to select numerical storage.
        let exact = components
            .iter()
            .any(|value| value.is_instance_of::<PythonExpression>())
            || (!components.is_empty()
                && components
                    .iter()
                    .all(|value| value.is_instance_of::<PyInt>()));
        if !exact {
            if let Ok(values) = value.extract::<Vec<f64>>() {
                return Ok(Self::Floats(values));
            }
            if let Ok(values) = value.extract::<Vec<Complex<f64>>>() {
                return Ok(Self::Complex(values));
            }
        }
        let values = value
            .extract::<Vec<ConvertibleToExpression>>()
            .map_err(|_| {
                PyTypeError::new_err(
                    "data must be a list of integers, floats, complex numbers, or Expressions",
                )
            })?;
        Ok(Self::Atoms(
            values
                .into_iter()
                .map(|value| value.to_expression().expr)
                .collect(),
        ))
    }
}

#[cfg(feature = "python_stubgen")]
impl_stub_type!(AtomsOrFloats = Vec<TensorElements>);

fn tensor_data_layout_error(error: TensorDataLayoutError) -> PyErr {
    match error {
        TensorDataLayoutError::Overflow => PyOverflowError::new_err(error.to_string()),
        TensorDataLayoutError::RankMismatch { .. }
        | TensorDataLayoutError::IndexOutOfBounds { .. }
        | TensorDataLayoutError::FlatIndexOutOfBounds { .. } => {
            PyIndexError::new_err(error.to_string())
        }
        TensorDataLayoutError::NonConcreteDimension | TensorDataLayoutError::DataLength { .. } => {
            exceptions::PyValueError::new_err(error.to_string())
        }
    }
}

/// Build the mapping from the public logical row-major layout to canonical storage.
fn tensor_data_layout(interface: &PartialStructure) -> PyResult<TensorDataLayout> {
    TensorDataLayout::from_canonicalized(interface).map_err(tensor_data_layout_error)
}

impl<'a, 'py> FromPyObject<'a, 'py> for TensorDataDescriptor {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let expression = value
            .extract::<PyRef<'_, TensorExpression>>()
            .map_err(|_| {
                PyTypeError::new_err("tensor data structure must be a TensorExpression")
            })?;
        let descriptor = TensorExpression::structured(&expression);
        let name = TensorExpression::descriptor_name(&expression).ok_or_else(|| {
            exceptions::PyValueError::new_err(
                "tensor data descriptor has no name; use expression.with_name(TensorName(...))",
            )
        })?;
        let args = TensorExpression::descriptor_args(&expression);
        Self::new(descriptor.clone(), name, args)
    }
}

impl TensorDataDescriptor {
    fn new(
        descriptor: SymbolicTensor<PartialStructure>,
        name: Symbol,
        args: Vec<Atom>,
    ) -> PyResult<Self> {
        let layout = tensor_data_layout(descriptor.structure())?;
        let owner = AbstractIndex::fresh_open_owner();
        let logical_axes = (0..layout.logical_shape().len()).collect::<Vec<_>>();
        let storage_axes = descriptor
            .structure()
            .layout()
            .logical_to_canonical(&logical_axes);
        let mut storage_indices = vec![0; storage_axes.len()];
        for (storage_index, &logical_axis) in storage_axes.iter().enumerate() {
            storage_indices[logical_axis] = storage_index;
        }
        let structure = OrderedStructure::new(
            descriptor
                .structure()
                .logical_slots()
                .into_iter()
                .zip(storage_indices)
                .map(|(slot, axis)| slot.rep().slot(AbstractIndex::Open { owner, axis }))
                .collect(),
        )
        .map_canonical(|structure| ShadowedStructure {
            structure,
            global_name: Some(name),
            additional_args: (!args.is_empty()).then_some(args.clone()),
        });
        Ok(Self {
            structure,
            descriptor,
            name,
            args,
        })
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for TensorDataDescriptor {
    fn type_output() -> TypeInfo {
        TensorExpression::type_output()
    }
}

impl Spensor {
    pub(crate) fn from_storage_with_descriptor(
        tensor: MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
        descriptor: SymbolicTensor<PartialStructure>,
        descriptor_name: Option<Symbol>,
        descriptor_args: Vec<Atom>,
    ) -> Self {
        Self {
            tensor,
            descriptor,
            descriptor_name,
            descriptor_args,
        }
    }

    pub(crate) fn scalar_with_descriptor(
        value: ConcreteOrParam<RealOrComplex<f64>>,
        descriptor: SymbolicTensor<PartialStructure>,
        descriptor_name: Option<Symbol>,
        descriptor_args: Vec<Atom>,
    ) -> Self {
        Self::from_storage_with_descriptor(
            ParamOrConcrete::new_scalar(value),
            descriptor,
            descriptor_name,
            descriptor_args,
        )
    }
}

impl Deref for Spensor {
    type Target = MixedTensor<f64, ShadowedStructure<AbstractIndex>>;

    fn deref(&self) -> &Self::Target {
        &self.tensor
    }
}

impl ModuleInit for Spensor {}

// #[gen_stub_pyclass_enum]

#[derive(FromPyObject)]
pub enum SliceOrIntOrExpanded<'a> {
    Slice(Bound<'a, PySlice>),
    Int(usize),
    Expanded(Vec<usize>),
    Selection(Bound<'a, PyAny>),
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for SliceOrIntOrExpanded<'_> {
    fn type_input() -> pyo3_stub_gen::TypeInfo {
        TypeInfo::builtin("slice") | usize::type_input() | TypeInfo::list_of::<usize>()
    }

    fn type_output() -> pyo3_stub_gen::TypeInfo {
        TypeInfo::builtin("slice") | usize::type_input() | TypeInfo::list_of::<usize>()
    }
}

#[derive(IntoPyObject)]
pub enum TensorElements {
    Real(Py<PyFloat>),
    Complex(Py<PyComplex>),
    Symbolica(PythonExpression),
}

impl TensorElements {
    fn validate_scalar(value: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Ok(tensor) = value.extract::<PyRef<'_, TensorExpression>>()
            && !tensor.interface().canonical().is_scalar()
        {
            return Err(PyTypeError::new_err("tensor components must be scalar"));
        }
        Ok(())
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for TensorElements {
    fn type_input() -> pyo3_stub_gen::TypeInfo {
        PythonExpression::type_input()
            | PyInt::type_input()
            | Complex::type_input()
            | PyFloat::type_input()
    }

    fn type_output() -> TypeInfo {
        PythonExpression::type_output() | Complex::type_output() | PyFloat::type_output()
    }
}

impl From<ConcreteOrParam<RealOrComplex<f64>>> for TensorElements {
    fn from(value: ConcreteOrParam<RealOrComplex<f64>>) -> Self {
        match value {
            ConcreteOrParam::Concrete(RealOrComplex::Real(f)) => {
                TensorElements::Real(Python::attach(|py| {
                    PyFloat::new(py, f).as_unbound().to_owned()
                }))
            }
            ConcreteOrParam::Concrete(RealOrComplex::Complex(c)) => {
                TensorElements::Complex(Python::attach(|py| {
                    PyComplex::from_doubles(py, c.re, c.im)
                        .as_unbound()
                        .to_owned()
                }))
            }
            ConcreteOrParam::Param(p) => TensorElements::Symbolica(PythonExpression::from(p)),
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl Spensor {
    /// Return the symbolic descriptor for these components.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A symbolic tensor with the same external axes.
    ///
    /// Notes
    /// -----
    /// The descriptor identifies the data but does not embed its components.
    /// Register the tensor in a TensorLibrary to evaluate that symbolic reference.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.expression().rank
    /// 2
    pub fn expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_known_parts(
            py,
            self.descriptor.expression().clone(),
            self.descriptor.structure().clone(),
            self.descriptor_name,
            self.descriptor_args.clone(),
        )
    }

    #[doc = python_doc!("Tensor.axes")]
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[Slot | Representation, ...]"))]
    fn axes(&self, py: Python<'_>) -> PyResult<Py<PyTuple>> {
        crate::metadata::axis_objects(py, self.descriptor.structure().logical_slots())
    }

    /// Canonical external signature of the whole tensor.
    ///
    /// Returns
    /// -------
    /// TensorStructure
    ///     Free axes in canonical order, independent of names and component layout.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.structure.rank
    /// 2
    #[getter]
    fn structure(&self) -> metadata::SpensoTensorStructure {
        crate::metadata::SpensoTensorStructure::from_interface(self.descriptor.structure())
    }

    /// Assign a data identity without rewriting the symbolic computation.
    ///
    /// Parameters
    /// ----------
    /// name : TensorName, str, Expression, or TensorExpression
    ///     New tensor name. A rank-zero atomic TensorExpression may also supply
    ///     scalar key arguments, for example TensorName("B")(7).
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A new value with the requested identity and the same axes and components.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> named = tensor.with_name("named_matrix")
    /// >>> named.rank
    /// 2
    fn with_name(&self, name: ConvertibleToSpensoName) -> Self {
        let ConvertibleToSpensoName(name, args) = name;
        let name = name.name;
        let storage = self
            .tensor
            .clone()
            .map_structure(|structure| ShadowedStructure {
                structure: structure.structure,
                global_name: Some(name),
                additional_args: (!args.is_empty()).then_some(args.clone()),
            });
        Self::from_storage_with_descriptor(storage, self.descriptor.clone(), Some(name), args)
    }

    /// Create an initially zero tensor using sparse component storage.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorExpression
    ///     Named descriptor with concrete dimensions in logical axis order.
    /// type_info : type
    ///     Component type: float or Expression. Expression storage preserves
    ///     exact symbolic values and accepts integer assignments without rounding.
    ///     Floating-point storage requires numerical assignments.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A zero tensor storing only explicitly populated entries.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.sparse(A, float)
    /// >>> tensor[1, 0] = 3.0
    /// >>> tensor[0, 0]
    /// 0.0
    #[staticmethod]
    pub fn sparse(
        structure: TensorDataDescriptor,
        type_info: Bound<'_, PyType>,
    ) -> PyResult<Spensor> {
        let TensorDataDescriptor {
            structure,
            descriptor,
            name,
            args,
        } = structure;
        if type_info.is_subclass_of::<PyFloat>()? {
            Ok(Spensor::from_storage_with_descriptor(
                structure
                    .map_canonical(|s| SparseTensor::<f64, _>::empty(s, 0.0).into())
                    .into_canonical(),
                descriptor,
                Some(name),
                args,
            ))
        } else if type_info.is_subclass_of::<PythonExpression>()? {
            Ok(Spensor::from_storage_with_descriptor(
                structure
                    .map_canonical(|s| {
                        ParamOrConcrete::Param(ParamTensor::from(SparseTensor::<Atom, _>::empty(
                            s,
                            Atom::Zero,
                        )))
                    })
                    .into_canonical(),
                descriptor,
                Some(name),
                args,
            ))
        } else {
            Err(PyTypeError::new_err("Only float type supported"))
        }
    }

    /// Store every component of a named tensor.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorExpression
    ///     Named descriptor with concrete dimensions and the desired logical
    ///     axis order. Use with_name() to give a composite expression a name.
    /// data : sequence of int, float, complex, or Expression
    ///     Flat component data in logical row-major order. The length must equal
    ///     the product of the dimensions. An explicit Expression selects symbolic
    ///     storage for the sequence, preserving exact rational and algebraic
    ///     values. Integer-only input is exact too. Float or complex sequences
    ///     without Expressions retain numerical storage.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A dense component tensor; the input descriptor remains symbolic.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.shape
    /// (2, 2)
    #[staticmethod]
    pub fn dense(structure: TensorDataDescriptor, data: AtomsOrFloats) -> PyResult<Spensor> {
        let TensorDataDescriptor {
            structure,
            descriptor,
            name,
            args,
        } = structure;
        let layout = tensor_data_layout(descriptor.structure())?;
        let storage_structure = structure.into_canonical();
        let dense = match data {
            AtomsOrFloats::Floats(f) => DenseTensor::<f64, _>::from_storage_data(
                layout
                    .reorder_to_storage(f)
                    .map_err(tensor_data_layout_error)?,
                storage_structure,
            )
            .map_err(|e| PyOverflowError::new_err(e.to_string()))?
            .into(),
            AtomsOrFloats::Atoms(a) => ParamOrConcrete::Param(ParamTensor::from(
                DenseTensor::<Atom, _>::from_storage_data(
                    layout
                        .reorder_to_storage(a)
                        .map_err(tensor_data_layout_error)?,
                    storage_structure,
                )
                .map_err(|e| PyOverflowError::new_err(e.to_string()))?,
            )),
            AtomsOrFloats::Complex(c) => {
                MixedTensor::Concrete(RealOrComplexTensor::Complex(DataTensor::Dense(
                    DenseTensor::<Complex<f64>, _>::from_storage_data(
                        layout
                            .reorder_to_storage(c)
                            .map_err(tensor_data_layout_error)?,
                        storage_structure,
                    )
                    .map_err(|e| PyOverflowError::new_err(e.to_string()))?,
                )))
            }
        };

        Ok(Spensor::from_storage_with_descriptor(
            dense,
            descriptor,
            Some(name),
            args,
        ))
    }
    /// Copy the components into dense storage.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     Independent dense tensor with the same values and logical axes.
    ///
    /// Notes
    /// -----
    /// Dense storage allocates space for every component.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.to_dense().storage
    /// 'dense'
    #[allow(clippy::wrong_self_convention)]
    fn to_dense(&self) -> Self {
        let mut result = self.clone();
        result.tensor = result.tensor.to_dense();
        result
    }

    /// Copy the components into sparse storage.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     Independent sparse tensor with the same values and logical axes.
    ///
    /// Notes
    /// -----
    /// Sparse storage retains nonzero entries and represents other entries by zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.to_sparse().storage
    /// 'sparse'
    #[allow(clippy::wrong_self_convention)]
    fn to_sparse(&self) -> Self {
        let mut result = self.clone();
        result.tensor = result.tensor.to_sparse();
        result
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> text = repr(tensor)
    fn __repr__(&self) -> String {
        format!("Tensor({})", display::format_concrete_tensor(self, false))
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> text = str(tensor)
    fn __str__(&self) -> String {
        display::format_concrete_tensor(self, false)
    }

    /// Produce compact plain-text tensor notation.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    ///
    /// Returns
    /// -------
    /// str
    ///     Readable tensor text, suitable for logs or terminal output.
    ///
    /// Notes
    /// -----
    /// Only ports layout and default spacing are supported in this source
    /// format. Use to_html() or to_svg() for other layouts and custom gaps.
    ///
    /// Components follow the current view's logical interface order (the order
    /// of axes), independently of the canonical structure signature. Large
    /// static displays show a bounded preview rather than every component.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> output = tensor.format_tensor()
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    fn format_tensor(
        &self,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_plain_source_settings(&settings)?;
        Ok(display::format_concrete_tensor_with_settings(
            self, &settings,
        ))
    }

    /// Produce static Typst math source.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    ///
    /// Returns
    /// -------
    /// str
    ///     Typst source without an enclosing document; no compilation is performed.
    ///
    /// Notes
    /// -----
    /// Only ports layout and default spacing are supported in this source
    /// format. Use to_html() or to_svg() for other layouts and custom gaps.
    ///
    /// Components follow the current view's logical interface order (the order
    /// of axes), independently of the canonical structure signature. Large
    /// static displays show a bounded preview rather than every component.
    ///
    /// Static Typst output remains available independently of the HTML explorer.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> output = tensor.to_typst()
    #[pyo3(signature = (show_dimensions = None, *, settings = None))]
    fn to_typst(
        &self,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::validate_typst_source_settings(&settings)?;
        Ok(display::concrete_tensor_to_typst(self, &settings))
    }

    /// Create a rich display value for a notebook.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    /// notation_source : str, optional
    ///     Custom Typst notation source used by the rich renderer. This is
    ///     Typst code, not a filename; omit it for the supplied tensor notation.
    ///
    /// Returns
    /// -------
    /// FormattedOutput
    ///     Symbolica display wrapper with text and available rich representations.
    ///
    /// Notes
    /// -----
    /// Components follow logical axis order. Large static displays show a
    /// bounded preview rather than every component.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> output = tensor.formatted()
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn formatted(
        &self,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PythonFormattedOutput {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::format_concrete_tensor_output_rich(py, self, &settings, notation_source.as_deref())
    }

    /// Render the tensor as an HTML fragment.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    /// notation_source : str, optional
    ///     Custom Typst notation source used by the rich renderer. This is
    ///     Typst code, not a filename; omit it for the supplied tensor notation.
    ///
    /// Returns
    /// -------
    /// str
    ///     Self-contained output fragment for an HTML-capable notebook or page.
    ///
    /// Notes
    /// -----
    /// Mathematical rendering uses the optional Typst runtime. The returned
    /// string is not automatically displayed; pass it to the notebook's HTML
    /// or SVG display facility.
    ///
    /// Components follow logical axis order. Large static displays show a
    /// bounded preview rather than every component.
    ///
    /// The default is a responsive component explorer with a Memory grid /
    /// Matrix toggle. Dark mode follows the surrounding page. The color scale
    /// measures stored component payload bytes, not total process memory. Use
    /// DisplaySettings(tensor_view="matrix") for static mathematical HTML.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> output = tensor.to_html()
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn to_html(
        &self,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::concrete_tensor_to_html(py, self, &settings, notation_source.as_deref())
    }

    /// Render static mathematical tensor notation to SVG.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override settings.show_dimensions for this call. None retains the
    ///     selected settings, whose default omits dimensions.
    /// settings : DisplaySettings, optional
    ///     Tensor notation and display choices. Defaults to DisplaySettings().
    /// notation_source : str, optional
    ///     Custom Typst notation source used by the rich renderer. This is
    ///     Typst code, not a filename; omit it for the supplied tensor notation.
    ///
    /// Returns
    /// -------
    /// str
    ///     SVG markup suitable for embedding or saving to a file.
    ///
    /// Notes
    /// -----
    /// Mathematical rendering uses the optional Typst runtime. The returned
    /// string is not automatically displayed; pass it to the notebook's HTML
    /// or SVG display facility.
    ///
    /// Components follow logical axis order. Large static displays show a
    /// bounded preview rather than every component.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> output = tensor.to_svg()
    #[pyo3(signature = (show_dimensions = None, *, settings = None, notation_source = None))]
    fn to_svg(
        &self,
        py: Python<'_>,
        show_dimensions: Option<bool>,
        settings: Option<PyRef<'_, display::DisplaySettings>>,
        notation_source: Option<String>,
    ) -> PyResult<String> {
        let settings = display::resolved_settings(show_dimensions, settings.as_deref());
        display::concrete_tensor_to_svg(py, self, &settings, notation_source.as_deref())
    }

    fn _repr_html_(&self, py: Python<'_>) -> Option<String> {
        display::concrete_tensor_to_html(py, self, &display::DisplaySettings::default(), None).ok()
    }

    fn _repr_latex_(&self) -> Option<String> {
        Some(display::concrete_tensor_to_latex(self, false))
    }

    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        let text = if cycle {
            "...".to_string()
        } else {
            display::format_concrete_tensor(self, false)
        };
        pretty.call_method1("text", (text,))?;
        Ok(())
    }

    /// Count all components in the tensor shape.
    ///
    /// Returns
    /// -------
    /// int
    ///     Product of the axis dimensions, including implicit zeros.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> len(tensor)
    /// 4
    fn __len__(&self) -> usize {
        self.size().unwrap()
    }

    #[doc = python_doc!("Tensor.__getitem__")]
    #[gen_stub(skip)]
    fn __getitem__(&self, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        let layout = tensor_data_layout(self.descriptor.structure())?;
        let size = layout.size();
        let get_owned_canonical = |index: usize| {
            self.get_owned_linear(index.into())
                .or_else(|| match &self.tensor {
                    ParamOrConcrete::Concrete(RealOrComplexTensor::Real(DataTensor::Sparse(
                        tensor,
                    ))) => Some(ConcreteOrParam::Concrete(RealOrComplex::Real(tensor.zero))),
                    ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(
                        DataTensor::Sparse(tensor),
                    )) => Some(ConcreteOrParam::Concrete(RealOrComplex::Complex(
                        tensor.zero,
                    ))),
                    ParamOrConcrete::Param(tensor) => match &tensor.tensor {
                        DataTensor::Sparse(tensor) => {
                            Some(ConcreteOrParam::Param(tensor.zero.clone()))
                        }
                        DataTensor::Dense(_) => None,
                    },
                    ParamOrConcrete::Concrete(
                        RealOrComplexTensor::Real(DataTensor::Dense(_))
                        | RealOrComplexTensor::Complex(DataTensor::Dense(_)),
                    ) => None,
                })
        };
        let get_owned_linear = |logical_index: usize| -> PyResult<_> {
            let canonical_index = layout
                .logical_flat_to_storage_flat(logical_index)
                .map_err(tensor_data_layout_error)?;
            get_owned_canonical(canonical_index)
                .ok_or_else(|| PyIndexError::new_err("flat index out of bounds"))
        };

        let out = match item {
            SliceOrIntOrExpanded::Int(i) => get_owned_linear(i)?,
            SliceOrIntOrExpanded::Expanded(idxs) => {
                let index = layout
                    .logical_expanded_to_storage_flat(&idxs)
                    .map_err(tensor_data_layout_error)?;
                get_owned_canonical(index)
                    .ok_or_else(|| PyIndexError::new_err("expanded index out of bounds"))?
            }
            SliceOrIntOrExpanded::Slice(s) => {
                let r = s.indices(size as isize)?;
                let slice = (0..r.slicelength)
                    .map(|i| {
                        get_owned_linear((r.start + i as isize * r.step) as usize)
                            .map(TensorElements::from)
                    })
                    .collect::<PyResult<Vec<_>>>()?;

                return Ok(
                    Python::attach(|py| slice.into_pyobject(py).map(|a| a.unbind()))?.into_any(),
                );
            }
            SliceOrIntOrExpanded::Selection(selection) => {
                return self.select_components(&selection);
            }
        };

        Python::attach(|py| {
            TensorElements::from(out)
                .into_pyobject(py)
                .map(|a| a.unbind())
        })
    }

    #[doc = python_doc!("Tensor.__setitem__")]
    #[gen_stub(skip)]
    fn __setitem__<'py>(
        &mut self,
        item: Bound<'py, PyAny>,
        value: Bound<'py, PyAny>,
    ) -> eyre::Result<()> {
        TensorElements::validate_scalar(&value)?;
        let value = if matches!(&self.tensor, ParamOrConcrete::Param(_)) {
            ConcreteOrParam::Param(
                value
                    .extract::<ConvertibleToExpression>()?
                    .to_expression()
                    .expr,
            )
        } else if let Ok(v) = value.extract::<PythonExpression>() {
            ConcreteOrParam::Param(v.expr)
        } else if let Ok(v) = value.extract::<f64>() {
            ConcreteOrParam::Concrete(RealOrComplex::Real(v))
        } else if let Ok(v) = value.extract::<Complex<f64>>() {
            ConcreteOrParam::Concrete(RealOrComplex::Complex(v))
        } else {
            return Err(eyre!(
                "value must be a float for a real tensor, complex for a complex tensor, or Expression for a parametric tensor"
            ));
        };

        let layout = tensor_data_layout(self.descriptor.structure())?;
        let coefficient_error = |error| {
            eyre!(
                "assigned coefficient kind must match tensor storage (float for real, complex for complex, Expression for parametric): {error}"
            )
        };
        if let Ok(flat_index) = item.extract::<isize>() {
            let flat_index = Self::normalized_component_index(flat_index, layout.size())?;
            let canonical = layout
                .logical_flat_to_storage_flat(flat_index)
                .map_err(tensor_data_layout_error)?;
            self.tensor
                .set_flat(canonical.into(), value)
                .map_err(coefficient_error)
        } else if let Ok(indices) = item.extract::<Vec<isize>>() {
            let shape = layout.logical_shape();
            if indices.len() != shape.len() {
                return Err(eyre!(
                    "expected {} coordinates, got {}",
                    shape.len(),
                    indices.len()
                ));
            }
            let expanded_idxs = indices
                .iter()
                .zip(shape)
                .map(|(&index, &dimension)| Self::normalized_component_index(index, dimension))
                .collect::<PyResult<Vec<_>>>()?;
            let canonical = layout
                .logical_expanded_to_storage_flat(&expanded_idxs)
                .map_err(tensor_data_layout_error)?;
            self.tensor
                .set_flat(canonical.into(), value)
                .map_err(coefficient_error)
        } else {
            Err(eyre!("Index must be an integer"))
        }
    }

    /// Prepare repeated numerical evaluation of symbolic component formulas.
    ///
    /// Parameters
    /// ----------
    /// constants : mapping of Expression to Expression
    ///     Fixed values that are exact rational real or complex numbers. Use
    ///     parameters for values that cannot be represented this way.
    /// funs : mapping
    ///     Function definitions keyed by (function_symbol, name, argument_symbols),
    ///     with Expression bodies. The name field is accepted but not used.
    /// params : sequence of Expression
    ///     Inputs in the exact order expected by every evaluation row.
    /// iterations : int, default 100
    ///     Horner-optimization iterations.
    /// n_cores : int, default 4
    ///     Number of cores for expression optimization.
    /// verbose : bool, default False
    ///     Print optimization progress.
    ///
    /// Returns
    /// -------
    /// TensorEvaluator
    ///     Optimized batch evaluator retaining this tensor's shape and identity.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> evaluator.evaluate([[2.0]])[0][1]
    /// 4.0
    #[pyo3(signature =
           (constants,
           funs,
           params,
           iterations = 100,
           n_cores = 4,
           verbose = false),
           )]
    pub fn evaluator(
        &self,
        constants: HashMap<PythonExpression, PythonExpression>,
        funs: HashMap<(PolyVariable, String, Vec<PolyVariable>), PythonExpression>,
        params: Vec<PythonExpression>,
        iterations: usize,
        n_cores: usize,
        verbose: bool,
    ) -> PyResult<SpensoExpressionEvaluator> {
        let mut fn_map = FunctionMap::new();

        for (k, v) in &constants {
            if let Ok(r) = SymComplex::<Rational>::try_from(v.expr.clone()) {
                fn_map
                    .add_aliases([(k.expr.clone(), Atom::num(r))])
                    .map_err(|e| {
                        exceptions::PyValueError::new_err(format!("Could not add constant: {}", e))
                    })?;
            } else {
                Err(exceptions::PyValueError::new_err(
                    "Constants must be rationals. If this is not possible, pass the value as a parameter",
                ))?
            }
        }

        for ((symbol, _rename, args), body) in &funs {
            let symbol = symbol
                .get_id()
                .ok_or(exceptions::PyValueError::new_err(format!(
                    "Bad function name {}",
                    symbol
                )))?;
            let args: Vec<_> = args
                .iter()
                .map(|x| {
                    x.get_id().ok_or(exceptions::PyValueError::new_err(format!(
                        "Bad function name {}",
                        symbol
                    )))
                })
                .collect::<Result<_, _>>()?;

            fn_map
                .add_function(symbol, args, body.expr.clone())
                .map_err(|e| {
                    exceptions::PyValueError::new_err(format!("Could not add function: {}", e))
                })?;
        }

        let settings = OptimizationSettings::new()
            .horner_iterations(iterations)
            .cores(n_cores)
            .verbose(verbose);

        let params: Vec<_> = params.iter().map(|x| x.expr.clone()).collect();

        use spenso::tensors::parametric::to_param::ToParam;
        let mut evaltensor = match &self.tensor {
            ParamOrConcrete::Param(s) => s.to_evaluation_tree(&fn_map, &params),
            ParamOrConcrete::Concrete(c) => {
                c.clone().to_param().to_evaluation_tree(&fn_map, &params)
            }
        }
        .map_err(|e| {
            exceptions::PyValueError::new_err(format!("Could not create evaluator: {e}"))
        })?;

        evaltensor.optimize_horner_scheme(&settings);

        evaltensor.common_subexpression_elimination();
        let linear = evaltensor.linearize(&settings);
        // Decide from exact coefficients before converting to floating point.
        let real_coefficients = Cell::new(true);
        let real = linear.clone().map_coeff(&|coefficient| {
            if !coefficient.im.is_zero() {
                real_coefficients.set(false);
            }
            coefficient.re.to_f64()
        });
        Ok(SpensoExpressionEvaluator {
            eval: real_coefficients.get().then_some(real),
            eval_complex: linear
                .clone()
                .map_coeff(&|x| Complex::new(x.re.to_f64(), x.im.to_f64())),
            eval_rat: linear,
            parameters: params,
            descriptor: self.descriptor.clone(),
            descriptor_name: self.descriptor_name,
            descriptor_args: self.descriptor_args.clone(),
        })
    }

    /// Extract the only component of a rank-zero tensor.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> scalar = sp.Tensor.dense(sp.TensorName("docs::scalar")(), [2.0])
    /// >>> value = scalar.scalar()
    ///
    /// Returns
    /// -------
    /// Expression
    ///     Scalar value converted to a Symbolica expression, including numeric data.
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     The tensor has external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorExpression
    /// >>> TensorExpression(3).to_tensor().scalar() == 3
    /// True
    fn scalar(&self) -> PyResult<PythonExpression> {
        self.tensor
            .clone()
            .scalar()
            .map(|r| PythonExpression { expr: r.into() })
            .ok_or_else(|| PyRuntimeError::new_err("No scalar found"))
    }

    /// Assign labels to the unresolved external axes.
    ///
    /// Parameters
    /// ----------
    /// *indices : int, str, Expression, Slot, or AUTO
    ///     One label per unresolved axis, in logical order. AUTO leaves the
    ///     corresponding axis unchanged. A Slot must have the matching
    ///     representation. Repeated compatible labels cause contraction.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     An indexed network retaining the component data; the original is unchanged.
    ///
    /// Notes
    /// -----
    /// Existing explicit labels are left in place. Use reindex to replace them.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = tensor.index("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern = None))]
    fn index(
        &self,
        indices: &Bound<'_, PyTuple>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor_reference(self.clone())?.index_network(indices, intern)
    }

    /// Assign labels to all external axes.
    ///
    /// Parameters
    /// ----------
    /// *indices : int, str, Expression, Slot, or AUTO
    ///     One label per external axis, in logical order. AUTO leaves the
    ///     corresponding axis unchanged. A Slot must have the matching
    ///     representation. Repeated compatible labels cause contraction.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     An indexed network retaining the component data; the original is unchanged.
    ///
    /// Notes
    /// -----
    /// Existing labels are replaced as well as unresolved axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = tensor.reindex("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern=None))]
    fn reindex(
        &self,
        indices: &Bound<'_, PyTuple>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor_reference(self.clone())?.reindex(indices, intern)
    }

    /// Rename selected external labels without changing the tensor rank.
    ///
    /// Parameters
    /// ----------
    /// mapping : dict
    ///     Old label to new label. A Slot key also restricts the representation.
    ///     Keys must identify existing external axes. Each axis may be renamed
    ///     once; renaming must preserve the external interface.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A new value with renamed axes; the original is unchanged.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = Tensor.dense(A("i", "j"), [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed.rename_indices({"i": "k"}).rank
    /// 2
    #[pyo3(signature = (mapping, *, intern=None))]
    fn rename_indices(
        &self,
        mapping: &Bound<'_, PyDict>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<Self> {
        SpensoNet::from_tensor_reference(self.clone())?
            .rename_indices(mapping, intern)?
            .result_tensor(None)
    }

    /// Reorder the external axes in logical component order.
    ///
    /// Parameters
    /// ----------
    /// axes : sequence of int
    ///     Each current axis position exactly once, in the desired new order.
    ///     For a matrix, [1, 0] exchanges the two axes.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A new tensor view with the reordered interface and correspondingly
    ///     reordered component access.
    ///
    /// Notes
    /// -----
    /// This permutes axes; it does not conjugate component values.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> transposed = tensor.permute_axes([1, 0])
    /// >>> transposed.shape
    /// (2, 2)
    fn permute_axes(&self, axes: Vec<usize>) -> PyResult<Self> {
        SpensoNet::from_tensor_reference(self.clone())?
            .permute_axes(axes)?
            .result_tensor(None)
    }

    /// Implement ``copy.copy`` using an independent Tensor copy.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     Copy of the component data and metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> import copy
    /// >>> copied = copy.copy(tensor)
    fn __copy__(&self) -> Self {
        self.clone()
    }

    /// Assign labels to the unresolved external axes.
    ///
    /// Parameters
    /// ----------
    /// *indices : int, str, Expression, Slot, or AUTO
    ///     One label per unresolved axis, in logical order. AUTO leaves the
    ///     corresponding axis unchanged. A Slot must have the matching
    ///     representation. Repeated compatible labels cause contraction.
    /// intern : {"indices", "flattened"} or None, optional
    ///     Encode compound index labels before checking the tensor structure.
    ///     "indices" retains reversible payloads; "flattened" creates readable names.
    ///     Both preserve index tags and leave scalar arguments and tensor heads intact.
    ///     None leaves labels unchanged and requires ordinary atomic index labels.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     An indexed network retaining the component data; the original is unchanged.
    ///
    /// Notes
    /// -----
    /// Equivalent to ``index(*indices)``.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = tensor("i", "j")
    /// >>> indexed.rank
    /// 2
    #[pyo3(signature = (*indices, intern = None))]
    fn __call__(
        &self,
        indices: &Bound<'_, PyTuple>,
        intern: Option<crate::simplification::Intern>,
    ) -> PyResult<SpensoNet> {
        self.index(indices, intern)
    }

    /// Negate every tensor component.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     Negated symbolic algebra retaining component data lazily.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> negated = -tensor
    fn __neg__(&self) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__neg__()
    }

    /// Add tensors with compatible external interfaces.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Scalar zero acts as the additive identity. Other scalars can be added
    /// only to rank-zero tensors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor + tensor
    fn __add__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.add_network(rhs.to_net(), false)
    }

    /// Implement reflected addition.
    ///
    /// Parameters
    /// ----------
    /// lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Scalar zero acts as the additive identity. Other scalars can be added
    /// only to rank-zero tensors.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor + tensor
    fn __radd__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .add_network(SpensoNet::from_tensor(self.clone())?, false)
    }

    /// Subtract tensors with compatible external interfaces.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Subtraction requires matching tensor axes; a nonzero scalar cannot
    /// be subtracted from a tensor with external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor - tensor
    fn __sub__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.add_network(rhs.to_net(), true)
    }

    /// Implement reflected subtraction.
    ///
    /// Parameters
    /// ----------
    /// lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Subtraction requires matching tensor axes; a nonzero scalar cannot
    /// be subtracted from a tensor with external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor - tensor
    fn __rsub__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .add_network(SpensoNet::from_tensor(self.clone())?, true)
    }

    /// Multiply tensors, contracting unambiguous compatible axes.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Matching explicit labels contract. Compatible unresolved axes are
    /// paired only when the choice is unambiguous. Use outer(), contract(),
    /// or compose() to make the intended pairing explicit.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor * 2
    fn __mul__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.multiply_network(rhs.to_net())
    }

    /// Implement reflected multiplication.
    ///
    /// Parameters
    /// ----------
    /// lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// Matching explicit labels contract. Compatible unresolved axes are
    /// paired only when the choice is unambiguous. Use outer(), contract(),
    /// or compose() to make the intended pairing explicit.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = 2 * tensor
    fn __rmul__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .multiply_network(SpensoNet::from_tensor(self.clone())?)
    }

    /// Divide tensor components by a scalar expression.
    ///
    /// Parameters
    /// ----------
    /// rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// The denominator must be scalar. For scalar divided by tensor, the
    /// tensor must also have rank zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor / 2
    fn __truediv__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__truediv__(rhs)
    }

    /// Implement reflected division.
    ///
    /// Parameters
    /// ----------
    /// lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    ///     Other operand. Component data are retained in a lazy network.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     New lazy calculation; operands are unchanged.
    ///
    /// Notes
    /// -----
    /// The denominator must be scalar. For scalar divided by tensor, the
    /// tensor must also have rank zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = 2 / TensorExpression(3).to_tensor()
    fn __rtruediv__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__rtruediv__(lhs)
    }

    /// Form an outer product without implicit contractions between the operands.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Tensor whose axes follow this tensor's axes.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy outer product retaining the operands' component data.
    ///
    /// Notes
    /// -----
    /// Use this when compatible unresolved axes should remain independent.
    /// Choose explicit distinct labels if you need to refer to them separately.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> product = tensor.outer(tensor)
    /// >>> product.rank
    /// 4
    fn outer(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.outer(rhs)
    }

    /// Contract one chosen pair of axes between two tensors.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Tensor to contract with this tensor.
    /// left, right : int
    ///     Zero-based axis positions in this tensor and rhs, respectively.
    ///     The representations must be compatible under contraction.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy contraction retaining component data.
    ///
    /// Notes
    /// -----
    /// Other axes remain external. Axis numbers refer to ``axes``.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> contracted = tensor.contract_ports(tensor, left=1, right=0)
    /// >>> contracted.rank
    /// 2
    #[pyo3(signature = (rhs, *, left, right))]
    fn contract_ports(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: usize,
        right: usize,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.contract_ports(rhs, left, right)
    }

    /// Multiply two tensors along explicitly selected matrix channels.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Next factor in the ordered matrix product.
    /// left, right : tuple of int and int
    ///     (input_axis, output_axis) in this tensor and rhs, respectively.
    ///     The left output contracts with the right input. Other axes are
    ///     retained as spectator axes.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy ordered matrix product retaining component data.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> product = tensor.compose(tensor, left=(0, 1), right=(0, 1))
    /// >>> product.rank
    /// 2
    #[pyo3(signature = (rhs, *, left, right))]
    fn compose(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: (usize, usize),
        right: (usize, usize),
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.compose(rhs, left, right)
    }

    /// Contract two rank-one tensors using their representation pairing.
    ///
    /// Parameters
    /// ----------
    /// rhs : TensorExpression, Tensor, or TensorNetwork
    ///     Rank-one tensor in a compatible space.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     A lazy scalar contraction, including the metric signs of the space.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorName, Tensor, TensorNetwork, Representation
    /// >>> vector = Tensor.dense(TensorName.vector("dot_v")(Representation.euc(2)), [1.0, 2.0])
    /// >>> vector.dot(vector).to_tensor().scalar() == 5
    /// True
    fn dot(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.dot(rhs)
    }

    /// Close a pair of matrix axes on this tensor.
    ///
    /// Parameters
    /// ----------
    /// channel : tuple of int and int, optional
    ///     (input_axis, output_axis) to contract. If omitted, the matrix channel
    ///     must be uniquely determined by the representations. Specify it when
    ///     more than one pairing is possible.
    ///
    /// Returns
    /// -------
    /// TensorNetwork
    ///     The traced tensor; spectator axes remain external.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.trace(channel=(0, 1)).is_scalar
    /// True
    #[pyo3(signature = (*, channel = None))]
    fn trace(&self, channel: Option<(usize, usize)>) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.trace(channel)
    }
}

/// Optimized numerical evaluation of a tensor's component formulas.
///
/// Create this object with ``Tensor.evaluator``. Each input row supplies values
/// for ``parameters`` and produces one Tensor with ``output_shape``. Real
/// coefficients support both real and complex evaluation; complex coefficients
/// require ``evaluate_complex``.
///
/// Examples
/// --------
/// >>> from symbolica import S
/// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
/// >>> x = S("x")
/// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
/// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
/// >>> evaluator.evaluate([[2.0], [3.0]])[1][1]
/// 9.0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "TensorEvaluator",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub struct SpensoExpressionEvaluator {
    pub eval_rat: LinearizedEvalTensor<SymComplex<Rational>, ShadowedStructure<AbstractIndex>>,
    pub eval: Option<LinearizedEvalTensor<f64, ShadowedStructure<AbstractIndex>>>,
    pub eval_complex: LinearizedEvalTensor<Complex<f64>, ShadowedStructure<AbstractIndex>>,
    parameters: Vec<Atom>,
    descriptor: SymbolicTensor<PartialStructure>,
    descriptor_name: Option<Symbol>,
    descriptor_args: Vec<Atom>,
}

impl SpensoExpressionEvaluator {
    fn validate_inputs<T>(expected: usize, inputs: &[Vec<T>]) -> PyResult<()> {
        if let Some((row, values)) = inputs
            .iter()
            .enumerate()
            .find(|(_, row)| row.len() != expected)
        {
            return Err(exceptions::PyValueError::new_err(format!(
                "input row {row} has {} parameters; expected {expected}",
                values.len()
            )));
        }
        Ok(())
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoExpressionEvaluator {
    /// Symbolic inputs in the required evaluation order.
    ///
    /// Returns
    /// -------
    /// list of Expression
    ///     One entry for each value in an input row.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> evaluator.parameters == [x]
    /// True
    #[getter]
    fn parameters(&self) -> Vec<PythonExpression> {
        self.parameters.iter().cloned().map(Into::into).collect()
    }

    /// Number of values required in each evaluation row.
    ///
    /// Returns
    /// -------
    /// int
    ///     Length of parameters.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> evaluator.input_size
    /// 1
    #[getter]
    fn input_size(&self) -> usize {
        self.parameters.len()
    }

    /// Logical component shape of each result tensor.
    ///
    /// Returns
    /// -------
    /// tuple of int
    ///     Shape excluding the batch dimension.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> evaluator.output_shape
    /// (2,)
    #[getter]
    fn output_shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            tensor_data_layout(self.descriptor.structure())?.logical_shape(),
        )
    }

    /// Whether a real-valued evaluation path is available.
    ///
    /// Returns
    /// -------
    /// bool
    ///     False when the formulas contain complex coefficients; use
    ///     evaluate_complex in that case.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> evaluator.supports_real
    /// True
    #[getter]
    fn supports_real(&self) -> bool {
        self.eval.is_some()
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> x = S("docs::x")
    /// >>> tensor = sp.Tensor.dense(A, [x, E("0"), E("0"), x + 1])
    /// >>> evaluator = tensor.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> text = repr(evaluator)
    fn __repr__(&self) -> String {
        format!(
            "TensorEvaluator({}, real={}, complex=True)",
            display::format_structured(&self.descriptor, false),
            if self.eval.is_some() { "True" } else { "False" }
        )
    }

    /// Evaluate a batch of real input rows.
    ///
    /// Parameters
    /// ----------
    /// inputs : sequence of sequences of float
    ///     One row per evaluation, with input_size entries in parameters order.
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     One real component tensor per input row, each with output_shape.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     A row has the wrong size, or the formulas contain complex
    ///     coefficients. Use evaluate_complex for those formulas.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> evaluator.evaluate([[2.0]])[0][1]
    /// 4.0
    fn evaluate(&mut self, inputs: Vec<Vec<f64>>) -> PyResult<Vec<Spensor>> {
        Self::validate_inputs(self.parameters.len(), &inputs)?;
        let eval = self.eval.as_mut().ok_or(exceptions::PyValueError::new_err(
            "Evaluator contains complex coefficients. Use evaluate_complex instead.",
        ))?;

        Ok(inputs
            .iter()
            .map(|s| {
                Spensor::from_storage_with_descriptor(
                    MixedTensor::Concrete(RealOrComplexTensor::Real(eval.evaluate(s))),
                    self.descriptor.clone(),
                    self.descriptor_name,
                    self.descriptor_args.clone(),
                )
            })
            .collect())
    }

    /// Evaluate a batch of complex input rows.
    ///
    /// Parameters
    /// ----------
    /// inputs : sequence of sequences of complex
    ///     One row per evaluation, with input_size entries in parameters order.
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     One complex component tensor per input row, each with output_shape.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     A row has the wrong size.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> evaluator.evaluate_complex([[2.0]])[0][1]
    /// (4+0j)
    fn evaluate_complex(&mut self, inputs: Vec<Vec<Complex<f64>>>) -> PyResult<Vec<Spensor>> {
        Self::validate_inputs(self.parameters.len(), &inputs)?;
        let eval = &mut self.eval_complex;

        Ok(inputs
            .iter()
            .map(|s| {
                Spensor::from_storage_with_descriptor(
                    MixedTensor::Concrete(RealOrComplexTensor::Complex(
                        eval.evaluate(s).map_data(|value| value),
                    )),
                    self.descriptor.clone(),
                    self.descriptor_name,
                    self.descriptor_args.clone(),
                )
            })
            .collect())
    }
}

#[cfg(feature = "native")]
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoExpressionEvaluator {
    /// Compile and load native C++ code for repeated tensor evaluation.
    ///
    /// Parameters
    /// ----------
    /// function_name : str
    ///     Exported C++ function name.
    /// filename : str
    ///     Path for generated complex-evaluation C++ source. When real evaluation
    ///     is supported, an additional file with suffix .real.cpp is written.
    /// library_name : str
    ///     Output library path. A separate .real library is generated when supported.
    /// inline_asm : str, default "default"
    ///     Assembly mode: "default", "x64", "aarch64", or "none".
    /// optimization_level : int, default 3
    ///     Compiler optimization level.
    /// compiler_path : str, optional
    ///     C++ compiler executable; omit to use Symbolica's default compiler.
    /// custom_header : str, optional
    ///     Additional C++ header source included in the generated code.
    ///
    /// Returns
    /// -------
    /// CompiledTensorEvaluator
    ///     Loaded native evaluator with the same input ordering and result shape.
    ///
    /// Notes
    /// -----
    /// This writes source and compiled library files to the supplied paths.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from tempfile import TemporaryDirectory
    /// >>> directory = TemporaryDirectory()
    /// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
    /// >>> compiled.input_size
    /// 1
    /// >>> directory.cleanup()
    #[pyo3(signature =
        (function_name,
        filename,
        library_name,
        // number_type,
        inline_asm = "default",
        optimization_level = 3,
        compiler_path = None,
        // compiler_flags = None,
        custom_header = None,
        // cuda_number_of_evaluations = 1,
        // cuda_block_size = 512
    ))]
    #[allow(clippy::too_many_arguments)]
    fn compile(
        &self,
        function_name: &str,
        filename: &str,
        library_name: &str,
        // number_type: &str,
        inline_asm: &str,
        optimization_level: u8,
        compiler_path: Option<&str>,
        // compiler_flags: Option<Vec<String>>,
        custom_header: Option<String>,
        // cuda_number_of_evaluations: usize,
        // cuda_block_size: usize,
        // py: Python<'_>,
    ) -> PyResult<SpensoCompiledExpressionEvaluator> {
        let mut options = CompileOptions::new().optimization_level(optimization_level as usize);

        if let Some(compiler_path) = compiler_path {
            options = options.compiler(compiler_path.to_string());
        }
        let inline_asm = match inline_asm.to_lowercase().as_str() {
            "default" => InlineASM::default(),
            "x64" => InlineASM::X64,
            "aarch64" => InlineASM::AArch64,
            "none" => InlineASM::None,
            _ => {
                return Err(exceptions::PyValueError::new_err(
                    "Invalid inline assembly type specified.",
                ));
            }
        };

        let real = self
            .eval
            .as_ref()
            .map(|eval| -> PyResult<_> {
                eval.export_cpp::<f64>(
                    &format!("{filename}.real.cpp"),
                    &format!("{function_name}_real"),
                    ExportSettings::new()
                        .include_header(true)
                        .inline_asm(inline_asm)
                        .custom_header(custom_header.clone()),
                )
                .map_err(|e| exceptions::PyValueError::new_err(format!("Export error: {e}")))?
                .compile(&format!("{library_name}.real"), options.clone())
                .map_err(|e| exceptions::PyValueError::new_err(format!("Compilation error: {e}")))?
                .load()
                .map_err(|e| {
                    exceptions::PyValueError::new_err(format!("Library loading error: {e}"))
                })
            })
            .transpose()?;

        Ok(SpensoCompiledExpressionEvaluator {
            real,
            parameters: self.parameters.clone(),
            eval: self
                .eval_complex
                .export_cpp::<Complex<f64>>(
                    filename,
                    function_name,
                    ExportSettings::new()
                        .include_header(true)
                        .inline_asm(inline_asm)
                        .custom_header(custom_header),
                )
                .map_err(|e| exceptions::PyValueError::new_err(format!("Export error: {}", e)))?
                .compile(library_name, options)
                .map_err(|e| {
                    exceptions::PyValueError::new_err(format!("Compilation error: {}", e))
                })?
                .load()
                .map_err(|e| {
                    exceptions::PyValueError::new_err(format!("Library loading error: {}", e))
                })?,
            descriptor: self.descriptor.clone(),
            descriptor_name: self.descriptor_name,
            descriptor_args: self.descriptor_args.clone(),
        })
    }
}

/// A tensor evaluator compiled into a loaded native library.
///
/// Create this object with ``TensorEvaluator.compile``. Its input ordering,
/// output shape, and real/complex evaluation rules match TensorEvaluator.
/// Compilation requires a C++ compiler and writes source and library files.
///
/// Examples
/// --------
/// >>> from symbolica import S
/// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
/// >>> x = S("x")
/// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
/// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
/// >>> from tempfile import TemporaryDirectory
/// >>> directory = TemporaryDirectory()
/// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
/// >>> compiled.evaluate([[2.0]])[0][1]
/// 4.0
/// >>> directory.cleanup()
#[cfg(feature = "native")]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "CompiledTensorEvaluator",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub struct SpensoCompiledExpressionEvaluator {
    pub eval: EvalTensor<CompiledComplexEvaluatorSpenso, ShadowedStructure<AbstractIndex>>,
    real: Option<
        EvalTensor<symbolica::evaluate::CompiledRealEvaluator, ShadowedStructure<AbstractIndex>>,
    >,
    parameters: Vec<Atom>,
    descriptor: SymbolicTensor<PartialStructure>,
    descriptor_name: Option<Symbol>,
    descriptor_args: Vec<Atom>,
}

#[cfg(feature = "native")]
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoCompiledExpressionEvaluator {
    /// Symbolic inputs in the required evaluation order.
    ///
    /// Returns
    /// -------
    /// list of Expression
    ///     One entry for each value in an input row.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from tempfile import TemporaryDirectory
    /// >>> directory = TemporaryDirectory()
    /// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
    /// >>> compiled.parameters == [x]
    /// True
    /// >>> directory.cleanup()
    #[getter]
    fn parameters(&self) -> Vec<PythonExpression> {
        self.parameters.iter().cloned().map(Into::into).collect()
    }

    /// Number of values required in each evaluation row.
    ///
    /// Returns
    /// -------
    /// int
    ///     Length of parameters.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from tempfile import TemporaryDirectory
    /// >>> directory = TemporaryDirectory()
    /// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
    /// >>> compiled.input_size
    /// 1
    /// >>> directory.cleanup()
    #[getter]
    fn input_size(&self) -> usize {
        self.parameters.len()
    }

    /// Logical component shape of each result tensor.
    ///
    /// Returns
    /// -------
    /// tuple of int
    ///     Shape excluding the batch dimension.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from tempfile import TemporaryDirectory
    /// >>> directory = TemporaryDirectory()
    /// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
    /// >>> compiled.output_shape
    /// (2,)
    /// >>> directory.cleanup()
    #[getter]
    fn output_shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            tensor_data_layout(self.descriptor.structure())?.logical_shape(),
        )
    }

    /// Whether a real-valued evaluation path is available.
    ///
    /// Returns
    /// -------
    /// bool
    ///     False when the formulas contain complex coefficients; use
    ///     evaluate_complex in that case.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from tempfile import TemporaryDirectory
    /// >>> directory = TemporaryDirectory()
    /// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
    /// >>> compiled.supports_real
    /// True
    /// >>> directory.cleanup()
    #[getter]
    fn supports_real(&self) -> bool {
        self.real.is_some()
    }

    /// Evaluate a batch of real input rows.
    ///
    /// Parameters
    /// ----------
    /// inputs : sequence of sequences of float
    ///     One row per evaluation, with input_size entries in parameters order.
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     One real component tensor per input row, each with output_shape.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     A row has the wrong size, or the formulas contain complex
    ///     coefficients. Use evaluate_complex for those formulas.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from tempfile import TemporaryDirectory
    /// >>> directory = TemporaryDirectory()
    /// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
    /// >>> compiled.evaluate([[2.0]])[0][1]
    /// 4.0
    /// >>> directory.cleanup()
    fn evaluate(&mut self, inputs: Vec<Vec<f64>>) -> PyResult<Vec<Spensor>> {
        SpensoExpressionEvaluator::validate_inputs(self.parameters.len(), &inputs)?;
        let eval = self.real.as_mut().ok_or_else(|| {
            exceptions::PyValueError::new_err(
                "Evaluator contains complex coefficients. Use evaluate_complex instead.",
            )
        })?;
        Ok(inputs
            .iter()
            .map(|row| {
                Spensor::from_storage_with_descriptor(
                    MixedTensor::Concrete(RealOrComplexTensor::Real(eval.evaluate(row))),
                    self.descriptor.clone(),
                    self.descriptor_name,
                    self.descriptor_args.clone(),
                )
            })
            .collect())
    }

    /// Return a description of the tensor layout and available numerical kernels.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> x = S("docs::x")
    /// >>> tensor = sp.Tensor.dense(A, [x, E("0"), E("0"), x + 1])
    /// >>> evaluator = tensor.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from pathlib import Path
    /// >>> from tempfile import TemporaryDirectory
    /// >>> with TemporaryDirectory() as directory:
    /// ...     folder = Path(directory)
    /// ...     compiled = evaluator.compile("docs_eval", str(folder / "eval.cpp"), str(folder / "eval.so"), inline_asm="none", optimization_level=0)
    /// ...     results = compiled.evaluate([[2.0], [3.0]])
    /// ...     text = repr(compiled)
    fn __repr__(&self) -> String {
        format!(
            "CompiledTensorEvaluator({}, real={}, complex=True)",
            display::format_structured(&self.descriptor, false),
            if self.real.is_some() { "True" } else { "False" }
        )
    }

    /// Evaluate a batch of complex input rows.
    ///
    /// Parameters
    /// ----------
    /// inputs : sequence of sequences of complex
    ///     One row per evaluation, with input_size entries in parameters order.
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     One complex component tensor per input row, each with output_shape.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     A row has the wrong size.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.tensor import Tensor, TensorName, Representation
    /// >>> x = S("x")
    /// >>> values = Tensor.dense(TensorName.vector("eval_v")(Representation.euc(2)), [x, x**2])
    /// >>> evaluator = values.evaluator({}, {}, [x], iterations=1, n_cores=1)
    /// >>> from tempfile import TemporaryDirectory
    /// >>> directory = TemporaryDirectory()
    /// >>> compiled = evaluator.compile("eval", directory.name + "/eval.cpp", directory.name + "/eval", inline_asm="none")
    /// >>> compiled.evaluate_complex([[2.0]])[0][1]
    /// (4+0j)
    /// >>> directory.cleanup()
    fn evaluate_complex(&mut self, inputs: Vec<Vec<Complex<f64>>>) -> PyResult<Vec<Spensor>> {
        SpensoExpressionEvaluator::validate_inputs(self.parameters.len(), &inputs)?;
        Ok(inputs
            .iter()
            .map(|s| {
                Spensor::from_storage_with_descriptor(
                    MixedTensor::Concrete(RealOrComplexTensor::Complex(
                        self.eval.evaluate(s).map_data(|value| value),
                    )),
                    self.descriptor.clone(),
                    self.descriptor_name,
                    self.descriptor_args.clone(),
                )
            })
            .collect())
    }
}

#[cfg(feature = "python_stubgen")]
submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<crate::structure::SpensoRepresentation>,
        attrs: &[],
        getters: &[],
        setters: &[],
        file: file!(),
        line: line!(),
        column: column!(),
        methods: &[
            MethodInfo {
                name: "__call__",
                parameters: &[
                    ParameterInfo {
                        name: "aind",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: || i64::type_input() | String::type_input(),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: structure::SpensoSlot::type_output,
                doc: python_doc!("Representation.__call__"),
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__call__",
                parameters: &[
                    ParameterInfo {
                        name: "aind",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: PythonExpression::type_input
                    },
                ],
                r#type: MethodType::Instance,
                r#return: || PythonExpression::type_output()| structure::SpensoSlot::type_output(),
                doc: python_doc!("Representation.__call__"),
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            }
        ]
    }
}

// static NONE: LazyLock<String> = LazyLock::new(|| "None".to_string());

#[cfg(feature = "python_stubgen")]
submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<Spensor>,
        attrs: &[],
        getters: &[],
        setters: &[],
        file: file!(),
        line: line!(),
        column: column!(),
        methods: &[
            MethodInfo {
                name: "__iter__",
                parameters: &[],
                r#type: MethodType::Instance,
                r#return:||
                TypeInfo {
                    name: "typing.Iterator[Expression | float | complex]".into(),
                    import: std::collections::HashSet::new(),
                },
                doc: python_doc!("Tensor.__iter__"),
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__getitem__",
                parameters: &[
                    ParameterInfo {
                        name: "item",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: || TypeInfo::builtin("slice"),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: Vec::<TensorElements>::type_output,
                doc: python_doc!("Tensor.__getitem__"),
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__getitem__",
                parameters: &[
                    ParameterInfo {
                        name: "item",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: || Vec::<usize>::type_input()|usize::type_input()
                    },
                ],
                r#type: MethodType::Instance,
                r#return: TensorElements::type_output,
                doc: python_doc!("Tensor.__getitem__"),
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__setitem__",
                parameters: &[
                    ParameterInfo {
                        name: "item",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: || usize::type_input()|Vec::<usize>::type_input()

                    },
                    ParameterInfo {
                        name: "value",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: ||TensorElements::type_input()
                    },
                ],
                r#type: MethodType::Instance,
                r#return: TypeInfo::none,
                doc: python_doc!("Tensor.__setitem__"),
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
        ]
    }
}

/// Gather the unified Spenso and Idenso Python API registered by Spynso3.
#[cfg(feature = "python_stubgen")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    let mut info = pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica".to_owned(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")),
    )?;
    if let Some(module) = info.modules.get_mut("symbolica.community.tensor") {
        SpensoModule::prepare_stub_module(module);
    }
    Ok(info)
}

#[cfg(feature = "python_stubgen")]
impl SpensoModule {
    /// Render the same checked signature inventory for docs and community packages.
    pub fn stub_source(module: &pyo3_stub_gen::generate::Module) -> String {
        module
            .to_string()
            .lines()
            .map(|line| {
                // Legacy TypeVars require an unannotated assignment; the upstream
                // variable generator currently always emits an annotation.
                if let Some((name, value)) = line.split_once(": typing.TypeVar = ") {
                    return format!("{name} = {value}");
                }
                // These methods deliberately specialize Expression's scalar API. Keep
                // both checkers' narrowly scoped override annotations on the definition.
                if line.trim_start().starts_with("def ")
                    && let Some(prefix) =
                        line.strip_suffix("# type: ignore[override,invalid-method-override]")
                {
                    format!(
                        "{prefix}# type: ignore[override]  # ty: ignore[invalid-method-override]"
                    )
                } else {
                    line.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n")
            + "\n"
    }

    /// PyO3 simple enums are classes with constants, not Python `enum.Enum` subclasses.
    /// Preserve their actual surface instead of promising `.name`, `.value`, or iteration.
    pub fn prepare_stub_module(module: &mut pyo3_stub_gen::generate::Module) {
        use pyo3_stub_gen::generate::{ClassDef, MemberDef, MethodDef, MethodType};
        for (id, enumeration) in std::mem::take(&mut module.enum_) {
            let mut class = ClassDef {
                name: enumeration.name,
                doc: enumeration.doc,
                attrs: enumeration.attrs,
                getter_setters: Default::default(),
                methods: Default::default(),
                bases: Vec::new(),
                classes: Vec::new(),
                match_args: None,
                subclass: false,
            };
            for (name, doc) in enumeration.variants {
                class.attrs.push(MemberDef {
                    name,
                    doc,
                    r#type: TypeInfo {
                        name: format!("typing.ClassVar[{}]", enumeration.name),
                        import: ["typing".into()].into(),
                    },
                    default: None,
                    deprecated: None,
                });
            }
            for getter in enumeration.getters {
                let name = getter.name.to_owned();
                class.getter_setters.entry(name).or_default().0 = Some(getter);
            }
            for setter in enumeration.setters {
                let name = setter.name.to_owned();
                class.getter_setters.entry(name).or_default().1 = Some(setter);
            }
            for method in enumeration.methods {
                class
                    .methods
                    .entry(method.name.to_owned())
                    .or_default()
                    .push(method);
            }
            class.methods.insert(
                "__int__".to_owned(),
                vec![MethodDef {
                    name: "__int__",
                    parameters: Default::default(),
                    r#return: i64::type_output(),
                    doc: "Return the underlying integer discriminant.",
                    r#type: MethodType::Instance,
                    is_async: false,
                    deprecated: None,
                    type_ignored: None,
                    is_overload: false,
                }],
            );
            module.class.insert(id, class);
        }
        stubs::refine(module);
    }
}

static CITATIONS_USED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

#[inline]
pub(crate) fn record_usage() {
    use std::sync::atomic::Ordering;
    if !CITATIONS_USED.load(Ordering::Relaxed) {
        CITATIONS_USED.store(true, Ordering::Relaxed);
    }
}

impl SpensoModule {
    /// Record tensor operations performed through another community binding.
    #[inline]
    pub fn record_usage() {
        record_usage();
    }
}

#[cfg(test)]
mod tests {
    use pyo3::types::PyInt;

    use super::*;

    #[test]
    fn complex_tensor_accepts_only_complex_assignment() {
        let interface = PartialStructure::from_logical_slots([]);
        let mut tensor = Spensor::from_storage_with_descriptor(
            ParamOrConcrete::new_scalar(ConcreteOrParam::Concrete(RealOrComplex::Complex(
                Complex::new(1., 2.),
            ))),
            SymbolicTensor::new(Atom::Zero, interface),
            None,
            Vec::new(),
        );

        Python::initialize();
        Python::attach(|py| {
            tensor
                .__setitem__(
                    PyInt::new(py, 0).into_any(),
                    PyComplex::from_doubles(py, 3., 4.).into_any(),
                )
                .unwrap();
            let value = tensor.__getitem__(SliceOrIntOrExpanded::Int(0)).unwrap();
            assert_eq!(
                value.bind(py).extract::<Complex<f64>>().unwrap(),
                Complex::new(3., 4.)
            );

            let error = tensor
                .__setitem__(
                    PyInt::new(py, 0).into_any(),
                    PyFloat::new(py, 5.).into_any(),
                )
                .unwrap_err();
            assert!(error.to_string().contains("coefficient kind must match"));
        });
    }

    #[test]
    fn tensor_public_surface_has_runtime_docstrings() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let tensor_type = py.get_type::<Spensor>();
            for name in [
                "structure",
                "expression",
                "with_name",
                "sparse",
                "dense",
                "to_dense",
                "to_sparse",
                "format_tensor",
                "to_typst",
                "to_html",
                "formatted",
                "evaluator",
                "scalar",
                "index",
                "outer",
                "contract_ports",
                "compose",
                "dot",
                "trace",
                "__getitem__",
                "__setitem__",
            ] {
                let documentation = tensor_type
                    .getattr(name)?
                    .getattr("__doc__")?
                    .extract::<Option<String>>()?;
                assert!(
                    documentation.is_some_and(|documentation| !documentation.trim().is_empty()),
                    "Tensor.{name} is missing its Python docstring"
                );
            }
            for (name, contract) in [
                ("dense", "logical row-major order"),
                ("format_tensor", "logical interface order"),
                ("to_typst", "logical interface order"),
            ] {
                let documentation = tensor_type
                    .getattr(name)?
                    .getattr("__doc__")?
                    .extract::<String>()?;
                assert!(
                    documentation.contains(contract),
                    "Tensor.{name} is missing its logical-order contract"
                );
            }
            // CPython supplies generic slot documentation for special methods, so their
            // detailed contracts cannot be asserted through `__doc__` at runtime.
            Ok(())
        })
        .unwrap();
    }

    #[cfg(feature = "python_stubgen")]
    #[test]
    fn tensor_indexing_stubs_document_logical_order() {
        let methods = pyo3_stub_gen::inventory::iter::<PyMethodsInfo>
            .into_iter()
            .filter(|info| (info.struct_id)() == std::any::TypeId::of::<Spensor>())
            .flat_map(|info| info.methods);
        let getitem = methods
            .clone()
            .filter(|method| method.name == "__getitem__" && method.is_overload)
            .collect::<Vec<_>>();
        assert_eq!(getitem.len(), 2);
        assert!(getitem.iter().all(|method| {
            method.doc.contains("logical row-major")
                && method
                    .doc
                    .contains("canonical storage-axis order is not exposed")
        }));

        let setitem = methods
            .filter(|method| method.name == "__setitem__" && method.is_overload)
            .collect::<Vec<_>>();
        assert_eq!(setitem.len(), 1);
        assert!(setitem[0].doc.contains("logical row-major"));
        assert!(
            setitem[0]
                .doc
                .contains("canonical storage-axis order is not exposed")
        );
        assert!(setitem[0].doc.contains("complex for complex storage"));
    }
}
