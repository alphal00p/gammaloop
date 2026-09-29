use std::{cell::Cell, collections::HashMap, ops::Deref};

use eyre::eyre;

use broadcast::SpensoBroadcastFunction;
use library::{SpensorFunctionLibrary, SpensorLibrary};
use network::{ConvertibleToSpensoNet, SpensoNet};

use pyo3::{
    PyClass,
    exceptions::{self, PyIndexError, PyOverflowError, PyRuntimeError, PyTypeError},
    prelude::*,
    types::{PyComplex, PyDict, PyFloat, PySlice, PyTuple, PyType},
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

use symbolica::api::python::{PythonExpression, PythonFormattedOutput};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, TypeInfo, derive::*, impl_stub_type};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass_enum, gen_stub_pyfunction};

pub mod aliases;
pub mod broadcast;
mod data;
pub mod display;
pub mod expression;
pub mod library;
pub mod metadata;
pub mod network;
pub mod pattern;
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

/// Policy for Rayon operations that manipulate Symbolica expressions.
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

/// Configure whether Spenso may use Rayon for Symbolica operations.
///
/// `Auto` checks the Symbolica license once during this call. When licensed,
/// operations may use Rayon and apply workload heuristics where available;
/// when unlicensed, symbolic operations remain serial. The returned boolean
/// reports whether the resolved policy permits Rayon, not whether every
/// operation will use it.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> sp.set_symbolica_rayon_enabled(sp.SymbolicParallelism.Serial)
/// >>> sp.set_symbolica_rayon_enabled(sp.SymbolicParallelism.Auto)
///
/// Parameters
/// ----------
/// policy : SymbolicParallelism
///     SymbolicParallelism.Auto uses licensing and workload checks; Serial stays on the calling thread; Parallel explicitly requests Rayon execution.
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
        aliases::AliasedTensorExpression,
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
    ],
    registered_functions: [
        set_symbolica_rayon_enabled => "set_symbolica_rayon_enabled",
    ],
    registered_modules: [
        simplification => [
            "CanonicalizationError", "ColorSimplifySettings",
            "CookMode", "CookSettings", "CookSourceFilter", "CookTagFilter", "CookingError",
            "DiracAdjointError", "GammaChainOrdering",
            "GammaSimplifySettings", "NetworkToolingError",
            "AlgebraSettings", "ContractSettings", "ReductionStatus",
        ],
        display => [
            "DisplaySettings", "format_tensor", "to_typst", "to_html", "to_svg", "formatted",
            "load_math_font",
        ],
        pattern => ["PortPattern", "TensorPattern"],
        expression => [
            "TensorExpression", "_AutoIndex", "AUTO", "_",
            "as_tensor", "dot", "chain", "trace",
        ],
    ],
    returned_opaque_classes: [],
}

/// A tensor class that can be either dense or sparse with flexible data types.
///
/// The tensor can store data as floats, complex numbers, or symbolic expressions.
/// Tensors have an associated structure that defines their shape and index properties.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::A")(r, r)
/// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
/// >>> assert tensor[1, 0] == 3.0
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
        if let Ok(values) = value.extract::<Vec<f64>>() {
            Ok(Self::Floats(values))
        } else if let Ok(values) = value.extract::<Vec<Complex<f64>>>() {
            Ok(Self::Complex(values))
        } else if let Ok(values) = value.extract::<Vec<PythonExpression>>() {
            Ok(Self::Atoms(
                values.into_iter().map(|value| value.expr).collect(),
            ))
        } else {
            Err(PyTypeError::new_err(
                "data must be a list of floats, complex numbers, or Expressions",
            ))
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl_stub_type!(AtomsOrFloats = Vec<PythonExpression> | Vec<f64> | Vec<Complex<f64>>);

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
        value: f64,
        descriptor: SymbolicTensor<PartialStructure>,
        descriptor_name: Option<Symbol>,
        descriptor_args: Vec<Atom>,
    ) -> Self {
        Self::from_storage_with_descriptor(
            ParamOrConcrete::new_scalar(ConcreteOrParam::Concrete(RealOrComplex::Real(value))),
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

#[cfg(feature = "python_stubgen")]
impl PyStubType for TensorElements {
    fn type_input() -> pyo3_stub_gen::TypeInfo {
        PythonExpression::type_input() | Complex::type_input() | PyFloat::type_input()
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
    /// Return the symbolic descriptor of this tensor, independently of its component data.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> descriptor = tensor.expression()
    /// >>> assert descriptor.shape == (2, 2)
    pub fn expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_known_parts(
            py,
            self.descriptor.expression().clone(),
            self.descriptor.structure().clone(),
            self.descriptor_name,
            self.descriptor_args.clone(),
        )
    }

    /// Tensor identity, scalar arguments, and external ports in logical order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor.structure
    #[getter]
    fn structure(&self) -> metadata::SpensoTensorStructure {
        metadata::SpensoTensorStructure {
            interface: self.descriptor.structure().clone(),
            name: self.descriptor_name,
            arguments: self.descriptor_args.clone(),
        }
    }

    /// Return a copy with a new data identity while preserving its expression and interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> renamed = tensor.with_name("docs::B")
    ///
    /// Parameters
    /// ----------
    /// name : TensorName | str | Expression | TensorExpression
    ///     New data identity, as a TensorName, name string, symbol, or named TensorExpression. The tensor’s mathematical expression and logical interface are retained.
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

    /// Create a tensor with zero-valued components and sparse storage.
    /// The named structure determines the logical axes. Assign nonzero components with
    /// `tensor[index] = value`; their coefficient type must match the chosen storage class.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.sparse(A, float)
    /// >>> tensor[1, 1] = 5.0
    ///
    /// Parameters
    /// ----------
    /// structure : TensorExpression
    ///     Named TensorExpression defining the component dimensions, scalar arguments, and logical axis order.
    /// type_info : type
    ///     Storage coefficient class: float or symbolica.Expression. All sparse entries initially contain zero.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A new sparse tensor with all elements initially zero
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

    /// Create a tensor from a flat list of components in logical row-major order.
    /// The last axis varies fastest. The named structure fixes the complete shape;
    /// the number of values must equal the product of its dimensions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> assert tensor[1, 0] == 3.0
    ///
    /// Parameters
    /// ----------
    /// structure : TensorExpression
    ///     Named TensorExpression defining the component dimensions, scalar arguments, and logical axis order.
    /// data : Sequence[Expression] | Sequence[float] | Sequence[complex]
    ///     Flat sequence of coefficients in logical row-major order: the last axis varies fastest. Its length must equal the product of the axis dimensions.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A new dense tensor with the specified data
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
    /// Convert this tensor to dense storage format.
    ///
    /// Returns an independent dense tensor; the original storage is unchanged.
    /// This allocates memory for all tensor elements.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> converted = tensor.to_dense()
    /// >>> assert converted[:] == tensor[:]
    #[allow(clippy::wrong_self_convention)]
    fn to_dense(&self) -> Self {
        let mut result = self.clone();
        result.tensor = result.tensor.to_dense();
        result
    }

    /// Convert this tensor to sparse storage format.
    ///
    /// Returns an independent sparse tensor, only storing non-zero elements.
    /// This can save memory for tensors with many zero elements.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> converted = tensor.to_sparse()
    /// >>> assert converted[:] == tensor[:]
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

    /// Format this concrete tensor using compact Spenso notation.
    ///
    /// Values are read in logical interface order without transposing or densifying
    /// the stored tensor. Large results show a bounded edge preview.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> text = tensor.format_tensor(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
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

    /// Format this concrete tensor as Typst math source.
    ///
    /// Values are read in logical interface order without transposing or densifying
    /// the stored tensor. Large results show a bounded edge preview.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> typst_source = tensor.to_typst(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
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

    /// Build Symbolica's rich display wrapper for this concrete tensor.
    ///
    /// Its text, HTML, and LaTeX representations all use logical interface order.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> display_value = tensor.formatted(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    /// notation_source : str, optional
    ///     Optional custom Typst notation definitions. None uses the bundled tensor notation. Pass trusted Typst source, not untrusted user text.
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

    /// Display concrete components in a self-contained interactive HTML explorer.
    ///
    /// Select matrix slices in the view control, or pass
    /// `DisplaySettings(tensor_view="matrix")` for static mathematical HTML.
    /// Uses the optional Typst renderer for symbolic component formulas.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> html = tensor.to_html(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    /// notation_source : str, optional
    ///     Optional custom Typst notation definitions. None uses the bundled tensor notation. Pass trusted Typst source, not untrusted user text.
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

    /// Compile this concrete tensor to SVG with the optional Typst renderer.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> svg = tensor.to_svg(settings=sp.DisplaySettings(index_style="raw", show_dimensions=True))
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Override whether representation dimensions are displayed. None uses settings; the default display settings hide dimensions.
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    /// notation_source : str, optional
    ///     Optional custom Typst notation definitions. None uses the bundled tensor notation. Pass trusted Typst source, not untrusted user text.
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
        display::format_concrete_tensor_output(self, false).latex
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

    /// Return the total number of component positions, including implicit sparse zeros.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> size = len(tensor)
    fn __len__(&self) -> usize {
        self.size().unwrap()
    }

    /// Return tensor data in logical interface order.
    ///
    /// An integer is a flat logical row-major position. A list supplies one coordinate
    /// per slot in `structure.slots`; slices likewise traverse flat logical order.
    /// Canonical storage-axis order is never exposed through this API.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> first = tensor[0]
    /// >>> last = tensor[-1]
    /// >>> column = tensor[:, 1]
    /// >>> assert column == [2.0, 4.0]
    ///
    /// Parameters
    /// ----------
    /// item : slice
    ///     Flat position, slice, or coordinate selection. Negative indices count from the end. Coordinate tuples follow the logical axis order.
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

    /// Set tensor element(s) at the specified index or indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor[0] = 9.0
    /// >>> tensor[1, 1] = 5.0
    /// >>> assert tensor[1, 1] == 5.0
    ///
    /// Parameters
    /// ----------
    /// item : int | Sequence[int]
    ///     Flat integer position or coordinate list/tuple for one component. Negative coordinates count from the end.
    /// value : Expression | complex | float
    ///     New component value: float for real storage, complex for complex storage, or Expression for parametric storage.
    #[gen_stub(skip)]
    fn __setitem__<'py>(
        &mut self,
        item: Bound<'py, PyAny>,
        value: Bound<'py, PyAny>,
    ) -> eyre::Result<()> {
        let value = if let Ok(v) = value.extract::<PythonExpression>() {
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

    /// Create an optimized evaluator for symbolic tensor expressions.
    ///
    /// Compiles the symbolic expressions in this tensor into an optimized evaluation tree
    /// that can efficiently compute numerical values for different parameter inputs.
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
    /// >>> results = evaluator.evaluate([[2.0], [3.0]])
    ///
    /// Parameters
    /// ----------
    /// constants : Mapping[Expression, Expression]
    ///     Exact scalar substitutions to apply before constructing the numerical evaluator.
    /// funs : Mapping[tuple[Expression, str, Sequence[Expression]], Expression]
    ///     Definitions for symbolic functions, keyed by (function_symbol, function_name, argument_list), as accepted by Symbolica’s evaluator.
    /// params : Sequence[Expression]
    ///     Symbolica variables in the exact order expected in each numerical input row. Include every free numerical parameter that remains after substitutions.
    /// iterations : int
    ///     Number of optimization iterations when constructing the evaluator. More iterations can improve repeated evaluations at a higher setup cost.
    /// n_cores : int
    ///     Worker count used while optimizing the evaluator. Use 1 for serial setup.
    /// verbose : bool
    ///     Print evaluator-optimization progress when True.
    ///
    /// Returns
    /// -------
    /// TensorEvaluator
    ///     An optimized evaluator for efficient numerical evaluation
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

    /// Extract the scalar value from a rank-0 (scalar) tensor.
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
    ///     The scalar expression contained in this tensor
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     If the tensor is not a scalar
    fn scalar(&self) -> PyResult<PythonExpression> {
        self.tensor
            .clone()
            .scalar()
            .map(|r| PythonExpression { expr: r.into() })
            .ok_or_else(|| PyRuntimeError::new_err("No scalar found"))
    }

    /// Reference this tensor with new abstract indices and return a lazy network.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = tensor.index("i", "j")
    ///
    /// Parameters
    /// ----------
    /// indices : index label, Slot, or AUTO
    ///     One index per selected port, in logical axis order. Accepts integers, strings, expressions, typed Slots, and AUTO; AUTO leaves that port unchanged.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (*indices, cook_indices = None))]
    fn index(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor_reference(self.clone())?.index_network(indices, cook_indices)
    }

    /// Assign all external indices, producing a network if any contractions are requested.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = tensor("i", "j")
    /// >>> renamed = indexed.reindex("mu", sp.AUTO)
    ///
    /// Parameters
    /// ----------
    /// indices : index label, Slot, or AUTO
    ///     One index per selected port, in logical axis order. Accepts integers, strings, expressions, typed Slots, and AUTO; AUTO leaves that port unchanged.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (*indices, cook_indices=None))]
    fn reindex(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor_reference(self.clone())?.reindex(indices, cook_indices)
    }

    /// Rename external indices without changing component values or rank.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = tensor("i", "j")
    /// >>> swapped = indexed.rename_indices({"i": "j", "j": "i"})
    ///
    /// Parameters
    /// ----------
    /// mapping : dict[int | str | Expression | Slot, index label, Slot, or AUTO]
    ///     Old-to-new external index labels, applied simultaneously. Typed Slot keys disambiguate representations. Renaming cannot contract ports or capture dummy indices; use reindex() for contractions.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (mapping, *, cook_indices=None))]
    fn rename_indices(
        &self,
        mapping: &Bound<'_, PyDict>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<Self> {
        SpensoNet::from_tensor_reference(self.clone())?
            .rename_indices(mapping, cook_indices)?
            .result_tensor(None)
    }

    /// Return an independent tensor with permuted logical axes, preserving sparse storage.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> transposed = tensor.permute_axes([1, 0])
    ///
    /// Parameters
    /// ----------
    /// axes : Sequence[int]
    ///     Permutation of all axis positions, starting at zero. For a matrix, [1, 0] transposes the logical axes.
    fn permute_axes(&self, axes: Vec<usize>) -> PyResult<Self> {
        SpensoNet::from_tensor_reference(self.clone())?
            .permute_axes(axes)?
            .result_tensor(None)
    }

    /// Copy the value and its tensor metadata without mutating the original.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> import copy
    /// >>> independent = copy.copy(tensor)
    fn __copy__(&self) -> Self {
        self.clone()
    }

    /// Reference this tensor with new abstract indices and return a lazy network.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> indexed = tensor("i", "j")
    ///
    /// Parameters
    /// ----------
    /// indices : index label, Slot, or AUTO
    ///     One index per selected port, in logical axis order. Accepts integers, strings, expressions, typed Slots, and AUTO; AUTO leaves that port unchanged.
    /// cook_indices : CookSettings, optional
    ///     None leaves index payloads unchanged. CookSettings selects encoding, filters, and output tags; CookSettings.indices() flattens nested index payloads.
    #[pyo3(signature = (*indices, cook_indices = None))]
    fn __call__(
        &self,
        indices: &Bound<'_, PyTuple>,
        cook_indices: Option<&crate::simplification::PyCookSettings>,
    ) -> PyResult<SpensoNet> {
        self.index(indices, cook_indices)
    }

    /// Negate every tensor coefficient, preserving the external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> negative = -tensor
    fn __neg__(&self) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__neg__()
    }

    /// Add compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor + tensor
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    fn __add__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.add_network(rhs.to_net(), false)
    }

    /// Add compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor + tensor
    ///
    /// Parameters
    /// ----------
    /// lhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    fn __radd__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .add_network(SpensoNet::from_tensor(self.clone())?, false)
    }

    /// Subtract compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor - tensor
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    fn __sub__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.add_network(rhs.to_net(), true)
    }

    /// Subtract compatible tensors, preserving their logical external axes.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor - tensor
    ///
    /// Parameters
    /// ----------
    /// lhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    fn __rsub__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .add_network(SpensoNet::from_tensor(self.clone())?, true)
    }

    /// Multiply tensor operands, contracting repeated compatible explicit indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor * 2
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    fn __mul__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.multiply_network(rhs.to_net())
    }

    /// Multiply tensor operands, contracting repeated compatible explicit indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = 2 * tensor
    ///
    /// Parameters
    /// ----------
    /// lhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    fn __rmul__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .multiply_network(SpensoNet::from_tensor(self.clone())?)
    }

    /// Divide tensor coefficients by a scalar denominator.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor / 2
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Right operand, with a compatible tensor interface for addition or subtraction. A denominator must be scalar.
    fn __truediv__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__truediv__(rhs)
    }

    /// Divide a scalar by a rank-zero tensor, returning a TensorNetwork.
    /// Call `.to_tensor().scalar()` to obtain its component value.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> scalar = sp.Tensor.dense(sp.TensorName("docs::scalar")(), [2.0])
    /// >>> reciprocal = 1 / scalar
    ///
    /// Parameters
    /// ----------
    /// lhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Left operand. Addition and subtraction require compatible tensor interfaces.
    fn __rtruediv__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__rtruediv__(lhs)
    }

    /// Form a lazy outer product without contracting compatible ports.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor.outer(tensor)
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
    fn outer(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.outer(rhs)
    }

    /// Contract one selected pair of ordered interface positions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor.contract_ports(tensor, left=1, right=0)
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
    /// left : int
    ///     Zero-based axis position in this operand to contract.
    /// right : int
    ///     Zero-based axis position in rhs to contract with the selected left axis.
    #[pyo3(signature = (rhs, *, left, right))]
    fn contract_ports(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: usize,
        right: usize,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.contract_ports(rhs, left, right)
    }

    /// Compose two selected `(input, output)` matrix channels.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> result = tensor.compose(tensor, left=(0, 1), right=(0, 1))
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
    /// left : tuple[int, int]
    ///     (input_axis, output_axis) channel in this operand.
    /// right : tuple[int, int]
    ///     (input_axis, output_axis) channel in rhs. The selected connecting ports must have compatible representations.
    #[pyo3(signature = (rhs, *, left, right))]
    fn compose(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: (usize, usize),
        right: (usize, usize),
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.compose(rhs, left, right)
    }

    /// Contract two rank-one operands into the canonical dot form.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> vector = sp.Tensor.dense(sp.TensorName.vector("docs::p")(r), [1.0, 2.0])
    /// >>> norm_squared = vector.dot(vector).to_tensor().scalar()
    ///
    /// Parameters
    /// ----------
    /// rhs : Expression or scalar | TensorExpression | TensorNetwork | Tensor
    ///     Other tensor, network, or scalar operand. Concrete component data or a network makes the result a TensorNetwork.
    fn dot(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.dot(rhs)
    }

    /// Close `channel`, or the unique matrix channel when it is omitted.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> trace = tensor.trace(channel=(0, 1))
    ///
    /// Parameters
    /// ----------
    /// channel : tuple[int, int], optional
    ///     Pair (input_axis, output_axis) identifying the matrix channel to close. None selects the unique compatible channel; ambiguous tensors require an explicit pair.
    #[pyo3(signature = (*, channel = None))]
    fn trace(&self, channel: Option<(usize, usize)>) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.trace(channel)
    }
}

/// A reusable numerical evaluator for tensors with symbolic components.
/// Create one with `Tensor.evaluator()`. Input rows follow `.parameters`; evaluation
/// returns one Tensor per row, preserving the original logical axes.
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
/// >>> results = evaluator.evaluate([[2.0], [3.0]])
/// >>> assert len(results) == 2
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
    /// Parameters in the exact order expected by each evaluation input row.
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
    /// >>> result = evaluator.parameters
    #[getter]
    fn parameters(&self) -> Vec<PythonExpression> {
        self.parameters.iter().cloned().map(Into::into).collect()
    }

    /// Number of parameter values required in each input row.
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
    /// >>> result = evaluator.input_size
    #[getter]
    fn input_size(&self) -> usize {
        self.parameters.len()
    }

    /// Shape of each returned tensor in logical axis order.
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
    /// >>> result = evaluator.output_shape
    #[getter]
    fn output_shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            tensor_data_layout(self.descriptor.structure())?.logical_shape(),
        )
    }

    /// Whether the exact coefficients admit a real-valued evaluator.
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
    /// >>> result = evaluator.supports_real
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

    /// Evaluate a batch of real parameter rows and return one Tensor per row.
    /// Use `.evaluate_complex()` if the expression has complex coefficients, even when
    /// the supplied parameter values are real.
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
    /// >>> results = evaluator.evaluate([[2.0], [3.0]])
    /// >>> assert len(results) == 2
    ///
    /// Parameters
    /// ----------
    /// inputs : Sequence[Sequence[float]]
    ///     Batch of parameter rows, each in the same order as .parameters. Every row must have exactly .input_size entries; an empty batch returns an empty list.
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     List of evaluated tensors, one for each input parameter set
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the evaluator contains complex coefficients
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

    /// Evaluate a batch of complex parameter rows and return one Tensor per row.
    /// Real parameter values are accepted too; results use complex component storage.
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
    /// >>> results = evaluator.evaluate_complex([[2.0 + 1.0j], [3.0 + 0.0j]])
    /// >>> assert len(results) == 2
    ///
    /// Parameters
    /// ----------
    /// inputs : Sequence[Sequence[complex]]
    ///     Batch of parameter rows, each in the same order as .parameters. Every row must have exactly .input_size entries; an empty batch returns an empty list.
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
    /// Generate C++ code and compile a numerical evaluator for repeated calls.
    /// A C++ compiler must be installed. The returned evaluator accepts the same parameter
    /// order and batches as this evaluator; keep the compiled library file available.
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
    ///
    /// Parameters
    /// ----------
    /// function_name : str
    ///     Name of the generated C++ evaluation function; use a valid C++ identifier.
    /// filename : str
    ///     Path at which to write the generated C++ source file.
    /// library_name : str
    ///     Output path for the compiled shared library. Keep the file available while you use the compiled evaluator.
    /// inline_asm : str
    ///     Assembly backend: "default", "x64", "aarch64", or "none". Use "none" for portable C++ without inline assembly.
    /// optimization_level : int
    ///     C++ compiler optimization level from 0 to 3. Higher levels take longer to compile and can improve repeated evaluation speed.
    /// compiler_path : str, optional
    ///     C++ compiler executable. None lets the evaluator use its default compiler; the compiler must be available on PATH.
    /// custom_header : str, optional
    ///     Additional C++ declarations to include in the generated source, or None.
    ///
    /// Returns
    /// -------
    /// CompiledTensorEvaluator
    ///     A compiled evaluator for maximum performance evaluation
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

/// A numerical tensor evaluator backed by a compiled C++ kernel.
/// Create one with `TensorEvaluator.compile()`. Real coefficients support real and
/// complex inputs; complex coefficients require `.evaluate_complex()`.
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
    /// Parameters in the exact order expected by each evaluation input row.
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
    /// ...     value = compiled.parameters
    #[getter]
    fn parameters(&self) -> Vec<PythonExpression> {
        self.parameters.iter().cloned().map(Into::into).collect()
    }

    /// Number of parameter values required in each input row.
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
    /// ...     value = compiled.input_size
    #[getter]
    fn input_size(&self) -> usize {
        self.parameters.len()
    }

    /// Shape of each returned tensor in logical axis order.
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
    /// ...     value = compiled.output_shape
    #[getter]
    fn output_shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            tensor_data_layout(self.descriptor.structure())?.logical_shape(),
        )
    }

    /// Whether a real-valued native kernel was compiled.
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
    /// ...     value = compiled.supports_real
    #[getter]
    fn supports_real(&self) -> bool {
        self.real.is_some()
    }

    /// Evaluate real parameter rows using the compiled C++ kernel.
    /// Returns one Tensor per row. Use `.evaluate_complex()` when the evaluator contains
    /// complex coefficients.
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
    ///
    /// Parameters
    /// ----------
    /// inputs : Sequence[Sequence[float]]
    ///     Batch of parameter rows, each in the same order as .parameters. Every row must have exactly .input_size entries; an empty batch returns an empty list.
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

    /// Evaluate complex parameter rows using the compiled C++ kernel.
    /// Returns one Tensor per row, with the same logical axes as the original expression.
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
    /// ...     results = compiled.evaluate_complex([[2.0 + 1.0j], [3.0 + 0.0j]])
    ///
    /// Parameters
    /// ----------
    /// inputs : Sequence[Sequence[complex]]
    ///     Batch of parameter rows, each in the same order as .parameters. Every row must have exactly .input_size entries; an empty batch returns an empty list.
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     List of evaluated tensors, one for each input parameter set
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
                doc:r##"Create a slot from this representation, by specifying an index.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.mink(4)
>>> mu = r("mu")
>>> nu = r(1)

Parameters
----------
aind : int | str
    Index label as an integer, string, or Symbolica expression. Equal compatible labels represent contracted indices.

Returns
-------
Slot
    A new Slot object with the specified index"##,
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
                doc:r##"Create a slot from this representation, by specifying an index.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.mink(4)
>>> mu = r("mu")
>>> nu = r(1)

Parameters
----------
aind : int | str
    Index label as an integer, string, or Symbolica expression. Equal compatible labels represent contracted indices.

Returns
-------
Slot
    A new Slot object with the specified index"##,
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
                doc:r##"Iterate over every component in logical row-major order, including sparse zeros.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> entries = list(tensor)"##,
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
                doc:r##"Get tensor elements at the specified range of indices.

Slices traverse the flat logical row-major order shown by `Tensor.structure.slots`;
canonical storage-axis order is not exposed.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> first = tensor[0]
>>> last = tensor[-1]
>>> column = tensor[:, 1]
>>> assert column == [2.0, 4.0]

Parameters
----------
item : slice
    Flat position, slice, or coordinate selection. Negative indices count from the end. Coordinate tuples follow the logical axis order.

Returns
-------
list of float, complex, or Expression
    The tensor elements at the specified range"##,
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
                doc:r##"Get tensor elements at the specified range of indices.

Slices traverse the flat logical row-major order shown by `Tensor.structure.slots`;
canonical storage-axis order is not exposed.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> first = tensor[0]
>>> last = tensor[-1]
>>> column = tensor[:, 1]
>>> assert column == [2.0, 4.0]

Parameters
----------
item : slice
    Flat position, slice, or coordinate selection. Negative indices count from the end. Coordinate tuples follow the logical axis order.

Returns
-------
list of float, complex, or Expression
    The tensor elements at the specified range"##,
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
                doc:r##"Set tensor element at the specified index.

Integers are flat logical row-major positions. Coordinate lists follow
`Tensor.structure.slots`; canonical storage-axis order is not exposed.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> A = sp.TensorName("docs::A")(r, r)
>>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> tensor[0] = 9.0
>>> tensor[1, 1] = 5.0
>>> assert tensor[1, 1] == 5.0

Parameters
----------
item : int | Sequence[int]
    Flat integer position or coordinate list/tuple for one component. Negative coordinates count from the end.
value : Expression | complex | float
    New component value: float for real storage, complex for complex storage, or Expression for parametric storage."##,
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
                if let Some(value) = line.strip_prefix("_LibraryDefault: typing.TypeVar = ") {
                    return format!("_LibraryDefault = {value}");
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
