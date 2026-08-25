use std::{
    collections::HashMap,
    ops::Deref,
    sync::atomic::{AtomicUsize, Ordering},
};

use eyre::eyre;

use broadcast::SpensoBroadcastFunction;
use library::{SpensorFunctionLibrary, SpensorLibrary};
use network::{ConvertibleToSpensoNet, SpensoNet};

use pyo3::{
    PyClass,
    exceptions::{self, PyIndexError, PyOverflowError, PyRuntimeError, PyTypeError},
    prelude::*,
    types::{PyComplex, PyFloat, PySlice, PyTuple, PyType},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};

use spenso::{
    algebra::complex::{Complex, RealOrComplex, symbolica_traits::CompiledComplexEvaluatorSpenso},
    tensors::{
        data::{DenseTensor, GetTensorData, SetTensorData, SparseOrDense, SparseTensor},
        parametric::{
            ConcreteOrParam, EvalTensor, ParamOrConcrete, ParamTensor, atomcore::TensorAtomOps,
        },
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
    api::python::SymbolicaCommunityModule, domains::float::Complex as SymComplex, prelude::*,
};

use symbolica::api::python::{PythonExpression, PythonFormattedOutput};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, TypeInfo, define_stub_info_gatherer, derive::*, impl_stub_type};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass_enum, gen_stub_pyfunction};

pub mod broadcast;
mod composition;
pub mod display;
pub mod expression;
pub mod library;
pub mod network;
pub mod structure;

use composition::StructuredAtom;
use expression::TensorExpression;
use structure::ConvertibleToSpensoName;

static OPEN_OWNER: AtomicUsize = AtomicUsize::new(0);

pub(crate) fn fresh_open_owner() -> usize {
    OPEN_OWNER.fetch_add(1, Ordering::Relaxed)
}

trait ModuleInit: PyClass {
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<Self>()
    }
}

pub struct SpensoModule;

/// Policy for Rayon operations that manipulate Symbolica expressions.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(from_py_object, eq, eq_int, module = "symbolica.community.spenso")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymbolicParallelism {
    /// Permit Rayon when licensed and use workload heuristics where available.
    Auto,
    /// Keep symbolic operations on the calling thread.
    Serial,
    /// Force Rayon without `Auto`'s Symbolica license safety check.
    Parallel,
}

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
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
fn set_symbolica_rayon_enabled(policy: SymbolicParallelism) -> bool {
    spenso::symbolic_parallelism::set_symbolica_rayon_enabled(policy.into());
    spenso::symbolic_parallelism::symbolica_rayon_enabled()
}

impl SymbolicaCommunityModule for SpensoModule {
    fn get_name() -> String {
        "spenso".to_string()
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

pub(crate) fn initialize_spenso(m: &Bound<'_, PyModule>) -> PyResult<()> {
    use network::ExecutionMode;

    // m.add_function(?)?;
    SpensoNet::init(m)?;
    ExecutionMode::init(m)?;
    m.add_class::<SymbolicParallelism>()?;
    m.add_class::<SpensoExpressionEvaluator>()?;
    m.add_class::<SpensoCompiledExpressionEvaluator>()?;
    m.add_function(wrap_pyfunction!(set_symbolica_rayon_enabled, m)?)?;
    display::register(m)?;
    expression::register(m)?;
    Spensor::init(m)?;
    m.add_class::<structure::SpensoName>()?;
    m.add_class::<structure::SpensoSlot>()?;
    m.add_class::<structure::SpensoRepresentation>()?;
    SpensorLibrary::init(m)?;
    SpensorFunctionLibrary::init(m)?;
    SpensoBroadcastFunction::init(m)?;
    let exports = m
        .dict()
        .keys()
        .iter()
        .filter_map(|key| key.extract::<String>().ok())
        .filter(|name| {
            name != "initialize"
                && name != "initialize_module"
                && (name == "_" || !name.starts_with('_'))
        })
        .collect::<Vec<_>>();
    m.add("__all__", exports)?;
    Ok(())
}

/// A tensor class that can be either dense or sparse with flexible data types.
///
/// The tensor can store data as floats, complex numbers, or symbolic expressions.
/// Tensors have an associated structure that defines their shape and index properties.
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
/// >>> rep = Representation.euc(4)
/// >>> structure = TensorName.vector("v")(rep("mu"))
/// >>> data = [1.0, 2.0, 3.0, 4.0]
/// >>> tensor = Tensor.dense(structure, data)
/// >>> sparse_tensor = Tensor.sparse(structure, float)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(from_py_object, name = "Tensor", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct Spensor {
    pub(crate) tensor: MixedTensor<f64, ShadowedStructure<AbstractIndex>>,
    pub(crate) descriptor: StructuredAtom,
    pub(crate) descriptor_name: Option<Symbol>,
    pub(crate) descriptor_args: Vec<Atom>,
}

pub struct TensorDataDescriptor {
    structure: Canonicalized<ShadowedStructure<AbstractIndex>>,
    descriptor: StructuredAtom,
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
        let layout = tensor_data_layout(&descriptor.interface)?;
        let owner = fresh_open_owner();
        let logical_axes = (0..layout.logical_shape().len()).collect::<Vec<_>>();
        let storage_axes = descriptor
            .interface
            .layout()
            .logical_to_canonical(&logical_axes);
        let mut storage_indices = vec![0; storage_axes.len()];
        for (storage_index, &logical_axis) in storage_axes.iter().enumerate() {
            storage_indices[logical_axis] = storage_index;
        }
        let structure = OrderedStructure::new(
            descriptor
                .interface
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
        descriptor: StructuredAtom,
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
        descriptor: StructuredAtom,
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
#[pymethods]
impl Spensor {
    /// Return the exact structured expression describing this tensor's data.
    pub fn structure(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface_descriptor(
            py,
            self.descriptor.atom.clone(),
            self.descriptor.interface.clone(),
            self.descriptor_name,
            self.descriptor_args.clone(),
        )
    }

    /// Return a copy with a new data identity while preserving its expression and interface.
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

    #[staticmethod]
    /// Create a new sparse empty tensor with the given structure and data type.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorExpression
    ///     A named expression with unresolved, explicit, or mixed interface ports
    /// type_info : type
    ///     The data type - either `float` or `Expression` class
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A new sparse tensor with all elements initially zero
    ///
    /// Examples
    /// --------
    /// >>> import symbolica
    /// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
    /// >>> rep = Representation.euc(3)
    /// >>> structure = TensorName("T")(rep, rep)
    /// >>> sparse_float = Tensor.sparse(structure, float)
    /// >>> sparse_sym = Tensor.sparse(structure, symbolica.Expression)
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

    #[staticmethod]
    /// Create a new dense tensor with the given structure and data.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorExpression
    ///     A named expression with unresolved, explicit, or mixed interface ports
    /// data : list of float, complex, or Expression
    ///     Tensor values in logical row-major order: the last interface axis varies
    ///     fastest, regardless of the canonical order used by internal storage.
    ///     Representation-axis movement occurs only at this input boundary; later
    ///     index-only storage permutations do not change the public layout.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A new dense tensor with the specified data
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
    /// >>> rep = Representation.euc(2)
    /// >>> structure = TensorName("T")(rep("mu"), rep("nu"))
    /// >>> data = [1.0, 2.0, 3.0, 4.0]
    /// >>> tensor = Tensor.dense(structure, data)
    /// >>> x, y = S("x", "y")
    /// >>> sym_data = [x, y, x * y, x + y]
    /// >>> sym_tensor = Tensor.dense(structure, sym_data)
    pub fn dense(structure: TensorDataDescriptor, data: AtomsOrFloats) -> PyResult<Spensor> {
        let TensorDataDescriptor {
            structure,
            descriptor,
            name,
            args,
        } = structure;
        let layout = tensor_data_layout(&descriptor.interface)?;
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
    #[allow(clippy::wrong_self_convention)]
    /// Convert this tensor to dense storage format.
    ///
    /// Convert this tensor to dense storage format.
    ///
    /// Converts sparse tensors to dense format in-place. Dense tensors are unchanged.
    /// This allocates memory for all tensor elements.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
    /// >>> rep = Representation.euc(4)
    /// >>> structure = TensorName.vector("v")(rep("mu"))
    /// >>> tensor = Tensor.sparse(structure, float)
    /// >>> tensor[0] = 1.0
    /// >>> tensor.to_dense()
    fn to_dense(&mut self) {
        self.tensor = self.tensor.clone().to_dense();
    }

    #[allow(clippy::wrong_self_convention)]
    /// Convert this tensor to sparse storage format.
    ///
    /// Convert this tensor to sparse storage format.
    ///
    /// Converts dense tensors to sparse format in-place, only storing non-zero elements.
    /// This can save memory for tensors with many zero elements.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
    /// >>> rep = Representation.euc(2)
    /// >>> structure = TensorName("T")(rep("mu"), rep("nu"))
    /// >>> data = [1.0, 0.0, 0.0, 2.0]
    /// >>> tensor = Tensor.dense(structure, data)
    /// >>> tensor.to_sparse()
    fn to_sparse(&mut self) {
        self.tensor = self.tensor.clone().to_sparse();
    }

    fn __repr__(&self) -> String {
        format!("Tensor({})", display::format_concrete_tensor(self, false))
    }

    fn __str__(&self) -> String {
        display::format_concrete_tensor(self, false)
    }

    /// Format this concrete tensor using compact Spenso notation.
    ///
    /// Values are read in logical interface order without transposing or densifying
    /// the stored tensor. Large results show a bounded edge preview.
    #[pyo3(signature = (show_dimensions = false))]
    fn format_tensor(&self, show_dimensions: bool) -> String {
        display::format_concrete_tensor(self, show_dimensions)
    }

    /// Format this concrete tensor as Typst math source.
    ///
    /// Values are read in logical interface order without transposing or densifying
    /// the stored tensor. Large results show a bounded edge preview.
    #[pyo3(signature = (show_dimensions = false))]
    fn to_typst(&self, show_dimensions: bool) -> String {
        display::concrete_tensor_to_typst(self, show_dimensions)
    }

    /// Format this concrete tensor as an escaped HTML table.
    ///
    /// Values are shown in logical interface order. Higher-rank tensors are
    /// displayed as labeled matrix slices, while large sparse tensors use a
    /// coordinate/value table without being densified.
    ///
    /// Parameters
    /// ----------
    /// show_dimensions : bool, optional
    ///     Include representation dimensions in the tensor interface.
    ///
    /// Returns
    /// -------
    /// str
    ///     A self-contained HTML fragment suitable for notebook display.
    #[pyo3(signature = (show_dimensions = false))]
    fn to_html(&self, show_dimensions: bool) -> String {
        display::concrete_tensor_to_html(self, show_dimensions)
    }

    /// Build Symbolica's rich display wrapper for this concrete tensor.
    ///
    /// Its text, HTML, and LaTeX representations all use logical interface order.
    #[pyo3(signature = (show_dimensions = false))]
    fn formatted(&self, show_dimensions: bool) -> PythonFormattedOutput {
        display::format_concrete_tensor_output(self, show_dimensions)
    }

    fn _repr_html_(&self) -> String {
        self.to_html(false)
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

    fn __len__(&self) -> usize {
        self.size().unwrap()
    }

    /// Return tensor data in logical interface order.
    ///
    /// An integer is a flat logical row-major position. A list supplies one coordinate
    /// per slot in `structure().interface`; slices likewise traverse flat logical order.
    /// Canonical storage-axis order is never exposed through this API.
    fn __getitem__(&self, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        let layout = tensor_data_layout(&self.descriptor.interface)?;
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

                let start = if r.start < 0 {
                    (r.slicelength as isize + r.start) as usize
                } else {
                    r.start as usize
                };

                let end = if r.stop < 0 {
                    (r.slicelength as isize + r.stop) as usize
                } else {
                    r.stop as usize
                };

                let (range, step) = if r.step < 0 {
                    (end..start, -r.step as usize)
                } else {
                    (start..end, r.step as usize)
                };

                let slice = range
                    .step_by(step)
                    .map(|i| get_owned_linear(i).map(TensorElements::from))
                    .collect::<PyResult<Vec<_>>>()?;

                return Ok(
                    Python::attach(|py| slice.into_pyobject(py).map(|a| a.unbind()))?.into_any(),
                );
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
    /// Parameters
    /// ----------
    /// item : int or list of int
    ///     Logical index specification (int for flat row-major position, list of int
    ///     for coordinates following `structure().interface`)
    /// value : float, complex, or Expression
    ///     The value to set. Its coefficient kind must match the tensor: `float` for
    ///     real storage, `complex` for complex storage, or `Expression` for parametric storage.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
    /// >>> rep = Representation.euc(2)
    /// >>> structure = TensorName("T")(rep("mu"), rep("nu"))
    /// >>> tensor = Tensor.sparse(structure, float)
    /// >>> tensor[0] = 4.0
    /// >>> tensor[1, 1] = 1.0
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

        let layout = tensor_data_layout(&self.descriptor.interface)?;
        let coefficient_error = |error| {
            eyre!(
                "assigned coefficient kind must match tensor storage (float for real, complex for complex, Expression for parametric): {error}"
            )
        };
        if let Ok(flat_index) = item.extract::<usize>() {
            let canonical = layout
                .logical_flat_to_storage_flat(flat_index)
                .map_err(tensor_data_layout_error)?;
            self.tensor
                .set_flat(canonical.into(), value)
                .map_err(coefficient_error)
        } else if let Ok(expanded_idxs) = item.extract::<Vec<usize>>() {
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

    #[pyo3(signature =
           (constants,
           funs,
           params,
           iterations = 100,
           n_cores = 4,
           verbose = false),
           )]
    /// Create an optimized evaluator for symbolic tensor expressions.
    ///
    /// Create an optimized evaluator for symbolic tensor expressions.
    ///
    /// Compiles the symbolic expressions in this tensor into an optimized evaluation tree
    /// that can efficiently compute numerical values for different parameter inputs.
    ///
    /// Parameters
    /// ----------
    /// constants : dict
    ///     Dict mapping symbolic expressions to their constant numerical values
    /// funs : dict
    ///     Dict mapping function signatures to their symbolic definitions
    /// params : list of Expression
    ///     List of symbolic parameters that will be varied during evaluation
    /// iterations : int, optional
    ///     Number of optimization iterations for Horner scheme (default: 100)
    /// n_cores : int, optional
    ///     Number of CPU cores to use for optimization (default: 4)
    /// verbose : bool, optional
    ///     Whether to print optimization progress (default: False)
    ///
    /// Returns
    /// -------
    /// TensorEvaluator
    ///     An optimized evaluator for efficient numerical evaluation
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
    /// >>> x, y = S("x", "y")
    /// >>> rep = Representation.euc(2)
    /// >>> structure = TensorName.vector("v")(rep("mu"))
    /// >>> tensor = Tensor.dense(structure, [x * y, x + y])
    /// >>> evaluator = tensor.evaluator(constants={}, funs={}, params=[x, y], iterations=50)
    /// >>> results = evaluator.evaluate_complex([[1.0, 2.0], [3.0, 4.0]])
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

        let mut evaltensor = match &self.tensor {
            ParamOrConcrete::Param(s) => s.to_evaluation_tree(&fn_map, &params).map_err(|e| {
                exceptions::PyValueError::new_err(format!("Could not create evaluator: {}", e))
            })?,
            ParamOrConcrete::Concrete(_) => return Err(PyRuntimeError::new_err("not atom")),
        };

        evaltensor.optimize_horner_scheme(&settings);

        evaltensor.common_subexpression_elimination();
        let linear = evaltensor.linearize(&settings);
        Ok(SpensoExpressionEvaluator {
            eval: None,
            eval_complex: linear
                .clone()
                .map_coeff(&|x| Complex::new(x.re.to_f64(), x.im.to_f64())),
            eval_rat: linear,
            descriptor: self.descriptor.clone(),
            descriptor_name: self.descriptor_name,
            descriptor_args: self.descriptor_args.clone(),
        })
    }

    /// Extract the scalar value from a rank-0 (scalar) tensor.
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
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import Tensor, TensorName, Representation
    /// >>> rep = Representation.euc(2)
    /// >>> p = TensorName.vector("p")
    /// >>> structure = (p(rep("mu")) * p(rep("mu"))).with_name(TensorName("p2"))
    /// >>> scalar_tensor = Tensor.dense(structure, [1.0])
    /// >>> value = scalar_tensor.scalar()
    fn scalar(&self) -> PyResult<PythonExpression> {
        self.tensor
            .clone()
            .scalar()
            .map(|r| PythonExpression { expr: r.into() })
            .ok_or_else(|| PyRuntimeError::new_err("No scalar found"))
    }

    /// Reference this tensor with new abstract indices and return a lazy network.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn index(&self, indices: &Bound<'_, PyTuple>, cook_indices: bool) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor_reference(self.clone())?.index_network(indices, cook_indices)
    }

    /// Reference this tensor with new abstract indices and return a lazy network.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn __call__(&self, indices: &Bound<'_, PyTuple>, cook_indices: bool) -> PyResult<SpensoNet> {
        self.index(indices, cook_indices)
    }

    fn __neg__(&self) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__neg__()
    }

    fn __add__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.add_network(rhs.to_net(), false)
    }

    fn __radd__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .add_network(SpensoNet::from_tensor(self.clone())?, false)
    }

    fn __sub__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.add_network(rhs.to_net(), true)
    }

    fn __rsub__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .add_network(SpensoNet::from_tensor(self.clone())?, true)
    }

    fn __mul__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.multiply_network(rhs.to_net())
    }

    fn __rmul__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        lhs.to_net()
            .multiply_network(SpensoNet::from_tensor(self.clone())?)
    }

    fn __truediv__(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__truediv__(rhs)
    }

    fn __rtruediv__(&self, lhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.__rtruediv__(lhs)
    }

    /// Form a lazy outer product without contracting compatible ports.
    fn outer(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.outer(rhs)
    }

    /// Contract one selected pair of ordered interface positions.
    #[pyo3(signature = (rhs, *, left, right))]
    fn contract(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: usize,
        right: usize,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.contract(rhs, left, right)
    }

    /// Compose two selected `(input, output)` matrix channels.
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
    fn dot(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.dot(rhs)
    }

    /// Close `channel`, or the unique matrix channel when it is omitted.
    #[pyo3(signature = (*, channel = None))]
    fn trace(&self, channel: Option<(usize, usize)>) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.trace(channel)
    }
}

/// An optimized evaluator for symbolic tensor expressions.
///
/// An optimized evaluator for symbolic tensor expressions.
///
/// This class provides efficient numerical evaluation of symbolic tensor expressions
/// after optimization. It supports both real and complex-valued evaluations.
///
/// Create instances using the `Tensor.evaluator()` method rather than directly.
///
/// Examples
/// --------
/// >>> evaluator = my_tensor.evaluator(constants={}, funs={}, params=[x, y])
/// >>> results = evaluator.evaluate([[1.0, 2.0], [3.0, 4.0]])
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "TensorEvaluator",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub struct SpensoExpressionEvaluator {
    pub eval_rat: LinearizedEvalTensor<SymComplex<Rational>, ShadowedStructure<AbstractIndex>>,
    pub eval: Option<LinearizedEvalTensor<f64, ShadowedStructure<AbstractIndex>>>,
    pub eval_complex: LinearizedEvalTensor<Complex<f64>, ShadowedStructure<AbstractIndex>>,
    descriptor: StructuredAtom,
    descriptor_name: Option<Symbol>,
    descriptor_args: Vec<Atom>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl SpensoExpressionEvaluator {
    /// Evaluate the tensor expression for multiple real-valued parameter inputs.
    ///
    /// Parameters
    /// ----------
    /// inputs : list of list of float
    ///     List of parameter value lists, where each inner list contains
    ///     numerical values for all parameters in the same order as specified
    ///     when creating the evaluator
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
    ///
    /// Examples
    /// --------
    /// >>> results = evaluator.evaluate([[1.0, 2.0], [3.0, 4.0]])
    fn evaluate(&mut self, inputs: Vec<Vec<f64>>) -> PyResult<Vec<Spensor>> {
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

    /// Evaluate the expression for multiple inputs and return the results.
    fn evaluate_complex(&mut self, inputs: Vec<Vec<Complex<f64>>>) -> Vec<Spensor> {
        let eval = &mut self.eval_complex;

        inputs
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
            .collect()
    }

    /// Compile the evaluator to a shared library using C++ for maximum performance.
    ///
    /// Compile the evaluator to a shared library using C++ for maximum performance.
    ///
    /// Generates optimized C++ code with optional inline assembly and compiles it
    /// into a shared library that can be loaded for extremely fast evaluation.
    ///
    /// Parameters
    /// ----------
    /// function_name : str
    ///     Name for the generated C++ function
    /// filename : str
    ///     Path for the generated C++ source file
    /// library_name : str
    ///     Name for the compiled shared library
    /// inline_asm : str, optional
    ///     Type of inline assembly optimization ("default", "x64", "aarch64", "none")
    /// optimization_level : int, optional
    ///     Compiler optimization level 0-3 (default: 3)
    /// compiler_path : str, optional
    ///     Path to specific C++ compiler (default: system default)
    /// custom_header : str, optional
    ///     Additional C++ header code to include
    ///
    /// Returns
    /// -------
    /// CompiledTensorEvaluator
    ///     A compiled evaluator for maximum performance evaluation
    ///
    /// Examples
    /// --------
    /// >>> compiled = evaluator.compile(
    /// ...     function_name="fast_eval",
    /// ...     filename="tensor_eval.cpp",
    /// ...     library_name="tensor_lib",
    /// ...     optimization_level=3,
    /// ... )
    /// >>> results = compiled.evaluate_complex([[1.0, 2.0], [3.0, 4.0]])
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

        Ok(SpensoCompiledExpressionEvaluator {
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

/// A compiled and optimized evaluator for maximum performance tensor evaluation.
///
/// This class wraps a compiled C++ shared library for extremely fast numerical
/// evaluation of tensor expressions. It only supports complex-valued evaluation
/// as this is the most general case.
///
/// A compiled and optimized evaluator for maximum performance tensor evaluation.
///
/// This class wraps a compiled C++ shared library for extremely fast numerical
/// evaluation of tensor expressions. It only supports complex-valued evaluation
/// as this is the most general case.
///
/// Create instances using the `TensorEvaluator.compile()` method.
///
/// Examples
/// --------
/// >>> compiled = evaluator.compile("eval_func", "code.cpp", "lib")
/// >>> results = compiled.evaluate_complex(large_input_batch)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "CompiledTensorEvaluator",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub struct SpensoCompiledExpressionEvaluator {
    pub eval: EvalTensor<CompiledComplexEvaluatorSpenso, ShadowedStructure<AbstractIndex>>,
    descriptor: StructuredAtom,
    descriptor_name: Option<Symbol>,
    descriptor_args: Vec<Atom>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl SpensoCompiledExpressionEvaluator {
    /// Evaluate the tensor expression for multiple complex-valued parameter inputs.
    ///
    /// Evaluate the tensor expression for multiple complex-valued parameter inputs.
    ///
    /// Uses the compiled C++ code for maximum performance evaluation with complex numbers.
    ///
    /// Parameters
    /// ----------
    /// inputs : list of list of complex
    ///     List of parameter value lists, where each inner list contains
    ///     complex values for all parameters in the same order as specified
    ///     when creating the original evaluator
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     List of evaluated tensors, one for each input parameter set
    ///
    /// Examples
    /// --------
    /// >>> complex_inputs = [
    /// ...     [1.0+2.0j, 3.0+0.0j],
    /// ...     [0.0+1.0j, 2.0+1.0j]
    /// ... ]
    /// >>> results = compiled_evaluator.evaluate_complex(complex_inputs)
    fn evaluate_complex(&mut self, inputs: Vec<Vec<Complex<f64>>>) -> Vec<Spensor> {
        inputs
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
            .collect()
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
                        type_info: structure::ConvertibleToAbstractIndex::type_input,
                    },
                ],
                r#type: MethodType::Instance,
                r#return: structure::SpensoSlot::type_output,
                doc:r##"Create a slot from this representation, by specifying an index.

Parameters
----------
aind : int, str, or Symbol
    The index specification

Returns
-------
Slot
    A new Slot object with the specified index

Examples
--------
>>> from symbolica.community.spenso import Representation
>>> import symbolica as sp
>>> rep = Representation.euc(3)
>>> slot1 = rep('mu')
>>> slot2 = rep(1)
>>> slot3 = rep(sp.S('nu'))
"##,
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
                doc:r##"Create a slot or symbolic expression from this representation.

Parameters
----------
aind : Expression
    The index specification (Expression creates symbolic representation)

Returns
-------
Expression
    A symbolic expression representing this representation

Examples
--------
>>> from symbolica.community.spenso import Representation
>>> import symbolica as sp
>>> rep = Representation.euc(3)
>>> expr = rep(sp.E("cos(x)"))
"##,
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
                    name: "typing.Iterator[typing.Any]".into(),
                    import: std::collections::HashSet::new(),
                },
                doc:r##"Iterator"##,
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

Slices traverse the flat logical row-major order shown by `Tensor.structure().interface`;
canonical storage-axis order is not exposed.

Parameters
----------
item : slice
    Slice object defining the range of indices

Returns
-------
list of float, complex, or Expression
    The tensor elements at the specified range
"##,
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
                doc:r##"Get tensor element at the specified index or indices.

Integers are flat logical row-major positions. Coordinate lists follow
`Tensor.structure().interface`; canonical storage-axis order is not exposed.

Parameters
----------
item : int or list of int
    Index specification (int for flat index, list of int for coordinates)

Returns
-------
float, complex, or Expression
    The tensor element at the specified index
"##,
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
`Tensor.structure().interface`; canonical storage-axis order is not exposed.

Parameters
----------
item : int or list of int
    Index specification (int for flat index, list of int for coordinates)
value : float, complex, or Expression
    The value to set. Use float for real storage, complex for complex storage,
    or Expression for parametric storage.

Examples
--------
>>> from symbolica.community.spenso import Tensor, TensorName, Representation
>>> rep = Representation.euc(2)
>>> structure = TensorName("T")(rep("mu"), rep("nu"))
>>> tensor = Tensor.sparse(structure, float)
>>> tensor[0] = 1.0
>>> tensor[1, 1] = 2.0
"##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
        ]
    }
}

#[cfg(feature = "python_stubgen")]
define_stub_info_gatherer!(stub_info);

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
            StructuredAtom::new(Atom::Zero, interface),
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
                "contract",
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
