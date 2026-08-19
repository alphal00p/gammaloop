use std::{collections::HashMap, ops::Deref};

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
        HasStructure, OrderedStructure, PermutedStructure, ScalarTensor, TensorStructure,
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
    pub(crate) tensor: PermutedStructure<MixedTensor<f64, ShadowedStructure<AbstractIndex>>>,
    pub(crate) descriptor: StructuredAtom,
    pub(crate) descriptor_name: Option<Symbol>,
    pub(crate) descriptor_args: Vec<Atom>,
}

pub struct TensorDataDescriptor {
    structure: PermutedStructure<ShadowedStructure<AbstractIndex>>,
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

/// Maps the public row-major layout in logical interface order to Spenso's
/// canonically ordered storage axes.
pub(crate) struct LogicalDataLayout {
    logical_shape: Vec<usize>,
    logical_strides: Vec<usize>,
    canonical_axes: Vec<usize>,
    canonical_strides: Vec<usize>,
    size: usize,
}

impl LogicalDataLayout {
    pub(crate) fn from_interface(interface: &PartialStructure) -> PyResult<Self> {
        let logical_shape = interface
            .logical_slots()
            .into_iter()
            .map(|slot| {
                usize::try_from(slot.dim()).map_err(|_| {
                    exceptions::PyValueError::new_err(
                        "tensor data requires concrete representation dimensions",
                    )
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        let logical_axes = (0..logical_shape.len()).collect::<Vec<_>>();
        let representation_sorted = interface.rep_permutation.apply_slice(&logical_axes);
        let canonical_axes = interface
            .index_permutation
            .apply_slice(&representation_sorted);
        let canonical_shape = canonical_axes
            .iter()
            .map(|&logical_axis| logical_shape[logical_axis])
            .collect::<Vec<_>>();
        let logical_strides = Self::strides(&logical_shape)?;
        let canonical_strides = Self::strides(&canonical_shape)?;
        let size = logical_shape.iter().try_fold(1usize, |size, &dimension| {
            size.checked_mul(dimension).ok_or_else(|| {
                PyOverflowError::new_err("tensor dimensions overflow addressable storage")
            })
        })?;

        Ok(Self {
            logical_shape,
            logical_strides,
            canonical_axes,
            canonical_strides,
            size,
        })
    }

    fn strides(shape: &[usize]) -> PyResult<Vec<usize>> {
        let mut strides: Vec<usize> = vec![1; shape.len()];
        for axis in (0..shape.len().saturating_sub(1)).rev() {
            strides[axis] = strides[axis + 1]
                .checked_mul(shape[axis + 1])
                .ok_or_else(|| {
                    PyOverflowError::new_err("tensor dimensions overflow addressable storage")
                })?;
        }
        Ok(strides)
    }

    pub(crate) fn size(&self) -> usize {
        self.size
    }

    fn logical_expanded(&self, flat: usize) -> PyResult<Vec<usize>> {
        if flat >= self.size {
            return Err(PyIndexError::new_err("flat index out of bounds"));
        }
        let mut remainder = flat;
        Ok(self
            .logical_strides
            .iter()
            .map(|&stride| {
                let index = remainder / stride;
                remainder %= stride;
                index
            })
            .collect())
    }

    pub(crate) fn canonical_expanded(&self, logical: &[usize]) -> PyResult<Vec<usize>> {
        if logical.len() != self.logical_shape.len() {
            return Err(PyIndexError::new_err(format!(
                "expected {} tensor indices, got {}",
                self.logical_shape.len(),
                logical.len()
            )));
        }
        for (axis, (&index, &dimension)) in logical.iter().zip(&self.logical_shape).enumerate() {
            if index >= dimension {
                return Err(PyIndexError::new_err(format!(
                    "index {index} out of bounds for axis {axis} of size {dimension}"
                )));
            }
        }
        Ok(self
            .canonical_axes
            .iter()
            .map(|&logical_axis| logical[logical_axis])
            .collect())
    }

    pub(crate) fn canonical_flat(&self, logical_flat: usize) -> PyResult<usize> {
        let logical = self.logical_expanded(logical_flat)?;
        self.canonical_flat_from_expanded(&logical)
    }

    pub(crate) fn canonical_flat_from_expanded(&self, logical: &[usize]) -> PyResult<usize> {
        let canonical = self.canonical_expanded(logical)?;
        Ok(canonical
            .iter()
            .zip(&self.canonical_strides)
            .map(|(&index, &stride)| index * stride)
            .sum())
    }

    pub(crate) fn reorder_to_canonical<T>(&self, logical: Vec<T>) -> PyResult<Vec<T>> {
        if logical.len() != self.size {
            return Err(exceptions::PyValueError::new_err(format!(
                "tensor data has {} elements, but the logical interface requires {}",
                logical.len(),
                self.size
            )));
        }
        let mut canonical = std::iter::repeat_with(|| None)
            .take(self.size)
            .collect::<Vec<_>>();
        for (logical_flat, value) in logical.into_iter().enumerate() {
            canonical[self.canonical_flat(logical_flat)?] = Some(value);
        }
        Ok(canonical
            .into_iter()
            .map(|value| value.expect("axis permutations are bijective"))
            .collect())
    }
}

#[cfg(test)]
mod logical_data_layout_tests {
    use super::LogicalDataLayout;
    use spenso::structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{ExtendibleReps, RepName},
    };

    #[test]
    fn mixed_representation_data_stays_in_logical_row_major_order() {
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let interface = PartialStructure::from_logical_slots([
            mink.slot(PartialIndex::Explicit(AbstractIndex::Normal(0))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(1))),
        ]);
        let layout = LogicalDataLayout::from_interface(&interface).unwrap();

        assert_eq!(layout.logical_shape, vec![2, 3]);
        assert_eq!(layout.canonical_axes, vec![1, 0]);
        assert_eq!(
            layout.reorder_to_canonical(vec![0, 1, 2, 3, 4, 5]).unwrap(),
            vec![0, 3, 1, 4, 2, 5]
        );
        assert_eq!(
            (0..6)
                .map(|index| layout.canonical_flat(index).unwrap())
                .collect::<Vec<_>>(),
            vec![0, 2, 4, 1, 3, 5]
        );
    }

    #[test]
    fn explicit_index_sorting_is_part_of_the_layout_mapping() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let interface = PartialStructure::from_logical_slots([
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(8))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(3))),
        ]);
        let layout = LogicalDataLayout::from_interface(&interface).unwrap();

        assert_eq!(layout.canonical_axes, vec![1, 0]);
        assert_eq!(
            layout.reorder_to_canonical(vec![0, 1, 2, 3]).unwrap(),
            vec![0, 2, 1, 3]
        );
    }

    #[test]
    fn mixed_open_and_explicit_ports_keep_logical_data_order() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let interface = PartialStructure::from_logical_slots([
            euc.slot(PartialIndex::open(0)),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(3))),
        ]);
        let layout = LogicalDataLayout::from_interface(&interface).unwrap();

        assert_eq!(layout.canonical_axes, vec![1, 0]);
        assert_eq!(
            layout.reorder_to_canonical(vec![0, 1, 2, 3]).unwrap(),
            vec![0, 2, 1, 3]
        );
    }

    #[test]
    fn rank_zero_layout_contains_one_scalar_element() {
        let interface = PartialStructure::from_logical_slots([]);
        let layout = LogicalDataLayout::from_interface(&interface).unwrap();

        assert!(layout.logical_shape.is_empty());
        assert!(layout.canonical_axes.is_empty());
        assert_eq!(layout.size(), 1);
        assert_eq!(layout.canonical_flat(0).unwrap(), 0);
        assert_eq!(layout.reorder_to_canonical(vec![7]).unwrap(), vec![7]);
    }
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
        let layout = LogicalDataLayout::from_interface(&descriptor.interface)?;
        let mut storage_indices = vec![0; layout.canonical_axes.len()];
        for (storage_index, &logical_axis) in layout.canonical_axes.iter().enumerate() {
            storage_indices[logical_axis] = storage_index;
        }
        let structure = OrderedStructure::new(
            descriptor
                .interface
                .logical_slots()
                .into_iter()
                .zip(storage_indices)
                .map(|(slot, index)| slot.rep().slot(AbstractIndex::Normal(index)))
                .collect(),
        )
        .map_structure(|structure| ShadowedStructure {
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
        tensor: PermutedStructure<MixedTensor<f64, ShadowedStructure<AbstractIndex>>>,
        descriptor: StructuredAtom,
        descriptor_name: Option<Symbol>,
        descriptor_args: Vec<Atom>,
    ) -> Self {
        let tensor = PermutedStructure {
            structure: tensor.structure,
            rep_permutation: descriptor.interface.rep_permutation.clone(),
            index_permutation: descriptor.interface.index_permutation.clone(),
        };
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
            PermutedStructure::identity(ParamOrConcrete::new_scalar(ConcreteOrParam::Concrete(
                RealOrComplex::Real(value),
            ))),
            descriptor,
            descriptor_name,
            descriptor_args,
        )
    }
}

impl Deref for Spensor {
    type Target = MixedTensor<f64, ShadowedStructure<AbstractIndex>>;

    fn deref(&self) -> &Self::Target {
        &self.tensor.structure
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
            .structure
            .clone()
            .map_structure(|structure| ShadowedStructure {
                structure: structure.structure,
                global_name: Some(name),
                additional_args: (!args.is_empty()).then_some(args.clone()),
            });
        Self::from_storage_with_descriptor(
            PermutedStructure {
                structure: storage,
                rep_permutation: self.tensor.rep_permutation.clone(),
                index_permutation: self.tensor.index_permutation.clone(),
            },
            self.descriptor.clone(),
            Some(name),
            args,
        )
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
                structure.map_structure(|s| SparseTensor::<f64, _>::empty(s, 0.0).into()),
                descriptor,
                Some(name),
                args,
            ))
        } else if type_info.is_subclass_of::<PythonExpression>()? {
            Ok(Spensor::from_storage_with_descriptor(
                structure.map_structure(|s| {
                    ParamOrConcrete::Param(ParamTensor::from(SparseTensor::<Atom, _>::empty(
                        s,
                        Atom::Zero,
                    )))
                }),
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
    ///     The tensor data in row-major order
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
        let layout = LogicalDataLayout::from_interface(&descriptor.interface)?;
        let dense = match data {
            AtomsOrFloats::Floats(f) => DenseTensor::<f64, _>::from_data(
                layout.reorder_to_canonical(f)?,
                structure.structure,
            )
            .map_err(|e| PyOverflowError::new_err(e.to_string()))?
            .into(),
            AtomsOrFloats::Atoms(a) => ParamOrConcrete::Param(ParamTensor::from(
                DenseTensor::<Atom, _>::from_data(
                    layout.reorder_to_canonical(a)?,
                    structure.structure,
                )
                .map_err(|e| PyOverflowError::new_err(e.to_string()))?,
            )),
            AtomsOrFloats::Complex(c) => {
                MixedTensor::Concrete(RealOrComplexTensor::Complex(DataTensor::Dense(
                    DenseTensor::<Complex<f64>, _>::from_data(
                        layout.reorder_to_canonical(c)?,
                        structure.structure,
                    )
                    .map_err(|e| PyOverflowError::new_err(e.to_string()))?,
                )))
            }
        };

        let dense = PermutedStructure {
            structure: dense,
            rep_permutation: structure.rep_permutation,
            index_permutation: structure.index_permutation,
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
        self.tensor.structure = self.tensor.structure.clone().to_dense();
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
        self.tensor.structure = self.tensor.structure.clone().to_sparse();
    }

    fn __repr__(&self) -> String {
        format!("Tensor({})", display::format_concrete_tensor(self, false))
    }

    fn __str__(&self) -> String {
        display::format_concrete_tensor(self, false)
    }

    /// Format this concrete tensor using compact Spenso notation.
    #[pyo3(signature = (show_dimensions = false))]
    fn format_tensor(&self, show_dimensions: bool) -> String {
        display::format_concrete_tensor(self, show_dimensions)
    }

    /// Format this concrete tensor as Typst math source.
    #[pyo3(signature = (show_dimensions = false))]
    fn to_typst(&self, show_dimensions: bool) -> String {
        display::concrete_tensor_to_typst(self, show_dimensions)
    }

    /// Build Symbolica's rich display wrapper for this concrete tensor.
    #[pyo3(signature = (show_dimensions = false))]
    fn formatted(&self, show_dimensions: bool) -> PythonFormattedOutput {
        display::format_concrete_tensor_output(self, show_dimensions)
    }

    fn _repr_html_(&self) -> Option<String> {
        display::format_concrete_tensor_output(self, false).html
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

    fn __getitem__(&self, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        let layout = LogicalDataLayout::from_interface(&self.descriptor.interface)?;
        let size = layout.size();
        let get_owned_canonical = |index: usize| {
            self.get_owned_linear(index.into())
                .or_else(|| match &self.tensor.structure {
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
            let canonical_index = layout.canonical_flat(logical_index)?;
            get_owned_canonical(canonical_index)
                .ok_or_else(|| PyIndexError::new_err("flat index out of bounds"))
        };

        let out = match item {
            SliceOrIntOrExpanded::Int(i) => get_owned_linear(i)?,
            SliceOrIntOrExpanded::Expanded(idxs) => {
                let index = layout.canonical_flat_from_expanded(&idxs)?;
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
    ///     Index specification (int for flat index, list of int for coordinates)
    /// value : float, complex, or Expression
    ///     The value to set
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
        } else {
            return Err(eyre!("Value must be a PythonExpression or a float"));
        };

        let layout = LogicalDataLayout::from_interface(&self.descriptor.interface)?;
        if let Ok(flat_index) = item.extract::<usize>() {
            let canonical = layout.canonical_flat(flat_index)?;
            self.tensor.structure.set_flat(canonical.into(), value)
        } else if let Ok(expanded_idxs) = item.extract::<Vec<usize>>() {
            let canonical = layout.canonical_expanded(&expanded_idxs)?;
            self.tensor.structure.set(&canonical, value)
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

        let mut evaltensor = match &self.tensor.structure {
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
            .structure
            .clone()
            .scalar()
            .map(|r| PythonExpression { expr: r.into() })
            .ok_or_else(|| PyRuntimeError::new_err("No scalar found"))
    }

    /// Reference this tensor with new abstract indices and return a lazy network.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn index(
        &self,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor_reference(self.clone())?.index_network(py, indices, cook_indices)
    }

    /// Reference this tensor with new abstract indices and return a lazy network.
    #[pyo3(signature = (*indices, cook_indices = false))]
    fn __call__(
        &self,
        py: Python<'_>,
        indices: &Bound<'_, PyTuple>,
        cook_indices: bool,
    ) -> PyResult<SpensoNet> {
        self.index(py, indices, cook_indices)
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

    fn outer(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.outer(rhs)
    }

    #[pyo3(signature = (rhs, *, left, right))]
    fn contract(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: usize,
        right: usize,
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.contract(rhs, left, right)
    }

    #[pyo3(signature = (rhs, *, left, right))]
    fn compose(
        &self,
        rhs: ConvertibleToSpensoNet,
        left: (usize, usize),
        right: (usize, usize),
    ) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.compose(rhs, left, right)
    }

    fn dot(&self, rhs: ConvertibleToSpensoNet) -> PyResult<SpensoNet> {
        SpensoNet::from_tensor(self.clone())?.dot(rhs)
    }

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
                    PermutedStructure::identity(MixedTensor::Concrete(RealOrComplexTensor::Real(
                        eval.evaluate(s),
                    ))),
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
                    PermutedStructure::identity(MixedTensor::Concrete(
                        RealOrComplexTensor::Complex(eval.evaluate(s).map_data(|value| value)),
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
                    PermutedStructure::identity(MixedTensor::Concrete(
                        RealOrComplexTensor::Complex(self.eval.evaluate(s).map_data(|value| value)),
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

Parameters
----------
item : int or list of int
    Index specification (int for flat index, list of int for coordinates)
value : float, complex, or Expression
    The value to set

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
