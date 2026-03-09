use std::ops::Deref;

use anyhow::anyhow;

use pyo3::{
    exceptions::{PyIndexError, PyOverflowError, PyRuntimeError, PyTypeError},
    prelude::*,
    types::{PyFloat, PyType},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    generate::MethodType,
    impl_stub_type,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};
use spenso::{
    algebra::complex::{Complex, RealOrComplex},
    tensors::{
        data::{DenseTensor, GetTensorData, SetTensorData, SparseOrDense, SparseTensor},
        parametric::{ConcreteOrParam, ParamOrConcrete, ParamTensor},
    },
};

use crate::SliceOrIntOrExpanded;
use spenso::{
    network::library::symbolic::ExplicitKey,
    structure::{
        HasStructure, PermutedStructure, ScalarTensor, TensorStructure,
        abstract_index::AbstractIndex, permuted::Perm,
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, StorageTensor},
        parametric::MixedTensor,
    },
};
use symbolica::atom::Atom;

use symbolica::api::python::PythonExpression;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{define_stub_info_gatherer, derive::*};

use super::{
    ModuleInit, TensorElements,
    structure::{ConvertibleToIndexLess, SpensoStructure},
};

/// A library tensor class optimized for use in tensor libraries and networks.
///
/// Library tensors are similar to regular tensors but use explicit keys for efficient
/// lookup and storage in tensor libraries. They can be either dense or sparse and
/// store data as floats, complex numbers, or symbolic expressions.
///
/// LibraryTensors are designed for:
/// - Registration in TensorLibrary instances
/// - Use in tensor networks where structure reuse is important
/// - Efficient symbolic manipulation and pattern matching
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import LibraryTensor, TensorStructure, Representation
/// >>> rep = Representation.euc(3)
/// >>> structure = TensorStructure(rep, rep, name="T")
/// >>> data = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
/// >>> tensor = LibraryTensor.dense(structure, data)
/// >>> sparse_tensor = LibraryTensor.sparse(structure, float)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "LibraryTensor", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct LibrarySpensor {
    pub tensor: PermutedStructure<MixedTensor<f64, ExplicitKey<AbstractIndex>>>,
}

impl Deref for LibrarySpensor {
    type Target = MixedTensor<f64, ExplicitKey<AbstractIndex>>;

    fn deref(&self) -> &Self::Target {
        &self.tensor.structure
    }
}

impl ModuleInit for LibrarySpensor {}

pub enum AtomsOrFloats {
    Atoms(Vec<Atom>),
    Floats(Vec<f64>),
    Complex(Vec<Complex<f64>>),
}

impl<'a, 'py> FromPyObject<'a, 'py> for AtomsOrFloats {
    type Error = PyErr;

    fn extract(aind: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let aind = if let Ok(i) = aind.extract::<Vec<f64>>() {
            AtomsOrFloats::Floats(i)
        } else if let Ok(i) = aind.extract::<Vec<Complex<f64>>>() {
            AtomsOrFloats::Complex(i)
        } else if let Ok(i) = aind.extract::<Vec<PythonExpression>>() {
            AtomsOrFloats::Atoms(i.into_iter().map(|e| e.expr).collect())
        } else {
            return Err(PyTypeError::new_err(
                "Argument must be a list of floats, complex numbers, or Atoms",
            ));
        };
        Ok(aind)
    }
}

#[cfg(feature = "python_stubgen")]
impl_stub_type!(AtomsOrFloats = Vec<PythonExpression> | Vec<f64> | Vec<Complex<f64>>);

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl LibrarySpensor {
    pub fn structure(&self) -> SpensoStructure {
        SpensoStructure {
            structure: PermutedStructure {
                structure: self.tensor.structure.structure().clone(),
                rep_permutation: self.tensor.rep_permutation.clone(),
                index_permutation: self.tensor.index_permutation.clone(),
            },
        }
    }

    #[staticmethod]
    /// Create a new sparse empty library tensor with the given structure and data type.
    ///
    /// Creates a sparse tensor that initially contains no non-zero elements.
    /// Elements can be set individually using indexing operations.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorStructure, list of Representations, or list of int
    ///     The tensor structure defining shape and index properties
    /// type_info : type
    ///     The data type - either `float` or `Expression` class
    ///
    /// Returns
    /// -------
    /// LibraryTensor
    ///     A new sparse library tensor with all elements initially zero
    ///
    /// Examples
    /// --------
    /// >>> import symbolica as sp
    /// >>> from symbolica.community.spenso import LibraryTensor, TensorStructure, Representation
    /// >>> rep = Representation.euc(3)
    /// >>> structure = TensorStructure(rep, rep)
    /// >>> sparse_float = LibraryTensor.sparse(structure, float)
    /// >>> sparse_sym = LibraryTensor.sparse(structure, sp.Expression)
    /// >>> sparse_float[0, 0] = 1.0
    /// >>> sparse_float[1, 1] = 2.0
    pub fn sparse(
        structure: ConvertibleToIndexLess,
        type_info: Bound<'_, PyType>,
    ) -> PyResult<Self> {
        if type_info.is_subclass_of::<PyFloat>()? {
            Ok(Self {
                tensor: structure
                    .0
                    .structure
                    .map_structure(|s| SparseTensor::<f64, _>::empty(s, 0.0).into()),
            })
        } else if type_info.is_subclass_of::<PythonExpression>()? {
            Ok(Self {
                tensor: structure.0.structure.map_structure(|s| {
                    ParamOrConcrete::Param(ParamTensor::from(SparseTensor::<Atom, _>::empty(
                        s,
                        Atom::Zero,
                    )))
                }),
            })
        } else {
            Err(PyTypeError::new_err("Only float type supported"))
        }
    }

    #[staticmethod]
    /// Create a new dense library tensor with the given structure and data.
    ///
    /// Dense tensors store all elements explicitly in row-major order. The structure
    /// defines the tensor's shape and indexing properties.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorStructure, list of Representations, or list of int
    ///     The tensor structure defining shape and index properties
    /// data : list of float, complex, or Expression
    ///     The tensor data in row-major order
    ///
    /// Returns
    /// -------
    /// LibraryTensor
    ///     A new dense library tensor with the specified data
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.spenso import LibraryTensor, TensorStructure, Representation
    /// >>> rep = Representation.euc(2)
    /// >>> sigma = S("sigma")
    /// >>> structure = TensorStructure(rep, rep, name=sigma)
    /// >>> data = [0.0, 1.0, 1.0, 0.0]
    /// >>> tensor = LibraryTensor.dense(structure, data)
    /// >>> x, y = S("x", "y")
    /// >>> sym_data = [x, y, -y, x]
    /// >>> sym_tensor = LibraryTensor.dense(structure, sym_data)
    pub fn dense(structure: ConvertibleToIndexLess, data: AtomsOrFloats) -> PyResult<Self> {
        let dense = match data {
            AtomsOrFloats::Floats(f) => {
                DenseTensor::<f64, _>::from_data(f, structure.0.structure.structure)
                    .map_err(|e| PyOverflowError::new_err(e.to_string()))?
                    .into()
            }
            AtomsOrFloats::Atoms(a) => ParamOrConcrete::Param(ParamTensor::from(
                DenseTensor::<Atom, _>::from_data(a, structure.0.structure.structure)
                    .map_err(|e| PyOverflowError::new_err(e.to_string()))?,
            )),
            AtomsOrFloats::Complex(c) => {
                MixedTensor::Concrete(RealOrComplexTensor::Complex(DataTensor::Dense(
                    DenseTensor::<Complex<f64>, _>::from_data(c, structure.0.structure.structure)
                        .map_err(|e| PyOverflowError::new_err(e.to_string()))?,
                )))
            }
        };

        let dense = PermutedStructure {
            structure: dense,
            rep_permutation: structure.0.structure.rep_permutation,
            index_permutation: structure.0.structure.index_permutation,
        };

        Ok(Self {
            tensor: dense.permute_inds_wrapped(),
        })
    }
    #[staticmethod]
    /// Create a scalar library tensor with value 1.0.
    ///
    /// Returns
    /// -------
    /// LibraryTensor
    ///     A scalar library tensor containing the value 1.0
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import LibraryTensor
    /// >>> one = LibraryTensor.one()
    pub fn one() -> Self {
        Self {
            tensor: PermutedStructure::identity(ParamOrConcrete::new_scalar(
                ConcreteOrParam::Concrete(RealOrComplex::Real(1.)),
            )),
        }
    }

    #[staticmethod]
    /// Create a scalar library tensor with value 0.0.
    ///
    /// Returns
    /// -------
    /// LibraryTensor
    ///     A scalar library tensor containing the value 0.0
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import LibraryTensor
    /// >>> zero = LibraryTensor.zero()
    pub fn zero() -> Self {
        Self {
            tensor: PermutedStructure::identity(ParamOrConcrete::new_scalar(
                ConcreteOrParam::Concrete(RealOrComplex::Real(2.)),
            )),
        }
    }

    #[allow(clippy::wrong_self_convention)]
    /// Convert this library tensor to dense storage format.
    ///
    /// Converts sparse tensors to dense format in-place. Dense tensors are unchanged.
    /// This allocates memory for all tensor elements.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import LibraryTensor, TensorStructure, Representation
    ///
    /// rep = Representation.cof(2)
    /// structure = TensorStructure([rep, rep])
    /// tensor = LibraryTensor.sparse(structure, float)
    /// tensor[0, 0] = 1.0
    /// tensor.to_dense()  # Now stores all 4 elements explicitly
    /// ```
    fn to_dense(&mut self) {
        self.tensor.structure = self.tensor.structure.clone().to_dense();
    }

    #[allow(clippy::wrong_self_convention)]
    /// Convert this library tensor to sparse storage format.
    ///
    /// Converts dense tensors to sparse format in-place, only storing non-zero elements.
    /// This can save memory for tensors with many zero elements.
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import LibraryTensor, TensorStructure, Representation
    ///
    /// rep = Representation.euc(2)
    /// structure = TensorStructure(rep, rep)
    /// data = [1.0, 0.0, 0.0, 2.0]
    /// tensor = LibraryTensor.dense(structure, data)
    /// tensor.to_sparse()  # Now only stores 2 non-zero elements
    /// ```
    fn to_sparse(&mut self) {
        self.tensor.structure = self.tensor.structure.clone().to_sparse();
    }

    fn __repr__(&self) -> String {
        format!("Spensor(\n{})", self.tensor)
    }

    fn __str__(&self) -> String {
        format!("{}", self.tensor)
    }

    fn __len__(&self) -> usize {
        self.size().unwrap()
    }

    fn __getitem__(&self, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        let out = match item {
            SliceOrIntOrExpanded::Int(i) => self
                .get_owned_linear(i.into())
                .ok_or(PyIndexError::new_err("flat index out of bounds"))?,
            SliceOrIntOrExpanded::Expanded(idxs) => self
                .get_owned(&idxs)
                .map_err(|s| PyIndexError::new_err(s.to_string()))?,
            SliceOrIntOrExpanded::Slice(s) => {
                let r = s.indices(self.size().unwrap() as isize)?;

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

                let slice: Option<Vec<TensorElements>> = range
                    .step_by(step)
                    .map(|i| self.get_owned_linear(i.into()).map(TensorElements::from))
                    .collect();

                if let Some(slice) = slice {
                    return Ok(
                        Python::attach(|py| slice.into_pyobject(py).map(|a| a.unbind()))?
                            .into_any(),
                    );
                } else {
                    return Err(PyIndexError::new_err("slice out of bounds"));
                }
            }
        };

        Python::attach(|py| {
            TensorElements::from(out)
                .into_pyobject(py)
                .map(|a| a.unbind())
        })
    }

    /// Set library tensor element(s) at the specified index or indices.
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
    /// >>> from symbolica.community.spenso import LibraryTensor, TensorStructure, Representation
    /// >>> rep = Representation.euc(2)
    /// >>> structure = TensorStructure(rep, rep)
    /// >>> tensor = LibraryTensor.sparse(structure, float)
    /// >>> tensor[0] = 1.0
    /// >>> tensor[1, 1] = 2.0
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

        if let Ok(flat_index) = item.extract::<usize>() {
            self.tensor.structure.set_flat(flat_index.into(), value)
        } else if let Ok(expanded_idxs) = item.extract::<Vec<usize>>() {
            self.tensor.structure.set(&expanded_idxs, value)
        } else {
            Err(eyre!("Index must be an integer"))
        }
    }

    /// Extract the scalar value from a rank-0 (scalar) library tensor.
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
    /// >>> from symbolica.community.spenso import LibraryTensor
    /// >>> scalar_tensor = LibraryTensor.one()
    /// >>> value = scalar_tensor.scalar()
    fn scalar(&self) -> PyResult<PythonExpression> {
        self.tensor
            .structure
            .clone()
            .scalar()
            .map(|r| PythonExpression { expr: r.into() })
            .ok_or_else(|| PyRuntimeError::new_err("No scalar found"))
    }
}

impl From<DataTensor<f64, ExplicitKey<AbstractIndex>>> for LibrarySpensor {
    fn from(value: DataTensor<f64, ExplicitKey<AbstractIndex>>) -> Self {
        LibrarySpensor {
            tensor: PermutedStructure::identity(MixedTensor::Concrete(RealOrComplexTensor::Real(
                value,
            ))),
        }
    }
}

impl From<DataTensor<Complex<f64>, ExplicitKey<AbstractIndex>>> for LibrarySpensor {
    fn from(value: DataTensor<Complex<f64>, ExplicitKey<AbstractIndex>>) -> Self {
        LibrarySpensor {
            tensor: PermutedStructure::identity(MixedTensor::Concrete(
                RealOrComplexTensor::Complex(value.map_data(|c| c)),
            )),
        }
    }
}

#[cfg(feature = "python_stubgen")]
submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<LibrarySpensor>,
        attrs: &[],
        getters: &[],
        setters: &[],
        file: file!(),
        line: line!(),
        column: column!(),
        methods: &[
            MethodInfo {
                name: "__getitem__",
                parameters: &[
                    ParameterInfo {
                        name: "item",
                        kind: ParameterKind::PositionalOrKeyword,
                        default:ParameterDefault::None,
                        type_info: || TypeInfo::builtin("slice"),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: Vec::<TensorElements>::type_output,
                doc:r##"Get library tensor elements at the specified range of indices.

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
                        default:ParameterDefault::None,
                        type_info: ||  Vec::<usize>::type_input()|usize::type_input()
                    },
                ],
                r#type: MethodType::Instance,
                r#return: TensorElements::type_output,
                doc:r##"Get library tensor element at the specified index or indices.

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
                        default:ParameterDefault::None,
                        type_info: ||Vec::<usize>::type_input()|usize::type_input()
                    },
                    ParameterInfo {
                        name: "value",
                        kind: ParameterKind::PositionalOrKeyword,
                        default:ParameterDefault::None,
                        type_info: ||TensorElements::type_input()
                    },
                ],
                r#type: MethodType::Instance,
                r#return: TypeInfo::none,
                doc:r##"Set library tensor element(s) at the specified index or indices.

Parameters
----------
item : int or list of int
    Index specification (int for flat index, list of int for coordinates)
value : float, complex, or Expression
    The value to set

Examples
--------
>>> from symbolica.community.spenso import LibraryTensor, TensorStructure, Representation
>>> rep = Representation.euc(2)
>>> structure = TensorStructure(rep, rep)
>>> tensor = LibraryTensor.sparse(structure, float)
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
