//! Python component access in logical axis order.
use crate::{
    AtomsOrFloats, SliceOrIntOrExpanded, Spensor, TensorDataDescriptor, TensorElements,
    tensor_data_layout,
};
use pyo3::{
    exceptions::{PyIndexError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PyComplex, PyFloat, PyList, PySlice, PyTuple, PyType},
};
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pymethods;
#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;
use spenso::{
    algebra::complex::Complex,
    structure::TensorStructure,
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, StorageTensor},
        parametric::{ParamOrConcrete, ParamTensor},
    },
};
use symbolica::api::python::{ConvertibleToExpression, PythonExpression};

enum AxisSelection {
    Index(usize),
    Slice(Vec<usize>),
}

// Sparse contractions only visit stored entries. A map that changes the implicit
// zero must therefore materialize those entries before any later contraction.
fn mapped_storage<T: Clone, S: TensorStructure + Clone>(
    data: DataTensor<T, S>,
    is_zero: impl Fn(&T) -> bool,
) -> DataTensor<T, S> {
    match data {
        DataTensor::Sparse(sparse) if !is_zero(&sparse.zero) => {
            DataTensor::Dense(sparse.to_dense())
        }
        data => data,
    }
}

impl Spensor {
    pub(crate) fn normalized_component_index(index: isize, size: usize) -> PyResult<usize> {
        let normalized = if index < 0 {
            size as isize + index
        } else {
            index
        };
        if normalized < 0 || normalized as usize >= size {
            Err(PyIndexError::new_err(format!(
                "index {index} is outside an axis of size {size}"
            )))
        } else {
            Ok(normalized as usize)
        }
    }

    pub(crate) fn select_components(&self, selection: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let py = selection.py();
        let layout = tensor_data_layout(self.descriptor.structure())?;
        if let Ok(index) = selection.extract::<isize>() {
            return self.__getitem__(SliceOrIntOrExpanded::Int(Self::normalized_component_index(
                index,
                layout.size(),
            )?));
        }
        let selectors = selection.extract::<Vec<Bound<'_, PyAny>>>().map_err(|_| {
            PyTypeError::new_err("index must be an integer, slice, or coordinate tuple")
        })?;
        let shape = layout.logical_shape();
        if selectors.len() != shape.len() {
            return Err(PyIndexError::new_err(format!(
                "expected {} coordinate selectors, got {}",
                shape.len(),
                selectors.len()
            )));
        }
        let axes = selectors
            .iter()
            .zip(shape)
            .map(|(selector, &dimension)| {
                if let Ok(index) = selector.extract::<isize>() {
                    Ok(AxisSelection::Index(Self::normalized_component_index(
                        index, dimension,
                    )?))
                } else if let Ok(slice) = selector.cast::<PySlice>() {
                    let slice = slice.indices(dimension as isize)?;
                    Ok(AxisSelection::Slice(
                        (0..slice.slicelength)
                            .map(|i| (slice.start + i as isize * slice.step) as usize)
                            .collect(),
                    ))
                } else {
                    Err(PyTypeError::new_err(
                        "coordinate selectors must be integers or slices",
                    ))
                }
            })
            .collect::<PyResult<Vec<_>>>()?;
        self.selected_axis(py, &axes, &mut Vec::new())
    }

    fn selected_axis(
        &self,
        py: Python<'_>,
        axes: &[AxisSelection],
        coordinates: &mut Vec<usize>,
    ) -> PyResult<Py<PyAny>> {
        let Some((first, rest)) = axes.split_first() else {
            return self.__getitem__(SliceOrIntOrExpanded::Expanded(coordinates.clone()));
        };
        match first {
            AxisSelection::Index(index) => {
                coordinates.push(*index);
                let value = self.selected_axis(py, rest, coordinates);
                coordinates.pop();
                value
            }
            AxisSelection::Slice(indices) => {
                let result = PyList::empty(py);
                for &index in indices {
                    coordinates.push(index);
                    result.append(self.selected_axis(py, rest, coordinates)?)?;
                    coordinates.pop();
                }
                Ok(result.unbind().into_any())
            }
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage, on_success)]
#[pymethods]
impl Spensor {
    /// Python type used for the stored components.
    ///
    /// Returns
    /// -------
    /// type
    ///     One of float, complex, or Expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.dtype is float
    /// True
    #[getter]
    #[gen_stub(override_return_type(type_repr = "type[float] | type[complex] | type[Expression]"))]
    fn dtype<'py>(&self, py: Python<'py>) -> Bound<'py, PyType> {
        match &self.tensor {
            ParamOrConcrete::Concrete(RealOrComplexTensor::Real(_)) => py.get_type::<PyFloat>(),
            ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(_)) => {
                py.get_type::<PyComplex>()
            }
            ParamOrConcrete::Param(_) => py.get_type::<PythonExpression>(),
        }
    }

    /// How components are stored.
    ///
    /// Returns
    /// -------
    /// str
    ///     "dense" for all entries, or "sparse" for populated entries.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.storage
    /// 'dense'
    #[getter]
    #[gen_stub(override_return_type(type_repr="typing.Literal['dense', 'sparse']", imports=("typing")))]
    fn storage(&self) -> &'static str {
        let sparse = match &self.tensor {
            ParamOrConcrete::Concrete(RealOrComplexTensor::Real(t)) => {
                matches!(t, DataTensor::Sparse(_))
            }
            ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(t)) => {
                matches!(t, DataTensor::Sparse(_))
            }
            ParamOrConcrete::Param(t) => matches!(t.tensor, DataTensor::Sparse(_)),
        };
        if sparse { "sparse" } else { "dense" }
    }

    /// Number of external tensor axes.
    ///
    /// Returns
    /// -------
    /// int
    ///     Zero for a scalar.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.rank
    /// 2
    #[getter]
    fn rank(&self) -> usize {
        self.descriptor.rank()
    }

    /// Whether the tensor has no external axes.
    ///
    /// Returns
    /// -------
    /// bool
    ///     True exactly when rank is zero.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.is_scalar
    /// False
    #[getter]
    fn is_scalar(&self) -> bool {
        self.descriptor.is_scalar()
    }

    /// Dimensions in logical axis order.
    ///
    /// Returns
    /// -------
    /// tuple of int
    ///     Axis sizes in the order of axes, the current component view.
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
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int, ...]"))]
    fn shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            tensor_data_layout(self.descriptor.structure())?.logical_shape(),
        )
    }

    /// Copy this tensor's components and metadata.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     An independent copy. Changes to it do not alter the original.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> copied = tensor.copy()
    /// >>> copied.shape == tensor.shape
    /// True
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Apply a scalar function to component values.
    ///
    /// Parameters
    /// ----------
    /// callback : callable
    ///     Function receiving one component and returning its replacement.
    ///     The traversal order is unspecified.
    /// dtype : type, optional
    ///     Output component type: float, complex, or Expression. Defaults to
    ///     the current type; callback results must convert to the selected type.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     An independent tensor with the same axes and mapped components.
    ///
    /// Notes
    /// -----
    /// Sparse storage remains sparse when the implicit zero maps to zero.
    /// If it maps to a nonzero value, the result becomes dense. The implicit
    /// zero is mapped once; avoid relying on callback invocation counts.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.map_components(lambda value: value * 2)[1, 0]
    /// 6.0
    #[pyo3(signature = (callback, *, dtype=None))]
    fn map_components(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Callable[[Expression | float | complex], Expression | float | complex]", imports=("typing")))]
        callback: &Bound<'_, PyAny>,
        dtype: Option<&Bound<'_, PyType>>,
    ) -> PyResult<Self> {
        let original = self.dtype(py);
        let dtype = dtype.unwrap_or(&original);
        if !dtype.is(py.get_type::<PyFloat>())
            && !dtype.is(py.get_type::<PyComplex>())
            && !dtype.is(py.get_type::<PythonExpression>())
        {
            return Err(PyTypeError::new_err(
                "dtype must be float, complex, or Expression",
            ));
        }
        let values = match &self.tensor {
            ParamOrConcrete::Concrete(RealOrComplexTensor::Real(t)) => {
                t.map_data_ref_result(|value| callback.call1((*value,)).map(Bound::unbind))?
            }
            ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(t)) => {
                t.map_data_ref_result(|value| {
                    callback
                        .call1((PyComplex::from_doubles(py, value.re, value.im),))
                        .map(Bound::unbind)
                })?
            }
            ParamOrConcrete::Param(t) => t.tensor.map_data_ref_result(|value| {
                callback
                    .call1((PythonExpression::from(value.clone()),))
                    .map(Bound::unbind)
            })?,
        };
        let tensor = if dtype.is(py.get_type::<PyFloat>()) {
            ParamOrConcrete::Concrete(RealOrComplexTensor::Real(mapped_storage(
                values.map_data_ref_result(|v| {
                    TensorElements::validate_scalar(v.bind(py))?;
                    v.extract::<f64>(py)
                })?,
                |v| *v == 0.0,
            )))
        } else if dtype.is(py.get_type::<PyComplex>()) {
            ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(mapped_storage(
                values.map_data_ref_result(|v| {
                    TensorElements::validate_scalar(v.bind(py))?;
                    v.extract::<Complex<f64>>(py)
                })?,
                |v| v.re == 0.0 && v.im == 0.0,
            )))
        } else {
            ParamOrConcrete::Param(ParamTensor::from(mapped_storage(
                values.map_data_ref_result(|v| {
                    TensorElements::validate_scalar(v.bind(py))?;
                    v.extract::<ConvertibleToExpression>(py)
                        .map(|e| e.to_expression().expr)
                })?,
                |v| v.as_view().is_zero(),
            )))
        };
        Ok(Self::from_storage_with_descriptor(
            tensor,
            self.descriptor.clone(),
            self.descriptor_name,
            self.descriptor_args.clone(),
        ))
    }

    /// Copy numerical components into a NumPy array.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     Array with the tensor's logical shape and dtype float64 or complex128.
    ///
    /// Raises
    /// ------
    /// TypeError
    ///     The tensor stores symbolic Expressions. Evaluate them first.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> tensor.to_numpy().tolist()
    /// [[1.0, 2.0], [3.0, 4.0]]
    #[gen_stub(override_return_type(type_repr="numpy.typing.NDArray[numpy.float64] | numpy.typing.NDArray[numpy.complex128]", imports=("numpy", "numpy.typing")))]
    fn to_numpy(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let dtype = match &self.tensor {
            ParamOrConcrete::Concrete(RealOrComplexTensor::Real(_)) => "float64",
            ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(_)) => "complex128",
            ParamOrConcrete::Param(_) => {
                return Err(PyTypeError::new_err(
                    "to_numpy requires numeric components; evaluate the symbolic tensor first",
                ));
            }
        };
        let values = self.__getitem__(SliceOrIntOrExpanded::Slice(PySlice::full(py)))?;
        let array = PyModule::import(py, "numpy")?.call_method1("array", (values, dtype))?;
        Ok(array.call_method1("reshape", (self.shape(py)?,))?.unbind())
    }

    /// Copy a numerical array into a named tensor.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorExpression
    ///     Named descriptor with concrete dimensions matching the array shape.
    /// array : numpy.typing.ArrayLike
    ///     Real or complex array-like data. Axes follow the descriptor's current axes.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     Dense tensor with float or complex components.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.from_numpy(A, [[1.0, 2.0], [3.0, 4.0]])
    /// >>> tensor[1, 0]
    /// 3.0
    #[staticmethod]
    fn from_numpy(
        structure: TensorDataDescriptor,
        #[gen_stub(override_type(type_repr="numpy.typing.ArrayLike", imports=("numpy.typing")))]
        array: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let array = PyModule::import(array.py(), "numpy")?.call_method1("asarray", (array,))?;
        let shape = array.getattr("shape")?.extract::<Vec<usize>>()?;
        let expected = tensor_data_layout(structure.descriptor.structure())?;
        if shape != expected.logical_shape() {
            return Err(PyValueError::new_err(format!(
                "array shape {shape:?} does not match tensor shape {:?}",
                expected.logical_shape()
            )));
        }
        let kind = array
            .getattr("dtype")?
            .getattr("kind")?
            .extract::<String>()?;
        let values = array
            .call_method1("ravel", ("C",))?
            .call_method0("tolist")?;
        let data = match kind.as_str() {
            "b" | "i" | "u" | "f" => AtomsOrFloats::Floats(values.extract()?),
            "c" => AtomsOrFloats::Complex(values.extract()?),
            _ => {
                return Err(PyTypeError::new_err(
                    "from_numpy requires a real or complex numeric array",
                ));
            }
        };
        Self::dense(structure, data)
    }
}
