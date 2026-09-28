//! Python component access in logical axis order.
use crate::{
    AtomsOrFloats, SliceOrIntOrExpanded, Spensor, TensorDataDescriptor, tensor_data_layout,
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
#[pymethods]
impl Spensor {
    /// Python component type: float, complex, or Expression.
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

    /// Component storage format, independent of tensor rank or expression type.
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

    /// Number of external axes.
    #[getter]
    fn rank(&self) -> usize {
        self.descriptor.rank()
    }

    /// Whether the tensor has no external axes.
    #[getter]
    fn is_scalar(&self) -> bool {
        self.descriptor.is_scalar()
    }

    /// Dimensions in logical axis order.
    #[getter]
    #[gen_stub(override_return_type(type_repr = "tuple[int, ...]"))]
    fn shape<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            tensor_data_layout(self.descriptor.structure())?.logical_shape(),
        )
    }

    /// Independent copy of component storage and metadata.
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Apply a scalar callback to components, preserving sparse storage and logical axes.
    ///
    /// The default dtype is unchanged; select float, complex, or Expression explicitly
    /// to convert it. Sparse defaults are mapped once; if zero maps to nonzero,
    /// the result becomes dense so later contractions include every component.
    /// Callback traversal order is unspecified. The input tensor is never modified.
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
                values.map_data_ref_result(|v| v.extract::<f64>(py))?,
                |v| *v == 0.0,
            )))
        } else if dtype.is(py.get_type::<PyComplex>()) {
            ParamOrConcrete::Concrete(RealOrComplexTensor::Complex(mapped_storage(
                values.map_data_ref_result(|v| v.extract::<Complex<f64>>(py))?,
                |v| v.re == 0.0 && v.im == 0.0,
            )))
        } else {
            ParamOrConcrete::Param(ParamTensor::from(mapped_storage(
                values.map_data_ref_result(|v| {
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

    /// Copy numeric components into a NumPy array in logical axis order.
    /// Imports NumPy only when called. Symbolic tensors must first be evaluated.
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

    /// Copy a numeric NumPy array into a named tensor, checking the full logical shape.
    /// Real numeric arrays become float64 storage; complex arrays become complex128.
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
