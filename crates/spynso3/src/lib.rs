use std::{collections::HashMap, ops::Deref};

use eyre::eyre;

use library::SpensorLibrary;
use library_tensor::AtomsOrFloats;
use network::SpensoNet;

use pyo3::{
    exceptions::{self, PyIndexError, PyOverflowError, PyRuntimeError, PyTypeError},
    prelude::*,
    types::{PyComplex, PyFloat, PySlice, PyType},
    PyClass,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};

use spenso::{
    algebra::complex::{symbolica_traits::CompiledComplexEvaluatorSpenso, Complex, RealOrComplex},
    tensors::{
        data::{DenseTensor, GetTensorData, SetTensorData, SparseOrDense, SparseTensor},
        parametric::{
            atomcore::TensorAtomOps, ConcreteOrParam, EvalTensor, ParamOrConcrete, ParamTensor,
        },
    },
};

use spenso::{
    network::parsing::ShadowedStructure,
    structure::{
        abstract_index::AbstractIndex, permuted::Perm, HasStructure, PermutedStructure,
        ScalarTensor, TensorStructure,
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, StorageTensor},
        parametric::{LinearizedEvalTensor, MixedTensor},
    },
};
use structure::{ConvertibleToStructure, SpensoIndices};
use symbolica::{
    api::python::SymbolicaCommunityModule,
    atom::Atom,
    domains::{float::Complex as SymComplex, rational::Rational},
    evaluate::{CompileOptions, ExportSettings, FunctionMap, InlineASM, OptimizationSettings},
    poly::PolyVariable,
};

use symbolica::api::python::PythonExpression;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{define_stub_info_gatherer, derive::*, PyStubType, TypeInfo};

pub mod library;
pub mod library_tensor;
pub mod network;
pub mod structure;

trait ModuleInit: PyClass {
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<Self>()
    }
}

pub struct SpensoModule;

impl SymbolicaCommunityModule for SpensoModule {
    fn get_name() -> String {
        "spenso".to_string()
    }

    fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
        initialize_spenso(m)
    }

    fn initialize(_py: Python) -> PyResult<()> {
        idenso::representations::initialize();
        Ok(())
    }
}

pub(crate) fn initialize_spenso(m: &Bound<'_, PyModule>) -> PyResult<()> {
    use library_tensor::LibrarySpensor;
    use network::ExecutionMode;

    // m.add_function(?)?;
    SpensoNet::init(m)?;
    ExecutionMode::init(m)?;
    Spensor::init(m)?;
    LibrarySpensor::init(m)?;
    SpensoIndices::init(m)?;
    SpensorLibrary::init(m)?;
    Ok(())
}

/// A tensor class that can be either dense or sparse with flexible data types.
///
/// The tensor can store data as floats, complex numbers, or symbolic expressions.
/// Tensors have an associated structure that defines their shape and index properties.
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import Tensor, TensorIndices, Representation
/// >>> structure = TensorIndices(Representation.euc(4)(1))
/// >>> data = [1.0, 2.0, 3.0, 4.0]
/// >>> tensor = Tensor.dense(structure, data)
/// >>> sparse_tensor = Tensor.sparse(structure, float)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Tensor", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct Spensor {
    tensor: PermutedStructure<MixedTensor<f64, ShadowedStructure<AbstractIndex>>>,
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
    pub fn structure(&self) -> SpensoIndices {
        SpensoIndices {
            structure: PermutedStructure {
                structure: self.tensor.structure.structure().clone(),
                rep_permutation: self.tensor.rep_permutation.clone(),
                index_permutation: self.tensor.index_permutation.clone(),
            },
        }
    }

    #[staticmethod]
    /// Create a new sparse empty tensor with the given structure and data type.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorIndices or list of Slots
    ///     The tensor structure defining shape and index properties
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
    /// >>> from symbolica.community.spenso import Tensor, TensorIndices, Representation as R
    /// >>> structure = TensorIndices(R.euc(3)(1), R.euc(3)(2))
    /// >>> sparse_float = Tensor.sparse(structure, float)
    /// >>> sparse_sym = Tensor.sparse(structure, symbolica.Expression)
    pub fn sparse(
        structure: ConvertibleToStructure,
        type_info: Bound<'_, PyType>,
    ) -> PyResult<Spensor> {
        if type_info.is_subclass_of::<PyFloat>()? {
            Ok(Spensor {
                tensor: structure
                    .0
                    .structure
                    .map_structure(|s| SparseTensor::<f64, _>::empty(s, 0.0).into()),
            })
        } else if type_info.is_subclass_of::<PythonExpression>()? {
            Ok(Spensor {
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
    /// Create a new dense tensor with the given structure and data.
    ///
    /// Parameters
    /// ----------
    /// structure : TensorIndices or list of Slots
    ///     The tensor structure defining shape and index properties
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
    /// >>> from symbolica.community.spenso import Tensor, TensorIndices, Representation as R
    /// >>> structure = TensorIndices(R.euc(2)(1), R.euc(2)(2))
    /// >>> data = [1.0, 2.0, 3.0, 4.0]
    /// >>> tensor = Tensor.dense(structure, data)
    /// >>> x, y = S("x", "y")
    /// >>> sym_data = [x, y, x * y, x + y]
    /// >>> sym_tensor = Tensor.dense(structure, sym_data)
    pub fn dense(structure: ConvertibleToStructure, data: AtomsOrFloats) -> PyResult<Spensor> {
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

        Ok(Spensor {
            tensor: dense.permute_inds_wrapped(),
        })
    }
    #[staticmethod]
    /// Create a scalar tensor with value 1.0.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A scalar tensor containing the value 1.0
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import Tensor
    /// >>> one = Tensor.one()
    pub fn one() -> Spensor {
        Spensor {
            tensor: PermutedStructure::identity(ParamOrConcrete::new_scalar(
                ConcreteOrParam::Concrete(RealOrComplex::Real(1.)),
            )),
        }
    }

    #[staticmethod]
    /// Create a scalar tensor with value 0.0.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     A scalar tensor containing the value 0.0
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import Tensor
    /// >>> zero = Tensor.zero()
    pub fn zero() -> Spensor {
        Spensor {
            tensor: PermutedStructure::identity(ParamOrConcrete::new_scalar(
                ConcreteOrParam::Concrete(RealOrComplex::Real(2.)),
            )),
        }
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
    /// >>> from symbolica.community.spenso import Tensor, TensorIndices, Representation as R
    /// >>> structure = TensorIndices(R.euc(4)(2))
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
    /// >>> from symbolica.community.spenso import Tensor, TensorIndices, Representation as R
    /// >>> structure = TensorIndices(R.euc(2)(2), R.euc(2)(1))
    /// >>> data = [1.0, 0.0, 0.0, 2.0]
    /// >>> tensor = Tensor.dense(structure, data)
    /// >>> tensor.to_sparse()
    fn to_sparse(&mut self) {
        self.tensor.structure = self.tensor.structure.clone().to_sparse();
    }

    fn __repr__(&self) -> String {
        format!("Spensor(\n{})", self.tensor)
    }

    fn __str__(&self) -> String {
        format!("{}", self.tensor.structure)
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
    /// >>> from symbolica.community.spenso import Tensor, TensorIndices, Representation as R
    /// >>> structure = TensorIndices(R.euc(2)(2), R.euc(2)(1))
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

        if let Ok(flat_index) = item.extract::<usize>() {
            self.tensor.structure.set_flat(flat_index.into(), value)
        } else if let Ok(expanded_idxs) = item.extract::<Vec<usize>>() {
            self.tensor.structure.set(&expanded_idxs, value)
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
    /// >>> from symbolica.community.spenso import Tensor, TensorIndices, Representation as R
    /// >>> x, y = S("x", "y")
    /// >>> structure = TensorIndices(R.euc(2)(1))
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
            if let Ok(r) = v.expr.clone().try_into() {
                fn_map.add_constant(k.expr.clone(), r);
            } else {
                Err(exceptions::PyValueError::new_err(
                    "Constants must be rationals. If this is not possible, pass the value as a parameter",
                ))?
            }
        }

        for ((symbol, rename, args), body) in &funs {
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
                .add_function(symbol, rename.clone(), args, body.expr.clone())
                .map_err(|e| {
                    exceptions::PyValueError::new_err(format!("Could not add function: {}", e))
                })?;
        }

        let settings = OptimizationSettings {
            horner_iterations: iterations,
            n_cores,
            verbose,
            ..OptimizationSettings::default()
        };

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
    /// >>> from symbolica.community.spenso import Tensor
    /// >>> scalar_tensor = Tensor.one()
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

impl From<DataTensor<f64, ShadowedStructure<AbstractIndex>>> for Spensor {
    fn from(value: DataTensor<f64, ShadowedStructure<AbstractIndex>>) -> Self {
        Spensor {
            tensor: PermutedStructure::identity(MixedTensor::Concrete(RealOrComplexTensor::Real(
                value,
            ))),
        }
    }
}

impl From<DataTensor<Complex<f64>, ShadowedStructure<AbstractIndex>>> for Spensor {
    fn from(value: DataTensor<Complex<f64>, ShadowedStructure<AbstractIndex>>) -> Self {
        Spensor {
            tensor: PermutedStructure::identity(MixedTensor::Concrete(
                RealOrComplexTensor::Complex(value.map_data(|c| c)),
            )),
        }
    }
}
impl From<MixedTensor<f64, ShadowedStructure<AbstractIndex>>> for Spensor {
    fn from(value: MixedTensor<f64, ShadowedStructure<AbstractIndex>>) -> Self {
        Spensor {
            tensor: PermutedStructure::identity(value),
        }
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
#[pyclass(name = "TensorEvaluator", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct SpensoExpressionEvaluator {
    pub eval_rat: LinearizedEvalTensor<SymComplex<Rational>, ShadowedStructure<AbstractIndex>>,
    pub eval: Option<LinearizedEvalTensor<f64, ShadowedStructure<AbstractIndex>>>,
    pub eval_complex: LinearizedEvalTensor<Complex<f64>, ShadowedStructure<AbstractIndex>>,
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

        Ok(inputs.iter().map(|s| eval.evaluate(s).into()).collect())
    }

    /// Evaluate the expression for multiple inputs and return the results.
    fn evaluate_complex(&mut self, inputs: Vec<Vec<Complex<f64>>>) -> Vec<Spensor> {
        let eval = &mut self.eval_complex;

        inputs.iter().map(|s| eval.evaluate(s).into()).collect()
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
        let mut options = CompileOptions {
            optimization_level: optimization_level as usize,
            ..Default::default()
        };

        if let Some(compiler_path) = compiler_path {
            options.compiler = compiler_path.to_string();
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
                    ExportSettings {
                        include_header: true,
                        inline_asm,
                        custom_header,
                        // ..Default::default()
                    },
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
    name = "CompiledTensorEvaluator",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub struct SpensoCompiledExpressionEvaluator {
    pub eval: EvalTensor<CompiledComplexEvaluatorSpenso, ShadowedStructure<AbstractIndex>>,
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
            .map(|s| self.eval.evaluate(s).into())
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
>>> from symbolica.community.spenso import Tensor, TensorIndices, Representation
>>> rep = Representation.euc(2)
>>> structure = TensorIndices(rep(1), rep(2))
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
