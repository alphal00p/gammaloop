use std::{cell::Cell, collections::HashMap};

use eyre::eyre;
use idenso::IndexTooling;
use pyo3::{
    FromPyObject, PyErr, exceptions,
    prelude::*,
    types::{PyAny, PyComplex},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
};

use spenso::algebra::complex::{Complex, RealOrComplex};
use spenso::network::library::function_lib::{PanicMissingConcrete, SymbolLib};
use spenso::network::library::{FunctionLibraryError, function_lib::INBUILTS};
use spenso::network::parsing::ShadowedStructure;
use spenso::tensors::complex::RealOrComplexTensor;
use spenso::tensors::data::StorageTensor;
use spenso::{
    network::library::symbolic::{ExplicitKey, TensorLibrary},
    structure::{
        Canonicalized, HasStructure, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialStructure, PartialStructureExt},
        slot::IsAbstractSlot,
    },
    tensors::parametric::MixedTensor,
};
use symbolica::atom::{Atom, DefaultNamespace, SymbolBuilder};
use symbolica::{
    api::python::PythonExpression,
    atom::{AtomView, Symbol},
    function,
};

use crate::{
    Spensor,
    broadcast::SpensoBroadcastFunction,
    expression::{TensorExpression, value_to_structured_atom},
    structure::SpensoName,
};

use super::ModuleInit;

/// A library for registering and managing tensor templates and structures.
///
/// The TensorLibrary provides a centralized registry for tensor definitions that can be
/// reused across tensor networks and expressions. It manages tensor structures with their
/// associated names and can resolve symbolic references to registered tensors.
///
/// ```python
/// from symbolica.community.spenso import Tensor, TensorLibrary, TensorName, Representation
///
/// lib = TensorLibrary()
/// rep = Representation.euc(3)
/// name = TensorName("my_tensor")
/// structure = name(rep, rep)
/// tensor = Tensor.dense(structure, [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
/// lib.register(tensor)
/// tensor_ref = lib[name]
/// ```
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorLibrary", module = "symbolica.community.spenso")]
pub struct SpensorLibrary {
    pub(crate) library: TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>,
    /// Logical interfaces retained for references registered through Python.
    references: HashMap<ExplicitKey<AbstractIndex>, PartialStructure>,
}

/// A registry of elementwise tensor functions used during network evaluation.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorFunctionLibrary", module = "symbolica.community.spenso")]
pub struct SpensorFunctionLibrary {
    pub(crate) library:
        SymbolLib<RealOrComplexTensor<f64, ShadowedStructure<AbstractIndex>>, PanicMissingConcrete>,
}

impl ModuleInit for SpensorLibrary {}
impl ModuleInit for SpensorFunctionLibrary {}

pub enum ConvertibleToSymbol {
    Name(String),
    Symbol(PythonExpression),
}

impl ConvertibleToSymbol {
    fn symbol(&self) -> eyre::Result<Symbol> {
        match self {
            ConvertibleToSymbol::Name(name) => {
                let namespace = DefaultNamespace {
                    namespace: "spenso_python".into(),
                    data: "",
                    file: "".into(),
                    line: 0,
                };
                Ok(SymbolBuilder::new(namespace.attach_namespace(name))
                    .build()
                    .map_err(|error| eyre!(error))?)
            }
            ConvertibleToSymbol::Symbol(symbol) => {
                if let AtomView::Var(a) = symbol.as_view() {
                    Ok(a.get_symbol())
                } else {
                    Err(eyre::eyre!("Symbol is not a variable"))
                }
            }
        }
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToSymbol {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.extract::<String>() {
            Ok(ConvertibleToSymbol::Name(a))
        } else if let Ok(num) = ob.extract::<PythonExpression>() {
            Ok(ConvertibleToSymbol::Symbol(num))
        } else if let Ok(num) = ob.extract::<SpensoName>() {
            Ok(ConvertibleToSymbol::Symbol(PythonExpression {
                expr: Atom::var(num.name),
            }))
        } else {
            Err(exceptions::PyTypeError::new_err(
                "Cannot convert to expression",
            ))
        }
    }
}
#[allow(clippy::new_without_default)]
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl SpensorFunctionLibrary {
    #[new]
    /// Create a new empty tensor function library.
    ///
    /// Initializes an empty library ready for registering tensor functions.
    ///
    /// Returns
    /// -------
    /// TensorFunctionLibrary
    ///     A new empty function library
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorFunctionLibrary
    /// >>> lib = TensorFunctionLibrary()
    pub fn new() -> Self {
        let mut a = Self {
            library: PanicMissingConcrete::new_lib(),
        };

        a.library.insert(INBUILTS.conj, |a| match a {
            RealOrComplexTensor::Complex(c) => {
                RealOrComplexTensor::Complex(c.map_data(|x| x.conj()))
            }
            RealOrComplexTensor::Real(r) => RealOrComplexTensor::Real(r),
        });
        a.library
            .insert_scalar_fallible(INBUILTS.conj, |scalar| Ok(scalar.spenso_conj()));

        a
    }

    /// Register an elementwise callback for concrete tensor execution.
    ///
    /// The callback receives a Python `float` or `complex` matching each input value and
    /// must return either type. The result uses complex storage if any returned value is
    /// complex; otherwise it uses real storage.
    pub fn register(&mut self, function: &SpensoBroadcastFunction, callback: Py<PyAny>) {
        let name = function.name;
        let scalar_callback = Python::attach(|py| callback.clone_ref(py));
        self.library.insert_scalar_fallible(name, move |scalar| {
            let Ok(value) = symbolica::domains::float::Complex::<f64>::try_from(scalar.as_view())
            else {
                return Ok(function!(name, scalar));
            };
            Python::attach(|py| -> PyResult<Atom> {
                let result = if value.im == 0.0 {
                    scalar_callback.call1(py, (value.re,))?
                } else {
                    let value = PyComplex::from_doubles(py, value.re, value.im);
                    scalar_callback.call1(py, (value,))?
                };
                if let Ok(value) = result.extract::<f64>(py) {
                    Ok(Atom::num(value))
                } else {
                    let value = result.extract::<Complex<f64>>(py)?;
                    Ok(Atom::num(value.re) + Atom::num(value.im) * Atom::i())
                }
            })
            .map_err(|error| {
                FunctionLibraryError::Other(eyre!(
                    "broadcast function `{name}` failed during concrete scalar execution: {error}"
                ))
            })
        });
        self.library.insert_fallible(name, move |tensor| {
            Python::attach(|py| -> PyResult<_> {
                let saw_output = Cell::new(false);
                let output_is_complex = Cell::new(false);
                let extract = |result: Py<PyAny>| -> PyResult<RealOrComplex<f64>> {
                    saw_output.set(true);
                    if let Ok(value) = result.extract::<f64>(py) {
                        Ok(RealOrComplex::Real(value))
                    } else {
                        let value = result.extract::<Complex<f64>>(py)?;
                        output_is_complex.set(true);
                        Ok(RealOrComplex::Complex(value))
                    }
                };
                let (tensor, input_is_complex) = match tensor {
                    RealOrComplexTensor::Real(tensor) => (
                        tensor.map_data_ref_result(|value| extract(callback.call1(py, (*value,))?)),
                        false,
                    ),
                    RealOrComplexTensor::Complex(tensor) => (
                        tensor.map_data_ref_result(|value| {
                            let value = PyComplex::from_doubles(py, value.re, value.im);
                            extract(callback.call1(py, (value,))?)
                        }),
                        true,
                    ),
                };
                let tensor = tensor?;
                let output_is_complex =
                    output_is_complex.get() || (!saw_output.get() && input_is_complex);

                if output_is_complex {
                    Ok(RealOrComplexTensor::Complex(
                        tensor.map_data(RealOrComplex::to_complex),
                    ))
                } else {
                    Ok(RealOrComplexTensor::Real(tensor.map_data(|value| {
                        let RealOrComplex::Real(value) = value else {
                            unreachable!("complex callback output was detected before conversion")
                        };
                        value
                    })))
                }
            })
            .map_err(|error| {
                FunctionLibraryError::Other(eyre!(
                    "broadcast function `{name}` failed during concrete execution: {error}"
                ))
            })
        });
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ConvertibleToSymbol {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        SpensoName::type_input() | PythonExpression::type_input() | String::type_input()
    }
}

struct ExactLibraryReference {
    key: Canonicalized<ExplicitKey<AbstractIndex>>,
    interface: PartialStructure,
    name: Symbol,
    args: Vec<Atom>,
}

impl ExactLibraryReference {
    fn from_expression(expression: &PyRef<'_, TensorExpression>) -> PyResult<Self> {
        let descriptor = TensorExpression::structured(expression);
        let name = TensorExpression::descriptor_name(expression).ok_or_else(|| {
            exceptions::PyValueError::new_err("exact tensor library lookup requires a name")
        })?;
        let args = TensorExpression::descriptor_args(expression);
        Self::new(descriptor.interface, name, args)
    }

    fn new(interface: PartialStructure, name: Symbol, args: Vec<Atom>) -> PyResult<Self> {
        let rank = interface.canonical().order();
        let open = interface.open_positions();
        if rank != 0 && open.len() != rank {
            let explicit = (0..rank)
                .filter(|position| !open.contains(position))
                .collect::<Vec<_>>();
            return Err(exceptions::PyValueError::new_err(format!(
                "tensor library keys require a fully unresolved interface; explicit ports are {explicit:?}"
            )));
        }
        let key = ExplicitKey::from_iter(
            interface.logical_slots().into_iter().map(|slot| slot.rep()),
            name,
            (!args.is_empty()).then_some(args.clone()),
        );
        Ok(Self {
            key,
            interface,
            name,
            args,
        })
    }

    fn signature(&self) -> String {
        let arguments = self
            .args
            .iter()
            .map(ToString::to_string)
            .chain(
                self.interface
                    .logical_slots()
                    .into_iter()
                    .map(|slot| slot.rep().to_string()),
            )
            .collect::<Vec<_>>()
            .join(", ");
        format!("{}({arguments})", self.name)
    }
}

enum LibraryReference {
    Symbol(Symbol),
    Exact(Box<ExactLibraryReference>),
}

pub struct ConvertibleToLibraryReference(LibraryReference);

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToLibraryReference {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(expression) = value.extract::<PyRef<'_, TensorExpression>>() {
            return ExactLibraryReference::from_expression(&expression)
                .map(Box::new)
                .map(LibraryReference::Exact)
                .map(Self);
        }
        value
            .extract::<ConvertibleToSymbol>()?
            .symbol()
            .map(LibraryReference::Symbol)
            .map(Self)
            .map_err(|error| exceptions::PyTypeError::new_err(error.to_string()))
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ConvertibleToLibraryReference {
    fn type_input() -> pyo3_stub_gen::TypeInfo {
        TensorExpression::type_input() | ConvertibleToSymbol::type_input()
    }

    fn type_output() -> pyo3_stub_gen::TypeInfo {
        Self::type_input()
    }
}

fn tensor_reference(
    py: Python<'_>,
    name: Symbol,
    args: Vec<Atom>,
    interface: PartialStructure,
) -> PyResult<Py<TensorExpression>> {
    let structure = ExplicitKey::from_iter(
        interface.logical_slots().into_iter().map(|slot| slot.rep()),
        name,
        (!args.is_empty()).then_some(args.clone()),
    );
    let value = value_to_structured_atom(&structure)?;
    TensorExpression::from_atom_interface_descriptor(py, value.atom, interface, Some(name), args)
}

#[allow(clippy::new_without_default)]
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl SpensorLibrary {
    #[new]
    /// Create a new empty tensor library.
    ///
    /// Initializes an empty library ready for registering tensor structures.
    /// The library automatically manages internal tensor IDs.
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     A new empty tensor library
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorLibrary
    /// >>> lib = TensorLibrary()
    pub fn new() -> Self {
        let mut a = Self {
            library: TensorLibrary::new(),
            references: HashMap::new(),
        };
        a.library.update_ids();
        a
    }

    #[staticmethod]
    /// Create a new empty tensor library (static method).
    ///
    /// Initializes an empty library ready for registering tensor structures.
    /// The library automatically manages internal tensor IDs.
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     A new empty tensor library
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorLibrary
    /// >>> lib = TensorLibrary.construct()
    pub fn construct() -> Self {
        Self::new()
    }

    /// Register a tensor in the library.
    ///
    /// Adds a tensor template to the library that can be referenced by name
    /// in tensor networks and symbolic expressions. The tensor must have a name
    /// set in its structure.
    ///
    /// Parameters
    /// ----------
    /// tensor : Tensor
    ///     A named tensor with a fully unresolved interface, or a named scalar
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import Tensor, TensorLibrary, TensorName, Representation
    /// >>> lib = TensorLibrary()
    /// >>> rep = Representation.euc(3)
    /// >>> name = TensorName("my_tensor")
    /// >>> structure = name(rep, rep)
    /// >>> tensor = Tensor.dense(structure, [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    /// >>> lib.register(tensor)
    /// >>> tensor_ref = lib[name]
    pub fn register(&mut self, tensor: PyRef<'_, Spensor>) -> PyResult<()> {
        let name = tensor.descriptor_name.ok_or_else(|| {
            exceptions::PyValueError::new_err(
                "TensorLibrary.register() requires a named Tensor; use tensor.with_name(...) first",
            )
        })?;
        let reference = ExactLibraryReference::new(
            tensor.descriptor.interface.clone(),
            name,
            tensor.descriptor_args.clone(),
        )?;
        let storage = tensor
            .tensor
            .clone()
            .map_structure(|_| reference.key.canonical().clone());
        self.library
            .insert_explicit(reference.key.clone().map_canonical(|_| storage));
        self.references
            .insert(reference.key.into_canonical(), reference.interface);
        Ok(())
    }

    /// Retrieve a registered tensor structure by name.
    ///
    /// Looks up a previously registered tensor by its name and returns
    /// a reference structure that can be used to create new tensor instances.
    ///
    /// Parameters
    /// ----------
    /// key : TensorExpression, TensorName, Expression, or str
    ///     An exact unresolved tensor signature, or a symbol-only convenience key
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     An atomic reference with the requested exact interface, or the registered
    ///     logical interface for a symbol-only lookup
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     If the tensor name is not found in the library
    ///
    /// Examples
    /// --------
    /// >>> exact = lib[TensorName("T")(1, Representation.euc(3))]
    /// >>> unique_by_name = lib["T"]
    pub fn __getitem__(
        &self,
        py: Python<'_>,
        key: ConvertibleToLibraryReference,
    ) -> PyResult<Py<TensorExpression>> {
        match key.0 {
            LibraryReference::Exact(reference) => {
                self.library
                    .get(reference.key.canonical())
                    .map_err(|error| exceptions::PyKeyError::new_err(error.to_string()))?;
                tensor_reference(py, reference.name, reference.args, reference.interface)
            }
            LibraryReference::Symbol(symbol) => {
                let key = self.library.get_key_from_name(symbol).map_err(|error| {
                    let mut variants = self
                        .references
                        .iter()
                        .filter(|(key, _)| key.global_name == Some(symbol))
                        .map(|(key, interface)| ExactLibraryReference {
                            key: Canonicalized::identity(key.clone()),
                            interface: interface.clone(),
                            name: symbol,
                            args: key.additional_args.clone().unwrap_or_default(),
                        })
                        .map(|reference| reference.signature())
                        .collect::<Vec<_>>();
                    variants.sort();
                    let message = if variants.len() > 1 {
                        format!(
                            "tensor name `{symbol}` is ambiguous; registered signatures: {}",
                            variants.join(", ")
                        )
                    } else {
                        error.to_string()
                    };
                    exceptions::PyKeyError::new_err(message)
                })?;
                if let Some(interface) = self.references.get(&key) {
                    return tensor_reference(
                        py,
                        symbol,
                        key.additional_args.clone().unwrap_or_default(),
                        interface.clone(),
                    );
                }
                let structure = self
                    .library
                    .get(&key)
                    .map_err(|error| exceptions::PyKeyError::new_err(error.to_string()))?
                    .into_owned()
                    .map_canonical(|tensor| tensor.structure().clone());
                TensorExpression::from_structure(py, &structure)
            }
        }
    }

    #[staticmethod]
    /// Create a library pre-loaded with High Energy Physics tensor definitions.
    ///
    /// Returns a library containing standard HEP tensors such as gamma matrices,
    /// color generators, metric tensors, and other commonly used structures in
    /// particle physics calculations. They are floating point tensors with f64 precision.
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     A TensorLibrary pre-populated with HEP tensor definitions
    ///
    /// Examples
    /// --------
    /// >>> import symbolica
    /// >>> from symbolica.community.spenso import TensorLibrary, TensorName
    /// >>> hep_lib = TensorLibrary.hep_lib()
    /// >>> gamma_structure = hep_lib[symbolica.S("spenso::gamma")]
    pub fn hep_lib() -> Self {
        Self {
            library: spenso_hep_lib::hep_lib(1., 0.),
            references: HashMap::new(),
        }
    }

    #[staticmethod]
    /// Create a library pre-loaded with High Energy Physics tensor definitions.
    ///
    /// Returns a library containing standard HEP tensors such as gamma matrices,
    /// color generators, metric tensors, and other commonly used structures in
    /// particle physics calculations. They are tensors with atom numeric entries.
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     A TensorLibrary pre-populated with HEP tensor definitions
    ///
    /// Examples
    /// --------
    /// >>> import symbolica
    /// >>> from symbolica.community.spenso import TensorLibrary, TensorName
    /// >>> hep_lib = TensorLibrary.hep_lib_atom()
    /// >>> gamma_structure = hep_lib[symbolica.S("spenso::gamma")]
    pub fn hep_lib_atom() -> Self {
        Self {
            library: spenso_hep_lib::hep_lib_atom(),
            references: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use super::*;
    use idenso::{color::CS, representations::initialize};
    use pyo3::types::PyFloat;
    use spenso::network::{ExecutionResult, Sequential, SmallestDegree};
    use spenso::structure::{
        OrderedStructure,
        dimension::Dimension,
        partial::{OpenPortId, PartialIndex},
        representation::{ExtendibleReps, RepName},
    };
    use spenso::tensors::data::{DataTensor, DenseTensor, SparseTensor};
    use spenso_hep_lib::HEP_LIB;
    use symbolica::{parse, symbol};

    #[test]
    fn builtin_scalar_conjugation_execution_preserves_exact_coefficients() {
        initialize();
        let library = SpensorFunctionLibrary::new();
        let mut network =
            crate::network::ParsingNet::from_scalar(parse!("1/3 + 2i/7")).fun(INBUILTS.conj);
        network
            .execute::<Sequential, SmallestDegree, _, _, _>(&*HEP_LIB, &library.library)
            .unwrap();
        let ExecutionResult::Val(conjugated) = network.result_scalar().unwrap() else {
            panic!("conjugation should produce a scalar value")
        };

        assert_eq!(conjugated.into_owned(), parse!("1/3 - 2i/7"));
    }

    #[test]
    fn tensor_callbacks_choose_numeric_kind_from_their_outputs() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let structure = OrderedStructure::new(vec![
                ExtendibleReps::EUCLIDEAN
                    .new_rep(Dimension::Concrete(2))
                    .slot(AbstractIndex::Normal(0)),
            ])
            .map_canonical(|structure| ShadowedStructure {
                structure,
                global_name: None,
                additional_args: None,
            })
            .into_canonical();
            let mut library = SpensorFunctionLibrary::new();

            let promote = SpensoBroadcastFunction {
                name: symbol!("spynso_callback_promote"),
            };
            let callback = CString::new("lambda value: 1j * value").unwrap();
            library.register(&promote, py.eval(callback.as_c_str(), None, None)?.unbind());
            let input = RealOrComplexTensor::Real(DataTensor::Dense(DenseTensor {
                data: vec![1.0, 2.0],
                structure: structure.clone(),
            }));
            let output = library.library.functions[&promote.name](input).unwrap();
            let RealOrComplexTensor::Complex(DataTensor::Dense(output)) = output else {
                panic!("a complex callback result must promote real tensor storage")
            };
            assert_eq!(
                output.data,
                vec![Complex::new(0.0, 1.0), Complex::new(0.0, 2.0)]
            );

            let demote = SpensoBroadcastFunction {
                name: symbol!("spynso_callback_demote"),
            };
            let callback = CString::new("lambda value: abs(value)").unwrap();
            library.register(&demote, py.eval(callback.as_c_str(), None, None)?.unbind());
            let input = RealOrComplexTensor::Complex(DataTensor::Dense(DenseTensor {
                data: vec![Complex::new(3.0, 4.0), Complex::new(-5.0, 12.0)],
                structure: structure.clone(),
            }));
            let output = library.library.functions[&demote.name](input).unwrap();
            let RealOrComplexTensor::Real(DataTensor::Dense(output)) = output else {
                panic!("real callback results must demote complex tensor storage")
            };
            assert_eq!(output.data, vec![5.0, 13.0]);

            let map_default = SpensoBroadcastFunction {
                name: symbol!("spynso_callback_map_sparse_default"),
            };
            let callback = CString::new("lambda value: 1j * (value + 1)").unwrap();
            library.register(
                &map_default,
                py.eval(callback.as_c_str(), None, None)?.unbind(),
            );
            let input = RealOrComplexTensor::Real(DataTensor::Sparse(SparseTensor {
                elements: HashMap::new(),
                zero: 2.0,
                structure,
            }));
            let output = library.library.functions[&map_default.name](input).unwrap();
            let RealOrComplexTensor::Complex(DataTensor::Sparse(output)) = output else {
                panic!("a complex sparse default must promote tensor storage")
            };
            assert_eq!(output.zero, Complex::new(0.0, 3.0));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn exact_library_reference_preserves_full_signature_and_logical_reps() {
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let name = symbol!("spynso_exact_library_reference");
        let interface = PartialStructure::from_logical_slots([
            mink.slot(PartialIndex::Open(OpenPortId(0))),
            euc.slot(PartialIndex::Open(OpenPortId(1))),
        ]);
        let reference = ExactLibraryReference::new(interface, name, vec![Atom::num(7)]).unwrap();

        assert_eq!(reference.key.canonical().global_name, Some(name));
        assert_eq!(
            reference.key.canonical().additional_args,
            Some(vec![Atom::num(7)])
        );
        let canonical = reference.key.canonical().reps();
        assert_eq!(
            reference.key.layout().canonical_to_logical(&canonical),
            vec![mink, euc]
        );
        assert!(reference.signature().contains("7"));
    }

    #[test]
    fn typed_structure_constant_uses_the_existing_library_key() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let factory = py
                .get_type::<TensorExpression>()
                .call_method1("f", (8,))?
                .extract::<Py<TensorExpression>>()?;
            let reference = ExactLibraryReference::from_expression(&factory.bind(py).borrow())?;
            assert_eq!(reference.key, CS.f_strct::<AbstractIndex>(8));
            assert!(reference.args.is_empty());

            let tensor = Py::new(
                py,
                Spensor::sparse(
                    factory.bind(py).as_any().extract()?,
                    py.get_type::<PyFloat>(),
                )?,
            )?;
            let mut library = SpensorLibrary::new();
            library.register(tensor.bind(py).borrow())?;

            let stored = library.__getitem__(
                py,
                ConvertibleToLibraryReference(LibraryReference::Exact(Box::new(reference))),
            )?;
            assert_ne!(stored.bind(py).borrow().as_super().expr, Atom::Zero);
            let indexed = stored.bind(py).call1(("a", "b", "c"))?;
            assert_eq!(indexed.getattr("rank")?.extract::<usize>()?, 3);
            assert_ne!(
                indexed
                    .call_method0("to_expression")?
                    .extract::<PythonExpression>()?
                    .expr,
                Atom::Zero
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn exact_lookup_materializes_data_in_the_requested_logical_order() {
        initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
            let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
            let name = spenso::network::tags::SPENSO_TAG
                .tensor_symbol("spynso_exact_lookup_logical_order");
            let registered = tensor_reference(
                py,
                name,
                Vec::new(),
                PartialStructure::from_logical_slots([
                    mink.slot(PartialIndex::Open(OpenPortId(0))),
                    euc.slot(PartialIndex::Open(OpenPortId(1))),
                ]),
            )?;
            let tensor = Py::new(
                py,
                Spensor::dense(
                    registered.bind(py).as_any().extract()?,
                    crate::AtomsOrFloats::Floats(vec![0., 1., 2., 3., 4., 5.]),
                )?,
            )?;
            let mut library = SpensorLibrary::new();
            library.register(tensor.bind(py).borrow())?;

            let requested = tensor_reference(
                py,
                name,
                Vec::new(),
                PartialStructure::from_logical_slots([
                    euc.slot(PartialIndex::Open(OpenPortId(0))),
                    mink.slot(PartialIndex::Open(OpenPortId(1))),
                ]),
            )?;
            let exact = library.__getitem__(py, requested.bind(py).extract()?)?;
            assert_eq!(
                exact
                    .bind(py)
                    .borrow()
                    .interface
                    .logical_slots()
                    .into_iter()
                    .map(|slot| slot.rep())
                    .collect::<Vec<_>>(),
                vec![euc, mink]
            );

            let by_name = library.__getitem__(
                py,
                ConvertibleToLibraryReference(LibraryReference::Symbol(name)),
            )?;
            assert_eq!(
                by_name
                    .bind(py)
                    .borrow()
                    .interface
                    .logical_slots()
                    .into_iter()
                    .map(|slot| slot.rep())
                    .collect::<Vec<_>>(),
                vec![mink, euc]
            );

            let library = Py::new(py, library)?;
            let indexed = exact.bind(py).call1(("i", "j"))?;
            let network = indexed.call_method1("to_network", (library.clone_ref(py),))?;
            network.call_method1("execute", (library.clone_ref(py),))?;
            let result = network.call_method1("result_tensor", (library,))?;
            let values = (0..6)
                .map(|index| result.get_item(index)?.extract::<f64>())
                .collect::<PyResult<Vec<_>>>()?;
            assert_eq!(values, vec![0., 3., 1., 4., 2., 5.]);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn library_reference_rejects_non_scalar_mixed_interfaces() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let interface = PartialStructure::from_logical_slots([
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(0))),
            euc.slot(PartialIndex::Open(OpenPortId(0))),
        ]);

        assert!(
            ExactLibraryReference::new(
                interface,
                symbol!("spynso_mixed_library_reference"),
                Vec::new(),
            )
            .is_err()
        );
    }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);
