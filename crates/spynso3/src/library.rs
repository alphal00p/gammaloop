use std::collections::HashMap;

use eyre::eyre;
use pyo3::{
    FromPyObject, PyErr, exceptions,
    prelude::*,
    types::{PyAny, PyComplex},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType,
    derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods},
};

use spenso::algebra::complex::Complex;
use spenso::network::library::function_lib::{PanicMissingConcrete, SymbolLib};
use spenso::network::library::{FunctionLibraryError, function_lib::INBUILTS};
use spenso::network::parsing::ShadowedStructure;
use spenso::tensors::complex::RealOrComplexTensor;
use spenso::tensors::data::StorageTensor;
use spenso::{
    network::library::symbolic::{ExplicitKey, TensorLibrary},
    structure::{
        HasStructure, PermutedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialStructure, PartialStructureExt},
        slot::IsAbstractSlot,
    },
    tensors::parametric::MixedTensor,
};
use symbolica::atom::{Atom, DefaultNamespace, FunctionBuilder, SymbolBuilder};
use symbolica::{
    api::python::PythonExpression,
    atom::{AtomView, Symbol},
    function,
};

use crate::{
    Spensor, broadcast::SpensoBroadcastFunction, expression::TensorExpression,
    structure::SpensoName,
};

use super::structure::SpensoStructure;

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
        a.library.insert_scalar_fallible(INBUILTS.conj, |scalar| {
            let Ok(value) = symbolica::domains::float::Complex::<f64>::try_from(scalar.as_view())
            else {
                return Ok(INBUILTS.conj(scalar));
            };
            Ok(Atom::num(value.re) - Atom::num(value.im) * Atom::i())
        });

        a
    }

    /// Register an elementwise callback for concrete tensor execution.
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
            Python::attach(|py| match tensor {
                RealOrComplexTensor::Real(tensor) => tensor
                    .map_data_ref_result(|value| callback.call1(py, (*value,))?.extract::<f64>(py))
                    .map(RealOrComplexTensor::Real),
                RealOrComplexTensor::Complex(tensor) => tensor
                    .map_data_ref_result(|value| {
                        let value = PyComplex::from_doubles(py, value.re, value.im);
                        callback.call1(py, (value,))?.extract::<Complex<f64>>(py)
                    })
                    .map(RealOrComplexTensor::Complex),
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
    key: PermutedStructure<ExplicitKey<AbstractIndex>>,
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
        let rank = interface.structure.order();
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
    Exact(ExactLibraryReference),
}

pub struct ConvertibleToLibraryReference(LibraryReference);

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToLibraryReference {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(expression) = value.extract::<PyRef<'_, TensorExpression>>() {
            return ExactLibraryReference::from_expression(&expression)
                .map(LibraryReference::Exact)
                .map(Self);
        }
        if let Ok(structure) = value.extract::<SpensoStructure>() {
            let expression = TensorExpression::from_structure(value.py(), &structure)?;
            let expression = expression.bind(value.py()).borrow();
            return ExactLibraryReference::from_expression(&expression)
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
    let atom = FunctionBuilder::new(name)
        .add_args(&args)
        .add_args(
            interface
                .logical_slots()
                .into_iter()
                .map(|slot| slot.rep().to_symbolic([])),
        )
        .finish();
    TensorExpression::from_atom_interface_descriptor(py, atom, interface, Some(name), args)
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
            .structure
            .clone()
            .map_structure(|_| reference.key.structure.clone());
        self.library.insert_explicit(PermutedStructure {
            structure: storage,
            rep_permutation: reference.key.rep_permutation.clone(),
            index_permutation: reference.key.index_permutation.clone(),
        });
        self.references
            .insert(reference.key.structure, reference.interface);
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
    ///     An atomic reference with the registered tensor's logical interface
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
                    .get(&reference.key.structure)
                    .map_err(|error| exceptions::PyKeyError::new_err(error.to_string()))?;
                let interface = self
                    .references
                    .get(&reference.key.structure)
                    .cloned()
                    .unwrap_or(reference.interface);
                tensor_reference(py, reference.name, reference.args, interface)
            }
            LibraryReference::Symbol(symbol) => {
                let key = self.library.get_key_from_name(symbol).map_err(|error| {
                    let mut variants = self
                        .references
                        .iter()
                        .filter(|(key, _)| key.global_name == Some(symbol))
                        .map(|(key, interface)| ExactLibraryReference {
                            key: PermutedStructure::identity(key.clone()),
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
                TensorExpression::from_structure(
                    py,
                    &SpensoStructure {
                        structure: PermutedStructure::identity(key),
                    },
                )
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

/// Enumeration for different tensor namespaces in physics.
///
/// Provides categorization for different types of tensor operations and structures
/// commonly used in theoretical physics calculations.
///
/// # Variants:
/// - Weyl: Tensors related to Weyl spinors and chiral representations
/// - Algebra: Tensors related to algebraic structures and Lie algebras
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(from_py_object, eq, eq_int, module = "symbolica.community.spenso")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TensorNamespace {
    Weyl,
    Algebra,
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::structure::{
        dimension::Dimension,
        partial::{OpenPortId, PartialIndex},
        representation::{ExtendibleReps, RepName},
    };
    use symbolica::symbol;

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

        assert_eq!(reference.key.structure.global_name, Some(name));
        assert_eq!(
            reference.key.structure.additional_args,
            Some(vec![Atom::num(7)])
        );
        let canonical = reference.key.structure.reps();
        assert_eq!(
            reference.key.rep_permutation.apply_slice_inv(&canonical),
            vec![mink, euc]
        );
        assert!(reference.signature().contains("7"));
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
