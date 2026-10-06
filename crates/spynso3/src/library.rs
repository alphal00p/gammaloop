use std::{cell::Cell, collections::HashMap};

use eyre::eyre;
use idenso::{IndexTooling, dirac::AGS, tensor::SymbolicTensor};
use pyo3::{
    FromPyObject, PyErr, exceptions,
    prelude::*,
    types::{PyAny, PyComplex, PyIterator, PyList},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType,
    derive::{gen_stub_pyclass, gen_stub_pymethods},
};

use spenso::algebra::complex::{Complex, RealOrComplex};
use spenso::network::library::function_lib::{PanicMissingConcrete, SymbolLib};
use spenso::network::library::{FunctionLibraryError, LibraryError, function_lib::INBUILTS};
use spenso::network::parsing::ShadowedStructure;
use spenso::tensors::complex::RealOrComplexTensor;
use spenso::tensors::data::StorageTensor;
use spenso::{
    network::library::symbolic::{ExplicitKey, TensorLibrary},
    structure::{
        Canonicalized, HasStructure, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
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
    Spensor, TensorDataDescriptor, broadcast::SpensoBroadcastFunction,
    expression::TensorExpression, structure::SpensoName,
};

use super::ModuleInit;

/// Reusable component definitions indexed by tensor name and signature.
///
/// Register a named Tensor, then use its symbolic expression in calculations
/// with this library. A full signature includes the name, scalar arguments, and
/// representations. A name alone is a shortcut when it identifies one entry.
///
/// Lookup returns independent component data. Enumeration lists stored entries;
/// dimension-dependent factories, such as metrics, are resolved on demand.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import Tensor
/// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
/// >>> from symbolica.community.tensor import TensorLibrary
/// >>> library = TensorLibrary()
/// >>> library.register(tensor)
/// >>> library[A.name][1, 0]
/// 3.0
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorLibrary", module = "symbolica.community.tensor")]
pub struct SpensorLibrary {
    pub(crate) library: TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>,
    /// Python-facing interfaces retained for registered and preloaded references.
    references: HashMap<ExplicitKey<AbstractIndex>, PartialStructure>,
}

/// Numerical implementations of unary elementwise tensor functions.
///
/// Use with network execution to evaluate BroadcastFunction calls on numeric
/// components. Complex conjugation is included by default. Symbolic function
/// calls whose arguments are not numeric can remain symbolic.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import TensorFunctionLibrary, BroadcastFunction
/// >>> functions = TensorFunctionLibrary()
/// >>> functions.register(BroadcastFunction("square"), lambda x: x * x)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorFunctionLibrary", module = "symbolica.community.tensor")]
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
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage, on_success)]
#[pymethods]
impl SpensorFunctionLibrary {
    /// Return a summary of the registered elementwise tensor functions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> functions = sp.TensorFunctionLibrary()
    /// >>> text = repr(functions)
    fn __repr__(&self) -> String {
        let mut names = self
            .library
            .functions
            .keys()
            .map(|name| name.to_string())
            .collect::<Vec<_>>();
        names.sort();
        format!("TensorFunctionLibrary(functions={names:?})")
    }

    /// Create a numerical function library with built-in conjugation.
    ///
    /// Returns
    /// -------
    /// TensorFunctionLibrary
    ///     Mutable library for BroadcastFunction implementations.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorFunctionLibrary
    /// >>> functions = TensorFunctionLibrary()
    #[new]
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

    /// Register a numerical implementation for an elementwise function.
    ///
    /// Parameters
    /// ----------
    /// function : BroadcastFunction
    ///     Function to implement; an existing implementation is replaced.
    /// callback : callable
    ///     Receives one float or complex value and returns a float or complex
    ///     value. Results determine the numerical output type.
    ///
    /// Returns
    /// -------
    /// None
    ///     Updates this function library.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorFunctionLibrary, BroadcastFunction
    /// >>> square = BroadcastFunction("numeric_square")
    /// >>> functions = TensorFunctionLibrary()
    /// >>> functions.register(square, lambda value: value * value)
    /// >>> square(tensor).to_tensor(function_library=functions)[0, 1]
    /// 4.0
    pub fn register(
        &mut self,
        function: &SpensoBroadcastFunction,
        #[gen_stub(override_type(
            type_repr = "typing.Callable[[float | complex], float | complex]",
            imports = ("typing",)
        ))]
        callback: Py<PyAny>,
    ) {
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
        Self::new(descriptor.structure().clone(), name, args)
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
    let value = SymbolicTensor::from_signature(&structure)
        .map_err(|error| exceptions::PyRuntimeError::new_err(error.to_string()))?;
    TensorExpression::from_known_parts(py, value.into_expression(), interface, Some(name), args)
}

/// Make a Python-facing interface directly from a key's canonical storage
/// order. Symbol-only lookup returns this interface; exact lookup validates
/// that gamma callers supplied the same order.
fn storage_interface(key: &ExplicitKey<AbstractIndex>) -> PartialStructure {
    PartialStructure::from_logical_slots(
        key.external_reps_iter()
            .enumerate()
            .map(|(position, representation)| representation.slot(PartialIndex::open(position))),
    )
}

fn storage_hep_references() -> HashMap<ExplicitKey<AbstractIndex>, PartialStructure> {
    [
        AGS.gamma_strct(4),
        AGS.gamma_adj_strct(4),
        AGS.gamma_conj_strct(4),
    ]
    .into_iter()
    .map(|key| {
        let interface = storage_interface(key.canonical());
        (key.into_canonical(), interface)
    })
    .collect()
}

impl SpensorLibrary {
    fn stored_reference(&self, key: &ExplicitKey<AbstractIndex>) -> ExactLibraryReference {
        ExactLibraryReference {
            key: Canonicalized::identity(key.clone()),
            interface: self
                .references
                .get(key)
                .cloned()
                .unwrap_or_else(|| storage_interface(key)),
            name: key.global_name.expect("library entries are named"),
            args: key.additional_args.clone().unwrap_or_default(),
        }
    }

    fn stored_references(&self) -> Vec<ExactLibraryReference> {
        let mut references = self
            .library
            .explicit_keys()
            .map(|key| self.stored_reference(key))
            .collect::<Vec<_>>();
        references.sort_by_cached_key(|reference| {
            (reference.name.get_name().to_owned(), reference.signature())
        });
        references
    }

    fn resolve_reference(&self, key: LibraryReference) -> PyResult<Option<ExactLibraryReference>> {
        let reference = match key {
            LibraryReference::Exact(reference) => {
                if [AGS.gamma, AGS.gammaadj, AGS.gammaconj].contains(&reference.name)
                    && let Some(registered) = self.references.get(reference.key.canonical())
                    && registered.logical_slots() != reference.interface.logical_slots()
                {
                    return Err(exceptions::PyKeyError::new_err(format!(
                        "tensor library signature '{}' does not use the registered storage-order interface",
                        reference.signature(),
                    )));
                }
                *reference
            }
            LibraryReference::Symbol(symbol) => {
                let key = match self.library.get_key_from_name(symbol) {
                    Ok(key) => key,
                    Err(LibraryError::InvalidKey | LibraryError::NotFound(_)) => return Ok(None),
                    Err(LibraryError::MultipleKeys(_)) => {
                        let variants = self
                            .stored_references()
                            .into_iter()
                            .filter(|reference| reference.name == symbol)
                            .map(|reference| reference.signature())
                            .collect::<Vec<_>>();
                        return Err(exceptions::PyKeyError::new_err(format!(
                            "tensor name `{symbol}` is ambiguous; registered signatures: {}",
                            variants.join(", ")
                        )));
                    }
                };
                self.stored_reference(&key)
            }
        };
        // Apply the same concrete-dimension validation as data lookup without
        // allocating component data when checking factory membership.
        crate::tensor_data_layout(&reference.interface)?;
        Ok(self
            .library
            .contains_key(reference.key.canonical())
            .then_some(reference))
    }

    fn stored_tensor(&self, reference: ExactLibraryReference) -> PyResult<Spensor> {
        // Rebuild the public logical layout while retaining canonical component storage.
        let key = ExplicitKey::from_iter(
            reference
                .interface
                .logical_slots()
                .into_iter()
                .map(|slot| slot.rep()),
            reference.name,
            (!reference.args.is_empty()).then_some(reference.args.clone()),
        );
        let descriptor = SymbolicTensor::from_signature(&key)
            .map_err(|error| exceptions::PyRuntimeError::new_err(error.to_string()))?;
        let descriptor =
            SymbolicTensor::checked_parts(descriptor.into_expression(), reference.interface)
                .map_err(|error| exceptions::PyRuntimeError::new_err(error.to_string()))?;
        // Validate concrete dimensions before invoking dimension-dependent factories.
        let descriptor = TensorDataDescriptor::new(descriptor, reference.name, reference.args)?;
        let tensor = self
            .library
            .get_storage(key.canonical())
            .map_err(|error| exceptions::PyKeyError::new_err(error.to_string()))?
            .into_owned()
            .map_structure(|_| descriptor.structure.canonical().clone());
        Ok(Spensor::from_storage_with_descriptor(
            tensor,
            descriptor.descriptor,
            Some(descriptor.name),
            descriptor.args,
        ))
    }
}

#[allow(clippy::new_without_default)]
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage, on_success)]
#[pymethods]
impl SpensorLibrary {
    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> library = sp.TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> text = repr(library)
    fn __repr__(&self) -> String {
        format!(
            "TensorLibrary(explicit={}, generic={})",
            self.library.explicit_len(),
            self.library.generic_len()
        )
    }

    /// Display stored tensors with on-demand component and construction details.
    ///
    /// Parameters
    /// ----------
    /// settings : DisplaySettings, optional
    ///     Index and tensor presentation choices; defaults to DisplaySettings().
    ///
    /// Returns
    /// -------
    /// str
    ///     HTML fragment for display in a notebook or page.
    ///
    /// Notes
    /// -----
    /// Stored entries can be selected to inspect the same component explorer
    /// as a standalone Tensor. On-demand factories are listed separately.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> value = library
    /// >>> html = value.to_html()
    #[pyo3(signature = (*, settings=None))]
    fn to_html(
        &self,
        py: Python<'_>,
        settings: Option<crate::display::DisplaySettings>,
    ) -> PyResult<String> {
        crate::display::library::to_html(py, self, &settings.unwrap_or_default())
    }

    fn _repr_html_(&self, py: Python<'_>) -> Option<String> {
        self.to_html(py, None).ok()
    }

    /// Create an empty tensor component library.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> library = sp.TensorLibrary()
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     Library ready to accept named tensors through register().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> len(library)
    /// 0
    #[new]
    pub fn new() -> Self {
        let mut a = Self {
            library: TensorLibrary::new(),
            references: HashMap::new(),
        };
        a.library.update_ids();
        a
    }

    /// Create an empty tensor component library.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> library = sp.TensorLibrary.construct()
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     Library ready to accept named tensors through register().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary.construct()
    /// >>> len(library)
    /// 0
    #[staticmethod]
    pub fn construct() -> Self {
        Self::new()
    }

    /// Store an independent copy of a named tensor's components.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::A")(r, r)
    /// >>> tensor = sp.Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> library = sp.TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> stored = library[A]
    ///
    /// Parameters
    /// ----------
    /// tensor : Tensor
    ///     Tensor with a name, scalar key arguments, concrete dimensions, and
    ///     fully unresolved external axes (no explicit index labels).
    ///     Use tensor.with_name(...) if the tensor is unnamed.
    ///
    /// Returns
    /// -------
    /// None
    ///     Adds or replaces the entry for this exact signature.
    ///
    /// Notes
    /// -----
    /// Changing the original tensor afterward does not change this library.
    /// Register a modified tensor again to update its stored definition.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> len(library)
    /// 1
    pub fn register(&mut self, tensor: PyRef<'_, Spensor>) -> PyResult<()> {
        let name = tensor.descriptor_name.ok_or_else(|| {
            exceptions::PyValueError::new_err(
                "TensorLibrary.register() requires a named Tensor; use tensor.with_name(...) first",
            )
        })?;
        let reference = ExactLibraryReference::new(
            tensor.descriptor.structure().clone(),
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

    /// Retrieve a tensor by exact signature or unambiguous name.
    ///
    /// Parameters
    /// ----------
    /// key : TensorExpression, TensorName, Expression, or str
    ///     A fully unresolved tensor expression selects an exact signature.
    ///     A TensorName, symbol, or string selects a name only when exactly one
    ///     stored signature has that name. Exact signatures also select factories.
    ///     Strings without a namespace use "spenso_python"; prefer TensorName
    ///     or a fully qualified string when the symbol uses another namespace.
    ///
    /// Returns
    /// -------
    /// Tensor
    ///     An independent copy in the requested logical axis order. Call its
    ///     expression() method for a symbolic reference.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     The key is missing or the name matches more than one stored signature.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> stored = library[A]
    /// >>> stored[0, 1]
    /// 2.0
    pub fn __getitem__(
        &self,
        py: Python<'_>,
        key: ConvertibleToLibraryReference,
    ) -> PyResult<Py<Spensor>> {
        let signature = match &key.0 {
            LibraryReference::Symbol(symbol) => symbol.to_string(),
            LibraryReference::Exact(reference) => reference.signature(),
        };
        let reference = self.resolve_reference(key.0)?.ok_or_else(|| {
            exceptions::PyKeyError::new_err(format!(
                "tensor library has no entry for `{signature}`"
            ))
        })?;
        Py::new(py, self.stored_tensor(reference)?)
    }

    /// Check whether a stored signature or factory can resolve a key.
    ///
    /// Parameters
    /// ----------
    /// key : TensorExpression, TensorName, Expression, or str
    ///     A fully unresolved tensor expression selects an exact signature.
    ///     A TensorName, symbol, or string selects a name only when exactly one
    ///     stored signature has that name. Exact signatures also select factories.
    ///     Strings without a namespace use "spenso_python"; prefer TensorName
    ///     or a fully qualified string when the symbol uses another namespace.
    ///
    /// Returns
    /// -------
    /// bool
    ///     False for a missing key. Factories are checked without constructing
    ///     their component data.
    ///
    /// Notes
    /// -----
    /// An ambiguous name or invalid key raises the same error as lookup.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> A in library
    /// True
    fn __contains__(&self, key: ConvertibleToLibraryReference) -> PyResult<bool> {
        Ok(self.resolve_reference(key.0)?.is_some())
    }

    /// Retrieve a tensor, returning a default when its signature is absent.
    ///
    /// Parameters
    /// ----------
    /// key : TensorExpression, TensorName, Expression, or str
    ///     A fully unresolved tensor expression selects an exact signature.
    ///     A TensorName, symbol, or string selects a name only when exactly one
    ///     stored signature has that name. Exact signatures also select factories.
    ///     Strings without a namespace use "spenso_python"; prefer TensorName
    ///     or a fully qualified string when the symbol uses another namespace.
    /// default : object, default None
    ///     Value returned only when the signature is missing.
    ///
    /// Returns
    /// -------
    /// Tensor or object
    ///     Independent tensor snapshot, or the supplied default.
    ///
    /// Notes
    /// -----
    /// Ambiguous names and invalid keys still raise; they do not select the default.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> library.get("missing_tensor") is None
    /// True
    #[pyo3(signature = (key, default=None))]
    fn get(
        &self,
        py: Python<'_>,
        key: ConvertibleToLibraryReference,
        default: Option<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        match self.resolve_reference(key.0)? {
            Some(reference) => Ok(Py::new(py, self.stored_tensor(reference)?)?.into_any()),
            None => Ok(default.unwrap_or_else(|| py.None())),
        }
    }

    /// Count stored tensor definitions.
    ///
    /// Returns
    /// -------
    /// int
    ///     Number of stored entries; excludes on-demand factories.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> len(library)
    /// 1
    fn __len__(&self) -> usize {
        self.library.explicit_len()
    }

    /// List the full unresolved signatures of stored tensors.
    ///
    /// Returns
    /// -------
    /// list of TensorExpression
    ///     Snapshot in name/signature order, excluding on-demand factories.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> len(library.keys())
    /// 1
    pub fn keys(&self, py: Python<'_>) -> PyResult<Vec<Py<TensorExpression>>> {
        self.stored_references()
            .into_iter()
            .map(|reference| {
                tensor_reference(py, reference.name, reference.args, reference.interface)
            })
            .collect()
    }

    /// List independent copies of the stored tensors.
    ///
    /// Returns
    /// -------
    /// list of Tensor
    ///     Snapshot in the same order as keys(), excluding on-demand factories.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> len(library.values())
    /// 1
    pub fn values(&self, py: Python<'_>) -> PyResult<Vec<Py<Spensor>>> {
        self.stored_references()
            .into_iter()
            .map(|reference| Py::new(py, self.stored_tensor(reference)?))
            .collect()
    }

    /// List stored tensor signatures together with their component data.
    ///
    /// Returns
    /// -------
    /// list of (TensorExpression, Tensor)
    ///     Independent snapshots in the same order as keys().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> len(library.items())
    /// 1
    pub fn items(&self, py: Python<'_>) -> PyResult<Vec<(Py<TensorExpression>, Py<Spensor>)>> {
        self.stored_references()
            .into_iter()
            .map(|reference| {
                let key = tensor_reference(
                    py,
                    reference.name,
                    reference.args.clone(),
                    reference.interface.clone(),
                )?;
                let value = Py::new(py, self.stored_tensor(reference)?)?;
                Ok((key, value))
            })
            .collect()
    }

    /// Iterate over a snapshot of stored tensor signatures.
    ///
    /// Returns
    /// -------
    /// iterator of TensorExpression
    ///     Same entries and ordering as keys().
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import TensorLibrary
    /// >>> library = TensorLibrary()
    /// >>> library.register(tensor)
    /// >>> len(list(library))
    /// 1
    #[gen_stub(override_return_type(type_repr = "typing.Iterator[TensorExpression]", imports = ("typing",)))]
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        PyList::new(py, self.keys(py)?)?.try_iter()
    }

    /// Load standard Dirac and SU(3) tensors with floating-point components.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> library = sp.TensorLibrary.hep_lib()
    /// >>> gamma = sp.TensorExpression.dirac_gamma(4)
    /// >>> components = gamma.to_tensor(library)
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     Independent library with four-dimensional Dirac matrices in the Weyl
    ///     basis, SU(3) color tensors, and dimension-dependent metric factories.
    ///
    /// Notes
    /// -----
    /// Values use double-precision real or complex numbers. This explicitly
    /// selects numerical components; omitted HEP libraries and hep_lib_atom()
    /// use exact Symbolica expressions.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorLibrary, TensorExpression
    /// >>> library = TensorLibrary.hep_lib()
    /// >>> library[TensorExpression.dirac_gamma(4)].shape
    /// (4, 4, 4)
    #[staticmethod]
    pub fn hep_lib() -> Self {
        Self {
            library: spenso_hep_lib::hep_lib_su3(),
            references: storage_hep_references(),
        }
    }

    /// Load standard Dirac and SU(3) tensors with exact symbolic components.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> library = sp.TensorLibrary.hep_lib_atom()
    /// >>> gamma = sp.TensorExpression.dirac_gamma(4)
    /// >>> components = gamma.to_tensor(library)
    ///
    /// Returns
    /// -------
    /// TensorLibrary
    ///     Independent library with four-dimensional Dirac matrices in the Weyl
    ///     basis, SU(3) color tensors, and dimension-dependent metric factories.
    ///
    /// Notes
    /// -----
    /// Components are Symbolica Expressions, retaining exact rational and
    /// algebraic constants for symbolic component calculations. This is the
    /// same component convention used when tensor-network evaluation omits a
    /// library. Use hep_lib() to request floating-point HEP data explicitly.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import TensorLibrary, TensorExpression
    /// >>> library = TensorLibrary.hep_lib_atom()
    /// >>> library[TensorExpression.dirac_gamma(4)].shape
    /// (4, 4, 4)
    #[staticmethod]
    pub fn hep_lib_atom() -> Self {
        Self {
            library: spenso_hep_lib::hep_lib_atom(),
            references: storage_hep_references(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use super::*;
    use idenso::{
        color::CS,
        representations::{Bispinor, initialize},
    };
    use pyo3::types::PyFloat;
    use spenso::network::{ExecutionResult, Sequential, SmallestDegree};
    use spenso::structure::{
        OrderedStructure,
        dimension::Dimension,
        partial::{OpenPortId, PartialIndex},
        representation::{ExtendibleReps, LibraryRep, Minkowski, RepName},
        slot::Slot,
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
                .call_method1("color_f", (8,))?
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
            assert_ne!(
                stored.bind(py).borrow().descriptor.expression(),
                &Atom::Zero
            );
            let expression = stored.bind(py).borrow().expression(py)?;
            let indexed = expression.bind(py).call1(("a", "b", "c"))?;
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
                    .descriptor
                    .structure()
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
                    .descriptor
                    .structure()
                    .logical_slots()
                    .into_iter()
                    .map(|slot| slot.rep())
                    .collect::<Vec<_>>(),
                vec![mink, euc]
            );

            let library = Py::new(py, library)?;
            let values = (0..6)
                .map(|index| exact.bind(py).get_item(index)?.extract::<f64>())
                .collect::<PyResult<Vec<_>>>()?;
            assert_eq!(values, vec![0., 3., 1., 4., 2., 5.]);
            let network = exact.bind(py).call1(("i", "j"))?;
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

    #[test]
    fn hep_gamma_family_references_require_storage_order() {
        idenso::representations::initialize();
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let minkowski = LibraryRep::from(Minkowski {}).new_rep(4);
            let bispinor = Bispinor {}.new_rep(4).cast::<LibraryRep>();
            let expected_representations = vec![bispinor, bispinor, minkowski];
            let storage_interface = PartialStructure::from_logical_slots([
                bispinor.slot(PartialIndex::open(0)),
                bispinor.slot(PartialIndex::open(1)),
                minkowski.slot(PartialIndex::open(2)),
            ]);
            let lorentz_first_interface = PartialStructure::from_logical_slots([
                minkowski.slot(PartialIndex::open(0)),
                bispinor.slot(PartialIndex::open(1)),
                bispinor.slot(PartialIndex::open(2)),
            ]);

            for library in [SpensorLibrary::hep_lib(), SpensorLibrary::hep_lib_atom()] {
                for name in [AGS.gamma, AGS.gammaadj, AGS.gammaconj] {
                    for key in [
                        LibraryReference::Symbol(name),
                        LibraryReference::Exact(Box::new(ExactLibraryReference::new(
                            storage_interface.clone(),
                            name,
                            Vec::new(),
                        )?)),
                    ] {
                        let reference =
                            library.__getitem__(py, ConvertibleToLibraryReference(key))?;
                        let reference = reference.bind(py).borrow();
                        let AtomView::Fun(function) = reference.descriptor.expression().as_view()
                        else {
                            panic!("a gamma library reference must remain an atomic tensor")
                        };
                        assert_eq!(function.get_symbol(), name);
                        assert_eq!(
                            function
                                .iter()
                                .map(|argument| {
                                    Slot::<LibraryRep, AbstractIndex>::try_from(argument)
                                        .unwrap()
                                        .rep()
                                })
                                .collect::<Vec<_>>(),
                            expected_representations
                        );
                        assert_eq!(
                            reference.descriptor.structure().logical_slots(),
                            storage_interface.logical_slots()
                        );
                    }

                    let lorentz_first = ExactLibraryReference::new(
                        lorentz_first_interface.clone(),
                        name,
                        Vec::new(),
                    )?;
                    let error = library
                        .__getitem__(
                            py,
                            ConvertibleToLibraryReference(LibraryReference::Exact(Box::new(
                                lorentz_first,
                            ))),
                        )
                        .unwrap_err();
                    assert!(error.to_string().contains("storage-order"));
                }
            }
            Ok(())
        })
        .unwrap();
    }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);
