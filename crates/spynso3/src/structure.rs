use std::ops::Deref;

use idenso::representations::{ColorAdjoint, ColorFundamental, ColorSextet};
use itertools::Itertools;

#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

use pyo3::{
    exceptions::{self, PyIndexError, PyRuntimeError, PyTypeError, PyValueError},
    prelude::*,
    pybacked::PyBackedStr,
    types::{PyList, PyTuple},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    TypeInfo,
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};
use spenso::structure::slot::DualSlotTo;
use spenso::{
    network::{
        library::symbolic::{ETS, ExplicitKey},
        parsing::ShadowedStructure,
    },
    structure::{
        HasName, IndexLess, OrderedStructure, PermutedStructure, TensorStructure, ToSymbolic,
        abstract_index::AbstractIndex,
        dimension::Dimension,
        permuted::Perm,
        representation::{
            Euclidean, ExtendibleReps, LibraryRep, Minkowski, RepName, Representation,
        },
        slot::{IsAbstractSlot, Slot},
    },
    tensors::symbolic::SymbolicTensor,
};
use symbolica::{
    api::python::PythonTransformer,
    atom::{
        Atom, AtomView, DefaultNamespace, FunctionBuilder, NamespacedSymbol, Symbol,
        SymbolAttribute,
    },
    state::Workspace,
    symbol,
    transformer::{Transformer, TransformerState},
    utils::Settable,
};

use symbolica::api::python::{ConvertibleToExpression, PythonExpression};

use thiserror::Error;

use idenso::{
    IndexTooling, color::CS, gamma::AGS, metric::PermuteWithMetric, representations::Bispinor,
};

use super::{ModuleInit, SliceOrIntOrExpanded};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, derive::*, impl_stub_type};

pub struct ConvertibleToSpensoName(pub SpensoName);

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToSpensoName {
    type Error = PyErr;

    fn extract(structure: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(structure) = structure.extract::<SpensoName>() {
            Ok(ConvertibleToSpensoName(structure))
        } else if let Ok(s) = structure.extract::<String>() {
            Ok(ConvertibleToSpensoName(SpensoName::symbol_shorthand(
                s, None, None, None, None, None,
            )?))
        } else if let Ok(s) = structure.extract::<PythonExpression>() {
            if let AtomView::Var(a) = s.as_view() {
                Ok(ConvertibleToSpensoName(SpensoName {
                    name: a.get_symbol(),
                }))
            } else {
                Err(PyTypeError::new_err(
                    "Tensor name cannot be built from non-variable expressions",
                ))
            }
        } else {
            Err(PyTypeError::new_err("Invalid input type for tensor name"))
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ConvertibleToSpensoName {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        SpensoName::type_output() | String::type_output() | PythonExpression::type_output()
    }
}
pub enum SpensoSlotOrArgOrRep {
    Slot(SpensoSlot),
    Arg(PythonExpression),
    Rep(SpensoRepresentation),
}

impl<'a, 'py> FromPyObject<'a, 'py> for SpensoSlotOrArgOrRep {
    type Error = PyErr;

    fn extract(structure: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(structure) = structure.extract::<SpensoSlot>() {
            Ok(SpensoSlotOrArgOrRep::Slot(structure))
        } else if let Ok(s) = structure.extract::<SpensoRepresentation>() {
            Ok(SpensoSlotOrArgOrRep::Rep(s))
        } else if let Ok(s) = structure.extract::<ConvertibleToExpression>() {
            Ok(SpensoSlotOrArgOrRep::Arg(s.to_expression()))
        } else {
            Err(PyTypeError::new_err(
                "Invalid input type for tensor slot, representation, or argument",
            ))
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for SpensoSlotOrArgOrRep {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        SpensoIndices::type_output() | <Vec<SpensoSlot>>::type_output()
    }
}

/// A symbolic name for tensor functions and structures.
///
/// TensorName represents named tensor functions that can be called with indices and arguments
/// to create tensor structures. Names can have various mathematical properties like symmetry,
/// antisymmetry, and custom normalization or printing behavior.
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import TensorName, Slot, Representation
/// >>> T = TensorName("T")
/// >>> symmetric_T = TensorName("S", is_symmetric=True)
/// >>> antisymmetric_T = TensorName("A", is_antisymmetric=True)
/// >>> rep = Representation.cof(3)
/// >>> mu = rep('mu')
/// >>> nu = rep('nu')
/// >>> tensor_structure = T(mu, nu)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorName", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct SpensoName {
    pub name: Symbol,
    // pub args: Vec<Atom>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoName {
    #[new]
    #[pyo3(signature = (name,is_symmetric=None,is_antisymmetric=None,is_cyclesymmetric=None,is_linear=None,custom_normalization=None))]
    /// Create a new tensor name with optional mathematical properties.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     The string name for the tensor function
    /// is_symmetric : bool, optional
    ///     If True, tensor is symmetric under index permutation
    /// is_antisymmetric : bool, optional
    ///     If True, tensor is antisymmetric under index permutation
    /// is_cyclesymmetric : bool, optional
    ///     If True, tensor is symmetric under cyclic permutations
    /// is_linear : bool, optional
    ///     If True, tensor is linear in its arguments
    /// custom_normalization : Transformer, optional
    ///     Custom normalization function (advanced)
    ///
    /// Returns
    /// -------
    /// TensorName
    ///     A new TensorName with the specified properties
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorName
    /// >>> T = TensorName("T")
    /// >>> g = TensorName("g", is_symmetric=True)
    /// >>> F = TensorName("F", is_antisymmetric=True)
    /// >>> D = TensorName("D", is_linear=True)
    fn symbol_shorthand(
        name: String,
        is_symmetric: Option<bool>,
        is_antisymmetric: Option<bool>,
        is_cyclesymmetric: Option<bool>,
        is_linear: Option<bool>,
        custom_normalization: Option<PythonTransformer>,
        // custom_print: Option<PyObject>,
    ) -> PyResult<Self> {
        let namespace = DefaultNamespace {
            namespace: "spenso_python".into(),
            data: "",
            file: "".into(),
            line: 0,
        };
        if is_symmetric.is_none()
            && is_antisymmetric.is_none()
            && is_cyclesymmetric.is_none()
            && is_linear.is_none()
            && custom_normalization.is_none()
        // && custom_print.is_none()
        {
            let id = Symbol::new(namespace.attach_namespace(&name))
                .build()
                .map_err(|e| exceptions::PyTypeError::new_err(e.to_string()))?;

            return Ok(SpensoName {
                name: id,
                // args: vec![],
            });
        }

        let count = (is_symmetric == Some(true)) as u8
            + (is_antisymmetric == Some(true)) as u8
            + (is_cyclesymmetric == Some(true)) as u8;

        if count > 1 {
            Err(exceptions::PyValueError::new_err(
                "Function cannot be both symmetric, antisymmetric or cyclesymmetric",
            ))?;
        }

        let mut opts = vec![];

        if let Some(true) = is_symmetric {
            opts.push(SymbolAttribute::Symmetric);
        }

        if let Some(true) = is_antisymmetric {
            opts.push(SymbolAttribute::Antisymmetric);
        }

        if let Some(true) = is_cyclesymmetric {
            opts.push(SymbolAttribute::Cyclesymmetric);
        }

        if let Some(true) = is_linear {
            opts.push(SymbolAttribute::Linear);
        }

        let name = namespace.attach_namespace(&name);

        let mut symbol = Symbol::new(name).with_attributes(opts);

        if let Some(f) = custom_normalization {
            symbol = symbol.with_normalization_function(Box::new(
                move |input: AtomView<'_>, out: &mut Settable<Atom>| {
                    let _ = Workspace::get_local()
                        .with(|ws| {
                            Transformer::execute_chain(
                                input,
                                &f.chain,
                                ws,
                                &TransformerState::default(),
                                out,
                            )
                        })
                        .unwrap();
                },
            ))
        }

        // if let Some(f) = custom_print {
        //     symbol = symbol.with_print_function(Box::new(
        //         move |input: AtomView<'_>, opts: &PrintOptions| {
        //             Python::with_gil(|py| {
        //                 let kwargs = opts.into_py_dict(py).unwrap();
        //                 f.call(
        //                     py,
        //                     (PythonExpression::from(input.to_owned()),),
        //                     Some(&kwargs),
        //                 )
        //                 .unwrap()
        //                 .extract::<Option<String>>(py)
        //                 .unwrap()
        //             })
        //         },
        //     ))
        // }

        let symbol = symbol
            .build()
            .map_err(|e| exceptions::PyTypeError::new_err(e.to_string()))?;

        Ok(SpensoName {
            name: symbol,
            // args: vec![],
        })
    }

    /// Call the tensor name with arguments to create tensor structures.
    ///
    /// Accepts a mix of slots (for indexed tensors), representations (for indexless tensors),
    /// and symbolic expressions (for additional arguments). Cannot mix slots and representations.
    ///
    /// Parameters
    /// ----------
    /// *args : Slot, Representation, or Expression
    ///     Mixed arguments (Slots for TensorIndices, Representations for TensorStructure, or Expressions for additional arguments)
    ///
    /// Returns
    /// -------
    /// TensorIndices or TensorStructure
    ///     Either TensorIndices (if slots provided) or TensorStructure (if representations provided)
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorName, Slot, Representation
    /// >>> import symbolica as sp
    /// >>> T = TensorName("T")
    /// >>> rep = Representation.euc(3)
    /// >>> mu = rep("mu")
    /// >>> nu = rep("nu")
    /// >>> indexed_tensor = T(mu, nu)
    /// >>> structure_tensor = T(rep, rep)
    /// >>> x = sp.S("x")
    /// >>> tensor_with_args = T(mu, nu, x)
    #[pyo3(signature = (*args))]
    #[gen_stub(skip)]
    fn __call__(&self, args: &Bound<'_, PyTuple>) -> PyResult<PossiblyIndexed> {
        let mut add_args: Vec<Atom> = Vec::new();
        let mut slots: Vec<_> = Vec::new();
        let mut reps: Vec<_> = Vec::new();

        for arg_bound in args.iter() {
            let convertible = arg_bound.extract::<SpensoSlotOrArgOrRep>()?;

            match convertible {
                SpensoSlotOrArgOrRep::Arg(expr) => add_args.push(expr.expr),
                SpensoSlotOrArgOrRep::Slot(slot) => slots.push(slot.slot),
                SpensoSlotOrArgOrRep::Rep(rep) => reps.push(rep.representation),
            }
        }

        let add_args = if add_args.is_empty() {
            None
        } else {
            Some(add_args)
        };

        if slots.is_empty() && reps.is_empty() {
            Err(exceptions::PyValueError::new_err(
                "No slots or representations provided",
            ))
        } else if reps.is_empty() {
            Ok(PossiblyIndexed::Indexed(SpensoIndices {
                structure: ShadowedStructure::<AbstractIndex>::from_iter(
                    slots, self.name, add_args,
                ),
            }))
        } else if slots.is_empty() {
            Ok(PossiblyIndexed::Unindexed(SpensoStructure {
                structure: ExplicitKey::from_iter(reps, self.name, add_args),
            }))
        } else {
            Err(exceptions::PyValueError::new_err(
                "Cannot generate structure with both slots and representations",
            ))
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.name)
    }

    fn __str__(&self) -> String {
        format!("{}", self.name)
    }

    /// Convert the tensor name to a symbolic expression.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     A symbolic Expression representing this tensor name
    ///
    /// Examples
    /// --------
    /// >>> T = TensorName("T")
    /// >>> expr = T.to_expression()
    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(Atom::var(self.name))
    }

    /// Predefined metric tensor name.
    #[staticmethod]
    fn g() -> SpensoName {
        SpensoName { name: ETS.metric }
    }

    /// Predefined musical isomorphism tensor name. This enables dualizing self dual indices.
    #[staticmethod]
    fn flat() -> SpensoName {
        SpensoName { name: ETS.flat }
    }

    /// Predefined gamma matrix name.
    #[staticmethod]
    fn gamma() -> SpensoName {
        SpensoName { name: AGS.gamma }
    }

    /// Predefined gamma5 matrix name.
    #[staticmethod]
    fn gamma5() -> SpensoName {
        SpensoName { name: AGS.gamma5 }
    }

    /// Predefined left chiral projector name.
    #[staticmethod]
    fn projm() -> SpensoName {
        SpensoName { name: AGS.projm }
    }

    /// Predefined right chiral projector name.
    #[staticmethod]
    fn projp() -> SpensoName {
        SpensoName { name: AGS.projp }
    }

    /// Predefined sigma matrix name.
    #[staticmethod]
    fn sigma() -> SpensoName {
        SpensoName { name: AGS.sigma }
    }

    /// Predefined color structure constant name.
    #[staticmethod]
    fn f() -> SpensoName {
        SpensoName { name: CS.f }
    }

    /// Predefined color generator name.
    #[staticmethod]
    fn t() -> SpensoName {
        SpensoName { name: CS.t }
    }
}

/// A tensor structure with abstract indices for symbolic tensor operations.
///
/// TensorIndices represents the index structure of tensors with named abstract indices
/// that can be contracted, manipulated symbolically, and converted to expressions.
/// It maintains both the representation structure and index assignments.
///
/// Examples
/// --------
/// >>> from symbolica.community.spenso import TensorIndices, Representation, TensorName
/// >>> rep = Representation.euc(3)
/// >>> mu = rep('mu')
/// >>> nu = rep('nu')
/// >>> indices = TensorIndices(mu, nu)
/// >>> T = TensorName("T")
/// >>> named_indices = T(mu, nu)
/// >>> expr = named_indices.to_expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorIndices", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct SpensoIndices {
    pub structure: PermutedStructure<ShadowedStructure<AbstractIndex>>,
}

impl Deref for SpensoIndices {
    type Target = ShadowedStructure<AbstractIndex>;

    fn deref(&self) -> &Self::Target {
        &self.structure.structure
    }
}

impl From<ShadowedStructure<AbstractIndex>> for SpensoIndices {
    fn from(value: ShadowedStructure<AbstractIndex>) -> Self {
        SpensoIndices {
            structure: PermutedStructure::identity(value),
        }
    }
}

impl ModuleInit for SpensoIndices {
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_class::<SpensoIndices>()?;
        m.add_class::<SpensoName>()?;
        m.add_class::<SpensoSlot>()?;
        m.add_class::<SpensoStructure>()?;
        m.add_class::<SpensoRepresentation>()?;
        Ok(())
    }
}

pub enum ArithmeticStructure {
    Convertible(ConvertibleToExpression),
    Structure(SpensoIndices),
    Expression(PythonExpression),
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ArithmeticStructure {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        ConvertibleToExpression::type_output()
            | SpensoIndices::type_output()
            | PythonExpression::type_output()
    }
}

impl ArithmeticStructure {
    pub fn to_expression(self) -> PyResult<PythonExpression> {
        match self {
            ArithmeticStructure::Convertible(expr) => Ok(expr.to_expression()),
            ArithmeticStructure::Structure(indices) => indices.to_expression(),
            ArithmeticStructure::Expression(expr) => Ok(expr),
        }
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for ArithmeticStructure {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(ob) = ob.extract::<ConvertibleToExpression>() {
            Ok(ArithmeticStructure::Convertible(ob))
        } else if let Ok(ob) = ob.extract::<SpensoIndices>() {
            Ok(ArithmeticStructure::Structure(ob))
        } else {
            Err(exceptions::PyTypeError::new_err(
                "Only convertible expressions and spenso indices can be used",
            ))
        }
    }
}

pub struct ConvertibleToStructure(pub SpensoIndices);

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToStructure {
    type Error = PyErr;

    fn extract(structure: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(structure) = structure.extract::<SpensoIndices>() {
            Ok(ConvertibleToStructure(structure))
        } else if let Ok(s) = structure.extract::<Vec<SpensoSlot>>() {
            Ok(ConvertibleToStructure(SpensoIndices {
                structure: PermutedStructure::<OrderedStructure>::from_iter(
                    s.into_iter().map(|s| s.slot),
                )
                .map_structure(Into::into),
            }))
        } else {
            Err(PyTypeError::new_err(
                "Internal tensor structure can only be build from TensorIndices or lists of Slots",
            ))
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ConvertibleToStructure {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        SpensoIndices::type_output() | <Vec<SpensoSlot>>::type_output()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoIndices {
    /// Create tensor structure from slots and optional arguments.
    ///
    /// Parameters
    /// ----------
    /// *additional_args : Slot or Expression
    ///     Mixed arguments (Slot objects and Expressions for additional arguments)
    /// name : TensorName, optional
    ///     Optional tensor name to assign to the structure
    ///
    /// Returns
    /// -------
    /// TensorIndices
    ///     A new TensorIndices object
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> from symbolica.community.spenso import TensorIndices, Representation, TensorName
    /// >>> rep = Representation.euc(3)
    /// >>> mu = rep('mu')
    /// >>> nu = rep('nu')
    /// >>> structure = TensorIndices(mu, nu)
    /// >>> x = S('x')
    /// >>> structure_with_args = TensorIndices(mu, nu, x)
    /// >>> T = TensorName("T")
    /// >>> named_structure = TensorIndices(mu, nu, name=T)

    #[new]
    #[pyo3(signature =
           (
           *slots,name=None))]
    #[gen_stub(skip)]
    pub fn from_list(
        slots: &Bound<'_, PyTuple>,
        name: Option<ConvertibleToSpensoName>,
    ) -> PyResult<Self> {
        let mut args = Vec::new();
        let mut true_slots = Vec::new();
        for a in slots {
            if let Ok(s) = a.extract::<SpensoSlot>() {
                true_slots.push(s.slot);
            } else if let Ok(arg) = a.extract::<PythonExpression>() {
                args.push(arg.expr);
            } else {
                return Err(exceptions::PyTypeError::new_err(
                    "Only slots and expressions can be used",
                ));
            }
        }

        let args = if args.is_empty() { None } else { Some(args) };
        let mut a: PermutedStructure<ShadowedStructure<AbstractIndex>> =
            PermutedStructure::<OrderedStructure>::from_iter(true_slots).map_structure(Into::into);
        if let Some(name) = name {
            a.structure.set_name(name.0.name);
        };
        a.structure.additional_args = args;

        Ok(SpensoIndices { structure: a })
    }

    /// Set the tensor name for this structure.
    ///
    /// Parameters
    /// ----------
    /// name : TensorName
    ///     The tensor name to assign
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorIndices, TensorName, Representation
    /// >>> rep = Representation.euc(3)
    /// >>> structure = TensorIndices(rep('mu'), rep('nu'))
    /// >>> T = TensorName("T")
    /// >>> structure.set_name(T)
    fn set_name(&mut self, name: ConvertibleToSpensoName) {
        self.structure.structure.set_name(name.0.name);
    }

    /// Get the tensor name of this structure.
    ///
    /// Returns
    /// -------
    /// TensorName or None
    ///     The tensor name if set, None otherwise
    ///
    /// Examples
    /// --------
    /// >>> name = structure.get_name()
    fn get_name(&self) -> Option<SpensoName> {
        self.structure
            .structure
            .name()
            .map(|a| SpensoName { name: a })
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.structure)
    }

    fn __str__(&self) -> String {
        if let Some(structure) = SymbolicTensor::from_named(&self.structure.structure) {
            let atom = PermutedStructure {
                index_permutation: self.structure.index_permutation.clone(),
                rep_permutation: self.structure.rep_permutation.clone(),
                structure,
            }
            .permute_inds()
            .expression;

            format!("{}", atom)
        } else {
            assert!(self.structure.index_permutation.is_identity());
            assert!(self.structure.rep_permutation.is_identity());
            let args = self
                .structure
                .structure
                .external_structure_iter()
                .map(|r| r.to_atom())
                .join(",");

            format!("({})", args.trim_end())
        }
    }

    /// Convert the tensor indices to a symbolic expression.
    ///
    /// Creates a symbolic representation of the tensor with its indices that can be
    /// used in algebraic manipulations and pattern matching.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     A symbolic Expression representing this indexed tensor
    ///
    /// Raises
    /// ------
    /// RuntimeError
    ///     If the tensor structure has no name
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.spenso import TensorName, Representation
    /// >>> T = TensorName("T")
    /// >>> rep = Representation.euc(3)
    /// >>> mu = rep('mu')
    /// >>> nu = rep('nu')
    /// >>> indices = T(mu, nu)
    /// >>> expr = indices.to_expression()
    fn to_expression(&self) -> PyResult<PythonExpression> {
        if self.structure.structure.name().is_none() {
            return Err(PyRuntimeError::new_err("No name"));
        }

        let atom = PermutedStructure {
            index_permutation: self.structure.index_permutation.clone(),
            rep_permutation: self.structure.rep_permutation.clone(),
            structure: self.structure.structure.clone(),
        }
        .permute_with_metric();

        Ok(atom.into())
    }

    fn __len__(&self) -> usize {
        self.structure.structure.size().unwrap()
    }

    #[gen_stub(skip)]
    fn __getitem__(&self, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        match item {
            SliceOrIntOrExpanded::Int(i) => {
                let out: Vec<_> = self
                    .expanded_index(i.into())
                    .map_err(|s| PyIndexError::new_err(s.to_string()))?
                    .into();

                Ok(Python::attach(|py| out.into_pyobject(py).map(|a| a.unbind()))?.into_any())
            }
            SliceOrIntOrExpanded::Expanded(idxs) => {
                let out: usize = self
                    .flat_index(&idxs)
                    .map_err(|s| PyIndexError::new_err(s.to_string()))?
                    .into();

                Ok(Python::attach(|py| out.into_pyobject(py).map(|a| a.unbind()))?.into_any())
            }
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

                let slice: Result<Vec<Vec<usize>>, _> = range
                    .step_by(step)
                    .map(|i| self.expanded_index(i.into()).map(Vec::<usize>::from))
                    .collect();

                match slice {
                    Ok(slice) => {
                        Ok(
                            Python::attach(|py| slice.into_pyobject(py).map(|a| a.unbind()))?
                                .into_any(),
                        )
                    }
                    Err(e) => Err(PyIndexError::new_err(e.to_string())),
                }
            }
        }
    }

    /// Add this expression to `other`, returning the result.
    pub fn __add__(&self, rhs: ArithmeticStructure) -> PyResult<PythonExpression> {
        let rhs = rhs.to_expression()?;
        Ok((self.to_expression()?.expr.as_ref() + rhs.expr.as_ref()).into())
    }

    /// Add this expression to `other`, returning the result.
    pub fn __radd__(&self, rhs: ArithmeticStructure) -> PyResult<PythonExpression> {
        self.__add__(rhs)
    }

    /// Subtract `other` from this expression, returning the result.
    pub fn __sub__(&self, rhs: ArithmeticStructure) -> PyResult<PythonExpression> {
        let rhs = rhs.to_expression()?.__neg__()?;
        self.__add__(ArithmeticStructure::Expression(rhs))
    }

    /// Subtract this expression from `other`, returning the result.
    pub fn __rsub__(&self, rhs: ArithmeticStructure) -> PyResult<PythonExpression> {
        let s = self.to_expression()?.__neg__()?.expr;

        let r = rhs.to_expression()?.expr;
        Ok((r + s).into())
    }

    /// Add this expression to `other`, returning the result.
    pub fn __mul__(&self, rhs: ArithmeticStructure) -> PyResult<PythonExpression> {
        let rhs = rhs.to_expression()?;
        Ok((self.to_expression()?.expr.as_ref() * rhs.expr.as_ref()).into())
    }

    /// Add this expression to `other`, returning the result.
    pub fn __rmul__(&self, rhs: ArithmeticStructure) -> PyResult<PythonExpression> {
        self.__mul__(rhs)
    }
}

/// A tensor structure without abstract indices, defined purely by representations.
///
/// TensorStructure represents the shape and representation structure of tensors
/// without specific index assignments. It's used for defining tensor templates
/// in libraries and for creating indexless tensor computations.
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import TensorStructure, Representation, TensorName
///
/// # Create from representations
/// rep = Representation.euc(3)
/// structure = TensorStructure(rep, rep)  # 3x3 matrix structure
///
/// # With name for library registration
/// T = TensorName("T")
/// named_structure = TensorStructure(rep, rep, name=T)
///
/// # Use to create indexed tensor
/// indices = structure.index('mu', 'nu')  # Assign specific indices
///
/// # Create symbolic expression
/// expr = structure.symbolic('a', 'b')  # T(a, b)
/// ```
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "TensorStructure", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct SpensoStructure {
    pub structure: PermutedStructure<ExplicitKey<AbstractIndex>>,
}

impl Deref for SpensoStructure {
    type Target = ExplicitKey<AbstractIndex>;

    fn deref(&self) -> &Self::Target {
        &self.structure.structure
    }
}

pub struct ConvertibleToIndexLess(pub SpensoStructure);

impl From<ExplicitKey<AbstractIndex>> for SpensoStructure {
    fn from(value: ExplicitKey<AbstractIndex>) -> Self {
        SpensoStructure {
            structure: PermutedStructure::identity(value),
        }
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToIndexLess {
    type Error = PyErr;

    fn extract(structure: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(structure) = structure.extract::<SpensoStructure>() {
            Ok(ConvertibleToIndexLess(structure))
        } else if let Ok(s) = structure.extract::<Vec<SpensoRepresentation>>() {
            Ok(ConvertibleToIndexLess(SpensoStructure {
                structure: PermutedStructure::<IndexLess>::from_iter(
                    s.into_iter().map(|s| s.representation),
                )
                .map_structure(Into::into),
            }))
        } else if let Ok(s) = structure.extract::<Vec<usize>>() {
            Ok(ConvertibleToIndexLess(SpensoStructure {
                structure: PermutedStructure::<IndexLess>::from_iter(
                    s.into_iter().map(|s| ExtendibleReps::EUCLIDEAN.new_rep(s)),
                )
                .map_structure(Into::into),
            }))
        } else {
            Err(PyTypeError::new_err(
                "Internal tensor structure can only be build from TensorStructure or lists of Representations or Integers",
            ))
        }
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ConvertibleToIndexLess {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        SpensoStructure::type_output()
            | <Vec<SpensoRepresentation>>::type_output()
            | <Vec<usize>>::type_output()
    }
}
#[derive(Error, Debug)]
pub enum SpensoError {
    #[error("Must have a name to register")]
    NoName,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoStructure {
    #[new]
    /// Create tensor structure from representations and optional arguments.
    ///
    /// # Parameters:
    /// - *reps_and_additional_args: Mixed arguments (Representation objects and Expressions for additional arguments)
    /// - name: Optional tensor name to assign to the structure
    ///
    /// # Examples:
    /// ```python
    /// from symbolica import S
    /// from symbolica.community.spenso import TensorStructure, Representation, TensorName
    ///
    /// # Create from representations
    /// rep = Representation.euc(3)
    /// structure = TensorStructure(rep, rep)  # 3x3 tensor
    ///
    /// # With additional arguments
    /// x = S('x')
    /// structure_with_args = TensorStructure(rep, rep, x)
    ///
    /// # With name
    /// T = TensorName("T")
    /// named_structure = TensorStructure(rep, rep, name=T)
    /// ```
    #[pyo3(signature =(*reps_and_additional_args,name=None))]
    #[gen_stub(skip)]
    pub fn from_list(
        reps_and_additional_args: &Bound<'_, PyTuple>,
        name: Option<ConvertibleToSpensoName>,
    ) -> PyResult<Self> {
        let mut args = Vec::new();
        let mut actual_slots = Vec::new();
        for a in reps_and_additional_args {
            if let Ok(s) = a.extract::<SpensoRepresentation>() {
                actual_slots.push(s.representation);
            } else if let Ok(arg) = a.extract::<PythonExpression>() {
                args.push(arg.expr);
            } else {
                return Err(exceptions::PyTypeError::new_err(
                    "Only slots and expressions can be used",
                ));
            }
        }

        let args = if args.is_empty() { None } else { Some(args) };

        let mut a: PermutedStructure<ExplicitKey<AbstractIndex>> =
            PermutedStructure::<IndexLess>::from_iter(actual_slots).map_structure(Into::into);
        if let Some(name) = name {
            a.structure.set_name(name.0.name);
        };
        a.structure.additional_args = args;

        Ok(SpensoStructure { structure: a })
    }

    fn set_name(&mut self, name: ConvertibleToSpensoName) {
        self.structure.structure.set_name(name.0.name);
    }

    fn get_name(&self) -> Option<SpensoName> {
        self.structure
            .structure
            .name()
            .map(|a| SpensoName { name: a })
    }

    fn __repr__(&self) -> String {
        format!(
            "{}",
            self.to_symbolic(Some(self.structure.rep_permutation.clone()))
                .unwrap()
        )
    }

    fn __str__(&self) -> String {
        let slot = self
            .external_reps()
            .into_iter()
            .map(|r| r.to_symbolic([]))
            .join(",");

        match (self.name(), self.args()) {
            (Some(name), Some(args)) => {
                let args = args.iter().join(",");
                format!("{}({})[{}]", name, args, slot)
            }
            (Some(name), None) => {
                format!("{}[{}]", name, slot)
            }
            (None, Some(args)) => {
                let args = args.iter().join(",");
                format!("({})[{}]", args, slot)
            }
            (None, None) => {
                format!("[{}]", slot)
            }
        }
    }

    fn __len__(&self) -> usize {
        self.size().unwrap()
    }

    #[gen_stub(skip)]
    fn __getitem__(&self, item: SliceOrIntOrExpanded) -> PyResult<Py<PyAny>> {
        match item {
            SliceOrIntOrExpanded::Int(i) => {
                let out: Vec<_> = self
                    .expanded_index(i.into())
                    .map_err(|s| PyIndexError::new_err(s.to_string()))?
                    .into();

                Ok(Python::attach(|py| out.into_pyobject(py).map(|a| a.unbind()))?.into_any())
            }
            SliceOrIntOrExpanded::Expanded(idxs) => {
                let out: usize = self
                    .flat_index(&idxs)
                    .map_err(|s| PyIndexError::new_err(s.to_string()))?
                    .into();

                Ok(Python::attach(|py| out.into_pyobject(py).map(|a| a.unbind()))?.into_any())
            }
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

                let slice: Result<Vec<Vec<usize>>, _> = range
                    .step_by(step)
                    .map(|i| self.expanded_index(i.into()).map(Vec::<usize>::from))
                    .collect();

                match slice {
                    Ok(slice) => {
                        Ok(
                            Python::attach(|py| slice.into_pyobject(py).map(|a| a.unbind()))?
                                .into_any(),
                        )
                    }
                    Err(e) => Err(PyIndexError::new_err(e.to_string())),
                }
            }
        }
    }

    // #[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
    /// Convenience method for creating symbolic expressions.
    ///
    /// This is a shorthand for calling `symbolic(*args, extra_args=extra_args)`.
    /// Creates a symbolic Expression representing this tensor structure.
    ///
    /// # Parameters:
    /// - *args: Positional arguments (indices and additional args)
    /// - extra_args: Optional list of additional non-tensorial arguments
    ///
    /// # Examples:
    /// ```python
    /// structure = TensorStructure(rep, rep, name="T")
    /// expr = structure('mu', 'nu')  # Same as structure.symbolic('mu', 'nu')
    /// ```
    #[pyo3(signature = (*args, extra_args=None))]
    #[gen_stub(skip)]
    fn __call__(
        &self,
        // #[gen_stub(override_type(type_repr = "int | float | complex | str |"))]
        args: &Bound<'_, PyTuple>,
        extra_args: Option<&Bound<'_, PyList>>,
    ) -> PyResult<PythonExpression> {
        // Directly delegate to symbolic, passing relevant arguments through
        self.symbolic(args, extra_args)
    }

    /// Create a symbolic expression representing this tensor structure.
    ///
    /// Builds a symbolic tensor expression with the specified indices. Arguments can be
    /// separated using a semicolon (';') to distinguish between additional arguments
    /// and tensor indices.
    ///
    /// # Parameters:
    /// - *args: Positional arguments (int, str, Symbol, Expression for indices, ';' for separator)
    /// - extra_args: Optional list of additional non-tensorial arguments
    ///
    /// # Examples:
    /// ```python
    /// import symbolica as sp
    /// from symbolica.community.spenso import TensorStructure, Representation, TensorName
    ///
    /// rep = Representation.euc(3)
    /// T = TensorName("T")
    /// structure = TensorStructure([rep, rep], name=T)
    ///
    /// # Basic usage
    /// expr = structure.symbolic('mu', 'nu')  # T(mu, nu)
    ///
    /// # With additional arguments
    /// x = sp.S('x')
    /// expr = structure.symbolic(x, ';', 'mu', 'nu')  # T(x; mu, nu)
    ///
    /// # Using extra_args parameter
    /// expr = structure.symbolic('mu', 'nu', extra_args=[x])  # T(x; mu, nu)
    /// ```
    #[pyo3(signature = (*args, extra_args=None))]
    #[gen_stub(skip)]
    fn symbolic(
        &self,
        args: &Bound<'_, PyTuple>,
        extra_args: Option<&Bound<'_, PyList>>,
    ) -> PyResult<PythonExpression> {
        // Use helper to parse arguments
        let (final_additional_args, potential_indices) =
            self.parse_args_for_indexing(args, extra_args)?;

        // --- Generate Symbolic Expression ---
        let name = self.name().ok_or_else(|| {
            PyRuntimeError::new_err("Cannot create symbolic atom: structure has no name")
        })?;

        let index_atoms: Vec<Atom> = potential_indices
            .iter()
            .map(|item| {
                match item {
                    // potential_indices now only contains Aind or Atom
                    ConvertibleToAbstractIndex::Aind(idx) => (*idx).into(),
                    ConvertibleToAbstractIndex::Atom(expr) => expr.expr.clone(),
                    ConvertibleToAbstractIndex::Separator => unreachable!(), // Helper ensures this
                }
            })
            .collect();

        if self.order() != index_atoms.len() {
            return Err(PyValueError::new_err(format!(
                "Number of index atoms {} does not match structure order {}",
                index_atoms.len(),
                self.order()
            )));
        }

        let slots_atoms = self
            .external_reps_iter()
            .zip(index_atoms)
            .map(|(rep, ind_atom)| rep.to_symbolic([ind_atom]))
            .collect::<Vec<_>>();

        let value_builder = FunctionBuilder::new(name);
        let final_expr = value_builder
            .add_args(&final_additional_args)
            .add_args(&slots_atoms)
            .finish();

        Ok(PythonExpression::from(final_expr))
    }

    /// Create an indexed tensor (TensorIndices) from this structure.
    ///
    /// Converts this structure template into a concrete indexed tensor by assigning
    /// specific abstract indices to each representation slot.
    ///
    /// # Parameters:
    /// - *args: Positional arguments (indices and ';' separator between additional args and indices)
    /// - extra_args: Optional list of additional non-tensorial arguments
    /// - cook_indices: If True, attempt to convert expressions to valid indices
    ///
    /// # Examples:
    /// ```python
    /// import symbolica as sp
    /// from symbolica.community.spenso import TensorStructure, Representation, TensorName
    ///
    /// rep = Representation.cof(3)
    /// T = TensorName("T")
    /// structure = TensorStructure([rep, rep], name=T)
    ///
    /// # Create indexed tensor
    /// indices = structure.index('mu', 'nu')  # T with indices mu, nu
    ///
    /// # With additional arguments
    /// x = sp.S('x')
    /// indices = structure.index(x, ';', 'mu', 'nu')  # T(x; mu, nu)
    /// ```
    #[pyo3(signature = (*args, extra_args=None, cook_indices=false))]
    #[gen_stub(skip)]
    fn index(
        &self,
        args: &Bound<'_, PyTuple>,
        extra_args: Option<&Bound<'_, PyList>>,
        cook_indices: bool,
    ) -> PyResult<SpensoIndices> {
        // Use helper to parse arguments
        let (final_additional_args, potential_indices) =
            self.parse_args_for_indexing(args, extra_args)?;

        // --- Resolve Indices (No change in this logic) ---
        let mut resolved_indices: Vec<AbstractIndex> = Vec::new();
        for item in potential_indices {
            // potential_indices now only contains Aind or Atom
            match item {
                ConvertibleToAbstractIndex::Aind(idx) => {
                    resolved_indices.push(idx);
                }
                ConvertibleToAbstractIndex::Atom(expr) => {
                    let converted_atom: Result<AbstractIndex, _> = if cook_indices {
                        expr.expr.cook_indices().as_view().try_into()
                    } else {
                        expr.expr.as_view().try_into()
                    };
                    match converted_atom {
                        Ok(idx) => resolved_indices.push(idx),
                        Err(e) => {
                            let cook_msg = if cook_indices {
                                ""
                            } else {
                                " Try setting cook_indices=True."
                            };
                            return Err(exceptions::PyValueError::new_err(format!(
                                "Cannot convert argument '{}' to an AbstractIndex: {}. Ensure it's a valid index type or cookable.{}",
                                expr.expr, e, cook_msg
                            )));
                        }
                    }
                }
                ConvertibleToAbstractIndex::Separator => unreachable!(), // Helper ensures this
            }
        }

        let mut structure_clone = self.structure.structure.clone();
        structure_clone.additional_args = if final_additional_args.is_empty() {
            None
        } else {
            Some(final_additional_args)
        };
        match structure_clone.reindex(&resolved_indices) {
            Ok(indexed_structure) => Ok(SpensoIndices {
                structure: indexed_structure,
            }),
            Err(e) => Err(PyValueError::new_err(format!(
                "Failed to create TensorIndices: {}",
                e
            ))),
        }
    }
}

impl SpensoStructure {
    fn parse_args_for_indexing(
        &self,
        args: &Bound<'_, PyTuple>,
        extra_args_opt: Option<&Bound<'_, PyList>>,
    ) -> PyResult<(Vec<Atom>, Vec<ConvertibleToAbstractIndex>)> {
        let mut pre_separator_args: Vec<ConvertibleToAbstractIndex> = Vec::new();
        let mut post_separator_args: Vec<ConvertibleToAbstractIndex> = Vec::new();
        let mut separator_found = false;

        for arg_bound in args.iter() {
            let convertible = arg_bound.extract::<ConvertibleToAbstractIndex>()?;

            match convertible {
                ConvertibleToAbstractIndex::Separator => {
                    if separator_found {
                        return Err(exceptions::PyValueError::new_err(
                            "Separator token ';' used more than once.",
                        ));
                    }

                    separator_found = true;
                    pre_separator_args.append(&mut post_separator_args);
                }
                item => {
                    post_separator_args.push(item);
                }
            }
        }

        let mut final_additional_args = self.args().unwrap_or_default();
        for item in pre_separator_args {
            match item {
                ConvertibleToAbstractIndex::Aind(idx) => final_additional_args.push(idx.into()),
                ConvertibleToAbstractIndex::Atom(expr) => {
                    final_additional_args.push(expr.expr.clone())
                }
                ConvertibleToAbstractIndex::Separator => unreachable!(),
            }
        }
        if let Some(extra_args_list) = extra_args_opt {
            for item_bound in extra_args_list.iter() {
                let expr = item_bound.extract::<PythonExpression>()?;
                final_additional_args.push(expr.expr);
            }
        }

        Ok((final_additional_args, post_separator_args))
    }
}
#[derive(IntoPyObject)]
pub enum PossiblyIndexed {
    Unindexed(SpensoStructure),
    Indexed(SpensoIndices),
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for PossiblyIndexed {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        SpensoStructure::type_output() | SpensoIndices::type_output()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Representation", module = "symbolica.community.spenso")]
#[derive(Clone)]
/// A representation in the sense of group representation theory for tensor indices.
///
/// Representations define the transformation properties of tensor indices under group operations.
/// They specify the dimension and duality structure, determining which indices can contract.
///
/// Key concepts:
/// - **Self-dual**: Indices can contract with other indices of the same representation
/// - **Dualizable**: Indices can only contract with their dual representation
/// - **Dimension**: Size of the representation space
///
/// Predefined representations are available as class methods:
/// - `Representation.euc(d)`: Euclidean space (self-dual)
/// - `Representation.mink(d)`: Minkowski space (self-dual)
/// - `Representation.bis(d)`: Bispinor (self-dual)
/// - `Representation.cof(d)`: Color fundamental (dualizable)
/// - `Representation.coad(d)`: Color adjoint (self-dual)
/// - `Representation.cos(d)`: Color sextet (dualizable)
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import Representation
///
/// # Standard representations
/// euclidean = Representation.euc(4)      # 4D Euclidean
/// lorentz = Representation.mink(4)       # 4D Minkowski
/// color = Representation.cof(3)          # SU(3) fundamental
/// adjoint = Representation.coad(8)       # SU(3) adjoint
///
/// # Custom representation
/// custom = Representation("MyRep", 5, is_self_dual=True)
///
/// # Create slots with indices
/// mu_slot = euclidean('mu')              # Euclidean index μ
/// a_slot = color('a')                    # Color index a
///
/// # Generate metric tensors
/// metric = euclidean.g('mu', 'nu')       # g_μν
/// ```
///
pub struct SpensoRepresentation {
    pub representation: Representation<LibraryRep>,
}

pub enum ConvertibleToAbstractIndex {
    Aind(AbstractIndex),
    Atom(PythonExpression),
    Separator,
}

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToAbstractIndex {
    type Error = PyErr;

    fn extract(aind: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let aind = if let Ok(i) = aind.extract::<char>() {
            if i == ';' {
                ConvertibleToAbstractIndex::Separator
            } else {
                let mut tmp = [0u8; 4];
                let name = i.encode_utf8(&mut tmp);
                ConvertibleToAbstractIndex::Aind(AbstractIndex::Symbol(symbol!(&name).into()))
            }
        } else if let Ok(i) = aind.extract::<isize>() {
            ConvertibleToAbstractIndex::Aind(i.into())
        } else if let Ok(expr) = aind.extract::<PythonExpression>() {
            match expr.expr.as_view() {
                AtomView::Var(v) => {
                    ConvertibleToAbstractIndex::Aind(AbstractIndex::Symbol(v.get_symbol().into()))
                }
                _ => ConvertibleToAbstractIndex::Atom(expr),
            }
        } else if let Ok(s) = aind.extract::<PyBackedStr>() {
            let id = symbol!(&s);
            ConvertibleToAbstractIndex::Aind(AbstractIndex::Symbol(id.into()))
        } else {
            return Err(PyTypeError::new_err(
                "Argument must be convertible to an index (int, str, Symbol), an Expression,, or the separator ';'",
            ));
        };

        Ok(aind)
    }
}

#[cfg(feature = "python_stubgen")]
impl_stub_type!(ConvertibleToAbstractIndex = isize | Symbol | PyBackedStr);

pub struct ConvertibleToDimension(Dimension);

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToDimension {
    type Error = PyErr;

    fn extract(dimension: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let dim = if let Ok(i) = dimension.extract::<usize>() {
            Dimension::from(i)
        } else if let Ok(expr) = dimension.extract::<PythonExpression>() {
            let id = match expr.expr.as_view() {
                AtomView::Var(v) => v.get_symbol(),
                _ => {
                    return Err(exceptions::PyTypeError::new_err(
                        "Only symbols can be abstract indices",
                    ));
                }
            };
            Dimension::from(id)
        } else if let Ok(s) = dimension.extract::<PyBackedStr>() {
            let ns = "spenso_python";
            let id = Symbol::new(NamespacedSymbol {
                symbol: format!("{}::{}", ns, s).into(),
                namespace: ns.into(),
                file: file!().into(),
                line: line!() as usize,
            })
            .build()
            .unwrap();

            Dimension::from(id)
        } else {
            return Err(PyTypeError::new_err(
                "dimension must be an non-zero integer or a symbol",
            ));
        };
        Ok(ConvertibleToDimension(dim))
    }
}

#[cfg(feature = "python_stubgen")]
impl_stub_type!(ConvertibleToDimension = usize | PythonExpression | PyBackedStr);

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[pymethods]
impl SpensoRepresentation {
    #[new]
    #[pyo3(signature =(name,dimension,is_self_dual=true))]
    /// Create and register a new representation with specified properties.
    ///
    /// # Parameters:
    /// - name: String name for the representation
    /// - dimension: Size of the representation (int or symbolic)
    /// - is_self_dual: If True, creates self-dual representation; if False, creates dualizable pair
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import Representation
    /// import symbolica as sp
    ///
    /// # Self-dual representation (indices contract with themselves)
    /// euclidean = Representation("Euclidean", 4, is_self_dual=True)
    ///
    /// # Dualizable representation (needs dual partner for contraction)
    /// vector_up = Representation("VectorUp", 4, is_self_dual=False)
    ///
    /// # Symbolic dimension
    /// n = sp.S('n')
    /// general = Representation("General", n, is_self_dual=True)
    /// ```
    pub fn register_new(
        name: String,
        dimension: ConvertibleToDimension,
        is_self_dual: bool,
    ) -> PyResult<Self> {
        let dim = dimension.0;

        let rep = if is_self_dual {
            LibraryRep::new_self_dual(&name).unwrap().new_rep(dim)
        } else {
            LibraryRep::new_dual(&name).unwrap().new_rep(dim)
        };
        Ok(SpensoRepresentation {
            representation: rep,
        })
    }

    fn dual(&self) -> Self {
        Self {
            representation: self.representation.dual(),
        }
    }

    /// Create a slot or symbolic expression from this representation.
    ///
    /// # Parameters:
    /// - aind: The index specification (Abstract index for Slot, Expression for symbolic representation)
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import Representation
    /// import symbolica as sp
    ///
    /// rep = Representation.euc(3)
    ///
    /// # Create slots with different index types
    /// slot1 = rep('mu')        # String index
    /// slot2 = rep(1)           # Integer index
    /// slot3 = rep(sp.S('nu'))  # Symbolic index
    ///
    /// # Create symbolic expression
    /// x = sp.S('x')
    /// sym_rep = rep(x)         # Symbolic representation
    /// ```
    #[gen_stub(skip)]
    fn __call__(&self, py: Python<'_>, aind: ConvertibleToAbstractIndex) -> PyResult<Py<PyAny>> {
        match aind {
            ConvertibleToAbstractIndex::Separator => {
                Err(PyValueError::new_err("separator cannot be an index"))
            }
            ConvertibleToAbstractIndex::Aind(aind) => Ok(SpensoSlot {
                slot: self.representation.slot(aind),
            }
            .into_pyobject(py)
            .map(|a| a.unbind())?
            .into_any()),
            ConvertibleToAbstractIndex::Atom(a) => {
                let a: PythonExpression = self.representation.to_symbolic([a.expr]).into();

                Ok(a.into_pyobject(py).map(|a| a.unbind())?.into_any())
            }
        }
    }

    /// Create a metric tensor for this representation.
    ///
    /// # Parameters:
    /// - i: First index
    /// - j: Second index
    ///
    /// # Examples:
    /// ```python
    /// rep = Representation.mink(4)
    /// metric = rep.g('mu', 'nu')  # Minkowski metric g_μν
    /// ```
    fn g(
        &self,
        i: ConvertibleToAbstractIndex,
        j: ConvertibleToAbstractIndex,
    ) -> PyResult<SpensoIndices> {
        match (i, j) {
            (ConvertibleToAbstractIndex::Aind(i), ConvertibleToAbstractIndex::Aind(j)) => {
                let structure = ShadowedStructure::<AbstractIndex>::from_iter(
                    [self.representation.slot(i), self.representation.slot(j)],
                    ETS.metric,
                    None,
                );

                Ok(SpensoIndices { structure })
            }
            _ => Err(PyValueError::new_err("indices must be abstract indices")),
        }
    }

    /// Create a musical isomorphism tensor for this representation.
    ///
    /// # Parameters:
    /// - i: First index
    /// - j: Second index
    ///
    /// # Examples:
    /// ```python
    /// rep = Representation.mink(4)
    /// flat = rep.flat('mu', 'nu')  # Flat isomorphism ♭_μν
    /// ```
    fn flat(
        &self,
        i: ConvertibleToAbstractIndex,
        j: ConvertibleToAbstractIndex,
    ) -> PyResult<SpensoIndices> {
        match (i, j) {
            (ConvertibleToAbstractIndex::Aind(i), ConvertibleToAbstractIndex::Aind(j)) => {
                let structure = ShadowedStructure::<AbstractIndex>::from_iter(
                    [self.representation.slot(i), self.representation.slot(j)],
                    ETS.flat,
                    None,
                );

                Ok(SpensoIndices { structure })
            }
            _ => Err(PyValueError::new_err("indices must be abstract indices")),
        }
    }

    /// Create an identity tensor for this representation.
    ///
    /// # Parameters:
    /// - i: First index
    /// - j: Second index
    ///
    /// # Examples:
    /// ```python
    /// rep = Representation.cof(3)
    /// identity = rep.id('a', 'b')  # Color identity δ_ab
    /// ```
    fn id(
        &self,
        i: ConvertibleToAbstractIndex,
        j: ConvertibleToAbstractIndex,
    ) -> PyResult<SpensoIndices> {
        match (i, j) {
            (ConvertibleToAbstractIndex::Aind(i), ConvertibleToAbstractIndex::Aind(j)) => {
                let structure = ShadowedStructure::<AbstractIndex>::from_iter(
                    [
                        self.representation.slot(i).dual(),
                        self.representation.slot(j),
                    ],
                    ETS.metric,
                    None,
                );

                Ok(SpensoIndices { structure })
            }
            _ => Err(PyValueError::new_err("indices must be abstract indices")),
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.representation)
    }

    fn __str__(&self) -> String {
        format!("{}", self.representation.to_symbolic([]))
    }

    /// Convert the representation to a symbolic expression.
    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(self.representation.to_symbolic([]))
    }

    /// Create a bispinor representation.
    ///
    /// # Parameters:
    /// - dimension: The dimension of the bispinor space
    #[staticmethod]
    fn bis(dimension: ConvertibleToDimension) -> Self {
        let dim = dimension.0;
        let rep = Bispinor {}.new_rep(dim).cast();
        Self {
            representation: rep,
        }
    }

    /// Create a Euclidean space representation.
    ///
    /// # Parameters:
    /// - dimension: The dimension of the Euclidean space
    #[staticmethod]
    fn euc(dimension: ConvertibleToDimension) -> Self {
        let dim = dimension.0;
        let rep = Euclidean {}.new_rep(dim).cast();
        Self {
            representation: rep,
        }
    }

    /// Create a Minkowski space representation.
    ///
    /// # Parameters:
    /// - dimension: The dimension of the Minkowski space
    #[staticmethod]
    fn mink(dimension: ConvertibleToDimension) -> Self {
        let dim = dimension.0;
        let rep = Minkowski {}.new_rep(dim).cast();
        Self {
            representation: rep,
        }
    }

    /// Create a color fundamental representation.
    ///
    /// # Parameters:
    /// - dimension: The dimension of the color group (e.g., 3 for SU(3))
    #[staticmethod]
    fn cof(dimension: ConvertibleToDimension) -> Self {
        let dim = dimension.0;
        let rep = ColorFundamental {}.new_rep(dim).cast();
        Self {
            representation: rep,
        }
    }

    /// Create a color adjoint representation.
    ///
    /// # Parameters:
    /// - dimension: The dimension of the adjoint representation (e.g., 8 for SU(3))
    #[staticmethod]
    fn coad(dimension: ConvertibleToDimension) -> Self {
        let dim = dimension.0;
        let rep = ColorAdjoint {}.new_rep(dim).cast();
        Self {
            representation: rep,
        }
    }

    /// Create a color sextet representation.
    ///
    /// # Parameters:
    /// - dimension: The dimension of the sextet representation (e.g., 6 for SU(3))
    #[staticmethod]
    fn cos(dimension: ConvertibleToDimension) -> Self {
        let dim = dimension.0;
        let rep = ColorSextet {}.new_rep(dim).cast();
        Self {
            representation: rep,
        }
    }
}

/// A tensor index slot combining a representation with an abstract index.
///
/// Slots are the building blocks for tensor structures, pairing a representation
/// (which defines transformation properties) with an abstract index identifier.
/// Slots with matching representations and indices can be contracted.
///
/// # Examples:
/// ```python
/// from symbolica.community.spenso import Slot, Representation
/// import symbolica as sp
///
/// # Create representation and slots
/// rep = Representation.euc(3)
/// slot1 = rep('mu')           # Slot with string index
/// slot2 = rep(1)              # Slot with integer index
/// slot3 = rep(sp.S('nu'))     # Slot with symbolic index
///
/// # Create custom slot
/// custom_slot = Slot("MyRep", 4, 'alpha', dual=False)
///
/// # Use in tensor structures
/// from symbolica.community.spenso import TensorIndices
/// tensor_structure = TensorIndices(slot1, slot2)
/// ```
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Slot", module = "symbolica.community.spenso")]
#[derive(Clone)]
pub struct SpensoSlot {
    pub slot: Slot<LibraryRep>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl SpensoSlot {
    fn __repr__(&self) -> String {
        format!("{:?}", self.slot)
    }

    fn __str__(&self) -> String {
        format!("{}", self.slot.to_atom())
    }

    fn dual(&self) -> Self {
        SpensoSlot {
            slot: self.slot.dual(),
        }
    }

    #[new]
    #[pyo3(signature =(name,dimension,aind,dual=false))]
    /// Create a new slot with a custom representation and index.
    ///
    /// # Parameters:
    /// - name: String name for the representation
    /// - dimension: Size of the representation space
    /// - aind: The abstract index (int, str, or Symbol)
    /// - dual: If True, creates dualizable representation; if False, self-dual
    ///
    /// # Examples:
    /// ```python
    /// from symbolica.community.spenso import Slot
    /// import symbolica as sp
    ///
    /// # Self-dual slot
    /// euclidean_slot = Slot("Euclidean", 4, 'mu', dual=False)
    ///
    /// # Dualizable slot
    /// vector_slot = Slot("Vector", 4, 'nu', dual=True)
    ///
    /// # With symbolic index
    /// sym_index = sp.S('alpha')
    /// symbolic_slot = Slot("Custom", 3, sym_index, dual=False)
    /// ```
    pub fn register_new(
        name: String,
        dimension: usize,
        aind: ConvertibleToAbstractIndex,
        dual: bool,
    ) -> PyResult<Self> {
        let rep = if dual {
            LibraryRep::new_dual(&name).unwrap().new_rep(dimension)
        } else {
            LibraryRep::new_self_dual(&name).unwrap().new_rep(dimension)
        };

        match aind {
            ConvertibleToAbstractIndex::Aind(a) => {
                let mut slot = rep.slot(a);
                if dual {
                    slot = slot.dual();
                }
                Ok(SpensoSlot { slot })
            }
            ConvertibleToAbstractIndex::Atom(a) => match a.expr.as_view() {
                AtomView::Var(v) => {
                    let slot = rep.slot(AbstractIndex::Symbol(v.get_symbol().into()));
                    Ok(SpensoSlot { slot })
                }
                _ => Err(exceptions::PyTypeError::new_err(
                    "Only symbols can be abstract indices",
                )),
            },
            ConvertibleToAbstractIndex::Separator => {
                Err(PyValueError::new_err("separator cannot be an index"))
            }
        }
    }

    /// Convert the slot to a symbolic expression.
    ///
    /// # Examples:
    /// ```python
    /// rep = Representation.euc(3)
    /// slot = rep('mu')
    /// expr = slot.to_expression()  # Symbolic representation of the slot
    /// ```
    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(self.slot.to_atom())
    }
}

#[cfg(feature = "python_stubgen")]
static EMPTY: fn() -> String = || "[]".into();

#[cfg(feature = "python_stubgen")]
static FALSE: fn() -> String = || "False".to_string();

#[cfg(feature = "python_stubgen")]
submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<SpensoName>,
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
                        name: "args",//: Slot | Expression | int | str | float | complex",
                        kind:ParameterKind::VarPositional,
                        default:ParameterDefault::None,
                        type_info: || SpensoSlot::type_input()| ConvertibleToExpression::type_input(),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: SpensoIndices::type_output,
                doc:r##"Call the tensor name with arguments to create tensor structures.

Accepts a mix of slots and symbolic expressions (for additional arguments).

Parameters
----------
*args : Slot or Expression
    Slot objects and Expressions for additional arguments

Returns
-------
TensorIndices
    A new TensorIndices object

Examples
--------
>>> from symbolica.community.spenso import TensorName, Slot, Representation
>>> import symbolica as sp
>>> T = TensorName("T")
>>> rep = Representation.euc(3)
>>> mu = rep("mu")
>>> nu = rep("nu")
>>> indexed_tensor = T(mu, nu)
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
                        name: "args",//: Representation | Expression",
                        kind:ParameterKind::VarPositional,
                        default:ParameterDefault::None,
                        type_info: || SpensoRepresentation::type_input()| PythonExpression::type_input(),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: SpensoStructure::type_output,
                doc:r##"Call the tensor name with arguments to create a TensorStructure.

Accepts a mix of representations and symbolic expressions (for additional arguments).

Parameters
----------
*args : Representation or Expression
    Representation objects and Expressions for additional arguments

Returns
-------
TensorStructure
    A new TensorStructure object

Examples
--------
>>> from symbolica.community.spenso import TensorName, Slot, Representation
>>> import symbolica as sp
>>> T = TensorName("T")
>>> rep = Representation.euc(3)
>>> structure_tensor = T(rep, rep)
"##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            }
        ]
    }
}

#[cfg(feature = "python_stubgen")]
static NONE_ARG: fn() -> String = || "None".into();

#[cfg(feature = "python_stubgen")]
submit! {
PyMethodsInfo {
    struct_id: std::any::TypeId::of::<SpensoIndices>,
    attrs: &[],
    getters: &[],
    setters: &[],
    file: file!(),
    line: line!(),
    column: column!(),
    methods: &[
        MethodInfo{
            name:"__new__",
            parameters:&[
                ParameterInfo{
                    name:"slots",//: Slot | Expression | int | str | float | complex",
                    type_info:ConvertibleToStructure::type_input,
                    kind:ParameterKind::VarPositional,
                    default:ParameterDefault::None,
                },
                ParameterInfo{
                    name:"name",
                    kind:ParameterKind::KeywordOnly,
                    default:ParameterDefault::Expr(NONE_ARG),
                    type_info:||ConvertibleToSpensoName::type_input() | TypeInfo::none(),
                },
            ]
            ,r#type: MethodType::New,
            r#return: SpensoIndices::type_output,
            doc:r##"Create tensor structure from slots and optional arguments.

Parameters
----------
*additional_args : Slot or Expression
    Mixed arguments (Slot objects and Expressions for additional arguments)
name : TensorName, optional
    Optional tensor name to assign to the structure

Returns
-------
TensorIndices
    A new TensorIndices object

Examples
--------
>>> from symbolica import S
>>> from symbolica.community.spenso import TensorIndices, Representation, TensorName
>>> rep = Representation.euc(3)
>>> mu = rep('mu')
>>> nu = rep('nu')
>>> structure = TensorIndices(mu, nu)
>>> x = S('x')
>>> structure_with_args = TensorIndices(mu, nu, x)
>>> T = TensorName("T")
>>> named_structure = TensorIndices(mu, nu, name=T)
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
                    default: ParameterDefault::None,
                    kind: ParameterKind::PositionalOrKeyword,
                    type_info: || TypeInfo::builtin("slice"),
                },
            ],
            r#type: MethodType::Instance,
            r#return: Vec::<crate::TensorElements>::type_output,
            doc:r##"Get expanded indices at the specified range of flattened indices.

Parameters
----------
item : slice
    Slice object defining the range of indices

Returns
-------
list of list of int
    List of expanded indices
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
                    default: ParameterDefault::None,
                    kind: ParameterKind::PositionalOrKeyword,
                    type_info:  Vec::<usize>::type_input
                },
            ],
            r#type: MethodType::Instance,
            r#return: crate::TensorElements::type_output,
            doc:r##"Get flattened index associated to this expanded index.

Parameters
----------
item : list of int
    Multi-dimensional index coordinates

Returns
-------
int
    The flat index
"##,
            is_async: false,
            deprecated: None,
            type_ignored: None,
            is_overload: true,
        }, MethodInfo {
            name: "__getitem__",
            parameters: &[
                ParameterInfo {
                    name: "item",
                    default: ParameterDefault::None,
                    kind: ParameterKind::PositionalOrKeyword,
                    type_info: usize::type_input,

                },
            ],
            r#type: MethodType::Instance,
            r#return: crate::TensorElements::type_output,
            doc:r##"Get expanded index associated to this flat index.

Parameters
----------
item : int
    Flat index into the tensor

Returns
-------
list of int
    Multi-dimensional index coordinates
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
submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<SpensoStructure>,
        attrs: &[],
        getters: &[],
        setters: &[],
        file: file!(),
        line: line!(),
        column: column!(),
        methods: &[
            MethodInfo{
                name:"__new__",
                parameters:&[
                    ParameterInfo{
                        name:"reps_and_additional_args",//: Representation |  Expression | int | str | float | complex",
                        type_info:ConvertibleToStructure::type_input,
                        kind:ParameterKind::VarPositional
                        ,default:ParameterDefault::None,
                    },
                    ParameterInfo{
                        name:"name",
                        kind:ParameterKind::KeywordOnly,
                        default:ParameterDefault::Expr(NONE_ARG),
                        type_info:||ConvertibleToSpensoName::type_input() | TypeInfo::none()
                    },
                ]
                ,r#type: MethodType::New,
                r#return: SpensoStructure::type_output,
                doc:r##"Construct a new TensorStructure with the given representations.

Parameters
----------
*reps_and_additional_args : Representation or Expression
    Mixed arguments (Representation objects and Expressions for additional arguments)
name : TensorName, optional
    Optional tensor name to assign to the structure

Returns
-------
TensorStructure
    A new TensorStructure object

Examples
--------
>>> from symbolica import S
>>> from symbolica.community.spenso import TensorStructure, Representation, TensorName
>>> rep = Representation.euc(3)
>>> structure = TensorStructure(rep, rep)
>>> x = S('x')
>>> structure_with_args = TensorStructure(rep, rep, x)
>>> T = TensorName("T")
>>> named_structure = TensorStructure(rep, rep, name=T)
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
                        default: ParameterDefault::None,
                        kind: ParameterKind::PositionalOrKeyword,
                        type_info: || TypeInfo::builtin("slice"),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: Vec::<crate::TensorElements>::type_output,
                doc:r##"Get expanded indices at the specified range of flattened indices.

Parameters
----------
item : slice
    Slice object defining the range of indices

Returns
-------
list of list of int
    List of expanded indices
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
                        default: ParameterDefault::None,
                        kind: ParameterKind::PositionalOrKeyword,
                        type_info:  Vec::<usize>::type_input
                    },
                ],
                r#type: MethodType::Instance,
                r#return: crate::TensorElements::type_output,
                doc:r##"Get flattened index associated to this expanded index.

Parameters
----------
item : list of int
    Multi-dimensional index coordinates

Returns
-------
int
    The flat index
"##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            }, MethodInfo {
                name: "__getitem__",
                parameters: &[
                    ParameterInfo {
                        name: "item",
                        default: ParameterDefault::None,
                        kind: ParameterKind::PositionalOrKeyword,
                        type_info:  usize::type_input
                    },
                ],
                r#type: MethodType::Instance,
                r#return: crate::TensorElements::type_output,
                doc:r##"Get expanded index associated to this flat index.

Parameters
----------
item : int
    Flat index into the tensor

Returns
-------
list of int
    Multi-dimensional index coordinates
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
                        name: "args",//: builtins.int | Expression | str",
                        kind:ParameterKind::VarPositional,
                        default:ParameterDefault::None,
                        type_info:|| ConvertibleToAbstractIndex::type_input(),
                    },
                    ParameterInfo {
                        name: "extra_args",
                        kind:ParameterKind::KeywordOnly,
                        default:ParameterDefault::Expr(EMPTY),
                        type_info: Vec::<ConvertibleToExpression>::type_input,
                    },
                ],
                r#type: MethodType::Instance,
                r#return: PythonExpression::type_output,
                doc:r##"Convenience method for creating symbolic expressions.

This is a shorthand for calling `symbolic(*args, extra_args=extra_args)`.
Creates a symbolic Expression representing this tensor structure.

Parameters
----------
*args : int, str, Symbol, or Expression
    Positional arguments (indices and additional args)
extra_args : list of Expression, optional
    Optional list of additional non-tensorial arguments

Returns
-------
Expression
    A symbolic Expression representing the tensor

Examples
--------
>>> structure = TensorStructure(rep, rep, name="T")
>>> expr = structure('mu', 'nu')
"##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "symbolic",
                parameters: &[
                    ParameterInfo {

                        name: "args",//: builtins.int | Expression | str",
                        kind:ParameterKind::VarPositional,
                        default:ParameterDefault::None,
                        type_info:|| ConvertibleToAbstractIndex::type_input(),
                    },
                    ParameterInfo {

                        name: "extra_args",
                        kind:ParameterKind::KeywordOnly,
                        default:ParameterDefault::Expr(EMPTY),
                        type_info: Vec::<ConvertibleToExpression>::type_input,
                    },
                ],
                r#type: MethodType::Instance,
                r#return: PythonExpression::type_output,
                doc:r##"Create a symbolic expression representing this tensor structure.

Builds a symbolic tensor expression with the specified indices. Arguments can be
separated using a semicolon (';') to distinguish between additional arguments
and tensor indices.

Parameters
----------
*args : int, str, Symbol, Expression, or ';'
    Positional arguments (int, str, Symbol, Expression for indices, ';' for separator)
extra_args : list of Expression, optional
    Optional list of additional non-tensorial arguments

Returns
-------
Expression
    A symbolic Expression representing the tensor with indices

Examples
--------
>>> import symbolica as sp
>>> from symbolica.community.spenso import TensorStructure, Representation, TensorName
>>> rep = Representation.euc(3)
>>> T = TensorName("T")
>>> structure = TensorStructure([rep, rep], name=T)
>>> expr = structure.symbolic('mu', 'nu')
>>> x = sp.S('x')
>>> expr = structure.symbolic(x, ';', 'mu', 'nu')
>>> expr = structure.symbolic('mu', 'nu', extra_args=[x])
"##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },

            MethodInfo {
                name: "index",
                parameters: &[
                    ParameterInfo {
                        name: "args",//: builtins.int | Expression | str",
                        kind:ParameterKind::VarPositional,
                        default:ParameterDefault::None,
                        type_info: ConvertibleToAbstractIndex::type_input,
                    },
                    ParameterInfo {
                        name: "extra_args",
                        kind:ParameterKind::KeywordOnly,
                        default:ParameterDefault::Expr(EMPTY),
                        type_info: std::vec::Vec::<PythonExpression>::type_input,
                    },
                    ParameterInfo {
                        name: "cook_indices",
                        kind:ParameterKind::KeywordOnly,
                        default:ParameterDefault::Expr(FALSE),
                        type_info: bool::type_input
                    }
                ],
                r#type: MethodType::Instance,
                r#return: SpensoIndices::type_output,
                doc:r##"Create an indexed tensor (TensorIndices) from this structure.

Converts this structure template into a concrete indexed tensor by assigning
specific abstract indices to each representation slot.

Parameters
----------
*args : int, str, Symbol, Expression, or ';'
    Positional arguments (indices and ';' separator between additional args and indices)
extra_args : list of Expression, optional
    Optional list of additional non-tensorial arguments
cook_indices : bool, optional
    If True, attempt to convert expressions to valid indices

Returns
-------
TensorIndices
    A TensorIndices object with concrete index assignments

Examples
--------
>>> import symbolica as sp
>>> from symbolica.community.spenso import TensorStructure, Representation, TensorName
>>> rep = Representation.cof(3)
>>> T = TensorName("T")
>>> structure = TensorStructure([rep, rep], name=T)
>>> indices = structure.index('mu', 'nu')
>>> x = sp.S('x')
>>> indices = structure.index(x, ';', 'mu', 'nu')
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
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);
