use idenso::representations::{ColorAdjoint, ColorFundamental, ColorSextet};

#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

use pyo3::{
    PyTypeInfo,
    exceptions::{self, PyTypeError, PyValueError},
    prelude::*,
    pybacked::PyBackedStr,
    types::{PyAny, PyTuple},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};
use spenso::structure::slot::DualSlotTo;
use spenso::{
    network::{library::symbolic::ETS, parsing::ShadowedStructure, tags::SPENSO_TAG},
    structure::{
        PermutedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{Euclidean, LibraryRep, Minkowski, RepName, Representation},
        slot::{IsAbstractSlot, Slot},
    },
};
use symbolica::{
    api::python::{PythonNormalization, PythonUserData},
    atom::{
        Atom, AtomView, DefaultNamespace, FunctionBuilder, NamespacedSymbol, Symbol, SymbolBuilder,
    },
    symbol,
};

use symbolica::api::python::{ConvertibleToExpression, PythonExpression};

use idenso::{color::CS, dirac::AGS, representations::Bispinor};

use super::expression::{TensorExpression, validate_color_structure_ports};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, derive::*, impl_stub_type};

pub struct ConvertibleToSpensoName(pub SpensoName, pub Vec<Atom>);

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToSpensoName {
    type Error = PyErr;

    fn extract(structure: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(structure) = structure.extract::<SpensoName>() {
            Ok(ConvertibleToSpensoName(structure, Vec::new()))
        } else if let Ok(expression) = structure.extract::<PyRef<'_, TensorExpression>>() {
            if !expression.interface.structure.is_scalar() {
                return Err(PyTypeError::new_err(
                    "a TensorExpression used as a name must have rank zero",
                ));
            }
            let AtomView::Fun(function) = expression.as_super().expr.as_view() else {
                return Err(PyTypeError::new_err(
                    "a TensorExpression used as a name must be an atomic tensor call",
                ));
            };
            if !function.get_symbol().has_tag(&SPENSO_TAG.tensor) {
                return Err(PyTypeError::new_err(
                    "a TensorExpression used as a name must be an atomic tensor call",
                ));
            }
            let mut args = Vec::new();
            for argument in function.iter() {
                if Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok()
                    || Representation::<LibraryRep>::try_from(argument).is_ok()
                {
                    return Err(PyTypeError::new_err(
                        "a TensorExpression used as a name cannot have structural ports",
                    ));
                }
                args.push(argument.to_owned());
            }
            Ok(ConvertibleToSpensoName(
                SpensoName {
                    name: function.get_symbol(),
                },
                args,
            ))
        } else if let Ok(s) = structure.extract::<String>() {
            Ok(ConvertibleToSpensoName(
                SpensoName::symbol_shorthand(
                    structure.py(),
                    s,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                )?,
                Vec::new(),
            ))
        } else if let Ok(s) = structure.extract::<PythonExpression>() {
            if let AtomView::Var(a) = s.as_view() {
                Ok(ConvertibleToSpensoName(
                    SpensoName {
                        name: a.get_symbol(),
                    },
                    Vec::new(),
                ))
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
        SpensoName::type_output()
            | String::type_output()
            | PythonExpression::type_output()
            | TensorExpression::type_output()
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
        } else if let Ok(s) = structure.extract::<PyRef<'_, TensorExpression>>() {
            if !s.interface.structure.is_scalar() {
                return Err(PyTypeError::new_err(
                    "tensor key arguments must be scalar expressions",
                ));
            }
            Ok(SpensoSlotOrArgOrRep::Arg(PythonExpression {
                expr: s.as_super().expr.clone(),
            }))
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
        SpensoSlot::type_output()
            | SpensoRepresentation::type_output()
            | ConvertibleToExpression::type_output()
    }
}

/// A symbolic name for tensor expressions.
///
/// TensorName represents named tensor functions that can be called with scalar arguments, slots,
/// and representations to create tensor expressions. Names can have various mathematical properties like symmetry,
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
/// >>> tensor_expression = T(mu, nu)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "TensorName",
    module = "symbolica.community.spenso"
)]
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
    #[pyo3(signature = (name, *, rank=None, is_symmetric=None, is_antisymmetric=None, is_cyclesymmetric=None, is_linear=None, is_flat=None, is_scalar=None, is_real=None, is_integer=None, is_positive=None, tags=None, aliases=None, normalization=None, print=None, derivative=None, series=None, eval=None, data=None))]
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
    /// rank : int, optional
    ///     The declared rank. Only rank one has a dedicated construction invariant.
    /// tags : list[str], optional
    ///     Extra Symbolica tags. The Spenso tensor tag is always included.
    /// normalization, print, derivative, series, eval, data : optional
    ///     Symbolica symbol callbacks and metadata.
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
    #[allow(clippy::too_many_arguments)]
    fn symbol_shorthand(
        py: Python<'_>,
        name: String,
        rank: Option<usize>,
        is_symmetric: Option<bool>,
        is_antisymmetric: Option<bool>,
        is_cyclesymmetric: Option<bool>,
        is_linear: Option<bool>,
        is_flat: Option<bool>,
        is_scalar: Option<bool>,
        is_real: Option<bool>,
        is_integer: Option<bool>,
        is_positive: Option<bool>,
        tags: Option<Vec<String>>,
        aliases: Option<Vec<String>>,
        normalization: Option<PythonNormalization>,
        print: Option<Py<PyAny>>,
        derivative: Option<Py<PyAny>>,
        series: Option<Py<PyAny>>,
        eval: Option<Py<PyAny>>,
        data: Option<PythonUserData>,
    ) -> PyResult<Self> {
        let rank_one = match rank {
            None => false,
            Some(1) => true,
            Some(rank) => {
                return Err(PyValueError::new_err(format!(
                    "TensorName only supports rank=None or rank=1, got rank={rank}"
                )));
            }
        };

        let mut tags = tags.unwrap_or_default();
        if !tags.iter().any(|tag| tag == &SPENSO_TAG.tensor) {
            tags.push(SPENSO_TAG.tensor.clone());
        }
        if rank_one && !tags.iter().any(|tag| tag == &SPENSO_TAG.rank1) {
            tags.push(SPENSO_TAG.rank1.clone());
        }

        let namespace = DefaultNamespace {
            namespace: "spenso_python".into(),
            data: "",
            file: "".into(),
            line: 0,
        };
        let name = namespace.attach_namespace(&name).symbol.to_string();
        let names = PyTuple::new(py, [name])?;
        let expression_type = PythonExpression::type_object(py);

        // Symbolica's public constructor owns all callback adapters, including
        // normalization and custom printing, so their behavior stays identical to S().
        let expression = PythonExpression::symbol(
            &expression_type,
            py,
            &names,
            is_symmetric,
            is_antisymmetric,
            is_cyclesymmetric,
            is_linear,
            is_flat,
            is_scalar,
            is_real,
            is_integer,
            is_positive,
            Some(tags),
            aliases,
            normalization,
            print,
            derivative,
            series,
            eval,
            data,
        )?
        .extract::<PythonExpression>(py)?;

        let AtomView::Var(name) = expression.as_view() else {
            unreachable!("a single Symbolica symbol constructor result is a variable")
        };
        Ok(SpensoName {
            name: name.get_symbol(),
        })
    }

    /// Create a rank-one tensor name.
    #[staticmethod]
    #[pyo3(signature = (name, *, is_symmetric=None, is_antisymmetric=None, is_cyclesymmetric=None, is_linear=None, is_flat=None, is_scalar=None, is_real=None, is_integer=None, is_positive=None, tags=None, aliases=None, normalization=None, print=None, derivative=None, series=None, eval=None, data=None))]
    #[allow(clippy::too_many_arguments)]
    fn vector(
        py: Python<'_>,
        name: String,
        is_symmetric: Option<bool>,
        is_antisymmetric: Option<bool>,
        is_cyclesymmetric: Option<bool>,
        is_linear: Option<bool>,
        is_flat: Option<bool>,
        is_scalar: Option<bool>,
        is_real: Option<bool>,
        is_integer: Option<bool>,
        is_positive: Option<bool>,
        tags: Option<Vec<String>>,
        aliases: Option<Vec<String>>,
        normalization: Option<PythonNormalization>,
        print: Option<Py<PyAny>>,
        derivative: Option<Py<PyAny>>,
        series: Option<Py<PyAny>>,
        eval: Option<Py<PyAny>>,
        data: Option<PythonUserData>,
    ) -> PyResult<Self> {
        Self::symbol_shorthand(
            py,
            name,
            Some(1),
            is_symmetric,
            is_antisymmetric,
            is_cyclesymmetric,
            is_linear,
            is_flat,
            is_scalar,
            is_real,
            is_integer,
            is_positive,
            tags,
            aliases,
            normalization,
            print,
            derivative,
            series,
            eval,
            data,
        )
    }

    /// Call the tensor name with scalar key arguments followed by structural ports.
    ///
    /// Slots become explicit ports and representations become unresolved ports. They may be
    /// mixed in one call, but every scalar key argument must precede the first structural port.
    ///
    /// Parameters
    /// ----------
    /// *args : Slot, Representation, or Expression
    ///     Scalar expressions followed by Slot and/or Representation ports
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A structured expression, including for calls with no structural ports
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
    /// >>> tensor_with_args = T(x, mu, nu)
    #[pyo3(signature = (*args))]
    #[gen_stub(skip)]
    fn __call__(
        &self,
        py: Python<'_>,
        args: &Bound<'_, PyTuple>,
    ) -> PyResult<Py<TensorExpression>> {
        let mut scalar_args = Vec::new();
        let mut ports = Vec::new();
        let mut port_atoms = Vec::new();
        let rank_one = self.name.has_tag(&SPENSO_TAG.rank1);
        let mut structural_seen = false;
        let mut next_open = 0;

        for arg_bound in args.iter() {
            let convertible = arg_bound.extract::<SpensoSlotOrArgOrRep>()?;

            match convertible {
                SpensoSlotOrArgOrRep::Arg(expr) => {
                    if structural_seen {
                        return Err(PyValueError::new_err(
                            "tensor scalar arguments must precede every slot or representation",
                        ));
                    }
                    scalar_args.push(expr.expr);
                }
                SpensoSlotOrArgOrRep::Slot(slot) => {
                    structural_seen = true;
                    ports.push(
                        slot.slot
                            .rep()
                            .slot(PartialIndex::Explicit(slot.slot.aind())),
                    );
                    port_atoms.push(slot.slot.to_atom());
                }
                SpensoSlotOrArgOrRep::Rep(rep) => {
                    structural_seen = true;
                    ports.push(rep.representation.slot(PartialIndex::open(next_open)));
                    port_atoms.push(rep.representation.to_symbolic([]));
                    next_open += 1;
                }
            }
        }

        if rank_one && ports.len() != 1 {
            return Err(PyValueError::new_err(format!(
                "rank-one tensor names require exactly one slot or representation, got {}",
                ports.len()
            )));
        }
        validate_color_structure_ports(
            self.name,
            args.len(),
            &ports.iter().map(|port| port.rep()).collect::<Vec<_>>(),
        )?;

        let atom = FunctionBuilder::new(self.name)
            .add_args(&scalar_args)
            .add_args(&port_atoms)
            .finish();
        TensorExpression::from_atom_interface_descriptor(
            py,
            atom,
            PartialStructure::from_logical_slots(ports),
            Some(self.name),
            scalar_args,
        )
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

    /// Check whether this tensor name carries `tag`.
    fn has_tag(&self, tag: &str) -> bool {
        self.name.has_tag(tag)
            || (!tag.contains("::")
                && self
                    .name
                    .get_tags()
                    .iter()
                    .any(|candidate| candidate.strip_prefix("python::") == Some(tag)))
    }

    /// Return all Symbolica tags carried by this tensor name.
    fn get_tags(&self) -> Vec<String> {
        self.name.get_tags().to_vec()
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

/// Internal indexed structure used while constructing a `TensorExpression`.
#[derive(Clone)]
pub(crate) struct SpensoIndices {
    pub(crate) structure: PermutedStructure<ShadowedStructure<AbstractIndex>>,
}

pub enum ArithmeticStructure {
    Tensor(Py<TensorExpression>),
    Convertible(ConvertibleToExpression),
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ArithmeticStructure {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        ConvertibleToExpression::type_output() | TensorExpression::type_output()
    }
}

impl ArithmeticStructure {
    pub fn to_expression(self) -> PyResult<PythonExpression> {
        match self {
            ArithmeticStructure::Tensor(expression) => Python::attach(|py| {
                let expression = expression.bind(py).borrow();
                Ok(PythonExpression {
                    expr: TensorExpression::materialized_atom(&expression)?,
                })
            }),
            ArithmeticStructure::Convertible(expr) => Ok(expr.to_expression()),
        }
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for ArithmeticStructure {
    type Error = PyErr;

    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(ob) = ob.extract::<Py<TensorExpression>>() {
            Ok(ArithmeticStructure::Tensor(ob))
        } else if let Ok(ob) = ob.extract::<ConvertibleToExpression>() {
            Ok(ArithmeticStructure::Convertible(ob))
        } else {
            Err(exceptions::PyTypeError::new_err(
                "expected a TensorExpression or Symbolica expression",
            ))
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    eq,
    name = "Representation",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, PartialEq, Eq)]
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
            let id = SymbolBuilder::new(NamespacedSymbol {
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

struct ConvertibleToInvariantDegree(Atom);

impl ConvertibleToInvariantDegree {
    fn quadratic() -> Self {
        Self(Atom::num(2))
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToInvariantDegree {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        value
            .extract::<ConvertibleToExpression>()
            .map(|degree| Self(degree.to_expression().expr))
    }
}

impl<'py> IntoPyObject<'py> for ConvertibleToInvariantDegree {
    type Target = PythonExpression;
    type Output = Bound<'py, PythonExpression>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        PythonExpression::from(self.0).into_pyobject(py)
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ConvertibleToInvariantDegree {
    fn type_input() -> pyo3_stub_gen::TypeInfo {
        ConvertibleToExpression::type_input()
    }

    fn type_output() -> pyo3_stub_gen::TypeInfo {
        PythonExpression::type_output()
    }
}

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

    /// Return the dimension carried by this representation.
    #[getter]
    fn dimension(&self) -> PythonExpression {
        self.representation.dim.to_symbolic().into()
    }

    /// Build the degree-k Casimir eigenvalue for this representation.
    #[pyo3(
        signature = (degree = ConvertibleToInvariantDegree::quadratic()),
        text_signature = "($self, degree=2)"
    )]
    fn casimir(&self, degree: ConvertibleToInvariantDegree) -> PythonExpression {
        CS.cas(degree.0, self.representation.to_symbolic([])).into()
    }

    /// Build the degree-k Dynkin index for this representation.
    #[pyo3(
        signature = (degree = ConvertibleToInvariantDegree::quadratic()),
        text_signature = "($self, degree=2)"
    )]
    fn dynkin_index(&self, degree: ConvertibleToInvariantDegree) -> PythonExpression {
        CS.idx(degree.0, self.representation.to_symbolic([])).into()
    }

    /// Build a degree-k Gram invariant with another representation.
    ///
    /// If `other` is omitted, both sides use this representation.
    #[pyo3(signature = (degree, other = None))]
    fn gram(
        &self,
        degree: ConvertibleToInvariantDegree,
        other: Option<&SpensoRepresentation>,
    ) -> PythonExpression {
        let other = other.unwrap_or(self);
        CS.gram(
            degree.0,
            self.representation.to_symbolic([]),
            other.representation.to_symbolic([]),
        )
        .into()
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
        py: Python<'_>,
        i: ConvertibleToAbstractIndex,
        j: ConvertibleToAbstractIndex,
    ) -> PyResult<Py<TensorExpression>> {
        match (i, j) {
            (ConvertibleToAbstractIndex::Aind(i), ConvertibleToAbstractIndex::Aind(j)) => {
                let structure = ShadowedStructure::<AbstractIndex>::from_iter(
                    [self.representation.slot(i), self.representation.slot(j)],
                    ETS.metric,
                    None,
                );

                TensorExpression::from_indices(py, &SpensoIndices { structure })
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
        py: Python<'_>,
        i: ConvertibleToAbstractIndex,
        j: ConvertibleToAbstractIndex,
    ) -> PyResult<Py<TensorExpression>> {
        match (i, j) {
            (ConvertibleToAbstractIndex::Aind(i), ConvertibleToAbstractIndex::Aind(j)) => {
                let structure = ShadowedStructure::<AbstractIndex>::from_iter(
                    [self.representation.slot(i), self.representation.slot(j)],
                    ETS.flat,
                    None,
                );

                TensorExpression::from_indices(py, &SpensoIndices { structure })
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
        py: Python<'_>,
        i: ConvertibleToAbstractIndex,
        j: ConvertibleToAbstractIndex,
    ) -> PyResult<Py<TensorExpression>> {
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

                TensorExpression::from_indices(py, &SpensoIndices { structure })
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
/// from symbolica.community.spenso import TensorName
/// tensor_expression = TensorName("T")(slot1, slot2)
/// ```
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    eq,
    name = "Slot",
    module = "symbolica.community.spenso"
)]
#[derive(Clone, PartialEq, Eq)]
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
                        name: "args",
                        kind:ParameterKind::VarPositional,
                        default:ParameterDefault::None,
                        type_info: || SpensoSlot::type_input()
                            | SpensoRepresentation::type_input()
                            | ConvertibleToExpression::type_input(),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: TensorExpression::type_output,
                doc:r##"Call the tensor name with arguments to create a tensor expression.

Accepts scalar key expressions followed by any mix of slots and representations.

Parameters
----------
*args : Expression, Slot, or Representation
    Scalar key expressions followed by structural ports

Returns
-------
TensorExpression
    A structured expression with explicit and/or unresolved ports

Examples
--------
>>> from symbolica.community.spenso import TensorName, Slot, Representation
>>> import symbolica as sp
>>> T = TensorName("T")
>>> rep = Representation.euc(3)
>>> mu = rep("mu")
>>> tensor = T(mu, rep)
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
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::structure::representation::ExtendibleReps;

    #[test]
    fn tensor_name_builds_scalar_and_mixed_structured_expressions() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4));
            let name = SpensoName {
                name: SPENSO_TAG.tensor_symbol("python_tensor_name_construction"),
            };
            let arguments = PyTuple::new(
                py,
                [
                    PythonExpression::from(Atom::num(11))
                        .into_pyobject(py)?
                        .unbind()
                        .into_any(),
                    Py::new(
                        py,
                        SpensoSlot {
                            slot: representation.slot(AbstractIndex::Normal(41)),
                        },
                    )?
                    .into_any(),
                    Py::new(py, SpensoRepresentation { representation })?.into_any(),
                ],
            )?;
            let expression = name.__call__(py, &arguments)?;
            let expression_ref = expression.bind(py).borrow();
            assert_eq!(expression_ref.interface.structure.order(), 2);
            assert!(matches!(
                expression_ref.interface.logical_slots()[0].aind,
                PartialIndex::Explicit(AbstractIndex::Normal(41))
            ));
            assert!(matches!(
                expression_ref.interface.logical_slots()[1].aind,
                PartialIndex::Open(_)
            ));
            assert_eq!(expression_ref.name_args, vec![Atom::num(11)]);
            drop(expression_ref);

            let scalar_arguments = PyTuple::new(
                py,
                [PythonExpression::from(Atom::num(13))
                    .into_pyobject(py)?
                    .unbind()
                    .into_any()],
            )?;
            let scalar = name.__call__(py, &scalar_arguments)?;
            assert!(scalar.bind(py).borrow().interface.structure.is_scalar());
            let descriptor = scalar
                .bind(py)
                .as_any()
                .extract::<ConvertibleToSpensoName>()?;
            assert_eq!(descriptor.0.name, name.name);
            assert_eq!(descriptor.1, vec![Atom::num(13)]);

            let invalid_arguments = PyTuple::new(
                py,
                [
                    Py::new(py, SpensoRepresentation { representation })?.into_any(),
                    PythonExpression::from(Atom::num(17))
                        .into_pyobject(py)?
                        .unbind()
                        .into_any(),
                ],
            )?;
            assert!(name.__call__(py, &invalid_arguments).is_err());

            let vector = SpensoName {
                name: spenso::vector_symbol!("python_vector_name_construction"),
            };
            assert!(vector.__call__(py, &PyTuple::empty(py)).is_err());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn representation_builds_scalar_invariants_from_its_metadata() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation_type = py.get_type::<SpensoRepresentation>();
            let fundamental = representation_type.call_method1("cof", (3,))?;
            let adjoint = representation_type.call_method1("coad", (8,))?;
            let euclidean = representation_type.call_method1("euc", (4,))?;

            let fundamental_representation = fundamental
                .extract::<SpensoRepresentation>()?
                .representation;
            let adjoint_representation = adjoint.extract::<SpensoRepresentation>()?.representation;
            let euclidean_representation =
                euclidean.extract::<SpensoRepresentation>()?.representation;

            let inspect_signature = py.import("inspect")?.getattr("signature")?;
            for (method, expected) in [
                ("casimir", "(degree=2)"),
                ("dynkin_index", "(degree=2)"),
                ("gram", "(degree, other=None)"),
            ] {
                let signature = inspect_signature
                    .call1((fundamental.getattr(method)?,))?
                    .str()?
                    .to_str()?
                    .to_owned();
                assert_eq!(signature, expected);
            }

            let dimension = fundamental
                .getattr("dimension")?
                .extract::<PythonExpression>()?;
            assert_eq!(dimension.expr, Atom::num(3));

            let casimir = fundamental.call_method0("casimir")?;
            assert!(!casimir.is_instance_of::<TensorExpression>());
            assert_eq!(
                casimir.extract::<PythonExpression>()?.expr,
                CS.cas(Atom::num(2), fundamental_representation.to_symbolic([]))
            );

            let degree = Atom::var(symbol!("python_color_invariant_degree"));
            let dynkin_index = fundamental
                .call_method1("dynkin_index", (PythonExpression::from(degree.clone()),))?;
            assert!(!dynkin_index.is_instance_of::<TensorExpression>());
            assert_eq!(
                dynkin_index.extract::<PythonExpression>()?.expr,
                CS.idx(degree.clone(), fundamental_representation.to_symbolic([]))
            );

            let self_gram = fundamental.call_method1("gram", (3,))?;
            assert!(!self_gram.is_instance_of::<TensorExpression>());
            assert_eq!(
                self_gram.extract::<PythonExpression>()?.expr,
                CS.gram(
                    Atom::num(3),
                    fundamental_representation.to_symbolic([]),
                    fundamental_representation.to_symbolic([]),
                )
            );

            let mixed_gram = fundamental
                .call_method1("gram", (PythonExpression::from(degree.clone()), &adjoint))?;
            assert_eq!(
                mixed_gram.extract::<PythonExpression>()?.expr,
                CS.gram(
                    degree,
                    fundamental_representation.to_symbolic([]),
                    adjoint_representation.to_symbolic([]),
                )
            );

            assert_eq!(
                euclidean
                    .call_method0("casimir")?
                    .extract::<PythonExpression>()?
                    .expr,
                CS.cas(Atom::num(2), euclidean_representation.to_symbolic([]))
            );
            Ok(())
        })
        .unwrap();
    }
}
