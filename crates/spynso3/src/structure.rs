use idenso::representations::{ColorAdjoint, ColorFundamental, ColorSextet};

#[cfg(not(feature = "python_stubgen"))]
use pyo3_stub_gen_derive::remove_gen_stub;

use pyo3::{
    PyTypeInfo,
    exceptions::{self, PyTypeError, PyValueError},
    prelude::*,
    pybacked::PyBackedStr,
    types::{PyAny, PyDict, PyTuple},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};
use spenso::structure::slot::DualSlotTo;
use spenso::{
    network::{
        library::symbolic::ETS,
        parsing::ShadowedStructure,
        tags::{SPENSO_TAG, TENSOR_PRINT_CALLBACK_TAG, TENSOR_PRINT_HEAD_PREFIX},
    },
    structure::{
        Canonicalized, TensorStructure,
        abstract_index::AbstractIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{Euclidean, LibraryRep, Minkowski, RepName, Representation},
        slot::{IsAbstractSlot, Slot},
    },
};
use symbolica::{
    api::python::{PythonNormalization, PythonUserData, get_namespace},
    atom::{
        Atom, AtomView, DefaultNamespace, FunctionBuilder, NamespacedSymbol, Symbol, SymbolBuilder,
    },
    symbol,
};

use symbolica::api::python::{ConvertibleToExpression, PythonExpression};

use idenso::{
    color::CS, dirac::AGS, representations::Bispinor, tensor::inference::InterfaceInference,
};

use super::{ModuleInit, expression::TensorExpression};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{PyStubType, derive::*, impl_stub_type};

pub struct ConvertibleToSpensoName(pub SpensoName, pub Vec<Atom>);

impl<'a, 'py> FromPyObject<'a, 'py> for ConvertibleToSpensoName {
    type Error = PyErr;

    fn extract(structure: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        if let Ok(structure) = structure.extract::<SpensoName>() {
            Ok(ConvertibleToSpensoName(structure, Vec::new()))
        } else if let Ok(expression) = structure.extract::<PyRef<'_, TensorExpression>>() {
            if !expression.interface().canonical().is_scalar() {
                return Err(PyTypeError::new_err(
                    "a TensorExpression used as a name must have rank zero",
                ));
            }
            let AtomView::Fun(function) = expression.atom().as_view() else {
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
            if !s.interface().canonical().is_scalar() {
                return Err(PyTypeError::new_err(
                    "tensor key arguments must be scalar expressions",
                ));
            }
            Ok(SpensoSlotOrArgOrRep::Arg(PythonExpression {
                expr: s.atom().clone(),
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

/// Register or look up a symbolic tensor head. Calling it with Slots or
/// Representations constructs a TensorExpression. Scalar arguments describe the
/// tensor’s identity and are kept separate from its external axes.
///
/// With only a name, an existing symbol keeps its attributes, callbacks, and printers.
/// Use a fresh namespaced name when registering a different tensor convention.
///
/// The predefined accessors such as `TensorName.dirac_gamma()` and `TensorName.color_t()`
/// return raw fixed heads for introspection and pattern construction; those
/// reserved names cannot be called directly. Use the matching
/// `TensorExpression` factory for concrete tensors or `TensorPattern` shortcut
/// for rewrite patterns.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> r = sp.Representation.euc(2)
/// >>> A = sp.TensorName("docs::symmetric_A", is_symmetric=True)
/// >>> tensor = A(r("i"), r("j"))
///
/// Parameters
/// ----------
/// name : str
///     Registered function name. Use a namespace such as "my_model::A" to avoid collisions; reusing a name also reuses its Symbolica identity.
/// rank : int, optional
///     Optional rank declaration. Rank one enables vector-specific construction; omit it for a general tensor head.
/// is_symmetric : bool, optional
///     Declare symmetry under permutations of function arguments. None retains the existing or default Symbolica attribute.
/// is_antisymmetric : bool, optional
///     Declare antisymmetry under permutations of function arguments. None retains the existing or default Symbolica attribute.
/// is_cyclesymmetric : bool, optional
///     Declare invariance under cyclic argument permutations. None retains the existing or default Symbolica attribute.
/// is_linear : bool, optional
///     Declare a linear function so Symbolica can distribute it over supported sums. None retains the existing or default attribute.
/// is_flat : bool, optional
///     Flatten nested applications of the same function. None retains the existing or default Symbolica attribute.
/// is_scalar : bool, optional
///     Mark the registered function as scalar for Symbolica. This flag does not supply tensor component data.
/// is_real : bool, optional
///     Declare real-valued output for symbolic simplification. None retains the existing or default attribute.
/// is_integer : bool, optional
///     Declare integer-valued output for symbolic simplification. None retains the existing or default attribute.
/// is_positive : bool, optional
///     Declare positive output for symbolic simplification. None retains the existing or default attribute.
/// tags : Sequence[str], optional
///     Additional Symbolica tags. The required Spenso tensor or broadcast tag is attached automatically.
/// aliases : Sequence[str], optional
///     Alternative names registered for the same Symbolica symbol.
/// normalization : Transformer | Callable[[Expression], Expression], optional
///     Optional Transformer or callable that normalizes a Symbolica function application. The callable receives and returns an Expression.
/// print : dict[str, str] | Callable[..., str | None], optional
///     A mapping with "typst", "latex", and/or "plain" keys customizes the tensor name only, for example {"latex": r"\bar{J}"}. Values omit math delimiters. A Symbolica print callback instead controls the full display; returning None selects standard tensor notation.
/// derivative : Any, optional
///     Optional Symbolica differentiation callback for the registered function; it follows Expression.symbol’s callback contract.
/// series : Any, optional
///     Optional Symbolica series-expansion callback for the registered function; it follows Expression.symbol’s callback contract.
/// eval : Any, optional
///     Optional Symbolica numerical-evaluation callback for the registered function. Concrete broadcasts additionally use TensorFunctionLibrary.
/// data : Expression or scalar | str | dict | list | bytes, optional
///     Optional Symbolica metadata attached to the symbol, such as an expression, text, or structured data.
///
/// Returns
/// -------
/// TensorName
///     A new TensorName with the specified properties
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "TensorName",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub struct SpensoName {
    pub name: Symbol,
    // pub args: Vec<Atom>,
}

impl ModuleInit for SpensoName {}

impl SpensoName {
    fn prepare_print(
        py: Python<'_>,
        print: Option<Py<PyAny>>,
        tags: &mut Vec<String>,
    ) -> PyResult<Option<Py<PyAny>>> {
        let Some(print) = print else {
            return Ok(None);
        };
        if let Ok(mapping) = print.bind(py).cast::<PyDict>() {
            let sources = mapping.extract::<std::collections::BTreeMap<String, String>>()?;
            for (backend, source) in &sources {
                if !matches!(backend.as_str(), "plain" | "latex" | "typst") {
                    return Err(PyValueError::new_err(
                        "print mapping keys must be 'plain', 'latex', or 'typst'",
                    ));
                }
                tags.push(format!("{TENSOR_PRINT_HEAD_PREFIX}{backend}:{source}"));
            }
            Ok(None)
        } else if print.bind(py).is_callable() {
            tags.push(TENSOR_PRINT_CALLBACK_TAG.to_owned());
            Ok(Some(print))
        } else {
            Err(PyTypeError::new_err(
                "print must be a mapping of backend names to strings, a callable, or None",
            ))
        }
    }

    /// Whether this name declares a rank-one tensor.
    pub fn is_vector(&self) -> bool {
        self.name.has_tag(&SPENSO_TAG.rank1)
    }

    pub(crate) fn builtin_factory(&self) -> Option<&'static str> {
        if self.name == ETS.metric {
            Some("g")
        } else if self.name == ETS.flat {
            Some("flat")
        } else if self.name == AGS.gamma {
            Some("dirac_gamma")
        } else if self.name == AGS.gamma5 {
            Some("gamma5")
        } else if self.name == AGS.projm {
            Some("projm")
        } else if self.name == AGS.projp {
            Some("projp")
        } else if self.name == AGS.sigma {
            Some("sigma")
        } else if self.name == CS.f {
            Some("color_f")
        } else if self.name == CS.t {
            Some("color_t")
        } else {
            None
        }
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoName {
    fn _repr_html_(&self, py: Python<'_>) -> Option<String> {
        crate::display::atom_to_html(
            py,
            &self.to_expression().expr,
            &crate::display::DisplaySettings::default(),
            None,
        )
        .ok()
    }

    fn _repr_latex_(&self) -> String {
        crate::display::atom_to_latex(&self.to_expression().expr, false, None)
    }

    /// Register or look up a symbolic tensor head. Calling it with Slots or
    /// Representations constructs a TensorExpression. Scalar arguments describe the
    /// tensor’s identity and are kept separate from its external axes.
    ///
    /// With only a name, an existing symbol keeps its attributes, callbacks, and printers.
    /// Use a fresh namespaced name when registering a different tensor convention.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> A = sp.TensorName("docs::symmetric_A", is_symmetric=True)
    /// >>> tensor = A(r("i"), r("j"))
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Registered function name. Use a namespace such as "my_model::A" to avoid collisions; reusing a name also reuses its Symbolica identity.
    /// rank : int, optional
    ///     Optional rank declaration. Rank one enables vector-specific construction; omit it for a general tensor head.
    /// is_symmetric : bool, optional
    ///     Declare symmetry under permutations of function arguments. None retains the existing or default Symbolica attribute.
    /// is_antisymmetric : bool, optional
    ///     Declare antisymmetry under permutations of function arguments. None retains the existing or default Symbolica attribute.
    /// is_cyclesymmetric : bool, optional
    ///     Declare invariance under cyclic argument permutations. None retains the existing or default Symbolica attribute.
    /// is_linear : bool, optional
    ///     Declare a linear function so Symbolica can distribute it over supported sums. None retains the existing or default attribute.
    /// is_flat : bool, optional
    ///     Flatten nested applications of the same function. None retains the existing or default Symbolica attribute.
    /// is_scalar : bool, optional
    ///     Mark the registered function as scalar for Symbolica. This flag does not supply tensor component data.
    /// is_real : bool, optional
    ///     Declare real-valued output for symbolic simplification. None retains the existing or default attribute.
    /// is_integer : bool, optional
    ///     Declare integer-valued output for symbolic simplification. None retains the existing or default attribute.
    /// is_positive : bool, optional
    ///     Declare positive output for symbolic simplification. None retains the existing or default attribute.
    /// tags : Sequence[str], optional
    ///     Additional Symbolica tags. The required Spenso tensor or broadcast tag is attached automatically.
    /// aliases : Sequence[str], optional
    ///     Alternative names registered for the same Symbolica symbol.
    /// normalization : Transformer | Callable[[Expression], Expression], optional
    ///     Optional Transformer or callable that normalizes a Symbolica function application. The callable receives and returns an Expression.
    /// print : dict[str, str] | Callable[..., str | None], optional
    ///     A mapping with "typst", "latex", and/or "plain" keys customizes the tensor name only, for example {"latex": r"\bar{J}"}. Values omit math delimiters. A Symbolica print callback instead controls the full display; returning None selects standard tensor notation.
    /// derivative : Any, optional
    ///     Optional Symbolica differentiation callback for the registered function; it follows Expression.symbol’s callback contract.
    /// series : Any, optional
    ///     Optional Symbolica series-expansion callback for the registered function; it follows Expression.symbol’s callback contract.
    /// eval : Any, optional
    ///     Optional Symbolica numerical-evaluation callback for the registered function. Concrete broadcasts additionally use TensorFunctionLibrary.
    /// data : Expression or scalar | str | dict | list | bytes, optional
    ///     Optional Symbolica metadata attached to the symbol, such as an expression, text, or structured data.
    ///
    /// Returns
    /// -------
    /// TensorName
    ///     A new TensorName with the specified properties
    #[new]
    #[pyo3(signature = (name, *, rank=None, is_symmetric=None, is_antisymmetric=None, is_cyclesymmetric=None, is_linear=None, is_flat=None, is_scalar=None, is_real=None, is_integer=None, is_positive=None, tags=None, aliases=None, normalization=None, print=None, derivative=None, series=None, eval=None, data=None))]
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
        #[gen_stub(override_type(
            type_repr = "typing.Optional[symbolica.core.Transformer | typing.Callable[[symbolica.core.Expression], symbolica.core.Expression]]",
            imports = ("typing", "symbolica.core")
        ))]
        normalization: Option<PythonNormalization>,
        #[gen_stub(override_type(
            type_repr = "typing.Optional[dict[str, str] | typing.Callable[..., str | None]]",
            imports = ("typing",)
        ))]
        print: Option<Py<PyAny>>,
        derivative: Option<Py<PyAny>>,
        series: Option<Py<PyAny>>,
        eval: Option<Py<PyAny>>,
        data: Option<PythonUserData>,
    ) -> PyResult<Self> {
        if rank.is_none()
            && [
                is_symmetric,
                is_antisymmetric,
                is_cyclesymmetric,
                is_linear,
                is_flat,
                is_scalar,
                is_real,
                is_integer,
                is_positive,
            ]
            .iter()
            .all(Option::is_none)
            && tags.is_none()
            && aliases.is_none()
            && normalization.is_none()
            && print.is_none()
            && derivative.is_none()
            && series.is_none()
            && eval.is_none()
            && data.is_none()
        {
            let namespace = DefaultNamespace {
                namespace: get_namespace(py)?.into(),
                data: "",
                file: "".into(),
                line: 0,
            };
            if let Some(name) = Symbol::get_symbol(namespace.attach_namespace(&name))
                && name.has_tag(&SPENSO_TAG.tensor)
            {
                return Ok(Self { name });
            }
        }
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

        let print = Self::prepare_print(py, print, &mut tags)?;
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

    /// Create a rank-one tensor name; ``print`` accepts the same mappings and callbacks as ``TensorName``.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> p = sp.TensorName.vector("docs::p")
    /// >>> vector = p(r("mu"))
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Registered function name. Use a namespace such as "my_model::A" to avoid collisions; reusing a name also reuses its Symbolica identity.
    /// is_symmetric : bool, optional
    ///     Declare symmetry under permutations of function arguments. None retains the existing or default Symbolica attribute.
    /// is_antisymmetric : bool, optional
    ///     Declare antisymmetry under permutations of function arguments. None retains the existing or default Symbolica attribute.
    /// is_cyclesymmetric : bool, optional
    ///     Declare invariance under cyclic argument permutations. None retains the existing or default Symbolica attribute.
    /// is_linear : bool, optional
    ///     Declare a linear function so Symbolica can distribute it over supported sums. None retains the existing or default attribute.
    /// is_flat : bool, optional
    ///     Flatten nested applications of the same function. None retains the existing or default Symbolica attribute.
    /// is_scalar : bool, optional
    ///     Mark the registered function as scalar for Symbolica. This flag does not supply tensor component data.
    /// is_real : bool, optional
    ///     Declare real-valued output for symbolic simplification. None retains the existing or default attribute.
    /// is_integer : bool, optional
    ///     Declare integer-valued output for symbolic simplification. None retains the existing or default attribute.
    /// is_positive : bool, optional
    ///     Declare positive output for symbolic simplification. None retains the existing or default attribute.
    /// tags : Sequence[str], optional
    ///     Additional Symbolica tags. The required Spenso tensor or broadcast tag is attached automatically.
    /// aliases : Sequence[str], optional
    ///     Alternative names registered for the same Symbolica symbol.
    /// normalization : Transformer | Callable[[Expression], Expression], optional
    ///     Optional Transformer or callable that normalizes a Symbolica function application. The callable receives and returns an Expression.
    /// print : dict[str, str] | Callable[..., str | None], optional
    ///     A mapping with "typst", "latex", and/or "plain" keys customizes the tensor name only, for example {"latex": r"\bar{J}"}. Values omit math delimiters. A Symbolica print callback instead controls the full display; returning None selects standard tensor notation.
    /// derivative : Any, optional
    ///     Optional Symbolica differentiation callback for the registered function; it follows Expression.symbol’s callback contract.
    /// series : Any, optional
    ///     Optional Symbolica series-expansion callback for the registered function; it follows Expression.symbol’s callback contract.
    /// eval : Any, optional
    ///     Optional Symbolica numerical-evaluation callback for the registered function. Concrete broadcasts additionally use TensorFunctionLibrary.
    /// data : Expression or scalar | str | dict | list | bytes, optional
    ///     Optional Symbolica metadata attached to the symbol, such as an expression, text, or structured data.
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
        #[gen_stub(override_type(
            type_repr = "typing.Optional[symbolica.core.Transformer | typing.Callable[[symbolica.core.Expression], symbolica.core.Expression]]",
            imports = ("typing", "symbolica.core")
        ))]
        normalization: Option<PythonNormalization>,
        #[gen_stub(override_type(
            type_repr = "typing.Optional[dict[str, str] | typing.Callable[..., str | None]]",
            imports = ("typing",)
        ))]
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

    /// Build a symbolic tensor from scalar arguments and ordered ports.
    /// A Slot fixes an index label; a Representation leaves a port unresolved. A rank-one
    /// TensorExpression can bind a port. Repeated compatible explicit indices contract.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> name = sp.TensorName("docs::A")
    /// >>> tensor = name(r("i"), r)
    /// >>> indexed = tensor("j")
    ///
    /// Parameters
    /// ----------
    /// args : Slot | Representation | Expression or scalar
    ///     Scalar arguments followed by structural ports, in logical order. Slots set explicit labels, Representations leave ports unresolved, and a rank-one tensor can bind a port.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     A structured expression with explicit and/or unresolved ports
    #[pyo3(signature = (*args))]
    #[gen_stub(skip)]
    fn __call__(
        &self,
        py: Python<'_>,
        args: &Bound<'_, PyTuple>,
    ) -> PyResult<Py<TensorExpression>> {
        if let Some(factory) = self.builtin_factory() {
            return Err(PyTypeError::new_err(format!(
                "predefined tensor `{factory}` has a fixed structure; use \
                 TensorExpression.{factory}(...) for concrete tensors or \
                 TensorPattern.{factory}(...) for patterns"
            )));
        }

        let mut scalar_args = Vec::new();
        let mut ports = Vec::new();
        let mut port_atoms = Vec::new();
        let rank_one = self.name.has_tag(&SPENSO_TAG.rank1);
        let mut structural_seen = false;
        let mut next_open = 0;
        let mut bound_ports = false;
        let mut inference = InterfaceInference::default();

        for arg_bound in args.iter() {
            let convertible = arg_bound.extract::<SpensoSlotOrArgOrRep>()?;

            match convertible {
                SpensoSlotOrArgOrRep::Arg(expr) => {
                    if !rank_one && inference.compact_vector_port(expr.expr.as_view()) {
                        structural_seen = true;
                        bound_ports = true;
                        port_atoms.push(expr.expr);
                        continue;
                    }
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
        let atom = FunctionBuilder::new(self.name)
            .add_args(&scalar_args)
            .add_args(&port_atoms)
            .finish();
        let interface = PartialStructure::from_logical_slots(ports);
        if bound_ports
            || !scalar_args.is_empty()
            || self.name.get_normalization_function().is_some()
            || self.name.get_evaluation_info().is_some()
        {
            idenso::tensor::SymbolicTensor::validate_interface(&atom, &interface)
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
        }
        let (name, arguments) = if bound_ports {
            // A bound vertex is a contraction, not the stored tensor with its
            // surviving signature. Keep the bound operands in the expression.
            (None, Vec::new())
        } else {
            (Some(self.name), scalar_args)
        };
        TensorExpression::from_known_parts(py, atom, interface, name, arguments)
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.TensorName("docs::A")
    /// >>> text = repr(name)
    fn __repr__(&self) -> String {
        format!("TensorName({:?})", self.name.get_name())
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.TensorName("docs::A")
    /// >>> text = str(name)
    fn __str__(&self) -> String {
        format!("{}", self.name)
    }

    /// Convert the tensor name to a symbolic expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.TensorName("docs::A")
    /// >>> symbolic_head = name.to_expression()
    ///
    /// Returns
    /// -------
    /// Expression
    ///     A symbolic Expression representing this tensor name
    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(Atom::var(self.name))
    }

    /// Check whether this tensor name carries `tag`.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.TensorName("docs::A")
    /// >>> tags = name.get_tags()
    /// >>> assert all(name.has_tag(tag) for tag in tags)
    ///
    /// Parameters
    /// ----------
    /// tag : str
    ///     Tag name to test. Fully qualified names are unambiguous; unqualified names also match tags in the Python namespace.
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
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> name = sp.TensorName("docs::A")
    /// >>> tags = name.get_tags()
    fn get_tags(&self) -> Vec<String> {
        self.name.get_tags().to_vec()
    }

    /// Predefined metric tensor name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.g()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn g() -> SpensoName {
        SpensoName { name: ETS.metric }
    }

    /// Predefined musical isomorphism tensor name. This enables dualizing self dual indices.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.flat()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn flat() -> SpensoName {
        SpensoName { name: ETS.flat }
    }

    /// Predefined gamma matrix name for introspection and pattern construction.
    ///
    /// The matching typed factories use storage order: bispinor-in,
    /// bispinor-out, then Minkowski.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.gamma()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn dirac_gamma() -> SpensoName {
        SpensoName { name: AGS.gamma }
    }

    /// Predefined gamma5 matrix name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.gamma5()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn gamma5() -> SpensoName {
        SpensoName { name: AGS.gamma5 }
    }

    /// Predefined left chiral projector name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.projm()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn projm() -> SpensoName {
        SpensoName { name: AGS.projm }
    }

    /// Predefined right chiral projector name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.projp()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn projp() -> SpensoName {
        SpensoName { name: AGS.projp }
    }

    /// Predefined sigma matrix name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.sigma()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn sigma() -> SpensoName {
        SpensoName { name: AGS.sigma }
    }

    /// Predefined color structure constant name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.f()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn color_f() -> SpensoName {
        SpensoName { name: CS.f }
    }

    /// Predefined color generator name.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> head = sp.TensorName.t()
    /// >>> tags = head.get_tags()
    #[staticmethod]
    fn color_t() -> SpensoName {
        SpensoName { name: CS.t }
    }
}

/// Internal indexed structure used while constructing a `TensorExpression`.
#[derive(Clone)]
pub(crate) struct SpensoIndices {
    pub(crate) structure: Canonicalized<ShadowedStructure<AbstractIndex>>,
}

pub enum ArithmeticStructure {
    Tensor(Py<TensorExpression>),
    Convertible(ConvertibleToExpression),
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for ArithmeticStructure {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        ConvertibleToExpression::type_input() | TensorExpression::type_input()
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
    module = "symbolica.community.tensor"
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
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> space = sp.Representation("docs::custom_space", 3)
/// >>> slot = space("i")
///
/// Parameters
/// ----------
/// name : str
///     Name of the custom representation. Use a namespace to distinguish it from built-in representations.
/// dimension : int | Expression | str
///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
/// is_self_dual : bool
///     Whether the representation pairs with itself. If False, .dual() identifies the distinct dual representation.
pub struct SpensoRepresentation {
    pub representation: Representation<LibraryRep>,
}

impl ModuleInit for SpensoRepresentation {}

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

pub struct ConvertibleToDimension(pub(crate) Dimension);

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
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoRepresentation {
    fn _repr_html_(&self) -> String {
        self.to_html()
    }

    /// Compact name, dimension, duality, and canonical metric metadata.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> html = r.to_html()
    fn to_html(&self) -> String {
        crate::display::metadata::representation(self.representation)
    }

    /// Dimension-independent representation identity, including duality.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> result = r.name
    #[getter]
    fn name(&self) -> crate::metadata::SpensoRepresentationName {
        crate::metadata::SpensoRepresentationName {
            rep: self.representation.rep,
        }
    }

    fn _repr_latex_(&self) -> String {
        crate::display::atom_to_latex(&self.to_expression().expr, false, None)
    }

    /// Create and register a new representation with specified properties.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> space = sp.Representation("docs::custom_space", 3)
    /// >>> slot = space("i")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Name of the custom representation. Use a namespace to distinguish it from built-in representations.
    /// dimension : int | Expression | str
    ///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
    /// is_self_dual : bool
    ///     Whether the representation pairs with itself. If False, .dual() identifies the distinct dual representation.
    #[new]
    #[pyo3(signature =(name,dimension,is_self_dual=true))]
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

    /// Return the representation paired with this one under index contraction.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> fundamental = sp.Representation.cof(3)
    /// >>> antifundamental = fundamental.dual()
    fn dual(&self) -> Self {
        Self {
            representation: self.representation.dual(),
        }
    }

    /// Return the dimension carried by this representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> result = r.dimension
    #[getter]
    fn dimension(&self) -> PythonExpression {
        self.representation.dim.to_symbolic().into()
    }

    /// Build the degree-k Casimir eigenvalue for this representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> fundamental = sp.Representation.cof(3)
    /// >>> invariant = fundamental.casimir()
    ///
    /// Parameters
    /// ----------
    /// degree : Expression or scalar
    ///     Degree of the representation invariant. For Casimir and Dynkin indices, the default 2 selects the quadratic invariant.
    #[pyo3(
        signature = (degree = ConvertibleToInvariantDegree::quadratic()),
        text_signature = "($self, degree=2)"
    )]
    fn casimir(&self, degree: ConvertibleToInvariantDegree) -> PythonExpression {
        CS.cas(degree.0, self.representation.to_symbolic([])).into()
    }

    /// Build the degree-k Dynkin index for this representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> fundamental = sp.Representation.cof(3)
    /// >>> invariant = fundamental.dynkin_index()
    ///
    /// Parameters
    /// ----------
    /// degree : Expression or scalar
    ///     Degree of the representation invariant. For Casimir and Dynkin indices, the default 2 selects the quadratic invariant.
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
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> fundamental = sp.Representation.cof(3)
    /// >>> invariant = fundamental.gram(2)
    ///
    /// Parameters
    /// ----------
    /// degree : Expression or scalar
    ///     Degree of the representation invariant. For Casimir and Dynkin indices, the default 2 selects the quadratic invariant.
    /// other : Representation, optional
    ///     Second representation in the Gram invariant. None uses this representation.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> mu = r("mu")
    /// >>> nu = r(1)
    ///
    /// Parameters
    /// ----------
    /// aind : int | str
    ///     Index label as an integer, string, or Symbolica expression. Equal compatible labels represent contracted indices.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> tensor = r.g("mu", "nu")
    ///
    /// Parameters
    /// ----------
    /// i : int | Expression | str
    ///     First index label, given as an integer, string, or Symbolica expression.
    /// j : int | Expression | str
    ///     Second index label; a repeated compatible label requests contraction.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> tensor = r.flat("mu", "nu")
    ///
    /// Parameters
    /// ----------
    /// i : int | Expression | str
    ///     First index label, given as an integer, string, or Symbolica expression.
    /// j : int | Expression | str
    ///     Second index label; a repeated compatible label requests contraction.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> tensor = r.id("mu", "nu")
    ///
    /// Parameters
    /// ----------
    /// i : int | Expression | str
    ///     First index label, given as an integer, string, or Symbolica expression.
    /// j : int | Expression | str
    ///     Second index label; a repeated compatible label requests contraction.
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

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> text = repr(r)
    fn __repr__(&self) -> String {
        format!("Representation({})", self.representation.to_symbolic([]))
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> text = str(r)
    fn __str__(&self) -> String {
        format!("{}", self.representation.to_symbolic([]))
    }

    /// Convert the representation to a symbolic expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.mink(4)
    /// >>> symbolic_space = r.to_expression()
    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(self.representation.to_symbolic([]))
    }

    /// Create a bispinor representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> space = sp.Representation.bis(4)
    /// >>> slot = space("i")
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | str
    ///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> space = sp.Representation.euc(3)
    /// >>> slot = space("i")
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | str
    ///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> space = sp.Representation.mink(4)
    /// >>> slot = space("i")
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | str
    ///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> space = sp.Representation.cof(3)
    /// >>> slot = space("i")
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | str
    ///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> space = sp.Representation.coad(8)
    /// >>> slot = space("i")
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | str
    ///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> space = sp.Representation.cos(6)
    /// >>> slot = space("i")
    ///
    /// Parameters
    /// ----------
    /// dimension : int | Expression | str
    ///     Positive concrete dimension or an atomic symbolic dimension. Built-in constructors reuse the registered representation identity.
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
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> slot = sp.Slot("docs::slot_space", 3, "i")
///
/// Parameters
/// ----------
/// name : str
///     Name of the custom representation. Use a namespace to distinguish it from built-in representations.
/// dimension : int
///     Concrete dimension of the custom representation.
/// aind : int | Expression | str
///     Index label as an integer, string, or Symbolica expression. Equal compatible labels represent contracted indices.
/// dual : bool
///     Use the dual orientation of the named representation when True.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    eq,
    name = "Slot",
    module = "symbolica.community.tensor"
)]
#[derive(Clone, PartialEq, Eq)]
pub struct SpensoSlot {
    pub slot: Slot<LibraryRep>,
}

impl ModuleInit for SpensoSlot {}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoSlot {
    /// Hash the immutable representation and index label for dictionaries and sets.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> labels = {r("i"): "free index"}
    /// >>> assert labels[r("i")] == "free index"
    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.slot.hash(&mut hasher);
        hasher.finish()
    }

    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.to_html(py, None)
    }

    /// Display the index label beside its representation and exact stored index.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Representation.mink(4)("mu")
    /// >>> html = slot.to_html()
    ///
    /// Parameters
    /// ----------
    /// settings : DisplaySettings, optional
    ///     Optional DisplaySettings controlling notation, indices, and component layout. None uses the default presentation.
    #[pyo3(signature = (*, settings=None))]
    fn to_html(
        &self,
        py: Python<'_>,
        settings: Option<crate::display::DisplaySettings>,
    ) -> PyResult<String> {
        crate::display::metadata::slot(py, self.slot, &settings.unwrap_or_default())
    }

    /// Exact symbolic index, independently of its presentation label.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Representation.mink(4)("mu")
    /// >>> result = slot.index
    #[getter]
    fn index(&self) -> PythonExpression {
        Atom::from(self.slot.aind).into()
    }

    fn _repr_latex_(&self) -> String {
        crate::display::atom_to_latex(&self.to_expression().expr, false, None)
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Representation.mink(4)("mu")
    /// >>> text = repr(slot)
    fn __repr__(&self) -> String {
        format!("Slot({})", self.slot.to_atom())
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Representation.mink(4)("mu")
    /// >>> text = str(slot)
    fn __str__(&self) -> String {
        format!("{}", self.slot.to_atom())
    }

    /// The representation carried by this slot, including its dimension and duality.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Representation.mink(4)("mu")
    /// >>> result = slot.representation
    #[getter]
    fn representation(&self) -> SpensoRepresentation {
        SpensoRepresentation {
            representation: self.slot.rep(),
        }
    }

    /// Return this slot with its representation replaced by the dual representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Representation.cof(3)("i")
    /// >>> dual_slot = slot.dual()
    fn dual(&self) -> Self {
        SpensoSlot {
            slot: self.slot.dual(),
        }
    }

    /// Create a new slot with a custom representation and index.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Slot("docs::slot_space", 3, "i")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Name of the custom representation. Use a namespace to distinguish it from built-in representations.
    /// dimension : int
    ///     Concrete dimension of the custom representation.
    /// aind : int | Expression | str
    ///     Index label as an integer, string, or Symbolica expression. Equal compatible labels represent contracted indices.
    /// dual : bool
    ///     Use the dual orientation of the named representation when True.
    #[new]
    #[pyo3(signature =(name,dimension,aind,dual=false))]
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
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> slot = sp.Representation.mink(4)("mu")
    /// >>> symbolic_slot = slot.to_expression()
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
                doc:r##"Build a symbolic tensor from scalar arguments and ordered ports.
A Slot fixes an index label; a Representation leaves a port unresolved. A rank-one
TensorExpression can bind a port. Repeated compatible explicit indices contract.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> name = sp.TensorName("docs::A")
>>> tensor = name(r("i"), r)
>>> indexed = tensor("j")

Parameters
----------
args : Slot | Representation | Expression or scalar
    Scalar arguments followed by structural ports, in logical order. Slots set explicit labels, Representations leave ports unresolved, and a rank-one tensor can bind a port.

Returns
-------
TensorExpression
    A structured expression with explicit and/or unresolved ports"##,
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
    fn tensor_name_bound_vectors_are_composites_and_metadata_stays_scalar() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let representation = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
            let tensor = SpensoName {
                name: spenso::tensor_symbol!("python_bound_tensor"),
            };
            let slot = || {
                Py::new(
                    py,
                    SpensoSlot {
                        slot: representation.slot(AbstractIndex::Normal(9)),
                    },
                )
                .map(Py::into_any)
            };
            let p = FunctionBuilder::new(spenso::vector_symbol!("python_bound_p"))
                .add_arg(representation.to_symbolic([]))
                .finish();
            let q = FunctionBuilder::new(spenso::vector_symbol!("python_bound_q"))
                .add_arg(representation.to_symbolic([]))
                .finish();
            let mut expressions = Vec::new();
            for vector in [&p, &q] {
                for bound_first in [true, false] {
                    let bound = PythonExpression::from(vector.clone())
                        .into_pyobject(py)?
                        .unbind()
                        .into_any();
                    let args = if bound_first {
                        vec![bound, slot()?]
                    } else {
                        vec![slot()?, bound]
                    };
                    let value = tensor.__call__(py, &PyTuple::new(py, args)?)?;
                    let value = value.bind(py).borrow();
                    assert_eq!(value.interface().canonical().order(), 1);
                    assert_eq!(
                        value.interface().logical_slots()[0].aind,
                        PartialIndex::Explicit(AbstractIndex::Normal(9))
                    );
                    assert!(value.name.is_none());
                    assert!(value.name_args.is_empty());
                    expressions.push(value.atom().clone());
                }
            }
            assert_ne!(expressions[0], expressions[2]);
            let nested_tensor = FunctionBuilder::new(tensor.name)
                .add_arg(representation.to_symbolic([]))
                .add_arg(representation.to_symbolic([]))
                .finish();
            let invalid = PyTuple::new(
                py,
                [
                    PythonExpression::from(nested_tensor.clone())
                        .into_pyobject(py)?
                        .unbind()
                        .into_any(),
                    slot()?,
                ],
            )?;
            assert!(tensor.__call__(py, &invalid).is_err());
            let malformed =
                FunctionBuilder::new(spenso::vector_symbol!("python_bound_missing_port"))
                    .add_arg(Atom::num(3))
                    .finish();
            let malformed = PyTuple::new(
                py,
                [
                    PythonExpression::from(malformed)
                        .into_pyobject(py)?
                        .unbind()
                        .into_any(),
                    slot()?,
                ],
            )?;
            assert!(tensor.__call__(py, &malformed).is_err());
            let metadata = FunctionBuilder::new(symbol!("python_bound_scalar"; Scalar))
                .add_arg(nested_tensor)
                .finish();
            let unbound = PyTuple::new(
                py,
                [
                    PythonExpression::from(metadata.clone())
                        .into_pyobject(py)?
                        .unbind()
                        .into_any(),
                    slot()?,
                ],
            )?;
            let value = tensor.__call__(py, &unbound)?;
            let value = value.bind(py).borrow();
            assert_eq!(value.name, Some(tensor.name));
            assert_eq!(value.name_args, vec![metadata]);
            Ok(())
        })
        .unwrap();
    }

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
            assert_eq!(expression_ref.interface().canonical().order(), 2);
            assert!(matches!(
                expression_ref.interface().logical_slots()[0].aind,
                PartialIndex::Explicit(AbstractIndex::Normal(41))
            ));
            assert!(matches!(
                expression_ref.interface().logical_slots()[1].aind,
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
            assert!(scalar.bind(py).borrow().interface().canonical().is_scalar());
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

            for (name, factory) in [
                (ETS.metric, "g"),
                (ETS.flat, "flat"),
                (AGS.gamma, "dirac_gamma"),
                (AGS.gamma5, "gamma5"),
                (AGS.projm, "projm"),
                (AGS.projp, "projp"),
                (AGS.sigma, "sigma"),
                (CS.f, "color_f"),
                (CS.t, "color_t"),
            ] {
                let error = SpensoName { name }
                    .__call__(py, &PyTuple::empty(py))
                    .expect_err("predefined tensors require their typed factory");
                assert!(error.is_instance_of::<PyTypeError>(py));
                assert!(
                    error
                        .to_string()
                        .contains(&format!("TensorExpression.{factory}"))
                );
                assert!(
                    error
                        .to_string()
                        .contains(&format!("TensorPattern.{factory}"))
                );
            }
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
