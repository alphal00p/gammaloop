use pyo3::{
    PyTypeInfo,
    prelude::*,
    types::{PyAny, PyTuple},
};

use spenso::network::{library::function_lib::INBUILTS, tags::SPENSO_TAG};
use symbolica::{
    api::python::{PythonExpression, PythonNormalization, PythonUserData},
    atom::{Atom, AtomView, DefaultNamespace, FunctionBuilder, Symbol},
};

use crate::{
    ModuleInit, Spensor,
    expression::{TensorExpression, TensorOperand},
    network::SpensoNet,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType,
    generate::MethodType,
    inventory::submit,
    type_info::{MethodInfo, ParameterDefault, ParameterInfo, ParameterKind, PyMethodsInfo},
};
#[cfg(feature = "python_stubgen")]
use symbolica::api::python::ConvertibleToExpression;

/// A unary Symbolica function whose action is broadcast over tensor entries.
///
/// Examples
/// --------
/// >>> from symbolica.community import tensor as sp
/// >>> square = sp.BroadcastFunction("docs::square")
/// >>> expression = square(sp.TensorExpression(2))
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
/// print : Any, optional
///     Optional Symbolica print callback. See Expression.symbol for the complete callback signature.
/// derivative : Any, optional
///     Optional Symbolica differentiation callback for the registered function; it follows Expression.symbol’s callback contract.
/// series : Any, optional
///     Optional Symbolica series-expansion callback for the registered function; it follows Expression.symbol’s callback contract.
/// eval : Any, optional
///     Optional Symbolica numerical-evaluation callback for the registered function. Concrete broadcasts additionally use TensorFunctionLibrary.
/// data : Expression or scalar | str | dict | list | bytes, optional
///     Optional Symbolica metadata attached to the symbol, such as an expression, text, or structured data.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "BroadcastFunction",
    module = "symbolica.community.tensor"
)]
#[derive(Clone)]
pub struct SpensoBroadcastFunction {
    pub(crate) name: Symbol,
}

impl ModuleInit for SpensoBroadcastFunction {}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl SpensoBroadcastFunction {
    /// Register a Symbolica function for elementwise application to tensor entries.
    ///
    /// The function is always tagged as a Spenso broadcast and cannot also be a tensor head.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> square = sp.BroadcastFunction("docs::square")
    /// >>> expression = square(sp.TensorExpression(2))
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
    /// print : Any, optional
    ///     Optional Symbolica print callback. See Expression.symbol for the complete callback signature.
    /// derivative : Any, optional
    ///     Optional Symbolica differentiation callback for the registered function; it follows Expression.symbol’s callback contract.
    /// series : Any, optional
    ///     Optional Symbolica series-expansion callback for the registered function; it follows Expression.symbol’s callback contract.
    /// eval : Any, optional
    ///     Optional Symbolica numerical-evaluation callback for the registered function. Concrete broadcasts additionally use TensorFunctionLibrary.
    /// data : Expression or scalar | str | dict | list | bytes, optional
    ///     Optional Symbolica metadata attached to the symbol, such as an expression, text, or structured data.
    #[new]
    #[pyo3(signature = (name, *, is_symmetric=None, is_antisymmetric=None, is_cyclesymmetric=None, is_linear=None, is_flat=None, is_scalar=None, is_real=None, is_integer=None, is_positive=None, tags=None, aliases=None, normalization=None, print=None, derivative=None, series=None, eval=None, data=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
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
        print: Option<Py<PyAny>>,
        derivative: Option<Py<PyAny>>,
        series: Option<Py<PyAny>>,
        eval: Option<Py<PyAny>>,
        data: Option<PythonUserData>,
    ) -> PyResult<Self> {
        let mut tags = tags.unwrap_or_default();
        tags.retain(|tag| tag != &SPENSO_TAG.tensor);
        if !tags.iter().any(|tag| tag == &SPENSO_TAG.broadcast) {
            tags.push(SPENSO_TAG.broadcast.clone());
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
        Ok(Self {
            name: name.get_symbol(),
        })
    }

    /// Return Spenso's registered complex-conjugation broadcast function.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> from symbolica import E, S
    /// >>> conjugate = sp.BroadcastFunction.conj()
    /// >>> result = conjugate(E("1 + 2i"))
    #[staticmethod]
    fn conj() -> Self {
        Self {
            name: INBUILTS.conj,
        }
    }

    /// Apply this unary function elementwise, preserving a structured interface.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> r = sp.Representation.euc(2)
    /// >>> conjugate = sp.BroadcastFunction.conj()
    /// >>> tensor = sp.Tensor.dense(sp.TensorName.vector("docs::p")(r), [1 + 2j, 3 + 4j])
    /// >>> result = conjugate(tensor).to_tensor()
    ///
    /// Parameters
    /// ----------
    /// arg : Tensor | TensorNetwork
    ///     Scalar, TensorExpression, Tensor, or TensorNetwork to transform elementwise. Component tensors are evaluated lazily as networks.
    #[gen_stub(skip)]
    fn __call__(&self, py: Python<'_>, arg: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(network) = arg.extract::<SpensoNet>() {
            return network.broadcast_function(self.name).and_then(|value| {
                value
                    .into_pyobject(py)
                    .map(|value| value.unbind().into_any())
            });
        }
        if let Ok(tensor) = arg.extract::<Spensor>() {
            return SpensoNet::from_tensor(tensor)?
                .broadcast_function(self.name)
                .and_then(|value| {
                    value
                        .into_pyobject(py)
                        .map(|value| value.unbind().into_any())
                });
        }
        match TensorOperand::extract(arg)? {
            TensorOperand::Structured(value) => {
                let (expression, interface) = value.into_parts();
                let atom = FunctionBuilder::new(self.name).add_arg(expression).finish();
                TensorExpression::from_known_parts(py, atom, interface, None, Vec::new())
                    .map(Py::into_any)
            }
            TensorOperand::Scalar(arg) => {
                let expression =
                    PythonExpression::from(FunctionBuilder::new(self.name).add_arg(arg).finish());
                Ok(expression.into_pyobject(py)?.unbind().into_any())
            }
        }
    }

    /// Return a readable object description for inspection.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> conjugate = sp.BroadcastFunction.conj()
    /// >>> text = repr(conjugate)
    fn __repr__(&self) -> String {
        format!("BroadcastFunction({:?})", self.name)
    }

    /// Return a readable text representation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> conjugate = sp.BroadcastFunction.conj()
    /// >>> text = str(conjugate)
    fn __str__(&self) -> String {
        self.name.to_string()
    }

    /// Return the registered function as a Symbolica variable expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> conjugate = sp.BroadcastFunction.conj()
    /// >>> symbolic_head = conjugate.to_expression()
    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(Atom::var(self.name))
    }

    /// Whether the registered function carries `tag`.
    ///
    /// Unqualified names also match tags registered in Symbolica's Python namespace.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> conjugate = sp.BroadcastFunction.conj()
    /// >>> tags = conjugate.get_tags()
    /// >>> assert all(conjugate.has_tag(tag) for tag in tags)
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

    /// Return the fully qualified tags attached to the registered function.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community import tensor as sp
    /// >>> conjugate = sp.BroadcastFunction.conj()
    /// >>> tags = conjugate.get_tags()
    fn get_tags(&self) -> Vec<String> {
        self.name.get_tags().to_vec()
    }
}

#[cfg(feature = "python_stubgen")]
submit! {
    PyMethodsInfo {
        struct_id: std::any::TypeId::of::<SpensoBroadcastFunction>,
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
                        name: "arg",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: || Spensor::type_input() | SpensoNet::type_input(),
                    },
                ],
                r#type: MethodType::Instance,
                r#return: SpensoNet::type_output,
                doc: r##"Apply this unary function lazily to a concrete tensor or tensor network.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> conjugate = sp.BroadcastFunction.conj()
>>> tensor = sp.Tensor.dense(sp.TensorName.vector("docs::p")(r), [1 + 2j, 3 + 4j])
>>> result = conjugate(tensor).to_tensor()

Parameters
----------
arg : Tensor | TensorNetwork
    Scalar, TensorExpression, Tensor, or TensorNetwork to transform elementwise. Component tensors are evaluated lazily as networks."##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__call__",
                parameters: &[
                    ParameterInfo {
                        name: "arg",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: TensorExpression::type_input,
                    },
                ],
                r#type: MethodType::Instance,
                r#return: TensorExpression::type_output,
                doc: r##"Apply this unary function lazily to a concrete tensor or tensor network.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> conjugate = sp.BroadcastFunction.conj()
>>> tensor = sp.Tensor.dense(sp.TensorName.vector("docs::p")(r), [1 + 2j, 3 + 4j])
>>> result = conjugate(tensor).to_tensor()

Parameters
----------
arg : Tensor | TensorNetwork
    Scalar, TensorExpression, Tensor, or TensorNetwork to transform elementwise. Component tensors are evaluated lazily as networks."##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
            MethodInfo {
                name: "__call__",
                parameters: &[
                    ParameterInfo {
                        name: "arg",
                        kind: ParameterKind::PositionalOrKeyword,
                        default: ParameterDefault::None,
                        type_info: ConvertibleToExpression::type_input,
                    },
                ],
                r#type: MethodType::Instance,
                r#return: PythonExpression::type_output,
                doc: r##"Apply this unary function lazily to a concrete tensor or tensor network.

Examples
--------
>>> from symbolica.community import tensor as sp
>>> r = sp.Representation.euc(2)
>>> conjugate = sp.BroadcastFunction.conj()
>>> tensor = sp.Tensor.dense(sp.TensorName.vector("docs::p")(r), [1 + 2j, 3 + 4j])
>>> result = conjugate(tensor).to_tensor()

Parameters
----------
arg : Tensor | TensorNetwork
    Scalar, TensorExpression, Tensor, or TensorNetwork to transform elementwise. Component tensors are evaluated lazily as networks."##,
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
        ],
    }
}
