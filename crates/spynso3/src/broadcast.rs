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

/// A unary function applied independently to each tensor component.
///
/// Calling it on a TensorExpression preserves the tensor interface. Calling it
/// on Tensor or TensorNetwork creates a lazy network. Scalar inputs produce
/// ordinary Symbolica Expressions. Register numerical behavior separately in
/// a TensorFunctionLibrary.
///
/// Examples
/// --------
/// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
/// >>> space = Representation.euc(2)
/// >>> A = TensorName("M")(space, space)
/// >>> from symbolica.community.tensor import Tensor
/// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
/// >>> from symbolica.community.tensor import BroadcastFunction
/// >>> conjugated = BroadcastFunction.conj()(tensor)
/// >>> conjugated.shape
/// (2, 2)
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
    /// Register a Symbolica function for elementwise tensor application.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Registered unary function name. It receives the broadcast tag and
    ///     cannot also be a tensor head.
    /// is_symmetric, is_antisymmetric, is_cyclesymmetric : bool, optional
    ///     Symmetry of the function arguments under all permutations, signed
    ///     permutations, or cyclic rotations. These affect scalar arguments as
    ///     well as index arguments. Repeated arguments in an antisymmetric call
    ///     make that call zero.
    /// is_linear : bool, optional
    ///     Distribute the function over sums in its arguments.
    /// is_flat : bool, optional
    ///     Flatten nested calls with the same head.
    /// is_scalar : bool, optional
    ///     Declare calls scalar for Symbolica's algebra. The broadcast acts on scalar component values.
    /// is_real, is_integer, is_positive : bool, optional
    ///     Assumptions used by Symbolica for the registered symbol.
    /// tags : sequence of str, optional
    ///     Additional Symbolica tags; Spenso's required tags are included automatically.
    /// aliases : sequence of str, optional
    ///     Additional names for the same Symbolica symbol.
    /// normalization : Transformer or callable, optional
    ///     Normalize a newly constructed call. A callable receives an Expression
    ///     and returns its normalized Expression; tensor interfaces must remain valid.
    /// print : callable, optional
    ///     Symbolica custom-print callback for the complete function call.
    ///     Returning None selects standard printing.
    /// derivative, series, eval : callable, optional
    ///     Symbolica callbacks for differentiation, series expansion, and numerical
    ///     evaluation. Their arguments and results follow ``symbolica.S``.
    /// data : object, optional
    ///     Symbolica user data attached to the symbol, such as a dict, list, or bytes.
    ///
    /// Returns
    /// -------
    /// BroadcastFunction
    ///     Callable wrapper with tensor-aware result types.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import BroadcastFunction, TensorExpression
    /// >>> function = BroadcastFunction("elementwise")
    /// >>> applied = function(TensorExpression(2))
    /// >>> applied.is_scalar
    /// True
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

    /// Return the built-in elementwise complex-conjugation function.
    ///
    /// Returns
    /// -------
    /// BroadcastFunction
    ///     Registered conjugation, with a built-in numerical implementation.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
    /// >>> space = Representation.euc(2)
    /// >>> A = TensorName("M")(space, space)
    /// >>> from symbolica.community.tensor import Tensor
    /// >>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
    /// >>> from symbolica.community.tensor import BroadcastFunction
    /// >>> BroadcastFunction.conj()(tensor).to_tensor()[0, 1]
    /// 2.0
    #[staticmethod]
    fn conj() -> Self {
        Self {
            name: INBUILTS.conj,
        }
    }

    #[doc = python_doc!("BroadcastFunction.__call__")]
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

    /// Return this function name as a Symbolica symbol.
    ///
    /// Returns
    /// -------
    /// Expression
    ///     The bare function head, without arguments or tensor structure.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import BroadcastFunction
    /// >>> name = BroadcastFunction("tagged_BroadcastFunction", tags=["example::example"])
    /// >>> head = name.to_expression()
    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(Atom::var(self.name))
    }

    /// Check whether the registered function carries a tag.
    ///
    /// Parameters
    /// ----------
    /// tag : str
    ///     Exact tag name. An unqualified name also matches its "python::" form.
    ///
    /// Returns
    /// -------
    /// bool
    ///     Whether the tag is present.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import BroadcastFunction
    /// >>> name = BroadcastFunction("tagged_BroadcastFunction", tags=["example::example"])
    /// >>> name.has_tag("example::example")
    /// True
    fn has_tag(&self, tag: &str) -> bool {
        self.name.has_tag(tag)
            || (!tag.contains("::")
                && self
                    .name
                    .get_tags()
                    .iter()
                    .any(|candidate| candidate.strip_prefix("python::") == Some(tag)))
    }

    /// List the tags attached to the registered function.
    ///
    /// Returns
    /// -------
    /// list of str
    ///     Fully qualified tags, including the tags required by Spenso.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica.community.tensor import BroadcastFunction
    /// >>> name = BroadcastFunction("tagged_BroadcastFunction", tags=["example::example"])
    /// >>> "example::example" in name.get_tags()
    /// True
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
                doc: python_doc!("BroadcastFunction.__call__"),
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
                doc: python_doc!("BroadcastFunction.__call__"),
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
                doc: python_doc!("BroadcastFunction.__call__"),
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
        ],
    }
}
