use pyo3::{
    PyTypeInfo,
    prelude::*,
    types::{PyAny, PyTuple},
};

use spenso::network::{library::function_lib::INBUILTS, tags::SPENSO_TAG};
use symbolica::{
    api::python::{PythonExpression, PythonTransformer, PythonUserData},
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
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    from_py_object,
    name = "BroadcastFunction",
    module = "symbolica.community.spenso"
)]
#[derive(Clone)]
pub struct SpensoBroadcastFunction {
    pub(crate) name: Symbol,
}

impl ModuleInit for SpensoBroadcastFunction {}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl SpensoBroadcastFunction {
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
        normalization: Option<PythonTransformer>,
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
    #[staticmethod]
    fn conj() -> Self {
        Self {
            name: INBUILTS.conj,
        }
    }

    /// Apply this unary function elementwise, preserving a structured interface.
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
                let atom = FunctionBuilder::new(self.name).add_arg(value.atom).finish();
                TensorExpression::from_atom_interface(py, atom, value.interface).map(Py::into_any)
            }
            TensorOperand::Scalar(arg) => {
                let expression =
                    PythonExpression::from(FunctionBuilder::new(self.name).add_arg(arg).finish());
                Ok(expression.into_pyobject(py)?.unbind().into_any())
            }
        }
    }

    fn __repr__(&self) -> String {
        format!("BroadcastFunction({:?})", self.name)
    }

    fn __str__(&self) -> String {
        self.name.to_string()
    }

    fn to_expression(&self) -> PythonExpression {
        PythonExpression::from(Atom::var(self.name))
    }

    fn has_tag(&self, tag: &str) -> bool {
        self.name.has_tag(tag)
            || (!tag.contains("::")
                && self
                    .name
                    .get_tags()
                    .iter()
                    .any(|candidate| candidate.strip_prefix("python::") == Some(tag)))
    }

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
                doc: "Apply this unary function lazily to a concrete tensor or tensor network.",
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
                doc: "Apply this unary function elementwise, preserving a structured interface.",
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
                doc: "Apply this unary function to an ordinary scalar expression.",
                is_async: false,
                deprecated: None,
                type_ignored: None,
                is_overload: true,
            },
        ],
    }
}
