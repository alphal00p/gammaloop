//! Shared typing refinements for the generated documentation and installed package.
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    generate::{Module, ParameterDefault, Parameters, VariableDef},
};
use symbolica::api::python::{ConvertibleToExpression, ConvertibleToReplaceWith};

fn parameters(
    parameters: &mut Parameters,
    mut refine: impl FnMut(&mut pyo3_stub_gen::generate::Parameter),
) {
    for parameter in parameters
        .positional_only
        .iter_mut()
        .chain(&mut parameters.positional_or_keyword)
        .chain(&mut parameters.keyword_only)
        .chain(&mut parameters.varargs)
        .chain(&mut parameters.varkw)
    {
        refine(parameter);
    }
}

fn named(name: &str) -> TypeInfo {
    TypeInfo {
        name: name.to_owned(),
        import: ["typing".into()].into(),
    }
}

const SYMBOLIC_FACTORS: &str = "TensorExpression | Expression | FactorProjector[TensorExpression]";
const CONCRETE_FACTORS: &str = "Tensor | TensorNetwork | FactorProjector[TensorNetwork]";
const ALL_FACTORS: &str = "TensorExpression | _ScalarInput | Tensor | TensorNetwork | FactorProjector[TensorExpression] | FactorProjector[TensorNetwork]";

pub(crate) fn refine(module: &mut Module) {
    module.doc = python_doc!("module").to_owned();
    let replacements = [
        (
            ConvertibleToReplaceWith::type_input().name,
            "_ReplacementInput",
        ),
        (ConvertibleToExpression::type_input().name, "_ScalarInput"),
        (
            "Float | builtins.int | builtins.float | builtins.str | decimal.Decimal".to_owned(),
            "_RealInput",
        ),
    ];
    let simplify = |info: &mut TypeInfo| {
        for (original, alias) in &replacements {
            info.name = info.name.replace(original, alias);
        }
        if info.name.ends_with(" | FactorProjector") {
            info.name
                .truncate(info.name.len() - " | FactorProjector".len());
            info.name
                .push_str(" | FactorProjector[TensorExpression] | FactorProjector[TensorNetwork]");
        }
    };
    let display_type = |name: &str| match name {
        "tensor_layout" => Some("typing.Literal['ports', 'schoonschip', 'call']"),
        "index_style" => Some("typing.Literal['alphabet', 'graph', 'raw']"),
        "component_style" => Some("typing.Literal['superscript', 'array']"),
        "invariant_style" => Some("typing.Literal['compact', 'explicit']"),
        "tensor_view" => Some("typing.Literal['interactive', 'matrix']"),
        _ => None,
    };
    for class in module.class.values_mut() {
        if class.name == "TensorExpression" {
            let overloads = class.methods.get_mut("__mul__").expect("tensor product");
            let mut extension = overloads[0].clone();
            extension.parameters.positional_or_keyword[0].type_info = TypeInfo {
                name: "symbolica.core._ExpressionProduct[_ProductT]".into(),
                import: ["symbolica.core".into()].into(),
            };
            extension.r#return = named("_ProductT");
            overloads.insert(0, extension);
        }
        if class.name == "FactorProjector" {
            class.bases = vec![named("typing.Generic[_Projected]")];
            for name in ["symmetric", "antisymmetric", "cyclic"] {
                let overloads = class
                    .methods
                    .get_mut(name)
                    .expect("factor projector constructor");
                if overloads.len() != 1 {
                    continue;
                }
                let mut symbolic = overloads[0].clone();
                symbolic.is_overload = true;
                symbolic.parameters.varargs.as_mut().unwrap().type_info = named(SYMBOLIC_FACTORS);
                symbolic.r#return = named("FactorProjector[TensorExpression]");
                let mut mixed = symbolic.clone();
                mixed.parameters.varargs.as_mut().unwrap().type_info = named(ALL_FACTORS);
                mixed.r#return =
                    named("FactorProjector[TensorExpression] | FactorProjector[TensorNetwork]");
                let mut concrete = mixed.clone();
                concrete.parameters.varargs.as_mut().unwrap().type_info = named(&format!(
                    "typing.Unpack[tuple[{CONCRETE_FACTORS}, *tuple[{ALL_FACTORS}, ...]]]"
                ));
                concrete.r#return = named("FactorProjector[TensorNetwork]");
                let mut second = concrete.clone();
                second.parameters.varargs.as_mut().unwrap().type_info = named(&format!(
                    "typing.Unpack[tuple[{ALL_FACTORS}, {CONCRETE_FACTORS}, *tuple[{ALL_FACTORS}, ...]]]"
                ));
                *overloads = vec![symbolic, concrete, second, mixed];
            }
        }
        for member in &mut class.attrs {
            simplify(&mut member.r#type);
        }
        for (name, (getter, setter)) in &mut class.getter_setters {
            for member in getter.iter_mut().chain(setter.iter_mut()) {
                simplify(&mut member.r#type);
                if matches!(class.name, "Tensor" | "TensorExpression" | "TensorNetwork")
                    && name == "axes"
                {
                    member.doc = python_doc!("Tensor.axes");
                }
                if class.name == "DisplaySettings"
                    && let Some(value) = display_type(name)
                {
                    member.r#type = named(value);
                }
                if matches!(class.name, "TensorEvaluator" | "CompiledTensorEvaluator")
                    && name == "output_shape"
                {
                    member.r#type = named("tuple[int, ...]");
                }
            }
        }
        for method in class.methods.values_mut().flatten() {
            // These established tensor operations intentionally specialize the
            // base Expression API: indexing, powers, and HEP display settings.
            // Keep their actual signatures instead of claiming base compatibility.
            if class.name == "TensorExpression"
                && matches!(
                    method.name,
                    "__call__" | "__pow__" | "to_typst" | "formatted" | "_repr_html_"
                )
            {
                method.type_ignored = Some(pyo3_stub_gen::type_info::IgnoreTarget::Specified(&[
                    "override",
                    "invalid-method-override",
                ]));
            }
            // PyO3 accepts macro doc literals; the stub derive only reads literal
            // attributes. Reuse the same docs for every manually specified overload.
            match (class.name, method.name) {
                ("Representation", "__call__") => {
                    method.doc = python_doc!("Representation.__call__")
                }
                ("TensorName", "__call__") => method.doc = python_doc!("TensorName.__call__"),
                ("BroadcastFunction", "__call__") => {
                    method.doc = python_doc!("BroadcastFunction.__call__")
                }
                ("Tensor", "__getitem__") => method.doc = python_doc!("Tensor.__getitem__"),
                ("Tensor", "__setitem__") => method.doc = python_doc!("Tensor.__setitem__"),
                ("Tensor", "__iter__") => method.doc = python_doc!("Tensor.__iter__"),
                ("TensorExpression", "__getitem__") => {
                    method.doc = python_doc!("TensorExpression.__getitem__")
                }
                ("TensorExpression", "__add__") => {
                    method.doc = python_doc!("TensorExpression.__add__")
                }
                ("TensorExpression", "__radd__") => {
                    method.doc = python_doc!("TensorExpression.__radd__")
                }
                ("TensorExpression", "__sub__") => {
                    method.doc = python_doc!("TensorExpression.__sub__")
                }
                ("TensorExpression", "__rsub__") => {
                    method.doc = python_doc!("TensorExpression.__rsub__")
                }
                ("TensorExpression", "__mul__") if method.r#return.name == "_ProductT" => {
                    method.doc = python_doc!("TensorExpression.__mul_extension__")
                }
                ("TensorExpression", "__mul__") => {
                    method.doc = python_doc!("TensorExpression.__mul__")
                }
                ("TensorExpression", "__rmul__") => {
                    method.doc = python_doc!("TensorExpression.__rmul__")
                }
                ("TensorExpression", "__truediv__") => {
                    method.doc = python_doc!("TensorExpression.__truediv__")
                }
                ("TensorExpression", "__rtruediv__") => {
                    method.doc = python_doc!("TensorExpression.__rtruediv__")
                }
                ("TensorExpression", "outer") => method.doc = python_doc!("TensorExpression.outer"),
                ("TensorExpression", "contract_ports") => {
                    method.doc = python_doc!("TensorExpression.contract_ports")
                }
                ("TensorExpression", "compose") => {
                    method.doc = python_doc!("TensorExpression.compose")
                }
                _ => {}
            }

            simplify(&mut method.r#return);
            parameters(&mut method.parameters, |parameter| {
                simplify(&mut parameter.type_info);
                if class.name == "TensorPattern" && matches!(parameter.name, "factors" | "indices")
                {
                    parameter.type_info = named("_ScalarInput");
                }
                if matches!(class.name, "TensorExpression" | "Tensor" | "TensorNetwork") {
                    if parameter.name == "indices" {
                        parameter.type_info = named("_IndexInput");
                    }
                    if method.name == "rename_indices" && parameter.name == "mapping" {
                        parameter.type_info =
                            named("dict[int | str | Expression | Slot, _IndexInput]");
                    }
                    if class.name == "Tensor" && parameter.name == "dtype" {
                        parameter.type_info = named(if method.name == "map_components" {
                            "type[float] | type[complex] | type[Expression] | None"
                        } else {
                            "type[float] | type[complex] | type[Expression]"
                        });
                    }
                }
                if class.name == "TensorExpression" {
                    match (method.name, parameter.name) {
                        ("__new__", "expression")
                        | ("__pow__", "exponent")
                        | ("__rpow__", "base") => {
                            parameter.type_info = named("TensorExpression | _ScalarInput");
                        }
                        ("replace", "pattern") => {
                            parameter.type_info =
                                named("TensorRule | list[TensorRule] | _ScalarInput");
                        }
                        _ => {}
                    }
                }
                if class.name == "TensorRule" && method.name == "__new__" {
                    parameter.type_info = match parameter.name {
                        "pattern" => named("TensorExpression | _ScalarInput | HeldExpression"),
                        "rhs" => named(
                            "TensorExpression | _ScalarInput | HeldExpression | typing.Callable[[dict[Expression, Expression]], TensorExpression | _ScalarInput]",
                        ),
                        _ => parameter.type_info.clone(),
                    };
                }
                if class.name == "DisplaySettings"
                    && let Some(value) = display_type(parameter.name)
                {
                    parameter.type_info = named(value);
                }
            });
        }
        if class.name == "Tensor"
            && let Some(overloads) = class.methods.get_mut("__getitem__")
            && !overloads
                .iter()
                .any(|method| method.r#return.name == "_Components")
        {
            let mut slice = overloads.last().expect("Tensor component overload").clone();
            slice.parameters.positional_or_keyword[0].type_info = named("tuple[int | slice, ...]");
            slice.r#return = named("_Components");
            slice.doc = python_doc!("Tensor.__getitem__");
            overloads.push(slice);
        }
        if class.name == "TensorLibrary"
            && let Some(overloads) = class.methods.get_mut("get")
            && overloads.len() == 1
        {
            let mut with_default = overloads[0].clone();
            with_default.is_overload = true;
            with_default.r#return = named("Tensor | _LibraryDefault");
            parameters(&mut with_default.parameters, |parameter| {
                if parameter.name == "default" {
                    parameter.type_info = named("_LibraryDefault");
                    parameter.default = ParameterDefault::None;
                }
            });
            let mut without_default = with_default.clone();
            without_default.r#return = named("Tensor | None");
            parameters(&mut without_default.parameters, |parameter| {
                if parameter.name == "default" {
                    parameter.type_info = named("None");
                    parameter.default = ParameterDefault::Expr("None".to_owned());
                }
            });
            *overloads = vec![without_default, with_default];
        }
    }
    for function in module.function.values_mut().flatten() {
        match function.name {
            "dot" => function.doc = python_doc!("dot"),
            "chain" => function.doc = python_doc!("chain"),
            "trace" => function.doc = python_doc!("trace"),
            _ => {}
        }
        simplify(&mut function.r#return);
        parameters(&mut function.parameters, |parameter| {
            simplify(&mut parameter.type_info);
            if matches!(function.name, "chain" | "trace") {
                if parameter.name == "factors" {
                    parameter.type_info = named(if function.r#return.name == "TensorExpression" {
                        SYMBOLIC_FACTORS
                    } else {
                        ALL_FACTORS
                    });
                } else if matches!(parameter.name, "factor" | "second") {
                    parameter.type_info = named(CONCRETE_FACTORS);
                }
            }
        });
    }
    // Module variables are emitted after actual imports by the upstream generator.
    // Quoted aliases support forward references without scanning Python docstrings.
    module.variables.insert(
        "_ProductT",
        VariableDef {
            name: "_ProductT",
            type_: named("typing.TypeVar"),
            default: Some("typing.TypeVar(\"_ProductT\")".to_owned()),
        },
    );
    module.variables.insert(
        "_Projected",
        VariableDef {
            name: "_Projected",
            type_: named("typing.TypeVar"),
            default: Some("typing.TypeVar(\"_Projected\", \"TensorExpression\", \"TensorNetwork\", covariant=True)".to_owned()),
        },
    );
    module.variables.insert(
        "_LibraryDefault",
        VariableDef {
            name: "_LibraryDefault",
            type_: named("typing.TypeVar"),
            default: Some("typing.TypeVar(\"_LibraryDefault\")".to_owned()),
        },
    );
    for (name, value) in [
        ("_RealInput", "Float | int | float | str | decimal.Decimal"),
        (
            "_ScalarInput",
            "Expression | int | float | complex | str | decimal.Decimal | Float | ComplexFloat | tuple[_RealInput, _RealInput]",
        ),
        (
            "_ReplacementInput",
            "_ScalarInput | HeldExpression | typing.Callable[[dict[Expression, Expression]], Expression]",
        ),
        ("_IndexInput", "int | str | Expression | Slot | _AutoIndex"),
        (
            "_Components",
            "Expression | float | complex | list[_Components]",
        ),
    ] {
        module.variables.insert(
            name,
            VariableDef {
                name,
                type_: named("typing.TypeAlias"),
                default: Some(format!("{value:?}")),
            },
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn public_overloads_keep_runtime_defaults_and_python_type_names() {
        let info = crate::stub_info().unwrap();
        let module = &info.modules["symbolica.community.tensor"];
        let library = module
            .class
            .values()
            .find(|class| class.name == "TensorLibrary")
            .unwrap();
        let get = &library.methods["get"];
        assert_eq!(get.len(), 2);
        let default = get[0]
            .parameters
            .positional_or_keyword
            .iter()
            .find(|parameter| parameter.name == "default")
            .unwrap();
        assert_eq!(default.type_info.name, "None");
        assert_eq!(default.default, ParameterDefault::Expr("None".to_owned()));
        let source = module.to_string();
        assert!(!source.contains("PythonExpression"));
        assert!(!source.contains("ConvertibleToExpression"));
    }
}
