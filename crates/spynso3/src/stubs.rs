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

pub(crate) fn refine(module: &mut Module) {
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
    };
    let display_type = |name: &str| match name {
        "tensor_layout" => Some("typing.Literal['ports', 'schoonschip', 'call']"),
        "index_style" => Some("typing.Literal['alphabet', 'graph', 'raw']"),
        "component_style" => Some("typing.Literal['superscript', 'array']"),
        "tensor_view" => Some("typing.Literal['interactive', 'matrix']"),
        _ => None,
    };
    for class in module.class.values_mut() {
        for member in &mut class.attrs {
            simplify(&mut member.r#type);
        }
        for (name, (getter, setter)) in &mut class.getter_setters {
            for member in getter.iter_mut().chain(setter.iter_mut()) {
                simplify(&mut member.r#type);
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
            simplify(&mut method.r#return);
            parameters(&mut method.parameters, |parameter| {
                simplify(&mut parameter.type_info);
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
            slice.doc =
                "Select components by logical coordinates, returning nested lists for sliced axes.";
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
            without_default
                .parameters
                .positional_or_keyword
                .retain(|parameter| parameter.name != "default");
            *overloads = vec![without_default, with_default];
        }
    }
    for function in module.function.values_mut().flatten() {
        simplify(&mut function.r#return);
        parameters(&mut function.parameters, |parameter| {
            simplify(&mut parameter.type_info)
        });
    }
    // Module variables are emitted after actual imports by the upstream generator.
    // Quoted aliases support forward references without scanning Python docstrings.
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
