use std::fs;

use pyo3::{
    prelude::*,
    types::{PyDict, PyModule},
};

const STUB_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/python/symbolica/community/feynkit/__init__.pyi"
);

fn main() -> pyo3_stub_gen::Result<()> {
    let info = feynkit_py::stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.feynkit")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "FeynKit did not contribute a symbolica.community.feynkit stub module",
            )
        })?;
    let stub = normalize_stub_source(&module.to_string());
    validate_stub(&stub)?;
    fs::write(STUB_PATH, stub)?;
    Ok(())
}

fn normalize_stub_source(source: &str) -> String {
    let mut lines = source
        .lines()
        .map(|line| line.trim_end_matches([' ', '\t']))
        .collect::<Vec<_>>();
    while lines.last().is_some_and(|line| line.is_empty()) {
        lines.pop();
    }
    lines.join("\n") + "\n"
}

fn validate_stub(source: &str) -> PyResult<()> {
    const DOCUMENTATION_AUDIT: &str = r#"
import ast

tree = ast.parse(source)
errors = []

def section_body(doc, heading):
    lines = doc.splitlines()
    try:
        start = lines.index(heading) + 2
    except ValueError:
        return None
    end = len(lines)
    for index in range(start, len(lines) - 1):
        underline = lines[index + 1].strip()
        if lines[index].strip() and len(underline) >= 3 and set(underline) == {"-"}:
            end = index
            break
    return [line.strip() for line in lines[start:end] if line.strip()]

for node in ast.walk(tree):
    if isinstance(node, ast.ClassDef):
        doc = ast.get_docstring(node) or ""
        if "Examples\n--------" not in doc:
            errors.append(f"class {node.name} has no Examples section")
        if "Examples\n--------" in doc and "Parameters\n----------" in doc:
            if doc.index("Examples\n--------") > doc.index("Parameters\n----------"):
                errors.append(f"class {node.name} puts Parameters before Examples")
        parameters = section_body(doc, "Parameters")
        if parameters is not None and (not parameters or parameters == ["None"]):
            errors.append(f"class {node.name} has an empty Parameters section")
        continue
    if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        continue

    doc = ast.get_docstring(node) or ""
    is_property = any(
        isinstance(decorator, ast.Name) and decorator.id == "property"
        for decorator in node.decorator_list
    )
    arguments = [
        *node.args.posonlyargs,
        *node.args.args,
        *node.args.kwonlyargs,
    ]
    arguments = [argument for argument in arguments if argument.arg not in {"self", "cls"}]
    if node.args.vararg is not None:
        arguments.append(node.args.vararg)
    if node.args.kwarg is not None:
        arguments.append(node.args.kwarg)
    has_examples = "Examples\n--------" in doc
    has_parameters = "Parameters\n----------" in doc
    parameters = section_body(doc, "Parameters")

    if not is_property and not has_examples:
        errors.append(f"callable {node.name} has no Examples section")
    if not is_property and arguments and not has_parameters:
        errors.append(f"callable {node.name} has undocumented parameters")
    if not arguments and has_parameters:
        errors.append(f"callable {node.name} has an empty Parameters section")
    if parameters is not None and (not parameters or parameters == ["None"]):
        errors.append(f"callable {node.name} has an empty Parameters section")
    if has_examples and has_parameters and doc.index("Examples\n--------") > doc.index("Parameters\n----------"):
        errors.append(f"callable {node.name} puts Parameters before Examples")
    if "isinstance(" in doc:
        errors.append(f"callable {node.name} uses a type-check-only example")

if errors:
    raise AssertionError("invalid FeynKit API documentation:\n" + "\n".join(errors))
"#;

    Python::initialize();
    Python::attach(|py| {
        PyModule::import(py, "ast")?
            .getattr("parse")?
            .call1((source,))?;
        let globals = PyDict::new(py);
        globals.set_item("source", source)?;
        PyModule::import(py, "builtins")?.getattr("exec")?.call1((
            DOCUMENTATION_AUDIT,
            &globals,
            &globals,
        ))?;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::normalize_stub_source;

    #[test]
    fn normalizes_trailing_whitespace_and_end_of_file() {
        assert_eq!(
            normalize_stub_source("class Example:  \n    ...\t\n   \n\n"),
            "class Example:\n    ...\n"
        );
    }
}
