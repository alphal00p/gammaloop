use std::collections::BTreeMap;

use alphal00p_docs_schema::{DocFormat, DocMember, DocMemberKind, DocText};
use eyre::{Result, ensure};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PythonDeclaration {
    pub(crate) name: String,
    pub(crate) signature: String,
    pub(crate) docs: String,
    pub(crate) kind: alphal00p_docs_schema::DocItemKind,
    pub(crate) line: u32,
    pub(crate) members: Vec<DocMember>,
}

pub(crate) fn python_declarations(source: &str) -> Result<Vec<PythonDeclaration>> {
    let lines = source.lines().collect::<Vec<_>>();
    let exports = python_exports(&lines);
    let mut declarations = vec![];
    for (index, line) in lines.iter().enumerate() {
        if indentation(line) != 0 {
            continue;
        }
        let trimmed = line.trim();
        let (kind, prefix) = if trimmed.starts_with("class ") {
            (alphal00p_docs_schema::DocItemKind::PythonClass, "class ")
        } else if trimmed.starts_with("def ") {
            (alphal00p_docs_schema::DocItemKind::PythonFunction, "def ")
        } else if trimmed.starts_with("async def ") {
            (
                alphal00p_docs_schema::DocItemKind::PythonFunction,
                "async def ",
            )
        } else {
            continue;
        };
        let name = declaration_name(trimmed, prefix);
        if name.is_empty() || (!exports.is_empty() && !exports.iter().any(|item| item == &name)) {
            continue;
        }
        let (signature, signature_end) = python_signature(&lines, index);
        let (docs, member_start) = python_docstring(&lines, signature_end + 1, 0);
        let members = if kind == alphal00p_docs_schema::DocItemKind::PythonClass {
            python_class_members(&lines, member_start, 0)
        } else {
            vec![]
        };
        declarations.push(PythonDeclaration {
            name,
            signature,
            docs,
            kind,
            line: index as u32 + 1,
            members,
        });
    }
    if !exports.is_empty() {
        declarations.sort_by_key(|declaration| {
            exports
                .iter()
                .position(|name| name == &declaration.name)
                .unwrap_or(usize::MAX)
        });
        for exported in exports {
            ensure!(
                declarations
                    .iter()
                    .any(|declaration| declaration.name == exported),
                "stub exports {exported} but has no top-level declaration"
            );
        }
    }
    Ok(declarations)
}

fn python_class_members(lines: &[&str], start: usize, class_indent: usize) -> Vec<DocMember> {
    let member_indent = class_indent + 4;
    let mut members = vec![];
    let mut decorators = vec![];
    let mut index = start;
    while index < lines.len() {
        let line = lines[index];
        let trimmed = line.trim();
        if trimmed.is_empty() {
            index += 1;
            continue;
        }
        let indent = indentation(line);
        if indent <= class_indent {
            break;
        }
        if indent != member_indent {
            index += 1;
            continue;
        }
        if trimmed.starts_with('@') {
            decorators.push(trimmed.to_owned());
            index += 1;
            continue;
        }
        let callable_prefix = if trimmed.starts_with("def ") {
            Some("def ")
        } else if trimmed.starts_with("async def ") {
            Some("async def ")
        } else {
            None
        };
        if let Some(prefix) = callable_prefix {
            let name = declaration_name(trimmed, prefix);
            let (signature, signature_end) = python_signature(lines, index);
            let (docs, next) = python_docstring(lines, signature_end + 1, member_indent);
            let mut member = DocMember::new(
                &name,
                python_callable_kind(&name, decorators.iter().map(String::as_str)),
            );
            member.signature = Some(signature.clone());
            if !docs.is_empty() {
                member.docs = Some(DocText::new(DocFormat::PythonDocstring, &docs));
            }
            member.members = python_parameters(&signature, &docs);
            members.push(member);
            decorators.clear();
            index = next.max(signature_end + 1);
            continue;
        }
        if let Some(mut field) = python_field(trimmed) {
            let (docs, next) = python_docstring(lines, index + 1, member_indent);
            if !docs.is_empty() {
                field.docs = Some(DocText::new(DocFormat::PythonDocstring, docs));
            }
            members.push(field);
            decorators.clear();
            index = next.max(index + 1);
            continue;
        }
        decorators.clear();
        index += 1;
    }
    members
}

fn python_callable_kind<'a>(
    name: &str,
    decorators: impl Iterator<Item = &'a str>,
) -> DocMemberKind {
    let decorators = decorators.collect::<Vec<_>>();
    if decorators
        .iter()
        .any(|decorator| decorator.ends_with(".setter"))
    {
        DocMemberKind::Setter
    } else if decorators
        .iter()
        .any(|decorator| *decorator == "@property" || decorator.ends_with(".getter"))
    {
        DocMemberKind::Getter
    } else if name == "__new__"
        || decorators
            .iter()
            .any(|decorator| matches!(*decorator, "@classmethod" | "@staticmethod"))
    {
        DocMemberKind::AssociatedFunction
    } else {
        DocMemberKind::Method
    }
}

fn python_field(line: &str) -> Option<DocMember> {
    if line.starts_with(['#', '@']) || line.starts_with("r\"\"\"") || line.starts_with("\"\"\"") {
        return None;
    }
    let (name, _) = split_top_level_once(line, ':')?;
    let name = name.trim();
    if !is_python_identifier(name) {
        return None;
    }
    let mut member = DocMember::new(name, DocMemberKind::Field);
    member.signature = Some(line.to_owned());
    member.default = split_top_level_once(line, '=')
        .map(|(_, default)| default.trim().to_owned())
        .filter(|default| !default.is_empty());
    Some(member)
}

fn python_parameters(signature: &str, docs: &str) -> Vec<DocMember> {
    let Some(opening) = signature.find('(') else {
        return vec![];
    };
    let Some(closing) = matching_delimiter(signature, opening, '(', ')') else {
        return vec![];
    };
    let parameter_docs = python_parameter_docs(docs);
    split_top_level(&signature[opening + 1..closing], ',')
        .into_iter()
        .filter_map(|parameter| {
            let parameter = parameter.trim();
            if parameter.is_empty() || matches!(parameter, "/" | "*") {
                return None;
            }
            let (without_default, default) = split_top_level_once(parameter, '=')
                .map_or((parameter, None), |(value, default)| {
                    (value.trim(), Some(default.trim().to_owned()))
                });
            let (raw_name, _) = split_top_level_once(without_default, ':')
                .map_or((without_default, ""), |(name, annotation)| {
                    (name.trim(), annotation.trim())
                });
            let name = raw_name.trim_start_matches('*').trim();
            if matches!(name, "self" | "cls") || !is_python_identifier(name) {
                return None;
            }
            let mut member = DocMember::new(name, DocMemberKind::Parameter);
            member.signature = Some(parameter.to_owned());
            member.default = default;
            if let Some(docs) = parameter_docs.get(name) {
                member.docs = Some(DocText::new(DocFormat::PythonDocstring, docs));
            }
            Some(member)
        })
        .collect()
}

fn python_parameter_docs(docs: &str) -> BTreeMap<String, String> {
    let lines = docs.lines().map(str::trim).collect::<Vec<_>>();
    let Some(mut index) = lines
        .windows(2)
        .position(|pair| pair[0] == "Parameters" && is_section_rule(pair[1]))
        .map(|index| index + 2)
    else {
        return BTreeMap::new();
    };
    let mut result = BTreeMap::new();
    while index < lines.len() {
        if index + 1 < lines.len() && is_section_rule(lines[index + 1]) {
            break;
        }
        let Some((names, _annotation)) = lines[index].split_once(" : ") else {
            index += 1;
            continue;
        };
        let names = names
            .split(',')
            .map(str::trim)
            .filter(|name| is_python_identifier(name))
            .collect::<Vec<_>>();
        if names.is_empty() {
            index += 1;
            continue;
        }
        index += 1;
        let mut description = vec![];
        while index < lines.len() {
            if lines[index].contains(" : ")
                || (index + 1 < lines.len() && is_section_rule(lines[index + 1]))
            {
                break;
            }
            if !lines[index].is_empty() {
                description.push(lines[index]);
            }
            index += 1;
        }
        let description = description.join(" ");
        if !description.is_empty() {
            for name in names {
                result.insert(name.to_owned(), description.clone());
            }
        }
    }
    result
}

fn is_section_rule(line: &str) -> bool {
    !line.is_empty() && line.chars().all(|character| character == '-')
}

fn python_exports(lines: &[&str]) -> Vec<String> {
    let mut exports = vec![];
    let mut inside = false;
    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with("__all__") && trimmed.contains('[') {
            inside = true;
        }
        if inside {
            if let Some(value) = trimmed
                .strip_prefix('"')
                .and_then(|value| value.split('"').next())
            {
                exports.push(value.to_owned());
            }
            if trimmed.contains(']') {
                break;
            }
        }
    }
    exports
}

fn python_signature(lines: &[&str], start: usize) -> (String, usize) {
    let mut depth = 0_i32;
    let mut signature = vec![];
    let mut end = start;
    for (index, line) in lines.iter().enumerate().skip(start) {
        let trimmed = line.trim();
        for character in trimmed.chars() {
            match character {
                '(' | '[' | '{' => depth += 1,
                ')' | ']' | '}' => depth -= 1,
                _ => {}
            }
        }
        end = index;
        signature.push(if depth <= 0 && trimmed.contains(':') {
            trimmed.strip_suffix(" ...").unwrap_or(trimmed)
        } else {
            trimmed
        });
        if depth <= 0 && trimmed.contains(':') {
            break;
        }
    }
    (signature.join(" "), end)
}

fn python_docstring(lines: &[&str], start: usize, parent_indent: usize) -> (String, usize) {
    let mut index = start;
    while index < lines.len() && lines[index].trim().is_empty() {
        index += 1;
    }
    if index >= lines.len() || indentation(lines[index]) <= parent_indent {
        return (String::new(), start);
    }
    let trimmed = lines[index].trim();
    let Some((delimiter, body)) = triple_quoted_body(trimmed) else {
        return (String::new(), start);
    };
    if let Some(body) = body.strip_suffix(delimiter) {
        return (body.trim().to_owned(), index + 1);
    }
    let mut docs = vec![];
    if !body.is_empty() {
        docs.push(body.trim());
    }
    index += 1;
    while index < lines.len() {
        let trimmed = lines[index].trim();
        if let Some(body) = trimmed.strip_suffix(delimiter) {
            if !body.is_empty() {
                docs.push(body.trim());
            }
            return (docs.join("\n").trim().to_owned(), index + 1);
        }
        docs.push(trimmed);
        index += 1;
    }
    (docs.join("\n").trim().to_owned(), index)
}

fn triple_quoted_body(line: &str) -> Option<(&'static str, &str)> {
    for (prefix, delimiter) in [("r\"\"\"", "\"\"\""), ("\"\"\"", "\"\"\"")] {
        if let Some(body) = line.strip_prefix(prefix) {
            return Some((delimiter, body));
        }
    }
    None
}

fn declaration_name(line: &str, prefix: &str) -> String {
    line[prefix.len()..]
        .split(['(', ':'])
        .next()
        .unwrap_or_default()
        .trim()
        .to_owned()
}

fn indentation(line: &str) -> usize {
    line.len() - line.trim_start().len()
}

fn is_python_identifier(value: &str) -> bool {
    let mut characters = value.chars();
    characters
        .next()
        .is_some_and(|character| character == '_' || character.is_ascii_alphabetic())
        && characters.all(|character| character == '_' || character.is_ascii_alphanumeric())
}

fn matching_delimiter(source: &str, opening: usize, open: char, close: char) -> Option<usize> {
    let mut depth = 0_u32;
    let mut quote = None;
    let mut escaped = false;
    for (offset, character) in source[opening..].char_indices() {
        if let Some(active_quote) = quote {
            if escaped {
                escaped = false;
            } else if character == '\\' {
                escaped = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        if matches!(character, '\'' | '"') {
            quote = Some(character);
        } else if character == open {
            depth += 1;
        } else if character == close {
            depth -= 1;
            if depth == 0 {
                return Some(opening + offset);
            }
        }
    }
    None
}

fn split_top_level(source: &str, delimiter: char) -> Vec<&str> {
    let mut result = vec![];
    let mut start = 0;
    let mut depth = 0_i32;
    let mut quote = None;
    let mut escaped = false;
    for (index, character) in source.char_indices() {
        if let Some(active_quote) = quote {
            if escaped {
                escaped = false;
            } else if character == '\\' {
                escaped = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => quote = Some(character),
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' => depth -= 1,
            _ if character == delimiter && depth == 0 => {
                result.push(&source[start..index]);
                start = index + character.len_utf8();
            }
            _ => {}
        }
    }
    result.push(&source[start..]);
    result
}

fn split_top_level_once(source: &str, delimiter: char) -> Option<(&str, &str)> {
    let mut depth = 0_i32;
    let mut quote = None;
    let mut escaped = false;
    for (index, character) in source.char_indices() {
        if let Some(active_quote) = quote {
            if escaped {
                escaped = false;
            } else if character == '\\' {
                escaped = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => quote = Some(character),
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' => depth -= 1,
            _ if character == delimiter && depth == 0 => {
                return Some((&source[..index], &source[index + character.len_utf8()..]));
            }
            _ => {}
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use alphal00p_docs_schema::{ApiLanguage, DocMemberKind};

    use super::*;
    use crate::{CatalogRequest, export_catalog, first_sentence};

    #[test]
    fn parses_python_exports_in_declared_order() {
        let stub = r#"
__all__ = [
    "run",
    "Engine",
]

class Engine:
    r"""
    Runs work. More detail.
    """

def run(value: int) -> int:
    r"""
    Return one value.
    """
"#;
        let declarations = python_declarations(stub).unwrap();
        assert_eq!(declarations[0].name, "run");
        assert_eq!(declarations[1].name, "Engine");
        assert_eq!(first_sentence(&declarations[1].docs), "Runs work.");
    }

    #[test]
    fn parses_methods_properties_fields_and_parameters() {
        let stub = r#"
class Engine:
    """Runs work."""
    retries: int = 2
    @property
    def status(self) -> str:
        """Current engine status."""
    @status.setter
    def status(self, value: str) -> None: ...
    @classmethod
    def open(cls, path: str, mode: str = "r") -> Engine:
        """
        Open an engine.

        Parameters
        ----------
        path : str
            Input path.
        mode : str
            File mode.
        """
"#;
        let declarations = python_declarations(stub).unwrap();
        let members = &declarations[0].members;
        assert_eq!(members[0].kind, DocMemberKind::Field);
        assert_eq!(members[0].default.as_deref(), Some("2"));
        assert_eq!(members[1].kind, DocMemberKind::Getter);
        assert_eq!(
            members[1].docs.as_ref().unwrap().body,
            "Current engine status."
        );
        assert_eq!(members[2].kind, DocMemberKind::Setter);
        assert_eq!(members[2].members[0].name, "value");
        assert_eq!(members[3].kind, DocMemberKind::AssociatedFunction);
        assert!(
            members[3]
                .signature
                .as_ref()
                .unwrap()
                .ends_with("-> Engine:")
        );
        assert_eq!(members[3].members[0].name, "path");
        assert_eq!(
            members[3].members[0].docs.as_ref().unwrap().body,
            "Input path."
        );
        assert_eq!(members[3].members[1].default.as_deref(), Some("\"r\""));
    }

    #[test]
    fn existing_snapshots_generate_deterministic_nested_catalogs() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        let components = [
            ("gammaloop-python", "gammaloop._gammaloop"),
            ("linnet-py", "linnet_py"),
            ("spynso3", "symbolica.community.spenso"),
            ("idenso-community", "symbolica.community.idenso"),
            ("vakint-community", "symbolica.community.vakint"),
        ];
        let mut catalogs = BTreeMap::new();
        for (component, module) in components {
            let request = CatalogRequest {
                product_id: "test".to_owned(),
                product_title: "Test".to_owned(),
                component_id: component.to_owned(),
                package: component.to_owned(),
                component_title: component.to_owned(),
                version: "0.1.0".to_owned(),
                language: ApiLanguage::Python,
                module: Some(module.to_owned()),
                features: vec![],
            };
            let stub = root
                .join("docs/api/python")
                .join(format!("{component}.pyi"));
            let first = export_catalog(&request, root, Some(&stub)).unwrap();
            let second = export_catalog(&request, root, Some(&stub)).unwrap();
            assert_eq!(
                serde_json::to_vec(&first).unwrap(),
                serde_json::to_vec(&second).unwrap()
            );
            catalogs.insert(component, first);
        }

        let vakint = &catalogs["vakint-community"].root.scopes["exports"].items["Vakint"];
        let canonical = vakint
            .members
            .iter()
            .find(|member| member.name == "to_canonical")
            .unwrap();
        assert_eq!(canonical.kind, DocMemberKind::Method);
        assert!(
            canonical
                .signature
                .as_ref()
                .unwrap()
                .ends_with("-> Expression:")
        );
        assert!(!canonical.docs.as_ref().unwrap().body.is_empty());
        assert_eq!(canonical.members[0].name, "integral_expression");
        assert_eq!(canonical.members[1].name, "short_form");
        assert_eq!(canonical.members[1].default.as_deref(), Some("None"));

        let linnet = &catalogs["linnet-py"].root.scopes["exports"].items["DotGraph"];
        assert!(linnet.members.iter().any(|member| {
            member.name == "global_data" && member.kind == DocMemberKind::Getter
        }));
        assert!(linnet.members.iter().any(|member| {
            member.name == "global_data" && member.kind == DocMemberKind::Setter
        }));
        assert!(linnet.members.iter().any(|member| {
            member.name == "from_string" && member.kind == DocMemberKind::AssociatedFunction
        }));

        let spenso = &catalogs["spynso3"].root.scopes["exports"].items["Tensor"];
        let sparse = spenso
            .members
            .iter()
            .find(|member| member.name == "sparse")
            .unwrap();
        assert_eq!(sparse.kind, DocMemberKind::AssociatedFunction);
        assert_eq!(
            sparse
                .members
                .iter()
                .map(|parameter| parameter.name.as_str())
                .collect::<Vec<_>>(),
            ["structure", "type_info"]
        );

        let gamma = &catalogs["gammaloop-python"].root.scopes["exports"].items["GammaLoopAPI"];
        let evaluate = gamma
            .members
            .iter()
            .find(|member| member.name == "evaluate_sample")
            .unwrap();
        assert!(
            evaluate
                .signature
                .as_ref()
                .unwrap()
                .ends_with("-> EvaluationResult:")
        );
        assert_eq!(evaluate.members[0].name, "point");
    }
}
