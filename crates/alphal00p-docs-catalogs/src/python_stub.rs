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
        let callable = if trimmed.starts_with("class ") {
            Some((alphal00p_docs_schema::DocItemKind::PythonClass, "class "))
        } else if trimmed.starts_with("def ") {
            Some((alphal00p_docs_schema::DocItemKind::PythonFunction, "def "))
        } else if trimmed.starts_with("async def ") {
            Some((
                alphal00p_docs_schema::DocItemKind::PythonFunction,
                "async def ",
            ))
        } else {
            None
        };

        if let Some((kind, prefix)) = callable {
            let name = declaration_name(trimmed, prefix);
            if name.is_empty() || (!exports.is_empty() && !exports.iter().any(|item| item == &name))
            {
                continue;
            }
            let (signature, signature_end) = python_signature(&lines, index);
            let (docs, member_start) = python_docstring(&lines, signature_end + 1, 0);
            let members = if kind == alphal00p_docs_schema::DocItemKind::PythonClass {
                python_class_members(
                    &lines,
                    member_start,
                    0,
                    &docs,
                    signature.contains("enum.Enum"),
                )
            } else {
                vec![]
            };
            declarations.push(PythonDeclaration {
                name,
                signature: python_display_signature(&signature),
                docs,
                kind,
                line: index as u32 + 1,
                members,
            });
            continue;
        }

        let Some((name, annotation)) = (!exports.is_empty())
            .then(|| split_top_level_once(trimmed, ':'))
            .flatten()
        else {
            continue;
        };
        let name = name.trim();
        if !is_python_identifier(name)
            || !exports.iter().any(|exported| exported == name)
            || declarations
                .iter()
                .any(|declaration| declaration.name == name)
        {
            continue;
        }
        let annotation = split_top_level_once(annotation, '=')
            .map_or(annotation, |(annotation, _)| annotation)
            .trim();
        if annotation.is_empty() {
            continue;
        }
        let annotation = python_display_signature(annotation);
        declarations.push(PythonDeclaration {
            name: name.to_owned(),
            signature: format!("{name}: {annotation}"),
            docs: format!("Exported constant of type `{annotation}`."),
            kind: alphal00p_docs_schema::DocItemKind::PythonConstant,
            line: index as u32 + 1,
            members: vec![],
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

fn python_class_members(
    lines: &[&str],
    start: usize,
    class_indent: usize,
    class_docs: &str,
    is_enum: bool,
) -> Vec<DocMember> {
    let member_indent = class_indent + 4;
    let variant_docs = python_named_docs(class_docs, "Variants", true);
    let mut members: Vec<DocMember> = vec![];
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
            let kind = python_callable_kind(&name, decorators.iter().map(String::as_str));
            let is_overload = decorators
                .iter()
                .any(|decorator| matches!(decorator.as_str(), "@overload" | "@typing.overload"));
            let mut member = DocMember::new(&name, kind);
            member.signature = Some(python_display_signature(&signature));
            if !docs.is_empty() {
                member.docs = Some(DocText::new(DocFormat::PythonDocstring, &docs));
            }
            member.members = python_parameters(&signature, &docs);
            if is_overload {
                member.kind = DocMemberKind::Overload;
                if let Some(callable) = members
                    .iter_mut()
                    .rev()
                    .find(|callable| callable.name == name && callable.kind == kind)
                {
                    callable.members.push(member);
                } else {
                    let mut callable = DocMember::new(name, kind);
                    callable.members.push(member);
                    members.push(callable);
                }
            } else if let Some(callable) = members.iter_mut().rev().find(|callable| {
                callable.name == name
                    && callable.kind == kind
                    && callable
                        .members
                        .iter()
                        .any(|member| member.kind == DocMemberKind::Overload)
            }) {
                callable.signature = member.signature;
                callable.docs = member.docs;
                callable.members.extend(member.members);
            } else {
                members.push(member);
            }
            decorators.clear();
            index = next.max(signature_end + 1);
            continue;
        }
        let field_kind = if is_enum {
            DocMemberKind::Variant
        } else {
            DocMemberKind::Field
        };
        if let Some(mut field) = python_field(trimmed, field_kind) {
            let (docs, next) = python_docstring(lines, index + 1, class_indent);
            if !docs.is_empty() {
                field.docs = Some(DocText::new(DocFormat::PythonDocstring, docs));
            } else if let Some(docs) = variant_docs.get(&field.name) {
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

fn python_field(line: &str, kind: DocMemberKind) -> Option<DocMember> {
    if line.starts_with(['#', '@']) || line.starts_with("r\"\"\"") || line.starts_with("\"\"\"") {
        return None;
    }
    let (name, _) = split_top_level_once(line, ':').or_else(|| {
        (kind == DocMemberKind::Variant)
            .then(|| split_top_level_once(line, '='))
            .flatten()
    })?;
    let name = name.trim();
    if !is_python_identifier(name) {
        return None;
    }
    let mut member = DocMember::new(name, kind);
    let default = split_top_level_once(line, '=')
        .map(|(_, default)| default.trim().to_owned())
        .filter(|default| !default.is_empty());
    if kind != DocMemberKind::Variant || default.as_deref() != Some("...") {
        member.signature = Some(line.to_owned());
        member.default = default;
    }
    Some(member)
}

fn python_parameters(signature: &str, docs: &str) -> Vec<DocMember> {
    let Some(opening) = signature.find('(') else {
        return vec![];
    };
    let Some(closing) = matching_delimiter(signature, opening, '(', ')') else {
        return vec![];
    };
    let mut parameter_docs = python_named_docs(docs, "Parameters", false);
    parameter_docs.extend(python_named_docs(docs, "Arguments", false));
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
            let (raw_name, annotation) = split_top_level_once(without_default, ':')
                .map_or((without_default, ""), |(name, annotation)| {
                    (name.trim(), annotation.trim())
                });
            let name = raw_name.trim_start_matches('*').trim();
            if matches!(name, "self" | "cls") || !is_python_identifier(name) {
                return None;
            }
            let mut member = DocMember::new(name, DocMemberKind::Parameter);
            member.signature =
                (!annotation.is_empty()).then(|| python_display_signature(annotation));
            member.default = default;
            if let Some(docs) = parameter_docs.get(name) {
                member.docs = Some(DocText::new(DocFormat::PythonDocstring, docs));
            }
            Some(member)
        })
        .collect()
}

fn python_display_signature(signature: &str) -> String {
    signature.replace("builtins.", "").replace("typing.", "")
}

fn python_named_docs(
    docs: &str,
    section: &str,
    inline_description: bool,
) -> BTreeMap<String, String> {
    let lines = docs.lines().map(str::trim).collect::<Vec<_>>();
    let numpy_start = lines
        .windows(2)
        .position(|pair| pair[0] == section && is_section_rule(pair[1]))
        .map(|index| index + 2);
    let markdown_start = lines
        .iter()
        .position(|line| {
            line.starts_with('#')
                && line.trim_start_matches('#').trim().trim_end_matches(':') == section
        })
        .map(|index| index + 1);
    let Some((mut index, markdown)) = numpy_start
        .map(|index| (index, false))
        .or_else(|| markdown_start.map(|index| (index, true)))
    else {
        return BTreeMap::new();
    };
    let mut result = BTreeMap::new();
    while index < lines.len() {
        if markdown && lines[index].starts_with('#')
            || !markdown && index + 1 < lines.len() && is_section_rule(lines[index + 1])
        {
            break;
        }
        let entry = if markdown {
            lines[index]
                .strip_prefix("- ")
                .and_then(|line| line.split_once(':'))
                .and_then(|(name, description)| {
                    let name = name.trim();
                    is_python_identifier(name).then(|| (vec![name], description.trim()))
                })
        } else {
            lines[index]
                .split_once(" : ")
                .and_then(|(names, annotation)| {
                    let names = names
                        .split(',')
                        .map(str::trim)
                        .filter(|name| is_python_identifier(name))
                        .collect::<Vec<_>>();
                    (!names.is_empty())
                        .then_some((names, if inline_description { annotation } else { "" }))
                })
        };
        let Some((names, initial_description)) = entry else {
            index += 1;
            continue;
        };
        index += 1;
        let mut description = (!initial_description.is_empty())
            .then_some(initial_description)
            .into_iter()
            .collect::<Vec<_>>();
        while index < lines.len() {
            if markdown && (lines[index].starts_with("- ") || lines[index].starts_with('#'))
                || !markdown
                    && (lines[index].contains(" : ")
                        || index + 1 < lines.len() && is_section_rule(lines[index + 1]))
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
    let doc_indent = indentation(lines[index]);
    let trimmed = lines[index].trim();
    let Some((delimiter, body)) = triple_quoted_body(trimmed) else {
        return (String::new(), start);
    };
    if let Some(body) = body.strip_suffix(delimiter) {
        return (deduplicate_adjacent_paragraphs(body.trim()), index + 1);
    }
    let mut docs = vec![];
    if !body.is_empty() {
        docs.push(body.trim().to_owned());
    }
    index += 1;
    while index < lines.len() {
        let line = strip_docstring_indent(lines[index], doc_indent);
        if let Some(body) = line.trim_end().strip_suffix(delimiter) {
            if !body.trim().is_empty() {
                docs.push(body.trim_end().to_owned());
            }
            return (finish_docstring(docs), index + 1);
        }
        docs.push(line.to_owned());
        index += 1;
    }
    (finish_docstring(docs), index)
}

fn strip_docstring_indent(line: &str, indent: usize) -> &str {
    if line.len() >= indent
        && line.as_bytes()[..indent]
            .iter()
            .all(u8::is_ascii_whitespace)
    {
        &line[indent..]
    } else {
        line.trim_start()
    }
}

fn finish_docstring(mut lines: Vec<String>) -> String {
    while lines.first().is_some_and(|line| line.trim().is_empty()) {
        lines.remove(0);
    }
    while lines.last().is_some_and(|line| line.trim().is_empty()) {
        lines.pop();
    }
    deduplicate_adjacent_paragraphs(&lines.join("\n"))
}

fn deduplicate_adjacent_paragraphs(docs: &str) -> String {
    let mut rendered = Vec::new();
    let mut previous_prose = None::<String>;
    let mut inside_fence = false;
    for paragraph in docs.split("\n\n") {
        let normalized = paragraph.trim();
        if normalized.is_empty() {
            continue;
        }
        let contains_fence = paragraph.lines().any(|line| line.trim().starts_with("```"));
        let protected = inside_fence || contains_fence;
        if !protected && previous_prose.as_deref() == Some(normalized) {
            continue;
        }
        rendered.push(paragraph.trim_end());
        previous_prose = (!protected).then(|| normalized.to_owned());
        if paragraph
            .lines()
            .filter(|line| line.trim().starts_with("```"))
            .count()
            % 2
            == 1
        {
            inside_fence = !inside_fence;
        }
    }
    rendered.join("\n\n")
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

    use alphal00p_docs_schema::{ApiLanguage, DocItemKind, DocMemberKind};

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
    fn parses_each_exported_annotated_constant_once() {
        let stub = r#"
__all__ = [
    "AUTO",
    "Engine",
]

AUTO: typing.Final[Auto] = ...
PRIVATE: int = 1

class Engine:
    """Runs work."""

AUTO: typing.Final[Auto]
"#;
        let declarations = python_declarations(stub).unwrap();
        assert_eq!(
            declarations
                .iter()
                .map(|declaration| declaration.name.as_str())
                .collect::<Vec<_>>(),
            ["AUTO", "Engine"]
        );
        assert_eq!(declarations[0].kind, DocItemKind::PythonConstant);
        assert_eq!(declarations[0].signature, "AUTO: Final[Auto]");
        assert_eq!(
            declarations[0].docs,
            "Exported constant of type `Final[Auto]`."
        );
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
    fn parses_markdown_parameter_and_argument_descriptions() {
        let stub = r#"
class Representation:
    def __new__(cls, name: str, dimension: int | Expression | str, is_self_dual: bool = True) -> Representation:
        r"""
        Create and register a new representation.

        # Parameters:
        - name: String name for the representation
        - dimension: Size of the representation
          as an integer or symbolic expression.
        - is_self_dual: If True, creates self-dual representation; if False, creates dualizable pair

        # Examples:
        ```python
        Representation("Euclidean", 4, is_self_dual=True)
        ```
        """

    def rename(self, value):
        """
        Return a renamed representation.

        # Arguments:
        - value: New representation name.
          This parameter deliberately has no annotation.
        """
"#;
        let declarations = python_declarations(stub).unwrap();
        let constructor = &declarations[0].members[0];
        assert_eq!(constructor.name, "__new__");
        assert_eq!(constructor.members[1].name, "dimension");
        assert_eq!(
            constructor.members[1].docs.as_ref().unwrap().body,
            "Size of the representation as an integer or symbolic expression."
        );
        let self_dual = &constructor.members[2];
        assert_eq!(self_dual.name, "is_self_dual");
        assert_eq!(self_dual.signature.as_deref(), Some("bool"));
        assert_eq!(self_dual.default.as_deref(), Some("True"));
        assert_eq!(
            self_dual.docs.as_ref().unwrap().body,
            "If True, creates self-dual representation; if False, creates dualizable pair"
        );

        let value = &declarations[0].members[1].members[0];
        assert_eq!(value.name, "value");
        assert!(value.signature.is_none());
        assert_eq!(
            value.docs.as_ref().unwrap().body,
            "New representation name. This parameter deliberately has no annotation."
        );
    }

    #[test]
    fn parses_enum_assignments_as_documented_variants() {
        let stub = r#"
class ExecutionMode(enum.Enum):
    """
    Execution modes.

    Variants
    --------
    Single : Execute one contraction at a time
    Scalar : Only contract scalar operations
    All : Execute every possible contraction
    """
    Single = ...
    Scalar = ...
    All = 3

class SymbolicParallelism(enum.Enum):
    """Symbolic parallelism policy."""
    Auto = ...
    """Select with workload heuristics."""
    Serial = ...
    """Remain on the calling thread."""

class LogLevel(enum.Enum):
    Off = ...
    """A level lower than all log levels."""
    Error = ...
    """The error log level."""
"#;
        let declarations = python_declarations(stub).unwrap();

        let execution = &declarations[0].members;
        assert_eq!(
            execution
                .iter()
                .map(|member| (member.name.as_str(), member.kind))
                .collect::<Vec<_>>(),
            [
                ("Single", DocMemberKind::Variant),
                ("Scalar", DocMemberKind::Variant),
                ("All", DocMemberKind::Variant),
            ]
        );
        assert!(execution[0].signature.is_none());
        assert!(execution[0].default.is_none());
        assert_eq!(execution[2].default.as_deref(), Some("3"));
        assert_eq!(
            execution[2].docs.as_ref().unwrap().body,
            "Execute every possible contraction"
        );
        assert_eq!(
            declarations[1].members[0].docs.as_ref().unwrap().body,
            "Select with workload heuristics."
        );
        assert_eq!(
            declarations[2].members[0].docs.as_ref().unwrap().body,
            "A level lower than all log levels."
        );
    }

    #[test]
    fn coalesces_python_overloads_under_one_callable() {
        let stub = r#"
class Engine:
    @typing.overload
    def run(self, value: int) -> int:
        """
        Run an integer.

        Parameters
        ----------
        value : int
            Integer input.
        """
    @overload
    def run(self, value: str, strict: bool = False) -> str:
        """Run a string."""
"#;
        let declarations = python_declarations(stub).unwrap();
        let runs = declarations[0]
            .members
            .iter()
            .filter(|member| member.name == "run")
            .collect::<Vec<_>>();
        assert_eq!(runs.len(), 1);
        let run = runs[0];
        assert_eq!(run.kind, DocMemberKind::Method);
        assert!(run.signature.is_none());
        assert_eq!(run.members.len(), 2);
        assert!(
            run.members
                .iter()
                .all(|member| member.kind == DocMemberKind::Overload)
        );
        assert_eq!(
            run.members[0].signature.as_deref(),
            Some("def run(self, value: int) -> int:")
        );
        assert_eq!(run.members[0].members[0].name, "value");
        assert_eq!(
            run.members[0].members[0].docs.as_ref().unwrap().body,
            "Integer input."
        );
        assert_eq!(run.members[1].members[1].name, "strict");
        assert_eq!(run.members[1].members[1].default.as_deref(), Some("False"));
    }

    #[test]
    fn display_signatures_hide_stub_generator_qualification() {
        let declarations = python_declarations(
            "class Result:\n    pass\n\nclass Engine:\n    def run(self, values: typing.Sequence[builtins.float]) -> typing.Optional[Result]: ...\n",
        )
        .unwrap();
        let run = &declarations[1].members[0];
        assert_eq!(
            run.signature.as_deref(),
            Some("def run(self, values: Sequence[float]) -> Optional[Result]:")
        );
        assert_eq!(run.members[0].signature.as_deref(), Some("Sequence[float]"));
    }

    #[test]
    fn exact_adjacent_prose_is_deduplicated_without_touching_distinct_or_code_blocks() {
        let declarations = python_declarations(
            r#"class Engine:
    r"""
    Reuse one compiled engine.

    Reuse one compiled engine.

    Reuse one compiled engine for each batch.

    ```python
    print("same")

    print("same")
    ```
    """
"#,
        )
        .unwrap();
        assert_eq!(
            declarations[0].docs,
            "Reuse one compiled engine.\n\nReuse one compiled engine for each batch.\n\n```python\nprint(\"same\")\n\nprint(\"same\")\n```"
        );
    }

    #[test]
    fn preserves_indentation_in_numpy_style_example_blocks() {
        let declarations = python_declarations(
            r#"class Engine:
    r"""
    Run an engine.

    Examples
    --------
    Iterate over grouped values::

        for name, values in groups.items():
            for value in values:
                print(name, value)
    """
"#,
        )
        .unwrap();

        assert!(declarations[0].docs.contains(
            "Iterate over grouped values::\n\n    for name, values in groups.items():\n        for value in values:\n            print(name, value)"
        ));
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

        let linnet = &catalogs["linnet-py"].root.scopes["exports"].items["Graph"];
        assert!(linnet.members.iter().any(|member| {
            member.name == "global_data" && member.kind == DocMemberKind::Getter
        }));
        assert!(linnet.members.iter().any(|member| {
            member.name == "global_data" && member.kind == DocMemberKind::Setter
        }));
        assert!(linnet.members.iter().any(|member| {
            member.name == "from_dot" && member.kind == DocMemberKind::AssociatedFunction
        }));
        let auto = &catalogs["linnet-py"].root.scopes["exports"].items["AUTO"];
        assert_eq!(auto.kind, DocItemKind::PythonConstant);
        assert_eq!(
            auto.docs.as_ref().unwrap().body,
            "Exported constant of type `Auto`."
        );

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

        let log_level = &catalogs["gammaloop-python"].root.scopes["exports"].items["LogLevel"];
        assert_eq!(
            log_level
                .members
                .iter()
                .filter(|member| member.kind == DocMemberKind::Variant)
                .map(|member| member.name.as_str())
                .collect::<Vec<_>>(),
            ["Off", "Error", "Warn", "Info", "Debug", "Trace"]
        );
        let execution = &catalogs["spynso3"].root.scopes["exports"].items["ExecutionMode"];
        assert_eq!(
            execution.members[0].docs.as_ref().unwrap().body,
            "Select one smallest-degree rewrite per step; without `n_steps`, continue until no work remains"
        );
        let parallelism = &catalogs["spynso3"].root.scopes["exports"].items["SymbolicParallelism"];
        assert_eq!(
            parallelism
                .members
                .iter()
                .map(|member| member.name.as_str())
                .collect::<Vec<_>>(),
            ["Auto", "Serial", "Parallel"]
        );
    }
}
