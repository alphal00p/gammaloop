use std::collections::{BTreeMap, BTreeSet};

use color_eyre::{
    eyre::{eyre, Result},
    Section,
};

pub type CommandEnvironment = BTreeMap<String, String>;

pub fn is_valid_placeholder_name(name: &str) -> bool {
    let mut chars = name.chars();
    let Some(first) = chars.next() else {
        return false;
    };
    matches!(first, 'A'..='Z' | 'a'..='z' | '_')
        && chars.all(|ch| matches!(ch, 'A'..='Z' | 'a'..='z' | '0'..='9' | '_'))
}

pub fn placeholder_names(input: &str) -> BTreeSet<String> {
    let mut names = BTreeSet::new();
    let mut index = 0usize;

    while let Some(relative_start) = input[index..].find('{') {
        let start = index + relative_start;
        if input[start..].starts_with("{{") {
            index = start + 2;
            continue;
        }

        let name_start = start + 1;
        let Some(relative_end) = input[name_start..].find('}') else {
            break;
        };
        let end = name_start + relative_end;
        let name = &input[name_start..end];
        if is_valid_placeholder_name(name) {
            names.insert(name.to_string());
        }
        index = end + 1;
    }

    names
}

pub fn contains_placeholder(input: &str) -> bool {
    !placeholder_names(input).is_empty()
}

pub fn expand(input: &str, environment: &CommandEnvironment) -> Result<String> {
    let mut expanded = String::with_capacity(input.len());
    let mut index = 0usize;

    while let Some(relative_start) = input[index..].find('{') {
        let start = index + relative_start;
        expanded.push_str(&input[index..start]);

        if input[start..].starts_with("{{") {
            if let Some(relative_end) = input[start + 2..].find("}}") {
                let inner_start = start + 2;
                let inner_end = inner_start + relative_end;
                let inner = &input[inner_start..inner_end];
                if is_valid_placeholder_name(inner) {
                    expanded.push('{');
                    expanded.push_str(inner);
                    expanded.push('}');
                    index = inner_end + 2;
                    continue;
                }
            }
            expanded.push('{');
            index = start + 1;
            continue;
        }

        let name_start = start + 1;
        let Some(relative_end) = input[name_start..].find('}') else {
            expanded.push_str(&input[start..]);
            index = input.len();
            break;
        };
        let end = name_start + relative_end;
        let name = &input[name_start..end];
        if !is_valid_placeholder_name(name) {
            expanded.push_str(&input[start..=end]);
            index = end + 1;
            continue;
        }

        let value = environment.get(name).ok_or_else(|| {
            let available = if environment.is_empty() {
                "no variables are defined".to_string()
            } else {
                format!(
                    "available variables: {}",
                    environment.keys().cloned().collect::<Vec<_>>().join(", ")
                )
            };
            eyre!("Missing command-block variable '{name}'")
                .with_note(|| format!("The placeholder {{{name}}} appears in '{input}'."))
                .with_note(|| available)
        })?;
        expanded.push_str(value);
        index = end + 1;
    }

    if index < input.len() {
        expanded.push_str(&input[index..]);
    }

    Ok(expanded)
}

#[cfg(test)]
mod tests {
    use super::{expand, placeholder_names, CommandEnvironment};

    #[test]
    fn placeholder_names_ignore_numeric_and_escaped_braces() {
        assert_eq!(
            placeholder_names("generate [{1}] {{2}} {process_name} {{literal}}"),
            ["process_name".to_string()].into_iter().collect()
        );
    }

    #[test]
    fn expand_replaces_valid_names_and_preserves_other_braces() {
        let environment = CommandEnvironment::from([
            ("process".to_string(), "aa_aa".to_string()),
            ("graph_id".to_string(), "0".to_string()),
        ]);

        assert_eq!(
            expand("run {process} -g {graph_id} [{1}] {{2}}", &environment).unwrap(),
            "run aa_aa -g 0 [{1}] {{2}}"
        );
    }

    #[test]
    fn expand_can_escape_named_placeholders() {
        let environment = CommandEnvironment::from([("process".to_string(), "aa_aa".to_string())]);

        assert_eq!(
            expand("literal {{process}} actual {process}", &environment).unwrap(),
            "literal {process} actual aa_aa"
        );
    }
}
