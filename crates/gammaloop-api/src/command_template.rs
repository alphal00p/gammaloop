use std::collections::{BTreeMap, BTreeSet};

use color_eyre::{
    eyre::{eyre, Result},
    Section,
};

pub type CommandEnvironment = BTreeMap<String, String>;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct PlaceholderSpec {
    pub name: String,
    pub default: Option<String>,
}

impl PlaceholderSpec {
    fn parse(inner: &str) -> Option<Self> {
        let (name, default) = inner
            .split_once(':')
            .map(|(name, default)| (name, Some(default)))
            .unwrap_or((inner, None));
        if !is_valid_placeholder_name(name) {
            return None;
        }
        Some(Self {
            name: name.to_string(),
            default: default.map(str::to_string),
        })
    }

    pub fn render(&self) -> String {
        match &self.default {
            Some(default) => format!("$({}:{default})", self.name),
            None => format!("$({})", self.name),
        }
    }
}

pub fn is_valid_placeholder_name(name: &str) -> bool {
    let mut chars = name.chars();
    let Some(first) = chars.next() else {
        return false;
    };
    matches!(first, 'A'..='Z' | 'a'..='z' | '_')
        && chars.all(|ch| matches!(ch, 'A'..='Z' | 'a'..='z' | '0'..='9' | '_'))
}

pub fn placeholder_specs(input: &str) -> BTreeSet<PlaceholderSpec> {
    let mut specs = BTreeSet::new();
    let mut index = 0usize;

    while let Some(relative_start) = input[index..].find("$(") {
        let start = index + relative_start;
        let spec_start = start + 2;
        let Some(relative_end) = input[spec_start..].find(')') else {
            break;
        };
        let end = spec_start + relative_end;
        if let Some(spec) = PlaceholderSpec::parse(&input[spec_start..end]) {
            specs.insert(spec);
        }
        index = end + 1;
    }

    specs
}

pub fn contains_placeholder(input: &str) -> bool {
    !placeholder_specs(input).is_empty()
}

pub fn expand(input: &str, environment: &CommandEnvironment) -> Result<String> {
    let mut expanded = String::with_capacity(input.len());
    let mut index = 0usize;

    while let Some(relative_start) = input[index..].find("$(") {
        let start = index + relative_start;
        expanded.push_str(&input[index..start]);

        let spec_start = start + 2;
        let Some(relative_end) = input[spec_start..].find(')') else {
            expanded.push_str(&input[start..]);
            index = input.len();
            break;
        };
        let end = spec_start + relative_end;
        let Some(spec) = PlaceholderSpec::parse(&input[spec_start..end]) else {
            expanded.push_str(&input[start..=end]);
            index = end + 1;
            continue;
        };

        let value = environment
            .get(&spec.name)
            .or(spec.default.as_ref())
            .ok_or_else(|| {
                let available = if environment.is_empty() {
                    "no variables are defined".to_string()
                } else {
                    format!(
                        "available variables: {}",
                        environment.keys().cloned().collect::<Vec<_>>().join(", ")
                    )
                };
                eyre!("Missing command-block variable '{}'", spec.name)
                    .with_note(|| {
                        format!("The placeholder {} appears in '{input}'.", spec.render())
                    })
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
    use super::{expand, placeholder_specs, CommandEnvironment};

    #[test]
    fn placeholder_specs_find_names_and_defaults() {
        let rendered = placeholder_specs("run $(process_name) -g $(graph_id:0) $(tag:a:b)")
            .into_iter()
            .map(|spec| spec.render())
            .collect::<Vec<_>>();

        assert_eq!(
            rendered,
            vec![
                "$(graph_id:0)".to_string(),
                "$(process_name)".to_string(),
                "$(tag:a:b)".to_string()
            ]
        );
    }

    #[test]
    fn placeholder_names_ignore_numeric_old_braces_and_invalid_forms() {
        let names =
            placeholder_specs("generate [{1}] {{2}} {process_name} $(process_name) $(4bad) $(")
                .into_iter()
                .map(|spec| spec.name)
                .collect::<Vec<_>>();
        assert_eq!(names, vec!["process_name".to_string()]);
    }

    #[test]
    fn expand_replaces_valid_names_uses_defaults_and_preserves_other_braces() {
        let environment = CommandEnvironment::from([
            ("process".to_string(), "aa_aa".to_string()),
            ("graph_id".to_string(), "0".to_string()),
        ]);

        assert_eq!(
            expand(
                "run $(process) -g $(graph_id:1) --tag $(tag:default) [{1}] {{2}}",
                &environment
            )
            .unwrap(),
            "run aa_aa -g 0 --tag default [{1}] {{2}}"
        );
    }

    #[test]
    fn expand_preserves_invalid_and_legacy_placeholder_forms() {
        let environment = CommandEnvironment::from([("process".to_string(), "aa_aa".to_string())]);

        assert_eq!(
            expand(
                "literal {process} invalid $(1process) actual $(process)",
                &environment
            )
            .unwrap(),
            "literal {process} invalid $(1process) actual aa_aa"
        );
    }

    #[test]
    fn expand_reports_missing_variables_without_defaults() {
        let error = expand("run $(process)", &CommandEnvironment::new()).unwrap_err();
        let error_text = format!("{error:?}");

        assert!(error_text.contains("Missing command-block variable 'process'"));
    }
}
