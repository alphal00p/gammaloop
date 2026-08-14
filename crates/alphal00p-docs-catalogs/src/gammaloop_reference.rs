//! GammaLoop CLI and settings metadata exported from the compiled application types.

use std::collections::BTreeSet;

use clap::{Arg, Command, CommandFactory};
use eyre::{Result, ensure};
use gammaloop_api::{CLISettings, OneShot};
use gammalooprs::{settings::RuntimeSettings, utils::serde_utils::ShowDefaultsGuard};
use schemars::schema_for;
use serde::Serialize;
use serde_json::Value;

use crate::generated::{
    CliArgument, CliCommand, GENERATED_REFERENCE_SCHEMA, GammaLoopReference, GeneratedSource,
    SettingReference,
};

pub fn export() -> Result<GammaLoopReference> {
    let mut command = OneShot::command();
    command.build();
    let mut commands = Vec::new();
    collect_commands(&command, "", &mut commands);

    let cli_schema = serde_json::to_value(schema_for!(CLISettings))?;
    let runtime_schema = serde_json::to_value(schema_for!(RuntimeSettings))?;
    let cli_defaults = settings_defaults(&CLISettings::default())?;
    let runtime_defaults = settings_defaults(&RuntimeSettings::default())?;
    let mut settings = Vec::new();
    collect_settings(
        "cli",
        &cli_schema,
        &cli_schema,
        Some(&cli_defaults),
        &mut settings,
        0,
    );
    collect_settings(
        "runtime",
        &runtime_schema,
        &runtime_schema,
        Some(&runtime_defaults),
        &mut settings,
        0,
    );

    ensure!(!commands.is_empty(), "Clap exported no GammaLoop commands");
    ensure!(
        !settings.is_empty(),
        "Schemars exported no GammaLoop settings"
    );
    ensure_unique(
        commands.iter().map(|command| command.path.as_str()),
        "command",
    )?;
    ensure_unique(
        settings.iter().map(|setting| setting.path.as_str()),
        "setting",
    )?;

    Ok(GammaLoopReference {
        schema_version: GENERATED_REFERENCE_SCHEMA,
        sources: vec![
            GeneratedSource {
                kind: "clap-command-factory".to_owned(),
                path: "crates/gammaloop-api/src/lib.rs".to_owned(),
                symbol: "gammaloop_api::OneShot::command".to_owned(),
            },
            GeneratedSource {
                kind: "schemars-and-serde-default".to_owned(),
                path: "crates/gammaloop-api/src/lib.rs".to_owned(),
                symbol: "gammaloop_api::CLISettings".to_owned(),
            },
            GeneratedSource {
                kind: "schemars-and-serde-default".to_owned(),
                path: "crates/gammalooprs/src/settings/mod.rs".to_owned(),
                symbol: "gammalooprs::settings::RuntimeSettings".to_owned(),
            },
        ],
        commands,
        settings,
    })
}

fn collect_commands(command: &Command, parent: &str, output: &mut Vec<CliCommand>) {
    let path = if parent.is_empty() {
        command.get_name().to_owned()
    } else {
        format!("{parent} {}", command.get_name())
    };
    output.push(CliCommand {
        path: path.clone(),
        name: command.get_name().to_owned(),
        about: command
            .get_long_about()
            .or_else(|| command.get_about())
            .map(ToString::to_string)
            .unwrap_or_default(),
        aliases: command.get_visible_aliases().map(str::to_owned).collect(),
        arguments: command.get_arguments().map(argument).collect(),
    });
    for child in command.get_subcommands() {
        collect_commands(child, &path, output);
    }
}

fn argument(argument: &Arg) -> CliArgument {
    CliArgument {
        id: argument.get_id().to_string(),
        long: argument.get_long().map(str::to_owned),
        short: argument.get_short(),
        value_names: argument
            .get_value_names()
            .into_iter()
            .flatten()
            .map(ToString::to_string)
            .collect(),
        help: argument
            .get_long_help()
            .or_else(|| argument.get_help())
            .map(ToString::to_string)
            .unwrap_or_default(),
        required: argument.is_required_set(),
        defaults: argument
            .get_default_values()
            .iter()
            .map(|value| value.to_string_lossy().into_owned())
            .collect(),
        possible_values: argument
            .get_possible_values()
            .into_iter()
            .filter(|value| !value.is_hide_set())
            .map(|value| value.get_name().to_owned())
            .collect(),
    }
}

fn settings_defaults(value: &impl Serialize) -> Result<Value> {
    let _show_defaults = ShowDefaultsGuard::new(true);
    Ok(serde_json::to_value(value)?)
}

fn collect_settings(
    prefix: &str,
    root: &Value,
    node: &Value,
    defaults: Option<&Value>,
    output: &mut Vec<SettingReference>,
    depth: usize,
) {
    if depth > 32 {
        return;
    }
    let Some(object) = object_schema(root, node) else {
        output.push(setting_reference(prefix, root, node, defaults, false));
        return;
    };
    let Some(properties) = object.get("properties").and_then(Value::as_object) else {
        output.push(setting_reference(prefix, root, node, defaults, false));
        return;
    };
    if properties.is_empty() {
        output.push(setting_reference(prefix, root, node, defaults, false));
        return;
    }
    let required = object
        .get("required")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(Value::as_str)
        .collect::<BTreeSet<_>>();
    for (name, child) in properties {
        let path = format!("{prefix}.{name}");
        let child_default = defaults.and_then(|value| value.get(name));
        if object_schema(root, child)
            .and_then(|value| value.get("properties"))
            .and_then(Value::as_object)
            .is_some_and(|properties| !properties.is_empty())
        {
            collect_settings(&path, root, child, child_default, output, depth + 1);
        } else {
            output.push(setting_reference(
                &path,
                root,
                child,
                child_default,
                required.contains(name.as_str()),
            ));
        }
    }
}

fn setting_reference(
    path: &str,
    root: &Value,
    node: &Value,
    default: Option<&Value>,
    required: bool,
) -> SettingReference {
    SettingReference {
        path: path.to_owned(),
        value_type: schema_type(root, node),
        description: schema_description(root, node),
        required,
        default: default
            .cloned()
            .or_else(|| resolved_schema(root, node).get("default").cloned()),
        possible_values: schema_values(root, node),
    }
}

fn object_schema<'a>(root: &'a Value, node: &'a Value) -> Option<&'a Value> {
    let node = resolved_schema(root, node);
    if node.get("properties").and_then(Value::as_object).is_some() {
        return Some(node);
    }
    for key in ["allOf", "anyOf", "oneOf"] {
        if let Some(values) = node.get(key).and_then(Value::as_array) {
            if let Some(object) = values.iter().find_map(|value| object_schema(root, value)) {
                return Some(object);
            }
        }
    }
    None
}

fn resolved_schema<'a>(root: &'a Value, node: &'a Value) -> &'a Value {
    let Some(reference) = node.get("$ref").and_then(Value::as_str) else {
        return node;
    };
    reference
        .strip_prefix('#')
        .and_then(|pointer| root.pointer(pointer))
        .unwrap_or(node)
}

fn schema_description(root: &Value, node: &Value) -> String {
    node.get("description")
        .and_then(Value::as_str)
        .or_else(|| {
            resolved_schema(root, node)
                .get("description")
                .and_then(Value::as_str)
        })
        .unwrap_or_default()
        .to_owned()
}

fn schema_type(root: &Value, node: &Value) -> String {
    let node = resolved_schema(root, node);
    if let Some(value) = node.get("type") {
        return match value {
            Value::String(value) => value.clone(),
            Value::Array(values) => values
                .iter()
                .filter_map(Value::as_str)
                .collect::<Vec<_>>()
                .join(" | "),
            _ => "unknown".to_owned(),
        };
    }
    for key in ["anyOf", "oneOf", "allOf"] {
        if let Some(values) = node.get(key).and_then(Value::as_array) {
            let kinds = values
                .iter()
                .map(|value| schema_type(root, value))
                .filter(|value| value != "unknown")
                .collect::<BTreeSet<_>>();
            if !kinds.is_empty() {
                return kinds.into_iter().collect::<Vec<_>>().join(" | ");
            }
        }
    }
    if node.get("enum").is_some() {
        "enum".to_owned()
    } else {
        "unknown".to_owned()
    }
}

fn schema_values(root: &Value, node: &Value) -> Vec<String> {
    let node = resolved_schema(root, node);
    let mut values = BTreeSet::new();
    if let Some(items) = node.get("enum").and_then(Value::as_array) {
        values.extend(items.iter().map(json_label));
    }
    if let Some(value) = node.get("const") {
        values.insert(json_label(value));
    }
    for key in ["anyOf", "oneOf", "allOf"] {
        if let Some(items) = node.get(key).and_then(Value::as_array) {
            for item in items {
                values.extend(schema_values(root, item));
            }
        }
    }
    values.into_iter().collect()
}

fn json_label(value: &Value) -> String {
    value
        .as_str()
        .map(str::to_owned)
        .unwrap_or_else(|| value.to_string())
}

fn ensure_unique<'a>(values: impl Iterator<Item = &'a str>, kind: &str) -> Result<()> {
    let mut seen = BTreeSet::new();
    for value in values {
        ensure!(seen.insert(value), "duplicate generated {kind} {value}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::export;

    #[test]
    fn compiled_cli_and_settings_are_exported() {
        let reference = export().unwrap();
        assert!(
            reference
                .commands
                .iter()
                .any(|command| command.path == "gammaLoop")
        );
        assert!(
            reference
                .settings
                .iter()
                .any(|setting| setting.path.starts_with("runtime.general."))
        );
        assert!(
            reference
                .settings
                .iter()
                .any(|setting| setting.default.is_some())
        );
    }
}
