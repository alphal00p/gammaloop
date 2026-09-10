//! GammaLoop CLI and settings metadata exported from the compiled application types.

use std::collections::{BTreeMap, BTreeSet};

use clap::{Arg, ArgAction, Command, CommandFactory};
use eyre::{Result, bail, ensure};
use gammaloop_api::{CLISettings, OneShot, commands::CliArgumentMetadata};
use gammalooprs::{settings::RuntimeSettings, utils::serde_utils::ShowDefaultsGuard};
use schemars::schema_for;
use serde::Serialize;
use serde_json::Value;

use crate::generated::{
    CliAlias, CliAliasKind, CliArgument, CliArgumentAction, CliCommand, CliValueArity,
    GENERATED_REFERENCE_SCHEMA, GammaLoopReference, GeneratedSource, SettingReference,
};

pub fn export() -> Result<GammaLoopReference> {
    let declared_command = OneShot::command();
    let mut command = declared_command.clone();
    command.build();
    let mut commands = Vec::new();
    collect_commands(&command, Some(&declared_command), "", &mut commands)?;

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
    if let Some(command) = commands
        .iter()
        .find(|command| !command.hidden && !command.generated_help && command.about.is_empty())
    {
        bail!("public command `{}` has no description", command.path);
    }
    if let Some((command, argument)) = commands.iter().find_map(|command| {
        command
            .arguments
            .iter()
            .find(|argument| !argument.hidden && argument.help.is_empty())
            .map(|argument| (command, argument))
    }) {
        bail!(
            "public argument `{}` on `{}` has no description",
            argument.id,
            command.path
        );
    }
    if let Some(setting) = settings
        .iter()
        .find(|setting| setting.description.is_empty())
    {
        bail!("setting `{}` has no description", setting.path);
    }

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

fn collect_commands(
    command: &Command,
    declared_command: Option<&Command>,
    parent: &str,
    output: &mut Vec<CliCommand>,
) -> Result<()> {
    let path = if parent.is_empty() {
        command.get_name().to_owned()
    } else {
        format!("{parent} {}", command.get_name())
    };
    let declared_arguments = declared_command
        .into_iter()
        .flat_map(Command::get_arguments)
        .map(|argument| argument.get_id().as_str())
        .collect::<BTreeSet<_>>();
    let visible_aliases = command.get_visible_aliases().collect::<BTreeSet<_>>();
    let visible_short_flag_aliases = command
        .get_visible_short_flag_aliases()
        .collect::<BTreeSet<_>>();
    let visible_long_flag_aliases = command
        .get_visible_long_flag_aliases()
        .collect::<BTreeSet<_>>();
    let aliases = command
        .get_all_aliases()
        .map(|alias| CliAlias {
            name: alias.to_owned(),
            kind: CliAliasKind::Name,
            visible: visible_aliases.contains(alias),
        })
        .chain(command.get_all_short_flag_aliases().map(|alias| CliAlias {
            name: alias.to_string(),
            kind: CliAliasKind::ShortFlag,
            visible: visible_short_flag_aliases.contains(&alias),
        }))
        .chain(command.get_all_long_flag_aliases().map(|alias| CliAlias {
            name: alias.to_owned(),
            kind: CliAliasKind::LongFlag,
            visible: visible_long_flag_aliases.contains(alias),
        }))
        .collect();
    output.push(CliCommand {
        path: path.clone(),
        name: command.get_name().to_owned(),
        usage: command.clone().render_usage().to_string(),
        about: command
            .get_long_about()
            .or_else(|| command.get_about())
            .map(ToString::to_string)
            .unwrap_or_default(),
        hidden: command.is_hide_set(),
        generated_help: declared_command.is_none(),
        short_flag: command.get_short_flag(),
        long_flag: command.get_long_flag().map(str::to_owned),
        aliases,
        arguments: command
            .get_arguments()
            .map(|argument| {
                argument_metadata(
                    command,
                    argument,
                    argument.is_global_set()
                        && !declared_arguments.contains(argument.get_id().as_str()),
                )
            })
            .collect::<Result<_>>()?,
    });
    for child in command.get_subcommands() {
        let declared_child = declared_command.and_then(|declared| {
            declared
                .get_subcommands()
                .find(|candidate| candidate.get_name() == child.get_name())
        });
        collect_commands(child, declared_child, &path, output)?;
    }
    Ok(())
}

fn argument_metadata(command: &Command, argument: &Arg, inherited: bool) -> Result<CliArgument> {
    let documentation = argument.get::<CliArgumentMetadata>();
    let action = match argument.get_action() {
        ArgAction::Set => CliArgumentAction::Set,
        ArgAction::Append => CliArgumentAction::Append,
        ArgAction::SetTrue => CliArgumentAction::SetTrue,
        ArgAction::SetFalse => CliArgumentAction::SetFalse,
        ArgAction::Count => CliArgumentAction::Count,
        ArgAction::Help => CliArgumentAction::Help,
        ArgAction::HelpShort => CliArgumentAction::HelpShort,
        ArgAction::HelpLong => CliArgumentAction::HelpLong,
        ArgAction::Version => CliArgumentAction::Version,
        action => bail!(
            "unsupported Clap action {action:?} for {}",
            argument.get_id()
        ),
    };
    let arity = argument.get_num_args().unwrap_or_else(|| {
        if argument.get_action().takes_values() {
            1.into()
        } else {
            0.into()
        }
    });
    let visible_long_aliases = argument
        .get_visible_aliases()
        .into_iter()
        .flatten()
        .collect::<BTreeSet<_>>();
    let visible_short_aliases = argument
        .get_visible_short_aliases()
        .into_iter()
        .flatten()
        .collect::<BTreeSet<_>>();
    let aliases = argument
        .get_all_aliases()
        .into_iter()
        .flatten()
        .map(|alias| CliAlias {
            name: alias.to_owned(),
            kind: CliAliasKind::LongFlag,
            visible: visible_long_aliases.contains(alias),
        })
        .chain(
            argument
                .get_all_short_aliases()
                .into_iter()
                .flatten()
                .map(|alias| CliAlias {
                    name: alias.to_string(),
                    kind: CliAliasKind::ShortFlag,
                    visible: visible_short_aliases.contains(&alias),
                }),
        )
        .collect();
    let mut conflicts_with = command
        .get_arg_conflicts_with(argument)
        .into_iter()
        .map(|conflict| conflict.get_id().to_string())
        .collect::<BTreeSet<_>>();
    for candidate in command.get_arguments() {
        if candidate.get_id() != argument.get_id()
            && (candidate.is_exclusive_set()
                || command
                    .get_arg_conflicts_with(candidate)
                    .into_iter()
                    .any(|conflict| conflict.get_id() == argument.get_id()))
        {
            conflicts_with.insert(candidate.get_id().to_string());
        }
    }
    for group in command.get_groups() {
        let arguments = group.get_args().collect::<Vec<_>>();
        let mut group = group.clone();
        if !group.is_multiple() && arguments.contains(&argument.get_id()) {
            conflicts_with.extend(
                arguments
                    .into_iter()
                    .filter(|id| *id != argument.get_id())
                    .map(ToString::to_string),
            );
        }
    }
    if argument.is_exclusive_set() {
        conflicts_with.extend(
            command
                .get_arguments()
                .filter(|candidate| candidate.get_id() != argument.get_id())
                .map(|candidate| candidate.get_id().to_string()),
        );
    }

    Ok(CliArgument {
        id: argument.get_id().to_string(),
        long: argument.get_long().map(str::to_owned),
        short: argument.get_short(),
        aliases,
        action,
        arity: CliValueArity {
            min: arity.min_values(),
            max: (arity.max_values() != usize::MAX).then_some(arity.max_values()),
        },
        takes_values: arity.takes_values(),
        value_required: arity.min_values() > 0,
        value_names: if arity.takes_values() {
            argument
                .get_value_names()
                .into_iter()
                .flatten()
                .map(ToString::to_string)
                .collect()
        } else {
            Vec::new()
        },
        help: argument
            .get_long_help()
            .or_else(|| argument.get_help())
            .map(ToString::to_string)
            .unwrap_or_default(),
        required: argument.is_required_set(),
        positional: argument.is_positional(),
        index: argument.get_index(),
        hidden: argument.is_hide_set(),
        global: argument.is_global_set(),
        inherited,
        exclusive: argument.is_exclusive_set(),
        require_equals: argument.is_require_equals_set(),
        value_delimiter: argument.get_value_delimiter(),
        value_terminator: argument.get_value_terminator().map(ToString::to_string),
        conflicts_with: conflicts_with.into_iter().collect(),
        requires: documentation
            .into_iter()
            .flat_map(|metadata| metadata.requires.iter().copied())
            .map(str::to_owned)
            .collect(),
        defaults: argument
            .get_default_values()
            .iter()
            .map(|value| value.to_string_lossy().into_owned())
            .collect(),
        default_missing_values: documentation
            .into_iter()
            .flat_map(|metadata| metadata.default_missing_values.iter().copied())
            .map(str::to_owned)
            .collect(),
        possible_values: argument
            .get_possible_values()
            .into_iter()
            .filter(|value| !value.is_hide_set())
            .map(|value| value.get_name().to_owned())
            .collect(),
    })
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
        if object_schema(root, child).is_some_and(|value| {
            value
                .get("properties")
                .and_then(Value::as_object)
                .is_some_and(|properties| !properties.is_empty())
        }) {
            collect_settings(&path, root, child, child_default, output, depth + 1);
        } else {
            output.push(setting_reference(
                &path,
                root,
                child,
                child_default,
                required.contains(name.as_str()),
            ));
            if let Some(items) = array_items_schema(root, child)
                && object_schema(root, items).is_some_and(|value| {
                    value
                        .get("properties")
                        .and_then(Value::as_object)
                        .is_some_and(|properties| !properties.is_empty())
                })
            {
                collect_settings(&format!("{path}[]"), root, items, None, output, depth + 1);
            }
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
    let description = match (schema_description(root, node), path) {
        // Internally tagged enums synthesize their `type` properties in Schemars, so
        // there is no Rust field on which a documentation comment can live.
        (description, "runtime.kinematics.externals.type") if description.is_empty() => {
            "Selects how external momenta and helicities are supplied; currently `constant`."
                .to_owned()
        }
        (description, "runtime.subtraction.integrated_ct_settings.range.type")
            if description.is_empty() =>
        {
            "Selects the integrated-counterterm radial domain: `infinite` with damping or `compact`."
                .to_owned()
        }
        (description, "runtime.stability.rotation_axis[].type") if description.is_empty() => {
            "Selects a Cartesian quarter-turn, no rotation, or explicit Euler angles.".to_owned()
        }
        (description, _) => description,
    };
    SettingReference {
        path: path.to_owned(),
        value_type: schema_type(root, node),
        description,
        required,
        default: default
            .cloned()
            .or_else(|| resolved_schema(root, node).get("default").cloned()),
        possible_values: schema_values(root, node),
    }
}

fn object_schema(root: &Value, node: &Value) -> Option<Value> {
    let node = resolved_schema(root, node);
    let mut components = Vec::new();
    if let Some(properties) = node.get("properties").and_then(Value::as_object) {
        let mut object = serde_json::Map::new();
        object.insert("properties".to_owned(), Value::Object(properties.clone()));
        if let Some(required) = node.get("required").and_then(Value::as_array) {
            object.insert("required".to_owned(), Value::Array(required.clone()));
        }
        components.push(Value::Object(object));
    }

    if let Some(values) = node.get("allOf").and_then(Value::as_array) {
        components.extend(values.iter().filter_map(|value| object_schema(root, value)));
    }
    for key in ["anyOf", "oneOf"] {
        let Some(values) = node.get(key).and_then(Value::as_array) else {
            continue;
        };
        let alternatives = values
            .iter()
            .map(|value| {
                object_schema(root, value).unwrap_or_else(|| {
                    serde_json::json!({
                        "properties": {},
                        "required": [],
                    })
                })
            })
            .collect::<Vec<_>>();
        let has_object = alternatives.iter().any(|value| {
            value
                .get("properties")
                .and_then(Value::as_object)
                .is_some_and(|properties| !properties.is_empty())
        });
        if has_object && let Some(object) = merge_object_schemas(alternatives, false) {
            components.push(object);
        }
    }
    merge_object_schemas(components, true)
}

fn merge_object_schemas(objects: Vec<Value>, require_every_component: bool) -> Option<Value> {
    let mut properties = BTreeMap::<String, Value>::new();
    let mut required = None::<BTreeSet<String>>;
    for object in objects {
        let object_required = object
            .get("required")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
            .map(str::to_owned)
            .collect::<BTreeSet<_>>();
        required = Some(match required {
            None => object_required,
            Some(mut required) if require_every_component => {
                required.extend(object_required);
                required
            }
            Some(mut required) => {
                required.retain(|name| object_required.contains(name));
                required
            }
        });
        for (name, schema) in object
            .get("properties")
            .and_then(Value::as_object)
            .into_iter()
            .flatten()
        {
            match properties.entry(name.clone()) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(schema.clone());
                }
                std::collections::btree_map::Entry::Occupied(mut entry)
                    if entry.get() != schema =>
                {
                    let previous = std::mem::take(entry.get_mut());
                    entry.insert(serde_json::json!({ "anyOf": [previous, schema] }));
                }
                std::collections::btree_map::Entry::Occupied(_) => {}
            }
        }
    }
    if properties.is_empty() {
        return None;
    }
    let mut object = serde_json::Map::new();
    object.insert(
        "properties".to_owned(),
        Value::Object(properties.into_iter().collect()),
    );
    if let Some(required) = required.filter(|required| !required.is_empty()) {
        object.insert(
            "required".to_owned(),
            Value::Array(required.into_iter().map(Value::String).collect()),
        );
    }
    Some(Value::Object(object))
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
    let description = node
        .get("description")
        .and_then(Value::as_str)
        .filter(|description| !description.trim().is_empty());
    if let Some(description) = description {
        return description.to_owned();
    }
    let node = resolved_schema(root, node);
    if let Some(description) = node
        .get("description")
        .and_then(Value::as_str)
        .filter(|description| !description.trim().is_empty())
    {
        return description.to_owned();
    }
    for key in ["allOf", "anyOf", "oneOf"] {
        if let Some(description) = node
            .get(key)
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .map(|value| schema_description(root, value))
            .find(|description| !description.is_empty())
        {
            return description;
        }
    }
    String::new()
}

fn schema_type(root: &Value, node: &Value) -> String {
    let reference_name = node
        .get("$ref")
        .and_then(Value::as_str)
        .and_then(|reference| reference.rsplit('/').next());
    let node = resolved_schema(root, node);
    if let Some(reference_name) = reference_name
        && (node.get("enum").is_some() || node.get("properties").is_some())
    {
        return reference_name.to_owned();
    }
    if let Some(value) = node.get("type") {
        return match value {
            Value::String(value) if value == "array" => {
                if let Some(items) = node.get("prefixItems").and_then(Value::as_array) {
                    return format!(
                        "tuple<{}>",
                        items
                            .iter()
                            .map(|item| schema_type(root, item))
                            .collect::<Vec<_>>()
                            .join(", ")
                    );
                }
                node.get("items").map_or_else(
                    || "array".to_owned(),
                    |items| format!("array<{}>", schema_type(root, items)),
                )
            }
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
    if let Some(items) = node.get("prefixItems").and_then(Value::as_array) {
        for (index, item) in items.iter().enumerate() {
            values.extend(
                schema_values(root, item)
                    .into_iter()
                    .map(|value| format!("item {}: {value}", index + 1)),
            );
        }
    } else if let Some(items) = node.get("items") {
        values.extend(schema_values(root, items));
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

fn array_items_schema<'a>(root: &'a Value, node: &'a Value) -> Option<&'a Value> {
    let node = resolved_schema(root, node);
    node.get("items").or_else(|| {
        ["anyOf", "oneOf", "allOf"].into_iter().find_map(|key| {
            node.get(key)
                .and_then(Value::as_array)
                .into_iter()
                .flatten()
                .find_map(|candidate| array_items_schema(root, candidate))
        })
    })
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
    use alphal00p_docs_schema::generated::{
        CliAliasKind, CliArgument, CliArgumentAction, CliCommand, GammaLoopReference,
        SettingReference,
    };
    use clap::{CommandFactory, error::ErrorKind};
    use gammaloop_api::OneShot;

    use super::export;

    fn command<'a>(reference: &'a GammaLoopReference, path: &str) -> &'a CliCommand {
        reference
            .commands
            .iter()
            .find(|command| command.path == path)
            .unwrap_or_else(|| panic!("missing exported command {path}"))
    }

    fn argument<'a>(
        reference: &'a GammaLoopReference,
        command_path: &str,
        id: &str,
    ) -> &'a CliArgument {
        command(reference, command_path)
            .arguments
            .iter()
            .find(|argument| argument.id == id)
            .unwrap_or_else(|| panic!("missing exported argument {command_path}:{id}"))
    }

    fn setting<'a>(reference: &'a GammaLoopReference, path: &str) -> &'a SettingReference {
        reference
            .settings
            .iter()
            .find(|setting| setting.path == path)
            .unwrap_or_else(|| panic!("missing exported setting {path}"))
    }

    fn assert_parses(arguments: &[&str]) {
        OneShot::command()
            .try_get_matches_from(arguments)
            .unwrap_or_else(|error| panic!("{arguments:?} did not parse: {error}"));
    }

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

        let execution_mode = reference
            .settings
            .iter()
            .find(|setting| setting.path == "cli.global.generation.evaluator.spenso_execution_mode")
            .expect("Spenso execution setting is exported");
        assert_eq!(
            execution_mode.value_type,
            "tuple<ExecutionMode, ContractionMode>"
        );
        assert_eq!(
            execution_mode.default.as_ref(),
            Some(&serde_json::json!(["Sequential", "MinResultRank"]))
        );
        assert!(
            execution_mode
                .possible_values
                .iter()
                .any(|value| value == "item 1: Parallel")
        );
        assert!(
            execution_mode
                .possible_values
                .iter()
                .any(|value| value == "item 2: MinResultRank")
        );

        let stability_precision = reference
            .settings
            .iter()
            .find(|setting| setting.path == "runtime.stability.levels[].precision")
            .expect("array item settings are exported");
        assert_eq!(stability_precision.value_type, "Precision");
        assert!(
            stability_precision
                .possible_values
                .iter()
                .any(|value| value == "Arb")
        );

        let rotation_type = setting(&reference, "runtime.stability.rotation_axis[].type");
        assert_eq!(
            rotation_type.possible_values,
            ["euler_angles", "none", "x", "y", "z"]
        );
        assert!(rotation_type.required);
        for (angle, description) in [
            ("alpha", "First Euler angle, in radians."),
            ("beta", "Second Euler angle, in radians."),
            ("gamma", "Third Euler angle, in radians."),
        ] {
            let angle = setting(
                &reference,
                &format!("runtime.stability.rotation_axis[].{angle}"),
            );
            assert_eq!(angle.value_type, "number");
            assert_eq!(angle.description, description);
            assert!(!angle.required);
        }

        let range_type = setting(
            &reference,
            "runtime.subtraction.integrated_ct_settings.range.type",
        );
        assert_eq!(range_type.possible_values, ["compact", "infinite"]);
        assert!(range_type.required);
        let damping_function = setting(
            &reference,
            "runtime.subtraction.integrated_ct_settings.range.h_function_settings.function",
        );
        assert_eq!(
            damping_function.description,
            "Functional form of the counterterm damping profile."
        );
    }

    #[test]
    fn clean_state_is_exported_as_a_valueless_flag() {
        let reference = export().unwrap();
        let root = command(&reference, "gammaLoop");
        let clean_state = argument(&reference, "gammaLoop", "clean_state");

        assert_eq!(
            root.usage,
            "Usage: gammaLoop [OPTIONS] [BOOT_COMMANDS_PATH] [COMMAND]"
        );
        assert_eq!(clean_state.long.as_deref(), Some("clean-state"));
        assert_eq!(clean_state.action, CliArgumentAction::SetTrue);
        assert_eq!(clean_state.arity.min, 0);
        assert_eq!(clean_state.arity.max, Some(0));
        assert!(!clean_state.takes_values);
        assert!(!clean_state.value_required);
        assert!(clean_state.value_names.is_empty());
        assert_parses(&["gammaLoop", "--clean-state"]);
        assert!(
            OneShot::command()
                .try_get_matches_from(["gammaLoop", "--clean-state=true"])
                .is_err()
        );
    }

    #[test]
    fn value_options_and_positionals_match_the_compiled_parser() {
        let reference = export().unwrap();
        let state_folder = argument(&reference, "gammaLoop", "state_folder");
        let boot_card = argument(&reference, "gammaLoop", "boot_commands_path");
        let with_uv = argument(&reference, "gammaLoop save dot", "with_uv");
        let tokens = argument(&reference, "gammaLoop generate amp", "tokens");

        assert_eq!(state_folder.action, CliArgumentAction::Set);
        assert_eq!(
            (state_folder.arity.min, state_folder.arity.max),
            (1, Some(1))
        );
        assert!(state_folder.takes_values);
        assert!(state_folder.value_required);
        assert!(!state_folder.positional);
        assert!(state_folder.default_missing_values.is_empty());
        assert_parses(&["gammaLoop", "--state-folder", "/tmp/docs-parity-state"]);

        assert!(boot_card.positional);
        assert_eq!(boot_card.index, Some(1));
        assert_eq!((boot_card.arity.min, boot_card.arity.max), (1, Some(1)));
        assert_parses(&["gammaLoop", "run-card.toml"]);

        assert_eq!(with_uv.action, CliArgumentAction::Set);
        assert_eq!((with_uv.arity.min, with_uv.arity.max), (0, Some(1)));
        assert!(with_uv.takes_values);
        assert!(!with_uv.value_required);
        assert_eq!(with_uv.default_missing_values, ["true"]);
        assert_parses(&["gammaLoop", "save", "dot", "--with-uv"]);
        assert_parses(&["gammaLoop", "save", "dot", "--with-uv=false"]);

        assert_eq!(tokens.action, CliArgumentAction::Append);
        assert_eq!((tokens.arity.min, tokens.arity.max), (1, None));
        assert!(tokens.positional);
        assert_parses(&["gammaLoop", "generate", "amp", "token"]);
    }

    #[test]
    fn command_and_argument_aliases_match_nested_parser_routes() {
        let reference = export().unwrap();
        let default_runtime = command(&reference, "gammaLoop display settings default-runtime");
        let amp = command(&reference, "gammaLoop generate amp");
        let clear = argument(
            &reference,
            "gammaLoop generate amp",
            "clear_existing_processes",
        );

        assert!(default_runtime.aliases.iter().any(|alias| {
            alias.name == "defaults" && alias.kind == CliAliasKind::Name && !alias.visible
        }));
        assert!(amp.aliases.iter().any(|alias| {
            alias.name == "amplitude" && alias.kind == CliAliasKind::Name && !alias.visible
        }));
        assert!(clear.aliases.iter().any(|alias| {
            alias.name == "clear" && alias.kind == CliAliasKind::LongFlag && !alias.visible
        }));

        assert_parses(&["gammaLoop", "display", "settings", "defaults"]);
        assert_parses(&["gammaLoop", "generate", "amplitude", "token"]);
        assert_parses(&["gammaLoop", "generate", "amplitude", "token", "--clear"]);
    }

    #[test]
    fn sample_evaluation_help_uses_the_runtime_coordinate_and_precision_contract() {
        let reference = export().unwrap();

        for command_path in ["gammaLoop inspect", "gammaLoop approach"] {
            let point = argument(&reference, command_path, "point");
            let use_arb_prec = argument(&reference, command_path, "use_arb_prec");
            let momentum_space = argument(&reference, command_path, "momentum_space");

            assert!(point.help.contains("(px, py, pz)"), "{}", point.help);
            assert!(
                use_arb_prec.help.contains("arbitrary-precision")
                    && use_arb_prec.help.contains("Arb"),
                "{}",
                use_arb_prec.help
            );
            assert!(
                momentum_space.help.contains("loop-momentum triplets"),
                "{}",
                momentum_space.help
            );
        }
    }

    #[test]
    fn conflicts_and_global_inheritance_match_the_compiled_parser() {
        let reference = export().unwrap();
        let python = argument(&reference, "gammaLoop save standalone", "python");
        let no_save_state = argument(&reference, "gammaLoop", "no_save_state");
        let keep_sources = argument(&reference, "gammaLoop generate", "keep_sources");
        let inherited_keep_sources = argument(&reference, "gammaLoop generate amp", "keep_sources");
        let orientation_id = argument(&reference, "gammaLoop inspect", "orientation_id");

        assert!(python.conflicts_with.iter().any(|id| id == "rust"));
        assert!(
            no_save_state
                .conflicts_with
                .iter()
                .any(|id| id == "override_state")
        );
        assert!(keep_sources.global);
        assert!(!keep_sources.inherited);
        assert!(inherited_keep_sources.global);
        assert!(inherited_keep_sources.inherited);
        assert_eq!(orientation_id.requires, ["graph_id"]);

        assert_parses(&["gammaLoop", "save", "standalone", "--python"]);
        let conflict = OneShot::command()
            .try_get_matches_from(["gammaLoop", "save", "standalone", "--python", "--rust"])
            .unwrap_err();
        assert_eq!(conflict.kind(), ErrorKind::ArgumentConflict);
        let group_conflict = OneShot::command()
            .try_get_matches_from(["gammaLoop", "--no-save-state", "--override-state"])
            .unwrap_err();
        assert_eq!(group_conflict.kind(), ErrorKind::ArgumentConflict);
        let missing_graph = OneShot::command()
            .try_get_matches_from(["gammaLoop", "inspect", "--orientation-id", "0"])
            .unwrap_err();
        assert_eq!(missing_graph.kind(), ErrorKind::MissingRequiredArgument);
        assert_parses(&[
            "gammaLoop",
            "inspect",
            "--graph-id",
            "1",
            "--orientation-id",
            "0",
        ]);
        assert_parses(&["gammaLoop", "generate", "amp", "token", "--keep-sources"]);
    }

    #[test]
    fn generated_help_tree_is_distinct_from_public_commands() {
        let reference = export().unwrap();

        assert!(!command(&reference, "gammaLoop generate amp").generated_help);
        assert!(command(&reference, "gammaLoop help").generated_help);
        assert!(command(&reference, "gammaLoop help generate amp").generated_help);
        assert!(command(&reference, "gammaLoop generate help amp").generated_help);
    }
}
