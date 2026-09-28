use color_eyre::Result;
use eyre::{eyre, Context};
use gammalooprs::utils::serde_utils::ShowDefaultsGuard;
use schemars::{schema_for, JsonSchema};
use serde::Serialize;
use serde_json::Value as JsonValue;

pub(crate) fn serialize_settings_with_defaults<T: Serialize>(
    settings: &T,
    context: &str,
) -> Result<JsonValue> {
    let _show_defaults_guard = ShowDefaultsGuard::new(true);
    serde_json::to_value(settings).context(format!("While serializing {context}"))
}

pub(crate) fn serialize_schema<T: JsonSchema>(context: &str) -> Result<JsonValue> {
    serde_json::to_value(schema_for!(T)).context(format!("While serializing schema for {context}"))
}

pub(crate) fn value_at_path<'a>(root: &'a JsonValue, path: &str) -> Result<&'a JsonValue> {
    let mut current = root;
    for segment in path.split('.') {
        if segment.is_empty() {
            continue;
        }
        current = match current {
            JsonValue::Object(map) => map.get(segment).ok_or_else(|| {
                let mut available: Vec<_> = map.keys().cloned().collect();
                available.sort();
                eyre!(
                    "Key segment '{}' not found while resolving '{}'. Available keys: {}",
                    segment,
                    path,
                    available.join(", ")
                )
            })?,
            JsonValue::Array(items) => {
                let index = segment.parse::<usize>().map_err(|_| {
                    eyre!(
                        "Array segment '{}' is not a valid index while resolving '{}'",
                        segment,
                        path
                    )
                })?;
                items.get(index).ok_or_else(|| {
                    eyre!(
                        "Array index {} out of bounds while resolving '{}', len={}",
                        index,
                        path,
                        items.len()
                    )
                })?
            }
            _ => {
                return Err(eyre!(
                    "Cannot descend into segment '{}' while resolving '{}': current value is {}",
                    segment,
                    path,
                    json_type_name(current)
                ));
            }
        };
    }
    Ok(current)
}

pub(crate) fn json_type_name(value: &JsonValue) -> &'static str {
    match value {
        JsonValue::Null => "null",
        JsonValue::Bool(_) => "boolean",
        JsonValue::Number(number) if number.is_i64() || number.is_u64() => "integer",
        JsonValue::Number(_) => "number",
        JsonValue::String(_) => "string",
        JsonValue::Array(_) => "array",
        JsonValue::Object(_) => "object",
    }
}

pub(crate) fn schema_at_path<'a>(root: &'a JsonValue, path: &str) -> Option<&'a JsonValue> {
    let mut current = root;
    for segment in path.split('.').filter(|segment| !segment.is_empty()) {
        current = schema_child_for_segment(root, current, segment)?;
    }
    resolve_schema_refs(root, current)
}

pub(crate) fn schema_enum_values(root: &JsonValue, node: &JsonValue) -> Vec<String> {
    let mut values = Vec::new();
    collect_schema_enum_values(root, node, &mut values);
    values
}

pub(crate) fn schema_value_hint(root: &JsonValue, node: &JsonValue) -> Option<String> {
    let description = schema_description(root, node);
    let type_hint = schema_type_hint(root, node);
    let example = schema_example_values(root, node).into_iter().next();

    match (description, type_hint, example) {
        (Some(description), _, Some(example)) => Some(format!("{description} Example: {example}")),
        (Some(description), _, None) => Some(description),
        (None, Some(type_hint), Some(example)) => Some(format!("{type_hint}; e.g. {example}")),
        (None, Some(type_hint), None) => Some(type_hint),
        (None, None, Some(example)) => Some(format!("e.g. {example}")),
        (None, None, None) => None,
    }
}

pub(crate) fn schema_object_property_names(root: &JsonValue, node: &JsonValue) -> Vec<String> {
    let mut values = Vec::new();
    collect_schema_object_property_names(root, node, &mut values);
    values
}

pub(crate) fn schema_is_object_container(root: &JsonValue, node: &JsonValue) -> bool {
    let Some(node) = resolve_schema_refs(root, node) else {
        return false;
    };

    if node
        .get("properties")
        .and_then(JsonValue::as_object)
        .is_some_and(|properties| !properties.is_empty())
    {
        return true;
    }

    match node.get("type") {
        Some(JsonValue::String(type_name)) if type_name == "object" => return true,
        Some(JsonValue::Array(type_names))
            if type_names
                .iter()
                .filter_map(JsonValue::as_str)
                .any(|type_name| type_name == "object") =>
        {
            return true;
        }
        _ => {}
    }

    for keyword in ["anyOf", "oneOf", "allOf"] {
        let Some(variants) = node.get(keyword).and_then(JsonValue::as_array) else {
            continue;
        };
        if variants
            .iter()
            .any(|variant| schema_is_object_container(root, variant))
        {
            return true;
        }
    }

    false
}

pub(crate) fn schema_example_values(root: &JsonValue, node: &JsonValue) -> Vec<String> {
    let mut values = Vec::new();
    collect_schema_examples(root, node, &mut values);
    values
}

fn schema_child_for_segment<'a>(
    root: &'a JsonValue,
    node: &'a JsonValue,
    segment: &str,
) -> Option<&'a JsonValue> {
    let node = resolve_schema_refs(root, node)?;

    if let Some(properties) = node.get("properties").and_then(JsonValue::as_object) {
        if let Some(child) = properties.get(segment) {
            return Some(child);
        }
    }

    if segment.parse::<usize>().is_ok() {
        if let Some(items) = node.get("items") {
            return Some(items);
        }
    }

    for keyword in ["anyOf", "oneOf", "allOf"] {
        if let Some(variants) = node.get(keyword).and_then(JsonValue::as_array) {
            for variant in variants {
                if let Some(child) = schema_child_for_segment(root, variant, segment) {
                    return Some(child);
                }
            }
        }
    }

    None
}

fn resolve_schema_refs<'a>(root: &'a JsonValue, mut node: &'a JsonValue) -> Option<&'a JsonValue> {
    loop {
        let Some(reference) = node.get("$ref").and_then(JsonValue::as_str) else {
            return Some(node);
        };
        let pointer = reference.strip_prefix('#')?;
        node = root.pointer(pointer)?;
    }
}

fn collect_schema_enum_values(root: &JsonValue, node: &JsonValue, values: &mut Vec<String>) {
    let Some(node) = resolve_schema_refs(root, node) else {
        return;
    };

    if let Some(enum_values) = node.get("enum").and_then(JsonValue::as_array) {
        for value in enum_values {
            let Some(value) = value.as_str() else {
                continue;
            };
            if !values.iter().any(|existing| existing == value) {
                values.push(value.to_string());
            }
        }
    }

    for keyword in ["anyOf", "oneOf", "allOf"] {
        let Some(variants) = node.get(keyword).and_then(JsonValue::as_array) else {
            continue;
        };
        for variant in variants {
            collect_schema_enum_values(root, variant, values);
        }
    }
}

fn collect_schema_object_property_names(
    root: &JsonValue,
    node: &JsonValue,
    values: &mut Vec<String>,
) {
    let Some(node) = resolve_schema_refs(root, node) else {
        return;
    };

    if let Some(properties) = node.get("properties").and_then(JsonValue::as_object) {
        for key in properties.keys() {
            if !values.iter().any(|existing| existing == key) {
                values.push(key.clone());
            }
        }
    }

    for keyword in ["anyOf", "oneOf", "allOf"] {
        let Some(variants) = node.get(keyword).and_then(JsonValue::as_array) else {
            continue;
        };
        for variant in variants {
            collect_schema_object_property_names(root, variant, values);
        }
    }
}

fn schema_description(root: &JsonValue, node: &JsonValue) -> Option<String> {
    let node = resolve_schema_refs(root, node)?;
    node.get("description")
        .and_then(JsonValue::as_str)
        .map(str::trim)
        .filter(|description| !description.is_empty())
        .map(str::to_string)
}

fn schema_type_hint(root: &JsonValue, node: &JsonValue) -> Option<String> {
    let mut type_names = Vec::new();
    collect_schema_type_names(root, node, &mut type_names);
    if type_names.is_empty() {
        return None;
    }

    let mut descriptions = type_names
        .into_iter()
        .map(schema_type_description)
        .collect::<Vec<_>>();
    descriptions.dedup();

    Some(match descriptions.as_slice() {
        [] => return None,
        [only] => format!("expects {only}"),
        [first, second] => format!("expects {first} or {second}"),
        _ => {
            let last = descriptions.pop().unwrap();
            format!("expects {} or {last}", descriptions.join(", "))
        }
    })
}

fn collect_schema_type_names(root: &JsonValue, node: &JsonValue, values: &mut Vec<&'static str>) {
    let Some(node) = resolve_schema_refs(root, node) else {
        return;
    };

    match node.get("type") {
        Some(JsonValue::String(name)) => push_unique_type_name(values, name),
        Some(JsonValue::Array(type_names)) => {
            for type_name in type_names {
                if let Some(name) = type_name.as_str() {
                    push_unique_type_name(values, name);
                }
            }
        }
        _ => {}
    }

    for keyword in ["anyOf", "oneOf", "allOf"] {
        let Some(variants) = node.get(keyword).and_then(JsonValue::as_array) else {
            continue;
        };
        for variant in variants {
            collect_schema_type_names(root, variant, values);
        }
    }
}

fn push_unique_type_name(values: &mut Vec<&'static str>, type_name: &str) {
    let normalized = match type_name {
        "boolean" => Some("boolean"),
        "integer" => Some("integer"),
        "number" => Some("number"),
        "string" => Some("string"),
        "array" => Some("array"),
        "object" => Some("object"),
        "null" => Some("null"),
        _ => None,
    };
    let Some(normalized) = normalized else {
        return;
    };
    if !values.contains(&normalized) {
        values.push(normalized);
    }
}

fn schema_type_description(type_name: &'static str) -> &'static str {
    match type_name {
        "array" => "an array",
        "boolean" => "a boolean",
        "integer" => "an integer",
        "null" => "null",
        "number" => "a number",
        "object" => "an object",
        "string" => "a string",
        _ => type_name,
    }
}

fn collect_schema_examples(root: &JsonValue, node: &JsonValue, values: &mut Vec<String>) {
    let Some(node) = resolve_schema_refs(root, node) else {
        return;
    };

    if let Some(examples) = node.get("examples").and_then(JsonValue::as_array) {
        for example in examples {
            push_unique_example(values, example);
        }
    }
    if let Some(example) = node.get("example") {
        push_unique_example(values, example);
    }

    for keyword in ["anyOf", "oneOf", "allOf"] {
        let Some(variants) = node.get(keyword).and_then(JsonValue::as_array) else {
            continue;
        };
        for variant in variants {
            collect_schema_examples(root, variant, values);
        }
    }
}

fn push_unique_example(values: &mut Vec<String>, example: &JsonValue) {
    let Ok(rendered) = serde_json::to_string(example) else {
        return;
    };
    if !values.contains(&rendered) {
        values.push(rendered);
    }
}
