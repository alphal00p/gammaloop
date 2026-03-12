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
