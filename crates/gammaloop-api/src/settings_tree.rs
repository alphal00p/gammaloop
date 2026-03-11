use std::sync::{atomic::Ordering, Mutex, MutexGuard};

use color_eyre::Result;
use eyre::{eyre, Context};
use gammalooprs::utils::serde_utils::SHOWDEFAULTS;
use serde::Serialize;
use serde_json::Value as JsonValue;

pub(crate) fn serialize_settings_with_defaults<T: Serialize>(
    settings: &T,
    context: &str,
) -> Result<JsonValue> {
    let _show_defaults_guard = ShowDefaultsGuard::new(true);
    serde_json::to_value(settings).context(format!("While serializing {context}"))
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

struct ShowDefaultsGuard {
    _lock: MutexGuard<'static, ()>,
    previous: bool,
}

impl ShowDefaultsGuard {
    fn new(show_defaults: bool) -> Self {
        static SHOWDEFAULTS_MUTEX: Mutex<()> = Mutex::new(());
        let lock = SHOWDEFAULTS_MUTEX
            .lock()
            .expect("SHOWDEFAULTS serialization mutex must not be poisoned");
        let previous = SHOWDEFAULTS.swap(show_defaults, Ordering::Relaxed);
        Self {
            _lock: lock,
            previous,
        }
    }
}

impl Drop for ShowDefaultsGuard {
    fn drop(&mut self) {
        SHOWDEFAULTS.store(self.previous, Ordering::Relaxed);
    }
}
