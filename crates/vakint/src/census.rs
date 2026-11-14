//! Diagnostics: with `VAKINT_LOWERING_CENSUS=<file>`, the stages that turn an
//! integrand into reduction input append one JSON line per call to `<file>`
//! (term counts, byte sizes, distinct structures). Off by default; the
//! checks below cost one relaxed load when off.

use std::sync::LazyLock;

use symbolica::atom::AtomView;

static PATH: LazyLock<Option<std::path::PathBuf>> =
    LazyLock::new(|| std::env::var_os("VAKINT_LOWERING_CENSUS").map(Into::into));

/// Whether the census is on.
pub fn enabled() -> bool {
    PATH.is_some()
}

/// The number of summands of an expression (1 for a non-sum, 0 for zero).
pub fn summands(view: AtomView<'_>) -> usize {
    match view {
        AtomView::Add(sum) => sum.get_nargs(),
        _ if view.is_zero() => 0,
        _ => 1,
    }
}

/// Append `{"stage": stage, ...fields}` to the census file.
pub fn record(stage: &str, fields: serde_json::Value) {
    use std::io::Write;
    static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let Some(path) = PATH.as_ref() else {
        return;
    };
    let mut fields = fields;
    if let serde_json::Value::Object(map) = &mut fields {
        map.insert("stage".into(), stage.into());
        map.insert(
            "thread".into(),
            format!("{:?}", std::thread::current().id()).into(),
        );
    }
    let _guard = LOCK.lock();
    if let Ok(mut file) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
    {
        let _ = writeln!(file, "{fields}");
    }
}
