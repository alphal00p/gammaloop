//! Python monitor keywords shared by the evaluation entry points.
use std::{path::PathBuf, time::Duration};

use pyo3::PyResult;

use crate::{RustRedMonitorDisplay, RustRedMonitorOptions};

/// Validate the per-call `monitor`, `progress_snapshot` and
/// `progress_interval_seconds` keywords.
pub(super) fn options(
    monitor: &str,
    progress_snapshot: Option<PathBuf>,
    interval: f64,
) -> PyResult<RustRedMonitorOptions> {
    let display = match monitor {
        "off" => RustRedMonitorDisplay::Off,
        "tui" => RustRedMonitorDisplay::Tui,
        "json" => RustRedMonitorDisplay::Json,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "monitor must be 'off', 'tui' or 'json'",
            ));
        }
    };
    let interval = Duration::try_from_secs_f64(interval).map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "progress_interval_seconds must be finite and nonnegative",
        )
    })?;
    let options = RustRedMonitorOptions {
        display,
        progress_snapshot,
        interval,
        ..Default::default()
    };
    options
        .validate()
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
    Ok(options)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monitor_keywords_are_validated() {
        assert!(!options("off", None, 1.0).unwrap().enabled());
        assert!(
            options("off", Some("latest.json".into()), 1.0)
                .unwrap()
                .enabled()
        );
        assert!(options("json", None, 1.0).is_err());
        assert!(options("verbose", None, 1.0).is_err());
        assert!(options("off", None, f64::NAN).is_err());
    }
}
