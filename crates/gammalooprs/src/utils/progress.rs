use std::{fmt, time::Duration};

use indicatif::{FormattedDuration, ProgressState, ProgressStyle};

const LONG_RUNNING_PROGRESS_TEMPLATE: &str = "[{elapsed_precise} | ETA: {eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}";

fn write_eta(eta: Duration, writer: &mut dyn fmt::Write) -> fmt::Result {
    // Indicatif saturates a floating-point ETA overflow to `u64::MAX` seconds.
    // Rendering that sentinel literally produces a meaningless 213503982334601-day ETA.
    if eta.as_secs() == u64::MAX {
        writer.write_str("--:--:--")
    } else {
        write!(writer, "{}", FormattedDuration(eta))
    }
}

fn write_safe_eta(state: &ProgressState, writer: &mut dyn fmt::Write) {
    let _ = write_eta(state.eta(), writer);
}

/// Standard style for bounded, long-running progress bars.
///
/// The custom ETA renderer keeps Indicatif's saturated-duration sentinel out of
/// terminal output while preserving its normal precise formatting.
pub fn long_running_progress_style() -> ProgressStyle {
    ProgressStyle::with_template(LONG_RUNNING_PROGRESS_TEMPLATE)
        .expect("long-running progress bar template should be valid")
        .with_key("eta_precise", write_safe_eta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn saturated_eta_is_rendered_as_unavailable() {
        let mut rendered = String::new();
        write_eta(Duration::new(u64::MAX, 0), &mut rendered).unwrap();
        assert_eq!(rendered, "--:--:--");
    }

    #[test]
    fn finite_eta_keeps_precise_duration_format() {
        let mut rendered = String::new();
        write_eta(Duration::from_secs(3_661), &mut rendered).unwrap();
        assert_eq!(rendered, "01:01:01");
    }
}
