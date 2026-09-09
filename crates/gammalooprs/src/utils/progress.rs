use std::{fmt, time::Duration};

use indicatif::{FormattedDuration, ProgressState, ProgressStyle};

const LONG_RUNNING_PROGRESS_TEMPLATE: &str = "[{elapsed_precise} | ETA: {eta_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({percent}%) {msg}";

fn write_eta(eta: Duration, writer: &mut dyn fmt::Write) -> fmt::Result {
    // Indicatif saturates a floating-point ETA overflow to `u64::MAX` seconds.
    // Rendering that sentinel literally produces a meaningless 213503982334601-day ETA.
    if eta.as_secs() == u64::MAX {
        writer.write_str("N/A")
    } else {
        write!(writer, "{}", FormattedDuration(eta))
    }
}

fn estimate_average_eta(
    elapsed: Duration,
    completed_steps: u64,
    total_steps: Option<u64>,
) -> Option<Duration> {
    let total_steps = total_steps?;
    if completed_steps == 0 {
        return None;
    }
    if completed_steps >= total_steps {
        return Some(Duration::ZERO);
    }

    let remaining_steps = total_steps - completed_steps;
    let eta_nanos = elapsed
        .as_nanos()
        .checked_mul(u128::from(remaining_steps))?
        / u128::from(completed_steps);
    let eta_seconds = eta_nanos / 1_000_000_000;
    if eta_seconds > u128::from(u64::MAX) {
        return None;
    }
    Some(Duration::from_secs(eta_seconds as u64))
}

fn write_eta_after_warmup(
    eta: Option<Duration>,
    completed_steps: u64,
    minimum_completed_steps: u64,
    is_finished: bool,
    writer: &mut dyn fmt::Write,
) -> fmt::Result {
    if !is_finished && completed_steps < minimum_completed_steps {
        writer.write_str("N/A")
    } else if let Some(eta) = eta {
        write_eta(eta, writer)
    } else {
        writer.write_str("N/A")
    }
}

fn write_safe_eta(state: &ProgressState, writer: &mut dyn fmt::Write) {
    let is_finished =
        state.is_finished() || state.len().is_some_and(|length| state.pos() >= length);
    let _ = write_eta_after_warmup(
        estimate_average_eta(state.elapsed(), state.pos(), state.len()),
        state.pos(),
        1,
        is_finished,
        writer,
    );
}

/// Standard style for bounded, long-running progress bars.
///
/// The custom ETA renderer uses average completed-work throughput and keeps
/// unavailable or overflowing estimates out of terminal output.
pub fn long_running_progress_style() -> ProgressStyle {
    ProgressStyle::with_template(LONG_RUNNING_PROGRESS_TEMPLATE)
        .expect("long-running progress bar template should be valid")
        .with_key("eta_precise", write_safe_eta)
}

/// Long-running progress style that withholds ETA until enough work units have
/// completed for the rate estimate to be meaningful.
pub fn long_running_progress_style_with_eta_warmup(minimum_completed_steps: u64) -> ProgressStyle {
    ProgressStyle::with_template(LONG_RUNNING_PROGRESS_TEMPLATE)
        .expect("long-running progress bar template should be valid")
        .with_key(
            "eta_precise",
            move |state: &ProgressState, writer: &mut dyn fmt::Write| {
                let is_finished =
                    state.is_finished() || state.len().is_some_and(|length| state.pos() >= length);
                let _ = write_eta_after_warmup(
                    estimate_average_eta(state.elapsed(), state.pos(), state.len()),
                    state.pos(),
                    minimum_completed_steps,
                    is_finished,
                    writer,
                );
            },
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn saturated_eta_is_rendered_as_unavailable() {
        let mut rendered = String::new();
        write_eta(Duration::new(u64::MAX, 0), &mut rendered).unwrap();
        assert_eq!(rendered, "N/A");
    }

    #[test]
    fn finite_eta_keeps_precise_duration_format() {
        let mut rendered = String::new();
        write_eta(Duration::from_secs(3_661), &mut rendered).unwrap();
        assert_eq!(rendered, "01:01:01");
    }

    #[test]
    fn eta_is_hidden_until_warmup_completes() {
        let mut rendered = String::new();
        write_eta_after_warmup(
            Some(Duration::from_secs(10_533)),
            1,
            10,
            false,
            &mut rendered,
        )
        .unwrap();
        assert_eq!(rendered, "N/A");

        rendered.clear();
        write_eta_after_warmup(Some(Duration::from_secs(600)), 10, 10, false, &mut rendered)
            .unwrap();
        assert_eq!(rendered, "00:10:00");
    }

    #[test]
    fn average_eta_uses_completed_work_instead_of_idle_ticks() {
        assert_eq!(
            estimate_average_eta(Duration::from_secs(655), 21, Some(332)),
            Some(Duration::from_secs(9_700)),
        );
        assert_eq!(
            estimate_average_eta(Duration::from_secs(655), 0, Some(332)),
            None
        );
    }
}
