//! Opt-in wall timing of the reduction-math boundaries on the calling thread.
//!
//! This measures the union of nested `MathSpan`s, not the enclosing operation.
//! Owner preparation is subtracted only when its cumulative native counters were
//! explicitly sampled. It is not inferred from observer events or UI timestamps.

use std::cell::{Cell, RefCell};
use std::rc::Rc;
use std::time::{Duration, Instant};

/// The operation's unchanged result and an independently qualified timing report.
#[derive(Debug)]
pub struct ReductionTimingCapture<T, E> {
    pub result: Result<T, E>,
    pub timing: ReductionMathTiming,
}

/// Main-thread wall intervals at the reduction-math boundary.
///
/// `reduction_math` is available only for complete, successful captures with
/// consistent counters. Neither this report nor a telemetry fault changes the
/// operation's result. Native owner wall counters may be added only when their
/// intervals do not overlap; concurrent preparation cannot qualify a net time.
/// Saved-application instrumentation therefore requires one preparation worker
/// for a qualified net duration and preserves the gross value otherwise.
#[derive(Clone, Debug)]
pub struct ReductionMathTiming {
    pub scope_id: &'static str,
    pub complete: bool,
    pub gross: Duration,
    pub measured_owner_preparation: Duration,
    pub reduction_math: Option<Duration>,
    pub interval_count: usize,
    pub incompleteness_reason: Option<&'static str>,
}

/// Capture reduction timing without changing the operation's error or retrying it.
///
/// Captures are local to the calling thread. Nested captures are independent and
/// restore their predecessor even when the operation unwinds. A surrounding math
/// span still includes its complete wall interval, including a nested capture.
/// If that nested capture measures owner preparation, the surrounding capture's
/// net time is unqualified: its caller may have sampled the same native counters,
/// so preparation deltas cannot safely be propagated between captures.
pub fn capture_reduction_timing<T, E>(
    operation: impl FnOnce() -> Result<T, E>,
) -> ReductionTimingCapture<T, E> {
    let state = Rc::new(CaptureState::default());
    let scope = CaptureScope::enter(Rc::clone(&state));
    let result = operation();
    drop(scope);
    let timing = state.report(result.is_ok());
    ReductionTimingCapture { result, timing }
}

thread_local! {
    static CURRENT_CAPTURE: RefCell<Option<Rc<CaptureState>>> = const { RefCell::new(None) };
}

#[derive(Clone, Copy, Default)]
struct Measurements {
    gross: Duration,
    preparation: Duration,
    intervals: usize,
    active_spans: usize,
    started: Option<Instant>,
    nested_preparation: bool,
    fault: Option<&'static str>,
}

#[derive(Default)]
struct CaptureState {
    measurements: Cell<Measurements>,
}

impl CaptureState {
    fn invalidate(&self, reason: &'static str) {
        let mut measured = self.measurements.get();
        measured.fault.get_or_insert(reason);
        self.measurements.set(measured);
    }

    fn report(&self, successful: bool) -> ReductionMathTiming {
        let measured = self.measurements.get();
        let net = measured.gross.checked_sub(measured.preparation);
        let reason = measured.fault.or_else(|| {
            if !successful {
                Some("operation returned an error")
            } else if measured.active_spans != 0 || measured.started.is_some() {
                Some("timing spans remained open")
            } else if net.is_none() {
                Some("owner preparation exceeds gross reduction time")
            } else {
                None
            }
        });
        ReductionMathTiming {
            scope_id: "vakint.rustred.reduction-math.v1",
            complete: reason.is_none(),
            gross: measured.gross,
            measured_owner_preparation: measured.preparation,
            reduction_math: if reason.is_none() { net } else { None },
            interval_count: measured.intervals,
            incompleteness_reason: reason,
        }
    }
}

struct CaptureScope {
    state: Rc<CaptureState>,
    previous: Option<Rc<CaptureState>>,
    overlaps_previous: bool,
    installed: bool,
}

impl CaptureScope {
    fn enter(state: Rc<CaptureState>) -> Self {
        let installed = CURRENT_CAPTURE.try_with(|slot| {
            slot.try_borrow_mut()
                .map(|mut current| current.replace(Rc::clone(&state)))
        });
        match installed {
            Ok(Ok(previous)) => {
                let overlaps_previous = previous
                    .as_ref()
                    .is_some_and(|state| state.measurements.get().active_spans != 0);
                Self {
                    state,
                    previous,
                    overlaps_previous,
                    installed: true,
                }
            }
            _ => {
                state.invalidate("thread-local timing capture unavailable");
                Self {
                    state,
                    previous: None,
                    overlaps_previous: false,
                    installed: false,
                }
            }
        }
    }
}

impl Drop for CaptureScope {
    fn drop(&mut self) {
        if self.installed {
            let nested = self.state.measurements.get();
            if (nested.preparation != Duration::ZERO || nested.nested_preparation)
                && let Some(previous) = &self.previous
            {
                let mut measured = previous.measurements.get();
                // Retain evidence through an inactive intermediate capture: a
                // more distant enclosing span may still overlap preparation.
                measured.nested_preparation = true;
                if self.overlaps_previous {
                    measured.fault.get_or_insert(
                        "nested owner preparation was not attributed to this capture",
                    );
                }
                previous.measurements.set(measured);
            }
            let previous = self.previous.take();
            let restored = CURRENT_CAPTURE.try_with(|slot| {
                if let Ok(mut current) = slot.try_borrow_mut() {
                    *current = previous;
                    true
                } else {
                    false
                }
            });
            if !matches!(restored, Ok(true)) {
                self.state
                    .invalidate("thread-local timing capture could not be restored");
            }
        }
    }
}

/// Counts only the outermost active interval in its capture.
///
/// With capture disabled this holds no state, allocates nothing and reads no
/// clock. Declare the guard before locals whose destruction belongs to the
/// measured boundary, so their natural drop remains inside the interval.
pub(crate) struct MathSpan {
    state: Option<Rc<CaptureState>>,
    entered: bool,
}

impl MathSpan {
    pub(crate) fn new() -> Self {
        let state = CURRENT_CAPTURE
            .try_with(|slot| slot.try_borrow().ok().and_then(|current| current.clone()))
            .ok()
            .flatten();
        let mut span = Self {
            state,
            entered: false,
        };
        if let Some(state) = &span.state {
            let mut measured = state.measurements.get();
            if let Some(depth) = measured.active_spans.checked_add(1) {
                if measured.active_spans == 0 {
                    measured.started = Some(Instant::now());
                }
                measured.active_spans = depth;
                span.entered = true;
                state.measurements.set(measured);
            } else {
                state.invalidate("timing span count overflow");
            }
        }
        span
    }

    pub(crate) fn is_active(&self) -> bool {
        self.state.is_some()
    }

    pub(crate) fn invalidate(&self, reason: &'static str) {
        if let Some(state) = &self.state {
            state.invalidate(reason);
        }
    }

    /// Add a measured cumulative wall-counter delta, never an estimated phase.
    /// Call sites must avoid recording overlapping before/after intervals.
    pub(crate) fn owner_preparation(&self, before: Duration, after: Duration) {
        if let Some(state) = &self.state {
            let mut measured = state.measurements.get();
            if let Some(total) = after
                .checked_sub(before)
                .and_then(|delta| measured.preparation.checked_add(delta))
            {
                measured.preparation = total;
                state.measurements.set(measured);
            } else {
                state.invalidate("owner preparation counter decreased or overflowed");
            }
        }
    }
}

impl Drop for MathSpan {
    fn drop(&mut self) {
        if !self.entered {
            return;
        }
        let Some(state) = &self.state else {
            return;
        };
        let mut measured = state.measurements.get();
        let Some(depth) = measured.active_spans.checked_sub(1) else {
            state.invalidate("timing span count underflow");
            return;
        };
        measured.active_spans = depth;
        if depth == 0 {
            let interval = measured
                .started
                .take()
                .and_then(|start| Instant::now().checked_duration_since(start));
            if let Some((gross, count)) = interval.and_then(|elapsed| {
                Some((
                    measured.gross.checked_add(elapsed)?,
                    measured.intervals.checked_add(1)?,
                ))
            }) {
                measured.gross = gross;
                measured.intervals = count;
            } else {
                measured
                    .fault
                    .get_or_insert("reduction timing interval missing or overflowed");
            }
        }
        state.measurements.set(measured);
    }
}

#[cfg(test)]
mod tests {
    use super::{MathSpan, capture_reduction_timing};
    use std::panic::{AssertUnwindSafe, catch_unwind};
    use std::time::Duration;

    #[test]
    fn disabled_span_has_no_state_and_does_not_affect_results() {
        let span = MathSpan::new();
        assert!(!span.is_active());
        span.owner_preparation(Duration::MAX, Duration::ZERO);
        span.invalidate("ignored without capture");
        drop(span);
        let captured = capture_reduction_timing(|| Ok::<_, ()>(42));
        assert_eq!(captured.result, Ok(42));
        assert!(captured.timing.complete);
        assert_eq!(captured.timing.interval_count, 0);
        assert_eq!(captured.timing.reduction_math, Some(Duration::ZERO));
    }

    #[test]
    fn nested_spans_count_only_the_outer_interval_and_include_later_drops() {
        struct DelayOnDrop;
        impl Drop for DelayOnDrop {
            fn drop(&mut self) {
                std::thread::sleep(Duration::from_millis(2));
            }
        }
        let captured = capture_reduction_timing(|| {
            let _outer = MathSpan::new();
            let _later_local = DelayOnDrop;
            let _inner = MathSpan::new();
            Ok::<_, ()>(7)
        });
        assert_eq!(captured.result, Ok(7));
        assert!(captured.timing.complete);
        assert_eq!(captured.timing.interval_count, 1);
        assert!(captured.timing.gross >= Duration::from_millis(2));
        assert_eq!(captured.timing.reduction_math, Some(captured.timing.gross));
    }

    #[test]
    fn nested_capture_is_independent_and_restores_previous_capture() {
        let outer = capture_reduction_timing(|| {
            drop(MathSpan::new());
            let inner = capture_reduction_timing(|| {
                let _span = MathSpan::new();
                Ok::<_, ()>(())
            });
            assert!(inner.timing.complete);
            assert_eq!(inner.timing.interval_count, 1);
            drop(MathSpan::new());
            Ok::<_, ()>(())
        });
        assert!(outer.timing.complete);
        assert_eq!(outer.timing.interval_count, 2);
        assert!(!MathSpan::new().is_active());
    }

    #[test]
    fn active_outer_span_survives_a_nested_capture() {
        let outer = capture_reduction_timing(|| {
            let _outer = MathSpan::new();
            let inner = capture_reduction_timing(|| {
                let _inner = MathSpan::new();
                Ok::<_, ()>(())
            });
            assert!(inner.timing.complete);
            assert_eq!(inner.timing.interval_count, 1);
            Ok::<_, ()>(())
        });
        assert!(outer.timing.complete);
        assert_eq!(outer.timing.interval_count, 1);
    }

    #[test]
    fn nested_owner_preparation_disqualifies_overlapping_outer_capture() {
        let outer = capture_reduction_timing(|| {
            let _outer = MathSpan::new();
            let inner = capture_reduction_timing(|| {
                let inner = MathSpan::new();
                inner.owner_preparation(Duration::ZERO, Duration::from_nanos(1));
                std::thread::sleep(Duration::from_millis(2));
                Ok::<_, ()>(())
            });
            assert!(inner.timing.complete);
            assert_eq!(
                inner.timing.measured_owner_preparation,
                Duration::from_nanos(1)
            );
            Ok::<_, ()>(23)
        });
        assert_eq!(outer.result, Ok(23));
        assert!(!outer.timing.complete);
        assert_eq!(outer.timing.reduction_math, None);
        assert_eq!(outer.timing.measured_owner_preparation, Duration::ZERO);
        assert_eq!(
            outer.timing.incompleteness_reason,
            Some("nested owner preparation was not attributed to this capture")
        );
    }

    #[test]
    fn nested_preparation_outside_outer_spans_does_not_disqualify_outer() {
        let outer = capture_reduction_timing(|| {
            drop(MathSpan::new());
            let inner = capture_reduction_timing(|| {
                let inner = MathSpan::new();
                inner.owner_preparation(Duration::ZERO, Duration::from_nanos(1));
                std::thread::sleep(Duration::from_millis(2));
                Ok::<_, ()>(())
            });
            assert!(inner.timing.complete);
            drop(MathSpan::new());
            Ok::<_, ()>(())
        });
        assert!(outer.timing.complete);
        assert_eq!(outer.timing.interval_count, 2);
        assert_eq!(outer.timing.measured_owner_preparation, Duration::ZERO);
    }

    #[test]
    fn nested_preparation_disqualification_propagates_through_deeper_captures() {
        let outer = capture_reduction_timing(|| {
            let _outer = MathSpan::new();
            let middle = capture_reduction_timing(|| {
                let middle = MathSpan::new();
                middle.invalidate("first diagnostic");
                let inner = capture_reduction_timing(|| {
                    let inner = MathSpan::new();
                    inner.owner_preparation(Duration::ZERO, Duration::from_nanos(1));
                    std::thread::sleep(Duration::from_millis(2));
                    Ok::<_, ()>(())
                });
                assert!(inner.timing.complete);
                Ok::<_, ()>(())
            });
            assert_eq!(
                middle.timing.incompleteness_reason,
                Some("first diagnostic")
            );
            assert_eq!(middle.timing.measured_owner_preparation, Duration::ZERO);
            Ok::<_, ()>(())
        });
        assert_eq!(outer.result, Ok(()));
        assert!(!outer.timing.complete);
        assert_eq!(outer.timing.reduction_math, None);
        assert_eq!(
            outer.timing.incompleteness_reason,
            Some("nested owner preparation was not attributed to this capture")
        );
    }

    #[test]
    fn nested_preparation_crosses_an_inactive_intermediate_capture() {
        let outer = capture_reduction_timing(|| {
            let _outer = MathSpan::new();
            let middle = capture_reduction_timing(|| {
                let inner = capture_reduction_timing(|| {
                    let inner = MathSpan::new();
                    inner.owner_preparation(Duration::ZERO, Duration::from_nanos(1));
                    std::thread::sleep(Duration::from_millis(2));
                    Ok::<_, ()>(())
                });
                assert!(inner.timing.complete);
                Ok::<_, ()>(())
            });
            assert!(middle.timing.complete);
            assert_eq!(middle.timing.gross, Duration::ZERO);
            Ok::<_, ()>(())
        });
        assert_eq!(outer.result, Ok(()));
        assert!(!outer.timing.complete);
        assert_eq!(outer.timing.reduction_math, None);
        assert_eq!(outer.timing.measured_owner_preparation, Duration::ZERO);
    }

    #[test]
    fn owner_preparation_uses_measured_deltas_and_subtracts_exactly() {
        let captured = capture_reduction_timing(|| {
            let outer = MathSpan::new();
            outer.owner_preparation(Duration::from_secs(1), Duration::from_secs(1));
            {
                let inner = MathSpan::new();
                inner.owner_preparation(
                    Duration::from_secs(1),
                    Duration::from_secs(1) + Duration::from_micros(1),
                );
            }
            std::thread::sleep(Duration::from_millis(2));
            Ok::<_, ()>(())
        });
        assert!(captured.timing.complete);
        assert_eq!(captured.timing.interval_count, 1);
        assert_eq!(
            captured.timing.measured_owner_preparation,
            Duration::from_micros(1)
        );
        assert_eq!(
            captured.timing.reduction_math,
            captured.timing.gross.checked_sub(Duration::from_micros(1))
        );
    }

    #[test]
    fn capture_does_not_activate_other_threads() {
        let captured = capture_reduction_timing(|| {
            let _span = MathSpan::new();
            std::thread::spawn(|| assert!(!MathSpan::new().is_active()))
                .join()
                .unwrap();
            Ok::<_, ()>(())
        });
        assert!(captured.timing.complete);
        assert_eq!(captured.timing.interval_count, 1);
    }

    #[test]
    fn operation_error_is_unchanged_and_disqualifies_net_time() {
        let captured = capture_reduction_timing(|| {
            let _span = MathSpan::new();
            Err::<(), _>("original error")
        });
        assert_eq!(captured.result, Err("original error"));
        assert!(!captured.timing.complete);
        assert_eq!(captured.timing.interval_count, 1);
        assert_eq!(captured.timing.reduction_math, None);
        assert_eq!(
            captured.timing.incompleteness_reason,
            Some("operation returned an error")
        );
    }

    #[test]
    fn unwinding_nested_capture_restores_outer_capture() {
        let outer = capture_reduction_timing(|| {
            drop(MathSpan::new());
            let unwind = catch_unwind(AssertUnwindSafe(|| {
                capture_reduction_timing(|| -> Result<(), ()> {
                    let _span = MathSpan::new();
                    panic!("original panic");
                })
            }));
            assert!(unwind.is_err());
            drop(MathSpan::new());
            Ok::<_, ()>(())
        });
        assert!(outer.timing.complete);
        assert_eq!(outer.timing.interval_count, 2);
        assert!(!MathSpan::new().is_active());
    }

    #[test]
    fn returned_open_span_disqualifies_capture_without_changing_result() {
        let captured = capture_reduction_timing(|| Ok::<_, ()>(MathSpan::new()));
        assert!(captured.result.is_ok());
        assert!(!captured.timing.complete);
        assert_eq!(captured.timing.reduction_math, None);
        assert_eq!(
            captured.timing.incompleteness_reason,
            Some("timing spans remained open")
        );
        drop(captured.result);
        assert!(!MathSpan::new().is_active());
    }

    #[test]
    fn invalid_preparation_counter_never_changes_the_operation_result() {
        let captured = capture_reduction_timing(|| {
            let span = MathSpan::new();
            span.owner_preparation(Duration::from_secs(1), Duration::ZERO);
            Ok::<_, ()>("original value")
        });
        assert_eq!(captured.result, Ok("original value"));
        assert!(!captured.timing.complete);
        assert_eq!(captured.timing.reduction_math, None);
    }

    #[test]
    fn preparation_overflow_is_reported_without_panicking() {
        let captured = capture_reduction_timing(|| {
            let span = MathSpan::new();
            let state = span.state.as_ref().unwrap();
            let mut measured = state.measurements.get();
            measured.preparation = Duration::MAX;
            state.measurements.set(measured);
            span.owner_preparation(Duration::ZERO, Duration::from_nanos(1));
            Ok::<_, ()>(())
        });
        assert_eq!(captured.result, Ok(()));
        assert!(!captured.timing.complete);
        assert_eq!(captured.timing.reduction_math, None);
    }

    #[test]
    fn preparation_larger_than_gross_is_not_clamped_to_zero() {
        let captured = capture_reduction_timing(|| {
            let span = MathSpan::new();
            span.owner_preparation(Duration::ZERO, Duration::MAX);
            Ok::<_, ()>(23)
        });
        assert_eq!(captured.result, Ok(23));
        assert!(!captured.timing.complete);
        assert_eq!(captured.timing.reduction_math, None);
        assert_eq!(
            captured.timing.incompleteness_reason,
            Some("owner preparation exceeds gross reduction time")
        );
    }

    #[test]
    fn explicit_invalidation_keeps_first_reason_and_original_result() {
        let captured = capture_reduction_timing(|| {
            let span = MathSpan::new();
            span.invalidate("preparation was concurrent");
            span.invalidate("later diagnostic");
            Ok::<_, ()>(31)
        });
        assert_eq!(captured.result, Ok(31));
        assert!(!captured.timing.complete);
        assert_eq!(captured.timing.reduction_math, None);
        assert_eq!(
            captured.timing.incompleteness_reason,
            Some("preparation was concurrent")
        );
    }
}
