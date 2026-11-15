//! Opt-in, thread-local wall-clock timings for nested calculation stages.
//!
//! Enable INFO events for the `vakint::stage_timing` tracing target. Exclusive
//! time subtracts direct children's inclusive time, including recursively
//! nested stages only once. It still includes waits and uninstrumented calls;
//! overlapping worker timings must not be summed as CPU time or process wall
//! time. Logging happens after the child's clock stops, so its overhead is
//! attributed to an enclosing stage's exclusive time.

use std::{
    cell::RefCell,
    marker::PhantomData,
    rc::Rc,
    time::{Duration, Instant},
};

thread_local! {
    static FRAMES: RefCell<FrameStack> = const { RefCell::new(FrameStack::new()) };
}

/// A lexical calculation stage, measured only when its tracing target is enabled.
///
/// Guards must stay on their originating thread and should be dropped in reverse
/// construction order. Out-of-order drops produce `nesting_valid = false`.
///
/// ```compile_fail
/// let stage = vakint::stage_timing::StageSpan::new("calculation");
/// std::thread::spawn(move || drop(stage));
/// ```
pub struct StageSpan {
    active: Option<ActiveSpan>,
    _not_send: PhantomData<Rc<()>>,
}

struct ActiveSpan {
    metadata: Metadata,
    started: Instant,
}

#[derive(Clone, Copy)]
struct Metadata {
    id: u64,
    stage: &'static str,
    depth: usize,
}

struct Frame {
    metadata: Metadata,
    parent: Option<u64>,
    children: Duration,
    valid: bool,
}

struct FrameStack {
    frames: Vec<Frame>,
    next_id: u64,
}

struct Timing {
    stage: &'static str,
    path: String,
    depth: usize,
    inclusive: Duration,
    exclusive: Duration,
    nesting_valid: bool,
}

impl StageSpan {
    pub fn new(stage: &'static str) -> Self {
        Self::start(
            stage,
            tracing::enabled!(target: "vakint::stage_timing", tracing::Level::INFO),
            Instant::now,
        )
    }

    fn start(stage: &'static str, enabled: bool, now: impl FnOnce() -> Instant) -> Self {
        let active = if enabled {
            FRAMES
                .try_with(|frames| {
                    let metadata = frames.try_borrow_mut().ok()?.enter(stage)?;
                    Some(ActiveSpan {
                        metadata,
                        started: now(),
                    })
                })
                .ok()
                .flatten()
        } else {
            // Do not touch TLS, allocate, or read a clock on the disabled path.
            None
        };
        Self {
            active,
            _not_send: PhantomData,
        }
    }
}

impl Drop for StageSpan {
    fn drop(&mut self) {
        let Some(active) = self.active.take() else {
            return;
        };
        // Profiling must not replace an operation's panic, including if a
        // subscriber panics while emitting a diagnostic during unwinding.
        let _ = std::panic::catch_unwind(|| {
            let inclusive = active.started.elapsed();
            let timing = FRAMES
                .try_with(|frames| {
                    Some(
                        frames
                            .try_borrow_mut()
                            .ok()?
                            .finish(active.metadata, inclusive),
                    )
                })
                .ok()
                .flatten();
            if let Some(timing) = timing {
                tracing::info!(
                    target: "vakint::stage_timing",
                    stage = timing.stage,
                    path = %timing.path,
                    inclusive_seconds = timing.inclusive.as_secs_f64(),
                    exclusive_seconds = timing.exclusive.as_secs_f64(),
                    depth = timing.depth,
                    nesting_valid = timing.nesting_valid,
                    "Completed calculation stage"
                );
            }
        });
    }
}

impl FrameStack {
    const fn new() -> Self {
        Self {
            frames: Vec::new(),
            next_id: 0,
        }
    }

    fn enter(&mut self, stage: &'static str) -> Option<Metadata> {
        self.next_id = self.next_id.checked_add(1)?;
        let metadata = Metadata {
            id: self.next_id,
            stage,
            depth: self.frames.len(),
        };
        let parent = self.frames.last();
        self.frames.push(Frame {
            metadata,
            parent: parent.map(|frame| frame.metadata.id),
            children: Duration::ZERO,
            valid: parent.is_none_or(|frame| frame.valid),
        });
        Some(metadata)
    }

    fn finish(&mut self, metadata: Metadata, inclusive: Duration) -> Timing {
        let Some(index) = self
            .frames
            .iter()
            .position(|frame| frame.metadata.id == metadata.id)
        else {
            return Timing {
                stage: metadata.stage,
                path: metadata.stage.into(),
                depth: metadata.depth,
                inclusive,
                exclusive: inclusive,
                nesting_valid: false,
            };
        };
        if index + 1 != self.frames.len() {
            // Remaining intervals overlap incorrectly after an out-of-order
            // drop. Keep their identities for cleanup, but invalidate timings.
            for frame in &mut self.frames {
                frame.valid = false;
            }
        }
        let path = self.frames[..=index]
            .iter()
            .map(|frame| frame.metadata.stage)
            .collect::<Vec<_>>()
            .join("/");
        let frame = self.frames.remove(index);
        let nesting_valid = frame.valid && frame.children <= inclusive;
        // Use the original parent identity: a removed parent must not cause
        // its child to be counted again against a surviving grandparent.
        if let Some(parent) = frame
            .parent
            .and_then(|id| self.frames.iter_mut().find(|frame| frame.metadata.id == id))
        {
            parent.children = parent.children.saturating_add(inclusive);
            parent.valid &= nesting_valid;
        }
        Timing {
            stage: metadata.stage,
            path,
            depth: metadata.depth,
            inclusive,
            exclusive: inclusive.saturating_sub(frame.children),
            nesting_valid,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nested_exclusive_times_subtract_each_child_once() {
        let mut stack = FrameStack::new();
        let outer = stack.enter("outer").unwrap();
        let child = stack.enter("child").unwrap();
        let grandchild = stack.enter("grandchild").unwrap();
        let grandchild = stack.finish(grandchild, Duration::from_secs(3));
        assert_eq!(grandchild.path, "outer/child/grandchild");
        assert_eq!(grandchild.depth, 2);
        assert_eq!(grandchild.exclusive, Duration::from_secs(3));
        let child = stack.finish(child, Duration::from_secs(7));
        assert_eq!(child.path, "outer/child");
        assert_eq!(child.exclusive, Duration::from_secs(4));
        let sibling = stack.enter("sibling").unwrap();
        let sibling = stack.finish(sibling, Duration::from_secs(2));
        let outer = stack.finish(outer, Duration::from_secs(12));
        assert_eq!(outer.path, "outer");
        assert_eq!(outer.depth, 0);
        assert_eq!(outer.exclusive, Duration::from_secs(3));
        assert!(stack.frames.is_empty());
        assert!(grandchild.nesting_valid && child.nesting_valid);
        assert!(sibling.nesting_valid && outer.nesting_valid);
        assert_eq!(
            grandchild.exclusive + child.exclusive + sibling.exclusive + outer.exclusive,
            outer.inclusive
        );
    }

    #[test]
    fn invalid_nesting_is_reported_and_recovers_without_double_counting() {
        let mut stack = FrameStack::new();
        let outer = stack.enter("outer").unwrap();
        let middle = stack.enter("middle").unwrap();
        let inner = stack.enter("inner").unwrap();
        assert!(!stack.finish(middle, Duration::from_secs(5)).nesting_valid);
        assert!(!stack.finish(inner, Duration::from_secs(3)).nesting_valid);
        let outer = stack.finish(outer, Duration::from_secs(10));
        assert!(!outer.nesting_valid);
        assert_eq!(outer.exclusive, Duration::from_secs(5));
        assert!(stack.frames.is_empty());
        let next = stack.enter("next").unwrap();
        let next = stack.finish(next, Duration::ZERO);
        assert!(next.nesting_valid);
        assert_eq!(next.depth, 0);
        assert_eq!(next.path, "next");
    }

    #[test]
    fn disabled_span_does_not_read_a_clock_or_change_the_stack() {
        FRAMES.with(|frames| {
            let stack = frames.borrow_mut();
            let guard = StageSpan::start("disabled", false, || panic!("clock was read"));
            assert!(guard.active.is_none());
            drop(guard);
            assert!(stack.frames.is_empty());
            assert_eq!(stack.next_id, 0);
        });
    }

    #[test]
    fn unwinding_cleans_the_stack_and_preserves_the_original_panic() {
        let result = std::panic::catch_unwind(|| {
            let _outer = StageSpan::start("outer", true, Instant::now);
            let _inner = StageSpan::start("inner", true, Instant::now);
            panic!("operation failed");
        });
        assert_eq!(
            result.unwrap_err().downcast_ref::<&str>(),
            Some(&"operation failed")
        );
        FRAMES.with(|frames| assert!(frames.borrow().frames.is_empty()));
    }
}
