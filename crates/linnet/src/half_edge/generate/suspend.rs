//! Cooperative checkpoints shared by the native search stages and the binding driver.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Step {
    Ready,
    Pending,
    Yield,
    Done,
}

/// A deterministic work allowance. Scalar topology search may run without a bound;
/// decoration and callback waves use finite allowances for cooperative cancellation.
pub(super) struct Work {
    remaining: Option<usize>,
}

impl Work {
    pub(super) fn unlimited() -> Self {
        Self { remaining: None }
    }

    pub(super) fn bounded(units: usize) -> Self {
        Self {
            remaining: Some(units),
        }
    }

    /// Reserve a unit before mutating an enumeration cursor. Exhaustion leaves the next
    /// attempt untouched, so resuming neither skips it nor repeats an earlier hook.
    pub(super) fn spend(&mut self) -> bool {
        match &mut self.remaining {
            None => true,
            Some(0) => false,
            Some(remaining) => {
                *remaining -= 1;
                true
            }
        }
    }
}
