//! Source-data precision diagnostics, distinct from working precision.
//!
//! Native float metadata is inspected before runtime resizing. Neither that
//! metadata nor a larger working precision is an error bound on a final result.

use std::collections::HashSet;

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::coefficient::CoefficientView;
use symbolica::id::{Condition, PatternRestriction};

use crate::VakintSettings;

/// Per-evaluation/table-pass deduplication; warnings never change arithmetic.
pub(crate) struct MasterPrecisionWarnings<'a> {
    settings: &'a VakintSettings,
    seen: HashSet<Atom>,
}

impl<'a> MasterPrecisionWarnings<'a> {
    pub(crate) fn new(settings: &'a VakintSettings) -> Self {
        Self {
            settings,
            seen: HashSet::new(),
        }
    }

    /// Return whether a new warning was emitted.
    /// Exact rationals and exact zero components impose no finite-data limit.
    pub(crate) fn check(&mut self, source: &Atom, value: &Atom) -> bool {
        if !self.seen.insert(source.clone()) {
            return false;
        }
        let mut available = u32::MAX;
        value.visitor(&mut |part| {
            if let AtomView::Num(number) = part
                && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
            {
                for component in [real, imaginary] {
                    let component = component.to_float();
                    if !component.as_raw().is_zero() {
                        available = available.min(component.prec());
                    }
                }
            }
            true
        });
        if available < self.settings.get_binary_precision() {
            log::warn!(
                "WARNING: master source {} provides only {available} bits of stored precision; {} decimal digits ({} working bits) were requested. Continuing at the requested working precision; this does not supply additional accurate source digits or guarantee final-result accuracy.",
                source.to_canonical_string(),
                self.settings.run_time_decimal_precision,
                self.settings.get_binary_precision(),
            );
            true
        } else {
            false
        }
    }

    /// Follow only reachable table records, including dependent constants.
    /// Inspect raw targets before any runtime precision resize; no substitutions
    /// or numerical evaluation of the user's expression are performed here.
    pub(crate) fn check_substitutions<'b>(
        &mut self,
        expression: AtomView<'b>,
        substitutions: impl IntoIterator<
            Item = (
                &'b Atom,
                &'b Atom,
                Option<&'b Condition<PatternRestriction>>,
            ),
        >,
    ) {
        let substitutions: Vec<_> = substitutions.into_iter().collect();
        let mut pending = vec![expression];
        while let Some(value) = pending.pop() {
            for &(source, target, condition) in &substitutions {
                if !self.seen.contains(source)
                    && value
                        .pattern_match(&source.to_pattern(), condition, None)
                        .next()
                        .is_some()
                {
                    self.check(source, target);
                    pending.push(target.as_view());
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;
