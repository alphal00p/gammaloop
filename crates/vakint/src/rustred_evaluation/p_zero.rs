//! Zero-sector pruning of the direct lowering ([`super::p_lowering`]).
//!
//! RustRed's census gives every sector `(presentation, support)` of a library
//! presentation a fate; a key whose sector is `Zero` (scaleless) reduces to
//! nothing: [`rustred`]'s engine admits it and its step is empty. The direct
//! lowering emits such keys like every other one. Here the census of the
//! lowering's top presentation becomes a bit table over the supports
//! ([`ZeroSectors`]) so that the lowering can drop them:
//!
//! - [`ZeroPrune::End`]: a lowered key in a zero sector is not emitted;
//! - [`ZeroPrune::Early`]: in addition, a partial product of the Horner slot
//!   substitution is dropped the moment a slot factor cancels a line and the
//!   remaining support is a zero sector.
//!
//! Why the early drop is exact. A key's power on a line slot is the request's
//! line power minus the slot's exponent; the support is the set of lines with
//! a positive power. Numerator factors only raise slot exponents, so every
//! descendant of a partial product has a support contained in the partial
//! product's support. The table is accepted only when the census holds every
//! support of the presentation and its zero set is closed under cancelling
//! lines (a subset of a zero support is zero, checked on the census itself,
//! [`ZeroSectors::from_census`]); then every descendant of a dropped term is a
//! key of a zero sector, and every key in a non-zero sector receives exactly
//! the contributions it receives without pruning (exact arithmetic: same keys,
//! same values). When the table is not accepted the oracle is unknown and
//! nothing is pruned.
//!
//! The switch is `VAKINT_P_ZERO_PRUNE` = `off` (default) | `end` | `early`
//! (`on`, `1`, `true` mean `early`) | `check`. Off is today's lowering, byte
//! for byte. `check` (a gate) prunes as `early` and also substitutes every
//! bucket without pruning: the keys outside zero sectors and their values
//! must be equal, else the lowering fails.

use std::cell::Cell;
use std::sync::OnceLock;

use rustred::propagator::PresentationId;
use rustred::propagator::census::PSectorFate;
use rustred::propagator::library::PLibrary;
use rustred::propagator::presentation::Slot;

/// Which zero-sector pruning the direct lowering applies.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum ZeroPrune {
    /// Every lowered key is emitted (today's lowering).
    #[default]
    Off,
    /// Lowered keys in zero sectors are not emitted.
    End,
    /// As [`Self::End`], and partial products of the slot substitution whose
    /// cancelled lines already make their sector zero are dropped with all
    /// their descendants.
    Early,
    /// [`Self::Early`], checked against the unpruned substitution of every
    /// bucket (twice the substitution work: a gate, not a production mode).
    Check,
}

impl ZeroPrune {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::End => "end",
            Self::Early => "early",
            Self::Check => "check",
        }
    }

    /// Partial products are dropped ([`Self::Early`] and [`Self::Check`]).
    pub(crate) fn early(self) -> bool {
        matches!(self, Self::Early | Self::Check)
    }

    fn parse(value: &str) -> Option<Self> {
        match value.trim().to_ascii_lowercase().as_str() {
            "" | "off" | "0" | "false" => Some(Self::Off),
            "end" => Some(Self::End),
            "early" | "on" | "1" | "true" => Some(Self::Early),
            "check" => Some(Self::Check),
            _ => None,
        }
    }
}

/// The pruning when `VAKINT_P_ZERO_PRUNE` is unset.
pub(crate) const DEFAULT_ZERO_PRUNE: ZeroPrune = ZeroPrune::Off;

/// The environment switch of the pruning (read once per process).
pub(crate) const ZERO_PRUNE_VARIABLE: &str = "VAKINT_P_ZERO_PRUNE";

thread_local! {
    /// Tests: a per-thread choice that takes precedence over the environment.
    static OVERRIDE: Cell<Option<ZeroPrune>> = const { Cell::new(None) };
}

/// The pruning of this lowering: a test override, else `VAKINT_P_ZERO_PRUNE`,
/// else [`DEFAULT_ZERO_PRUNE`].
pub(crate) fn zero_prune() -> ZeroPrune {
    if let Some(mode) = OVERRIDE.with(Cell::get) {
        return mode;
    }
    static FROM_ENVIRONMENT: OnceLock<ZeroPrune> = OnceLock::new();
    *FROM_ENVIRONMENT.get_or_init(|| match std::env::var(ZERO_PRUNE_VARIABLE) {
        Err(_) => DEFAULT_ZERO_PRUNE,
        Ok(value) => ZeroPrune::parse(&value).unwrap_or_else(|| {
            eprintln!(
                "vakint: {ZERO_PRUNE_VARIABLE}={value:?} is not off, end, early or check; using {}",
                DEFAULT_ZERO_PRUNE.as_str()
            );
            DEFAULT_ZERO_PRUNE
        }),
    })
}

/// Run `f` with the pruning `mode` on this thread (tests).
#[cfg(test)]
pub(crate) fn with_zero_prune<R>(mode: ZeroPrune, f: impl FnOnce() -> R) -> R {
    struct Restore(Option<ZeroPrune>);
    impl Drop for Restore {
        fn drop(&mut self) {
            OVERRIDE.with(|cell| cell.set(self.0));
        }
    }
    let _restore = Restore(OVERRIDE.with(|cell| cell.replace(Some(mode))));
    f()
}

/// The oracle has a table for at most this many lines (2^24 bits, 2 MiB).
pub(crate) const MAX_TABLE_LINES: usize = 24;

/// The zero sectors of one presentation whose first `lines` slots are its
/// lines (vakint's top presentations): bit `support` is set iff the census
/// gives `(presentation, support)` the fate `Zero`.
#[derive(Clone, Debug)]
pub(crate) struct ZeroSectors {
    pub lines: usize,
    /// Supports with the fate `Zero` (of `2^lines`).
    pub zero_count: usize,
    zero: Vec<u64>,
}

impl ZeroSectors {
    /// The table of `presentation` from `library.census()`, or why there is
    /// none (the oracle is then unknown): more than [`MAX_TABLE_LINES`] lines,
    /// lines that are not the first `lines` slots, ε-lines (forced axes), a
    /// support of the presentation missing from the census, or a zero set
    /// that is not closed under cancelling a line.
    pub(crate) fn from_census(
        library: &PLibrary,
        presentation: PresentationId,
        lines: usize,
    ) -> Result<Self, &'static str> {
        let p = library.presentation(presentation);
        if lines == 0 || lines > MAX_TABLE_LINES || lines > p.slots.len() {
            return Err("a line count outside the table range");
        }
        if p.slots
            .iter()
            .enumerate()
            .any(|(index, slot)| matches!(slot, Slot::Line(_)) != (index < lines))
        {
            return Err("the lines are not the first slots");
        }
        if p.forced_axes() != 0 {
            return Err("a presentation with forced (epsilon) lines");
        }
        let count = 1_usize << lines;
        let mut table = Self {
            lines,
            zero_count: 0,
            zero: vec![0; count.div_ceil(64)],
        };
        let mut held = 0_usize;
        for (&(_, support), fate) in library
            .census()
            .sectors
            .range((presentation, 0)..=(presentation, count as u64 - 1))
        {
            held += 1;
            if matches!(fate, PSectorFate::Zero) {
                table.zero[(support >> 6) as usize] |= 1 << (support & 63);
                table.zero_count += 1;
            }
        }
        if held != count {
            return Err("the census does not hold every support of the presentation");
        }
        if table.not_closed().is_some() {
            return Err("the zero sectors are not closed under cancelling lines");
        }
        Ok(table)
    }

    /// A zero support with a non-zero support one line smaller, if any: then
    /// the zero set is not closed under cancelling lines (closed under
    /// cancelling one line at a time is closed under cancelling any set).
    fn not_closed(&self) -> Option<(u64, u64)> {
        for support in 0..1_u64 << self.lines {
            if !self.is_zero(support) {
                continue;
            }
            let mut rest = support;
            while rest != 0 {
                let line = rest & rest.wrapping_neg();
                rest ^= line;
                if !self.is_zero(support ^ line) {
                    return Some((support, support ^ line));
                }
            }
        }
        None
    }

    /// Whether the sector `support` (bit `s` = line slot `s` has a positive
    /// power) is a zero sector.
    #[inline]
    pub(crate) fn is_zero(&self, support: u64) -> bool {
        (self.zero[(support >> 6) as usize] >> (support & 63)) & 1 == 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_switch_reads_its_values() {
        assert_eq!(ZeroPrune::parse("off"), Some(ZeroPrune::Off));
        assert_eq!(ZeroPrune::parse(""), Some(ZeroPrune::Off));
        assert_eq!(ZeroPrune::parse("0"), Some(ZeroPrune::Off));
        assert_eq!(ZeroPrune::parse("end"), Some(ZeroPrune::End));
        assert_eq!(ZeroPrune::parse("Early"), Some(ZeroPrune::Early));
        assert_eq!(ZeroPrune::parse("on"), Some(ZeroPrune::Early));
        assert_eq!(ZeroPrune::parse("1"), Some(ZeroPrune::Early));
        assert_eq!(ZeroPrune::parse("check"), Some(ZeroPrune::Check));
        assert_eq!(ZeroPrune::parse("yes please"), None);
        assert!(ZeroPrune::Check.early() && ZeroPrune::Early.early());
        assert!(!ZeroPrune::End.early() && !ZeroPrune::Off.early());
        assert_eq!(DEFAULT_ZERO_PRUNE, ZeroPrune::Off);
        let outside = zero_prune();
        with_zero_prune(ZeroPrune::Early, || {
            assert_eq!(zero_prune(), ZeroPrune::Early);
            with_zero_prune(ZeroPrune::End, || assert_eq!(zero_prune(), ZeroPrune::End));
            assert_eq!(zero_prune(), ZeroPrune::Early);
        });
        assert_eq!(zero_prune(), outside);
    }

    /// The closure test of a table: zero supports of a three-line
    /// presentation given by hand.
    #[test]
    fn a_zero_set_must_be_closed_under_cancelling_lines() {
        let table = |zero: &[u64]| {
            let mut table = ZeroSectors {
                lines: 3,
                zero_count: zero.len(),
                zero: vec![0],
            };
            for support in zero {
                table.zero[0] |= 1 << support;
            }
            table
        };
        // {}, {0}, {1}, {0,1} zero: closed; {0,1,2} etc. non-zero.
        let closed = table(&[0b000, 0b001, 0b010, 0b011]);
        assert_eq!(closed.not_closed(), None);
        assert!(closed.is_zero(0b011) && !closed.is_zero(0b111) && !closed.is_zero(0b100));
        // {0,1} zero but {1} not.
        let open = table(&[0b000, 0b001, 0b011]);
        assert_eq!(open.not_closed(), Some((0b011, 0b010)));
    }
}
