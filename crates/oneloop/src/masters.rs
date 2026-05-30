//! One-loop scalar master integrals and their analytic closed forms.

use symbolica::atom::Atom;
use symbolica::function;

use crate::error::OneLoopError;
use crate::symbols::S;

/// A one-loop scalar master integral, keyed by its propagator masses².
#[derive(Debug, Clone)]
pub enum MasterIntegral {
    Tadpole {
        m_sq: Atom,
    },
    Bubble {
        m1_sq: Atom,
        m2_sq: Atom,
    },
    Triangle {
        m1_sq: Atom,
        m2_sq: Atom,
        m3_sq: Atom,
    },
    Box {
        m1_sq: Atom,
        m2_sq: Atom,
        m3_sq: Atom,
        m4_sq: Atom,
    },
}

pub trait MasterBasis {
    /// Is `integral` an irreducible master in this basis
    fn is_master(&self, integral: &MasterIntegral) -> bool;

    /// The ε-dependent analytic closed form as a Symbolica `Atom`
    fn closed_form(&self, integral: &MasterIntegral) -> Result<Atom, OneLoopError>;
}

/// The one-loop master basis: tadpole, bubble, triangle, box.
pub struct OneLoopMasters;

impl MasterBasis for OneLoopMasters {
    fn is_master(&self, integral: &MasterIntegral) -> bool {
        matches!(
            integral,
            MasterIntegral::Tadpole { .. }
                | MasterIntegral::Bubble { .. }
                | MasterIntegral::Triangle { .. }
                | MasterIntegral::Box { .. }
        )
    }

    fn closed_form(&self, integral: &MasterIntegral) -> Result<Atom, OneLoopError> {
        match integral {
            // Massless bubble  B₀(p²; 0, 0) = 1/ε − ln(−p²/μ²) + 2
            MasterIntegral::Bubble { m1_sq, m2_sq }
                if *m1_sq == Atom::Zero && *m2_sq == Atom::Zero =>
            {
                let ep = Atom::var(S.ep);
                let psq = Atom::var(S.psq);
                let musq = Atom::var(S.musq);
                let log_term = function!(S.log, -&psq / &musq);
                Ok(Atom::num(1) / &ep - log_term + Atom::num(2))
            }

            other => Err(OneLoopError::MasterNotInLibrary {
                which: format!("{other:?}"),
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{MasterBasis, MasterIntegral, OneLoopMasters};
    use crate::error::OneLoopError;
    use symbolica::atom::Atom;

    #[test]
    fn massless_bubble_is_a_master() {
        crate::ensure_symbolica_license();
        let basis = OneLoopMasters;
        let m = MasterIntegral::Bubble {
            m1_sq: Atom::Zero,
            m2_sq: Atom::Zero,
        };
        assert!(basis.is_master(&m));
    }

    #[test]
    fn massless_bubble_closed_form_has_pole_and_log() {
        crate::ensure_symbolica_license();
        let basis = OneLoopMasters;
        let m = MasterIntegral::Bubble {
            m1_sq: Atom::Zero,
            m2_sq: Atom::Zero,
        };
        let cf = basis
            .closed_form(&m)
            .expect("massless bubble is in the library");
        let s = cf.to_string();
        assert!(s.contains("ep"), "must carry the 1/ε UV pole, got: {s}");
        assert!(s.contains("log"), "must carry the log, got: {s}");
    }

    #[test]
    fn missing_master_reports_which() {
        crate::ensure_symbolica_license();
        let basis = OneLoopMasters;
        let m = MasterIntegral::Triangle {
            m1_sq: Atom::num(1),
            m2_sq: Atom::num(1),
            m3_sq: Atom::num(1),
        };
        let err = basis.closed_form(&m).unwrap_err();
        assert!(matches!(err, OneLoopError::MasterNotInLibrary { .. }));
    }
}
