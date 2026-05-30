use symbolica::atom::Atom;
use symbolica::function;

use crate::symbols::S;

#[derive(Debug, Clone)]
pub enum MasterIntegral {
    Tadpole {
        m_sq: Atom,
    },
    Bubble {
        p_sq: Atom,
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
    fn is_master(&self, integral: &MasterIntegral) -> bool;
    fn symbol(&self, integral: &MasterIntegral) -> Atom;
}

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

    fn symbol(&self, integral: &MasterIntegral) -> Atom {
        match integral {
            MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => function!(S.b0, p_sq, m1_sq, m2_sq),
            _ => todo!("A0/C0/D0 symbolic forms land with their masters (M4)"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{MasterBasis, MasterIntegral, OneLoopMasters};
    use crate::symbols::S;
    use symbolica::atom::Atom;
    use symbolica::function;

    fn massless_bubble() -> MasterIntegral {
        MasterIntegral::Bubble {
            p_sq: Atom::var(S.psq),
            m1_sq: Atom::Zero,
            m2_sq: Atom::Zero,
        }
    }

    #[test]
    fn massless_bubble_is_a_master() {
        crate::ensure_symbolica_license();
        assert!(OneLoopMasters.is_master(&massless_bubble()));
    }

    #[test]
    fn bubble_maps_to_symbolic_b0() {
        crate::ensure_symbolica_license();
        let got = OneLoopMasters.symbol(&massless_bubble());
        let want = function!(S.b0, Atom::var(S.psq), Atom::Zero, Atom::Zero);
        assert_eq!(got, want);
    }
}
