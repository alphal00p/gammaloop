use symbolica::atom::Atom;
use symbolica::function;

use crate::symbols::S;

#[derive(Debug, Clone, PartialEq)]
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
        p1_sq: Atom,
        p2_sq: Atom,
        p12_sq: Atom,
        m1_sq: Atom,
        m2_sq: Atom,
        m3_sq: Atom,
    },
    Box {
        p1_sq: Atom,
        p2_sq: Atom,
        p3_sq: Atom,
        p4_sq: Atom,
        s: Atom,
        t: Atom,
        m1_sq: Atom,
        m2_sq: Atom,
        m3_sq: Atom,
        m4_sq: Atom,
    },
}

pub trait MasterBasis {
    fn symbol(&self, integral: &MasterIntegral) -> Atom;
}

pub struct OneLoopMasters;

impl MasterBasis for OneLoopMasters {
    fn symbol(&self, integral: &MasterIntegral) -> Atom {
        match integral {
            MasterIntegral::Tadpole { m_sq } => function!(S.a0, m_sq),
            MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => function!(S.b0, p_sq, m1_sq, m2_sq),
            MasterIntegral::Triangle {
                p1_sq,
                p2_sq,
                p12_sq,
                m1_sq,
                m2_sq,
                m3_sq,
            } => function!(S.c0, p1_sq, p2_sq, p12_sq, m1_sq, m2_sq, m3_sq),
            MasterIntegral::Box {
                p1_sq,
                p2_sq,
                p3_sq,
                p4_sq,
                s,
                t,
                m1_sq,
                m2_sq,
                m3_sq,
                m4_sq,
            } => function!(
                S.d0, p1_sq, p2_sq, p3_sq, p4_sq, s, t, m1_sq, m2_sq, m3_sq, m4_sq
            ),
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
    fn bubble_maps_to_symbolic_b0() {
        crate::ensure_symbolica_license();
        let got = OneLoopMasters.symbol(&massless_bubble());
        let want = function!(S.b0, Atom::var(S.psq), Atom::Zero, Atom::Zero);
        assert_eq!(got, want);
    }

    #[test]
    fn tadpole_maps_to_a0() {
        crate::ensure_symbolica_license();
        let m = MasterIntegral::Tadpole { m_sq: Atom::num(1) };
        assert_eq!(OneLoopMasters.symbol(&m), function!(S.a0, Atom::num(1)));
    }

    #[test]
    fn triangle_maps_to_c0() {
        crate::ensure_symbolica_license();
        let m = MasterIntegral::Triangle {
            p1_sq: Atom::num(1),
            p2_sq: Atom::num(2),
            p12_sq: Atom::num(3),
            m1_sq: Atom::Zero,
            m2_sq: Atom::Zero,
            m3_sq: Atom::Zero,
        };
        let want = function!(
            S.c0,
            Atom::num(1),
            Atom::num(2),
            Atom::num(3),
            Atom::Zero,
            Atom::Zero,
            Atom::Zero
        );
        assert_eq!(OneLoopMasters.symbol(&m), want);
    }

    #[test]
    fn box_maps_to_d0() {
        crate::ensure_symbolica_license();
        let m = MasterIntegral::Box {
            p1_sq: Atom::Zero,
            p2_sq: Atom::Zero,
            p3_sq: Atom::Zero,
            p4_sq: Atom::Zero,
            s: Atom::num(1),
            t: Atom::num(2),
            m1_sq: Atom::Zero,
            m2_sq: Atom::Zero,
            m3_sq: Atom::Zero,
            m4_sq: Atom::Zero,
        };
        let want = function!(
            S.d0,
            Atom::Zero,
            Atom::Zero,
            Atom::Zero,
            Atom::Zero,
            Atom::num(1),
            Atom::num(2),
            Atom::Zero,
            Atom::Zero,
            Atom::Zero,
            Atom::Zero
        );
        assert_eq!(OneLoopMasters.symbol(&m), want);
    }
}
