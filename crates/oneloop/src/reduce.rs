use symbolica::atom::Atom;

use crate::family::IntegralFamily;
use crate::masters::MasterIntegral;

pub struct Reduction {
    pub terms: Vec<(Atom, MasterIntegral)>,
}

pub fn reduce(family: &IntegralFamily) -> Reduction {
    match family.propagators.len() {
        2 => {
            let p_sq = family
                .kinematics
                .invariants
                .first()
                .cloned()
                .unwrap_or(Atom::Zero);
            let m1_sq = family.propagators[0].mass_sq.clone();
            let m2_sq = family.propagators[1].mass_sq.clone();
            Reduction {
                terms: vec![(Atom::num(1), MasterIntegral::Bubble { p_sq, m1_sq, m2_sq })],
            }
        }
        _ => todo!("only the one-loop bubble is reduced so far"),
    }
}

#[cfg(test)]
mod tests {
    use super::reduce;
    use crate::family::{Integral, IntegralFamily, Kinematics, Propagator};
    use crate::masters::MasterIntegral;
    use crate::symbols::S;
    use symbolica::atom::Atom;

    fn massless_bubble_family() -> IntegralFamily {
        IntegralFamily {
            propagators: vec![
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: Atom::Zero,
                },
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: Atom::Zero,
                },
            ],
            isps: vec![],
            kinematics: Kinematics {
                invariants: vec![Atom::var(S.psq)],
                masses_sq: vec![],
            },
            targets: vec![Integral {
                propagator_exponents: vec![1, 1],
                isp_exponents: vec![],
            }],
        }
    }

    #[test]
    fn scalar_bubble_reduces_to_unit_b0() {
        crate::ensure_symbolica_license();
        let reduction = reduce(&massless_bubble_family());
        assert_eq!(reduction.terms.len(), 1);
        let (coeff, master) = &reduction.terms[0];
        assert_eq!(*coeff, Atom::num(1));
        match master {
            MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => {
                assert_eq!(*p_sq, Atom::var(S.psq));
                assert_eq!(*m1_sq, Atom::Zero);
                assert_eq!(*m2_sq, Atom::Zero);
            }
            other => panic!("expected a bubble master, got {other:?}"),
        }
    }
}
