use symbolica::atom::Atom;

use crate::family::IntegralFamily;
use crate::masters::MasterIntegral;
use crate::symbols::S;

pub struct Reduction {
    pub terms: Vec<(Atom, MasterIntegral)>,
}

pub fn reduce(family: &IntegralFamily) -> Reduction {
    let inv = |i: usize| {
        family
            .kinematics
            .invariants
            .get(i)
            .cloned()
            .unwrap_or(Atom::Zero)
    };
    let mass = |i: usize| family.propagators[i].mass_sq.clone();

    let term = match family.propagators.len() {
        1 => {
            let m_sq = mass(0);
            let coeff = tadpole_coefficient(family.targets[0].propagator_exponents[0], &m_sq);
            (coeff, MasterIntegral::Tadpole { m_sq })
        }
        2 => (
            Atom::num(1),
            MasterIntegral::Bubble {
                p_sq: inv(0),
                m1_sq: mass(0),
                m2_sq: mass(1),
            },
        ),
        3 => (
            Atom::num(1),
            MasterIntegral::Triangle {
                p1_sq: inv(0),
                p2_sq: inv(1),
                p12_sq: inv(2),
                m1_sq: mass(0),
                m2_sq: mass(1),
                m3_sq: mass(2),
            },
        ),
        4 => (
            Atom::num(1),
            MasterIntegral::Box {
                p1_sq: inv(0),
                p2_sq: inv(1),
                p3_sq: inv(2),
                p4_sq: inv(3),
                s: inv(4),
                t: inv(5),
                m1_sq: mass(0),
                m2_sq: mass(1),
                m3_sq: mass(2),
                m4_sq: mass(3),
            },
        ),
        n => todo!("scalar reduction for {n}-propagator families"),
    };

    Reduction { terms: vec![term] }
}

fn tadpole_coefficient(exponent: i32, m_sq: &Atom) -> Atom {
    let d = Atom::var(S.d);
    let mut coeff = Atom::num(1);
    for k in 1..i64::from(exponent) {
        coeff = coeff * (&d - Atom::num(2 * k)) / (Atom::num(2 * k) * m_sq);
    }
    coeff
}

#[cfg(test)]
mod tests {
    use super::reduce;
    use crate::family::{Integral, IntegralFamily, Kinematics, Propagator};
    use crate::masters::MasterIntegral;
    use crate::symbols::S;
    use symbolica::atom::Atom;
    use symbolica::symbol;

    fn family(masses: Vec<Atom>, invariants: Vec<Atom>, exponents: Vec<i32>) -> IntegralFamily {
        IntegralFamily {
            propagators: masses
                .into_iter()
                .map(|mass_sq| Propagator {
                    momentum: Atom::Zero,
                    mass_sq,
                })
                .collect(),
            isps: vec![],
            kinematics: Kinematics {
                invariants,
                masses_sq: vec![],
            },
            targets: vec![Integral {
                propagator_exponents: exponents,
                isp_exponents: vec![],
            }],
        }
    }

    fn scalar_family(masses: Vec<Atom>, invariants: Vec<Atom>) -> IntegralFamily {
        let n = masses.len();
        family(masses, invariants, vec![1; n])
    }

    #[test]
    fn scalar_tadpole_reduces_to_unit_a0() {
        crate::ensure_symbolica_license();
        let r = reduce(&scalar_family(vec![Atom::num(1)], vec![]));
        assert_eq!(r.terms.len(), 1);
        let (coeff, master) = &r.terms[0];
        assert_eq!(*coeff, Atom::num(1));
        match master {
            MasterIntegral::Tadpole { m_sq } => assert_eq!(*m_sq, Atom::num(1)),
            other => panic!("expected a tadpole master, got {other:?}"),
        }
    }

    #[test]
    fn dotted_tadpole_reduces_with_recursion_coefficient() {
        crate::ensure_symbolica_license();
        let msq = Atom::var(symbol!("oneloop::msq"));
        let r = reduce(&family(vec![msq.clone()], vec![], vec![2]));
        assert_eq!(r.terms.len(), 1);
        let (coeff, master) = &r.terms[0];
        let d = Atom::var(S.d);
        let want = (&d - Atom::num(2)) / (Atom::num(2) * &msq);
        assert_eq!(*coeff, want);
        match master {
            MasterIntegral::Tadpole { m_sq } => assert_eq!(*m_sq, msq),
            other => panic!("expected a tadpole master, got {other:?}"),
        }
    }

    #[test]
    fn scalar_bubble_reduces_to_unit_b0() {
        crate::ensure_symbolica_license();
        let r = reduce(&scalar_family(
            vec![Atom::Zero, Atom::Zero],
            vec![Atom::var(S.psq)],
        ));
        assert_eq!(r.terms.len(), 1);
        let (coeff, master) = &r.terms[0];
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

    #[test]
    fn scalar_triangle_reduces_to_unit_c0() {
        crate::ensure_symbolica_license();
        let r = reduce(&scalar_family(
            vec![Atom::Zero, Atom::Zero, Atom::Zero],
            vec![Atom::num(1), Atom::num(2), Atom::num(3)],
        ));
        assert_eq!(r.terms.len(), 1);
        let (coeff, master) = &r.terms[0];
        assert_eq!(*coeff, Atom::num(1));
        match master {
            MasterIntegral::Triangle {
                p1_sq,
                p2_sq,
                p12_sq,
                ..
            } => {
                assert_eq!(*p1_sq, Atom::num(1));
                assert_eq!(*p2_sq, Atom::num(2));
                assert_eq!(*p12_sq, Atom::num(3));
            }
            other => panic!("expected a triangle master, got {other:?}"),
        }
    }

    #[test]
    fn scalar_box_reduces_to_unit_d0() {
        crate::ensure_symbolica_license();
        let r = reduce(&scalar_family(
            vec![Atom::Zero, Atom::Zero, Atom::Zero, Atom::Zero],
            vec![
                Atom::num(1),
                Atom::num(2),
                Atom::num(3),
                Atom::num(4),
                Atom::num(5),
                Atom::num(6),
            ],
        ));
        assert_eq!(r.terms.len(), 1);
        let (coeff, master) = &r.terms[0];
        assert_eq!(*coeff, Atom::num(1));
        match master {
            MasterIntegral::Box { s, t, .. } => {
                assert_eq!(*s, Atom::num(5));
                assert_eq!(*t, Atom::num(6));
            }
            other => panic!("expected a box master, got {other:?}"),
        }
    }
}
