use symbolica::atom::{Atom, AtomCore};
use symbolica::{function, symbol};

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
    let exponents = &family.targets[0].propagator_exponents;

    if family.numerator != Atom::num(1) && family.propagators.len() != 2 {
        todo!("numerator reduction for bubble only so far");
    }

    let terms = match family.propagators.len() {
        1 => {
            let m_sq = mass(0);
            vec![(
                tadpole_coefficient(exponents[0], &m_sq),
                MasterIntegral::Tadpole { m_sq },
            )]
        }
        2 => {
            let p_sq = inv(0);
            let m1_sq = mass(0);
            let m2_sq = mass(1);
            let (c_b0, c_a1, c_a2) = if family.numerator == Atom::num(1) {
                reduce_bubble(exponents[0], exponents[1], &p_sq, &m1_sq, &m2_sq)
            } else {
                bubble_numerator(
                    &family.numerator,
                    exponents[0],
                    exponents[1],
                    &p_sq,
                    &m1_sq,
                    &m2_sq,
                )
            };
            let mut terms = Vec::new();
            push_nonzero(
                &mut terms,
                c_b0,
                MasterIntegral::Bubble {
                    p_sq,
                    m1_sq: m1_sq.clone(),
                    m2_sq: m2_sq.clone(),
                },
            );
            push_nonzero(&mut terms, c_a1, MasterIntegral::Tadpole { m_sq: m1_sq });
            push_nonzero(&mut terms, c_a2, MasterIntegral::Tadpole { m_sq: m2_sq });
            terms
        }
        3 => vec![(
            Atom::num(1),
            MasterIntegral::Triangle {
                p1_sq: inv(0),
                p2_sq: inv(1),
                p12_sq: inv(2),
                m1_sq: mass(0),
                m2_sq: mass(1),
                m3_sq: mass(2),
            },
        )],
        4 => vec![(
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
        )],
        n => todo!("scalar reduction for {n}-propagator families"),
    };

    Reduction { terms }
}

fn tadpole_coefficient(exponent: i32, m_sq: &Atom) -> Atom {
    let d = Atom::var(S.d);
    let mut coeff = Atom::num(1);
    for k in 1..i64::from(exponent) {
        coeff = coeff * (&d - Atom::num(2 * k)) / (Atom::num(2 * k) * m_sq);
    }
    coeff
}

fn push_nonzero(terms: &mut Vec<(Atom, MasterIntegral)>, coeff: Atom, master: MasterIntegral) {
    if coeff != Atom::Zero {
        terms.push((coeff, master));
    }
}

// Reduce a bubble B(a1,a2) to (coeff of B0, coeff of A0(m1^2), coeff of A0(m2^2))
fn reduce_bubble(a1: i32, a2: i32, p_sq: &Atom, m1_sq: &Atom, m2_sq: &Atom) -> (Atom, Atom, Atom) {
    if a1 == 1 && a2 == 1 {
        return (Atom::num(1), Atom::Zero, Atom::Zero);
    }
    if a2 == 0 {
        return (Atom::Zero, tadpole_coefficient(a1, m1_sq), Atom::Zero);
    }
    if a1 == 0 {
        return (Atom::Zero, Atom::Zero, tadpole_coefficient(a2, m2_sq));
    }

    let d = Atom::var(S.d);
    let lam = kallen(p_sq, m1_sq, m2_sq);
    let pinch = (p_sq - m1_sq - m2_sq) / &lam;

    if a1 >= 2 {
        let den = Atom::num(i64::from(a1 - 1)) * &lam;
        let c = (&d + Atom::num(i64::from(1 - a1 - 2 * a2))) * m1_sq
            + (Atom::num(i64::from(3 * a1 - 3)) - &d) * m2_sq
            + (Atom::num(i64::from(a1 + 2 * a2 - 1)) - &d) * p_sq;
        combine(
            &(Atom::num(i64::from(2 * a2)) * m2_sq / &den),
            reduce_bubble(a1 - 2, a2 + 1, p_sq, m1_sq, m2_sq),
            &(c / &den),
            reduce_bubble(a1 - 1, a2, p_sq, m1_sq, m2_sq),
            &pinch,
            reduce_bubble(a1, a2 - 1, p_sq, m1_sq, m2_sq),
        )
    } else {
        let den = Atom::num(i64::from(a2 - 1)) * &lam;
        let c = (&d + Atom::num(i64::from(1 - a2 - 2 * a1))) * m2_sq
            + (Atom::num(i64::from(3 * a2 - 3)) - &d) * m1_sq
            + (Atom::num(i64::from(a2 + 2 * a1 - 1)) - &d) * p_sq;
        combine(
            &(Atom::num(i64::from(2 * a1)) * m1_sq / &den),
            reduce_bubble(a1 + 1, a2 - 2, p_sq, m1_sq, m2_sq),
            &(c / &den),
            reduce_bubble(a1, a2 - 1, p_sq, m1_sq, m2_sq),
            &pinch,
            reduce_bubble(a1 - 1, a2, p_sq, m1_sq, m2_sq),
        )
    }
}

fn combine(
    k1: &Atom,
    r1: (Atom, Atom, Atom),
    k2: &Atom,
    r2: (Atom, Atom, Atom),
    k3: &Atom,
    r3: (Atom, Atom, Atom),
) -> (Atom, Atom, Atom) {
    (
        k1 * &r1.0 + k2 * &r2.0 + k3 * &r3.0,
        k1 * &r1.1 + k2 * &r2.1 + k3 * &r3.1,
        k1 * &r1.2 + k2 * &r2.2 + k3 * &r3.2,
    )
}

fn kallen(p_sq: &Atom, m1_sq: &Atom, m2_sq: &Atom) -> Atom {
    p_sq * p_sq + m1_sq * m1_sq + m2_sq * m2_sq
        - Atom::num(2) * p_sq * m1_sq
        - Atom::num(2) * p_sq * m2_sq
        - Atom::num(2) * m1_sq * m2_sq
}

// Reduce a bubble with a numerator
fn bubble_numerator(
    numerator: &Atom,
    a1: i32,
    a2: i32,
    p_sq: &Atom,
    m1_sq: &Atom,
    m2_sq: &Atom,
) -> (Atom, Atom, Atom) {
    let k = Atom::var(S.k);
    let p = Atom::var(S.p);
    let den1 = symbol!("oneloop::den1");
    let den2 = symbol!("oneloop::den2");
    let d1 = Atom::var(den1);
    let d2 = Atom::var(den2);

    //   l^2 = D1 + m1^2 ;   l.p = (D2 - D1 - m1^2 - p^2 + m2^2)/2 ;   p^2 = p_sq
    let n = numerator
        .replace(function!(S.dot, &k, &k).to_pattern())
        .with(&d1 + m1_sq)
        .replace(function!(S.dot, &k, &p).to_pattern())
        .with((&d2 - &d1 - m1_sq - p_sq + m2_sq) / Atom::num(2))
        .replace(function!(S.dot, &p, &p).to_pattern())
        .with(p_sq.clone());

    // n is linear in D1, D2:  n = c_const + c_den1*D1 + c_den2*D2
    let c_den1 = n.derivative(den1);
    let c_den2 = n.derivative(den2);
    let c_const = n
        .replace(d1.to_pattern())
        .with(Atom::Zero)
        .replace(d2.to_pattern())
        .with(Atom::Zero);

    combine(
        &c_const,
        reduce_bubble(a1, a2, p_sq, m1_sq, m2_sq),
        &c_den1,
        reduce_bubble(a1 - 1, a2, p_sq, m1_sq, m2_sq),
        &c_den2,
        reduce_bubble(a1, a2 - 1, p_sq, m1_sq, m2_sq),
    )
}

#[cfg(test)]
mod tests {
    use super::reduce;
    use crate::family::{Integral, IntegralFamily, Kinematics, Propagator};
    use crate::masters::MasterIntegral;
    use crate::symbols::S;
    use symbolica::atom::Atom;
    use symbolica::function;
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
            numerator: Atom::num(1),
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
    fn dotted_bubble_reduces_to_a_bubble_and_two_tadpoles() {
        crate::ensure_symbolica_license();
        let psq = Atom::var(S.psq);
        let m1 = Atom::var(symbol!("oneloop::m1sq"));
        let m2 = Atom::var(symbol!("oneloop::m2sq"));
        let r = reduce(&family(vec![m1, m2], vec![psq], vec![3, 1]));
        assert_eq!(r.terms.len(), 3);
        assert!(matches!(r.terms[0].1, MasterIntegral::Bubble { .. }));
        assert!(matches!(r.terms[1].1, MasterIntegral::Tadpole { .. }));
        assert!(matches!(r.terms[2].1, MasterIntegral::Tadpole { .. }));
    }

    #[test]
    fn bubble_with_linear_numerator_reduces_to_masters() {
        crate::ensure_symbolica_license();
        let psq = Atom::var(S.psq);
        let m1 = Atom::var(symbol!("oneloop::m1sq"));
        let m2 = Atom::var(symbol!("oneloop::m2sq"));
        // numerator = l . p
        let fam = IntegralFamily {
            propagators: vec![
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: m1,
                },
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: m2,
                },
            ],
            isps: vec![],
            kinematics: Kinematics {
                invariants: vec![psq],
                masses_sq: vec![],
            },
            targets: vec![Integral {
                propagator_exponents: vec![1, 1],
                isp_exponents: vec![],
            }],
            numerator: function!(S.dot, Atom::var(S.k), Atom::var(S.p)),
        };
        let r = reduce(&fam);
        assert_eq!(r.terms.len(), 3);
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, MasterIntegral::Bubble { .. }))
        );
        assert_eq!(
            r.terms
                .iter()
                .filter(|(_, m)| matches!(m, MasterIntegral::Tadpole { .. }))
                .count(),
            2
        );
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
