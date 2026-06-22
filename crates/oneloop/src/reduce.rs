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
        3 => reduce_triangle(
            exponents[0],
            exponents[1],
            exponents[2],
            &inv(0),
            &inv(1),
            &inv(2),
            &mass(0),
            &mass(1),
            &mass(2),
        ),
        4 => reduce_box(
            exponents[0],
            exponents[1],
            exponents[2],
            exponents[3],
            &inv(0),
            &inv(1),
            &inv(2),
            &inv(3),
            &inv(4),
            &inv(5),
            &mass(0),
            &mass(1),
            &mass(2),
            &mass(3),
        ),
        n => todo!("scalar reduction for {n}-propagator families"),
    };

    Reduction { terms }
}

fn tadpole_coefficient(exponent: i32, m_sq: &Atom) -> Atom {
    if *m_sq == Atom::Zero {
        return Atom::Zero;
    }
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

fn cayley(
    s1: &Atom,
    s2: &Atom,
    s3: &Atom,
    m1: &Atom,
    m2: &Atom,
    m3: &Atom,
) -> (Atom, Atom, Atom, Atom) {
    let y11 = Atom::num(2) * m1;
    let y22 = Atom::num(2) * m2;
    let y33 = Atom::num(2) * m3;
    let y12 = m1 + m2 - s1;
    let y13 = m1 + m3 - s3;
    let y23 = m2 + m3 - s2;

    let m11 = &y22 * &y33 - &y23 * &y23;
    let m12 = &y12 * &y33 - &y13 * &y23;
    let m13 = &y12 * &y23 - &y13 * &y22;
    let det = &y11 * &m11 - &y12 * &m12 + &y13 * &m13;

    let r12 = &y13 * &y23 - &y12 * &y33;
    let r13 = &y12 * &y23 - &y13 * &y22;
    let r23 = &y12 * &y13 - &y11 * &y23;
    (det, r12, r13, r23)
}

#[allow(clippy::too_many_arguments)]
fn diag1(
    a1: i32,
    a2: i32,
    a3: i32,
    s1: &Atom,
    s2: &Atom,
    s3: &Atom,
    m1: &Atom,
    m2: &Atom,
    m3: &Atom,
) -> Atom {
    let p1 = Atom::num(2) * m1 * s2 + m2 * m2 - Atom::num(2) * m2 * m3 + m2 * s1
        - Atom::num(3) * m2 * s2
        - m2 * s3
        + m3 * m3
        - m3 * s1
        - Atom::num(3) * m3 * s2
        + m3 * s3
        - s1 * s2
        + Atom::num(2) * s2 * s2
        - s2 * s3;
    let p2 = m1 * m3 - m1 * m2 + Atom::num(3) * m1 * s2 + m2 * m3 + m2 * s1
        - m2 * s2
        - m3 * m3
        - Atom::num(3) * m3 * s1
        + Atom::num(2) * m3 * s3
        - s1 * s2
        + s2 * s2
        - Atom::num(2) * s2 * s3;
    let p3 =
        m1 * m2 - m1 * m3 + Atom::num(3) * m1 * s2 - m2 * m2 + m2 * m3 + Atom::num(2) * m2 * s1
            - Atom::num(3) * m2 * s3
            - m3 * s2
            + m3 * s3
            - Atom::num(2) * s1 * s2
            + s2 * s2
            - s2 * s3;
    let pd = m2 * s2 - Atom::num(2) * m1 * s2 - m2 * s1 + m2 * s3 + m3 * s1 + m3 * s2 - m3 * s3
        + s1 * s2
        - s2 * s2
        + s2 * s3;
    let pc = Atom::num(2) * m2 * m3 - Atom::num(2) * m1 * s2 - m2 * m2 - m2 * s1
        + Atom::num(3) * m2 * s2
        + m2 * s3
        - m3 * m3
        + m3 * s1
        + Atom::num(3) * m3 * s2
        - m3 * s3
        + s1 * s2
        - Atom::num(2) * s2 * s2
        + s2 * s3;
    Atom::num(i64::from(a1)) * &p1
        + Atom::num(i64::from(a2)) * &p2
        + Atom::num(i64::from(a3)) * &p3
        + Atom::var(S.d) * &pd
        + pc
}

fn promote_bubble(
    r: (Atom, Atom, Atom),
    p_sq: Atom,
    ma: Atom,
    mb: Atom,
) -> Vec<(Atom, MasterIntegral)> {
    let mut terms = Vec::new();
    push_nonzero(
        &mut terms,
        r.0,
        MasterIntegral::Bubble {
            p_sq,
            m1_sq: ma.clone(),
            m2_sq: mb.clone(),
        },
    );
    push_nonzero(&mut terms, r.1, MasterIntegral::Tadpole { m_sq: ma });
    push_nonzero(&mut terms, r.2, MasterIntegral::Tadpole { m_sq: mb });
    terms
}

fn add_scaled(
    acc: &mut Vec<(Atom, MasterIntegral)>,
    coeff: &Atom,
    child: Vec<(Atom, MasterIntegral)>,
) {
    for (c, m) in child {
        let scaled = coeff * &c;
        if scaled == Atom::Zero {
            continue;
        }
        match acc.iter_mut().find(|pair| pair.1 == m) {
            Some(slot) => slot.0 = &slot.0 + &scaled,
            None => acc.push((scaled, m)),
        }
    }
}

// Reduce a triangle C(a1,a2,a3) to C0 + pinched bubbles/tadpoles
#[allow(clippy::too_many_arguments)]
fn reduce_triangle(
    a1: i32,
    a2: i32,
    a3: i32,
    s1: &Atom,
    s2: &Atom,
    s3: &Atom,
    m1: &Atom,
    m2: &Atom,
    m3: &Atom,
) -> Vec<(Atom, MasterIntegral)> {
    if a1 == 1 && a2 == 1 && a3 == 1 {
        return vec![(
            Atom::num(1),
            MasterIntegral::Triangle {
                p1_sq: s1.clone(),
                p2_sq: s2.clone(),
                p12_sq: s3.clone(),
                m1_sq: m1.clone(),
                m2_sq: m2.clone(),
                m3_sq: m3.clone(),
            },
        )];
    }
    if a2 == 0 {
        return promote_bubble(
            reduce_bubble(a1, a3, s3, m1, m3),
            s3.clone(),
            m1.clone(),
            m3.clone(),
        );
    }
    if a1 == 0 {
        return promote_bubble(
            reduce_bubble(a2, a3, s2, m2, m3),
            s2.clone(),
            m2.clone(),
            m3.clone(),
        );
    }
    if a3 == 0 {
        return promote_bubble(
            reduce_bubble(a1, a2, s1, m1, m2),
            s1.clone(),
            m1.clone(),
            m2.clone(),
        );
    }

    let (det, r12, r13, r23) = cayley(s1, s2, s3, m1, m2, m3);
    let (den, table) = if a1 >= 2 && a1 >= a2 && a1 >= a3 {
        let k = kallen(s2, m2, m3);
        (
            Atom::num(i64::from(a1 - 1)) * &det,
            [
                (a1 - 2, a2, a3 + 1, Atom::num(i64::from(a3)) * &k),
                (a1 - 2, a2 + 1, a3, Atom::num(i64::from(a2)) * &k),
                (a1 - 1, a2 - 1, a3 + 1, Atom::num(i64::from(-a3)) * &r12),
                (a1 - 1, a2 + 1, a3 - 1, Atom::num(i64::from(-a2)) * &r13),
                (a1, a2 - 1, a3, Atom::num(i64::from(-(a1 - 1))) * &r12),
                (a1, a2, a3 - 1, Atom::num(i64::from(-(a1 - 1))) * &r13),
                (a1 - 1, a2, a3, diag1(a1, a2, a3, s1, s2, s3, m1, m2, m3)),
            ],
        )
    } else if a2 >= 2 && a2 >= a3 {
        let k = kallen(s3, m1, m3);
        (
            Atom::num(i64::from(a2 - 1)) * &det,
            [
                (a1 - 1, a2 - 1, a3 + 1, Atom::num(i64::from(-a3)) * &r12),
                (a1 - 1, a2, a3, Atom::num(i64::from(-(a2 - 1))) * &r12),
                (a1, a2 - 2, a3 + 1, Atom::num(i64::from(a3)) * &k),
                (a1 + 1, a2 - 2, a3, Atom::num(i64::from(a1)) * &k),
                (a1, a2, a3 - 1, Atom::num(i64::from(-(a2 - 1))) * &r23),
                (a1 + 1, a2 - 1, a3 - 1, Atom::num(i64::from(-a1)) * &r23),
                (a1, a2 - 1, a3, diag1(a2, a3, a1, s2, s3, s1, m2, m3, m1)),
            ],
        )
    } else {
        let k = kallen(s1, m1, m2);
        (
            Atom::num(i64::from(a3 - 1)) * &det,
            [
                (a1 - 1, a2, a3, Atom::num(i64::from(-(a3 - 1))) * &r13),
                (a1 - 1, a2 + 1, a3 - 1, Atom::num(i64::from(-a2)) * &r13),
                (a1, a2 - 1, a3, Atom::num(i64::from(-(a3 - 1))) * &r23),
                (a1 + 1, a2 - 1, a3 - 1, Atom::num(i64::from(-a1)) * &r23),
                (a1, a2 + 1, a3 - 2, Atom::num(i64::from(a2)) * &k),
                (a1 + 1, a2, a3 - 2, Atom::num(i64::from(a1)) * &k),
                (a1, a2, a3 - 1, diag1(a3, a1, a2, s3, s1, s2, m3, m1, m2)),
            ],
        )
    };

    let mut acc: Vec<(Atom, MasterIntegral)> = Vec::new();
    for (b1, b2, b3, num) in table {
        if b1 < 0 || b2 < 0 || b3 < 0 {
            continue;
        }
        let coeff = num / &den;
        add_scaled(
            &mut acc,
            &coeff,
            reduce_triangle(b1, b2, b3, s1, s2, s3, m1, m2, m3),
        );
    }
    acc
}

fn det3(m: &[[Atom; 3]; 3]) -> Atom {
    let c0 = &m[1][1] * &m[2][2] - &m[1][2] * &m[2][1];
    let c1 = &m[1][0] * &m[2][2] - &m[1][2] * &m[2][0];
    let c2 = &m[1][0] * &m[2][1] - &m[1][1] * &m[2][0];
    &m[0][0] * &c0 - &m[0][1] * &c1 + &m[0][2] * &c2
}

fn minor3(y: &[[Atom; 4]; 4], skip_row: usize, skip_col: usize) -> [[Atom; 3]; 3] {
    let rows: Vec<usize> = (0..4).filter(|&r| r != skip_row).collect();
    let cols: Vec<usize> = (0..4).filter(|&c| c != skip_col).collect();
    std::array::from_fn(|r| std::array::from_fn(|c| y[rows[r]][cols[c]].clone()))
}

#[allow(clippy::too_many_arguments)]
fn cayley4(
    p1: &Atom,
    p2: &Atom,
    p3: &Atom,
    p4: &Atom,
    s: &Atom,
    t: &Atom,
    m1: &Atom,
    m2: &Atom,
    m3: &Atom,
    m4: &Atom,
) -> (Atom, [[Atom; 4]; 4]) {
    let two = Atom::num(2);
    let y: [[Atom; 4]; 4] = [
        [&two * m1, m1 + m2 - p1, m1 + m3 - s, m1 + m4 - p4],
        [m1 + m2 - p1, &two * m2, m2 + m3 - p2, m2 + m4 - t],
        [m1 + m3 - s, m2 + m3 - p2, &two * m3, m3 + m4 - p3],
        [m1 + m4 - p4, m2 + m4 - t, m3 + m4 - p3, &two * m4],
    ];

    let adj: [[Atom; 4]; 4] = std::array::from_fn(|k| {
        std::array::from_fn(|i| {
            let minor = det3(&minor3(&y, i, k));
            if (i + k) % 2 == 0 {
                minor
            } else {
                Atom::num(-1) * &minor
            }
        })
    });

    let mut det = Atom::Zero;
    for j in 0..4 {
        det += &y[0][j] * &adj[j][0];
    }
    (det, adj)
}

// Reduce a box D(a1,a2,a3,a4) to D0 + pinched triangles
#[allow(clippy::too_many_arguments)]
fn reduce_box(
    a1: i32,
    a2: i32,
    a3: i32,
    a4: i32,
    p1: &Atom,
    p2: &Atom,
    p3: &Atom,
    p4: &Atom,
    s: &Atom,
    t: &Atom,
    m1: &Atom,
    m2: &Atom,
    m3: &Atom,
    m4: &Atom,
) -> Vec<(Atom, MasterIntegral)> {
    if a4 == 0 {
        return reduce_triangle(a1, a2, a3, p1, p2, s, m1, m2, m3);
    }
    if a1 == 0 {
        return reduce_triangle(a2, a3, a4, p2, p3, t, m2, m3, m4);
    }
    if a2 == 0 {
        return reduce_triangle(a1, a3, a4, s, p3, p4, m1, m3, m4);
    }
    if a3 == 0 {
        return reduce_triangle(a1, a2, a4, p1, t, p4, m1, m2, m4);
    }
    if a1 == 1 && a2 == 1 && a3 == 1 && a4 == 1 {
        return vec![(
            Atom::num(1),
            MasterIntegral::Box {
                p1_sq: p1.clone(),
                p2_sq: p2.clone(),
                p3_sq: p3.clone(),
                p4_sq: p4.clone(),
                s: s.clone(),
                t: t.clone(),
                m1_sq: m1.clone(),
                m2_sq: m2.clone(),
                m3_sq: m3.clone(),
                m4_sq: m4.clone(),
            },
        )];
    }

    let (det, adj) = cayley4(p1, p2, p3, p4, s, t, m1, m2, m3, m4);
    let b = [a1, a2, a3, a4];
    let mut k = 0;
    for i in 1..4 {
        if b[i] > b[k] {
            k = i;
        }
    }
    let mut a = b;
    a[k] -= 1;
    let total: i32 = a.iter().sum();
    let den = Atom::num(i64::from(a[k])) * &det;

    let d = Atom::var(S.d);
    let mut diag = Atom::Zero;
    for i in 0..4 {
        let factor = &d - Atom::num(i64::from(total + a[i]));
        diag += &adj[k][i] * &factor;
    }

    let mut acc: Vec<(Atom, MasterIntegral)> = Vec::new();
    add_scaled(
        &mut acc,
        &(diag / &den),
        reduce_box(a[0], a[1], a[2], a[3], p1, p2, p3, p4, s, t, m1, m2, m3, m4),
    );
    for i in 0..4 {
        for j in 0..4 {
            if i == j {
                continue;
            }
            let num = Atom::num(i64::from(-a[j])) * &adj[k][i];
            let mut ch = a;
            ch[j] += 1;
            ch[i] -= 1;
            add_scaled(
                &mut acc,
                &(num / &den),
                reduce_box(
                    ch[0], ch[1], ch[2], ch[3], p1, p2, p3, p4, s, t, m1, m2, m3, m4,
                ),
            );
        }
    }
    acc
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
    let q1 = Atom::var(S.q1);
    let den1 = symbol!("oneloop::den1");
    let den2 = symbol!("oneloop::den2");
    let d1 = Atom::var(den1);
    let d2 = Atom::var(den2);

    //   l^2 = D1 + m1^2 ;   l.q1 = (D2 - D1 - m1^2 - p_sq + m2^2)/2 ;   q1.q1 = p_sq
    let n = numerator
        .replace(function!(S.dot, &k, &k).to_pattern())
        .with(&d1 + m1_sq)
        .replace(function!(S.dot, &k, &q1).to_pattern())
        .with((&d2 - &d1 - m1_sq - p_sq + m2_sq) / Atom::num(2))
        .replace(function!(S.dot, &q1, &q1).to_pattern())
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

    // [p1, p2, p3, p4, s, t, m1, m2, m3, m4]
    fn box_syms() -> [Atom; 10] {
        [
            Atom::var(symbol!("oneloop::p1")),
            Atom::var(symbol!("oneloop::p2")),
            Atom::var(symbol!("oneloop::p3")),
            Atom::var(symbol!("oneloop::p4")),
            Atom::var(symbol!("oneloop::sinv")),
            Atom::var(symbol!("oneloop::tinv")),
            Atom::var(symbol!("oneloop::m1sq")),
            Atom::var(symbol!("oneloop::m2sq")),
            Atom::var(symbol!("oneloop::m3sq")),
            Atom::var(symbol!("oneloop::m4sq")),
        ]
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
    fn massless_dotted_tadpole_vanishes() {
        crate::ensure_symbolica_license();
        let r = reduce(&family(vec![Atom::Zero], vec![], vec![2]));
        assert_eq!(r.terms.len(), 1);
        assert_eq!(r.terms[0].0, Atom::Zero);
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
            numerator: function!(S.dot, Atom::var(S.k), Atom::var(S.q1)),
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
    fn dotted_triangle_reduces_to_masters() {
        crate::ensure_symbolica_license();
        let s1 = Atom::var(symbol!("oneloop::s1"));
        let s2 = Atom::var(symbol!("oneloop::s2"));
        let s3 = Atom::var(symbol!("oneloop::s3"));
        let m1 = Atom::var(symbol!("oneloop::m1sq"));
        let m2 = Atom::var(symbol!("oneloop::m2sq"));
        let m3 = Atom::var(symbol!("oneloop::m3sq"));
        let r = reduce(&family(vec![m1, m2, m3], vec![s1, s2, s3], vec![2, 2, 2]));
        let triangles = r
            .terms
            .iter()
            .filter(|(_, m)| matches!(m, MasterIntegral::Triangle { .. }))
            .count();
        let bubbles = r
            .terms
            .iter()
            .filter(|(_, m)| matches!(m, MasterIntegral::Bubble { .. }))
            .count();
        assert_eq!(triangles, 1);
        assert!(bubbles >= 1);
        assert!(r.terms.iter().all(|(_, m)| matches!(
            m,
            MasterIntegral::Triangle { .. }
                | MasterIntegral::Bubble { .. }
                | MasterIntegral::Tadpole { .. }
        )));
    }

    #[test]
    fn pinched_triangle_routes_to_the_right_bubble() {
        crate::ensure_symbolica_license();
        let s1 = Atom::var(symbol!("oneloop::s1"));
        let s2 = Atom::var(symbol!("oneloop::s2"));
        let s3 = Atom::var(symbol!("oneloop::s3"));
        let m1 = Atom::var(symbol!("oneloop::m1sq"));
        let m2 = Atom::var(symbol!("oneloop::m2sq"));
        let m3 = Atom::var(symbol!("oneloop::m3sq"));
        // third line pinched -> bubble of lines 1,2 carrying s1
        let r = reduce(&family(
            vec![m1.clone(), m2.clone(), m3],
            vec![s1.clone(), s2, s3],
            vec![1, 1, 0],
        ));
        assert!(r.terms.iter().any(|(_, m)| matches!(
            m,
            MasterIntegral::Bubble { p_sq, m1_sq, m2_sq }
                if *p_sq == s1 && *m1_sq == m1 && *m2_sq == m2
        )));
        assert!(
            r.terms
                .iter()
                .all(|(_, m)| !matches!(m, MasterIntegral::Triangle { .. }))
        );
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
    fn dotted_box_reduces_to_masters() {
        crate::ensure_symbolica_license();
        let p1 = Atom::var(symbol!("oneloop::p1"));
        let p2 = Atom::var(symbol!("oneloop::p2"));
        let p3 = Atom::var(symbol!("oneloop::p3"));
        let p4 = Atom::var(symbol!("oneloop::p4"));
        let s = Atom::var(symbol!("oneloop::sinv"));
        let t = Atom::var(symbol!("oneloop::tinv"));
        let m1 = Atom::var(symbol!("oneloop::m1sq"));
        let m2 = Atom::var(symbol!("oneloop::m2sq"));
        let m3 = Atom::var(symbol!("oneloop::m3sq"));
        let m4 = Atom::var(symbol!("oneloop::m4sq"));
        let r = reduce(&family(
            vec![m1, m2, m3, m4],
            vec![p1, p2, p3, p4, s, t],
            vec![2, 2, 1, 1],
        ));
        let boxes = r
            .terms
            .iter()
            .filter(|(_, m)| matches!(m, MasterIntegral::Box { .. }))
            .count();
        let triangles = r
            .terms
            .iter()
            .filter(|(_, m)| matches!(m, MasterIntegral::Triangle { .. }))
            .count();
        assert_eq!(boxes, 1);
        assert!(triangles >= 1);
        assert!(r.terms.iter().all(|(_, m)| matches!(
            m,
            MasterIntegral::Box { .. }
                | MasterIntegral::Triangle { .. }
                | MasterIntegral::Bubble { .. }
                | MasterIntegral::Tadpole { .. }
        )));
    }

    #[test]
    fn pinched_box_routes_to_a_triangle() {
        crate::ensure_symbolica_license();
        let [p1, p2, p3, p4, s, t, m1, m2, m3, m4] = box_syms();

        let r = reduce(&family(
            vec![m1.clone(), m2.clone(), m3.clone(), m4],
            vec![p1.clone(), p2.clone(), p3, p4, s.clone(), t],
            vec![1, 1, 1, 0],
        ));
        assert!(r.terms.iter().any(|(_, m)| matches!(
            m,
            MasterIntegral::Triangle { p1_sq, p2_sq, p12_sq, m1_sq, m2_sq, m3_sq }
                if *p1_sq == p1 && *p2_sq == p2 && *p12_sq == s
                    && *m1_sq == m1 && *m2_sq == m2 && *m3_sq == m3
        )));
        assert!(
            r.terms
                .iter()
                .all(|(_, m)| !matches!(m, MasterIntegral::Box { .. }))
        );
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
