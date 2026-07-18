use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::domains::integer::{IntegerRing, Z};
use symbolica::domains::rational::Q;
use symbolica::domains::rational_polynomial::{RationalPolynomial, RationalPolynomialField};
use symbolica::tensors::matrix::Matrix;
use symbolica::{function, symbol};

use crate::family::IntegralFamily;
use crate::masters::MasterIntegral;
use crate::symbols::S;

pub struct Reduction {
    pub terms: Vec<(Atom, MasterIntegral)>,
}

impl Reduction {
    /// Reduce every coefficient to lowest terms
    pub fn simplify(mut self) -> Self {
        for (coeff, _) in &mut self.terms {
            *coeff = coeff.cancel();
        }
        self
    }
}

/// Reduce a one-loop integral family to the scalar masters `A0`/`B0`/`C0`/`D0`.
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

    let terms = match family.propagators.len() {
        1 => {
            let m_sq = mass(0);
            if family.numerator == Atom::num(1) {
                vec![(
                    tadpole_coefficient(exponents[0], &m_sq),
                    MasterIntegral::Tadpole { m_sq },
                )]
            } else {
                tadpole_numerator(&family.numerator, exponents[0], &m_sq)
            }
        }
        2 => {
            let p_sq = inv(0);
            let m1_sq = mass(0);
            let m2_sq = mass(1);
            if family.numerator == Atom::num(1) {
                let m = [m1_sq, m2_sq];
                reduce_cayley(&modified_cayley(&m, &[p_sq]), &m, exponents)
            } else {
                let (c_b0, c_a1, c_a2) = bubble_numerator(
                    &family.numerator,
                    exponents[0],
                    exponents[1],
                    &p_sq,
                    &m1_sq,
                    &m2_sq,
                );
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
        }
        3 => {
            if family.numerator == Atom::num(1) {
                // Triangle invariants (s1,s2,s3)
                let masses: Vec<Atom> = (0..3).map(mass).collect();
                let lex = [inv(0), inv(1), inv(2)];
                reduce_cayley(&modified_cayley(&masses, &lex), &masses, exponents)
            } else {
                triangle_numerator(
                    &family.numerator,
                    exponents[0],
                    exponents[1],
                    exponents[2],
                    &inv(0),
                    &inv(1),
                    &inv(2),
                    &mass(0),
                    &mass(1),
                    &mass(2),
                )
            }
        }
        4 => {
            if family.numerator == Atom::num(1) {
                // Box invariants (p1,p2,p3,p4,s,t)
                let masses: Vec<Atom> = (0..4).map(mass).collect();
                let lex = [inv(0), inv(1), inv(2), inv(3), inv(4), inv(5)];
                reduce_cayley(&modified_cayley(&masses, &lex), &masses, exponents)
            } else {
                box_numerator(
                    &family.numerator,
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
                )
            }
        }
        n => {
            // N > 4: van Neerven-Vermaseren / FJT-Tarasov reduction to boxes.
            let expected = n * (n - 1) / 2;
            assert_eq!(
                family.kinematics.invariants.len(),
                expected,
                "an {n}-point family needs {expected} lexicographic pairwise invariants"
            );
            let masses: Vec<Atom> = (0..n).map(mass).collect();
            if family.numerator != Atom::num(1) {
                ngon_numerator(
                    &family.numerator,
                    exponents,
                    &family.kinematics.invariants,
                    &masses,
                )
            } else {
                let y = modified_cayley(&masses, &family.kinematics.invariants);
                reduce_cayley(&y, &masses, exponents)
            }
        }
    };

    Reduction { terms }
}

fn tadpole_coefficient(exponent: i32, m_sq: &Atom) -> Atom {
    if exponent <= 0 || *m_sq == Atom::Zero {
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
    if coeff == Atom::Zero {
        return;
    }
    match terms.iter_mut().find(|(_, m)| *m == master) {
        Some(slot) => slot.0 = &slot.0 + &coeff,
        None => terms.push((coeff, master)),
    }
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

// General symbolic determinant by cofactor expansion along the first row.
type Rp = RationalPolynomial<IntegerRing, u16>;

fn to_matrix(m: &[Vec<Atom>]) -> Matrix<RationalPolynomialField<IntegerRing, u16>> {
    let rows: Vec<Vec<Rp>> = m
        .iter()
        .map(|r| {
            r.iter()
                .map(|a| a.to_rational_polynomial::<_, _, u16>(&Q, &Z, None))
                .collect()
        })
        .collect();
    Matrix::from_nested_vec(rows, RationalPolynomialField::new(Z)).unwrap()
}

// Exact symbolic determinant
fn det(m: &[Vec<Atom>]) -> Atom {
    to_matrix(m).det().unwrap().to_expression()
}

// Inverse as an Atom matrix
fn matrix_inv(m: &[Vec<Atom>]) -> Vec<Vec<Atom>> {
    let inv = to_matrix(m).inv().unwrap();
    (0..inv.nrows())
        .map(|i| {
            (0..inv.ncols())
                .map(|j| inv[(i as u32, j as u32)].to_expression())
                .collect()
        })
        .collect()
}

// Modified Cayley matrix Y_ij = m_i^2 + m_j^2 - (r_i - r_j)^2 from the masses and the
// C(n,2) pairwise invariants (r_i - r_j)^2 in lexicographic (i<j) order.
fn modified_cayley(masses: &[Atom], pairwise: &[Atom]) -> Vec<Vec<Atom>> {
    let n = masses.len();
    let idx = |a: usize, b: usize| a * n - a * (a + 1) / 2 + (b - a - 1);
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        Atom::num(2) * &masses[i]
                    } else {
                        let (a, b) = (i.min(j), i.max(j));
                        &masses[i] + &masses[j] - &pairwise[idx(a, b)]
                    }
                })
                .collect()
        })
        .collect()
}

fn delete_row_col(m: &[Vec<Atom>], skip_row: usize, skip_col: usize) -> Vec<Vec<Atom>> {
    let n = m.len();
    let cols: Vec<usize> = (0..n).filter(|&c| c != skip_col).collect();
    (0..n)
        .filter(|&r| r != skip_row)
        .map(|r| cols.iter().map(|&c| m[r][c].clone()).collect())
        .collect()
}

// Rational null-space basis of `m`
fn nullspace(m: &[Vec<Atom>]) -> Vec<Vec<Atom>> {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    let mut a = m.to_vec();
    let mut pivot_of_col: Vec<Option<usize>> = vec![None; cols];
    let mut r = 0;
    for col in 0..cols {
        if r >= rows {
            break;
        }
        let Some(piv) = (r..rows).find(|&i| a[i][col] != Atom::Zero) else {
            continue;
        };
        a.swap(r, piv);
        let inv = Atom::num(1) / &a[r][col];
        for cell in a[r].iter_mut() {
            *cell = &*cell * &inv;
        }
        let pivot_row = a[r].clone();
        for (i, row) in a.iter_mut().enumerate() {
            if i != r && row[col] != Atom::Zero {
                let f = row[col].clone();
                for (cell, p) in row.iter_mut().zip(&pivot_row) {
                    *cell = &*cell - &(&f * p);
                }
            }
        }
        pivot_of_col[col] = Some(r);
        r += 1;
    }
    (0..cols)
        .filter(|&free| pivot_of_col[free].is_none())
        .map(|free| {
            let mut v = vec![Atom::Zero; cols];
            v[free] = Atom::num(1);
            for (col, &pr) in pivot_of_col.iter().enumerate() {
                if let Some(pr) = pr {
                    v[col] = -&a[pr][free];
                }
            }
            v
        })
        .collect()
}

// Degenerate reduction for det(Y) = 0
fn degenerate_coeffs(y: &[Vec<Atom>]) -> Vec<Atom> {
    let n = y.len();
    let mut yb = vec![vec![Atom::Zero; n + 1]; n + 1];
    for i in 0..n {
        for j in 0..n {
            yb[i][j] = y[i][j].clone();
        }
        yb[i][n] = Atom::num(1);
        yb[n][i] = Atom::num(1);
    }

    for v in nullspace(&yb) {
        if v[n] != Atom::Zero {
            return (0..n).map(|i| &v[i] / &v[n]).collect();
        }
    }
    panic!("degenerate Cayley reduction: exceptional kinematics");
}

// van Neerven-Vermaseren coefficients for I_N = sum_i c_i I_{N-1}^(i) (scalar, d=4):
// c_i = (YB^{-1})_{i0} / (YB^{-1})_{00}, YB the bordered Cayley matrix.
fn high_point_coeffs(y: &[Vec<Atom>]) -> Vec<Atom> {
    let det_y = det(y);
    if det_y == Atom::Zero {
        return degenerate_coeffs(y);
    }
    let n = y.len();
    let mut yb = vec![vec![Atom::Zero; n + 1]; n + 1];
    for i in 1..=n {
        yb[0][i] = Atom::num(1);
        yb[i][0] = Atom::num(1);
        for j in 1..=n {
            yb[i][j] = y[i - 1][j - 1].clone();
        }
    }
    (1..=n)
        .map(|col| {
            let minor = det(&delete_row_col(&yb, 0, col));
            let cofactor = if col % 2 == 0 { minor } else { -minor };
            cofactor / &det_y
        })
        .collect()
}

// Reduce an N-point (N>=4) from its Cayley matrix, masses, and powers.
fn emit_master(y: &[Vec<Atom>], masses: &[Atom]) -> Vec<(Atom, MasterIntegral)> {
    let pair = |a: usize, b: usize| (&masses[a] + &masses[b] - &y[a][b]).expand();
    let m = |i: usize| masses[i].clone();
    let master = match y.len() {
        2 => MasterIntegral::Bubble {
            p_sq: pair(0, 1),
            m1_sq: m(0),
            m2_sq: m(1),
        },
        3 => MasterIntegral::Triangle {
            p1_sq: pair(0, 1),
            p2_sq: pair(1, 2),
            p12_sq: pair(0, 2),
            m1_sq: m(0),
            m2_sq: m(1),
            m3_sq: m(2),
        },
        4 => MasterIntegral::Box {
            p1_sq: pair(0, 1),
            p2_sq: pair(1, 2),
            p3_sq: pair(2, 3),
            p4_sq: pair(0, 3),
            s: pair(0, 2),
            t: pair(1, 3),
            m1_sq: m(0),
            m2_sq: m(1),
            m3_sq: m(2),
            m4_sq: m(3),
        },
        _ => unreachable!("emit_master: scalar leaf must be a bubble/triangle/box"),
    };
    vec![(Atom::num(1), master)]
}

fn reduce_cayley(
    y: &[Vec<Atom>],
    masses: &[Atom],
    exponents: &[i32],
) -> Vec<(Atom, MasterIntegral)> {
    let drop = |i: usize, v: &[Atom]| -> Vec<Atom> {
        v.iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, x)| x.clone())
            .collect()
    };
    let n = y.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        let c = tadpole_coefficient(exponents[0], &masses[0]);
        return if c == Atom::Zero {
            Vec::new()
        } else {
            vec![(
                c,
                MasterIntegral::Tadpole {
                    m_sq: masses[0].clone(),
                },
            )]
        };
    }
    if exponents.contains(&0) {
        let z = exponents.iter().position(|&e| e == 0).unwrap();
        let sub_exp: Vec<i32> = exponents
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != z)
            .map(|(_, &e)| e)
            .collect();
        return reduce_cayley(&delete_row_col(y, z, z), &drop(z, masses), &sub_exp);
    }
    if exponents.iter().all(|&e| e == 1) {
        if n <= 4 {
            return emit_master(y, masses);
        }
        let ones = vec![1; n - 1];
        let mut terms = Vec::new();
        for (i, c_i) in high_point_coeffs(y).iter().enumerate() {
            add_scaled(
                &mut terms,
                c_i,
                reduce_cayley(&delete_row_col(y, i, i), &drop(i, masses), &ones),
            );
        }
        return terms;
    }

    let det_y = det(y);
    if det_y == Atom::Zero {
        let mut terms = Vec::new();
        for (l, cl) in degenerate_coeffs(y).iter().enumerate() {
            if exponents[l] == 0 {
                continue;
            }
            let mut a2 = exponents.to_vec();
            a2[l] -= 1;
            add_scaled(&mut terms, cl, reduce_cayley(y, masses, &a2));
        }
        return terms;
    }

    // dotted: FJT/Tarasov index-lowering
    let mut k = 0;
    for i in 1..n {
        if exponents[i] > exponents[k] {
            k = i;
        }
    }
    let mut a = exponents.to_vec();
    a[k] -= 1;
    let total: i32 = a.iter().sum();
    // inv[k][i] = adj[k][i]/det(y)
    let inv = matrix_inv(y);
    let den = Atom::num(i64::from(a[k]));

    let d = Atom::var(S.d);
    let mut diag = Atom::Zero;
    for i in 0..n {
        let factor = &d - Atom::num(i64::from(total + a[i]));
        diag += &inv[k][i] * &factor;
    }
    let mut acc = Vec::new();
    add_scaled(&mut acc, &(diag / &den), reduce_cayley(y, masses, &a));
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let num = Atom::num(i64::from(-a[j])) * &inv[k][i];
            let mut ch = a.clone();
            ch[j] += 1;
            ch[i] -= 1;
            add_scaled(&mut acc, &(num / &den), reduce_cayley(y, masses, &ch));
        }
    }
    acc
}

type DotMono = (Vec<i32>, Atom);

fn dot_ll() -> Atom {
    function!(S.dot, Atom::var(S.k), Atom::var(S.k))
}
fn dot_lq(j: usize) -> Atom {
    let q = symbol!(format!("oneloop::q{}", j + 1));
    function!(S.dot, Atom::var(S.k), Atom::var(q))
}

fn numerator_to_monos(numerator: &Atom, n_ext: usize) -> Vec<DotMono> {
    let xll = symbol!("oneloop::xll");
    let mut vars = vec![xll];
    let mut n = numerator
        .replace(dot_ll().to_pattern())
        .with(Atom::var(xll));
    for a in 0..n_ext {
        let xq = symbol!(format!("oneloop::xq{}", a + 1));
        vars.push(xq);
        n = n.replace(dot_lq(a).to_pattern()).with(Atom::var(xq));
    }
    extract_monomials(&n, &vars)
}

// Monomials of a polynomial in `vars`
fn extract_monomials(expr: &Atom, vars: &[Symbol]) -> Vec<DotMono> {
    fn rec(expr: &Atom, vars: &[Symbol], idx: usize, exps: &mut Vec<i32>, out: &mut Vec<DotMono>) {
        if idx == vars.len() {
            if *expr != Atom::Zero {
                out.push((exps.clone(), expr.clone()));
            }
            return;
        }
        let v = vars[idx];
        let vatom = Atom::var(v);
        // Peel off powers of `v` one at a time via Taylor coefficients:
        //   coeff of v^c = (1/c!) d^c/dv^c expr |_{v=0}
        let mut deriv = expr.clone();
        let mut factorial: i64 = 1;
        let mut c = 0i32;
        loop {
            // coefficient of v^c is (deriv with v->0) / factorial
            let at_zero = deriv.replace(vatom.to_pattern()).with(Atom::Zero);
            if at_zero != Atom::Zero {
                let coeff = at_zero / Atom::num(factorial);
                exps[idx] = c;
                rec(&coeff, vars, idx + 1, exps, out);
            }
            // step to next power
            let next = deriv.derivative(v);
            if next == Atom::Zero {
                break;
            }
            deriv = next;
            c += 1;
            assert!(
                c <= 20,
                "numerator degree in one variable exceeds the supported bound (20)"
            );
            factorial *= i64::from(c);
        }
        exps[idx] = 0;
    }
    let mut out = Vec::new();
    let mut exps = vec![0i32; vars.len()];
    rec(expr, vars, 0, &mut exps, &mut out);
    out
}

// A direction in external-momentum space: integer coefficients over (q1,...,q_{n_ext}).
type Dir = Vec<i32>;

// Base Gram closure: (a,b) in {0,1,2} -> q_a . q_b at a kinematic point.
type GramFn = Box<dyn Fn(usize, usize) -> Atom>;
// A reduction result: a list of (coefficient, master integral).
type Masters = Vec<(Atom, MasterIntegral)>;
// Scalar reducer for a topology: propagator exponents -> masters.
type ScalarFn = Box<dyn Fn(&[i32]) -> Masters>;
// Pinch routine: (kept lines, their exponents, residual dot-polynomial) -> masters.
type PinchFn = Box<dyn Fn(&[usize], &[i32], &[DotMono]) -> Masters>;

// A topology level for arbitrary-rank numerator reduction.
struct Topo {
    // Reducible directions (RSP span), as coefficient vectors over (q1,q2,q3).
    // dot(k, d) for d in this span is reducible; everything else is an ISP.
    rsp_dirs: Vec<Dir>,
    // Base Gram q_a . q_b for a,b in {0,1,2} (q1,q2,q3) at this kinematic point.
    base_gram: GramFn,
    // RSP rule for dot(k,k): linear in den symbols.
    rule_ll: Atom,
    // RSP rule for dot(k, d) for each reducible direction d (same order as rsp_dirs).
    rule_lq: Vec<Atom>,
    // den symbols, one per surviving line, in order.
    den_syms: Vec<Symbol>,
    // Propagator exponents a_i of the surviving lines.
    a: Vec<i32>,
    // r_i: shift of each surviving line as a direction over (q1,q2,q3).
    r_coeffs: Vec<Dir>,
    // Mass^2 of each surviving line.
    masses: Vec<Atom>,
    // Reduce a pure propagator-shift (all exponents) of this topology to masters.
    scalar: ScalarFn,
    // Reduce a residual numerator on the sub-topology that KEEPS the given lines/exponents.
    pinch: PinchFn,
}

impl Topo {
    // Gram dot of two directions u, v over the external momenta: sum u_a v_b (q_a.q_b).
    fn dir_dot(&self, u: &Dir, v: &Dir) -> Atom {
        dir_dot_with(&*self.base_gram, u, v)
    }
}

// Solve the symbolic Gram system G * alpha = rhs (G_{ij} = rsp_dirs[i].rsp_dirs[j]) by
// Cramer's rule with the exact symbolic determinant, so det cancels downstream.
fn gram_solve(topo: &Topo, rhs: &[Atom]) -> Vec<Atom> {
    let dirs = &topo.rsp_dirs;
    let n = dirs.len();
    if n == 0 {
        return vec![];
    }
    let g: Vec<Vec<Atom>> = (0..n)
        .map(|i| (0..n).map(|j| topo.dir_dot(&dirs[i], &dirs[j])).collect())
        .collect();
    // Solve G c = rhs
    let b: Vec<Vec<Atom>> = rhs.iter().map(|x| vec![x.clone()]).collect();
    let sol = to_matrix(&g).solve(&to_matrix(&b)).unwrap();
    (0..n).map(|i| sol[(i as u32, 0)].to_expression()).collect()
}

// All perfect matchings (pairings) of 0..m (m even).  Each pairing is a Vec of
// (i,j) index pairs.
fn pairings(items: &[usize]) -> Vec<Vec<(usize, usize)>> {
    if items.is_empty() {
        return vec![vec![]];
    }
    let first = items[0];
    let mut out = Vec::new();
    for i in 1..items.len() {
        let mut rest = Vec::with_capacity(items.len() - 2);
        for (j, &it) in items.iter().enumerate() {
            if j != 0 && j != i {
                rest.push(it);
            }
        }
        for mut sub in pairings(&rest) {
            sub.insert(0, (first, items[i]));
            out.push(sub);
        }
    }
    out
}

// Step 1: PV transverse average.  Split each loop dot into a reducible part
// (par = sum alpha_ab dot(k,q_b)) plus an ISP part xi; average even xi products over
// pairings of the transverse metric in n_t = d-|rsp| dims.
#[allow(clippy::assign_op_pattern)]
fn isp_project(topo: &Topo, monos: &[DotMono]) -> Vec<DotMono> {
    let d = Atom::var(S.d);
    let n_rsp = topo.rsp_dirs.len();
    let n_out = n_rsp + 1; // output DotMono length: l.l slot + one per reducible direction
    let n_ext = topo
        .rsp_dirs
        .first()
        .map(Vec::len)
        .or_else(|| monos.first().map(|m| m.0.len() - 1))
        .unwrap_or(0); // input external dimension
    let nt = &d - Atom::num(n_rsp as i64);
    let unit = |a: usize| -> Dir {
        let mut e = vec![0i32; n_ext];
        e[a] = 1;
        e
    };

    // Reducible-span coeffs of v: exact integer decomposition (finite even for a singular
    // Gram), else the oblique Gram inverse (only reached for genuine ISPs on full-rank subs).
    let alpha_of = |v: &Dir| -> Vec<Atom> {
        if let Some(c) = span_coeffs(&topo.rsp_dirs, v) {
            return c;
        }
        let rhs: Vec<Atom> = topo.rsp_dirs.iter().map(|r| topo.dir_dot(r, v)).collect();
        gram_solve(topo, &rhs)
    };
    // Transverse metric (u.v)_perp = u.v - (rsp.u)^T G^{-1} (rsp.v).
    let perp_metric = |u: &Dir, v: &Dir| -> Atom {
        let alpha = alpha_of(u);
        let mut par = Atom::Zero;
        for (i, r) in topo.rsp_dirs.iter().enumerate() {
            par += &alpha[i] * topo.dir_dot(r, v);
        }
        topo.dir_dot(u, v) - par
    };
    // par(v): coefficients on the reducible loop dots dot(k, rsp_dirs[j]).
    let par_of = |v: &Dir| -> Vec<Atom> { alpha_of(v) };
    // A unit direction is reducible iff it lies in the integer span of rsp_dirs.
    let is_reducible = |v: &Dir| -> bool { dir_in_span(&topo.rsp_dirs, v) };

    // Only build the transverse metric (needs the Gram inverse) when a real ISP is present.
    let any_isp = monos
        .iter()
        .any(|(exps, _)| (0..n_ext).any(|a| exps[a + 1] > 0 && !is_reducible(&unit(a))));

    // L2perp = dot(k,k) - sum_{i,j} dot(k,rsp_i) Ginv_ij dot(k,rsp_j), keyed in the
    // OUTPUT convention (slot j+1 = dot(k, rsp_dirs[j])).
    let l2perp: Vec<DotMono> = if !any_isp {
        Vec::new()
    } else {
        let mut ll_slot = vec![0i32; n_out];
        ll_slot[0] = 1;
        let mut out: Vec<DotMono> = vec![(ll_slot, Atom::num(1))];
        for i in 0..n_rsp {
            let mut e = vec![Atom::Zero; n_rsp];
            e[i] = Atom::num(1);
            let ginv_col = gram_solve(topo, &e); // G^{-1} column i
            for (j, gij) in ginv_col.iter().enumerate() {
                let mut ek = vec![0i32; n_out];
                ek[i + 1] += 1;
                ek[j + 1] += 1;
                out.push((ek, Atom::num(-1) * gij));
            }
        }
        out
    };

    let mut result: Vec<DotMono> = Vec::new();
    for (exps, coeff) in monos {
        let ll_pow = exps[0];
        // factors: the multiset of unit external directions from the q-dots.
        let mut factors: Vec<Dir> = Vec::new();
        for a in 0..n_ext {
            for _ in 0..exps[a + 1] {
                factors.push(unit(a));
            }
        }
        // base poly: ll^{ll_pow} (in output convention) times coeff.
        let mut base_e = vec![0i32; n_out];
        base_e[0] = ll_pow;
        let base_mono: Vec<DotMono> = vec![(base_e, coeff.clone())];
        // term state: (par-poly so far, xi-multiset of ISP directions).
        let mut terms: Vec<(Vec<DotMono>, Vec<Dir>)> = vec![(base_mono, Vec::new())];
        for fdir in &factors {
            let reducible = is_reducible(fdir);
            let par = par_of(fdir);
            let mut next: Vec<(Vec<DotMono>, Vec<Dir>)> = Vec::new();
            for (poly, xis) in &terms {
                // par branch
                let mut par_poly: Vec<DotMono> = Vec::new();
                for (j, c) in par.iter().enumerate() {
                    let mut ek = vec![0i32; n_out];
                    ek[j + 1] += 1;
                    par_poly.push((ek, c.clone()));
                }
                next.push((mono_mul(poly, &par_poly), xis.clone()));
                // xi branch (only if the direction is a genuine ISP)
                if !reducible {
                    let mut new_xis = xis.clone();
                    new_xis.push(fdir.clone());
                    next.push((poly.clone(), new_xis));
                }
            }
            terms = next;
        }
        for (poly, xis) in terms {
            let m = xis.len();
            if m % 2 == 1 {
                continue;
            }
            let kk = m / 2;
            let idxs: Vec<usize> = (0..m).collect();
            let mut s = Atom::Zero;
            for pr in pairings(&idxs) {
                let mut term = Atom::num(1);
                for (i, j) in pr {
                    term *= perp_metric(&xis[i], &xis[j]);
                }
                s += term;
            }
            if s == Atom::Zero {
                continue;
            }
            let mut norm = Atom::num(1);
            for jj in 0..kk {
                norm = norm * (&nt + Atom::num(2 * jj as i64));
            }
            let avg_coeff = s / &norm;
            let mut piece = poly.clone();
            for _ in 0..kk {
                piece = mono_mul(&piece, &l2perp);
            }
            for (e, c) in &piece {
                push_mono(&mut result, e.clone(), &avg_coeff * c);
            }
        }
    }
    result
}

fn dir_in_span(dirs: &[Dir], v: &Dir) -> bool {
    let to_i64 = |d: &Dir| -> Vec<i64> { d.iter().map(|&x| i64::from(x)).collect() };
    let mut basis: Vec<Vec<i64>> = Vec::new();
    for d in dirs {
        let mut row = to_i64(d);
        reduce_against(&mut row, &basis);
        if row.iter().any(|&x| x != 0) {
            basis.push(row);
        }
    }
    let mut row = to_i64(v);
    reduce_against(&mut row, &basis);
    row.iter().all(|&x| x == 0)
}

// Rational c with v = sum c_i dirs[i] (None if v is not in the coordinate span), by exact
// Gauss-Jordan.
fn span_coeffs(dirs: &[Dir], v: &Dir) -> Option<Vec<Atom>> {
    let n_ext = v.len();
    let n = dirs.len();
    if n == 0 {
        return v.iter().all(|&x| x == 0).then(Vec::new);
    }
    // Augmented rows = coordinates (n_ext), cols = n unknowns then rhs.
    let mut a: Vec<Vec<Atom>> = (0..n_ext)
        .map(|row| {
            let mut r: Vec<Atom> = dirs.iter().map(|d| Atom::num(i64::from(d[row]))).collect();
            r.push(Atom::num(i64::from(v[row])));
            r
        })
        .collect();
    let mut pivot_col = vec![usize::MAX; n_ext];
    let mut r = 0;
    for col in 0..n {
        let Some(piv) = (r..n_ext).find(|&i| a[i][col] != Atom::Zero) else {
            continue;
        };
        a.swap(r, piv);
        let pivot_row = a[r].clone();
        for (i, row) in a.iter_mut().enumerate() {
            if i != r && row[col] != Atom::Zero {
                let factor = &row[col] / &pivot_row[col];
                for (cell, p) in row.iter_mut().zip(&pivot_row) {
                    *cell = &*cell - &(&factor * p);
                }
            }
        }
        pivot_col[r] = col;
        r += 1;
        if r == n_ext {
            break;
        }
    }
    // Inconsistent row (all unknowns zero, rhs nonzero) => v not in the span.
    for row in &a {
        if row[..n].iter().all(|x| *x == Atom::Zero) && row[n] != Atom::Zero {
            return None;
        }
    }
    // Free (non-pivot) unknowns are set to zero; read each pivot off its row.
    let mut c = vec![Atom::Zero; n];
    for (row, &col) in a.iter().zip(&pivot_col) {
        if col != usize::MAX {
            c[col] = &row[n] / &row[col];
        }
    }
    Some(c)
}

// Clear each basis pivot from `row` via row <- b_piv*row - row_piv*b.
fn reduce_against(row: &mut [i64], basis: &[Vec<i64>]) {
    for b in basis {
        let piv = b.iter().position(|&x| x != 0).unwrap();
        if row[piv] != 0 {
            let (rp, bp) = (row[piv], b[piv]);
            for (rc, &bc) in row.iter_mut().zip(b) {
                *rc = bp * *rc - rp * bc;
            }
        }
    }
}

// Multiply two dot-polynomials.
fn mono_mul(a: &[DotMono], b: &[DotMono]) -> Vec<DotMono> {
    let mut out: Vec<DotMono> = Vec::new();
    for (ea, ca) in a {
        for (eb, cb) in b {
            let e: Vec<i32> = ea.iter().zip(eb).map(|(x, y)| x + y).collect();
            push_mono(&mut out, e, ca * cb);
        }
    }
    out
}

fn push_mono(acc: &mut Vec<DotMono>, e: Vec<i32>, c: Atom) {
    if c == Atom::Zero {
        return;
    }
    match acc.iter_mut().find(|(ee, _)| *ee == e) {
        Some(slot) => slot.1 = &slot.1 + &c,
        None => acc.push((e, c)),
    }
}

// Gram dot of two directions using a base-Gram closure over the external momenta.
#[allow(clippy::needless_range_loop)]
fn dir_dot_with(gram: &dyn Fn(usize, usize) -> Atom, u: &Dir, v: &Dir) -> Atom {
    let mut out = Atom::Zero;
    for a in 0..u.len() {
        for b in 0..v.len() {
            if u[a] != 0 && v[b] != 0 {
                out += Atom::num(i64::from(u[a]) * i64::from(v[b])) * gram(a, b);
            }
        }
    }
    out
}

// Clone a base-Gram closure into a fresh boxed closure
fn clone_gram(
    gram: &dyn Fn(usize, usize) -> Atom,
    n_ext: usize,
) -> Box<dyn Fn(usize, usize) -> Atom> {
    let entries: Vec<Vec<Atom>> = (0..n_ext)
        .map(|a| (0..n_ext).map(|b| gram(a, b)).collect())
        .collect();
    Box::new(move |a: usize, b: usize| entries[a][b].clone())
}

// Apply the loop-momentum shift l -> l - r to a dot-polynomial numerator.
fn shift_numerator(
    monos: &[DotMono],
    r: &Dir,
    gram: &dyn Fn(usize, usize) -> Atom,
) -> Vec<DotMono> {
    if r.iter().all(|&x| x == 0) {
        return monos.to_vec();
    }
    let n_ext = r.len();
    let n_dots = n_ext + 1;
    // dot(l,r) = sum_a r_a dot(l,q_a)
    let ll_shift: Vec<DotMono> = {
        let mut ll = vec![0i32; n_dots];
        ll[0] = 1;
        let mut v: Vec<DotMono> = vec![(ll, Atom::num(1))]; // dot(l,l)
        for a in 0..n_ext {
            if r[a] != 0 {
                let mut e = vec![0i32; n_dots];
                e[a + 1] = 1;
                push_mono(&mut v, e, Atom::num(-2 * i64::from(r[a])));
            }
        }
        push_mono(&mut v, vec![0i32; n_dots], dir_dot_with(gram, r, r));
        v
    };
    let lq_shift = |a: usize| -> Vec<DotMono> {
        let mut e = vec![0i32; n_dots];
        e[a + 1] = 1;
        let mut ra: Dir = vec![0i32; n_ext];
        ra[a] = 1;
        vec![
            (e, Atom::num(1)),
            (
                vec![0i32; n_dots],
                Atom::num(-1) * dir_dot_with(gram, r, &ra),
            ),
        ]
    };
    let mut out: Vec<DotMono> = Vec::new();
    for (e, c) in monos {
        // build the shifted product: ll_shift^{e0} * prod_a lq_shift(a)^{e_{a+1}}
        let mut prod: Vec<DotMono> = vec![(vec![0i32; n_dots], c.clone())];
        for _ in 0..e[0] {
            prod = mono_mul(&prod, &ll_shift);
        }
        for a in 0..n_ext {
            for _ in 0..e[a + 1] {
                prod = mono_mul(&prod, &lq_shift(a));
            }
        }
        for (pe, pc) in prod {
            push_mono(&mut out, pe, pc);
        }
    }
    out
}

// Re-express the inverse propagator D_i = dot(k,k) + 2 dot(k,r_i) + (r_i.r_i - m_i) as a dot-polynomial.
fn di_as_dots(r_coeff: &[i32], mass: &Atom, gram: &dyn Fn(usize, usize) -> Atom) -> Vec<DotMono> {
    let n_ext = r_coeff.len();
    let n_dots = n_ext + 1;
    let mut ll = vec![0i32; n_dots];
    ll[0] = 1;
    let mut out: Vec<DotMono> = vec![(ll, Atom::num(1))]; // dot(k,k)
    for j in 0..n_ext {
        if r_coeff[j] != 0 {
            let mut e = vec![0i32; n_dots];
            e[j + 1] = 1;
            push_mono(&mut out, e, Atom::num(2 * i64::from(r_coeff[j])));
        }
    }
    // r_i . r_i = sum_{a,b} c_a c_b (q_a.q_b)
    let mut rr = Atom::Zero;
    for a in 0..n_ext {
        for b in 0..n_ext {
            if r_coeff[a] != 0 && r_coeff[b] != 0 {
                rr += Atom::num(i64::from(r_coeff[a]) * i64::from(r_coeff[b])) * gram(a, b);
            }
        }
    }
    push_mono(&mut out, vec![0i32; n_dots], rr - mass);
    out
}

// The generic engine: reduce a dot-polynomial numerator on topology `topo`.
#[allow(clippy::needless_range_loop)]
fn reduce_num(topo: &Topo, numerator_monos: &[DotMono]) -> Vec<(Atom, MasterIntegral)> {
    // Step 1: ISP projection (no-op for the top topology).
    let projected = isp_project(topo, numerator_monos);
    if projected.is_empty() {
        return Vec::new();
    }

    // Step 2: RSP substitution (slot 0 -> rule_ll, slot j+1 -> rule_lq[j]) -> polynomial in den symbols.
    let n_lines = topo.den_syms.len();
    let mut n = Atom::Zero;
    for (e, c) in &projected {
        let mut term = c.clone();
        for _ in 0..e[0] {
            term = &term * &topo.rule_ll;
        }
        for (j, rule) in topo.rule_lq.iter().enumerate() {
            for _ in 0..e[j + 1] {
                term = &term * rule;
            }
        }
        n += term;
    }

    // Step 3: extract den-monomials.
    let den_monos = extract_monomials(&n, &topo.den_syms);

    // Step 4: route each monomial.
    let mut acc: Vec<(Atom, MasterIntegral)> = Vec::new();
    for (cexp, coeff) in den_monos {
        let mut b = topo.a.clone();
        for i in 0..n_lines {
            b[i] -= cexp[i];
        }
        if b.iter().all(|&x| x >= 0) {
            // pure propagator shift
            add_scaled(&mut acc, &coeff, (topo.scalar)(&b));
        } else {
            // pinch lines with b_i <= 0; residual prod D_i^{|b_i|} on survivors.
            let keep: Vec<usize> = (0..n_lines).filter(|&i| b[i] > 0).collect();
            let keep_exps: Vec<i32> = keep.iter().map(|&i| b[i]).collect();
            // build residual dot-polynomial (in this topology's external basis)
            let n_ext = topo.r_coeffs.first().map_or(0, |r| r.len());
            let mut residual: Vec<DotMono> = vec![(vec![0i32; n_ext + 1], Atom::num(1))];
            for i in 0..n_lines {
                if b[i] < 0 {
                    let di = di_as_dots(&topo.r_coeffs[i], &topo.masses[i], &*topo.base_gram);
                    for _ in 0..(-b[i]) {
                        residual = mono_mul(&residual, &di);
                    }
                }
            }
            add_scaled(&mut acc, &coeff, (topo.pinch)(&keep, &keep_exps, &residual));
        }
    }
    acc
}

// --- Topology builders -----------------------------------------------------

fn den_symbol(i: usize) -> Symbol {
    symbol!(format!("oneloop::den{}", i + 1))
}

// Reduce a tadpole numerator
fn tadpole_numerator(numerator: &Atom, a1: i32, m_sq: &Atom) -> Vec<(Atom, MasterIntegral)> {
    let stripped = numerator
        .replace(dot_ll().to_pattern())
        .with(Atom::var(symbol!("oneloop::tad_ll")));
    let has_k = stripped
        .replace(Atom::var(S.k).to_pattern())
        .with(Atom::Zero)
        != stripped;
    assert!(
        !has_k,
        "tadpole numerator must be a polynomial in dot(k,k) (no external momenta at one line)"
    );
    let z = Atom::Zero;
    reduce_num(
        &tadpole_topo(a1, m_sq, base_gram_box(&z, &z, &z, &z, &z, &z), 0),
        &numerator_to_monos(numerator, 0),
    )
}

// Reduce a bubble with an arbitrary-rank numerator.
fn bubble_numerator(
    numerator: &Atom,
    a1: i32,
    a2: i32,
    p_sq: &Atom,
    m1_sq: &Atom,
    m2_sq: &Atom,
) -> (Atom, Atom, Atom) {
    let terms = reduce_num(
        &bubble_topo(a1, a2, p_sq, m1_sq, m2_sq),
        &numerator_to_monos(numerator, 3),
    );
    let mut c_b0 = Atom::Zero;
    let mut c_a1 = Atom::Zero;
    let mut c_a2 = Atom::Zero;
    for (c, m) in terms {
        match m {
            MasterIntegral::Bubble { .. } => c_b0 = &c_b0 + &c,
            MasterIntegral::Tadpole { m_sq } if m_sq == *m1_sq => c_a1 = &c_a1 + &c,
            MasterIntegral::Tadpole { m_sq } if m_sq == *m2_sq => c_a2 = &c_a2 + &c,
            _ => unreachable!("bubble numerator produced an unexpected master"),
        }
    }
    (c_b0, c_a1, c_a2)
}

fn bubble_topo(a1: i32, a2: i32, p_sq: &Atom, m1_sq: &Atom, m2_sq: &Atom) -> Topo {
    // Top-level bubble: line 1 has zero shift, line 2 has shift q1.
    bubble_topo_general(
        a1,
        a2,
        vec![0, 0, 0],
        vec![1, 0, 0],
        m1_sq,
        m2_sq,
        p_sq.clone(),
        base_gram_box(
            p_sq,
            &Atom::Zero,
            &Atom::Zero,
            &Atom::Zero,
            &Atom::Zero,
            &Atom::Zero,
        ),
    )
}

// A general bubble sub-topology: lines r_a (reference) and r_b, with the parent Gram for ISP projection.
#[allow(clippy::too_many_arguments)]
fn bubble_topo_general(
    a1: i32,
    a2: i32,
    r_a: Dir,
    r_b: Dir,
    m_a: &Atom,
    m_b: &Atom,
    p_sq: Atom,
    base_gram: Box<dyn Fn(usize, usize) -> Atom>,
) -> Topo {
    let den1 = den_symbol(0);
    let den2 = den_symbol(1);
    let d1 = Atom::var(den1);
    let d2 = Atom::var(den2);
    let two = Atom::num(2);
    // Reducible direction w = r_b - r_a.
    let w: Dir = r_b.iter().zip(&r_a).map(|(b, a)| b - a).collect();
    // In the reference frame:
    //   l.l = D1 + m_a ;  l.w = (D2 - D1 - m_a + m_b - p_sq)/2
    let rule_ll = &d1 + m_a;
    let rule_w = (&d2 - &d1 - m_a + m_b - &p_sq) / &two;
    let (p1, m1a, m2a) = (p_sq.clone(), m_a.clone(), m_b.clone());

    let pinch_gram = clone_gram(&*base_gram, r_a.len());
    let pinch_r = vec![r_a.clone(), r_b.clone()];
    let pinch_masses = [m_a.clone(), m_b.clone()];
    Topo {
        rsp_dirs: vec![w],
        base_gram,
        rule_ll,
        rule_lq: vec![rule_w],
        den_syms: vec![den1, den2],
        a: vec![a1, a2],
        r_coeffs: vec![r_a, r_b],
        masses: vec![m_a.clone(), m_b.clone()],
        scalar: Box::new(move |b| {
            let m = [m1a.clone(), m2a.clone()];
            reduce_cayley(&modified_cayley(&m, std::slice::from_ref(&p1)), &m, b)
        }),
        // A bubble can pinch to a tadpole carrying a residual; route generically.
        pinch: Box::new(move |keep, exps, num| {
            pinch_to_subtopo(keep, exps, num, &pinch_r, &pinch_masses, &*pinch_gram)
        }),
    }
}

// Base Gram over (q1,q2,q3) from box invariants (triangle: p3=p4=t=0; bubble: only q1.q1).
fn base_gram_box(
    p1: &Atom,
    p2: &Atom,
    p3: &Atom,
    p4: &Atom,
    s: &Atom,
    t: &Atom,
) -> Box<dyn Fn(usize, usize) -> Atom> {
    let two = Atom::num(2);
    let q12 = (s - p1 - p2) / &two;
    let q23 = (t - p2 - p3) / &two;
    let q13 = (p4 - p1 - p2 - p3) / &two - &q12 - &q23;
    let (p1, p2, p3) = (p1.clone(), p2.clone(), p3.clone());
    Box::new(move |a: usize, b: usize| -> Atom {
        let (i, j) = if a <= b { (a, b) } else { (b, a) };
        match (i, j) {
            (0, 0) => p1.clone(),
            (1, 1) => p2.clone(),
            (2, 2) => p3.clone(),
            (0, 1) => q12.clone(),
            (1, 2) => q23.clone(),
            (0, 2) => q13.clone(),
            _ => unreachable!(),
        }
    })
}

// Base Gram q_a.q_b over n_ext external momenta from the C(n_ext+1,2)
fn base_gram_from_pairwise(n_ext: usize, pairwise: &[Atom]) -> Box<dyn Fn(usize, usize) -> Atom> {
    let n = n_ext + 1;
    let idx = |i: usize, j: usize| i * n - i * (i + 1) / 2 + (j - i - 1);
    let s = |i: usize, j: usize| -> Atom {
        if i == j {
            Atom::Zero
        } else {
            let (lo, hi) = (i.min(j), i.max(j));
            pairwise[idx(lo, hi)].clone()
        }
    };
    let gram: Vec<Vec<Atom>> = (0..n_ext)
        .map(|a| {
            (0..n_ext)
                .map(|b| {
                    if a == b {
                        s(a, a + 1)
                    } else {
                        (s(a + 1, b) + s(a, b + 1) - s(a + 1, b + 1) - s(a, b)) / Atom::num(2)
                    }
                })
                .collect()
        })
        .collect();
    Box::new(move |a: usize, b: usize| gram[a][b].clone())
}

// Reduce a triangle with an arbitrary-rank numerator.
#[allow(clippy::too_many_arguments)]
fn triangle_numerator(
    numerator: &Atom,
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
    reduce_num(
        &triangle_topo(a1, a2, a3, s1, s2, s3, m1, m2, m3),
        &numerator_to_monos(numerator, 3),
    )
}

#[allow(clippy::too_many_arguments)]
fn triangle_topo(
    a1: i32,
    a2: i32,
    a3: i32,
    s1: &Atom,
    s2: &Atom,
    s3: &Atom,
    m1: &Atom,
    m2: &Atom,
    m3: &Atom,
) -> Topo {
    triangle_topo_general(
        a1,
        a2,
        a3,
        [vec![0, 0, 0], vec![1, 0, 0], vec![1, 1, 0]],
        [m1.clone(), m2.clone(), m3.clone()],
        [s1.clone(), s2.clone(), s3.clone()],
        base_gram_box(s1, s2, &Atom::Zero, s3, s3, &Atom::Zero),
    )
}

// A general triangle sub-topology: lines r[0..3], invariants s, parent Gram for ISP projection.
#[allow(clippy::too_many_arguments)]
fn triangle_topo_general(
    a1: i32,
    a2: i32,
    a3: i32,
    r: [Dir; 3],
    masses: [Atom; 3],
    s: [Atom; 3],
    base_gram: Box<dyn Fn(usize, usize) -> Atom>,
) -> Topo {
    let den1 = den_symbol(0);
    let den2 = den_symbol(1);
    let den3 = den_symbol(2);
    let d1 = Atom::var(den1);
    let d2 = Atom::var(den2);
    let d3 = Atom::var(den3);
    let two = Atom::num(2);
    let (m1, m2, m3) = (masses[0].clone(), masses[1].clone(), masses[2].clone());
    let (s1, s2, s3) = (s[0].clone(), s[1].clone(), s[2].clone());
    // Reducible directions w1 = r[1]-r[0], w2 = r[2]-r[1]; RSP rules below.
    let w1: Dir = r[1].iter().zip(&r[0]).map(|(b, a)| b - a).collect();
    let w2: Dir = r[2].iter().zip(&r[1]).map(|(b, a)| b - a).collect();
    let rule_ll = &d1 + &m1;
    let rule_w1 = (&d2 - &d1 - &m1 + &m2 - &s1) / &two;
    let rule_w2 = (&d3 - &d2 - &m2 + &m3 + &s1 - &s3) / &two;
    let (ss1, ss2, ss3) = (s1.clone(), s2.clone(), s3.clone());
    let (sm1, sm2, sm3) = (m1.clone(), m2.clone(), m3.clone());
    let r_for_pinch = r.clone();
    let masses_for_pinch = [m1.clone(), m2.clone(), m3.clone()];
    // Pinch-to-bubble reuses this triangle's full Gram
    let gram_for_pinch = clone_gram(&*base_gram, r[0].len());
    Topo {
        rsp_dirs: vec![w1, w2],
        base_gram,
        rule_ll,
        rule_lq: vec![rule_w1, rule_w2],
        den_syms: vec![den1, den2, den3],
        a: vec![a1, a2, a3],
        r_coeffs: r.to_vec(),
        masses: vec![m1, m2, m3],
        scalar: Box::new(move |b| {
            let m = [sm1.clone(), sm2.clone(), sm3.clone()];
            let lex = [ss1.clone(), ss3.clone(), ss2.clone()];
            reduce_cayley(&modified_cayley(&m, &lex), &m, b)
        }),
        pinch: Box::new(move |keep, exps, num| {
            pinch_to_subtopo(
                keep,
                exps,
                num,
                &r_for_pinch,
                &masses_for_pinch,
                &*gram_for_pinch,
            )
        }),
    }
}

fn pinch_to_subtopo(
    keep: &[usize],
    exps: &[i32],
    num: &[DotMono],
    all_r: &[Dir],
    all_masses: &[Atom],
    base_gram: &dyn Fn(usize, usize) -> Atom,
) -> Vec<(Atom, MasterIntegral)> {
    if keep.is_empty() {
        // No denominator survives
        return Vec::new();
    }
    let r_ref = all_r[keep[0]].clone();
    let n_ext = r_ref.len();
    // Shift the residual into the child's reference frame
    let shifted = shift_numerator(num, &r_ref, base_gram);
    // Surviving line shift-directions relative to the reference, and masses.
    let dirs: Vec<Dir> = keep
        .iter()
        .map(|&i| all_r[i].iter().zip(&r_ref).map(|(x, y)| x - y).collect())
        .collect();
    let kmasses: Vec<Atom> = keep.iter().map(|&i| all_masses[i].clone()).collect();

    match keep.len() {
        k if k >= 5 => {
            let topo = ngon_topo(exps, &dirs, &kmasses, clone_gram(base_gram, n_ext));
            reduce_num(&topo, &shifted)
        }
        4 => {
            let topo = box_topo_general(
                [exps[0], exps[1], exps[2], exps[3]],
                [
                    dirs[0].clone(),
                    dirs[1].clone(),
                    dirs[2].clone(),
                    dirs[3].clone(),
                ],
                [
                    kmasses[0].clone(),
                    kmasses[1].clone(),
                    kmasses[2].clone(),
                    kmasses[3].clone(),
                ],
                clone_gram(base_gram, n_ext),
            );
            reduce_num(&topo, &shifted)
        }
        3 => {
            // Triangle invariants of the kept triple.
            let w_ab = dirs[1].clone();
            let w_bc: Dir = dirs[2].iter().zip(&dirs[1]).map(|(a, b)| a - b).collect();
            let w_ac = dirs[2].clone();
            let s1 = dir_dot_with(base_gram, &w_ab, &w_ab);
            let s2 = dir_dot_with(base_gram, &w_bc, &w_bc);
            let s3 = dir_dot_with(base_gram, &w_ac, &w_ac);
            let topo = triangle_topo_general(
                exps[0],
                exps[1],
                exps[2],
                [dirs[0].clone(), dirs[1].clone(), dirs[2].clone()],
                [kmasses[0].clone(), kmasses[1].clone(), kmasses[2].clone()],
                [s1, s2, s3],
                clone_gram(base_gram, n_ext),
            );
            reduce_num(&topo, &shifted)
        }
        2 => {
            let w = dirs[1].clone();
            let p_sq = dir_dot_with(base_gram, &w, &w);
            let topo = bubble_topo_general(
                exps[0],
                exps[1],
                vec![0i32; n_ext],
                w,
                &kmasses[0],
                &kmasses[1],
                p_sq,
                clone_gram(base_gram, n_ext),
            );
            reduce_num(&topo, &shifted)
        }
        1 => {
            let topo = tadpole_topo(exps[0], &kmasses[0], clone_gram(base_gram, n_ext), n_ext);
            reduce_num(&topo, &shifted)
        }
        _ => unreachable!("pinch keeps 1..=4 lines"),
    }
}

// A tadpole sub-topology: one line, no reducible directions; every loop dot is an ISP.
fn tadpole_topo(
    a1: i32,
    mass: &Atom,
    base_gram: Box<dyn Fn(usize, usize) -> Atom>,
    n_ext: usize,
) -> Topo {
    let den1 = den_symbol(0);
    let d1 = Atom::var(den1);
    let rule_ll = &d1 + mass;
    let m = mass.clone();
    Topo {
        rsp_dirs: vec![],
        base_gram,
        rule_ll,
        rule_lq: vec![],
        den_syms: vec![den1],
        a: vec![a1],
        r_coeffs: vec![vec![0i32; n_ext]],
        masses: vec![mass.clone()],
        scalar: Box::new(move |b| {
            // A0 with exponent <= 0 is a scaleless integral (no scale) -> 0.
            if b[0] <= 0 {
                return Vec::new();
            }
            let coeff = tadpole_coefficient(b[0], &m);
            if coeff == Atom::Zero {
                Vec::new()
            } else {
                vec![(coeff, MasterIntegral::Tadpole { m_sq: m.clone() })]
            }
        }),
        pinch: Box::new(|_, _, _| Vec::new()),
    }
}

// Reduce a box with an arbitrary-rank numerator.
#[allow(clippy::too_many_arguments)]
fn box_numerator(
    numerator: &Atom,
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
    reduce_num(
        &box_topo(a1, a2, a3, a4, p1, p2, p3, p4, s, t, m1, m2, m3, m4),
        &numerator_to_monos(numerator, 3),
    )
}

#[allow(clippy::too_many_arguments)]
fn box_topo(
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
) -> Topo {
    let den1 = den_symbol(0);
    let den2 = den_symbol(1);
    let den3 = den_symbol(2);
    let den4 = den_symbol(3);
    let d1 = Atom::var(den1);
    let d2 = Atom::var(den2);
    let d3 = Atom::var(den3);
    let d4 = Atom::var(den4);
    let two = Atom::num(2);
    // Reducible directions q1,q2,q3 (the full span); RSP rules below.
    let rule_ll = &d1 + m1;
    let rule_q1 = (&d2 - &d1 - m1 + m2 - p1) / &two;
    let rule_q2 = (&d3 - &d2 - m2 + m3 + p1 - s) / &two;
    let rule_q3 = (&d4 - &d3 - m3 + m4 - p4 + s) / &two;
    let base_gram = base_gram_box(p1, p2, p3, p4, s, t);
    let masses = [m1.clone(), m2.clone(), m3.clone(), m4.clone()];
    let r: [Dir; 4] = [vec![0, 0, 0], vec![1, 0, 0], vec![1, 1, 0], vec![1, 1, 1]];
    let (sp1, sp2, sp3, sp4) = (p1.clone(), p2.clone(), p3.clone(), p4.clone());
    let (ss, st) = (s.clone(), t.clone());
    let (sm1, sm2, sm3, sm4) = (m1.clone(), m2.clone(), m3.clone(), m4.clone());
    let pinch_gram = base_gram_box(p1, p2, p3, p4, s, t);
    let pinch_masses = masses.clone();
    Topo {
        rsp_dirs: vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]],
        base_gram,
        rule_ll,
        rule_lq: vec![rule_q1, rule_q2, rule_q3],
        den_syms: vec![den1, den2, den3, den4],
        a: vec![a1, a2, a3, a4],
        r_coeffs: r.to_vec(),
        masses: masses.to_vec(),
        scalar: Box::new(move |b| {
            let m = [sm1.clone(), sm2.clone(), sm3.clone(), sm4.clone()];
            let lex = [
                sp1.clone(),
                ss.clone(),
                sp4.clone(),
                sp2.clone(),
                st.clone(),
                sp3.clone(),
            ];
            reduce_cayley(&modified_cayley(&m, &lex), &m, b)
        }),
        pinch: Box::new(move |keep, exps, num| {
            pinch_to_subtopo(keep, exps, num, &r, &pinch_masses, &*pinch_gram)
        }),
    }
}

// A general box sub-topology: lines r[0..4], parent Gram for ISP projection and the six D0 invariants.
fn box_topo_general(
    exps: [i32; 4],
    r: [Dir; 4],
    masses: [Atom; 4],
    base_gram: Box<dyn Fn(usize, usize) -> Atom>,
) -> Topo {
    let dens: Vec<Symbol> = (0..4).map(den_symbol).collect();
    let d: Vec<Atom> = dens.iter().map(|&s| Atom::var(s)).collect();
    let two = Atom::num(2);
    let add = |u: &Dir, v: &Dir| -> Dir { u.iter().zip(v).map(|(a, b)| a + b).collect() };
    let sub = |u: &Dir, v: &Dir| -> Dir { u.iter().zip(v).map(|(a, b)| a - b).collect() };
    let w1 = sub(&r[1], &r[0]);
    let w2 = sub(&r[2], &r[1]);
    let w3 = sub(&r[3], &r[2]);
    // Box invariants (r_i - r_j)^2 from the parent Gram, in reduce_box's convention.
    let p1 = dir_dot_with(&*base_gram, &w1, &w1);
    let p2 = dir_dot_with(&*base_gram, &w2, &w2);
    let p3 = dir_dot_with(&*base_gram, &w3, &w3);
    let w12 = add(&w1, &w2);
    let w23 = add(&w2, &w3);
    let w123 = add(&w12, &w3);
    let s = dir_dot_with(&*base_gram, &w12, &w12);
    let t = dir_dot_with(&*base_gram, &w23, &w23);
    let p4 = dir_dot_with(&*base_gram, &w123, &w123);
    let (m1, m2, m3, m4) = (
        masses[0].clone(),
        masses[1].clone(),
        masses[2].clone(),
        masses[3].clone(),
    );
    // RSP rules, as in box_topo.
    let rule_ll = &d[0] + &m1;
    let rule_w1 = (&d[1] - &d[0] - &m1 + &m2 - &p1) / &two;
    let rule_w2 = (&d[2] - &d[1] - &m2 + &m3 + &p1 - &s) / &two;
    let rule_w3 = (&d[3] - &d[2] - &m3 + &m4 - &p4 + &s) / &two;
    let (sp1, sp2, sp3, sp4, ss, st) = (
        p1.clone(),
        p2.clone(),
        p3.clone(),
        p4.clone(),
        s.clone(),
        t.clone(),
    );
    let (sm1, sm2, sm3, sm4) = (m1.clone(), m2.clone(), m3.clone(), m4.clone());
    let r_for_pinch = r.clone();
    let masses_for_pinch = masses.clone();
    let gram_for_pinch = clone_gram(&*base_gram, r[0].len());
    Topo {
        rsp_dirs: vec![w1, w2, w3],
        base_gram,
        rule_ll,
        rule_lq: vec![rule_w1, rule_w2, rule_w3],
        den_syms: dens,
        a: exps.to_vec(),
        r_coeffs: r.to_vec(),
        masses: vec![m1, m2, m3, m4],
        scalar: Box::new(move |b| {
            let m = [sm1.clone(), sm2.clone(), sm3.clone(), sm4.clone()];
            let lex = [
                sp1.clone(),
                ss.clone(),
                sp4.clone(),
                sp2.clone(),
                st.clone(),
                sp3.clone(),
            ];
            reduce_cayley(&modified_cayley(&m, &lex), &m, b)
        }),
        pinch: Box::new(move |keep, exps, num| {
            pinch_to_subtopo(
                keep,
                exps,
                num,
                &r_for_pinch,
                &masses_for_pinch,
                &*gram_for_pinch,
            )
        }),
    }
}

// Reduce an N-gon (N >= 5) numerator
fn ngon_numerator(
    numerator: &Atom,
    exps: &[i32],
    invariants: &[Atom],
    masses: &[Atom],
) -> Vec<(Atom, MasterIntegral)> {
    let n = masses.len();
    let n_ext = n - 1;
    // Standard offsets r_i = q1 + ... + q_{i-1}, as directions over (q1..q_{n_ext}).
    let r_dirs: Vec<Dir> = (0..n)
        .map(|i| (0..n_ext).map(|k| (k < i) as i32).collect())
        .collect();
    reduce_num(
        &ngon_topo(
            exps,
            &r_dirs,
            masses,
            base_gram_from_pairwise(n_ext, invariants),
        ),
        &numerator_to_monos(numerator, n_ext),
    )
}

// A general K-point topology (K >= 5): RSP rules, reducible dirs, and Cayley invariants from the Gram.
fn ngon_topo(
    exps: &[i32],
    r_dirs: &[Dir],
    masses: &[Atom],
    base_gram: Box<dyn Fn(usize, usize) -> Atom>,
) -> Topo {
    let k = r_dirs.len();
    let n_ext = r_dirs[0].len();
    let dens: Vec<Symbol> = (0..k).map(den_symbol).collect();
    let d: Vec<Atom> = dens.iter().map(|&s| Atom::var(s)).collect();
    let two = Atom::num(2);
    let sub = |u: &Dir, v: &Dir| -> Dir { u.iter().zip(v).map(|(a, b)| a - b).collect() };
    // r_i^2 in the parent metric.
    let rr: Vec<Atom> = (0..k)
        .map(|i| dir_dot_with(&*base_gram, &r_dirs[i], &r_dirs[i]))
        .collect();
    let rule_ll = &d[0] + &masses[0];
    // l.q_a = [D_{a+1} - D_a - (r_{a+1}^2 - r_a^2) + m_{a+1} - m_a]/2, a = 0..k-1.
    let rule_lq: Vec<Atom> = (0..k - 1)
        .map(|a| (&d[a + 1] - &d[a] - (&rr[a + 1] - &rr[a]) + &masses[a + 1] - &masses[a]) / &two)
        .collect();
    // Reducible directions = the external momenta q_a = r_{a+1} - r_a.
    let rsp_dirs: Vec<Dir> = (0..k - 1)
        .map(|a| sub(&r_dirs[a + 1], &r_dirs[a]))
        .collect();

    let pairwise: Vec<Atom> = (0..k)
        .flat_map(|i| (i + 1..k).map(move |j| (i, j)))
        .map(|(i, j)| {
            let w = sub(&r_dirs[i], &r_dirs[j]);
            dir_dot_with(&*base_gram, &w, &w)
        })
        .collect();
    let smasses = masses.to_vec();
    let (r_for_pinch, pmasses) = (r_dirs.to_vec(), masses.to_vec());
    let gram_for_pinch = clone_gram(&*base_gram, n_ext);
    Topo {
        rsp_dirs,
        base_gram,
        rule_ll,
        rule_lq,
        den_syms: dens,
        a: exps.to_vec(),
        r_coeffs: r_dirs.to_vec(),
        masses: masses.to_vec(),
        scalar: Box::new(move |b| {
            reduce_cayley(&modified_cayley(&smasses, &pairwise), &smasses, b)
        }),
        pinch: Box::new(move |keep, exps, num| {
            pinch_to_subtopo(keep, exps, num, &r_for_pinch, &pmasses, &*gram_for_pinch)
        }),
    }
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

    #[test]
    fn det_matches_a_known_value() {
        crate::ensure_symbolica_license();
        // M = [[2,1,0],[1,3,1],[0,1,2]], det 8
        let m = vec![
            vec![Atom::num(2), Atom::num(1), Atom::num(0)],
            vec![Atom::num(1), Atom::num(3), Atom::num(1)],
            vec![Atom::num(0), Atom::num(1), Atom::num(2)],
        ];
        assert_eq!(super::det(&m), Atom::num(8));
    }

    #[test]
    fn scalar_pentagon_reduces_to_five_boxes() {
        crate::ensure_symbolica_license();
        let masses: Vec<Atom> = [1, 2, 3, 4, 5].iter().map(|&x| Atom::num(x)).collect();
        let invariants: Vec<Atom> = [3, 5, 7, 9, 4, 6, 8, 5, 7, 6]
            .iter()
            .map(|&x| Atom::num(x))
            .collect();
        let r = reduce(&scalar_family(masses, invariants));
        assert_eq!(r.terms.len(), 5);
        assert!(
            r.terms
                .iter()
                .all(|(_, m)| matches!(m, MasterIntegral::Box { .. }))
        );
        // van Neerven-Vermaseren coefficients
        let want_c = [
            Atom::num(-1088) / Atom::num(639),
            Atom::num(-212) / Atom::num(639),
            Atom::num(-44) / Atom::num(213),
            Atom::num(-191) / Atom::num(639),
            Atom::num(-341) / Atom::num(639),
        ];
        for ((coeff, _), want) in r.terms.iter().zip(&want_c) {
            assert_eq!(coeff, want);
        }
        // pinching the last propagator leaves the box on propagators 1..4
        assert_eq!(
            r.terms[4].1,
            MasterIntegral::Box {
                p1_sq: Atom::num(3),
                p2_sq: Atom::num(4),
                p3_sq: Atom::num(5),
                p4_sq: Atom::num(7),
                s: Atom::num(5),
                t: Atom::num(6),
                m1_sq: Atom::num(1),
                m2_sq: Atom::num(2),
                m3_sq: Atom::num(3),
                m4_sq: Atom::num(4),
            }
        );
    }

    #[test]
    fn scalar_hexagon_recurses_down_to_boxes() {
        crate::ensure_symbolica_license();
        let masses: Vec<Atom> = [1, 2, 3, 4, 5, 6].iter().map(|&x| Atom::num(x)).collect();
        let invariants: Vec<Atom> = [3, 5, 7, 9, 11, 4, 6, 8, 10, 5, 7, 9, 6, 8, 7]
            .iter()
            .map(|&x| Atom::num(x))
            .collect();
        let r = reduce(&scalar_family(masses, invariants));
        assert!(!r.terms.is_empty());
        assert!(
            r.terms
                .iter()
                .all(|(_, m)| matches!(m, MasterIntegral::Box { .. }))
        );
    }

    #[test]
    fn heptagon_numerator_reduces_to_finite_masters() {
        crate::ensure_symbolica_license();
        let masses: Vec<Atom> = (1..=7).map(Atom::num).collect();
        let invariants: Vec<Atom> = [
            4, -8, -20, -28, -12, -22, -8, -16, -28, -20, -26, -16, -38, -20, -16, -32, 4, -22,
            -10, -10, -16,
        ]
        .iter()
        .map(|&x| Atom::num(x))
        .collect();
        for num in [super::dot_lq(0), &super::dot_lq(0) * &super::dot_lq(5)] {
            let mut fam = family(masses.clone(), invariants.clone(), vec![1; 7]);
            fam.numerator = num;
            let r = reduce(&fam);
            assert!(!r.terms.is_empty());
            for (c, m) in &r.terms {
                assert!(c.to_string().is_ascii(), "non-finite heptagon coeff: {c}");
                assert!(matches!(
                    m,
                    MasterIntegral::Box { .. }
                        | MasterIntegral::Triangle { .. }
                        | MasterIntegral::Bubble { .. }
                        | MasterIntegral::Tadpole { .. }
                ));
            }
        }
    }

    #[test]
    #[ignore = "slow (~80s): dotted N>=7 FJT dimension-shift tree; verifies the gap-1 fix"]
    fn dotted_heptagon_reduces_to_finite_masters() {
        crate::ensure_symbolica_license();
        let masses: Vec<Atom> = (1..=7).map(Atom::num).collect();
        let invariants: Vec<Atom> = [
            4, -8, -20, -28, -12, -22, -8, -16, -28, -20, -26, -16, -38, -20, -16, -32, 4, -22,
            -10, -10, -16,
        ]
        .iter()
        .map(|&x| Atom::num(x))
        .collect();
        let r = reduce(&family(masses, invariants, vec![2, 1, 1, 1, 1, 1, 1]));
        assert!(!r.terms.is_empty());
        for (c, m) in &r.terms {
            assert!(
                c.to_string().is_ascii(),
                "non-finite dotted-heptagon coeff: {c}"
            );
            assert!(matches!(
                m,
                MasterIntegral::Box { .. }
                    | MasterIntegral::Triangle { .. }
                    | MasterIntegral::Bubble { .. }
                    | MasterIntegral::Tadpole { .. }
            ));
        }
    }

    #[test]
    fn base_gram_from_pairwise_matches_the_box_gram() {
        crate::ensure_symbolica_license();
        // physics box invariants (p1,p2,p3,p4,s,t)
        let n = |x| Atom::num(x);
        let phys = super::base_gram_box(&n(2), &n(3), &n(4), &n(5), &n(6), &n(7));
        let lex: Vec<Atom> = [2, 6, 5, 3, 7, 4].iter().map(|&x| n(x)).collect();
        let general = super::base_gram_from_pairwise(3, &lex);
        for a in 0..3 {
            for b in 0..3 {
                assert_eq!(general(a, b), phys(a, b), "q{a}.q{b}");
            }
        }
    }

    #[test]
    fn dotted_pentagon_reduces_to_valid_masters() {
        crate::ensure_symbolica_license();
        let masses: Vec<Atom> = [1, 2, 3, 4, 5].iter().map(|&x| Atom::num(x)).collect();
        let invariants: Vec<Atom> = [3, 5, 7, 9, 4, 6, 8, 5, 7, 6]
            .iter()
            .map(|&x| Atom::num(x))
            .collect();
        // pentagon with one squared propagator
        let r = reduce(&family(masses, invariants, vec![2, 1, 1, 1, 1]));
        assert!(!r.terms.is_empty());
        assert!(r.terms.iter().all(|(_, m)| matches!(
            m,
            MasterIntegral::Box { .. }
                | MasterIntegral::Triangle { .. }
                | MasterIntegral::Bubble { .. }
                | MasterIntegral::Tadpole { .. }
        )));
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, MasterIntegral::Triangle { .. }))
        );
    }

    #[test]
    fn pentagon_with_numerator_reduces_to_valid_masters() {
        crate::ensure_symbolica_license();
        let masses: Vec<Atom> = [1, 2, 3, 4, 5].iter().map(|&x| Atom::num(x)).collect();
        let invariants: Vec<Atom> = [3, 5, 7, 9, 4, 6, 8, 5, 7, 6]
            .iter()
            .map(|&x| Atom::num(x))
            .collect();
        let mut fam = family(masses, invariants, vec![1, 1, 1, 1, 1]);
        // numerator = l . q4
        fam.numerator = super::dot_lq(3);
        let r = reduce(&fam);
        assert!(!r.terms.is_empty());
        assert!(r.terms.iter().all(|(_, m)| matches!(
            m,
            MasterIntegral::Box { .. }
                | MasterIntegral::Triangle { .. }
                | MasterIntegral::Bubble { .. }
                | MasterIntegral::Tadpole { .. }
        )));
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, MasterIntegral::Box { .. }))
        );
    }

    #[test]
    fn hexagon_numerator_reduces_to_finite_masters() {
        crate::ensure_symbolica_license();
        let masses: Vec<Atom> = [1, 2, 3, 4, 5, 6].iter().map(|&x| Atom::num(x)).collect();
        let invariants: Vec<Atom> = [
            -4, -10, -20, -30, -30, -10, -24, -30, -22, -18, -38, -28, -34, -14, -18,
        ]
        .iter()
        .map(|&x| Atom::num(x))
        .collect();
        // rank-1 (l.q5, boxes only) and rank-2 (l.q1^2, pinches to pentagons w/ residual)
        for num in [super::dot_lq(4), super::dot_lq(0) * super::dot_lq(0)] {
            let mut fam = family(masses.clone(), invariants.clone(), vec![1; 6]);
            fam.numerator = num;
            let r = reduce(&fam);
            assert!(!r.terms.is_empty());
            for (c, m) in &r.terms {
                assert!(
                    c.to_string().is_ascii(),
                    "non-finite coeff (degenerate Gram): {c}"
                );
                assert!(matches!(
                    m,
                    MasterIntegral::Box { .. }
                        | MasterIntegral::Triangle { .. }
                        | MasterIntegral::Bubble { .. }
                        | MasterIntegral::Tadpole { .. }
                ));
            }
        }
    }

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
            kinematics: Kinematics { invariants },
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
    fn triangle_with_linear_numerator_reduces_to_masters() {
        crate::ensure_symbolica_license();
        let s1 = Atom::var(symbol!("oneloop::s1"));
        let s2 = Atom::var(symbol!("oneloop::s2"));
        let s3 = Atom::var(symbol!("oneloop::s3"));
        let m1 = Atom::var(symbol!("oneloop::m1sq"));
        let m2 = Atom::var(symbol!("oneloop::m2sq"));
        let m3 = Atom::var(symbol!("oneloop::m3sq"));
        // numerator = 2*(l.q1) - (l.q2) + (l.l)
        let numerator = Atom::num(2) * function!(S.dot, Atom::var(S.k), Atom::var(S.q1))
            - function!(S.dot, Atom::var(S.k), Atom::var(S.q2))
            + function!(S.dot, Atom::var(S.k), Atom::var(S.k));
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
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: m3,
                },
            ],
            isps: vec![],
            kinematics: Kinematics {
                invariants: vec![s1, s2, s3],
            },
            targets: vec![Integral {
                propagator_exponents: vec![1, 1, 1],
                isp_exponents: vec![],
            }],
            numerator,
        };
        let r = reduce(&fam);
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, MasterIntegral::Triangle { .. }))
        );
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, MasterIntegral::Bubble { .. }))
        );
        assert!(r.terms.iter().all(|(_, m)| matches!(
            m,
            MasterIntegral::Triangle { .. }
                | MasterIntegral::Bubble { .. }
                | MasterIntegral::Tadpole { .. }
        )));
    }

    #[test]
    fn scalar_triangle_reduces_to_unit_c0() {
        crate::ensure_symbolica_license();
        // Invariants are lexicographic pairwise (s01, s02, s12); the C0 legs come out in the
        // physical (adjacent, opposite) order (leg01, leg12, leg02) = (s01, s12, s02).
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
                assert_eq!(*p2_sq, Atom::num(3));
                assert_eq!(*p12_sq, Atom::num(2));
            }
            other => panic!("expected a triangle master, got {other:?}"),
        }
    }

    #[test]
    fn mixed_dotted_and_numerator_reduces_to_masters() {
        crate::ensure_symbolica_license();
        let v = |s: &str| Atom::var(symbol!(format!("oneloop::{s}")));
        let mut bub = family(vec![v("m1sq"), v("m2sq")], vec![v("psq")], vec![2, 1]);
        bub.numerator = super::dot_lq(0);
        let mut bx = family(
            vec![v("m1sq"), v("m2sq"), v("m3sq"), v("m4sq")],
            vec![v("p1"), v("p2"), v("p3"), v("p4"), v("sinv"), v("tinv")],
            vec![2, 1, 1, 1],
        );
        bx.numerator = super::dot_ll();
        for fam in [bub, bx] {
            let r = reduce(&fam);
            assert!(!r.terms.is_empty());
            for (c, m) in &r.terms {
                assert!(c.to_string().is_ascii(), "non-finite mixed coeff: {c}");
                assert!(matches!(
                    m,
                    MasterIntegral::Box { .. }
                        | MasterIntegral::Triangle { .. }
                        | MasterIntegral::Bubble { .. }
                        | MasterIntegral::Tadpole { .. }
                ));
            }
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
    fn box_with_linear_numerator_reduces_to_masters() {
        crate::ensure_symbolica_license();
        let [p1, p2, p3, p4, s, t, m1, m2, m3, m4] = box_syms();
        // numerator = (l.q1) + 3*(l.q3) - 2*(l.l)
        let numerator = function!(S.dot, Atom::var(S.k), Atom::var(S.q1))
            + Atom::num(3) * function!(S.dot, Atom::var(S.k), Atom::var(S.q3))
            - Atom::num(2) * function!(S.dot, Atom::var(S.k), Atom::var(S.k));
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
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: m3,
                },
                Propagator {
                    momentum: Atom::Zero,
                    mass_sq: m4,
                },
            ],
            isps: vec![],
            kinematics: Kinematics {
                invariants: vec![p1, p2, p3, p4, s, t],
            },
            targets: vec![Integral {
                propagator_exponents: vec![1, 1, 1, 1],
                isp_exponents: vec![],
            }],
            numerator,
        };
        let r = reduce(&fam);
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, MasterIntegral::Box { .. }))
        );
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, MasterIntegral::Triangle { .. }))
        );
        assert!(r.terms.iter().all(|(_, m)| matches!(
            m,
            MasterIntegral::Box { .. }
                | MasterIntegral::Triangle { .. }
                | MasterIntegral::Bubble { .. }
                | MasterIntegral::Tadpole { .. }
        )));
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
