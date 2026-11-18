//! Numerator profile of an R* remainder for the numerator-aware routing keys
//! (numerator-routing lane, `GAMMALOOP_ROOT_ROUTING=num-*`).
//!
//! The remainder numerator (numerator x insertions x closure, before routing)
//! depends on the routing only through its edge momenta `Q(e, ...)`: routed,
//! `Q(e)` becomes a sum of `l_e` loop and `r_e` boundary momenta. The profile
//! is a *count circuit* of the expression: every term counts 1 (coefficients
//! and non-momentum factors are 1), `Q(e)` is the variable `x_e`, `den(...)` and
//! negative powers are 1, and sums of products are kept as small polynomials
//! in the `x_e` (collapsed while their size stays below [`CAP`]). Evaluated at
//! `x_e = l_e + r_e y` it gives the number of terms of the routed, expanded
//! numerator per order `y^i` in the boundary momenta (an upper bound: no
//! cancellation, no symmetrisation). No value of the integrand depends on it.

use std::collections::BTreeMap;

use symbolica::atom::{AtomCore, AtomView};

use crate::utils::GS;

/// Largest polynomial kept collapsed; larger products stay circuit nodes.
const CAP: usize = 4096;

/// One monomial: sorted (edge, power) pairs.
type Mono = Vec<(u32, u32)>;

#[derive(Clone, Debug)]
enum Node {
    /// Σ coefficient x Π x_e^power.
    Poly(Vec<(Mono, f64)>),
    Mul(Vec<Node>),
    Add(Vec<Node>),
    Pow(Box<Node>, u32),
}

fn one() -> Node {
    Node::Poly(vec![(Vec::new(), 1.0)])
}

fn mono_mul(a: &Mono, b: &Mono) -> Mono {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        if j == b.len() || (i < a.len() && a[i].0 < b[j].0) {
            out.push(a[i]);
            i += 1;
        } else if i == a.len() || b[j].0 < a[i].0 {
            out.push(b[j]);
            j += 1;
        } else {
            out.push((a[i].0, a[i].1 + b[j].1));
            i += 1;
            j += 1;
        }
    }
    out
}

fn poly_mul(a: &[(Mono, f64)], b: &[(Mono, f64)]) -> Vec<(Mono, f64)> {
    let mut map: BTreeMap<Mono, f64> = BTreeMap::new();
    for (ma, ca) in a {
        for (mb, cb) in b {
            *map.entry(mono_mul(ma, mb)).or_insert(0.0) += ca * cb;
        }
    }
    map.into_iter().collect()
}

fn poly_add(parts: Vec<Vec<(Mono, f64)>>) -> Vec<(Mono, f64)> {
    let mut map: BTreeMap<Mono, f64> = BTreeMap::new();
    for part in parts {
        for (m, c) in part {
            *map.entry(m).or_insert(0.0) += c;
        }
    }
    map.into_iter().collect()
}

/// Multiply a list of nodes, collapsing adjacent polynomials while small.
fn make_mul(children: Vec<Node>) -> Node {
    let mut acc: Vec<(Mono, f64)> = vec![(Vec::new(), 1.0)];
    let mut rest = Vec::new();
    for child in children {
        match child {
            Node::Poly(p) => {
                if p.len() == 1 && p[0].0.is_empty() {
                    acc.iter_mut().for_each(|(_, c)| *c *= p[0].1);
                } else if acc.len() * p.len() <= CAP {
                    acc = poly_mul(&acc, &p);
                } else {
                    rest.push(Node::Poly(p));
                }
            }
            Node::Mul(inner) => rest.extend(inner),
            other => rest.push(other),
        }
    }
    if rest.is_empty() {
        return Node::Poly(acc);
    }
    if !(acc.len() == 1 && acc[0].0.is_empty() && acc[0].1 == 1.0) {
        rest.push(Node::Poly(acc));
    }
    if rest.len() == 1 {
        rest.pop().unwrap()
    } else {
        Node::Mul(rest)
    }
}

fn make_add(children: Vec<Node>) -> Node {
    let (polys, rest): (Vec<Node>, Vec<Node>) =
        children.into_iter().partition(|c| matches!(c, Node::Poly(_)));
    let merged = poly_add(
        polys
            .into_iter()
            .map(|p| match p {
                Node::Poly(p) => p,
                _ => unreachable!(),
            })
            .collect(),
    );
    if rest.is_empty() {
        return Node::Poly(merged);
    }
    let mut rest = rest;
    if !merged.is_empty() {
        rest.push(Node::Poly(merged));
    }
    if rest.len() == 1 { rest.pop().unwrap() } else { Node::Add(rest) }
}

fn build(view: AtomView<'_>) -> Node {
    match view {
        AtomView::Num(_) | AtomView::Var(_) => one(),
        AtomView::Fun(f) => {
            if f.get_symbol() == GS.emr_mom {
                match usize::try_from(f.get(0)) {
                    Ok(edge) => Node::Poly(vec![(vec![(edge as u32, 1)], 1.0)]),
                    Err(_) => one(),
                }
            } else if f.get_symbol() == GS.den {
                one()
            } else {
                make_mul(f.iter().map(build).collect())
            }
        }
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            match i64::try_from(exponent) {
                Ok(n) if n > 0 => {
                    let base = build(base);
                    match base {
                        Node::Poly(ref b) if b.len() == 1 && b[0].0.is_empty() => one(),
                        Node::Poly(b) => {
                            let mut acc = vec![(Vec::new(), 1.0)];
                            let mut ok = true;
                            for _ in 0..n {
                                if acc.len() * b.len() > CAP {
                                    ok = false;
                                    break;
                                }
                                acc = poly_mul(&acc, &b);
                            }
                            if ok {
                                Node::Poly(acc)
                            } else {
                                Node::Pow(Box::new(Node::Poly(b)), n as u32)
                            }
                        }
                        other => Node::Pow(Box::new(other), n as u32),
                    }
                }
                _ => one(),
            }
        }
        AtomView::Mul(m) => make_mul(m.iter().map(build).collect()),
        AtomView::Add(a) => make_add(a.iter().map(build).collect()),
    }
}

/// Truncated power series in the boundary scaling `y` (orders `0..=degree`).
type Series = Vec<f64>;

fn series_mul(a: &[f64], b: &[f64], degree: usize) -> Series {
    let mut out = vec![0.0; degree + 1];
    for (i, x) in a.iter().enumerate() {
        if *x == 0.0 {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if i + j > degree {
                break;
            }
            out[i + j] += x * y;
        }
    }
    out
}

/// The numerator profile of one remainder.
#[derive(Clone, Debug)]
pub(crate) struct NumeratorProfile {
    root: Node,
    /// Edges whose momentum appears, with their maximal power in one term.
    pub max_power: BTreeMap<u32, u32>,
    /// Maximal total momentum degree of one term.
    pub max_degree: u32,
    /// Circuit nodes and collapsed monomials (diagnostics).
    pub nodes: usize,
    pub monomials: usize,
}

fn max_powers(node: &Node) -> (BTreeMap<u32, u32>, u32) {
    match node {
        Node::Poly(p) => {
            let mut map = BTreeMap::new();
            let mut degree = 0;
            for (m, _) in p {
                let mut d = 0;
                for (e, k) in m {
                    let entry = map.entry(*e).or_insert(0);
                    *entry = (*entry).max(*k);
                    d += k;
                }
                degree = degree.max(d);
            }
            (map, degree)
        }
        Node::Mul(children) => {
            let mut map = BTreeMap::new();
            let mut degree = 0;
            for child in children {
                let (m, d) = max_powers(child);
                for (e, k) in m {
                    *map.entry(e).or_insert(0) += k;
                }
                degree += d;
            }
            (map, degree)
        }
        Node::Add(children) => {
            let mut map = BTreeMap::new();
            let mut degree = 0;
            for child in children {
                let (m, d) = max_powers(child);
                for (e, k) in m {
                    let entry = map.entry(e).or_insert(0);
                    *entry = (*entry).max(k);
                }
                degree = degree.max(d);
            }
            (map, degree)
        }
        Node::Pow(base, n) => {
            let (m, d) = max_powers(base);
            (m.into_iter().map(|(e, k)| (e, k * n)).collect(), d * n)
        }
    }
}

fn sizes(node: &Node) -> (usize, usize) {
    match node {
        Node::Poly(p) => (1, p.len()),
        Node::Mul(c) | Node::Add(c) => c
            .iter()
            .map(sizes)
            .fold((1, 0), |(a, b), (x, y)| (a + x, b + y)),
        Node::Pow(b, _) => {
            let (a, m) = sizes(b);
            (a + 1, m)
        }
    }
}

impl NumeratorProfile {
    pub(crate) fn new(view: AtomView<'_>) -> Self {
        let root = build(view);
        let (max_power, max_degree) = max_powers(&root);
        let (nodes, monomials) = sizes(&root);
        Self { root, max_power, max_degree, nodes, monomials }
    }

    /// Terms of the routed numerator per boundary order `0..=degree`, with
    /// `terms[e] = (l_e, r_e)` (an edge absent from `terms` counts as one
    /// boundary momentum).
    pub(crate) fn orders(&self, terms: &dyn Fn(u32) -> (f64, f64), degree: usize) -> Series {
        let mut cache: BTreeMap<(u32, u32), Series> = BTreeMap::new();
        eval(&self.root, terms, degree, &mut cache)
    }
}

fn eval(
    node: &Node,
    terms: &dyn Fn(u32) -> (f64, f64),
    degree: usize,
    cache: &mut BTreeMap<(u32, u32), Series>,
) -> Series {
    match node {
        Node::Poly(p) => {
            let mut total = vec![0.0; degree + 1];
            for (m, c) in p {
                let mut s = vec![0.0; degree + 1];
                s[0] = *c;
                for (e, k) in m {
                    let factor = cache.entry((*e, *k)).or_insert_with(|| {
                        let (l, r) = terms(*e);
                        let mut base = vec![0.0; degree + 1];
                        base[0] = l;
                        if degree >= 1 {
                            base[1] = r;
                        }
                        let mut acc = vec![0.0; degree + 1];
                        acc[0] = 1.0;
                        for _ in 0..*k {
                            acc = series_mul(&acc, &base, degree);
                        }
                        acc
                    });
                    s = series_mul(&s, factor, degree);
                }
                for (t, x) in total.iter_mut().zip(s) {
                    *t += x;
                }
            }
            total
        }
        Node::Mul(children) => {
            let mut acc = vec![0.0; degree + 1];
            acc[0] = 1.0;
            for child in children {
                acc = series_mul(&acc, &eval(child, terms, degree, cache), degree);
            }
            acc
        }
        Node::Add(children) => {
            let mut total = vec![0.0; degree + 1];
            for child in children {
                for (t, x) in total.iter_mut().zip(eval(child, terms, degree, cache)) {
                    *t += x;
                }
            }
            total
        }
        Node::Pow(base, n) => {
            let b = eval(base, terms, degree, cache);
            let mut acc = vec![0.0; degree + 1];
            acc[0] = 1.0;
            for _ in 0..*n {
                acc = series_mul(&acc, &b, degree);
            }
            acc
        }
    }
}

/// The den-product Taylor shape estimate: per den carrying `r > 0` boundary
/// and `l` loop momenta, `1 + (l r) y + ((l r)^2 + r^2) y^2 + ...`.
pub(crate) fn den_orders(dens: &[(f64, f64)], degree: usize) -> Series {
    let mut acc = vec![0.0; degree + 1];
    acc[0] = 1.0;
    for &(l, r) in dens {
        if r == 0.0 {
            continue;
        }
        let mut f = vec![0.0; degree + 1];
        f[0] = 1.0;
        let a = l * r + r;
        let mut p = 1.0;
        for (i, slot) in f.iter_mut().enumerate().skip(1) {
            p *= a;
            *slot = p + if i >= 2 { r * r } else { 0.0 };
        }
        acc = series_mul(&acc, &f, degree);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn series_and_dens() {
        assert_eq!(series_mul(&[1.0, 2.0, 0.0], &[3.0, 1.0, 0.0], 2), vec![3.0, 7.0, 2.0]);
        let d = den_orders(&[(2.0, 1.0)], 2);
        assert_eq!(d, vec![1.0, 3.0, 10.0]);
    }
}
