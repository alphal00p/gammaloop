//! rewrite a gammaloop *scalar* numerator into oneloop's dot convention.
//!
//! When gammaloop is asked for
//! `save dot --output-full-numerator --do-gamma-algebra --do-color-algebra`, the emitted numerator
//! is a polynomial in Minkowski dot products of the loop momentum `K(0,·)` and the external momenta
//! `P(i,·)`. A dot product `a·b` is written as a shared-index product
//! `a(mink(4,idx)) b(mink(4,idx))` — Einstein summation on the repeated `mink(4,idx)` — rather than
//! as a `dot(...)` function. The oneloop reducer instead consumes a polynomial in the symmetric-linear
//! [`crate::symbols`] `dot(k, q_i)` / `dot(k, k)`.
//!
//! Momentum map: loop `K(0,·)` -> `oneloop::k`; externals `P(0/1/2,·)` -> `oneloop::q1/q2/q3`.

use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::{function, symbol};

use crate::family::{Integral, IntegralFamily, Kinematics, Propagator};
use crate::symbols::S;

#[derive(Clone, Copy)]
pub struct GammaloopHeads {
    pub loop_mom: Symbol,
    pub external_mom: Symbol,
    pub index: Symbol,
    pub metric: Symbol,
}

fn bare_momentum(head: Symbol, id: i64, oneloop_sym: Atom, index: Symbol) -> (Atom, Atom) {
    let tensor = function!(head, Atom::num(id), function!(index, Atom::num(4)));
    (tensor, oneloop_sym)
}

/// The loop/external momenta in bare form
fn bare_momenta(heads: &GammaloopHeads) -> Vec<(Atom, Atom)> {
    vec![
        bare_momentum(heads.loop_mom, 0, Atom::var(S.k), heads.index),
        bare_momentum(heads.external_mom, 0, Atom::var(S.q1), heads.index),
        bare_momentum(heads.external_mom, 1, Atom::var(S.q2), heads.index),
        bare_momentum(heads.external_mom, 2, Atom::var(S.q3), heads.index),
    ]
}

/// A gammaloop momentum tensor `head(id, index(4, idx_))`
fn known_momentum(head: Symbol, id: i64, oneloop_sym: Atom, index: Symbol) -> (Atom, Atom) {
    let idx = Atom::var(symbol!("idx_"));
    let tensor = function!(head, Atom::num(id), function!(index, Atom::num(4), idx));
    (tensor, oneloop_sym)
}

/// The loop/external momenta the bridge currently recognizes.
///
/// TODO: only externals up to `q3` (≤4-point) are mapped; generalize to build `q{j+1}` dynamically
/// so pentagons and beyond (`P(3,·)`, …) are handled too.
fn known_momenta(heads: &GammaloopHeads) -> Vec<(Atom, Atom)> {
    vec![
        known_momentum(heads.loop_mom, 0, Atom::var(S.k), heads.index),
        known_momentum(heads.external_mom, 0, Atom::var(S.q1), heads.index),
        known_momentum(heads.external_mom, 1, Atom::var(S.q2), heads.index),
        known_momentum(heads.external_mom, 2, Atom::var(S.q3), heads.index),
    ]
}

/// Rewrite the shared-index momentum contractions of a gammaloop scalar numerator into oneloop
/// `dot(...)` form. Self-contractions `a·a` (which Symbolica stores as squares) and contractions
/// `a·b` between two distinct momenta (shared-index products) both become `dot(...)`.
pub fn numerator_to_dot_form(num: &Atom, heads: &GammaloopHeads) -> Atom {
    let moms = known_momenta(heads);
    let mut out = num.clone();
    // Self-contractions `a·a` appear as squares `mom(mink(4,i))^2` (the shared index makes the
    // two factors identical, so Symbolica collects them into a power).
    for (tensor, sym) in &moms {
        let square = tensor * tensor;
        let dot = function!(S.dot, sym.clone(), sym.clone());
        out = out.replace(square.to_pattern()).with(dot);
    }
    // Contractions `a·b` between two distinct momenta appear as shared-index products.
    for (i, (tensor_i, sym_i)) in moms.iter().enumerate() {
        for (j, (tensor_j, sym_j)) in moms.iter().enumerate() {
            if i == j {
                continue;
            }
            let contraction = tensor_i * tensor_j;
            let dot = function!(S.dot, sym_i.clone(), sym_j.clone());
            out = out.replace(contraction.to_pattern()).with(dot);
        }
    }
    // Metric-dot form `g(a, b)` produced by `simplify_metrics` (both self `a·a` and distinct `a·b`).
    let bare = bare_momenta(heads);
    for (tensor_i, sym_i) in &bare {
        for (tensor_j, sym_j) in &bare {
            let g = function!(heads.metric, tensor_i.clone(), tensor_j.clone());
            let dot = function!(S.dot, sym_i.clone(), sym_j.clone());
            out = out.replace(g.to_pattern()).with(dot);
        }
    }
    out
}

/// The external-momentum offset of a propagator, extracted from its gammaloop `lmb_rep`
pub fn external_offset_from_lmb_rep(lmb_rep: &Atom, heads: &GammaloopHeads) -> Atom {
    // A single wildcard that swallows the whole `mink(4, idx)` argument.
    let any_index = Atom::var(symbol!("midx_"));
    let mut offset = lmb_rep.clone();
    for l in 0..MAX_MOMENTUM_ID {
        let loop_mom = function!(heads.loop_mom, Atom::num(l), any_index.clone());
        offset = offset.replace(loop_mom.to_pattern()).with(Atom::Zero);
    }
    for j in 0..MAX_MOMENTUM_ID {
        let external = function!(heads.external_mom, Atom::num(j), any_index.clone());
        let q = Atom::var(symbol!(format!("oneloop::q{}", j + 1)));
        offset = offset.replace(external.to_pattern()).with(q);
    }
    offset
}

/// How many loop/external momentum ids the bridge recognizes (0..N).
const MAX_MOMENTUM_ID: i64 = 8;

pub struct GammaloopEdge {
    pub lmb_rep: Atom,
    pub mass_sq: Atom,
}

fn square_external_momentum(momentum: &Atom) -> Atom {
    let qs = [Atom::var(S.q1), Atom::var(S.q2), Atom::var(S.q3)];
    let mut out = (momentum * momentum).expand();
    // `q_a^2 -> dot(q_a, q_a)` (squares) then `q_a*q_b -> dot(q_a, q_b)`
    for qa in &qs {
        out = out
            .replace((qa * qa).to_pattern())
            .with(function!(S.dot, qa.clone(), qa.clone()));
    }
    for a in 0..qs.len() {
        for b in (a + 1)..qs.len() {
            out = out.replace((&qs[a] * &qs[b]).to_pattern()).with(function!(
                S.dot,
                qs[a].clone(),
                qs[b].clone()
            ));
        }
    }
    out
}

/// The `C(n,2)` pairwise invariants `(r_i - r_j)^2`
fn invariants_from_offsets(offsets: &[Atom]) -> Vec<Atom> {
    let mut invariants = Vec::new();
    for i in 0..offsets.len() {
        for j in (i + 1)..offsets.len() {
            invariants.push(square_external_momentum(&(&offsets[i] - &offsets[j])));
        }
    }
    invariants
}

/// Assemble a reducer [`IntegralFamily`] from a gammaloop one-loop numerator and its internal edges.
/// The numerator is translated to dot form; each edge contributes a massive propagator; and the
/// external kinematics are the pairwise invariants of the edges' external offsets.
pub fn family_from_gammaloop(
    numerator: &Atom,
    edges: &[GammaloopEdge],
    heads: &GammaloopHeads,
) -> IntegralFamily {
    let offsets: Vec<Atom> = edges
        .iter()
        .map(|e| external_offset_from_lmb_rep(&e.lmb_rep, heads))
        .collect();
    let n = edges.len();
    IntegralFamily {
        propagators: edges
            .iter()
            .map(|e| Propagator {
                momentum: Atom::Zero,
                mass_sq: e.mass_sq.clone(),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: invariants_from_offsets(&offsets),
        },
        targets: vec![Integral {
            propagator_exponents: vec![1; n],
            isp_exponents: vec![],
        }],
        numerator: numerator_to_dot_form(numerator, heads),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mink4(idx: i64) -> Atom {
        function!(symbol!("mink"), Atom::num(4), Atom::num(idx))
    }

    /// Standalone gammaloop heads for tests (self-consistent with the `K`/`P`/`mink` inputs built
    /// above; the real glue passes gammalooprs's `GS.loop_mom`/`GS.external_mom`/`spenso::mink`).
    fn heads() -> GammaloopHeads {
        GammaloopHeads {
            loop_mom: symbol!("K"),
            external_mom: symbol!("P"),
            index: symbol!("mink"),
            metric: symbol!("g"),
        }
    }

    #[test]
    fn contracts_metric_dot_form_into_dot() {
        crate::ensure_symbolica_license();
        // `simplify_metrics` can leave a scalar product as `g(a, b)` (bare momenta): map it too.
        let bare = |head: &str, id: i64| {
            function!(
                symbol!(head),
                Atom::num(id),
                function!(symbol!("mink"), Atom::num(4))
            )
        };
        let g = |a: Atom, b: Atom| function!(symbol!("g"), a, b);
        // g(K,K) = k·k ; g(K,P) = k·q1
        let input =
            &g(bare("K", 0), bare("K", 0)) + &(Atom::num(2) * g(bare("K", 0), bare("P", 0)));
        let got = numerator_to_dot_form(&input, &heads());
        let kk = function!(S.dot, Atom::var(S.k), Atom::var(S.k));
        let kq1 = function!(S.dot, Atom::var(S.k), Atom::var(S.q1));
        let want = (&kk + &(Atom::num(2) * &kq1)).expand();
        assert_eq!(got.expand(), want);
    }

    #[test]
    fn contracts_loop_external_product_into_dot() {
        crate::ensure_symbolica_license();

        let input = &function!(symbol!("K"), Atom::num(0), mink4(5))
            * &function!(symbol!("P"), Atom::num(0), mink4(5));
        let got = numerator_to_dot_form(&input, &heads());
        let want = function!(S.dot, Atom::var(S.k), Atom::var(S.q1));
        assert_eq!(got, want);
    }

    #[test]
    fn contracts_loop_self_square_into_dot_kk() {
        crate::ensure_symbolica_license();

        let kmom = function!(symbol!("K"), Atom::num(0), mink4(2));
        let input = &kmom * &kmom;
        let got = numerator_to_dot_form(&input, &heads());
        let want = function!(S.dot, Atom::var(S.k), Atom::var(S.k));
        assert_eq!(got, want);
    }

    #[test]
    fn contracts_mixed_rank2_numerator() {
        crate::ensure_symbolica_license();

        let k = |i: i64| function!(symbol!("K"), Atom::num(0), mink4(i));
        let p0 = |i: i64| function!(symbol!("P"), Atom::num(0), mink4(i));
        let p1 = |i: i64| function!(symbol!("P"), Atom::num(1), mink4(i));
        let input = Atom::num(2) * &k(1) * &p0(1) * &k(2) * &k(2) - &p0(3) * &p1(3);
        let got = numerator_to_dot_form(&input, &heads());
        let kk = function!(S.dot, Atom::var(S.k), Atom::var(S.k));
        let kq1 = function!(S.dot, Atom::var(S.k), Atom::var(S.q1));
        let q1q2 = function!(S.dot, Atom::var(S.q1), Atom::var(S.q2));
        let want = (Atom::num(2) * &kq1 * &kk - &q1q2).expand();
        assert_eq!(got.expand(), want);
    }

    #[test]
    fn external_offset_drops_loop_and_maps_externals() {
        crate::ensure_symbolica_license();
        let k = function!(symbol!("K"), Atom::num(0), mink4(9));
        let p0 = function!(symbol!("P"), Atom::num(0), mink4(9));
        assert_eq!(external_offset_from_lmb_rep(&k, &heads()), Atom::Zero);
        let offset = external_offset_from_lmb_rep(&(&k - &p0), &heads());
        assert_eq!(offset, Atom::num(-1) * Atom::var(S.q1));
    }

    #[test]
    fn bubble_invariants_from_offsets() {
        crate::ensure_symbolica_license();

        let offsets = vec![Atom::Zero, Atom::num(-1) * Atom::var(S.q1)];
        let got = invariants_from_offsets(&offsets);
        assert_eq!(
            got,
            vec![function!(S.dot, Atom::var(S.q1), Atom::var(S.q1))]
        );
    }

    #[test]
    fn scalar_massless_bubble_reduces_to_b0() {
        crate::ensure_symbolica_license();

        let k = function!(symbol!("K"), Atom::num(0), mink4(0));
        let p0 = function!(symbol!("P"), Atom::num(0), mink4(0));
        let edges = vec![
            GammaloopEdge {
                lmb_rep: k.clone(),
                mass_sq: Atom::Zero,
            },
            GammaloopEdge {
                lmb_rep: &k - &p0,
                mass_sq: Atom::Zero,
            },
        ];
        let fam = family_from_gammaloop(&Atom::num(1), &edges, &heads());
        let r = crate::reduce::reduce(&fam);
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, crate::masters::MasterIntegral::Bubble { .. })),
            "expected a B0 master, got {:?}",
            r.terms.iter().map(|(_, m)| m).collect::<Vec<_>>()
        );
    }

    #[test]
    fn rank1_bubble_numerator_reduces_through_the_bridge() {
        crate::ensure_symbolica_license();

        let k = |i: i64| function!(symbol!("K"), Atom::num(0), mink4(i));
        let p0 = |i: i64| function!(symbol!("P"), Atom::num(0), mink4(i));
        let numerator = &k(1) * &p0(1);
        let edges = vec![
            GammaloopEdge {
                lmb_rep: k(0),
                mass_sq: Atom::Zero,
            },
            GammaloopEdge {
                lmb_rep: &k(0) - &p0(0),
                mass_sq: Atom::Zero,
            },
        ];
        let fam = family_from_gammaloop(&numerator, &edges, &heads());
        assert_eq!(
            fam.numerator,
            function!(S.dot, Atom::var(S.k), Atom::var(S.q1))
        );
        let r = crate::reduce::reduce(&fam);
        assert!(!r.terms.is_empty());
        assert!(
            r.terms
                .iter()
                .any(|(_, m)| matches!(m, crate::masters::MasterIntegral::Bubble { .. })),
            "expected a B0 master, got {:?}",
            r.terms.iter().map(|(_, m)| m).collect::<Vec<_>>()
        );
    }
}
