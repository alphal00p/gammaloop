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

use symbolica::atom::{Atom, AtomCore};
use symbolica::{function, symbol};

use crate::symbols::S;

/// A gammaloop momentum tensor `H(id, mink(4, idx_))` 
fn known_momentum(head: &str, id: i64, oneloop_sym: Atom) -> (Atom, Atom) {
    let idx = Atom::var(symbol!("idx_"));
    let tensor = function!(
        symbol!(head),
        Atom::num(id),
        function!(symbol!("mink"), Atom::num(4), idx)
    );
    (tensor, oneloop_sym)
}

/// The loop/external momenta the bridge currently recognizes.
///
/// TODO: only externals up to `q3` (≤4-point) are mapped; generalize to build `q{j+1}` dynamically
/// so pentagons and beyond (`P(3,·)`, …) are handled too.
fn known_momenta() -> Vec<(Atom, Atom)> {
    vec![
        known_momentum("K", 0, Atom::var(S.k)),
        known_momentum("P", 0, Atom::var(S.q1)),
        known_momentum("P", 1, Atom::var(S.q2)),
        known_momentum("P", 2, Atom::var(S.q3)),
    ]
}

/// Rewrite the shared-index momentum contractions of a gammaloop scalar numerator into oneloop
/// `dot(...)` form. Self-contractions `a·a` (which Symbolica stores as squares) and contractions
/// `a·b` between two distinct momenta (shared-index products) both become `dot(...)`.
pub fn numerator_to_dot_form(num: &Atom) -> Atom {
    let moms = known_momenta();
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
    out
}

/// The external-momentum offset of a propagator, extracted from its gammaloop `lmb_rep`
pub fn external_offset_from_lmb_rep(lmb_rep: &Atom) -> Atom {
    // A single wildcard that swallows the whole `mink(4, idx)` argument.
    let any_index = Atom::var(symbol!("midx_"));
    let mut offset = lmb_rep.clone();
    for l in 0..MAX_MOMENTUM_ID {
        let loop_mom = function!(symbol!("K"), Atom::num(l), any_index.clone());
        offset = offset.replace(loop_mom.to_pattern()).with(Atom::Zero);
    }
    for j in 0..MAX_MOMENTUM_ID {
        let external = function!(symbol!("P"), Atom::num(j), any_index.clone());
        let q = Atom::var(symbol!(format!("oneloop::q{}", j + 1)));
        offset = offset.replace(external.to_pattern()).with(q);
    }
    offset
}

/// How many loop/external momentum ids the bridge recognizes (0..N).
const MAX_MOMENTUM_ID: i64 = 8;

#[cfg(test)]
mod tests {
    use super::*;

    fn mink4(idx: i64) -> Atom {
        function!(symbol!("mink"), Atom::num(4), Atom::num(idx))
    }

    #[test]
    fn contracts_loop_external_product_into_dot() {
        crate::ensure_symbolica_license();
        // gammaloop's `K(0,mink(4,5)) * P(0,mink(4,5))` is the scalar product k·q1.
        let input = &function!(symbol!("K"), Atom::num(0), mink4(5))
            * &function!(symbol!("P"), Atom::num(0), mink4(5));
        let got = numerator_to_dot_form(&input);
        let want = function!(S.dot, Atom::var(S.k), Atom::var(S.q1));
        assert_eq!(got, want);
    }

    #[test]
    fn contracts_loop_self_square_into_dot_kk() {
        crate::ensure_symbolica_license();
        // gammaloop's `K(0,mink(4,2))^2` is the scalar product k·k.
        let kmom = function!(symbol!("K"), Atom::num(0), mink4(2));
        let input = &kmom * &kmom;
        let got = numerator_to_dot_form(&input);
        let want = function!(S.dot, Atom::var(S.k), Atom::var(S.k));
        assert_eq!(got, want);
    }

    #[test]
    fn contracts_mixed_rank2_numerator() {
        crate::ensure_symbolica_license();
        // 2*(k·q1)*(k·k) - (q1·q2), written in gammaloop shared-index form.
        let k = |i: i64| function!(symbol!("K"), Atom::num(0), mink4(i));
        let p0 = |i: i64| function!(symbol!("P"), Atom::num(0), mink4(i));
        let p1 = |i: i64| function!(symbol!("P"), Atom::num(1), mink4(i));
        let input = Atom::num(2) * &k(1) * &p0(1) * &k(2) * &k(2) - &p0(3) * &p1(3);
        let got = numerator_to_dot_form(&input);
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
        // A bubble's two internal edges (the real lmb_reps `K(0,·)` and `-P(0,·)+K(0,·)`):
        assert_eq!(external_offset_from_lmb_rep(&k), Atom::Zero);
        let offset = external_offset_from_lmb_rep(&(&k - &p0));
        assert_eq!(offset, Atom::num(-1) * Atom::var(S.q1));
    }
}
