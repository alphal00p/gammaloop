//! Reduction of closed colour networks to group invariants.
//!
//! A product of colour factors whose indices are all contracted (a vacuum
//! colour graph) is rewritten with the identities of color.h (T. van
//! Ritbergen, A.N. Schellekens, J.A.M. Vermaseren, "Group theory factors for
//! Feynman diagrams", hep-ph/9802376): fundamental and adjoint traces lose
//! their contracted generator pairs, become symmetrized traces with
//! commutator terms, structure constants between symmetrized traces move by
//! the invariance of the symmetric tensors, and the remaining networks of
//! symmetric tensors are named. Subnetworks with two legs factorize (Schur's
//! lemma). Every intermediate network is canonically labelled, so equal
//! subnetworks are reduced once per thread.

mod atoms;
mod graph;
mod names;
mod poly;
mod reduce;
mod rules;
mod su_n;

use std::{collections::HashMap, sync::LazyLock};

use symbolica::atom::{Atom, AtomView, Symbol};

pub(crate) use reduce::Failure;

static SYMBOLS: LazyLock<HashMap<&'static str, Symbol>> = LazyLock::new(|| {
    names::template_names()
        .map(|name| (name, symbolica::symbol!(&format!("spenso::{name}"); Real, Scalar)))
        .collect()
});

/// The function symbol of a named invariant with more than two tensors.
pub(crate) fn invariant_symbol(name: &'static str) -> Symbol {
    SYMBOLS[name]
}

/// The name of a named-invariant function symbol.
pub(crate) fn invariant_name(symbol: Symbol) -> Option<&'static str> {
    SYMBOLS.iter().find(|(_, s)| **s == symbol).map(|(n, _)| *n)
}

/// SU(N) value (T_F = 1/2) of an invariant as (power of N, numerator,
/// denominator) terms; `adjoint` marks the adjoint arguments.
pub(crate) fn su_n_terms(name: &str, adjoint: &[bool]) -> Option<&'static [(i32, i64, i64)]> {
    let reps: Vec<graph::Rep> = adjoint.iter().map(|&a| if a { graph::Rep::A } else { graph::Rep::R }).collect();
    su_n::value(name, &reps)
}

/// A product of colour factors read as a closed network, not yet reduced.
pub(crate) struct ClosedNetwork {
    net: graph::Net,
    coefficient: poly::C,
    ctx: atoms::Context,
}

/// The closed network of a product of colour factors with all colour
/// indices summed, if the factors are a supported closed network.
pub(crate) fn parse_closed(factors: &[AtomView<'_>]) -> Option<ClosedNetwork> {
    let (net, coefficient, ctx) = atoms::parse(factors)?;
    Some(ClosedNetwork { net, coefficient, ctx })
}

impl ClosedNetwork {
    /// The group invariants of the network: `cas`, `idx`, `gram`,
    /// dimensions and the higher named invariants.
    pub(crate) fn reduce(self) -> Result<Atom, Failure> {
        reduce::trim_memo();
        poly::take_overflow();
        let mut reducer = reduce::Reducer::new();
        let value = reducer.evaluate(self.net).map(|v| v.scale(&self.coefficient).with_fundamental_closures());
        if poly::take_overflow() {
            // Values computed in this reduction (and memoized) are garbage.
            reduce::clear_memo();
            return Err(Failure::Overflow);
        }
        Ok(atoms::to_atom(&value?, &self.ctx))
    }
}

/// Reduce a product of colour factors with all colour indices contracted to
/// group invariants. `Ok(None)` if the factors are not a supported closed
/// network.
pub(crate) fn reduce_closed(factors: &[AtomView<'_>]) -> Result<Option<Atom>, Failure> {
    parse_closed(factors).map(ClosedNetwork::reduce).transpose()
}

/// Reduce a closed colour network given as a product (or a single factor).
/// Returns `None` when the expression is not a closed network of supported
/// colour factors, or when the reduction does not finish.
pub fn reduce_closed_color_network(expression: AtomView<'_>) -> Option<Atom> {
    let factors: Vec<AtomView<'_>> = match expression {
        AtomView::Mul(m) => m.iter().collect(),
        v => vec![v],
    };
    reduce_closed(&factors).ok().flatten()
}

/// As [`reduce_closed_color_network`], with the reason of a failure.
pub fn try_reduce_closed_color_network(expression: AtomView<'_>) -> Result<Option<Atom>, String> {
    let factors: Vec<AtomView<'_>> = match expression {
        AtomView::Mul(m) => m.iter().collect(),
        v => vec![v],
    };
    reduce_closed(&factors).map_err(|e| format!("{e:?}"))
}

/// The value of a color.h invariant `name(args)` (for example `d433` with
/// `[true, false, false]` for d433(A1, R1, R2)), after color.h's relations,
/// in the representations `fundamental` = cof(N) and `adjoint` = coad(N_A).
pub fn named_invariant(name: &str, adjoint_args: &[bool], fundamental: Atom, adjoint: Atom) -> Option<Atom> {
    let reps: Vec<graph::Rep> = adjoint_args.iter().map(|&a| if a { graph::Rep::A } else { graph::Rep::R }).collect();
    let inv = names::invariant(name, &reps)?;
    let value = names::value(inv);
    Some(atoms::to_atom(&value, &atoms::Context { fundamental: Some(fundamental), adjoint: Some(adjoint) }))
}

/// Drop this thread's cache of reduced networks.
pub fn clear_closed_network_cache() {
    reduce::clear_memo();
}

/// Entries and estimated bytes of this thread's cache of reduced networks.
/// Before a network is reduced the cache is dropped if it is over a fixed
/// per-thread budget (32 MiB), so it does not grow with the number of calls.
pub fn closed_network_cache_stats() -> (usize, usize) {
    reduce::memo_stats()
}

/// Networks reduced (cache misses) on this thread since the last call.
pub fn take_closed_network_visits() -> usize {
    reduce::take_visits()
}

#[cfg(test)]
mod tests {
    use symbolica::parse_lit;
    use symbolica_utils::AtomPrintExt;

    #[test]
    fn open_networks_are_not_closed() {
        crate::test_support::test_initialize();
        let open = parse_lit!(
            f(coad(dA, a), coad(dA, c), coad(dA, x)) * f(coad(dA, b), coad(dA, c), coad(dA, x)),
            default_namespace = "spenso"
        );
        assert_eq!(super::try_reduce_closed_color_network(open.as_view()), Ok(None), "{open}");
        let closed = parse_lit!(
            f(coad(dA, a), coad(dA, c), coad(dA, x)) * f(coad(dA, a), coad(dA, c), coad(dA, x)),
            default_namespace = "spenso"
        );
        let value = super::reduce_closed_color_network(closed.as_view()).unwrap();
        assert_eq!(value.to_bare_ordered_string(), "cas(2,coad(dA))*dA");
    }

    /// The memo stays within its budget across calls: once over it, the next
    /// top-level network starts from an empty cache. Values do not depend on
    /// the budget, and the byte count matches the entries.
    #[test]
    fn closed_network_cache_respects_its_budget() {
        use crate::color::CS;
        use crate::representations::{ColorAdjoint, ColorFundamental};
        use spenso::structure::representation::RepName;
        use symbolica::{atom::Atom, symbol};
        crate::test_support::test_initialize();
        let a = |i: String| {
            ColorAdjoint {}.to_symbolic([Atom::var(symbol!("spenso::dA")), Atom::var(symbol!(&format!("memo_budget::{i}")))])
        };
        let r = |i: String| {
            ColorFundamental {}.to_symbolic([Atom::var(symbol!("spenso::Nc")), Atom::var(symbol!(&format!("memo_budget::{i}")))])
        };
        let dual = |x: Atom| symbolica::function!(spenso::structure::abstract_index::AIND_SYMBOLS.dind, x);
        let mut networks = Vec::new();
        for n in 2..=5 {
            // Prism of n rungs of structure constants.
            let mut prism = Atom::num(1);
            for k in 0..n {
                let p = (k + n - 1) % n;
                prism = prism
                    * CS.structure_f(a(format!("t{p}")), a(format!("t{k}")), a(format!("r{k}")))
                    * CS.structure_f(a(format!("u{p}")), a(format!("u{k}")), a(format!("r{k}")));
            }
            networks.push(prism);
            // A quark loop with n gluons attached in order and in reverse.
            let mut loops = Atom::num(1);
            for k in 0..n {
                let (next, back) = ((k + 1) % n, (k + n - 1) % n);
                loops = loops
                    * CS.explicit_t(a(format!("g{k}")), r(format!("i{k}")), dual(r(format!("i{next}"))))
                    * CS.explicit_t(a(format!("g{k}")), r(format!("j{k}")), dual(r(format!("j{back}"))));
            }
            networks.push(loops);
        }
        let reduce = |x: &Atom| super::reduce_closed_color_network(x.as_view()).unwrap_or_else(|| panic!("{x}"));
        let consistent = || {
            let (_, bytes) = super::closed_network_cache_stats();
            assert_eq!(bytes, super::reduce::recount_memo_bytes());
            bytes
        };
        // Cold values and the cache each network fills on its own.
        let mut cold = Vec::new();
        for x in &networks {
            super::clear_closed_network_cache();
            assert_eq!(super::closed_network_cache_stats(), (0, 0));
            cold.push((reduce(x), consistent()));
        }
        assert!(cold.iter().filter(|(_, bytes)| *bytes > 0).count() >= 4, "{:?}", cold.iter().map(|c| c.1).collect::<Vec<_>>());
        let budget = cold.iter().map(|(_, b)| *b).max().unwrap();
        super::reduce::set_memo_budget(budget);
        super::clear_closed_network_cache();
        let mut trims = 0;
        for _ in 0..3 {
            for (x, (value, own)) in networks.iter().zip(&cold) {
                let before = consistent();
                assert_eq!(&reduce(x), value, "{x}");
                let after = consistent();
                if before > budget {
                    trims += 1;
                    assert!(after <= *own, "{after} > {own} after a trim");
                }
                assert!(after <= budget + *own, "{after} > {budget} + {own}");
            }
        }
        assert!(trims > 0);
        super::reduce::set_memo_budget(super::reduce::MEMO_BUDGET);
        super::clear_closed_network_cache();
    }

    /// Fixed components are not summed: (f^{123})^2 is one component.
    #[test]
    fn concrete_components_are_not_closed() {
        crate::test_support::test_initialize();
        let fixed = parse_lit!(
            f(coad(8, cind(1)), coad(8, cind(2)), coad(8, cind(3))) ^ 2,
            default_namespace = "spenso"
        );
        assert_eq!(super::try_reduce_closed_color_network(fixed.as_view()), Ok(None), "{fixed}");
    }

    /// Every relation applied to a named invariant (color.h's and the ones it
    /// does not apply: d4433b(R,A,R,R), d653(A,R,R), dddff(R,R,A)) agrees in
    /// SU(N) with the invariant's own SU(N) value, computed with SUn.prc from
    /// its defining network.
    #[test]
    fn relations_agree_with_su_n_values() {
        use crate::color::{CS, ColorSimplifier};
        use crate::representations::{ColorAdjoint, ColorFundamental};
        use spenso::structure::representation::RepName;
        use symbolica::atom::{Atom, AtomCore};
        crate::test_support::test_initialize();
        let n = Atom::var(CS.nc);
        let fundamental = ColorFundamental {}.to_symbolic([n.clone()]);
        let adjoint = ColorAdjoint {}.to_symbolic([n.clone().pow(Atom::num(2)) - Atom::num(1)]);
        let mut checked = 0;
        // Entries with a relation: the others compare the table with itself.
        let mut related = 0;
        for (name, reps, terms) in super::su_n::entries() {
            let entry = super::names::invariant(name, &reps[..]).unwrap_or_else(|| panic!("{name} {reps:?}"));
            if super::names::value(entry) != super::poly::Poly::var(super::poly::FIRST_NAMED + entry) {
                related += 1;
            }
            let adjoint_args: Vec<bool> = reps.iter().map(|r| *r == super::graph::Rep::A).collect();
            let value = super::named_invariant(name, &adjoint_args, fundamental.clone(), adjoint.clone())
                .unwrap_or_else(|| panic!("{name} {reps:?}"));
            let value = value.to_cof_dimension_invariants();
            let mut expected = Atom::Zero;
            for &(power, numerator, denominator) in terms {
                expected = expected + Atom::num(numerator) / Atom::num(denominator) * n.clone().pow(Atom::num(power as i64));
            }
            let difference = (value.clone() - &expected).expand().together().expand();
            assert!(difference.is_zero(), "{name} {reps:?}: {value} != {expected}");
            checked += 1;
        }
        assert_eq!(checked, super::su_n::entries().count());
        assert!(checked > 20, "{checked}");
        assert!(related >= 4, "{related} entries with a relation");
    }
}
