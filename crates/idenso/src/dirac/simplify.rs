use spenso::{
    chain,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    structure::representation::LibraryRep,
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    id::Context,
    utils::Settable,
};

use crate::{chain::Chain, metric::MetricSimplifier, representations::Bispinor};

use super::{AGS, id_atom};

/// Controls how open gamma chains are reordered during simplification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GammaChainOrdering {
    /// Only move a repeated gamma toward its mate, matching FORM's conservative
    /// strategy for avoiding unnecessary terms.
    RepeatedPairs,
    /// Canonicalize gamma order with adjacent Clifford swaps.
    Canonical,
}

/// Settings for the chain-based Dirac simplifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GammaSimplifySettings {
    /// Ordering strategy for open chains.
    pub chain_ordering: GammaChainOrdering,
    /// Whether closed chains should be evaluated as traces.
    pub evaluate_traces: bool,
}

impl Default for GammaSimplifySettings {
    fn default() -> Self {
        Self {
            chain_ordering: GammaChainOrdering::RepeatedPairs,
            evaluate_traces: true,
        }
    }
}

impl GammaSimplifySettings {
    /// FORM-like chain simplification with trace evaluation enabled.
    pub fn repeated_pairs() -> Self {
        Self::default()
    }

    /// Canonical chain ordering with trace evaluation enabled.
    pub fn canonical() -> Self {
        Self {
            chain_ordering: GammaChainOrdering::Canonical,
            ..Self::default()
        }
    }

    /// Leaves `trace(...)` nodes inert after chain collection and normalization.
    pub fn without_trace_evaluation(mut self) -> Self {
        self.evaluate_traces = false;
        self
    }
}

pub(super) fn simplify_dirac_chains_impl(expr: AtomView, settings: GammaSimplifySettings) -> Atom {
    let rep = Bispinor {}.into();
    let mut expr = chainify_fixpoint(expr.to_owned().expand().simplify_metrics(), rep)
        .collect_chains(rep)
        .expand()
        .simplify_metrics();

    loop {
        let next = match (settings.chain_ordering, settings.evaluate_traces) {
            (GammaChainOrdering::RepeatedPairs, true) => {
                expr.replace_map(&dirac_chain_rewrite::<false, true>)
            }
            (GammaChainOrdering::RepeatedPairs, false) => {
                expr.replace_map(&dirac_chain_rewrite::<false, false>)
            }
            (GammaChainOrdering::Canonical, true) => {
                expr.replace_map(&dirac_chain_rewrite::<true, true>)
            }
            (GammaChainOrdering::Canonical, false) => {
                expr.replace_map(&dirac_chain_rewrite::<true, false>)
            }
        }
        .expand()
        .simplify_metrics()
        .normalize_chains()
        .simplify_metrics();

        if next == expr {
            return next;
        }

        expr = next;
    }
}

fn chainify_fixpoint(mut expr: Atom, rep: LibraryRep) -> Atom {
    loop {
        let next = expr.chainify(rep);
        if next == expr {
            return next;
        }
        expr = next;
    }
}

// `replace_map` needs a function item with sufficiently general lifetimes.
// Const generics keep the settings static without four handwritten wrappers.
fn dirac_chain_rewrite<const CANONICAL: bool, const EVALUATE_TRACES: bool>(
    arg: AtomView,
    _context: &Context,
    out: &mut Settable<'_, Atom>,
) {
    let AtomView::Fun(f) = arg else {
        return;
    };

    if f.get_symbol() == T.chain {
        if let Some(rewritten) = simplify_chain_node(arg, chain_ordering::<CANONICAL>()) {
            **out = rewritten;
        }
    } else if EVALUATE_TRACES
        && f.get_symbol() == T.trace
        && let Some(rewritten) = simplify_trace_node(arg)
    {
        **out = rewritten;
    }
}

fn chain_ordering<const CANONICAL: bool>() -> GammaChainOrdering {
    if CANONICAL {
        GammaChainOrdering::Canonical
    } else {
        GammaChainOrdering::RepeatedPairs
    }
}

fn simplify_chain_node(chain: AtomView, chain_ordering: GammaChainOrdering) -> Option<Atom> {
    let AtomView::Fun(f) = chain else {
        return None;
    };

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    let [start, end, factors @ ..] = args.as_slice() else {
        return None;
    };

    if factors.is_empty() {
        return Some(id_atom(start.clone(), end.clone()));
    }

    contract_adjacent_gamma_pair(start, end, factors).or_else(|| match chain_ordering {
        GammaChainOrdering::RepeatedPairs => {
            bubble_repeated_gamma_towards_contraction(start, end, factors)
        }
        GammaChainOrdering::Canonical => canonicalize_gamma_chain_order(start, end, factors),
    })
}

fn contract_adjacent_gamma_pair(start: &Atom, end: &Atom, factors: &[Atom]) -> Option<Atom> {
    for i in 0..factors.len().saturating_sub(1) {
        let Some(mu) = gamma_factor_lorentz(factors[i].as_view()) else {
            continue;
        };
        let Some(nu) = gamma_factor_lorentz(factors[i + 1].as_view()) else {
            continue;
        };

        if mu != nu {
            continue;
        }

        let mut rest = Vec::with_capacity(factors.len() - 2);
        rest.extend_from_slice(&factors[..i]);
        rest.extend_from_slice(&factors[i + 2..]);

        return Some(function!(ETS.metric, mu, nu) * chain!(start, end; rest));
    }

    None
}

fn bubble_repeated_gamma_towards_contraction(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
) -> Option<Atom> {
    // Repeated-pair mode only pays the anticommutation cost when it exposes a
    // contraction in the next fixed-point step.
    let (_i, j) = shortest_repeated_gamma_pair(factors)?;
    if j == 0 {
        return None;
    }

    anticommute_adjacent_gamma_pair(start, end, factors, j - 1)
}

fn canonicalize_gamma_chain_order(start: &Atom, end: &Atom, factors: &[Atom]) -> Option<Atom> {
    for i in 0..factors.len().saturating_sub(1) {
        let Some(mu) = gamma_factor_lorentz(factors[i].as_view()) else {
            continue;
        };
        let Some(nu) = gamma_factor_lorentz(factors[i + 1].as_view()) else {
            continue;
        };

        if gamma_order_key(&mu) > gamma_order_key(&nu) {
            return anticommute_adjacent_gamma_pair(start, end, factors, i);
        }
    }

    None
}

fn anticommute_adjacent_gamma_pair(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    swap_at: usize,
) -> Option<Atom> {
    let mu = gamma_factor_lorentz(factors[swap_at].as_view())?;
    let nu = gamma_factor_lorentz(factors[swap_at + 1].as_view())?;

    let mut metric_rest = Vec::with_capacity(factors.len() - 2);
    metric_rest.extend_from_slice(&factors[..swap_at]);
    metric_rest.extend_from_slice(&factors[swap_at + 2..]);

    // Try to consume the generated metric immediately by rewriting another
    // gamma factor. This keeps intermediate expressions smaller.
    let metric_term = Atom::num(2) * generated_metric_chain_term(start, end, &mu, &nu, metric_rest);

    let mut swapped = factors.to_vec();
    swapped.swap(swap_at, swap_at + 1);

    Some(metric_term - chain!(start, end; swapped))
}

fn gamma_order_key(mu: &Atom) -> String {
    mu.to_canonical_string()
}

fn generated_metric_chain_term(
    start: &Atom,
    end: &Atom,
    mu: &Atom,
    nu: &Atom,
    factors: Vec<Atom>,
) -> Atom {
    let mut contracted = factors.clone();
    if replace_first_gamma_lorentz(&mut contracted, mu, nu)
        || replace_first_gamma_lorentz(&mut contracted, nu, mu)
    {
        chain!(start, end; contracted)
    } else {
        function!(ETS.metric, mu.clone(), nu.clone()) * chain!(start, end; factors)
    }
}

fn replace_first_gamma_lorentz(factors: &mut [Atom], old: &Atom, new: &Atom) -> bool {
    for factor in factors {
        let Some(replaced) = replace_gamma_lorentz(factor.as_view(), old, new) else {
            continue;
        };

        *factor = replaced;
        return true;
    }

    false
}

fn replace_gamma_lorentz(factor: AtomView, old: &Atom, new: &Atom) -> Option<Atom> {
    let AtomView::Fun(f) = factor else {
        return None;
    };

    if f.get_symbol() != AGS.gamma || f.get_nargs() != 3 {
        return None;
    }

    let args = f.iter().collect::<Vec<_>>();
    if args[2] != old.as_view() {
        return None;
    }

    Some(
        FunctionBuilder::new(AGS.gamma)
            .add_arg(args[0])
            .add_arg(args[1])
            .add_arg(new.as_view())
            .finish(),
    )
}

fn shortest_repeated_gamma_pair(factors: &[Atom]) -> Option<(usize, usize)> {
    let mut best = None;

    for i in 0..factors.len() {
        let Some(mu) = gamma_factor_lorentz(factors[i].as_view()) else {
            continue;
        };

        for (j, factor) in factors.iter().enumerate().skip(i + 1) {
            let Some(nu) = gamma_factor_lorentz(factor.as_view()) else {
                continue;
            };

            if mu == nu && best.is_none_or(|(a, b)| j - i < b - a) {
                best = Some((i, j));
            }
        }
    }

    best
}

fn simplify_trace_node(trace: AtomView) -> Option<Atom> {
    let AtomView::Fun(f) = trace else {
        return None;
    };

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    let [rep, factors @ ..] = args.as_slice() else {
        return None;
    };

    if factors.is_empty() {
        return Some(Atom::num(4));
    }

    if factors
        .iter()
        .any(|factor| gamma_factor_lorentz(factor.as_view()).is_none())
    {
        return None;
    }

    if factors.len() % 2 == 1 {
        return Some(Atom::Zero);
    }

    let first = gamma_factor_lorentz(factors[0].as_view())?;
    let mut sum = Atom::Zero;

    // Standard recursive even trace formula:
    // tr(g1...gn) = sum_i (-1)^i g(1,i) tr(g2...g_{i-1}g_{i+1}...gn).
    for i in 1..factors.len() {
        let mu_i = gamma_factor_lorentz(factors[i].as_view())?;
        let sign = if i % 2 == 1 { 1 } else { -1 };

        let mut rest = Vec::with_capacity(factors.len() - 2);
        rest.extend_from_slice(&factors[1..i]);
        rest.extend_from_slice(&factors[i + 1..]);

        let rest_trace = if rest.is_empty() {
            Atom::num(4)
        } else {
            trace!(rep; rest)
        };

        let term = function!(ETS.metric, first.clone(), mu_i) * rest_trace;
        if sign == 1 {
            sum += term;
        } else {
            sum -= term;
        }
    }

    Some(sum)
}

fn gamma_factor_lorentz(factor: AtomView) -> Option<Atom> {
    let AtomView::Fun(f) = factor else {
        return None;
    };

    if f.get_symbol() != AGS.gamma || f.get_nargs() != 3 {
        return None;
    }

    let args = f.iter().collect::<Vec<_>>();
    if !(is_chain_gamma_endpoint(args[0], T.chain_in)
        && is_chain_gamma_endpoint(args[1], T.chain_out)
        || is_chain_gamma_endpoint(args[0], T.chain_out)
            && is_chain_gamma_endpoint(args[1], T.chain_in))
    {
        return None;
    }

    Some(args[2].to_owned())
}

fn is_chain_gamma_endpoint(arg: AtomView, expected: Symbol) -> bool {
    matches!(arg, AtomView::Var(symbol) if symbol.get_symbol() == expected)
}
