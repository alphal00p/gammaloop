use std::sync::LazyLock;

use spenso::{
    chain,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    structure::representation::{Minkowski, RepName},
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Context,
    utils::Settable,
};

use crate::{
    W_,
    chain::Chain,
    metric::MetricSimplifier,
    representations::Bispinor,
    schoonschip::{Schoonschip, SchoonschipSettings},
};

use super::{AGS, id_atom};

static MINKOWSKI_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
    let minkowski = Minkowski {}.to_symbolic(std::iter::empty::<Atom>());
    let AtomView::Fun(f) = minkowski.as_view() else {
        unreachable!("Minkowski representations are symbolic functions")
    };

    f.get_symbol()
});

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DiracRuleDimension {
    AnyDimension,
    FourDimensional,
}

const ADJACENT_GAMMA_CONTRACTION: DiracRuleDimension = DiracRuleDimension::AnyDimension;
const GAMMA_ANTICOMMUTATION: DiracRuleDimension = DiracRuleDimension::AnyDimension;
const FOUR_DIM_CHISHOLM: DiracRuleDimension = DiracRuleDimension::FourDimensional;
const TRACE_GAMMA_RECURSION: DiracRuleDimension = DiracRuleDimension::AnyDimension;

#[derive(Debug, Clone, PartialEq, Eq)]
enum DiracFactor {
    Gamma {
        lorentz: Atom,
        dimension: Option<Atom>,
    },
    Gamma5,
    Gamma0,
    ProjectorPlus,
    ProjectorMinus,
    Other(Atom),
}

impl DiracFactor {
    fn parse(factor: AtomView) -> Self {
        let AtomView::Fun(f) = factor else {
            return Self::Other(factor.to_owned());
        };

        if f.get_symbol() == AGS.gamma && f.get_nargs() == 3 {
            let args = f.iter().collect::<Vec<_>>();
            if has_chain_endpoints(args[0], args[1]) {
                let lorentz = args[2].to_owned();
                let dimension = gamma_lorentz_dimension(lorentz.as_view());
                return Self::Gamma { lorentz, dimension };
            }
        }

        if f.get_nargs() == 2 {
            let args = f.iter().collect::<Vec<_>>();
            if has_chain_endpoints(args[0], args[1]) {
                return match f.get_symbol() {
                    symbol if symbol == AGS.gamma5 => Self::Gamma5,
                    symbol if symbol == AGS.gamma0 => Self::Gamma0,
                    symbol if symbol == AGS.projp => Self::ProjectorPlus,
                    symbol if symbol == AGS.projm => Self::ProjectorMinus,
                    _ => Self::Other(factor.to_owned()),
                };
            }
        }

        Self::Other(factor.to_owned())
    }

    fn gamma_lorentz(&self, rule_dimension: DiracRuleDimension) -> Option<&Atom> {
        match self {
            Self::Gamma { lorentz, dimension }
                if rule_dimension.allows_gamma_dimension(dimension.as_ref()) =>
            {
                Some(lorentz)
            }
            Self::Other(atom) => {
                let _ = atom;
                None
            }
            _ => None,
        }
    }

    fn gamma_dimension(&self) -> Option<&Atom> {
        match self {
            Self::Gamma { dimension, .. } => dimension.as_ref(),
            Self::Other(atom) => {
                let _ = atom;
                None
            }
            _ => None,
        }
    }
}

impl DiracRuleDimension {
    fn allows_gamma_dimension(self, dimension: Option<&Atom>) -> bool {
        match self {
            Self::AnyDimension => true,
            Self::FourDimensional => dimension.is_some_and(is_four_dimension),
        }
    }
}

pub(super) fn simplify_dirac_chains_impl(expr: AtomView, settings: GammaSimplifySettings) -> Atom {
    let rep = Bispinor {}.into();
    let mut expr = expr
        .to_owned()
        .expand()
        .simplify_metrics()
        .chainify(rep)
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

    let parsed_factors = factors
        .iter()
        .map(|factor| DiracFactor::parse(factor.as_view()))
        .collect::<Vec<_>>();

    contract_adjacent_gamma_pair(start, end, factors, &parsed_factors)
        .or_else(|| four_dim_chisholm_contraction(start, end, factors, &parsed_factors))
        .or_else(|| match chain_ordering {
            GammaChainOrdering::RepeatedPairs => {
                bubble_repeated_gamma_towards_contraction(start, end, factors, &parsed_factors)
            }
            GammaChainOrdering::Canonical => {
                canonicalize_gamma_chain_order(start, end, factors, &parsed_factors)
            }
        })
}

fn contract_adjacent_gamma_pair(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    parsed_factors: &[DiracFactor],
) -> Option<Atom> {
    for i in 0..factors.len().saturating_sub(1) {
        let Some((mu, nu)) = gamma_pair_lorentz_for(
            ADJACENT_GAMMA_CONTRACTION,
            &parsed_factors[i],
            &parsed_factors[i + 1],
        ) else {
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
    parsed_factors: &[DiracFactor],
) -> Option<Atom> {
    // Repeated-pair mode only pays the anticommutation cost when it exposes a
    // contraction in the next fixed-point step.
    let (_i, j) = shortest_repeated_gamma_pair(parsed_factors)?;
    if j == 0 {
        return None;
    }

    anticommute_adjacent_gamma_pair(start, end, factors, parsed_factors, j - 1)
}

fn four_dim_chisholm_contraction(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    parsed_factors: &[DiracFactor],
) -> Option<Atom> {
    let (left, right) = shortest_repeated_four_dim_gamma_pair(parsed_factors)?;
    let interior = &factors[left + 1..right];
    let parsed_interior = &parsed_factors[left + 1..right];
    let interior_lorentz = gamma_lorentz_sequence_for(FOUR_DIM_CHISHOLM, parsed_interior)?;

    // FORM and FeynCalc use direct 4D Chisholm identities for short interiors;
    // this avoids generating the larger anticommutation expansion first.
    match interior.len() {
        1 => Some(
            Atom::num(-2)
                * chain!(start, end; chain_factors(factors, left, right, interior.iter().cloned())),
        ),
        2 => {
            let mu = interior_lorentz[0].clone();
            let nu = interior_lorentz[1].clone();

            Some(
                Atom::num(4)
                    * function!(ETS.metric, mu, nu)
                    * chain!(start, end; chain_factors(factors, left, right, [])),
            )
        }
        3 => {
            let reversed = interior.iter().rev().cloned().collect::<Vec<_>>();
            Some(Atom::num(-2) * chain!(start, end; chain_factors(factors, left, right, reversed)))
        }
        4 => {
            let term_1 = [
                interior[2].clone(),
                interior[1].clone(),
                interior[0].clone(),
                interior[3].clone(),
            ];
            let term_2 = [
                interior[3].clone(),
                interior[0].clone(),
                interior[1].clone(),
                interior[2].clone(),
            ];

            Some(
                Atom::num(2) * chain!(start, end; chain_factors(factors, left, right, term_1))
                    + Atom::num(2)
                        * chain!(start, end; chain_factors(factors, left, right, term_2)),
            )
        }
        _ => None,
    }
}

fn chain_factors(
    factors: &[Atom],
    left: usize,
    right: usize,
    middle: impl IntoIterator<Item = Atom>,
) -> Vec<Atom> {
    let mut result = Vec::with_capacity(factors.len() - 2);
    result.extend_from_slice(&factors[..left]);
    result.extend(middle);
    result.extend_from_slice(&factors[right + 1..]);
    result
}

fn canonicalize_gamma_chain_order(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    parsed_factors: &[DiracFactor],
) -> Option<Atom> {
    for i in 0..factors.len().saturating_sub(1) {
        let Some((mu, nu)) = gamma_pair_lorentz_for(
            GAMMA_ANTICOMMUTATION,
            &parsed_factors[i],
            &parsed_factors[i + 1],
        ) else {
            continue;
        };

        if mu > nu {
            return anticommute_adjacent_gamma_pair(start, end, factors, parsed_factors, i);
        }
    }

    None
}

fn anticommute_adjacent_gamma_pair(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    parsed_factors: &[DiracFactor],
    swap_at: usize,
) -> Option<Atom> {
    let (mu, nu) = gamma_pair_lorentz_for(
        GAMMA_ANTICOMMUTATION,
        &parsed_factors[swap_at],
        &parsed_factors[swap_at + 1],
    )?;

    let mut metric_rest = Vec::with_capacity(factors.len() - 2);
    metric_rest.extend_from_slice(&factors[..swap_at]);
    metric_rest.extend_from_slice(&factors[swap_at + 2..]);

    // Run the local metric term through chain-aware Schoonschip before it can
    // swell the Clifford expansion.
    let metric_term = Atom::num(2) * generated_metric_chain_term(start, end, &mu, &nu, metric_rest);

    let mut swapped = factors.to_vec();
    swapped.swap(swap_at, swap_at + 1);

    Some(metric_term - chain!(start, end; swapped))
}

fn generated_metric_chain_term(
    start: &Atom,
    end: &Atom,
    mu: &Atom,
    nu: &Atom,
    factors: Vec<Atom>,
) -> Atom {
    (function!(ETS.metric, mu.clone(), nu.clone()) * chain!(start, end; factors))
        .schoonschip_with_settings(
            &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
        )
}

fn shortest_repeated_gamma_pair(factors: &[DiracFactor]) -> Option<(usize, usize)> {
    let mut best = None;

    for i in 0..factors.len() {
        if factors[i].gamma_lorentz(GAMMA_ANTICOMMUTATION).is_none() {
            continue;
        }

        for (j, factor) in factors.iter().enumerate().skip(i + 1) {
            let Some((mu, nu)) = gamma_pair_lorentz_for(GAMMA_ANTICOMMUTATION, &factors[i], factor)
            else {
                continue;
            };

            if mu == nu && best.is_none_or(|(a, b)| j - i < b - a) {
                best = Some((i, j));
            }
        }
    }

    best
}

fn shortest_repeated_four_dim_gamma_pair(factors: &[DiracFactor]) -> Option<(usize, usize)> {
    let mut best = None;

    for i in 0..factors.len() {
        let Some(mu) = factors[i].gamma_lorentz(FOUR_DIM_CHISHOLM) else {
            continue;
        };
        if !is_minkowski_slot(mu) {
            continue;
        }

        for (j, factor) in factors.iter().enumerate().skip(i + 1) {
            let Some(nu) = factor.gamma_lorentz(FOUR_DIM_CHISHOLM) else {
                continue;
            };

            if mu == nu && best.is_none_or(|(a, b)| j - i < b - a) {
                best = Some((i, j));
            }
        }
    }

    best
}

fn gamma_pair_lorentz_for(
    rule_dimension: DiracRuleDimension,
    left: &DiracFactor,
    right: &DiracFactor,
) -> Option<(Atom, Atom)> {
    let left_lorentz = left.gamma_lorentz(rule_dimension)?;
    let right_lorentz = right.gamma_lorentz(rule_dimension)?;
    if !compatible_gamma_dimensions(rule_dimension, left, right) {
        return None;
    }

    Some((left_lorentz.clone(), right_lorentz.clone()))
}

fn gamma_lorentz_sequence_for(
    rule_dimension: DiracRuleDimension,
    factors: &[DiracFactor],
) -> Option<Vec<Atom>> {
    let mut known_dimension = None;
    let mut lorentz = Vec::with_capacity(factors.len());

    for factor in factors {
        let gamma_lorentz = factor.gamma_lorentz(rule_dimension)?;
        if !merge_known_gamma_dimension(
            rule_dimension,
            &mut known_dimension,
            factor.gamma_dimension(),
        ) {
            return None;
        }
        lorentz.push(gamma_lorentz.clone());
    }

    Some(lorentz)
}

fn compatible_gamma_dimensions(
    rule_dimension: DiracRuleDimension,
    left: &DiracFactor,
    right: &DiracFactor,
) -> bool {
    match rule_dimension {
        DiracRuleDimension::AnyDimension => match (left.gamma_dimension(), right.gamma_dimension())
        {
            (Some(left), Some(right)) => left == right,
            _ => true,
        },
        DiracRuleDimension::FourDimensional => true,
    }
}

fn merge_known_gamma_dimension(
    rule_dimension: DiracRuleDimension,
    known_dimension: &mut Option<Atom>,
    dimension: Option<&Atom>,
) -> bool {
    if rule_dimension == DiracRuleDimension::FourDimensional {
        return true;
    }

    let Some(dimension) = dimension else {
        return true;
    };

    match known_dimension {
        Some(known_dimension) => known_dimension == dimension,
        None => {
            *known_dimension = Some(dimension.clone());
            true
        }
    }
}

fn gamma_lorentz_dimension(lorentz: AtomView) -> Option<Atom> {
    if let Some(dimension) = minkowski_dimension(lorentz) {
        return Some(dimension);
    }

    let AtomView::Fun(f) = lorentz else {
        return None;
    };

    f.iter().find_map(minkowski_dimension)
}

fn minkowski_dimension(atom: AtomView) -> Option<Atom> {
    let AtomView::Fun(f) = atom else {
        return None;
    };

    if f.get_symbol() != *MINKOWSKI_SYMBOL || f.get_nargs() == 0 {
        return None;
    }

    f.iter().next().map(|dimension| dimension.to_owned())
}

fn is_minkowski_slot(atom: &Atom) -> bool {
    let AtomView::Fun(f) = atom.as_view() else {
        return false;
    };

    f.get_symbol() == *MINKOWSKI_SYMBOL && f.get_nargs() == 2
}

fn is_four_dimension(dimension: &Atom) -> bool {
    *dimension == Atom::num(4)
}

fn has_chain_endpoints(left: AtomView, right: AtomView) -> bool {
    is_chain_endpoint(left, T.chain_in) && is_chain_endpoint(right, T.chain_out)
        || is_chain_endpoint(left, T.chain_out) && is_chain_endpoint(right, T.chain_in)
}

fn is_chain_endpoint(arg: AtomView, expected: Symbol) -> bool {
    matches!(arg, AtomView::Var(symbol) if symbol.get_symbol() == expected)
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
        return simplify_trace_terminal(trace);
    }

    let parsed_factors = factors
        .iter()
        .map(|factor| DiracFactor::parse(factor.as_view()))
        .collect::<Vec<_>>();
    let trace_lorentz = gamma_lorentz_sequence_for(TRACE_GAMMA_RECURSION, &parsed_factors)?;

    if factors.len() % 2 == 1 {
        return Some(Atom::Zero);
    }

    let first = trace_lorentz[0].clone();
    let mut sum = Atom::Zero;

    // Standard recursive even trace formula:
    // tr(g1...gn) = sum_i (-1)^i g(1,i) tr(g2...g_{i-1}g_{i+1}...gn).
    for i in 1..factors.len() {
        let mu_i = trace_lorentz[i].clone();
        let sign = if i % 2 == 1 { 1 } else { -1 };

        let mut rest = Vec::with_capacity(factors.len() - 2);
        rest.extend_from_slice(&factors[1..i]);
        rest.extend_from_slice(&factors[i + 1..]);

        let rest_trace = if rest.is_empty() {
            let terminal_trace = trace!(rep; std::iter::empty::<Atom>());
            simplify_trace_terminal(terminal_trace.as_view())?
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

fn simplify_trace_terminal(trace: AtomView) -> Option<Atom> {
    let trace = trace.to_owned();
    let simplified = trace
        .replace(function!(T.trace, T.rep_::<0, _>([W_.d_])).to_pattern())
        .with(Atom::var(W_.d_));
    (simplified != trace).then_some(simplified)
}
