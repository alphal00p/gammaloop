use std::sync::LazyLock;

use spenso::{
    chain,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    rep_,
    structure::representation::{LibraryRep, Minkowski, RepName},
    symbolica_atom::IntoAtom,
    tensors::parametric::atomcore::PatternReplacement,
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, representation::FunView},
    function,
    id::{Context, Replacement},
    utils::Settable,
};

use crate::{
    W_,
    chain::Chain,
    epsilon::{EpsilonSimplifier, epsilon4},
    metric::MetricSimplifier,
    representations::Bispinor,
    schoonschip::{Schoonschip, SchoonschipSettings},
};

use super::{AGS, id_atom};

static MINKOWSKI_SYMBOL: LazyLock<Symbol> =
    LazyLock::new(|| LibraryRep::from(Minkowski {}).symbol());

static BISPINOR_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| LibraryRep::from(Bispinor {}).symbol());

static EPSILON_DUMMY_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| symbolica::symbol!("sigma"));

static TRACE_TERMINALS: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    [Replacement::new(
        trace!(rep_!(0; W_.d_)).to_pattern(),
        Atom::var(W_.d_),
    )]
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
    /// Whether three 4D gammas may be expanded into the gamma5-epsilon basis.
    pub expand_three_gamma_epsilon: bool,
}

impl Default for GammaSimplifySettings {
    fn default() -> Self {
        Self {
            chain_ordering: GammaChainOrdering::RepeatedPairs,
            evaluate_traces: true,
            expand_three_gamma_epsilon: false,
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

    /// Enables the four-dimensional identity that rewrites three gammas into
    /// metric terms plus a gamma5-epsilon term.
    pub fn with_gamma5_epsilon_expansion(mut self) -> Self {
        self.expand_three_gamma_epsilon = true;
        self
    }

    fn rewrite_expression(&self, expr: Atom) -> Atom {
        expr.replace_map(|a, b, c| self.rewrite_node(a, b, c))
    }

    fn rewrite_node(&self, arg: AtomView, _context: &Context, out: &mut Settable<'_, Atom>) {
        let AtomView::Fun(f) = arg else {
            return;
        };

        let simplifier = DiracSimplifier::new(self);
        if f.get_symbol() == T.chain {
            if let Some(rewritten) = simplifier.simplify_chain_node(f) {
                **out = rewritten;
            }
        } else if self.evaluate_traces
            && f.get_symbol() == T.trace
            && let Some(rewritten) = simplifier.simplify_trace_node(f)
        {
            **out = rewritten;
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DiracRuleDimension {
    /// The rule is valid in arbitrary dimension.
    ///
    /// Pair and sequence helpers still compare dimensions when both gamma
    /// factors expose one. If a dimension cannot be inferred, the rule stays
    /// permissive and leaves the caller's structural index checks to decide
    /// whether a contraction or ordering step is valid.
    AnyDimension,
    /// The rule is intrinsically four-dimensional.
    ///
    /// Gamma factors are admitted only when their Minkowski index exposes a
    /// four-dimensional Minkowski representation. Once admitted, no additional
    /// pairwise dimension comparison is needed for that rule.
    FourDimensional,
}

const ADJACENT_GAMMA_CONTRACTION: DiracRuleDimension = DiracRuleDimension::AnyDimension;
const GAMMA_ANTICOMMUTATION: DiracRuleDimension = DiracRuleDimension::AnyDimension;
const FOUR_DIM_CHISHOLM: DiracRuleDimension = DiracRuleDimension::FourDimensional;
const FOUR_DIM_GAMMA5_ANTICOMMUTATION: DiracRuleDimension = DiracRuleDimension::FourDimensional;
const FOUR_DIM_THREE_GAMMA_EPSILON: DiracRuleDimension = DiracRuleDimension::FourDimensional;
const TRACE_GAMMA_RECURSION: DiracRuleDimension = DiracRuleDimension::AnyDimension;
const TRACE_GAMMA5_RECURSION: DiracRuleDimension = DiracRuleDimension::FourDimensional;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DiracFactor<'a> {
    /// A chain-local gamma factor, borrowing its Minkowski index and inferred
    /// Minkowski dimension from the source expression.
    Gamma {
        factor: AtomView<'a>,
        mink_index: AtomView<'a>,
        dimension: Option<AtomView<'a>>,
    },
    Gamma5(AtomView<'a>),
    Gamma0(AtomView<'a>),
    ProjectorPlus(AtomView<'a>),
    ProjectorMinus(AtomView<'a>),
    Other(AtomView<'a>),
}

impl<'a> DiracFactor<'a> {
    /// Classifies a factor inside a `chain(...)` or `trace(...)` without
    /// materializing owned atoms.
    ///
    /// Ordinary gammas are recognized as `gamma(in,out,mu)` or
    /// `gamma(out,in,mu)` and keep `mu` by view. Their dimension is inferred
    /// from `mu` when it is a Minkowski slot, or from the first visible
    /// Minkowski representation inside slash-like tensorial indices such as
    /// `P(1,mink(D))`.
    fn parse(factor: AtomView<'a>) -> Self {
        let AtomView::Fun(f) = factor else {
            return Self::Other(factor);
        };

        if f.get_symbol() == AGS.gamma && f.get_nargs() == 3 {
            let mut args = f.iter();
            let (Some(left), Some(right), Some(mink_index)) =
                (args.next(), args.next(), args.next())
            else {
                return Self::Other(factor);
            };
            if has_chain_endpoints(left, right) {
                let dimension = mink_slot_dimension(mink_index);
                return Self::Gamma {
                    factor,
                    mink_index,
                    dimension,
                };
            }
        }

        if f.get_nargs() == 2 {
            let mut args = f.iter();
            let (Some(left), Some(right)) = (args.next(), args.next()) else {
                return Self::Other(factor);
            };
            if has_chain_endpoints(left, right) {
                return match f.get_symbol() {
                    symbol if symbol == AGS.gamma5 => Self::Gamma5(factor),
                    symbol if symbol == AGS.gamma0 => Self::Gamma0(factor),
                    symbol if symbol == AGS.projp => Self::ProjectorPlus(factor),
                    symbol if symbol == AGS.projm => Self::ProjectorMinus(factor),
                    _ => Self::Other(factor),
                };
            }
        }

        Self::Other(factor)
    }

    /// Returns the original factor view used to build this parsed factor.
    fn as_view(self) -> AtomView<'a> {
        match self {
            Self::Gamma { factor, .. }
            | Self::Gamma5(factor)
            | Self::Gamma0(factor)
            | Self::ProjectorPlus(factor)
            | Self::ProjectorMinus(factor)
            | Self::Other(factor) => factor,
        }
    }

    /// Returns the Minkowski index when this factor is a gamma admitted by the
    /// requested rule dimension.
    fn gamma_mink_index(self, rule_dimension: DiracRuleDimension) -> Option<AtomView<'a>> {
        match self {
            Self::Gamma {
                mink_index,
                dimension,
                ..
            } if rule_dimension.allows_gamma_dimension(dimension) => Some(mink_index),
            _ => None,
        }
    }

    /// Returns the dimension inferred from this gamma's Minkowski index, if
    /// syntactically visible.
    fn gamma_dimension(self) -> Option<AtomView<'a>> {
        match self {
            Self::Gamma { dimension, .. } => dimension,
            _ => None,
        }
    }

    fn is_gamma5(self) -> bool {
        matches!(self, Self::Gamma5(_))
    }

    fn anticommutes_with_gamma5(self) -> bool {
        matches!(self, Self::Gamma { .. } | Self::Gamma0(_))
    }
}

impl DiracRuleDimension {
    /// Checks the per-factor dimension gate for a rule.
    ///
    /// Arbitrary-dimensional rules do not reject an individual gamma here; they
    /// rely on the pair/sequence compatibility checks below. Four-dimensional
    /// rules require an explicitly inferred dimension equal to integer `4`.
    fn allows_gamma_dimension(self, dimension: Option<AtomView<'_>>) -> bool {
        match self {
            Self::AnyDimension => true,
            Self::FourDimensional => dimension.is_some_and(is_four_dimension),
        }
    }

    /// Checks whether two gamma factors have dimensions compatible with this
    /// rule.
    ///
    /// For arbitrary-dimensional rules, visible dimensions must match on both
    /// sides. Missing dimensions are treated as unknown rather than
    /// contradictory; the caller still compares the actual index atoms when it
    /// needs a repeated-index contraction. Four-dimensional rules do not compare
    /// here because each gamma was already required to expose dimension `4`.
    fn gamma_compatible(self, left: DiracFactor<'_>, right: DiracFactor<'_>) -> bool {
        match self {
            Self::AnyDimension => match (left.gamma_dimension(), right.gamma_dimension()) {
                (Some(left), Some(right)) => left == right,
                _ => true,
            },
            Self::FourDimensional => true,
        }
    }

    /// Extends the known dimension for a gamma sequence under this rule.
    ///
    /// Arbitrary-dimensional rules reject mixed known dimensions, but tolerate
    /// unknown dimensions. Four-dimensional rules skip this because per-factor
    /// admission has already enforced dimension `4`.
    fn merge_known_gamma_dimension<'a>(
        self,
        known_dimension: &mut Option<AtomView<'a>>,
        dimension: Option<AtomView<'a>>,
    ) -> bool {
        if self == Self::FourDimensional {
            return true;
        }

        let Some(dimension) = dimension else {
            return true;
        };

        match known_dimension {
            Some(known_dimension) => *known_dimension == dimension,
            None => {
                *known_dimension = Some(dimension);
                true
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct DiracSimplifier<'settings> {
    settings: &'settings GammaSimplifySettings,
}

impl<'settings> DiracSimplifier<'settings> {
    pub(crate) fn new(settings: &'settings GammaSimplifySettings) -> Self {
        Self { settings }
    }

    pub(crate) fn simplify(self, expr: AtomView) -> Atom {
        let rep = Bispinor {}.into();
        let mut expr = expr
            .to_owned()
            .expand()
            .simplify_metrics()
            .simplify_epsilon()
            .chainify(rep)
            .collect_chains(rep)
            .expand()
            .simplify_metrics()
            .simplify_epsilon();

        loop {
            let next = self
                .settings
                .rewrite_expression(expr.clone())
                .expand()
                .simplify_metrics()
                .simplify_epsilon()
                .normalize_chains()
                .simplify_metrics()
                .simplify_epsilon();

            if next == expr {
                return next;
            }

            expr = next;
        }
    }

    fn simplify_chain_node(self, f: FunView) -> Option<Atom> {
        let args = f.iter().collect::<Vec<_>>();
        let [start, end, factors @ ..] = args.as_slice() else {
            return None;
        };

        let factors = factors
            .iter()
            .map(|factor| DiracFactor::parse(*factor))
            .collect::<Vec<_>>();

        if factors.is_empty() {
            return Some(id_atom(*start, *end));
        }

        Self::contract_adjacent_gamma_pair(*start, *end, &factors)
            .or_else(|| Self::contract_adjacent_special_dirac_pair(*start, *end, &factors))
            .or_else(|| Self::conjugate_special_factor_by_gamma0(*start, *end, &factors))
            .or_else(|| Self::four_dim_chisholm_contraction(*start, *end, &factors))
            .or_else(|| {
                self.settings
                    .expand_three_gamma_epsilon
                    .then(|| Self::four_dim_three_gamma_epsilon_expansion(*start, *end, &factors))
                    .flatten()
            })
            .or_else(|| Self::move_gamma5_right_of_gamma(*start, *end, &factors))
            .or_else(|| Self::move_projector_right_of_gamma(*start, *end, &factors))
            .or_else(|| match self.settings.chain_ordering {
                GammaChainOrdering::RepeatedPairs => {
                    Self::bubble_repeated_gamma_towards_contraction(*start, *end, &factors)
                }
                GammaChainOrdering::Canonical => {
                    Self::canonicalize_gamma_chain_order(*start, *end, &factors)
                }
            })
    }
}

impl DiracSimplifier<'_> {
    /// Contracts adjacent equal-index gammas:
    /// `...[gamma(mu), gamma(mu)]... -> g(mu, mu) * ...[...]...`.
    fn contract_adjacent_gamma_pair(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        for i in 0..factors.len().saturating_sub(1) {
            let Some((mu, nu)) =
                Self::mink_index_pair(ADJACENT_GAMMA_CONTRACTION, &factors[i], &factors[i + 1])
            else {
                continue;
            };

            if mu != nu {
                continue;
            }

            let mut rest = Vec::with_capacity(factors.len() - 2);
            Self::extend_factors(&mut rest, &factors[..i]);
            Self::extend_factors(&mut rest, &factors[i + 2..]);

            return Some(function!(ETS.metric, mu, nu) * chain!(start, end; rest));
        }

        None
    }

    /// Reduces adjacent special 4D factors:
    /// `gamma5 gamma5 -> 1`, `gamma0 gamma0 -> 1`, `P+ P+ -> P+`,
    /// `P- P- -> P-`, and `P+ P- -> 0`.
    fn contract_adjacent_special_dirac_pair(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        if !has_four_dimensional_spin_endpoints(start, end) {
            return None;
        }

        for i in 0..factors.len().saturating_sub(1) {
            let replacement = match (&factors[i], &factors[i + 1]) {
                (DiracFactor::Gamma5(_), DiracFactor::Gamma5(_))
                | (DiracFactor::Gamma0(_), DiracFactor::Gamma0(_)) => Vec::new(),
                (DiracFactor::ProjectorPlus(_), DiracFactor::ProjectorPlus(_)) => {
                    vec![endpoint_factor(AGS.projp)]
                }
                (DiracFactor::ProjectorMinus(_), DiracFactor::ProjectorMinus(_)) => {
                    vec![endpoint_factor(AGS.projm)]
                }
                (DiracFactor::ProjectorPlus(_), DiracFactor::ProjectorMinus(_))
                | (DiracFactor::ProjectorMinus(_), DiracFactor::ProjectorPlus(_)) => {
                    return Some(Atom::Zero);
                }
                _ => continue,
            };

            return Some(chain!(start, end; Self::chain_factors(
                factors,
                i,
                i + 1,
                replacement,
            )));
        }

        None
    }

    /// Conjugates a special 4D factor by `gamma0`:
    /// `gamma0 gamma5 gamma0 -> -gamma5`,
    /// `gamma0 P+ gamma0 -> P-`, and `gamma0 P- gamma0 -> P+`.
    fn conjugate_special_factor_by_gamma0(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        if !has_four_dimensional_spin_endpoints(start, end) {
            return None;
        }

        for i in 0..factors.len().saturating_sub(2) {
            if !matches!(factors[i], DiracFactor::Gamma0(_))
                || !matches!(factors[i + 2], DiracFactor::Gamma0(_))
            {
                continue;
            }

            let (sign, replacement) = match factors[i + 1] {
                DiracFactor::Gamma5(_) => (-1, gamma5_factor()),
                DiracFactor::ProjectorPlus(_) => (1, endpoint_factor(AGS.projm)),
                DiracFactor::ProjectorMinus(_) => (1, endpoint_factor(AGS.projp)),
                _ => continue,
            };

            let rewritten =
                chain!(start, end; Self::chain_factors(factors, i, i + 2, [replacement]));
            return Some(if sign == 1 {
                rewritten
            } else {
                Atom::num(sign) * rewritten
            });
        }

        None
    }

    /// Moves repeated gammas toward each other with
    /// `gamma(mu) gamma(nu) -> 2 g(mu,nu) - gamma(nu) gamma(mu)`.
    ///
    /// Repeated-pair mode only pays the anticommutation cost when it exposes a
    /// contraction in the next fixed-point step.
    fn bubble_repeated_gamma_towards_contraction(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        let (_i, j) = Self::shortest_repeated_gamma_pair(factors)?;
        if j == 0 {
            return None;
        }

        Self::anticommute_adjacent_gamma_pair(start, end, factors, j - 1)
    }

    /// Applies short 4D Chisholm contractions around repeated endpoint gammas:
    /// `gamma(mu) A gamma(mu)` is replaced directly for interiors of length
    /// one through four.
    fn four_dim_chisholm_contraction(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        let (left, right) = Self::shortest_repeated_four_dim_gamma_pair(factors)?;
        let parsed_interior = &factors[left + 1..right];
        let interior_mink_indices =
            Self::gamma_mink_index_sequence_for(FOUR_DIM_CHISHOLM, parsed_interior)?;

        match parsed_interior.len() {
            1 => Some(
                Atom::num(-2)
                    * chain!(start, end; Self::chain_factors(
                        factors,
                        left,
                        right,
                        parsed_interior.iter().copied().map(DiracFactor::as_view),
                    )),
            ),
            2 => {
                let mu = interior_mink_indices[0];
                let nu = interior_mink_indices[1];

                Some(
                    Atom::num(4)
                        * function!(ETS.metric, mu, nu)
                        * chain!(start, end; Self::chain_factors(
                            factors,
                            left,
                            right,
                            std::iter::empty::<Atom>(),
                        )),
                )
            }
            3 => {
                let reversed = parsed_interior
                    .iter()
                    .rev()
                    .copied()
                    .map(DiracFactor::as_view)
                    .collect::<Vec<_>>();
                Some(
                    Atom::num(-2)
                        * chain!(start, end; Self::chain_factors(
                            factors,
                            left,
                            right,
                            reversed,
                        )),
                )
            }
            4 => {
                let term_1 = [
                    parsed_interior[2].as_view(),
                    parsed_interior[1].as_view(),
                    parsed_interior[0].as_view(),
                    parsed_interior[3].as_view(),
                ];
                let term_2 = [
                    parsed_interior[3].as_view(),
                    parsed_interior[0].as_view(),
                    parsed_interior[1].as_view(),
                    parsed_interior[2].as_view(),
                ];

                Some(
                    Atom::num(2)
                        * chain!(start, end; Self::chain_factors(
                            factors,
                            left,
                            right,
                            term_1,
                        ))
                        + Atom::num(2)
                            * chain!(start, end; Self::chain_factors(
                                factors,
                                left,
                                right,
                                term_2,
                            )),
                )
            }
            _ => None,
        }
    }

    /// Expands three 4D gammas into metric terms plus an epsilon-gamma5 term:
    /// `gamma(mu) gamma(nu) gamma(rho) -> g(mu,nu) gamma(rho)
    /// - g(mu,rho) gamma(nu) + g(nu,rho) gamma(mu)
    /// - epsilon(mu,nu,rho,sigma) gamma(sigma) gamma5`.
    fn four_dim_three_gamma_epsilon_expansion(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        if !has_four_dimensional_spin_endpoints(start, end) {
            return None;
        }

        for i in 0..factors.len().saturating_sub(2) {
            let Some(mink_indices) = Self::gamma_mink_index_sequence_for(
                FOUR_DIM_THREE_GAMMA_EPSILON,
                &factors[i..i + 3],
            ) else {
                continue;
            };
            if !mink_indices.iter().copied().all(is_minkowski_slot) {
                continue;
            }

            let [mu, nu, rho] = mink_indices.as_slice() else {
                unreachable!("the window always contains three gamma factors")
            };
            let sigma = epsilon_dummy_minkowski_slot();

            // The chain normal form keeps gamma5 to the right of ordinary gammas.
            let epsilon_term = Atom::num(-1)
                * chain!(start, end; Self::chain_factors(
                    factors,
                    i,
                    i + 2,
                    [gamma_factor(sigma.clone()), gamma5_factor()],
                ))
                * epsilon4(*mu, *nu, *rho, &sigma);

            let metric_mu_nu = Self::generated_metric_chain_term(
                start,
                end,
                *mu,
                *nu,
                Self::chain_factors(factors, i, i + 2, [factors[i + 2].as_view()]),
            );
            let metric_mu_rho = Self::generated_metric_chain_term(
                start,
                end,
                *mu,
                *rho,
                Self::chain_factors(factors, i, i + 2, [factors[i + 1].as_view()]),
            );
            let metric_nu_rho = Self::generated_metric_chain_term(
                start,
                end,
                *nu,
                *rho,
                Self::chain_factors(factors, i, i + 2, [factors[i].as_view()]),
            );

            return Some(epsilon_term + metric_mu_nu - metric_mu_rho + metric_nu_rho);
        }

        None
    }

    fn chain_factors<M: IntoAtom>(
        factors: &[DiracFactor<'_>],
        left: usize,
        right: usize,
        middle: impl IntoIterator<Item = M>,
    ) -> Vec<Atom> {
        let mut result = Vec::with_capacity(factors.len() - 2);
        Self::extend_factors(&mut result, &factors[..left]);
        result.extend(middle.into_iter().map(IntoAtom::into_atom));
        Self::extend_factors(&mut result, &factors[right + 1..]);
        result
    }

    fn owned_factors(factors: &[DiracFactor<'_>]) -> Vec<Atom> {
        let mut result = Vec::with_capacity(factors.len());
        Self::extend_factors(&mut result, factors);
        result
    }

    fn extend_factors(result: &mut Vec<Atom>, factors: &[DiracFactor<'_>]) {
        result.extend(
            factors
                .iter()
                .copied()
                .map(DiracFactor::as_view)
                .map(IntoAtom::into_atom),
        );
    }

    /// Canonicalizes adjacent gamma order using the Clifford anticommutator:
    /// `gamma(mu) gamma(nu) -> 2 g(mu,nu) - gamma(nu) gamma(mu)`.
    fn canonicalize_gamma_chain_order(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        for i in 0..factors.len().saturating_sub(1) {
            let Some((mu, nu)) =
                Self::mink_index_pair(GAMMA_ANTICOMMUTATION, &factors[i], &factors[i + 1])
            else {
                continue;
            };

            if mu > nu {
                return Self::anticommute_adjacent_gamma_pair(start, end, factors, i);
            }
        }

        None
    }

    /// Moves `gamma5` to the right of a 4D gamma:
    /// `gamma5 gamma(mu) -> -gamma(mu) gamma5`.
    fn move_gamma5_right_of_gamma(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        if !has_four_dimensional_spin_endpoints(start, end) {
            return None;
        }

        for i in 0..factors.len().saturating_sub(1) {
            let DiracFactor::Gamma5(_) = factors[i] else {
                continue;
            };
            if factors[i + 1]
                .gamma_mink_index(FOUR_DIM_GAMMA5_ANTICOMMUTATION)
                .is_none()
            {
                continue;
            }

            let mut swapped = Self::owned_factors(factors);
            swapped.swap(i, i + 1);
            return Some(Atom::num(-1) * chain!(start, end; swapped));
        }

        None
    }

    /// Moves chiral projectors to the right of a 4D gamma:
    /// `P+ gamma(mu) -> gamma(mu) P-` and
    /// `P- gamma(mu) -> gamma(mu) P+`.
    fn move_projector_right_of_gamma(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        if !has_four_dimensional_spin_endpoints(start, end) {
            return None;
        }

        for i in 0..factors.len().saturating_sub(1) {
            let opposite_projector = match factors[i] {
                DiracFactor::ProjectorPlus(_) => endpoint_factor(AGS.projm),
                DiracFactor::ProjectorMinus(_) => endpoint_factor(AGS.projp),
                _ => continue,
            };
            if factors[i + 1]
                .gamma_mink_index(FOUR_DIM_GAMMA5_ANTICOMMUTATION)
                .is_none()
            {
                continue;
            }

            let mut moved = Self::owned_factors(factors);
            moved[i] = factors[i + 1].as_view().into_atom();
            moved[i + 1] = opposite_projector;
            return Some(chain!(start, end; moved));
        }

        None
    }

    fn anticommute_adjacent_gamma_pair(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
        swap_at: usize,
    ) -> Option<Atom> {
        let (mu, nu) = Self::mink_index_pair(
            GAMMA_ANTICOMMUTATION,
            &factors[swap_at],
            &factors[swap_at + 1],
        )?;

        let mut metric_rest = Vec::with_capacity(factors.len() - 2);
        Self::extend_factors(&mut metric_rest, &factors[..swap_at]);
        Self::extend_factors(&mut metric_rest, &factors[swap_at + 2..]);

        // Run the local metric term through chain-aware Schoonschip before it can
        // swell the Clifford expansion.
        let metric_term =
            Atom::num(2) * Self::generated_metric_chain_term(start, end, mu, nu, metric_rest);

        let mut swapped = Self::owned_factors(factors);
        swapped.swap(swap_at, swap_at + 1);

        Some(metric_term - chain!(start, end; swapped))
    }

    fn generated_metric_chain_term(
        start: AtomView<'_>,
        end: AtomView<'_>,
        mu: AtomView<'_>,
        nu: AtomView<'_>,
        factors: Vec<Atom>,
    ) -> Atom {
        (function!(ETS.metric, mu, nu) * chain!(start, end; factors)).schoonschip_with_settings(
            &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
        )
    }

    fn shortest_repeated_gamma_pair(factors: &[DiracFactor<'_>]) -> Option<(usize, usize)> {
        let mut best = None;

        for i in 0..factors.len() {
            if factors[i].gamma_mink_index(GAMMA_ANTICOMMUTATION).is_none() {
                continue;
            }

            for (j, factor) in factors.iter().enumerate().skip(i + 1) {
                let Some((mu, nu)) =
                    Self::mink_index_pair(GAMMA_ANTICOMMUTATION, &factors[i], factor)
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

    fn shortest_repeated_four_dim_gamma_pair(
        factors: &[DiracFactor<'_>],
    ) -> Option<(usize, usize)> {
        let mut best = None;

        for i in 0..factors.len() {
            let Some(mu) = factors[i].gamma_mink_index(FOUR_DIM_CHISHOLM) else {
                continue;
            };
            if !is_minkowski_slot(mu) {
                continue;
            }

            for (j, factor) in factors.iter().enumerate().skip(i + 1) {
                let Some(nu) = factor.gamma_mink_index(FOUR_DIM_CHISHOLM) else {
                    continue;
                };

                if mu == nu && best.is_none_or(|(a, b)| j - i < b - a) {
                    best = Some((i, j));
                }
            }
        }

        best
    }

    /// Asks for the Minkowski indices of a gamma pair, checking that their
    /// dimensions are compatible with the rule.
    fn mink_index_pair<'a>(
        rule_dimension: DiracRuleDimension,
        left: &DiracFactor<'a>,
        right: &DiracFactor<'a>,
    ) -> Option<(AtomView<'a>, AtomView<'a>)> {
        let left_mink_index = (*left).gamma_mink_index(rule_dimension)?;
        let right_mink_index = (*right).gamma_mink_index(rule_dimension)?;
        if !rule_dimension.gamma_compatible(*left, *right) {
            return None;
        }

        Some((left_mink_index, right_mink_index))
    }

    /// Returns the Minkowski indices of a gamma-only sequence accepted by a
    /// rule.
    ///
    /// For arbitrary-dimensional rules, all known dimensions in the sequence
    /// must agree. Unknown dimensions are accepted so symbolic slash arguments
    /// without a visible representation can still flow through dimension-generic
    /// identities. Four-dimensional rules have already filtered every gamma
    /// through `allows_gamma_dimension`.
    fn gamma_mink_index_sequence_for<'a>(
        rule_dimension: DiracRuleDimension,
        factors: &[DiracFactor<'a>],
    ) -> Option<Vec<AtomView<'a>>> {
        let mut known_dimension = None;
        let mut mink_indices = Vec::with_capacity(factors.len());

        for factor in factors {
            let gamma_mink_index = (*factor).gamma_mink_index(rule_dimension)?;
            if !rule_dimension
                .merge_known_gamma_dimension(&mut known_dimension, (*factor).gamma_dimension())
            {
                return None;
            }
            mink_indices.push(gamma_mink_index);
        }

        Some(mink_indices)
    }
}

/// Infers the Minkowski dimension carried by a slot, that could be schoonschipped
///
/// Direct slots use the tensor-slot convention `mink(dim,index)`, so the first
/// argument is the dimension. Slash-like indices can be tensorial, e.g.
/// `P(1,mink(D))`; for those, the first visible Minkowski representation supplies
/// the dimension.
fn mink_slot_dimension(mink_index: AtomView<'_>) -> Option<AtomView<'_>> {
    if let Some(dimension) = minkowski_dimension(mink_index) {
        return Some(dimension);
    }

    let AtomView::Fun(f) = mink_index else {
        return None;
    };

    f.iter().find_map(minkowski_dimension)
}

/// Returns the first argument of a Minkowski representation function.
///
/// Spenso tensor slots are encoded as a function headed by the representation
/// symbol; the first argument is always the dimension and the optional second
/// argument is the index.
fn minkowski_dimension(atom: AtomView<'_>) -> Option<AtomView<'_>> {
    let AtomView::Fun(f) = atom else {
        return None;
    };

    if f.get_symbol() != *MINKOWSKI_SYMBOL || f.get_nargs() == 0 {
        return None;
    }

    f.iter().next()
}

/// Checks for a concrete Minkowski slot, not a stripped representation.
fn is_minkowski_slot(atom: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = atom else {
        return false;
    };

    f.get_symbol() == *MINKOWSKI_SYMBOL && f.get_nargs() == 2
}

fn is_four_dimension(dimension: AtomView<'_>) -> bool {
    matches!(i64::try_from(dimension), Ok(4))
}

/// Checks that open-chain spin endpoints both live in four-dimensional
/// bispinor space.
fn has_four_dimensional_spin_endpoints(start: AtomView<'_>, end: AtomView<'_>) -> bool {
    let Some(start_dimension) = bispinor_dimension(start) else {
        return false;
    };
    let Some(end_dimension) = bispinor_dimension(end) else {
        return false;
    };

    start_dimension == end_dimension && is_four_dimension(start_dimension)
}

fn has_four_dimensional_trace_rep(rep: AtomView<'_>) -> bool {
    bispinor_dimension(rep).is_some_and(is_four_dimension)
}

/// Returns the first argument of a bispinor representation function, following
/// the same `rep(dim,index)` slot convention as Minkowski slots.
fn bispinor_dimension(atom: AtomView<'_>) -> Option<AtomView<'_>> {
    let AtomView::Fun(f) = atom else {
        return None;
    };

    if f.get_symbol() != *BISPINOR_SYMBOL || f.get_nargs() == 0 {
        return None;
    }

    f.iter().next()
}

fn endpoint_factor(symbol: Symbol) -> Atom {
    FunctionBuilder::new(symbol)
        .add_arg(Atom::var(T.chain_in))
        .add_arg(Atom::var(T.chain_out))
        .finish()
}

fn gamma5_factor() -> Atom {
    endpoint_factor(AGS.gamma5)
}

fn gamma_factor(mink_index: impl IntoAtom) -> Atom {
    FunctionBuilder::new(AGS.gamma)
        .add_arg(Atom::var(T.chain_in))
        .add_arg(Atom::var(T.chain_out))
        .add_arg(mink_index.into_atom())
        .finish()
}

fn epsilon_dummy_minkowski_slot() -> Atom {
    Minkowski {}.to_symbolic([Atom::num(4), Atom::var(*EPSILON_DUMMY_SYMBOL)])
}

fn has_chain_endpoints(left: AtomView, right: AtomView) -> bool {
    is_chain_endpoint(left, T.chain_in) && is_chain_endpoint(right, T.chain_out)
        || is_chain_endpoint(left, T.chain_out) && is_chain_endpoint(right, T.chain_in)
}

fn is_chain_endpoint(arg: AtomView, expected: Symbol) -> bool {
    matches!(arg, AtomView::Var(symbol) if symbol.get_symbol() == expected)
}

impl DiracSimplifier<'_> {
    fn simplify_trace_node(self, f: FunView) -> Option<Atom> {
        let args = f.iter().collect::<Vec<_>>();
        let [rep, factors @ ..] = args.as_slice() else {
            return None;
        };

        if factors.is_empty() {
            return Self::simplify_trace_terminal(f.as_view());
        }

        let factors = factors
            .iter()
            .map(|factor| DiracFactor::parse(*factor))
            .collect::<Vec<_>>();

        if let Some(rewritten) = Self::simplify_special_trace_pair(*rep, &factors) {
            return Some(rewritten);
        }

        if let Some(rewritten) = Self::simplify_gamma5_trace_node(*rep, &factors) {
            return Some(rewritten);
        }

        let trace_mink_indices =
            Self::gamma_mink_index_sequence_for(TRACE_GAMMA_RECURSION, &factors)?;

        if factors.len() % 2 == 1 {
            return Some(Atom::Zero);
        }

        let first = trace_mink_indices[0];
        let mut sum = Atom::Zero;

        // Standard recursive even trace formula:
        // tr(g1...gn) = sum_i (-1)^i g(1,i) tr(g2...g_{i-1}g_{i+1}...gn).
        for i in 1..factors.len() {
            let mu_i = trace_mink_indices[i];
            let sign = if i % 2 == 1 { 1 } else { -1 };

            let mut rest = Vec::with_capacity(factors.len() - 2);
            Self::extend_factors(&mut rest, &factors[1..i]);
            Self::extend_factors(&mut rest, &factors[i + 1..]);

            let rest_trace = if rest.is_empty() {
                let terminal_trace = trace!(*rep; std::iter::empty::<Atom>());
                Self::simplify_trace_terminal(terminal_trace.as_view())?
            } else {
                trace!(*rep; rest)
            };

            let term = function!(ETS.metric, first, mu_i) * rest_trace;
            if sign == 1 {
                sum += term;
            } else {
                sum -= term;
            }
        }

        Some(sum)
    }

    fn simplify_special_trace_pair(rep: AtomView<'_>, factors: &[DiracFactor<'_>]) -> Option<Atom> {
        if !has_four_dimensional_trace_rep(rep) {
            return None;
        }

        for i in 0..factors.len().saturating_sub(1) {
            match (&factors[i], &factors[i + 1]) {
                (DiracFactor::Gamma5(_), DiracFactor::Gamma5(_))
                | (DiracFactor::Gamma0(_), DiracFactor::Gamma0(_)) => {
                    let mut rest = Vec::with_capacity(factors.len() - 2);
                    Self::extend_factors(&mut rest, &factors[..i]);
                    Self::extend_factors(&mut rest, &factors[i + 2..]);
                    return Some(Self::trace_or_terminal(rep, rest));
                }
                _ => {}
            }
        }

        None
    }

    fn simplify_gamma5_trace_node(rep: AtomView<'_>, factors: &[DiracFactor<'_>]) -> Option<Atom> {
        if !has_four_dimensional_trace_rep(rep) {
            return None;
        }

        let gamma5_positions = factors
            .iter()
            .enumerate()
            .filter_map(|(i, factor)| factor.is_gamma5().then_some(i))
            .collect::<Vec<_>>();

        match gamma5_positions.len() {
            0 => None,
            1 => Self::simplify_single_gamma5_trace(rep, factors, gamma5_positions[0]),
            _ => Self::reduce_gamma5_trace_pair(rep, factors, gamma5_positions[0]),
        }
    }

    fn simplify_single_gamma5_trace(
        rep: AtomView<'_>,
        factors: &[DiracFactor<'_>],
        gamma5_position: usize,
    ) -> Option<Atom> {
        let parsed_after_gamma5 = Self::cyclic_without_position(factors, gamma5_position);
        let mink_indices =
            Self::gamma_mink_index_sequence_for(TRACE_GAMMA5_RECURSION, &parsed_after_gamma5)?;

        if mink_indices.len() < 4 || mink_indices.len() % 2 == 1 {
            return Some(Atom::Zero);
        }

        if mink_indices.len() == 4 {
            return Some(
                Atom::num(4)
                    * epsilon4(
                        mink_indices[0],
                        mink_indices[1],
                        mink_indices[2],
                        mink_indices[3],
                    ),
            );
        }

        let first = mink_indices[0];
        let mut sum = Atom::Zero;

        for i in 1..parsed_after_gamma5.len() {
            let mut rest = Vec::with_capacity(parsed_after_gamma5.len() - 1);
            rest.push(gamma5_factor());
            Self::extend_factors(&mut rest, &parsed_after_gamma5[1..i]);
            Self::extend_factors(&mut rest, &parsed_after_gamma5[i + 1..]);

            let term =
                function!(ETS.metric, first, mink_indices[i]) * Self::trace_or_terminal(rep, rest);

            if i % 2 == 1 {
                sum += term;
            } else {
                sum -= term;
            }
        }

        Some(sum)
    }

    fn reduce_gamma5_trace_pair(
        rep: AtomView<'_>,
        factors: &[DiracFactor<'_>],
        first_gamma5_position: usize,
    ) -> Option<Atom> {
        let factors = Self::cyclic_from_position(factors, first_gamma5_position);
        let second_gamma5_position = factors
            .iter()
            .enumerate()
            .skip(1)
            .find_map(|(i, factor)| factor.is_gamma5().then_some(i))?;

        let crossing_count = factors[1..second_gamma5_position]
            .iter()
            .filter(|factor| Self::factor_anticommutes_with_gamma5(factor))
            .count();
        if factors[1..second_gamma5_position]
            .iter()
            .any(|factor| !Self::factor_anticommutes_with_gamma5(factor))
        {
            return None;
        }

        let mut rest = Vec::with_capacity(factors.len() - 2);
        Self::extend_factors(&mut rest, &factors[1..second_gamma5_position]);
        Self::extend_factors(&mut rest, &factors[second_gamma5_position + 1..]);

        let reduced = Self::trace_or_terminal(rep, rest);
        Some(if crossing_count % 2 == 0 {
            reduced
        } else {
            Atom::num(-1) * reduced
        })
    }

    fn factor_anticommutes_with_gamma5(factor: &DiracFactor<'_>) -> bool {
        factor.anticommutes_with_gamma5()
    }

    fn cyclic_without_position<T: Clone>(items: &[T], position: usize) -> Vec<T> {
        let mut result = Vec::with_capacity(items.len().saturating_sub(1));
        result.extend_from_slice(&items[position + 1..]);
        result.extend_from_slice(&items[..position]);
        result
    }

    fn cyclic_from_position<T: Clone>(items: &[T], position: usize) -> Vec<T> {
        let mut result = Vec::with_capacity(items.len());
        result.extend_from_slice(&items[position..]);
        result.extend_from_slice(&items[..position]);
        result
    }

    fn trace_or_terminal(rep: AtomView<'_>, factors: Vec<Atom>) -> Atom {
        if factors.is_empty() {
            let terminal_trace = trace!(rep; std::iter::empty::<Atom>());
            Self::simplify_trace_terminal(terminal_trace.as_view()).unwrap_or(terminal_trace)
        } else {
            trace!(rep; factors)
        }
    }

    fn simplify_trace_terminal(trace: AtomView) -> Option<Atom> {
        let trace = trace.to_owned();
        let simplified = trace.replace_multiple_repeat(TRACE_TERMINALS.as_ref());
        (simplified != trace).then_some(simplified)
    }
}
