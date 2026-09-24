use std::sync::LazyLock;

use spenso::{
    chain, g,
    network::tags::SPENSO_TAG as T,
    rep_,
    shadowing::{self, IntoAtom},
    structure::{
        abstract_index::AbstractIndex,
        representation::{LibraryRep, Minkowski, RepName},
        slot::{DummyAind, ParseableAind, SlotMatch, SlotMatcher},
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, representation::FunView},
    coefficient::CoefficientView,
    id::{Context, Replacement},
    utils::Settable,
};
use symbolica_utils::PatternReplacement;

use crate::{
    W_,
    dirac::GammaSimplifier,
    epsilon::{EpsilonSimplifier, epsilon4},
    representations::Bispinor,
    shorthands::schoonschip::{Schoonschip, SchoonschipSettings, SimplificationCandidates},
};

use super::{AGS, id_atom};

mod trace_kernel;

// gamma(a) gamma(b) gamma(c): the metric part of the 4D decomposition.
// The remaining term is -epsilon(a,b,c,d) gamma(d) gamma5.
const THREE_GAMMA_METRIC_TERMS: [(i32, [usize; 2], usize); 3] =
    [(1, [0, 1], 2), (-1, [0, 2], 1), (1, [1, 2], 0)];

#[derive(Clone, Copy)]
enum DiracWord<'a> {
    Chain(AtomView<'a>, AtomView<'a>),
    Trace(AtomView<'a>),
}

impl DiracWord<'_> {
    fn build(self, factors: Vec<Atom>) -> Atom {
        match self {
            Self::Chain(start, end) => chain!(start, end; factors),
            Self::Trace(rep) => DiracSimplifier::trace_or_terminal(rep, factors),
        }
    }

    fn is_four_dimensional(self) -> bool {
        match self {
            Self::Chain(start, end) => has_four_dimensional_spin_endpoints(start, end),
            Self::Trace(rep) => has_four_dimensional_trace_rep(rep),
        }
    }
}

static MINKOWSKI_SYMBOL: LazyLock<Symbol> =
    LazyLock::new(|| LibraryRep::from(Minkowski {}).symbol());

static BISPINOR_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| LibraryRep::from(Bispinor {}).symbol());

static EPSILON_DUMMY_SYMBOL: LazyLock<Symbol> = LazyLock::new(|| symbolica::symbol!("sigma"));

static TRACE_TERMINALS: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    [Replacement::new(
        // Empty spin trace: Tr_rep(1) -> dim(rep).
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
    /// Fully expand each evaluated trace body, preserving surrounding factors.
    /// Has no effect when `evaluate_traces` is false.
    pub expand_traces: bool,
    /// Whether three 4D gammas may be expanded into the gamma5-epsilon basis.
    pub expand_three_gamma_epsilon: bool,
}

impl Default for GammaSimplifySettings {
    fn default() -> Self {
        Self {
            chain_ordering: GammaChainOrdering::RepeatedPairs,
            evaluate_traces: true,
            expand_traces: false,
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

    /// Request a fully expanded result for each evaluated trace body.
    pub fn with_expanded_traces(mut self) -> Self {
        self.expand_traces = true;
        self
    }

    /// Enables the four-dimensional identity that rewrites three gammas into
    /// metric terms plus a gamma5-epsilon term.
    pub fn with_gamma5_epsilon_expansion(mut self) -> Self {
        self.expand_three_gamma_epsilon = true;
        self
    }

    fn rewrite_expression(&self, expr: Atom) -> Atom {
        // Empty chains and traces still have identities, even without gamma
        // factors. Only the absence of both eligible heads makes this a no-op.
        let mut candidate = false;
        expr.visitor(&mut |node| {
            if candidate {
                return false;
            }
            if let AtomView::Fun(function) = node {
                let head = function.get_symbol();
                candidate = head == T.chain || (self.evaluate_traces && head == T.trace);
            }
            !candidate
        });
        if !candidate {
            return expr;
        }
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
            // Complete the original trace body, including coefficient factors,
            // before expanding it. Nested callback-created traces retain the
            // same expansion setting; surrounding spectators stay outside.
            **out = if self.expand_traces && rewritten.as_view() != arg {
                simplifier
                    .simplify_complete::<false>(rewritten.as_view())
                    .expand()
            } else {
                rewritten
            };
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
    ChargeConjugation(AtomView<'a>),
    ProjectorPlus(AtomView<'a>),
    ProjectorMinus(AtomView<'a>),
    Other(AtomView<'a>),
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
struct DiracFactorKinds {
    has_gamma5: bool,
    has_gamma0: bool,
    has_charge_conjugation: bool,
    has_projector: bool,
}

impl DiracFactorKinds {
    fn observe(&mut self, factor: DiracFactor<'_>) {
        match factor {
            DiracFactor::Gamma5(_) => self.has_gamma5 = true,
            DiracFactor::Gamma0(_) => self.has_gamma0 = true,
            DiracFactor::ChargeConjugation(_) => self.has_charge_conjugation = true,
            DiracFactor::ProjectorPlus(_) | DiracFactor::ProjectorMinus(_) => {
                self.has_projector = true;
            }
            DiracFactor::Gamma { .. } | DiracFactor::Other(_) => {}
        }
    }

    fn has_non_pure_gamma(self) -> bool {
        self.has_gamma5 || self.has_gamma0 || self.has_charge_conjugation || self.has_projector
    }

    fn has_special_trace_pair(self) -> bool {
        self.has_gamma5 || self.has_gamma0 || self.has_charge_conjugation
    }

    fn has_conjugation_candidate(self) -> bool {
        self.has_charge_conjugation || (self.has_gamma0 && (self.has_gamma5 || self.has_projector))
    }
}

impl<'a> DiracFactor<'a> {
    /// Classifies a factor inside a `chain(...)` or `trace(...)` without
    /// materializing owned atoms.
    ///
    /// Ordinary gammas are recognized as `gamma(in,out,mu)` and keep `mu`
    /// by view. Reversed matrix endpoints remain opaque to ordinary Clifford
    /// and special-matrix rules, which cannot safely act on mixed
    /// ordinary/transposed words. An explicit four-dimensional C sandwich may
    /// flip their orientation using the charge-conjugation identities.
    /// Their dimension is inferred
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
            if has_forward_chain_endpoints(left, right) {
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
            // Antisymmetric normalization puts every C in the same endpoint
            // order, extracting any transposition sign into the chain scalar.
            if f.get_symbol() == AGS.charge_conjugation
                && (has_forward_chain_endpoints(left, right)
                    || has_forward_chain_endpoints(right, left))
            {
                return Self::ChargeConjugation(factor);
            }
            if has_forward_chain_endpoints(left, right) {
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
            | Self::ChargeConjugation(factor)
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
        self.simplify_complete::<true>(expr)
    }

    // Local trace-body completion must not grant terminal admission to a
    // subtrace whose explicit slots can still meet its surrounding factors.
    fn simplify_complete<const TERMINAL_CONTEXT: bool>(self, expr: AtomView<'_>) -> Atom {
        if TERMINAL_CONTEXT
            && self.settings.evaluate_traces
            && let Some(result) = if self.settings.expand_traces {
                // Exact callback-free words already satisfy terminal closure,
                // including mixed free slots and compact vector components.
                self.evaluate_terminal_trace::<true>(expr)
                    .or_else(|| self.evaluate_terminal_trace::<false>(expr))
            } else {
                self.evaluate_terminal_trace::<false>(expr)
            }
        {
            return result;
        }
        let mut expr = expr.to_owned();
        // A trace can emit metrics connected to surviving words. Contract those
        // before each rewrite, so later iterations expose their repeated indices.
        let metric_settings = SchoonschipSettings::default().with_chain_like_functions();
        let bispinor: LibraryRep = Bispinor {}.into();
        let heads = [
            bispinor.symbol(),
            T.chain,
            T.bracket,
            T.trace,
            *crate::epsilon::EPSILON_SYMBOL,
        ];
        let [bispinor, chain, bracket, trace, epsilon_head] = heads.map(|head| head.get_id());
        let mut observed = Some(SimplificationCandidates::scan(expr.as_view(), heads));
        if TERMINAL_CONTEXT
            && self.settings.evaluate_traces
            && observed.as_ref().is_some_and(|candidates| candidates.symbols[3])
            && let AtomView::Mul(product) = expr.as_view()
            && let Some(trace) = product.iter().find(
                |factor| matches!(factor, AtomView::Fun(function) if function.get_symbol_id() == trace),
            )
            && product.iter().all(|factor| {
                if factor == trace {
                    return true;
                }
                // Exact function-free scalars have no tensor slots or callbacks.
                // Approximate coefficients retain the original evaluation order.
                let mut scalar = true;
                factor.visitor(&mut |node| {
                    scalar &= match node {
                        AtomView::Fun(_) => false,
                        AtomView::Num(_) => Self::exact_scalar_leaf(node),
                        _ => true,
                    };
                    scalar
                });
                scalar
            })
            && let Some(result) = self.evaluate_terminal_trace::<true>(trace)
        {
            // The terminal kernel emits each free explicit slot once per
            // term. Canonical vector components have no callbacks, and the exact
            // function-free spectators cannot introduce tensor cleanup work.
            return Atom::mul_many(product.iter().map(|factor| {
                if factor == trace { result.as_view() } else { factor }
            }));
        }

        loop {
            // Metric contraction can close a chain or connect separate chains.
            // Include their collection in the fixed point of the complete pass.
            // Close metric-linked chains before rewriting, including inert traces.
            let candidates = observed
                .take()
                .unwrap_or_else(|| SimplificationCandidates::scan(expr.as_view(), heads));
            let normalized = if candidates.normalized() {
                expr.clone()
            } else {
                expr.schoonschip_with_settings(&metric_settings)
            };
            // Share the absence check across the remaining passes. Large trace
            // results need no chain collection or Dirac rewrite, and ordinary
            // generic-dimensional traces contain no epsilon either.
            let (mut collect, mut rewrite, mut epsilon) = (false, false, false);
            if candidates.complete && normalized == expr {
                let [has_bispinor, has_chain, has_bracket, has_trace, has_epsilon] =
                    candidates.symbols;
                collect = has_bispinor || has_chain || has_bracket;
                rewrite = has_chain || (self.settings.evaluate_traces && has_trace);
                epsilon = has_epsilon;
            } else {
                normalized.visitor(&mut |node| {
                    let head = match node {
                        AtomView::Fun(function) => function.get_symbol_id(),
                        AtomView::Var(variable) => variable.get_symbol_id(),
                        _ => return true,
                    };
                    collect |= head == bispinor || head == chain || head == bracket;
                    rewrite |= head == chain || (self.settings.evaluate_traces && head == trace);
                    epsilon |= head == epsilon_head;
                    !(collect && rewrite && epsilon)
                });
            }
            let mut next = if collect {
                normalized.collect_gamma_chains()
            } else {
                normalized.clone()
            };
            // A preceding pass can introduce work absent from the initial scan:
            // collection closes traces, and a trace rewrite can emit epsilons.
            if rewrite || next != normalized {
                next = self.settings.rewrite_expression(next);
            }
            // Inspect the whole rebuilt expression: unchanged outer factors
            // may now contract with the trace output or callback result.
            let rewritten = next != normalized;
            if rewritten {
                let complete = SimplificationCandidates::scan(next.as_view(), heads);
                if complete.finished() {
                    return next;
                }
                // Carry this scan only across exactly unchanged cleanup.
                observed = Some(complete);
            }
            if epsilon || rewritten {
                let updated = next.simplify_epsilon();
                if updated != next {
                    observed = None;
                }
                next = updated;
            }
            // Schoonschip already normalized dots. Repeat that walk only when a
            // later operation changed its result.
            if next != normalized {
                let updated = next.normalize_dots();
                if updated != next {
                    observed = None;
                }
                next = updated;
            }

            if next == expr {
                return next;
            }

            expr = next;
        }
    }

    // Only exact rational coefficients may be regrouped by the scalar-context
    // trace recurrence. Other numeric domains keep the original schedule.
    fn exact_scalar_leaf(value: AtomView<'_>) -> bool {
        match value {
            AtomView::Var(_) => true,
            AtomView::Num(number) => matches!(
                number.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            ),
            _ => false,
        }
    }

    /// Evaluate standalone traces after reducing summed indices within the
    /// word. The surviving explicit indices occur once per monomial and compact
    /// arguments produce scalar dots, so no outer tensor contraction remains.
    /// A scalar-product context additionally requires callback-free ordinary
    /// words with exact leaf metadata. Free slots occur once per term; wildcard
    /// metadata may retain inert closed diagonals. Exact scalar spectators cannot
    /// introduce further tensor cleanup.
    fn evaluate_terminal_trace<const SCALAR_CONTEXT: bool>(
        self,
        expr: AtomView<'_>,
    ) -> Option<Atom> {
        let output = if self.settings.expand_traces {
            trace_kernel::TraceOutput::Expanded
        } else {
            trace_kernel::TraceOutput::Factored
        };
        let AtomView::Fun(f) = expr else {
            return None;
        };
        let (rep, factors) = shadowing::trace_parts(f)?;
        if !SCALAR_CONTEXT && self.settings.expand_traces {
            // Scalar parameters and representation metadata may contain nested
            // traces. The strict scalar-context admission already rejects them.
            for argument in f.iter() {
                let mut nested_trace = false;
                argument.visitor(&mut |node| {
                    nested_trace |=
                        matches!(node, AtomView::Fun(inner) if inner.get_symbol() == T.trace);
                    !nested_trace
                });
                if nested_trace {
                    return None;
                }
            }
        }
        if SCALAR_CONTEXT
            && !matches!(rep, AtomView::Fun(spin)
                if spin.get_symbol() == *BISPINOR_SYMBOL
                    && spin.get_nargs() == 1
                    && spin.iter().next().is_some_and(Self::exact_scalar_leaf))
        {
            return None;
        }
        let mut factors = factors
            .iter()
            .copied()
            .map(DiracFactor::parse)
            .collect::<Vec<_>>();
        let gamma5_position = factors.iter().position(|factor| factor.is_gamma5());
        if SCALAR_CONTEXT && gamma5_position.is_some() {
            return None;
        }
        // Trace cyclicity moves one gamma5 to the front without a sign. Any
        // further gamma5 or non-gamma factor is rejected by the sequence parser.
        if let Some(position) = gamma5_position {
            factors = Self::cyclic_without_position(&factors, position);
        }
        let axial = gamma5_position.is_some();
        let rule = if axial {
            if !has_four_dimensional_trace_rep(rep) {
                return None;
            }
            FOUR_DIM_CHISHOLM
        } else {
            TRACE_GAMMA_RECURSION
        };
        let indices = Self::gamma_mink_index_sequence_for(rule, &factors)?;
        let mut slots = SlotMatcher::default();
        let mut repeated = 0;
        let mut compact_count = 0;
        for (position, &index) in indices.iter().enumerate() {
            if is_minkowski_slot(index) {
                let SlotMatch::Explicit(slot) = slots.classify(index) else {
                    return None;
                };
                if slot.is_concrete_index()
                    || (SCALAR_CONTEXT
                        && (!Self::exact_scalar_leaf(slot.dimension())
                            || !matches!(slot.index(), AtomView::Var(_))))
                {
                    return None;
                }
                // More than two occurrences do not define an Einstein sum.
                // Preserve the complete pass's existing ordered substitutions.
                let occurrences = indices[..position].iter().filter(|&&i| i == index).count();
                if occurrences > usize::from(!axial) {
                    return None;
                }
                repeated += occurrences;
            } else {
                compact_count += 1;
                let AtomView::Fun(vector) = index else {
                    return None;
                };
                if axial
                    || (SCALAR_CONTEXT
                        && (vector.get_symbol().get_normalization_function().is_some()
                            || !vector
                                .iter()
                                .take(vector.get_nargs().saturating_sub(1))
                                .all(Self::exact_scalar_leaf)))
                    || !vector.get_symbol().has_tag(&T.rank1)
                    || !slots
                        .vector_argument(vector)
                        .and_then(|argument| slots.compact_representation(argument))
                        .is_some_and(|rep| {
                            rep.head() == *MINKOWSKI_SYMBOL
                                && rep.is_base()
                                && (!SCALAR_CONTEXT || Self::exact_scalar_leaf(rep.dimension()))
                        })
                {
                    return None;
                }
            }
        }
        if axial && indices.len() > 14 {
            return None;
        }
        if indices.len() % 2 == 1 || (axial && indices.len() < 4) {
            return Some(Atom::Zero);
        }
        let compact = !axial && indices.iter().all(|&index| !is_minkowski_slot(index));
        let result = if !compact
            && repeated == 0
            && (axial
                || (has_four_dimensional_trace_rep(rep)
                    && factors
                        .iter()
                        .all(|factor| factor.gamma_dimension().is_some_and(is_four_dimension))))
        {
            trace_kernel::evaluate(&indices, axial, output)?
        } else {
            let terminal = trace!(rep; std::iter::empty::<Atom>());
            let trace_unit = Self::simplify_trace_terminal(terminal.as_view())?;
            trace_kernel::evaluate_generic(
                &indices,
                trace_unit.as_view(),
                // Distinct explicit arguments need no factored reductions.
                // Expanded emission uses the canonical admission contract.
                self.settings.expand_traces
                    || !(SCALAR_CONTEXT
                        && repeated == 0
                        && indices.iter().all(|&index| is_minkowski_slot(index))),
                output,
            )
        };
        let free_explicit = indices.len() - compact_count - 2 * repeated;
        Some(
            if !SCALAR_CONTEXT
                && self.settings.expand_traces
                && compact_count > 0
                && free_explicit > 0
            {
                // Broad standalone admission permits component callbacks, which
                // can create sums or nested traces. Strict context excludes them.
                self.simplify_complete::<false>(result.as_view()).expand()
            } else {
                result
            },
        )
    }

    fn simplify_chain_node(self, f: FunView) -> Option<Atom> {
        let args = f.iter().collect::<Vec<_>>();
        let [start, end, factors @ ..] = args.as_slice() else {
            return None;
        };

        let mut factor_kinds = DiracFactorKinds::default();
        let factors = factors
            .iter()
            .map(|factor| {
                let factor = DiracFactor::parse(*factor);
                factor_kinds.observe(factor);
                factor
            })
            .collect::<Vec<_>>();

        if factors.is_empty() {
            return Some(id_atom(*start, *end));
        }

        Self::contract_adjacent_gamma_pair(DiracWord::Chain(*start, *end), &factors)
            .or_else(|| {
                factor_kinds
                    .has_non_pure_gamma()
                    .then(|| Self::contract_adjacent_special_dirac_pair(*start, *end, &factors))
                    .flatten()
            })
            .or_else(|| {
                (factor_kinds.has_conjugation_candidate()
                    && has_four_dimensional_spin_endpoints(*start, *end))
                .then(|| Self::conjugate_special_dirac_factor(&factors))
                .flatten()
                .map(|(sign, factors)| Atom::num(sign) * chain!(*start, *end; factors))
            })
            .or_else(|| {
                Self::four_dim_chisholm_contraction(DiracWord::Chain(*start, *end), &factors)
            })
            .or_else(|| {
                self.settings
                    .expand_three_gamma_epsilon
                    .then(|| {
                        Self::four_dim_three_gamma_epsilon_expansion(
                            DiracWord::Chain(*start, *end),
                            &factors,
                        )
                    })
                    .flatten()
            })
            .or_else(|| {
                factor_kinds
                    .has_gamma5
                    .then(|| Self::move_gamma5_right_of_gamma(*start, *end, &factors))
                    .flatten()
            })
            .or_else(|| {
                factor_kinds
                    .has_projector
                    .then(|| Self::move_projector_right_of_gamma(*start, *end, &factors))
                    .flatten()
            })
            .or_else(|| match self.settings.chain_ordering {
                GammaChainOrdering::RepeatedPairs => {
                    Self::bubble_repeated_gamma_towards_contraction(
                        DiracWord::Chain(*start, *end),
                        &factors,
                    )
                }
                GammaChainOrdering::Canonical => {
                    Self::canonicalize_gamma_chain_order(DiracWord::Chain(*start, *end), &factors)
                }
            })
    }
}

impl DiracSimplifier<'_> {
    /// Contracts adjacent equal-index gammas:
    /// `...[gamma(mu), gamma(mu)]... -> g(mu, mu) * ...[...]...`.
    fn contract_adjacent_gamma_pair(
        word: DiracWord<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        for (i, pair) in factors.windows(2).enumerate() {
            let [left, right] = pair else {
                unreachable!("windows(2) always yields pairs")
            };
            let Some((mu, nu)) = Self::mink_index_pair(ADJACENT_GAMMA_CONTRACTION, left, right)
            else {
                continue;
            };

            if mu != nu {
                continue;
            }

            let mut rest = Vec::with_capacity(factors.len() - 2);
            Self::extend_factors(&mut rest, &factors[..i]);
            Self::extend_factors(&mut rest, &factors[i + 2..]);

            return Some(g!(mu, nu) * word.build(rest));
        }

        None
    }

    /// Reduces adjacent special 4D factors:
    /// `gamma5 gamma5 -> 1`, `gamma0 gamma0 -> 1`, `P+ P+ -> P+`,
    /// `P- P- -> P-`, `P+ P- -> 0`, and `C C -> -1`.
    fn contract_adjacent_special_dirac_pair(
        start: AtomView<'_>,
        end: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        if !has_four_dimensional_spin_endpoints(start, end) {
            return None;
        }

        for (i, pair) in factors.windows(2).enumerate() {
            let [left, right] = pair else {
                unreachable!("windows(2) always yields pairs")
            };
            let sign = if matches!(
                (left, right),
                (
                    DiracFactor::ChargeConjugation(_),
                    DiracFactor::ChargeConjugation(_)
                )
            ) {
                -1
            } else {
                1
            };
            let replacement = match (left, right) {
                (DiracFactor::ChargeConjugation(_), DiracFactor::ChargeConjugation(_))
                | (DiracFactor::Gamma5(_), DiracFactor::Gamma5(_))
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

            return Some(
                Atom::num(sign)
                    * chain!(start, end; Self::chain_factors(
                        factors,
                        i,
                        i + 1,
                        replacement,
                    )),
            );
        }

        None
    }

    /// Conjugates special 4D factors by `gamma0` or charge conjugation:
    /// `gamma0 gamma5 gamma0 -> -gamma5`,
    /// `gamma0 P+ gamma0 -> P-`, and `gamma0 P- gamma0 -> P+`.
    /// With C^-1 = -C, `C gamma(mu) C = gamma(mu)^T` and
    /// `C gamma(mu)^T C = gamma(mu)`. For a word, insert C^-1 C between
    /// factors: `C A1...An C = -product(sigma_i) A1^T...An^T`, with
    /// sigma = -1 for gamma/gamma0 and +1 for gamma5/projectors. Preserve
    /// factor order: this is a product of transposes, not the transposed word.
    /// Symbolic-D gamma sandwiches stay opaque; their charge-conjugation
    /// convention has not been specified.
    fn conjugate_special_dirac_factor(factors: &[DiracFactor<'_>]) -> Option<(i64, Vec<Atom>)> {
        for (i, left) in factors.iter().enumerate() {
            if matches!(left, DiracFactor::Gamma0(_)) {
                let [middle, DiracFactor::Gamma0(_), ..] = &factors[i + 1..] else {
                    continue;
                };
                let (sign, replacement) = match *middle {
                    DiracFactor::Gamma5(_) => (-1, gamma5_factor()),
                    DiracFactor::ProjectorPlus(_) => (1, endpoint_factor(AGS.projm)),
                    DiracFactor::ProjectorMinus(_) => (1, endpoint_factor(AGS.projp)),
                    _ => continue,
                };
                return Some((sign, Self::chain_factors(factors, i, i + 2, [replacement])));
            }
            if !matches!(left, DiracFactor::ChargeConjugation(_)) {
                continue;
            }
            let Some(end) = factors
                .iter()
                .enumerate()
                .skip(i + 1)
                .find_map(|(index, factor)| {
                    matches!(factor, DiracFactor::ChargeConjugation(_)).then_some(index)
                })
            else {
                continue;
            };
            let mut sign = -1;
            let mut transposes = Vec::with_capacity(end - i - 1);
            for middle in &factors[i + 1..end] {
                let AtomView::Fun(matrix) = middle.as_view() else {
                    break;
                };
                let args = matrix.iter().collect::<Vec<_>>();
                let [row, column, rest @ ..] = args.as_slice() else {
                    break;
                };
                if !has_forward_chain_endpoints(*row, *column)
                    && !has_forward_chain_endpoints(*column, *row)
                {
                    break;
                }
                let sigma = match (matrix.get_symbol(), rest) {
                    (symbol, [mu])
                        if symbol == AGS.gamma
                            && mink_slot_dimension(*mu).is_some_and(is_four_dimension) =>
                    {
                        -1
                    }
                    (symbol, []) if symbol == AGS.gamma0 => -1,
                    (symbol, [])
                        if symbol == AGS.gamma5 || symbol == AGS.projp || symbol == AGS.projm =>
                    {
                        1
                    }
                    _ => break,
                };
                sign *= sigma;
                transposes.push(
                    FunctionBuilder::new(matrix.get_symbol())
                        .add_arg(*column)
                        .add_arg(*row)
                        .add_args(rest.iter().copied())
                        .finish(),
                );
            }
            // A single unsupported factor invalidates the whole sandwich.
            if transposes.len() == end - i - 1 {
                return Some((sign, Self::chain_factors(factors, i, end, transposes)));
            }
        }
        None
    }

    /// Moves repeated gammas toward each other with
    /// `gamma(mu) gamma(nu) -> 2 g(mu,nu) - gamma(nu) gamma(mu)`.
    ///
    /// Repeated-pair mode only pays the anticommutation cost when it exposes a
    /// contraction in the next fixed-point step.
    fn bubble_repeated_gamma_towards_contraction(
        word: DiracWord<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        let (_i, j) = Self::shortest_repeated_gamma_pair(factors)?;
        if j == 0 {
            return None;
        }

        Self::anticommute_adjacent_gamma_pair(word, factors, j - 1)
    }

    /// Applies the 4D Chisholm contractions around repeated endpoint gammas:
    /// odd interiors reverse with factor -2; even interiors give two words.
    /// The two-gamma interior has the shorter metric terminal.
    fn four_dim_chisholm_contraction(
        word: DiracWord<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        let (left, right) = Self::repeated_four_dim_gamma_pair(word, factors)?;
        let parsed_interior = &factors[left + 1..right];
        let interior_mink_indices =
            Self::gamma_mink_index_sequence_for(FOUR_DIM_CHISHOLM, parsed_interior)?;

        match parsed_interior.len() {
            0 => None,
            2 => Some(
                Atom::num(4)
                    * g!(interior_mink_indices[0], interior_mink_indices[1])
                    * word.build(Self::chain_factors(
                        factors,
                        left,
                        right,
                        std::iter::empty::<Atom>(),
                    )),
            ),
            n if n % 2 == 1 => Some(
                Atom::num(-2)
                    * word.build(Self::chain_factors(
                        factors,
                        left,
                        right,
                        parsed_interior
                            .iter()
                            .rev()
                            .copied()
                            .map(DiracFactor::as_view),
                    )),
            ),
            n => {
                let last = parsed_interior[n - 1].as_view();
                let rest = &parsed_interior[..n - 1];
                let first = rest
                    .iter()
                    .rev()
                    .copied()
                    .map(DiracFactor::as_view)
                    .chain([last]);
                let second = [last]
                    .into_iter()
                    .chain(rest.iter().copied().map(DiracFactor::as_view));
                Some(
                    Atom::num(2) * word.build(Self::chain_factors(factors, left, right, first))
                        + Atom::num(2)
                            * word.build(Self::chain_factors(factors, left, right, second)),
                )
            }
        }
    }

    /// Expands three 4D gammas into metric terms plus an epsilon-gamma5 term:
    /// `gamma(mu) gamma(nu) gamma(rho) -> g(mu,nu) gamma(rho)
    /// - g(mu,rho) gamma(nu) + g(nu,rho) gamma(mu)
    /// - epsilon(mu,nu,rho,sigma) gamma(sigma) gamma5`.
    fn four_dim_three_gamma_epsilon_expansion(
        word: DiracWord<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        if !word.is_four_dimensional() {
            return None;
        }

        for (i, triple) in factors.windows(3).enumerate() {
            let Some(mink_indices) =
                Self::gamma_mink_index_sequence_for(FOUR_DIM_THREE_GAMMA_EPSILON, triple)
            else {
                continue;
            };
            if matches!(word, DiracWord::Chain(..))
                && !mink_indices.iter().copied().all(is_minkowski_slot)
            {
                continue;
            }

            let [mu, nu, rho] = mink_indices.as_slice() else {
                unreachable!("the window always contains three gamma factors")
            };
            let sigma = match word {
                DiracWord::Chain(..) => epsilon_dummy_minkowski_slot(),
                // A longer axial trace may need several reductions before
                // its epsilon tensors contract; each needs its own dummy.
                DiracWord::Trace(_) => Minkowski {}
                    .new_rep(4)
                    .pattern(AbstractIndex::new_dummy().to_atom()),
            };

            // The chain normal form keeps gamma5 to the right of ordinary gammas.
            let epsilon_term = Atom::num(-1)
                * word.build(Self::chain_factors(
                    factors,
                    i,
                    i + 2,
                    [gamma_factor(sigma.clone()), gamma5_factor()],
                ))
                * epsilon4(*mu, *nu, *rho, &sigma);

            let metric_terms = THREE_GAMMA_METRIC_TERMS.map(|(coefficient, [a, b], remaining)| {
                Atom::num(coefficient)
                    * Self::generated_metric_word_term(
                        word,
                        mink_indices[a],
                        mink_indices[b],
                        Self::chain_factors(factors, i, i + 2, [triple[remaining].as_view()]),
                    )
            });
            return Some(epsilon_term + Atom::add_many(metric_terms));
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
        word: DiracWord<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
        for (i, pair) in factors.windows(2).enumerate() {
            let [left, right] = pair else {
                unreachable!("windows(2) always yields pairs")
            };
            let Some((mu, nu)) = Self::mink_index_pair(GAMMA_ANTICOMMUTATION, left, right) else {
                continue;
            };

            if mu > nu {
                return Self::anticommute_adjacent_gamma_pair(word, factors, i);
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

        for (i, pair) in factors.windows(2).enumerate() {
            let [left, right] = pair else {
                unreachable!("windows(2) always yields pairs")
            };
            let DiracFactor::Gamma5(_) = *left else {
                continue;
            };
            if (*right)
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

        for (i, pair) in factors.windows(2).enumerate() {
            let [left, right] = pair else {
                unreachable!("windows(2) always yields pairs")
            };
            let opposite_projector = match *left {
                DiracFactor::ProjectorPlus(_) => endpoint_factor(AGS.projm),
                DiracFactor::ProjectorMinus(_) => endpoint_factor(AGS.projp),
                _ => continue,
            };
            if (*right)
                .gamma_mink_index(FOUR_DIM_GAMMA5_ANTICOMMUTATION)
                .is_none()
            {
                continue;
            }

            let mut moved = Self::owned_factors(factors);
            moved[i] = (*right).as_view().into_atom();
            moved[i + 1] = opposite_projector;
            return Some(chain!(start, end; moved));
        }

        None
    }

    fn anticommute_adjacent_gamma_pair(
        word: DiracWord<'_>,
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
            Atom::num(2) * Self::generated_metric_word_term(word, mu, nu, metric_rest);

        let mut swapped = Self::owned_factors(factors);
        swapped.swap(swap_at, swap_at + 1);

        Some(metric_term - word.build(swapped))
    }

    fn generated_metric_word_term(
        word: DiracWord<'_>,
        mu: AtomView<'_>,
        nu: AtomView<'_>,
        factors: Vec<Atom>,
    ) -> Atom {
        (g!(mu, nu) * word.build(factors)).schoonschip_with_settings(
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

    fn repeated_four_dim_gamma_pair(
        word: DiracWord<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<(usize, usize)> {
        let mut best = None;
        let priority = |(left, right)| match word {
            DiracWord::Trace(_) => Self::trace_pair_priority(factors, left, right - left),
            DiracWord::Chain(..) => (0, right - left),
        };

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

                if mu == nu {
                    let score = priority((i, j));
                    if best.is_none_or(|(_, best_score)| score < best_score) {
                        best = Some(((i, j), score));
                    }
                }
            }
        }

        best.map(|(pair, _)| pair)
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

fn has_forward_chain_endpoints(left: AtomView, right: AtomView) -> bool {
    is_chain_endpoint(left, T.chain_in) && is_chain_endpoint(right, T.chain_out)
}

fn is_chain_endpoint(arg: AtomView, expected: Symbol) -> bool {
    matches!(arg, AtomView::Var(symbol) if symbol.get_symbol() == expected)
}

impl DiracSimplifier<'_> {
    /// Evaluates a Dirac trace.
    ///
    /// The pass first handles four-dimensional special factors (`gamma0`,
    /// `gamma5`, chiral projectors, charge conjugation), then reduces repeated ordinary gammas using the shared chain identities.
    /// Generic dimensions build the factored pairing polynomial directly.
    /// Four-dimensional words of length at most fourteen use shorter generated
    /// kernels; longer 4D words recurse until those kernels apply.
    fn simplify_trace_node(self, f: FunView) -> Option<Atom> {
        let (rep, factors) = shadowing::trace_parts(f)?;

        if factors.is_empty() {
            return Self::simplify_trace_terminal(f.as_view());
        }

        let mut factor_kinds = DiracFactorKinds::default();
        let factors = factors
            .iter()
            .map(|factor| {
                let factor = DiracFactor::parse(*factor);
                factor_kinds.observe(factor);
                factor
            })
            .collect::<Vec<_>>();

        if factor_kinds.has_projector
            && has_four_dimensional_trace_rep(rep)
            && factors.iter().all(|factor| match factor {
                DiracFactor::Gamma { dimension, .. } => dimension.is_some_and(is_four_dimension),
                DiracFactor::ChargeConjugation(_) | DiracFactor::Other(_) => false,
                _ => true,
            })
        {
            for (position, factor) in factors.iter().enumerate() {
                let sign = match factor {
                    DiracFactor::ProjectorPlus(_) => 1,
                    DiracFactor::ProjectorMinus(_) => -1,
                    _ => continue,
                };
                // P± = (1 ± gamma5)/2. Expand one projector per pass so the
                // existing fixed point and gamma5 trace rules own the reduction.
                let mut rest = Vec::with_capacity(factors.len());
                Self::extend_factors(&mut rest, &factors[..position]);
                Self::extend_factors(&mut rest, &factors[position + 1..]);
                let ordinary = Self::trace_or_terminal(rep, rest.clone());
                rest.insert(position, gamma5_factor());
                let axial = Self::trace_or_terminal(rep, rest);
                return Some((ordinary + Atom::num(sign) * axial) / 2);
            }
        }

        if factor_kinds.has_special_trace_pair()
            && let Some(rewritten) = Self::simplify_special_trace_pair(rep, &factors)
        {
            return Some(rewritten);
        }

        if factor_kinds.has_conjugation_candidate()
            && has_four_dimensional_trace_rep(rep)
            && let Some((sign, factors)) = Self::conjugate_special_dirac_factor(&factors)
        {
            return Some(Atom::num(sign) * Self::trace_or_terminal(rep, factors));
        }

        if factor_kinds.has_gamma5
            && let Some(rewritten) = Self::simplify_gamma5_trace_node(rep, &factors)
        {
            return Some(rewritten);
        }

        let trace_mink_indices =
            Self::gamma_mink_index_sequence_for(TRACE_GAMMA_RECURSION, &factors)?;

        if factors.len() % 2 == 1 {
            // Tr(gamma(mu_1)...gamma(mu_{2n+1})) -> 0.
            return Some(Atom::Zero);
        }

        // Cyclicity lets the preferred repeated pair cross the stored boundary.
        // Prefer a one-word Chisholm reduction before a branching contraction.
        let rotated = Self::rotate_trace_to_repeated_pair(&factors);
        let word = DiracWord::Trace(rep);
        if let Some(reduced) = Self::contract_adjacent_gamma_pair(word, &rotated)
            .or_else(|| Self::four_dim_chisholm_contraction(word, &rotated))
            .or_else(|| Self::bubble_repeated_gamma_towards_contraction(word, &rotated))
        {
            return Some(reduced);
        }
        let four_dimensional = has_four_dimensional_trace_rep(rep)
            && Self::gamma_mink_index_sequence_for(FOUR_DIM_CHISHOLM, &factors).is_some();
        if four_dimensional
            && let Some(result) = trace_kernel::evaluate(
                &trace_mink_indices,
                false,
                trace_kernel::TraceOutput::Expanded,
            )
        {
            return Some(result);
        }

        if !four_dimensional {
            let terminal = trace!(rep; std::iter::empty::<Atom>());
            let trace_unit = Self::simplify_trace_terminal(terminal.as_view())?;
            return Some(trace_kernel::evaluate_generic(
                &trace_mink_indices,
                trace_unit.as_view(),
                false,
                trace_kernel::TraceOutput::Factored,
            ));
        }

        let first = trace_mink_indices[0];
        let mut sum = Atom::Zero;

        // Standard recursive even trace formula, with one-based positions:
        // Tr(g1...gn) =
        //   sum_{k=2..n} (-1)^k g(1,k) Tr(g2...g_{k-1}g_{k+1}...gn).
        for i in 1..factors.len() {
            let mu_i = trace_mink_indices[i];
            let sign = if i % 2 == 1 { 1 } else { -1 };

            let mut rest = Vec::with_capacity(factors.len() - 2);
            Self::extend_factors(&mut rest, &factors[1..i]);
            Self::extend_factors(&mut rest, &factors[i + 1..]);

            let rest_trace = if rest.is_empty() {
                let terminal_trace = trace!(rep; std::iter::empty::<Atom>());
                Self::simplify_trace_terminal(terminal_trace.as_view())?
            } else {
                trace!(rep; rest)
            };

            let term = g!(first, mu_i) * rest_trace;
            if sign == 1 {
                sum += term;
            } else {
                sum -= term;
            }
        }

        Some(sum)
    }

    /// Adjacent pairs remain cheapest. Among the rest, an explicit 4D pair
    /// with an odd gamma-only interior reduces to one word instead of two.
    /// Applicable even interiors also precede arcs crossing gamma5 or another
    /// unsupported factor. Rotation and contraction use this same ordering.
    fn trace_pair_priority(
        factors: &[DiracFactor<'_>],
        start: usize,
        distance: usize,
    ) -> (u8, usize) {
        let priority = if distance == 1 {
            0
        } else if factors[start]
            .gamma_mink_index(FOUR_DIM_CHISHOLM)
            .is_some_and(is_minkowski_slot)
            && factors
                .iter()
                .cycle()
                .skip(start + 1)
                .take(distance - 1)
                .all(|factor| factor.gamma_mink_index(FOUR_DIM_CHISHOLM).is_some())
        {
            if distance.is_multiple_of(2) { 1 } else { 2 }
        } else {
            3
        };
        (priority, distance)
    }

    fn rotate_trace_to_repeated_pair<'a>(factors: &[DiracFactor<'a>]) -> Vec<DiracFactor<'a>> {
        let mut best = None;
        for (i, left) in factors.iter().enumerate() {
            for (j, right) in factors.iter().enumerate().skip(i + 1) {
                if Self::mink_index_pair(GAMMA_ANTICOMMUTATION, left, right)
                    .is_some_and(|(a, b)| a == b)
                {
                    let distance = j - i;
                    for candidate in [(i, distance), (j, factors.len() - distance)] {
                        let score = Self::trace_pair_priority(factors, candidate.0, candidate.1);
                        if best.is_none_or(|(_, best_score)| score < best_score) {
                            best = Some((candidate, score));
                        }
                    }
                }
            }
        }
        let start = best.map_or(0, |((start, _), _)| start);
        Self::cyclic_from_position(factors, start)
    }

    fn simplify_special_trace_pair(rep: AtomView<'_>, factors: &[DiracFactor<'_>]) -> Option<Atom> {
        if !has_four_dimensional_trace_rep(rep) {
            return None;
        }

        if matches!(factors, [DiracFactor::ChargeConjugation(_)]) {
            return Some(Atom::Zero);
        }

        for (i, pair) in factors.windows(2).enumerate() {
            let [left, right] = pair else {
                unreachable!("windows(2) always yields pairs")
            };
            match (left, right) {
                (DiracFactor::ChargeConjugation(_), DiracFactor::ChargeConjugation(_)) => {
                    let mut rest = Vec::with_capacity(factors.len() - 2);
                    Self::extend_factors(&mut rest, &factors[..i]);
                    Self::extend_factors(&mut rest, &factors[i + 2..]);
                    return Some(-Self::trace_or_terminal(rep, rest));
                }
                (DiracFactor::Gamma5(_), DiracFactor::Gamma5(_))
                | (DiracFactor::Gamma0(_), DiracFactor::Gamma0(_)) => {
                    // Four-dimensional involutions inside a trace:
                    // Tr(... gamma5 gamma5 ...) -> Tr(...)
                    // Tr(... gamma0 gamma0 ...) -> Tr(...).
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
        // Use trace cyclicity to put gamma5 first, then evaluate
        // Tr(gamma5 gamma(mu_1)...gamma(mu_n)).
        let parsed_after_gamma5 = Self::cyclic_without_position(factors, gamma5_position);
        let mink_indices =
            Self::gamma_mink_index_sequence_for(TRACE_GAMMA5_RECURSION, &parsed_after_gamma5)?;

        if mink_indices.len() < 4 || mink_indices.len() % 2 == 1 {
            // With one gamma5, traces with fewer than four gammas or an odd
            // number of ordinary gammas vanish.
            return Some(Atom::Zero);
        }

        let rotated = Self::rotate_trace_to_repeated_pair(factors);
        let word = DiracWord::Trace(rep);
        if let Some(reduced) = Self::contract_adjacent_gamma_pair(word, &rotated)
            .or_else(|| Self::four_dim_chisholm_contraction(word, &rotated))
            .or_else(|| Self::bubble_repeated_gamma_towards_contraction(word, &rotated))
        {
            return Some(reduced);
        }
        if let Some(result) =
            trace_kernel::evaluate(&mink_indices, true, trace_kernel::TraceOutput::Expanded)
        {
            return Some(result);
        }

        // The ordinary pairing recurrence does not apply with one gamma5.
        // Reduce a triple using the same 4D identity as the short kernels.
        Self::four_dim_three_gamma_epsilon_expansion(word, factors)
    }

    fn reduce_gamma5_trace_pair(
        rep: AtomView<'_>,
        factors: &[DiracFactor<'_>],
        first_gamma5_position: usize,
    ) -> Option<Atom> {
        // Rotate one gamma5 to the front. If every intervening factor
        // anticommutes with gamma5, move the second gamma5 next to it:
        // Tr(gamma5 A_1...A_m gamma5 B) -> (-1)^m Tr(A_1...A_m B).
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
        // Tr_rep(1) -> dim(rep).
        let trace = trace.to_owned();
        let simplified = trace.replace_multiple_repeat(TRACE_TERMINALS.as_ref());
        (simplified != trace).then_some(simplified)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{gamma, gamma5, test_support::test_initialize};
    use spenso::slot;

    #[test]
    fn exact_scalar_trace_products_preserve_the_terminal_result() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let (a, b) = (
            slot!(r.mink_d, a).into_atom(),
            slot!(r.mink_d, b).into_atom(),
        );
        let p = momenta(&r.mink_d.to_symbolic([]));
        let input = trace!(&spin; [&a, &p[0], &b, &p[1], &a, &p[2], &b, &p[3]]
            .map(|argument| gamma!(argument)));
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        let decorated = &spectator * &input;
        let expected = (&spectator * input.simplify_gamma()).simplify_gamma();
        let output = decorated.simplify_gamma();
        assert_eq!(output, expected);
        assert!(matches!(output.as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == spectator.as_view())));
        assert_eq!(output.simplify_gamma(), output);
        assert_eq!(
            decorated
                .simplify_gamma_with(GammaSimplifySettings::default().without_trace_evaluation()),
            decorated
        );
    }

    #[test]
    fn terminal_trace_products_are_closed_with_wildcard_and_tagged_scalar_leaves() {
        test_initialize();
        let vector = T.rank_one_tensor_symbol("idenso::terminal_closure_vector");
        let dimension = Atom::var(symbolica::symbol!("terminal_closure_dimension"));
        let four = Atom::num(4);
        let index = Atom::var(symbolica::symbol!("terminal_closure_index"));
        let wildcard = Atom::var(symbolica::symbol!("terminal_closure_wild_"));
        let metadata = Atom::var(vector);
        // These are scalar variables, even when their symbols also name tensor
        // functions. Keep the entire scalar numerator as one product factor.
        let spectator = (symbolica::parse_lit!(x + y)
            + Atom::add_many(
                [
                    T.chain,
                    T.trace,
                    T.bracket,
                    *crate::epsilon::EPSILON_SYMBOL,
                    vector,
                ]
                .map(Atom::var),
            ))
        .pow(8);
        for (name, [dim, unit, label, parameter], repeated) in [
            ("index", [&dimension, &four, &wildcard, &metadata], true),
            ("dimension", [&wildcard, &four, &index, &metadata], true),
            ("spin", [&dimension, &wildcard, &index, &metadata], true),
            ("metadata", [&dimension, &four, &index, &wildcard], true),
            (
                "free index",
                [&dimension, &four, &wildcard, &metadata],
                false,
            ),
        ] {
            let compact = symbolica::function!(*MINKOWSKI_SYMBOL, dim);
            let a = symbolica::function!(*MINKOWSKI_SYMBOL, dim, label);
            let b = if repeated {
                a.clone()
            } else {
                symbolica::function!(
                    *MINKOWSKI_SYMBOL,
                    dim,
                    symbolica::symbol!("terminal_closure_other_index")
                )
            };
            let p = symbolica::function!(vector, parameter, &compact);
            let q = spenso::q!(&compact);
            let input = trace!(symbolica::function!(*BISPINOR_SYMBOL, unit);
                [&a, &p, &b, &q].map(|argument| gamma!(argument)));
            let admitted = DiracSimplifier::new(&GammaSimplifySettings::default())
                .evaluate_terminal_trace::<true>(input.as_view())
                .unwrap_or_else(|| panic!("{name}: expected scalar-context admission"));
            let standalone = input.simplify_gamma();
            assert_eq!(admitted, standalone, "{name}: terminal result");
            let raw = &spectator * &standalone;
            // Wildcards may keep a diagonal metric unevaluated. Its eliminated
            // dummy has no partner, so even that raw output needs no cleanup.
            assert_eq!(raw.simplify_gamma(), raw, "{name}: full cleanup");
            let contextual = (&spectator * &input).simplify_gamma();
            assert_eq!(contextual, raw, "{name}: contextual result");
            assert_eq!(contextual.simplify_gamma(), contextual, "{name}: rerun");
            assert!(
                matches!(contextual.as_view(), AtomView::Mul(product)
                if product.iter().any(|factor| factor == spectator.as_view())),
                "{name}: factored scalar numerator"
            );
        }
    }

    #[test]
    fn scalar_trace_products_preserve_component_callback_order() {
        let r = test_initialize();
        fn component(value: AtomView<'_>, out: &mut Settable<Atom>) {
            if let AtomView::Fun(vector) = value
                && let Some(AtomView::Fun(slot)) = vector.iter().last()
                && slot.get_symbol() == *MINKOWSKI_SYMBOL
                && slot.get_nargs() == 2
            {
                **out = Atom::num(1);
            }
        }
        let p = spenso::vector_symbol!("idenso::scalar_context_callback_p", norm = component);
        let q = spenso::vector_symbol!("idenso::scalar_context_callback_q", norm = component);
        let compact = r.mink_d.to_symbolic([]);
        let p = symbolica::function!(p, &compact);
        let q = symbolica::function!(q, &compact);
        let a = slot!(r.mink_d, a).into_atom();
        let input =
            trace!(r.bis4.to_symbolic([]); [&a, &p, &a, &q].map(|argument| gamma!(argument)));
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        let dimension = mink_slot_dimension(a.as_view()).unwrap();
        // Component callbacks are applied by the existing full trace route.
        let expected = &spectator * (Atom::num(8) - Atom::num(4) * dimension * g!(&p, &q));
        let output = (&spectator * &input).simplify_gamma();
        assert_eq!(output, expected);
        assert_eq!(output.simplify_gamma(), output);
        assert!(((&output / &spectator) - input.simplify_gamma()).expand() != Atom::Zero);
    }

    #[test]
    fn scalar_trace_products_preserve_rounded_operation_order() {
        let r = test_initialize();
        let dimension = r.mink_d.to_symbolic([]).as_view().to_owned();
        let d = minkowski_dimension(dimension.as_view()).unwrap().to_owned();
        let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
        let four = Atom::num(4);
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        let inert = symbolica::function!(symbolica::symbol!("scalar_context_inert"), Atom::num(0));
        let word = |dimension: &Atom, unit: &Atom| {
            let compact = symbolica::function!(*MINKOWSKI_SYMBOL, dimension);
            let p = momenta(&compact);
            let a = symbolica::function!(
                *MINKOWSKI_SYMBOL,
                dimension,
                symbolica::symbol!("scalar_context_a")
            );
            let b = symbolica::function!(
                *MINKOWSKI_SYMBOL,
                dimension,
                symbolica::symbol!("scalar_context_b")
            );
            trace!(symbolica::function!(*BISPINOR_SYMBOL, unit);
                [&a, &p[0], &b, &p[1], &a, &p[2], &b, &p[3]].map(|argument| gamma!(argument)))
        };
        // These independently guard both compact-only metadata positions;
        // an explicit slot elsewhere must not mask a missing admission check.
        let compact = r.mink_d.to_symbolic([]);
        let p = momenta(&compact);
        let metadata = symbolica::function!(
            T.rank_one_tensor_symbol("idenso::rounded_compact_metadata"),
            &rounded,
            &compact
        );
        let leading_metadata = trace!(r.bis4.to_symbolic([]);
            [&metadata, &p[0]].map(|argument| gamma!(argument)));
        let rounded_compact = symbolica::function!(*MINKOWSKI_SYMBOL, &rounded);
        let p = momenta(&rounded_compact);
        let compact_dimension = trace!(r.bis4.to_symbolic([]);
            [&p[0], &p[1]].map(|argument| gamma!(argument)));
        for input in [leading_metadata, compact_dimension] {
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<true>(input.as_view())
                    .is_none()
            );
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .is_some()
            );
        }
        for (name, input) in [
            ("trace unit", &spectator * word(&d, &rounded)),
            ("dimension", &spectator * word(&rounded, &four)),
            ("spectator", &rounded * &spectator * word(&d, &four)),
        ] {
            // An inert function keeps the established full-route ordering.
            // Compare exact floating coefficients on both calls, without a tolerance.
            let reference = (&inert * &input).simplify_gamma() / &inert;
            let output = input.simplify_gamma();
            assert_eq!(output, reference, "{name}: first call");
            assert_eq!(
                output.simplify_gamma(),
                reference.simplify_gamma(),
                "{name}: rerun"
            );
        }
    }

    #[test]
    fn cleanup_callbacks_can_introduce_new_traces() {
        let reps = test_initialize();
        let compact = reps.mink4.to_symbolic([]);
        let spin = reps.bis4.to_symbolic([]);
        let mink = *MINKOWSKI_SYMBOL;
        let inner = trace!(
            &spin,
            gamma!(slot!(reps.mink4, rho)),
            gamma!(slot!(reps.mink4, sigma))
        );
        let emitted = inner.clone();
        let inner_vector = spenso::vector_symbol!(
            "idenso::cleanup_inner_vector",
            norm = move |value, out| {
                if let AtomView::Fun(function) = value
                    && let Some(AtomView::Fun(slot)) = function.iter().last()
                    && slot.get_symbol() == mink
                    && slot.get_nargs() == 1
                {
                    **out = emitted.clone();
                }
            }
        );
        let other = T.rank_one_tensor_symbol("idenso::cleanup_other_vector");
        let nested = symbolica::function!(inner_vector, symbolica::function!(other, &compact));
        let component = nested.clone();
        let outer_vector = spenso::vector_symbol!(
            "idenso::cleanup_outer_vector",
            norm = move |value, out| {
                if let AtomView::Fun(function) = value
                    && let Some(AtomView::Fun(slot)) = function.iter().last()
                    && slot.get_symbol() == mink
                    && slot.get_nargs() == 2
                {
                    **out = component.clone();
                }
            }
        );
        // The outer trace emits a nested vector. Dot normalization then
        // constructs its compact component, whose callback introduces a trace.
        assert!(!nested.contains_symbol(T.trace));
        assert!(nested.normalize_dots().contains_symbol(T.trace));
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        let input = &spectator
            * trace!(
                &spin,
                gamma!(symbolica::function!(outer_vector, &compact)),
                gamma!(slot!(reps.mink4, mu))
            );
        let expected = &spectator
            * 4
            * g!(
                inner.simplify_gamma(),
                symbolica::function!(other, &compact)
            );
        let output = input.simplify_gamma();
        assert!(!output.contains_symbol(T.trace));
        assert_eq!(output, expected);
        assert_eq!(output.simplify_gamma(), output);
    }

    #[test]
    fn rewrite_guard_retains_empty_trace_and_non_gamma_chain_identities() {
        let r = test_initialize();
        let trace = trace!(r.bis4.to_symbolic([]); std::iter::empty::<Atom>());
        let start = slot!(r.bis4, a).into_atom();
        let end = slot!(r.bis4, b).into_atom();
        let coefficient = symbolica::parse_lit!((x + y) ^ 8);
        let settings = GammaSimplifySettings::default();
        assert_eq!(
            settings.rewrite_expression(&coefficient * &trace),
            &coefficient * Atom::num(4)
        );
        assert_eq!(
            settings
                .without_trace_evaluation()
                .rewrite_expression(&coefficient * &trace),
            &coefficient * trace
        );
        for chain in [
            chain!(&start, &end),
            chain!(&start, &end, gamma5!(), gamma5!()),
        ] {
            assert_eq!(
                settings.rewrite_expression(settings.rewrite_expression(&coefficient * chain)),
                &coefficient * id_atom(start.as_view(), end.as_view())
            );
        }
        let terminal = &coefficient * g!(slot!(r.mink4, mu), slot!(r.mink4, nu));
        assert_eq!(settings.rewrite_expression(terminal.clone()), terminal);
    }

    #[test]
    fn free_axial_shortcut_matches_complete_pass_at_each_gamma5_position() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let slots: Vec<_> = (0..14)
            .map(|index| r.mink4.pattern(Atom::num(index)))
            .collect();
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        for length in 0..=14 {
            let positions: Vec<_> = if length <= 8 {
                (0..=length).collect()
            } else {
                vec![0, length / 2, length]
            };
            for position in positions {
                let mut factors: Vec<_> = slots[..length].iter().map(|slot| gamma!(slot)).collect();
                factors.insert(position, gamma5!());
                let input = trace!(&spin; factors);
                let shortcut = DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .expect("a distinct-index short axial trace is terminal");
                // A factorized scalar spectator forces the full existing
                // pipeline, without expanding or changing the trace word.
                let fallback = &spectator * &input;
                assert!(
                    DiracSimplifier::new(&GammaSimplifySettings::default())
                        .evaluate_terminal_trace::<false>(fallback.as_view())
                        .is_none()
                );
                assert_eq!(
                    fallback.simplify_gamma(),
                    &spectator * shortcut.expand(),
                    "length {length}, gamma5 position {position}"
                );
                assert_eq!(input.simplify_gamma(), shortcut);
                assert_eq!(
                    input.simplify_gamma_with(
                        GammaSimplifySettings::default().without_trace_evaluation()
                    ),
                    input
                );
            }
        }
    }

    #[test]
    fn free_axial_shortcut_leaves_contractions_and_unsupported_words_to_full_pass() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let mu = slot!(r.mink4, mu).into_atom();
        let nu = slot!(r.mink4, nu).into_atom();
        let repeated = trace!(&spin, gamma5!(), gamma!(&mu), gamma!(&mu));
        let two_gamma5 = trace!(&spin, gamma5!(), gamma5!(), gamma!(&mu), gamma!(&nu));
        let slash = spenso::p!(r.mink4.to_symbolic([]));
        let compact = trace!(&spin, gamma5!(), gamma!(slash));
        let generic_dimension = trace!(&spin, gamma5!(), gamma!(slot!(r.mink_d, mu)));
        let long = trace!(&spin;
            [gamma5!()].into_iter().chain((0..16).map(|index| gamma!(r.mink4.pattern(Atom::num(index)))))
        );
        for input in [repeated, two_gamma5, compact, generic_dimension, long] {
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .is_none(),
                "{input}"
            );
        }
    }

    #[test]
    fn generic_trace_kernel_preserves_spin_dimension_and_spectator_factorization() {
        let r = test_initialize();
        let slots: Vec<_> = (0..4)
            .map(|index| r.mink_d.pattern(Atom::num(index)))
            .collect();
        let pairing = g!(&slots[0], &slots[1]) * g!(&slots[2], &slots[3])
            - g!(&slots[0], &slots[2]) * g!(&slots[1], &slots[3])
            + g!(&slots[0], &slots[3]) * g!(&slots[1], &slots[2]);
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        for spin in [r.bis4.to_symbolic([]), r.bis_d.to_symbolic([])] {
            let input = trace!(&spin; slots.iter().map(|slot| gamma!(slot)));
            let unit = bispinor_dimension(spin.as_view()).unwrap();
            let expected = pairing.as_view() * unit;
            let result = input.simplify_gamma();
            assert!((&result - expected).expand().is_zero());
            assert_eq!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view()),
                Some(result.clone())
            );
            assert_eq!(
                (&spectator * input).simplify_gamma(),
                &spectator * &result,
                "the outer scalar numerator remains factored"
            );
            assert_eq!(result.simplify_gamma(), result);
        }
    }

    #[test]
    fn generic_trace_dispatch_preserves_dimension_and_axial_boundaries() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let a = slot!(r.mink_d, a).into_atom();
        let b = slot!(r.mink4, b).into_atom();
        for input in [
            trace!(&spin, gamma!(&a), gamma!(&b)),
            trace!(&spin, gamma5!(), gamma!(&a)),
        ] {
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .is_none()
            );
            assert_eq!(input.simplify_gamma(), input);
        }
        let odd = trace!(&spin;
            (0..5).map(|index| gamma!(r.mink_d.pattern(Atom::num(index))))
        );
        assert_eq!(odd.simplify_gamma(), Atom::Zero);
        assert_eq!(
            odd.simplify_gamma_with(GammaSimplifySettings::default().without_trace_evaluation()),
            odd
        );
        let symbolic_spin = trace!(r.bis_d.to_symbolic([]);
            (0..10).map(|index| gamma!(r.mink4.pattern(Atom::num(index))))
        );
        // The 4D kernel fixes Tr(1)=4. A symbolic spin dimension must instead
        // multiply the dimension-generic pairing formula.
        assert_eq!(symbolic_spin.simplify_gamma().expand().nterms(), 945);
    }

    #[test]
    fn generic_trace_kernel_keeps_external_metric_and_vector_contractions_complete() {
        use spenso::network::parsing::AtomStructureExt;

        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let slots: Vec<_> = (0..10)
            .map(|index| r.mink_d.pattern(Atom::num(index)))
            .collect();
        let input = g!(&slots[0], &slots[1]) * trace!(&spin; slots.iter().map(|slot| gamma!(slot)));
        let result = input.simplify_gamma();
        let tail = trace!(&spin; slots[2..].iter().map(|slot| gamma!(slot)));
        let expected = tail.simplify_gamma() * mink_slot_dimension(slots[0].as_view()).unwrap();
        assert!((&result - expected).expand().is_zero());
        assert!(!result.has_repeated_explicit_indices());

        let p = spenso::p!(r.mink_d.to_symbolic([]));
        let q = spenso::q!(r.mink_d.to_symbolic([]));
        let input = spenso::p!(&slots[0])
            * spenso::q!(&slots[1])
            * trace!(&spin; slots[..4].iter().map(|slot| gamma!(slot)));
        let expected: Atom = 4
            * (g!(&p, &q) * g!(&slots[2], &slots[3]) - g!(&p, &slots[2]) * g!(&q, &slots[3])
                + g!(&p, &slots[3]) * g!(&q, &slots[2]));
        let result = input.simplify_gamma();
        assert!((&result - expected).expand().is_zero());
        assert!(!result.has_repeated_explicit_indices());
    }

    #[test]
    fn compact_trace_words_preserve_scalar_recurrences_and_spin_dimension() {
        let r = test_initialize();
        for representation in [r.mink4.to_symbolic([]), r.mink_d.to_symbolic([])] {
            let p = spenso::p!(&representation);
            let q = spenso::q!(&representation);
            let pp = g!(&p, &p);
            let qq = g!(&q, &q);
            let pq = g!(&p, &q);
            for spin in [r.bis4.to_symbolic([]), r.bis_d.to_symbolic([])] {
                let unit = bispinor_dimension(spin.as_view()).unwrap();
                let paired = trace!(&spin;
                    [&p, &p].into_iter().chain(std::iter::repeat_n(&q, 12)).map(|p| gamma!(p))
                );
                assert_eq!(paired.simplify_gamma(), unit.to_owned() * &pp * qq.pow(6));

                // For A=p/ q/, A²-2(p.q)A+p²q²=0. This independent scalar
                // recurrence certifies all seven alternating pairs exactly.
                let mut previous = unit.to_owned();
                let mut expected = unit * pq.as_view();
                for pairs in 1..=7 {
                    let input = trace!(&spin;
                        (0..2 * pairs).map(|i| gamma!(if i % 2 == 0 { &p } else { &q }))
                    );
                    let result = input.simplify_gamma();
                    assert!((&result - &expected).expand().is_zero());
                    assert_eq!(result.simplify_gamma(), result);
                    assert!(
                        DiracSimplifier::new(&GammaSimplifySettings::default())
                            .evaluate_terminal_trace::<false>(input.as_view())
                            .is_some()
                    );
                    let next = Atom::num(2) * &pq * &expected - &pp * &qq * &previous;
                    previous = expected;
                    expected = next;
                }
            }
        }
    }

    #[test]
    fn compact_trace_shortcut_keeps_metadata_and_rejects_noncanonical_arguments() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let rep = r.mink4.to_symbolic([]);
        let head = T.rank_one_tensor_symbol("idenso::compact_trace_metadata");
        let vector = |arguments: &[Atom]| FunctionBuilder::new(head).add_args(arguments).finish();
        let p = vector(&[Atom::num(1), rep.clone()]);
        let q = vector(&[Atom::num(2), rep.clone()]);
        let input = trace!(&spin, gamma!(&p), gamma!(&q));
        assert_eq!(input.simplify_gamma(), Atom::num(4) * g!(&p, &q));
        let odd = trace!(&spin, gamma!(&p), gamma!(&q), gamma!(&p));
        assert!(odd.simplify_gamma().is_zero());

        let unknown = FunctionBuilder::new(symbolica::symbol!("compact_trace_unknown"))
            .add_arg(&rep)
            .finish();
        let malformed = vector(&[rep.clone(), rep.clone()]);
        let indexed = spenso::p!(r.mink4.pattern(Atom::num(91)));
        for argument in [unknown, malformed, indexed] {
            let input = trace!(&spin, gamma!(&p), gamma!(argument));
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .is_none()
            );
        }
        let mixed = trace!(
            &spin,
            gamma!(&p),
            gamma!(spenso::q!(r.mink_d.to_symbolic([])))
        );
        assert!(
            DiracSimplifier::new(&GammaSimplifySettings::default())
                .evaluate_terminal_trace::<false>(mixed.as_view())
                .is_none()
        );
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        assert!(
            DiracSimplifier::new(&GammaSimplifySettings::default())
                .evaluate_terminal_trace::<false>((&spectator * &input).as_view())
                .is_none()
        );
        assert_eq!(
            input.simplify_gamma_with(GammaSimplifySettings::default().without_trace_evaluation()),
            input
        );
    }

    fn momenta(rep: &Atom) -> [Atom; 8] {
        std::array::from_fn(|i| {
            FunctionBuilder::new(T.rank_one_tensor_symbol(&format!("idenso::trace_order::p{i}")))
                .add_arg(rep)
                .finish()
        })
    }

    #[test]
    fn trace_prefers_an_odd_interior_over_a_shorter_branching_pair() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let (a, b) = (slot!(r.mink4, a).into_atom(), slot!(r.mink4, b).into_atom());
        let p = momenta(&r.mink4.to_symbolic([]));
        let factors = [
            &a, &p[0], &b, &p[1], &p[2], &a, &p[3], &p[4], &b, &p[5], &p[6], &p[7],
        ]
        .map(|index| gamma!(index));
        let expected_first = Atom::num(-2)
            * trace!(&spin; [
                &a, &p[0], &p[4], &p[3], &a, &p[2], &p[1], &p[5], &p[6], &p[7],
            ].map(|index| gamma!(index)));
        let parsed = factors
            .iter()
            .map(|factor| DiracFactor::parse(factor.as_view()))
            .collect::<Vec<_>>();
        let rotated = DiracSimplifier::rotate_trace_to_repeated_pair(&parsed);
        assert_eq!(
            DiracSimplifier::four_dim_chisholm_contraction(
                DiracWord::Trace(spin.as_view()),
                &rotated,
            ),
            Some(expected_first.clone()),
        );
        // Rotation and the actual selector must agree on the one-word route.
        let input = trace!(&spin; factors.iter());
        assert!(
            (input.simplify_gamma() - expected_first.simplify_gamma())
                .expand()
                .is_zero()
        );
        // Open chains retain their existing shortest-pair ordering.
        let (start, end) = (slot!(r.bis4, i).into_atom(), slot!(r.bis4, j).into_atom());
        assert_eq!(
            DiracSimplifier::repeated_four_dim_gamma_pair(
                DiracWord::Chain(start.as_view(), end.as_view()),
                &parsed,
            ),
            Some((0, 5)),
        );
    }

    #[test]
    fn contracted_terminal_trace_preserves_dimensions_and_free_slots() {
        use spenso::structure::dimension::Dimension;

        let r = test_initialize();
        for dimension in [
            Atom::num(4),
            Atom::num(6),
            Atom::var(symbolica::symbol!("terminal_trace_D")),
        ] {
            let mink = Minkowski {}.new_rep(Dimension::try_from(dimension.as_view()).unwrap());
            let slots: Vec<_> = (0..3).map(|i| mink.pattern(Atom::num(i))).collect();
            for spin in [r.bis4.to_symbolic([]), r.bis_d.to_symbolic([])] {
                let unit = bispinor_dimension(spin.as_view()).unwrap();
                for (word, coefficient) in [
                    ([0, 1, 2, 0], dimension.clone()),
                    ([0, 1, 0, 2], Atom::num(2) - &dimension),
                ] {
                    for start in 0..4 {
                        let input = trace!(&spin; word.iter().cycle().skip(start).take(4)
                            .map(|&i| gamma!(&slots[i])));
                        let expected = coefficient.as_view() * unit * g!(&slots[1], &slots[2]);
                        let result = DiracSimplifier::new(&GammaSimplifySettings::default())
                            .evaluate_terminal_trace::<false>(input.as_view())
                            .expect("the summed pair reduces without emitting tensor metrics");
                        assert_eq!(result, expected, "dimension {dimension}, rotation {start}");
                        assert_eq!(input.simplify_gamma(), result);
                    }
                }
            }
        }
    }

    #[test]
    fn terminal_trace_contracts_branching_words_but_rejects_ambiguous_indices() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        for dimension in [4, 6] {
            let mink = Minkowski {}.new_rep(dimension);
            let a = mink.pattern(Atom::num(10));
            let b = mink.pattern(Atom::num(11));
            let p = momenta(&mink.to_symbolic([]));
            for word in [
                vec![&a, &p[0], &p[1], &a, &p[2], &p[3], &p[4], &p[5]],
                vec![&a, &p[0], &b, &p[1], &a, &p[2], &b, &p[3]],
            ] {
                let input = trace!(&spin; word.into_iter().map(|index| gamma!(index)));
                let result = DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .unwrap();
                result.visitor(&mut |node| {
                    assert!(node != a.as_view() && node != b.as_view());
                    true
                });
            }
            let component = mink.pattern(symbolica::function!(
                spenso::structure::abstract_index::AIND_SYMBOLS.cind,
                Atom::num(0)
            ));
            let fixed_component = trace!(&spin; [&component, &p[0], &p[1], &component, &p[2], &p[3]]
                .map(|index| gamma!(index)));
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(fixed_component.as_view())
                    .is_none(),
                "fixed components do not define Einstein sums"
            );
            let ambiguous = trace!(&spin; [&a, &p[0], &a, &p[1], &a, &p[2]]
                .map(|index| gamma!(index)));
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(ambiguous.as_view())
                    .is_none()
            );
            let input = trace!(&spin; [&a, &p[0], &p[1], &p[2], &a, &p[3], &p[4], &p[5]]
                .map(|index| gamma!(index)));
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .is_some()
            );
        }
    }

    #[test]
    fn branching_terminal_trace_preserves_free_slots_and_scalar_spectators() {
        let r = test_initialize();
        let spin = r.bis_d.to_symbolic([]);
        let slots: Vec<_> = (0..5).map(|i| r.mink_d.pattern(Atom::num(i))).collect();
        let dimension = mink_slot_dimension(slots[0].as_view()).unwrap();
        let unit = bispinor_dimension(spin.as_view()).unwrap();
        let input = trace!(&spin; [0, 1, 2, 0, 3, 4].map(|i| gamma!(&slots[i])));
        let reduced = trace!(&spin; slots[1..].iter().map(|index| gamma!(index))).simplify_gamma();
        let expected = Atom::num(4) * unit * g!(&slots[1], &slots[2]) * g!(&slots[3], &slots[4])
            + (dimension.to_owned() - Atom::num(4)) * reduced;
        let result = input.simplify_gamma();
        assert!((&result - &expected).expand().is_zero());
        result.visitor(&mut |node| {
            assert_ne!(node, slots[0].as_view(), "the summed index must be absent");
            true
        });
        let spectator = symbolica::parse_lit!((x + y) ^ 8);
        let decorated = &spectator * &input;
        assert!(
            DiracSimplifier::new(&GammaSimplifySettings::default())
                .evaluate_terminal_trace::<false>(decorated.as_view())
                .is_none()
        );
        let complete = decorated.simplify_gamma();
        // Remove only the spectator for the standalone polynomial certificate;
        // the compound scalar numerator itself remains unexpanded.
        let body = complete
            .replace(spectator.to_pattern())
            .with(Atom::num(1).to_pattern());
        assert_eq!(complete, &spectator * &body);
        assert!((body - expected).expand().is_zero());
    }

    #[test]
    fn trace_keeps_cyclic_adjacent_pairs_ahead_of_odd_interiors() {
        let r = test_initialize();
        let (a, b) = (slot!(r.mink4, a).into_atom(), slot!(r.mink4, b).into_atom());
        let p = momenta(&r.mink4.to_symbolic([]));
        let factors =
            [&a, &p[0], &b, &p[1], &p[2], &p[3], &b, &p[4], &p[5], &a].map(|index| gamma!(index));
        let parsed = factors
            .iter()
            .map(|factor| DiracFactor::parse(factor.as_view()))
            .collect::<Vec<_>>();
        let rotated = DiracSimplifier::rotate_trace_to_repeated_pair(&parsed);
        assert_eq!(rotated[0].as_view(), factors[9].as_view());
        assert_eq!(rotated[1].as_view(), factors[0].as_view());
    }

    #[test]
    fn trace_rotates_even_only_pairs_across_the_boundary() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let (a, b) = (slot!(r.mink4, a).into_atom(), slot!(r.mink4, b).into_atom());
        let p = momenta(&r.mink4.to_symbolic([]));
        let factors = [&b, &a, &p[1], &p[2], &p[3], &p[4], &p[5], &p[6], &a, &p[7]]
            .map(|index| gamma!(index));
        let parsed = factors
            .iter()
            .map(|factor| DiracFactor::parse(factor.as_view()))
            .collect::<Vec<_>>();
        let rotated = DiracSimplifier::rotate_trace_to_repeated_pair(&parsed);
        let expected = 4 * g!(&p[7], &b) * trace!(&spin; p[1..7].iter().map(|index| gamma!(index)));
        assert_eq!(
            DiracSimplifier::four_dim_chisholm_contraction(
                DiracWord::Trace(spin.as_view()),
                &rotated,
            ),
            Some(expected),
        );
    }

    #[test]
    fn trace_pair_priority_preserves_symbolic_dimensions_and_slashes() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        for (rep, a, b) in [
            (
                r.mink_d.to_symbolic([]),
                slot!(r.mink_d, a).into_atom(),
                slot!(r.mink_d, b).into_atom(),
            ),
            (
                r.mink4.to_symbolic([]),
                spenso::p!(r.mink4.to_symbolic([])),
                spenso::q!(r.mink4.to_symbolic([])),
            ),
        ] {
            let p = momenta(&rep);
            let factors = [
                &a, &p[0], &b, &p[1], &p[2], &a, &p[3], &p[4], &b, &p[5], &p[6], &p[7],
            ]
            .map(|index| gamma!(index));
            let parsed = factors
                .iter()
                .map(|factor| DiracFactor::parse(factor.as_view()))
                .collect::<Vec<_>>();
            let rotated = DiracSimplifier::rotate_trace_to_repeated_pair(&parsed);
            assert_eq!(rotated[0].as_view(), factors[0].as_view());
            assert_eq!(
                DiracSimplifier::repeated_four_dim_gamma_pair(
                    DiracWord::Trace(spin.as_view()),
                    &rotated,
                ),
                None,
            );
        }
    }

    #[test]
    fn axial_trace_can_contract_a_cyclic_odd_interior_without_crossing_gamma5() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let (a, b) = (slot!(r.mink4, a).into_atom(), slot!(r.mink4, b).into_atom());
        let p = momenta(&r.mink4.to_symbolic([]));
        let factors = [
            gamma!(&a),
            gamma!(&p[0]),
            gamma5!(),
            gamma!(&b),
            gamma!(&p[1]),
            gamma!(&p[2]),
            gamma!(&b),
            gamma!(&p[3]),
            gamma!(&p[4]),
            gamma!(&a),
            gamma!(&p[5]),
        ];
        let parsed = factors
            .iter()
            .map(|factor| DiracFactor::parse(factor.as_view()))
            .collect::<Vec<_>>();
        let rotated = DiracSimplifier::rotate_trace_to_repeated_pair(&parsed);
        assert_eq!(
            DiracSimplifier::repeated_four_dim_gamma_pair(
                DiracWord::Trace(spin.as_view()),
                &rotated,
            ),
            Some((0, 2)),
        );
        assert!(
            DiracSimplifier::four_dim_chisholm_contraction(
                DiracWord::Trace(spin.as_view()),
                &rotated,
            )
            .is_some()
        );
    }

    #[test]
    fn external_metrics_contract_before_trace_expansion() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let (a, b, c, d) = (
            slot!(r.mink4, a).into_atom(),
            slot!(r.mink4, b).into_atom(),
            slot!(r.mink4, c).into_atom(),
            slot!(r.mink4, d).into_atom(),
        );
        let p = momenta(&r.mink4.to_symbolic([]));
        let input = g!(&a, &c)
            * g!(&b, &d)
            * trace!(&spin; [
                &a, &p[0], &b, &p[1], &p[2], &c, &p[3], &p[4], &d, &p[5], &p[6], &p[7],
            ].map(|index| gamma!(index)));
        let contracted = input
            .schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions());
        let result = input.simplify_gamma();
        assert!((&result - contracted.simplify_gamma()).expand().is_zero());
        assert_eq!(result.simplify_gamma(), result);
    }

    #[test]
    fn initial_metric_contraction_preserves_inert_traces_and_open_chains() {
        let r = test_initialize();
        for rep in [&r.mink4, &r.mink_d] {
            let (a, b, c) = (
                slot!(rep, a).into_atom(),
                slot!(rep, b).into_atom(),
                slot!(rep, c).into_atom(),
            );
            let spin = r.bis4.to_symbolic([]);
            let input = g!(&a, &b) * trace!(&spin; [gamma!(&b), gamma!(&c)]);
            let expected = trace!(&spin; [gamma!(&a), gamma!(&c)]);
            let settings = GammaSimplifySettings::default().without_trace_evaluation();
            assert_eq!(input.simplify_gamma_with(settings), expected);

            let (start, end) = (slot!(r.bis4, i).into_atom(), slot!(r.bis4, j).into_atom());
            let input = g!(&a, &b) * chain!(&start, &end; [gamma!(&b), gamma!(&c)]);
            let expected = chain!(&start, &end; [gamma!(&a), gamma!(&c)]);
            assert_eq!(input.simplify_gamma(), expected);
        }
    }

    #[test]
    fn chisholm_metric_contraction_reaches_an_idempotent_result() {
        let r = test_initialize();
        let (start, end) = (slot!(r.bis4, i).into_atom(), slot!(r.bis4, j).into_atom());
        let (a, b, c, d) = (
            slot!(r.mink4, a).into_atom(),
            slot!(r.mink4, b).into_atom(),
            slot!(r.mink4, c).into_atom(),
            slot!(r.mink4, d).into_atom(),
        );
        let input = chain!(&start, &end;
            [&a, &b, &c, &a, &b, &d].map(|index| gamma!(index))
        );
        let expected = 4 * chain!(&start, &end; [gamma!(&c), gamma!(&d)]);
        let result = input.simplify_gamma();
        assert_eq!(result, expected);
        assert_eq!(result.simplify_gamma(), result);
    }

    #[test]
    fn evaluated_trace_metrics_contract_into_surviving_chains() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let (start, end) = (slot!(r.bis4, i).into_atom(), slot!(r.bis4, j).into_atom());
        for rep in [&r.mink4, &r.mink_d] {
            let (a, b, c) = (
                slot!(rep, a).into_atom(),
                slot!(rep, b).into_atom(),
                slot!(rep, c).into_atom(),
            );
            let input = trace!(&spin; [gamma!(&a), gamma!(&b)])
                * chain!(&start, &end; [gamma!(&a), gamma!(&c)]);
            let expected = 4 * chain!(&start, &end; [gamma!(&b), gamma!(&c)]);
            let result = input.simplify_gamma();
            assert_eq!(result, expected);
            assert_eq!(result.simplify_gamma(), result);
        }
    }

    #[test]
    fn axial_trace_prefers_an_applicable_arc_over_a_shorter_gamma5_crossing() {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let a = slot!(r.mink4, a).into_atom();
        let p = momenta(&r.mink4.to_symbolic([]));
        let factors = [
            gamma5!(),
            gamma!(&a),
            gamma!(&p[0]),
            gamma!(&p[1]),
            gamma!(&p[2]),
            gamma!(&p[3]),
            gamma!(&p[4]),
            gamma!(&p[5]),
            gamma!(&a),
            gamma!(&p[6]),
            gamma!(&p[7]),
        ];
        let parsed = factors
            .iter()
            .map(|factor| DiracFactor::parse(factor.as_view()))
            .collect::<Vec<_>>();
        let rotated = DiracSimplifier::rotate_trace_to_repeated_pair(&parsed);
        assert_eq!(
            DiracSimplifier::repeated_four_dim_gamma_pair(
                DiracWord::Trace(spin.as_view()),
                &rotated,
            ),
            Some((0, 7)),
        );
        let reduced = DiracSimplifier::four_dim_chisholm_contraction(
            DiracWord::Trace(spin.as_view()),
            &rotated,
        )
        .unwrap();
        assert_eq!(reduced.nterms(), 2);
    }
}
