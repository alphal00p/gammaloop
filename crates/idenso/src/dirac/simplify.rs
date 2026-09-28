#[cfg(test)]
use crate::shorthands::schoonschip::Schoonschip;
use std::{cell::RefCell, sync::LazyLock};

use spenso::{
    chain, g,
    network::tags::SPENSO_TAG as T,
    rep_,
    shadowing::{self, IntoAtom, TensorCollectFilter},
    structure::{
        abstract_index::AbstractIndex,
        partial::PartialStructure,
        representation::{LibraryRep, Minkowski, RepName},
        slot::{DummyAind, ParseableAind, SlotMatch, SlotMatcher},
    },
    trace,
};
use symbolica::{
    atom::{
        Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol, representation::FunView,
    },
    coefficient::CoefficientView,
    id::{Context, Replacement},
    utils::Settable,
};
use symbolica_utils::PatternReplacement;

use crate::{
    W_,
    epsilon::epsilon4,
    representations::Bispinor,
    shorthands::{
        bracket::BracketNormalizer,
        chain::Chain,
        schoonschip::{SchoonschipSettings, SchoonschipWithSettings},
    },
    tensor::{SymbolicTensor, aliases::Definition, inference::TensorInferenceError},
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

/// Select the algebraic output while retaining the shared tensor alias registry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GammaOutput {
    /// Apply Clifford identities and the configured trace evaluation.
    #[default]
    Reduced,
    /// Collect spinor words without reducing their gamma algebra.
    Chains,
}

/// Settings for the chain-based Dirac simplifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GammaSimplifySettings {
    /// Select reduced gamma algebra or collected spinor chains.
    pub output: GammaOutput,
    /// Rewrite conjugated Dirac matrices into gamma-zero sandwiches first.
    pub conjugate: bool,
    /// Apply gamma-zero factoring and repeated gamma-zero identities first.
    pub gamma0: bool,
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
            output: GammaOutput::Reduced,
            conjugate: false,
            gamma0: false,
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

    fn rewrite_expression(&self, expr: Atom, output: trace_kernel::TraceOutput<'_>) -> Atom {
        // The shared scheduler or the full raw pass already observed an
        // eligible domain. Do not repeat its head scan on the unchanged input.
        expr.replace_map(|a, b, c| self.rewrite_node(a, b, c, output))
    }

    fn rewrite_node(
        &self,
        arg: AtomView,
        _context: &Context,
        out: &mut Settable<'_, Atom>,
        output: trace_kernel::TraceOutput<'_>,
    ) {
        let AtomView::Fun(f) = arg else {
            return;
        };

        let simplifier = DiracSimplifier {
            settings: self,
            output,
        };
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

/// Borrowed inputs admitted by the terminal trace evaluator.
struct TerminalTrace<'a> {
    representation: AtomView<'a>,
    indices: Vec<AtomView<'a>>,
    axial: bool,
    repeated: usize,
    four_dimensional: bool,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct DiracSimplifier<'settings> {
    settings: &'settings GammaSimplifySettings,
    output: trace_kernel::TraceOutput<'settings>,
}

impl SymbolicTensor<PartialStructure> {
    pub(crate) fn simplify_gamma_parts(
        &self,
        registry: &[Definition],
        settings: GammaSimplifySettings,
    ) -> Result<(Self, Vec<Definition>), TensorInferenceError> {
        let representation = Bispinor {}.into();
        let standalone_trace = matches!(self.expression.as_view(), AtomView::Fun(function)
            if function.get_symbol() == T.trace);
        self.collect_with_map(
            None,
            registry,
            |atom| TensorCollectFilter::Reps([representation]).matches(atom),
            |selected, _registry, complete| {
                let definitions = RefCell::new(trace_kernel::TraceDefinitions::default());
                let mut scalar_callback = false;
                let mut component_callback = false;
                selected.expression.visitor(&mut |node| {
                    if let AtomView::Fun(function) = node {
                        let symbol = function.get_symbol();
                        let callback = symbol.get_normalization_function().is_some()
                            || symbol.get_evaluation_info().is_some();
                        scalar_callback |= symbol.is_scalar() && callback;
                        component_callback |= symbol.has_tag(&T.rank1) && callback;
                    }
                    true
                });
                let contextual_callback = !standalone_trace && (scalar_callback || component_callback);
                if contextual_callback && !complete {
                    return Ok((selected, Vec::new()));
                }
                // A scalar normalizer observes the actual child result, as in
                // the former local Factored rewrite. Literal aliases cannot be
                // substituted for that callback input. Ordinary cyclic trace
                // wrappers retain the aliased emission path.
                let output = if scalar_callback || contextual_callback {
                    trace_kernel::TraceOutput::Factored
                } else {
                    trace_kernel::TraceOutput::Aliased(&definitions)
                };
                let prepared =
                    DiracSimplifier::new(&settings).prepare(selected.expression.as_view());
                let expression = BracketNormalizer::normalize(prepared.as_view())
                    .chainify(representation)
                    .join_chains(representation);
                let expression = BracketNormalizer::normalize(expression.as_view());
                // Finish joining the selected spinor factors before Clifford
                // reduction can turn an open prefix into a sum of chains. Such
                // a sum can hide the cycle closed by a later selected factor.
                let expression = if settings.output == GammaOutput::Chains || !complete {
                    expression
                } else if contextual_callback {
                    // Component normalizers can temporarily remove one port
                    // before the next local trace identity removes its partner.
                    // Finish this original selected atom context before typed
                    // publication. The shared scheduler remains the only owner
                    // of contractions and cross-domain fixed points.
                    let limit = crate::tensor::simplification::SimplifySettings::default().max_passes;
                    let mut current = expression;
                    let mut finished = false;
                    for _ in 0..limit {
                        let next = settings.rewrite_expression(current.clone(), output);
                        if next == current {
                            finished = true;
                            break;
                        }
                        current = next;
                    }
                    if !finished {
                        return Err(TensorInferenceError::Invalid(format!(
                            "callback-sensitive gamma kernel did not stabilize within {limit} passes"
                        )));
                    }
                    current
                } else {
                    let simplifier = DiracSimplifier { settings: &settings, output };
                    let terminal = settings.evaluate_traces.then(|| {
                        if standalone_trace {
                            simplifier.evaluate_terminal_trace::<false>(expression.as_view())
                        } else {
                            simplifier.evaluate_terminal_trace::<true>(expression.as_view())
                        }
                    }).flatten();
                    terminal.unwrap_or_else(|| settings.rewrite_expression(expression, output))
                };
                let definitions = definitions.into_inner().into_definitions()?;
                Ok((selected.with_rewritten_expression(expression)?, definitions))
            },
        )
    }
}

impl<'settings> DiracSimplifier<'settings> {
    pub(crate) fn new(settings: &'settings GammaSimplifySettings) -> Self {
        Self {
            settings,
            output: trace_kernel::TraceOutput::Factored,
        }
    }

    fn prepare<'a>(self, expr: AtomView<'a>) -> AtomOrView<'a> {
        let mut result = AtomOrView::View(expr);
        if self.settings.conjugate {
            result = Self::conjugate_matrices::<AbstractIndex>(result.as_view()).into();
        }
        if self.settings.gamma0 {
            result = Self::factor_gamma_zero(result.as_view()).into();
        }
        result
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

    /// Recognize one trace multiplied by the terminal evaluator's exact scalar spectators.
    fn scalar_trace_factor(expr: AtomView<'_>, trace: u32) -> Option<AtomView<'_>> {
        let AtomView::Mul(product) = expr else {
            return None;
        };
        let trace = product.iter().find(
            |factor| matches!(factor, AtomView::Fun(function) if function.get_symbol_id() == trace),
        )?;
        product
            .iter()
            .all(|factor| {
                if factor == trace {
                    return true;
                }
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
            .then_some(trace)
    }

    /// Borrow the admitted trace's cyclic wrapper and Minkowski arguments.
    /// The tensor owner still checks callbacks, compact leaves and its declared interface.
    pub(crate) fn terminal_trace_interface_inputs(
        expr: AtomView<'_>,
    ) -> Option<(Option<AtomView<'_>>, Vec<AtomView<'_>>)> {
        // Retain the existing rejection of nested trace metadata even when the
        // requested algebra operation leaves its result factored.
        let settings = GammaSimplifySettings::default();
        let simplifier = DiracSimplifier::new(&settings);
        let trace = match expr {
            AtomView::Fun(_) => expr,
            // Exact spectators cannot change tensor ports. Their arithmetic
            // emission still follows the evaluator's stricter contextual rules.
            AtomView::Mul(_) => Self::scalar_trace_factor(expr, T.trace.get_id())?,
            _ => return None,
        };
        let word = simplifier.parse_terminal_trace::<false>(trace)?;
        // General trace evaluation also admits formal spin dimensions. That
        // admission does not establish this stronger source-only interface proof.
        if !has_four_dimensional_trace_rep(word.representation) {
            return None;
        }
        let AtomView::Fun(function) = trace else {
            unreachable!("terminal traces are functions");
        };
        let cyclic = if function.get_nargs() == 2 {
            function.iter().nth(1).filter(|argument| {
                matches!(argument, AtomView::Fun(wrapper) if wrapper.get_symbol() == *shadowing::CYCLIC)
            })
        } else {
            None
        };
        Some((cyclic, word.indices))
    }

    fn parse_terminal_trace<'a, const SCALAR_CONTEXT: bool>(
        self,
        expr: AtomView<'a>,
    ) -> Option<TerminalTrace<'a>> {
        let AtomView::Fun(f) = expr else {
            return None;
        };
        let (rep, factors) = shadowing::trace_parts(f)?;
        if !SCALAR_CONTEXT {
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
        Some(TerminalTrace {
            representation: rep,
            indices,
            axial,
            repeated,
            four_dimensional: has_four_dimensional_trace_rep(rep)
                && factors
                    .iter()
                    .all(|factor| factor.gamma_dimension().is_some_and(is_four_dimension)),
        })
    }

    /// Evaluate standalone traces after reducing summed indices within the
    /// word. The surviving explicit indices occur once per monomial and compact
    /// arguments produce scalar dots, so no outer tensor contraction remains.
    /// A scalar-product context additionally requires callback-free words
    /// with exact leaf metadata. Free slots occur once per term; wildcard
    /// metadata may retain inert closed diagonals. Exact scalar spectators cannot
    /// introduce further tensor cleanup.
    pub(crate) fn evaluate_terminal_trace<const SCALAR_CONTEXT: bool>(
        self,
        expr: AtomView<'_>,
    ) -> Option<Atom> {
        let TerminalTrace {
            representation: rep,
            indices,
            axial,
            repeated,
            four_dimensional,
        } = self.parse_terminal_trace::<SCALAR_CONTEXT>(expr)?;
        let output = self.output;
        if indices.len() % 2 == 1 || (axial && indices.len() < 4) {
            return Some(Atom::Zero);
        }
        let compact = !axial && indices.iter().all(|&index| !is_minkowski_slot(index));
        let result = if !compact && repeated == 0 && (axial || four_dimensional) {
            trace_kernel::evaluate(&indices, axial, output)?
        } else {
            let terminal = trace!(rep; std::iter::empty::<Atom>());
            let trace_unit = Self::simplify_trace_terminal(terminal.as_view())?;
            trace_kernel::evaluate_generic(
                &indices,
                trace_unit.as_view(),
                !(SCALAR_CONTEXT
                    && repeated == 0
                    && indices.iter().all(|&index| is_minkowski_slot(index))),
                output,
            )
        };
        Some(result)
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
        let metric_term = g!(mu, nu) * word.build(factors);
        SchoonschipWithSettings {
            settings: &SchoonschipSettings::default().with_chain_like_functions(),
        }
        .run(metric_term.as_view(), &mut Vec::new())
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

        // Tr(Aᵀ ... Zᵀ) = Tr(Z ... A). Normalize only uniformly reversed
        // gamma/gamma5 words here; mixed words and open chains retain their
        // existing opaque treatment. Normalize the local factor
        // views, so surrounding scalar callbacks see the evaluated result in
        // this same pass rather than an intermediate forward trace.
        let transposed =
            (matches!(rep, AtomView::Fun(spin) if spin.get_symbol() == *BISPINOR_SYMBOL)
                && factors.iter().all(|factor| {
                    matches!(factor, AtomView::Fun(word)
                    if ((word.get_symbol() == AGS.gamma && word.get_nargs() == 3)
                        || (word.get_symbol() == AGS.gamma5 && word.get_nargs() == 2))
                        && has_forward_chain_endpoints(
                            word.iter().nth(1).unwrap(), word.iter().next().unwrap()))
                }))
            .then(|| {
                factors
                    .iter()
                    .rev()
                    .map(|factor| {
                        let AtomView::Fun(word) = factor else {
                            unreachable!()
                        };
                        if word.get_symbol() == AGS.gamma {
                            gamma_factor(word.iter().nth(2).unwrap())
                        } else {
                            gamma5_factor()
                        }
                    })
                    .collect::<Vec<_>>()
            });
        let factors = transposed
            .as_ref()
            .map(|word| word.iter().map(Atom::as_view).collect())
            .unwrap_or(factors);

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
            && let Some(rewritten) = self.simplify_gamma5_trace_node(rep, &factors)
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
            && let Some(result) = trace_kernel::evaluate(&trace_mink_indices, false, self.output)
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
                self.output,
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

    fn simplify_gamma5_trace_node(
        self,
        rep: AtomView<'_>,
        factors: &[DiracFactor<'_>],
    ) -> Option<Atom> {
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
            1 => self.simplify_single_gamma5_trace(rep, factors, gamma5_positions[0]),
            _ => Self::reduce_gamma5_trace_pair(rep, factors, gamma5_positions[0]),
        }
    }

    fn simplify_single_gamma5_trace(
        self,
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
        if let Some(result) = trace_kernel::evaluate(&mink_indices, true, self.output) {
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

    mod contracted_trace;

    #[test]
    fn uniformly_reversed_gamma_traces_evaluate_in_the_same_pass() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let reverse_gamma = |argument: &Atom| {
            FunctionBuilder::new(AGS.gamma)
                .add_arg(Atom::var(T.chain_out))
                .add_arg(Atom::var(T.chain_in))
                .add_arg(argument)
                .finish()
        };
        for representation in [&reps.mink4, &reps.mink_d] {
            let slots = [97101, 97102, 97103, 97104]
                .map(|index| representation.to_symbolic([Atom::num(index)]));
            let compact = representation.to_symbolic([]);
            let p = momenta(&compact);
            for arguments in [
                slots.to_vec(),
                p[..4].to_vec(),
                vec![
                    slots[0].clone(),
                    p[0].clone(),
                    slots[1].clone(),
                    p[1].clone(),
                ],
            ] {
                for length in [2, 3, 4] {
                    let arguments = &arguments[..length];
                    let input = trace!(&spin; arguments.iter().map(reverse_gamma));
                    let forward = trace!(&spin; arguments.iter().rev().map(gamma_factor));
                    let settings = GammaSimplifySettings::default();
                    let expected =
                        settings.rewrite_expression(forward, trace_kernel::TraceOutput::Factored);
                    assert_eq!(
                        settings
                            .rewrite_expression(input.clone(), trace_kernel::TraceOutput::Factored),
                        expected
                    );
                    if length == 3 {
                        assert!(expected.is_zero());
                    }
                    let source = SymbolicTensor::<PartialStructure>::infer(input).unwrap();
                    let result = source.simplify_gamma(settings).unwrap();
                    assert_eq!(result.root().structure, source.structure);
                    assert_eq!(result.resolved().unwrap().into_expression(), expected);
                    assert_eq!(
                        result
                            .simplify_gamma(settings)
                            .unwrap()
                            .resolved()
                            .unwrap()
                            .into_expression(),
                        expected
                    );
                }
            }
        }
        let a = reps.mink4.to_symbolic([Atom::num(97105)]);
        let b = reps.mink4.to_symbolic([Atom::num(97106)]);
        let input = trace!(&spin; [&a, &b].map(reverse_gamma));
        assert_eq!(
            GammaSimplifySettings::default()
                .rewrite_expression(input, trace_kernel::TraceOutput::Factored),
            Atom::num(4) * g!(&a, &b)
        );
    }

    #[test]
    fn reversed_gamma_trace_normalization_keeps_mixed_words_and_open_chains_opaque() {
        let reps = test_initialize();
        let a = reps.mink4.to_symbolic([Atom::num(97201)]);
        let b = reps.mink4.to_symbolic([Atom::num(97202)]);
        let reverse_gamma = |argument: &Atom| {
            FunctionBuilder::new(AGS.gamma)
                .add_arg(Atom::var(T.chain_out))
                .add_arg(Atom::var(T.chain_in))
                .add_arg(argument)
                .finish()
        };
        let mixed = trace!(reps.bis4.to_symbolic([]); [reverse_gamma(&a), gamma_factor(&b)]);
        let start = reps.bis4.to_symbolic([Atom::num(97203)]);
        let end = reps.bis4.to_symbolic([Atom::num(97204)]);
        let open = chain!(&start, &end; [&a, &b].map(reverse_gamma));
        for input in [mixed, open] {
            assert_eq!(
                GammaSimplifySettings::default()
                    .rewrite_expression(input.clone(), trace_kernel::TraceOutput::Factored),
                input
            );
        }
    }

    #[test]
    fn reversed_gamma_trace_callbacks_observe_only_the_same_pass_result() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let [a, b, c, d] =
            [97301, 97302, 97303, 97304].map(|index| reps.mink4.to_symbolic([Atom::num(index)]));
        let reversed = trace!(&spin; [&a, &b].map(|argument| {
            FunctionBuilder::new(AGS.gamma)
                .add_arg(Atom::var(T.chain_out))
                .add_arg(Atom::var(T.chain_in))
                .add_arg(argument)
                .finish()
        }));
        let forward = trace!(&spin; [&b, &a].map(gamma_factor));
        let calls = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let recorded = std::sync::Arc::clone(&calls);
        let callback = symbolica::symbol!("idenso::reversed_trace_scalar_callback"; Scalar;
            norm = move |node, output| {
                recorded.lock().unwrap().push(node.to_owned());
                if !node.contains_symbol(T.trace) { **output = Atom::num(7); }
            }
        );
        let mut results = Vec::new();
        for inner in [forward, reversed] {
            let input =
                symbolica::function!(callback, inner) * trace!(&spin; [&c, &d].map(gamma_factor));
            let source = SymbolicTensor::<PartialStructure>::infer(input).unwrap();
            calls.lock().unwrap().clear();
            let output = source
                .simplify_gamma(GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
            results.push((output, std::mem::take(&mut *calls.lock().unwrap())));
        }
        assert_eq!(results[0], results[1]);
        assert_eq!(results[0].0, Atom::num(28) * g!(&c, &d));
        assert_eq!(results[0].1.len(), 1);
        assert!(!results[0].1[0].contains_symbol(T.trace));
    }

    #[test]
    fn reversed_gamma_trace_keeps_component_callback_results_and_order() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let calls = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let recorded = std::sync::Arc::clone(&calls);
        let vector = spenso::vector_symbol!(
            "idenso::reversed_trace_component_callback",
            norm = move |value, output| {
                if let AtomView::Fun(vector) = value
                    && let Some(AtomView::Fun(slot)) = vector.iter().last()
                    && slot.get_symbol() == *MINKOWSKI_SYMBOL
                    && slot.get_nargs() == 2
                {
                    recorded.lock().unwrap().push(value.to_owned());
                    **output = Atom::num(1);
                }
            }
        );
        let compact = reps.mink_d.to_symbolic([]);
        let p = symbolica::function!(vector, Atom::num(1), &compact);
        let q = symbolica::function!(vector, Atom::num(2), &compact);
        let a = reps.mink_d.to_symbolic([Atom::num(97401)]);
        let arguments = [&a, &p, &a, &q];
        let forward = trace!(&spin; arguments.iter().rev().map(|argument| gamma_factor(*argument)));
        let reversed = trace!(&spin; arguments.map(|argument| {
            FunctionBuilder::new(AGS.gamma)
                .add_arg(Atom::var(T.chain_out))
                .add_arg(Atom::var(T.chain_in))
                .add_arg(argument)
                .finish()
        }));
        let spectator = symbolica::parse_lit!((reversed_trace_x + reversed_trace_y) ^ 8);
        let mut results = Vec::new();
        for input in [forward, reversed] {
            let source = SymbolicTensor::<PartialStructure>::infer(&spectator * input).unwrap();
            calls.lock().unwrap().clear();
            let result = source
                .simplify_gamma(GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
            results.push((result, std::mem::take(&mut *calls.lock().unwrap())));
        }
        assert_eq!(results[0], results[1]);
        assert!(!results[0].1.is_empty());
        let dimension = mink_slot_dimension(a.as_view()).unwrap();
        assert_eq!(
            results[0].0,
            spectator * (Atom::num(8) - Atom::num(4) * dimension * g!(&p, &q))
        );
    }

    #[test]
    fn uniformly_reversed_axial_traces_preserve_epsilon_and_pair_signs() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let slots =
            [97501, 97502, 97503, 97504].map(|index| reps.mink4.to_symbolic([Atom::num(index)]));
        let reverse = |factor: &Atom| {
            let AtomView::Fun(word) = factor.as_view() else {
                unreachable!()
            };
            let mut result = FunctionBuilder::new(word.get_symbol())
                .add_arg(Atom::var(T.chain_out))
                .add_arg(Atom::var(T.chain_in));
            for argument in word.iter().skip(2) {
                result = result.add_arg(argument);
            }
            result.finish()
        };
        let settings = GammaSimplifySettings::default();
        let compact = momenta(&reps.mink4.to_symbolic([]));
        for arguments in [slots.to_vec(), compact[..4].to_vec()] {
            for length in [2, 3, 4] {
                for position in 0..=length {
                    for paired in [false, true] {
                        let mut word = arguments[..length]
                            .iter()
                            .map(gamma_factor)
                            .collect::<Vec<_>>();
                        word.insert(position, gamma5_factor());
                        if paired {
                            word.push(gamma5_factor());
                        }
                        let input = trace!(&spin; word.iter().map(reverse));
                        let forward = trace!(&spin; word.iter().rev());
                        assert_eq!(
                            settings.rewrite_expression(
                                input.clone(),
                                trace_kernel::TraceOutput::Factored
                            ),
                            settings.rewrite_expression(
                                forward.clone(),
                                trace_kernel::TraceOutput::Factored
                            )
                        );
                        let expected = SymbolicTensor::<PartialStructure>::infer(forward)
                            .unwrap()
                            .simplify_gamma(settings)
                            .unwrap()
                            .resolved()
                            .unwrap()
                            .into_expression();
                        let source = SymbolicTensor::<PartialStructure>::infer(input).unwrap();
                        let output = source.simplify_gamma(settings).unwrap();
                        assert_eq!(output.root().structure, source.structure);
                        assert_eq!(output.resolved().unwrap().into_expression(), expected);
                        assert_eq!(
                            output
                                .simplify_gamma(settings)
                                .unwrap()
                                .resolved()
                                .unwrap()
                                .into_expression(),
                            expected
                        );
                        if length == 3 || (!paired && length == 2) {
                            assert!(expected.is_zero());
                        }
                        if !paired && length == 4 && arguments == slots {
                            let epsilon = FunctionBuilder::new(*crate::epsilon::EPSILON_SYMBOL)
                                .add_arg(&slots[0])
                                .add_arg(&slots[1])
                                .add_arg(&slots[2])
                                .add_arg(&slots[3])
                                .finish();
                            assert_eq!(
                                expected,
                                Atom::num(if position % 2 == 0 { 4 } else { -4 }) * epsilon
                            );
                        }
                        if paired && length == 2 {
                            assert_eq!(
                                expected,
                                Atom::num(if position % 2 == 0 { 4 } else { -4 })
                                    * g!(&arguments[0], &arguments[1])
                            );
                        }
                    }
                }
            }
        }
        // The existing gamma5 owner is four-dimensional, even after orientation normalization.
        let d = reps.mink_d.to_symbolic([Atom::num(97510)]);
        let generic = trace!(&spin; [gamma5_factor(), gamma_factor(&d)].iter().map(reverse));
        assert_eq!(
            settings.rewrite_expression(generic.clone(), trace_kernel::TraceOutput::Factored),
            generic
        );
        let mixed = trace!(&spin; [reverse(&gamma5_factor()), gamma_factor(&slots[0])]);
        let open = chain!(reps.bis4.to_symbolic([Atom::num(97511)]), reps.bis4.to_symbolic([Atom::num(97512)]); [gamma5_factor(), gamma_factor(&slots[0])].iter().map(reverse));
        for input in [mixed, open] {
            assert_eq!(
                settings.rewrite_expression(input.clone(), trace_kernel::TraceOutput::Factored),
                input
            );
        }
    }

    #[test]
    fn reversed_axial_trace_callbacks_observe_only_the_evaluated_result() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let arguments =
            [97601, 97602, 97603].map(|index| reps.mink4.to_symbolic([Atom::num(index)]));
        let reversed_factors = arguments
            .iter()
            .map(|argument| {
                FunctionBuilder::new(AGS.gamma)
                    .add_arg(Atom::var(T.chain_out))
                    .add_arg(Atom::var(T.chain_in))
                    .add_arg(argument)
                    .finish()
            })
            .chain(std::iter::once(
                FunctionBuilder::new(AGS.gamma5)
                    .add_arg(Atom::var(T.chain_out))
                    .add_arg(Atom::var(T.chain_in))
                    .finish(),
            ))
            .collect::<Vec<_>>();
        let reverse = trace!(&spin; &reversed_factors);
        let forward = trace!(&spin; std::iter::once(gamma5_factor()).chain(arguments.iter().rev().map(gamma_factor)));
        let calls = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let recorded = std::sync::Arc::clone(&calls);
        let callback = symbolica::symbol!("idenso::reversed_axial_trace_scalar_callback"; Scalar;
            norm = move |node, output| {
                recorded.lock().unwrap().push(node.to_owned());
                if !node.contains_symbol(T.trace) { **output = Atom::num(7); }
            }
        );
        let settings = GammaSimplifySettings::default();
        let mut results = Vec::new();
        for inner in [forward, reverse] {
            let input = symbolica::function!(callback, inner);
            calls.lock().unwrap().clear();
            let output = settings.rewrite_expression(input, trace_kernel::TraceOutput::Factored);
            results.push((output, std::mem::take(&mut *calls.lock().unwrap())));
        }
        assert_eq!(results[0], results[1]);
        assert_eq!(results[0].0, Atom::num(7));
        assert_eq!(results[0].1.len(), 1);
        assert!(!results[0].1[0].contains_symbol(T.trace));
    }

    #[test]
    fn typed_trace_definitions_keep_free_ports_and_factored_spectators() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let spectator = symbolica::parse_lit!((typed_trace_x + 1) * (typed_trace_y + 1));
        for (representation, length) in [(&reps.mink_d, 6), (&reps.mink4, 12)] {
            let slots = (0..length)
                .map(|index| representation.to_symbolic([Atom::num(93100 + index)]))
                .collect::<Vec<_>>();
            let input = &spectator * trace!(&spin; slots.iter().map(|slot| gamma!(slot)));
            let source = SymbolicTensor::<PartialStructure>::infer(input.clone()).unwrap();
            let result = source
                .simplify_gamma(GammaSimplifySettings::default())
                .unwrap();
            assert_eq!(result.root().structure, source.structure);
            assert!(
                result
                    .aliases()
                    .unwrap()
                    .iter()
                    .any(|(handle, _)| !handle.is_scalar())
            );
            let root = result.root();
            let AtomView::Mul(product) = root.expression.as_view() else {
                panic!("trace result must retain the surrounding product");
            };
            let AtomView::Mul(spectators) = spectator.as_view() else {
                unreachable!();
            };
            for spectator in spectators.iter() {
                assert!(product.iter().any(|factor| factor == spectator));
            }
            assert!(
                (result.expanded().unwrap().expression
                    - GammaSimplifySettings::default()
                        .rewrite_expression(input.clone(), trace_kernel::TraceOutput::Expanded))
                .expand()
                .is_zero()
            );
        }
    }

    #[test]
    fn typed_odd_trace_retains_its_zero_interface() {
        let reps = test_initialize();
        let input = trace!(reps.bis4.to_symbolic([]); (0..3).map(|index|
            gamma!(reps.mink_d.to_symbolic([Atom::num(93200 + index)]))));
        let source = SymbolicTensor::<PartialStructure>::infer(input).unwrap();
        let result = source
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap();
        assert!(result.root().expression.is_zero());
        assert_eq!(result.root().structure, source.structure);
    }

    #[test]
    fn production_odd_trace_closes_before_clifford_reduction() {
        test_initialize();
        T.rank_one_tensor_symbol("gammalooprs::Q");
        // One existing top-level summand of the captured GL16 aa -> aa
        // numerator: nine gamma matrices form one closed spinor cycle.
        let input = Atom::parse(
            include_str!("../../tests/fixtures/aa_aa_2l_gl16_odd_trace_term.sym"),
            "spenso",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let spectator = symbolica::parse_lit!((gamma_completion_x + gamma_completion_y) ^ 30);
        let source = SymbolicTensor::infer(&spectator * input).unwrap();
        let settings = GammaSimplifySettings::default();
        let chains = source
            .simplify_gamma(GammaSimplifySettings {
                output: GammaOutput::Chains,
                ..settings
            })
            .unwrap();
        assert!(
            matches!(chains.root().expression.as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == spectator.as_view()))
        );
        let expected = chains.simplify_gamma(settings).unwrap();
        assert!(expected.resolved().unwrap().expression.is_zero());
        let result = source.simplify_gamma(settings).unwrap();
        assert!(result.resolved().unwrap().expression.is_zero());
        assert_eq!(result.root().structure, source.structure);
        assert_eq!(result.resolved().unwrap(), expected.resolved().unwrap());
        let completed = result
            .simplify_color(crate::color::ColorSimplifySettings::default())
            .unwrap()
            .contract(Default::default())
            .unwrap();
        assert!(completed.resolved().unwrap().expression.is_zero());
        assert_eq!(completed.root().structure, source.structure);
    }

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
        let expected = crate::tensor::SymbolicTensor::infer(
            (&spectator
                * crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression())
            .as_atom_view()
            .to_owned(),
        )
        .unwrap()
        .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
        let output = crate::tensor::SymbolicTensor::infer((decorated).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        assert_eq!(output, expected);
        assert!(matches!(output.as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == spectator.as_view())));
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((output).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            output
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((decorated).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(GammaSimplifySettings::default().without_trace_evaluation())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
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
            let settings = GammaSimplifySettings::default();
            let standalone = DiracSimplifier::new(&settings)
                .evaluate_terminal_trace::<false>(input.as_view())
                .unwrap();
            assert_eq!(admitted, standalone, "{name}: terminal result");
            let raw = &spectator * &standalone;
            assert_eq!(
                settings.rewrite_expression(raw.clone(), trace_kernel::TraceOutput::Factored),
                raw,
                "{name}: local cleanup"
            );
            let decorated = &spectator * &input;
            let trace = DiracSimplifier::scalar_trace_factor(decorated.as_view(), T.trace.get_id())
                .expect("exact scalar spectators retain the terminal word");
            let contextual = &spectator
                * DiracSimplifier::new(&settings)
                    .evaluate_terminal_trace::<true>(trace)
                    .unwrap();
            assert_eq!(contextual, raw, "{name}: contextual result");
            assert_eq!(
                settings
                    .rewrite_expression(contextual.clone(), trace_kernel::TraceOutput::Factored),
                contextual,
                "{name}: rerun"
            );
            // Bare tagged tensor heads are raw scalar variables in this kernel
            // fixture, but deliberately fail the shared typed admission boundary.
            assert!(SymbolicTensor::<PartialStructure>::infer(decorated).is_err());
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
        let output =
            crate::tensor::SymbolicTensor::infer((&spectator * &input).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
        assert_eq!(output, expected);
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((output).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            output
        );
        assert!(
            ((&output / &spectator)
                - crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression())
            .expand()
                != Atom::Zero
        );
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
        // Rounded representation metadata belongs to this local algebra
        // regression, not typed admission. Compare the same trace-node pass
        // with and without an inert spectator, then compare the next pass too.
        let settings = GammaSimplifySettings::default();
        let rewrite = |expression| {
            settings.rewrite_expression(expression, trace_kernel::TraceOutput::Factored)
        };
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
            // Compare exact floating coefficients after the same local pass,
            // with and without a surrounding inert function, without a tolerance.
            let reference = rewrite(&inert * &input) / &inert;
            let output = rewrite(input);
            assert_eq!(output, reference, "{name}: first call");
            assert_eq!(rewrite(output), rewrite(reference), "{name}: rerun");
        }
    }

    #[test]
    fn scalar_callback_observes_factored_trace_results_in_the_shared_pipeline() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let [a, b, c, d] =
            [95101, 95102, 95103, 95104].map(|index| reps.mink4.to_symbolic([Atom::num(index)]));
        let inner = trace!(&spin, gamma!(&a), gamma!(&b));
        let outer = trace!(&spin, gamma!(&c), gamma!(&d));
        let calls = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let recorded = std::sync::Arc::clone(&calls);
        let callback = symbolica::symbol!("shared_gamma_scalar_trace_callback"; Scalar;
            norm = move |node, output| {
                recorded.lock().unwrap().push(node.to_owned());
                if !node.contains_symbol(T.trace) {
                    **output = Atom::num(7);
                }
            }
        );
        let input = symbolica::function!(callback, &inner) * &outer;
        let source = SymbolicTensor::<PartialStructure>::infer(input.clone()).unwrap();
        let settings = GammaSimplifySettings::default();
        calls.lock().unwrap().clear();
        // This is the original local kernel output mode, not a second call
        // through the new aliased scheduler.
        let expected = settings.rewrite_expression(input, trace_kernel::TraceOutput::Factored);
        let expected_calls = std::mem::take(&mut *calls.lock().unwrap());
        assert_eq!(expected, Atom::num(28) * g!(&c, &d));
        assert_eq!(expected_calls.len(), 1);
        let AtomView::Fun(observed) = expected_calls[0].as_view() else {
            panic!()
        };
        assert_eq!(
            observed.iter().next().unwrap(),
            (Atom::num(4) * g!(&a, &b)).as_view()
        );
        let actual = source
            .simplify_gamma(settings)
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        assert_eq!(actual, expected);
        assert_eq!(*calls.lock().unwrap(), expected_calls);
    }

    #[test]
    fn contextual_callback_kernel_reports_nonstabilization_at_shared_default_bound() {
        use std::sync::{
            Arc, OnceLock,
            atomic::{AtomicUsize, Ordering},
        };
        let r = test_initialize();
        let word = trace!(
            r.bis4.to_symbolic([]),
            gamma!(slot!(r.mink4, 95301)),
            gamma!(slot!(r.mink4, 95301))
        );
        let heads = Arc::new(OnceLock::<[Symbol; 2]>::new());
        let calls = Arc::new(AtomicUsize::new(0));
        let mut symbols = Vec::new();
        for (position, name) in ["gamma_callback_cycle_left", "gamma_callback_cycle_right"]
            .into_iter()
            .enumerate()
        {
            let heads = Arc::clone(&heads);
            let calls = Arc::clone(&calls);
            let word = word.clone();
            symbols.push(
                symbolica::symbol!(name; Scalar; norm = move |node, output| {
                    if !node.contains_symbol(T.trace) {
                        calls.fetch_add(1, Ordering::Relaxed);
                        **output = symbolica::function!(heads.get().unwrap()[1 - position], &word);
                    }
                }),
            );
        }
        heads.set([symbols[0], symbols[1]]).unwrap();
        let source =
            SymbolicTensor::<PartialStructure>::infer(symbolica::function!(symbols[0], word))
                .unwrap();
        calls.store(0, Ordering::Relaxed);
        let error = source
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("callback-sensitive gamma kernel did not stabilize")
        );
        assert_eq!(
            calls.load(Ordering::Relaxed),
            crate::tensor::simplification::SimplifySettings::default().max_passes
        );
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
                crate::tensor::SymbolicTensor::infer((inner).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                symbolica::function!(other, &compact)
            );
        // Retain the original raw primitive oracle: one gamma rewrite emits
        // the nested vector, dot cleanup invokes its compact callback, and a
        // second gamma rewrite evaluates the newly introduced trace.
        let settings = GammaSimplifySettings::default();
        let output_kind = trace_kernel::TraceOutput::Factored;
        let rewritten = settings.rewrite_expression(input.clone(), output_kind);
        let cleaned = rewritten.normalize_dots();
        assert!(cleaned.contains_symbol(T.trace));
        let output = settings.rewrite_expression(cleaned, output_kind);
        assert!(!output.contains_symbol(T.trace));
        assert_eq!(output, expected);
        assert_eq!(
            settings.rewrite_expression(output.normalize_dots(), output_kind),
            output
        );

        // This raw callback fixture loses mu and introduces rho and sigma in
        // its place. Its final metric operand is a rank-two metric expression,
        // not a slot or compact vector. The typed owner must reject that change.
        use spenso::structure::partial::PartialStructureExt;
        let source = crate::tensor::SymbolicTensor::infer(input).unwrap();
        assert_eq!(source.structure.logical_slots().len(), 1);
        let inner = crate::tensor::SymbolicTensor::infer(inner).unwrap();
        assert_eq!(inner.structure.logical_slots().len(), 2);
        assert!(crate::tensor::SymbolicTensor::infer(output).is_err());
        assert!(source.simplify_gamma(settings).is_err());
    }

    #[test]
    fn rewrite_guard_retains_empty_trace_and_non_gamma_chain_identities() {
        let r = test_initialize();
        let trace = trace!(r.bis4.to_symbolic([]); std::iter::empty::<Atom>());
        let start = slot!(r.bis4, a).into_atom();
        let end = slot!(r.bis4, b).into_atom();
        let coefficient = symbolica::parse_lit!((x + y) ^ 8);
        let settings = GammaSimplifySettings::default();
        let output = trace_kernel::TraceOutput::Factored;
        assert_eq!(
            settings.rewrite_expression(&coefficient * &trace, output),
            &coefficient * Atom::num(4)
        );
        assert_eq!(
            settings
                .without_trace_evaluation()
                .rewrite_expression(&coefficient * &trace, output),
            &coefficient * trace
        );
        for chain in [
            chain!(&start, &end),
            chain!(&start, &end, gamma5!(), gamma5!()),
        ] {
            assert_eq!(
                settings.rewrite_expression(
                    settings.rewrite_expression(&coefficient * chain, output),
                    output
                ),
                &coefficient * id_atom(start.as_view(), end.as_view())
            );
        }
        let terminal = &coefficient * g!(slot!(r.mink4, mu), slot!(r.mink4, nu));
        assert_eq!(
            settings.rewrite_expression(terminal.clone(), output),
            terminal
        );
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
                let complete = crate::tensor::SymbolicTensor::infer(fallback)
                    .unwrap()
                    .simplify_gamma(GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression();
                // Alias definitions and the direct recurrence may group the
                // same standalone Clifford polynomial differently. Keep the
                // scalar spectator factored, and compare the trace exactly.
                assert!(
                    complete.is_zero()
                        || matches!(complete.as_view(), AtomView::Mul(product)
                    if product.iter().any(|factor| factor == spectator.as_view()))
                );
                assert!(
                    (&complete / &spectator - &shortcut).expand().is_zero(),
                    "length {length}, gamma5 position {position}"
                );
                let standalone = crate::tensor::SymbolicTensor::infer(input.clone())
                    .unwrap()
                    .simplify_gamma(GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression();
                assert!((standalone - &shortcut).expand().is_zero());
                assert_eq!(
                    crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                        .unwrap()
                        .simplify_gamma(GammaSimplifySettings::default().without_trace_evaluation())
                        .unwrap()
                        .resolved()
                        .unwrap()
                        .into_expression(),
                    input
                );
            }
        }
    }

    #[test]
    fn scalar_axial_traces_match_the_full_pass_at_every_gamma5_position() {
        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let slots = (0..14)
            .map(|position| {
                reps.mink4.pattern(symbolica::symbol!(format!(
                    "contextual_axial_mu_{position}"
                )))
            })
            .collect::<Vec<_>>();
        let scalar = symbolica::parse_lit!((contextual_axial_x + contextual_axial_y) ^ 3);
        let spectators = [
            Atom::num(-3),
            Atom::num((2, 7)),
            scalar.clone(),
            -&scalar / 5,
        ];
        let inert = symbolica::function!(symbolica::symbol!("contextual_axial_inert"), 0);
        for length in 2..=14 {
            for position in 0..=length {
                let mut factors = slots[..length]
                    .iter()
                    .map(|slot| gamma!(slot))
                    .collect::<Vec<_>>();
                factors.insert(position, gamma5!());
                let trace = trace!(&spin; factors);
                {
                    let settings = GammaSimplifySettings::default();
                    let admitted = DiracSimplifier::new(&settings)
                        .evaluate_terminal_trace::<true>(trace.as_view())
                        .expect("a free axial word with symbolic slots is terminal");
                    // A function spectator forces the established full pass.
                    // Exact coefficients can then be applied without regrouping
                    // approximate arithmetic or invoking normalization callbacks.
                    let reference = crate::tensor::SymbolicTensor::infer(
                        (&inert * &trace).as_atom_view().to_owned(),
                    )
                    .unwrap()
                    .simplify_gamma(settings)
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression()
                        / &inert;
                    // Compare the same Clifford polynomial independently of
                    // how the alias path groups its metric/epsilon factors.
                    assert!(
                        (&admitted - &reference).expand().is_zero(),
                        "length {length}, gamma5 position {position}"
                    );
                    for spectator in &spectators {
                        let decorated = spectator * &trace;
                        let result = crate::tensor::SymbolicTensor::infer(
                            (decorated).as_atom_view().to_owned(),
                        )
                        .unwrap()
                        .simplify_gamma(settings)
                        .unwrap()
                        .resolved()
                        .unwrap()
                        .into_expression();
                        assert!(
                            (&result / spectator - &reference).expand().is_zero(),
                            "length {length}, gamma5 position {position}, spectator {spectator}"
                        );
                        assert_eq!(
                            crate::tensor::SymbolicTensor::infer(
                                (result).as_atom_view().to_owned()
                            )
                            .unwrap()
                            .simplify_gamma(settings)
                            .unwrap()
                            .resolved()
                            .unwrap()
                            .into_expression(),
                            result
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn scalar_axial_shortcut_retains_unsupported_and_callback_boundaries() {
        use std::sync::{Arc, Mutex};

        let reps = test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let slots = (0..16)
            .map(|position| {
                reps.mink4.pattern(symbolica::symbol!(format!(
                    "contextual_axial_guard_mu_{position}"
                )))
            })
            .collect::<Vec<_>>();
        let trace = trace!(&spin; std::iter::once(gamma5!())
            .chain(slots[..4].iter().map(|slot| gamma!(slot))));
        let unsupported = [
            trace!(&spin, gamma5!(), gamma!(&slots[0]), gamma!(&slots[0])),
            trace!(
                &spin,
                gamma5!(),
                gamma5!(),
                gamma!(&slots[0]),
                gamma!(&slots[1])
            ),
            trace!(
                &spin,
                gamma5!(),
                gamma!(spenso::p!(reps.mink4.to_symbolic([])))
            ),
            trace!(
                &spin,
                gamma5!(),
                gamma!(slot!(reps.mink_d, contextual_axial_guard_mu))
            ),
            trace!(&spin; std::iter::once(gamma5!())
                .chain(slots.iter().map(|slot| gamma!(slot)))),
        ];
        for input in unsupported {
            assert!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<true>(input.as_view())
                    .is_none(),
                "{input}"
            );
        }
        let inert = symbolica::function!(symbolica::symbol!("contextual_axial_guard_inert"), 0);
        let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
        let rounded_input = &rounded * &trace;
        assert!(
            DiracSimplifier::scalar_trace_factor(rounded_input.as_view(), T.trace.get_id())
                .is_none()
        );
        {
            let settings = GammaSimplifySettings::default();
            assert_eq!(
                crate::tensor::SymbolicTensor::infer((rounded_input).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(settings)
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                crate::tensor::SymbolicTensor::infer(
                    (&inert * &rounded_input).as_atom_view().to_owned()
                )
                .unwrap()
                .simplify_gamma(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression()
                    / &inert
            );
        }

        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = Arc::clone(&calls);
        let callback = symbolica::symbol!("contextual_axial_scalar_callback"; Scalar;
        norm = move |node, output| {
            recorded.lock().unwrap().push(node.to_owned());
            if !node.contains_symbol(T.trace) {
                **output = Atom::num(7);
            }
        });
        let spectator = symbolica::function!(
            callback,
            trace!(&spin, gamma!(&slots[4]), gamma!(&slots[5]))
        );
        let input = spectator * &trace;
        assert!(DiracSimplifier::scalar_trace_factor(input.as_view(), T.trace.get_id()).is_none());
        {
            let settings = GammaSimplifySettings::default();
            calls.lock().unwrap().clear();
            let reference =
                crate::tensor::SymbolicTensor::infer((&inert * &input).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(settings)
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression()
                    / &inert;
            let expected_calls = std::mem::take(&mut *calls.lock().unwrap());
            assert!(!expected_calls.is_empty());
            let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(settings)
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
            assert_eq!(result, reference);
            assert_eq!(*calls.lock().unwrap(), expected_calls);
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
            let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
            assert!((&result - expected).expand().is_zero());
            assert_eq!(
                DiracSimplifier::new(&GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view()),
                Some(result.clone())
            );
            assert_eq!(
                crate::tensor::SymbolicTensor::infer(
                    (&spectator * input).as_atom_view().to_owned()
                )
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
                &spectator * &result,
                "the outer scalar numerator remains factored"
            );
            assert_eq!(
                crate::tensor::SymbolicTensor::infer((result).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                result
            );
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
            assert_eq!(
                crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                input
            );
        }
        let odd = trace!(&spin;
            (0..5).map(|index| gamma!(r.mink_d.pattern(Atom::num(index))))
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((odd).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            Atom::Zero
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((odd).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(GammaSimplifySettings::default().without_trace_evaluation())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            odd
        );
        let symbolic_spin = trace!(r.bis_d.to_symbolic([]);
            (0..10).map(|index| gamma!(r.mink4.pattern(Atom::num(index))))
        );
        // The 4D kernel fixes Tr(1)=4. A symbolic spin dimension must instead
        // multiply the dimension-generic pairing formula.
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((symbolic_spin).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression()
                .expand()
                .nterms(),
            945
        );
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
        let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        let tail = trace!(&spin; slots[2..].iter().map(|slot| gamma!(slot)));
        let expected = crate::tensor::SymbolicTensor::infer((tail).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
            * mink_slot_dimension(slots[0].as_view()).unwrap();
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
        let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
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
                assert_eq!(
                    crate::tensor::SymbolicTensor::infer((paired).as_atom_view().to_owned())
                        .unwrap()
                        .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                        .unwrap()
                        .resolved()
                        .unwrap()
                        .into_expression(),
                    unit.to_owned() * &pp * qq.pow(6)
                );

                // For A=p/ q/, A²-2(p.q)A+p²q²=0. This independent scalar
                // recurrence certifies all seven alternating pairs exactly.
                let mut previous = unit.to_owned();
                let mut expected = unit * pq.as_view();
                for pairs in 1..=7 {
                    let input = trace!(&spin;
                        (0..2 * pairs).map(|i| gamma!(if i % 2 == 0 { &p } else { &q }))
                    );
                    let result =
                        crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                            .unwrap()
                            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                            .unwrap()
                            .resolved()
                            .unwrap()
                            .into_expression();
                    assert!((&result - &expected).expand().is_zero());
                    assert_eq!(
                        crate::tensor::SymbolicTensor::infer((result).as_atom_view().to_owned())
                            .unwrap()
                            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                            .unwrap()
                            .resolved()
                            .unwrap()
                            .into_expression(),
                        result
                    );
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
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            Atom::num(4) * g!(&p, &q)
        );
        let odd = trace!(&spin, gamma!(&p), gamma!(&q), gamma!(&p));
        assert!(
            crate::tensor::SymbolicTensor::infer((odd).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression()
                .is_zero()
        );

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
            crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(GammaSimplifySettings::default().without_trace_evaluation())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
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
            (crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression()
                - crate::tensor::SymbolicTensor::infer((expected_first).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression())
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
                        assert_eq!(
                            crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                                .unwrap()
                                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                                .unwrap()
                                .resolved()
                                .unwrap()
                                .into_expression(),
                            result
                        );
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
        let reduced = crate::tensor::SymbolicTensor::infer(
            (trace!(&spin; slots[1..].iter().map(|index| gamma!(index))))
                .as_atom_view()
                .to_owned(),
        )
        .unwrap()
        .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
        let expected = Atom::num(4) * unit * g!(&slots[1], &slots[2]) * g!(&slots[3], &slots[4])
            + (dimension.to_owned() - Atom::num(4)) * reduced;
        let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
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
        let complete = crate::tensor::SymbolicTensor::infer((decorated).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
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
        let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        assert!(
            (&result
                - crate::tensor::SymbolicTensor::infer((contracted).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression())
            .expand()
            .is_zero()
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((result).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            result
        );
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
            assert_eq!(
                crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(settings)
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                expected
            );

            let (start, end) = (slot!(r.bis4, i).into_atom(), slot!(r.bis4, j).into_atom());
            let input = g!(&a, &b) * chain!(&start, &end; [gamma!(&b), gamma!(&c)]);
            let expected = chain!(&start, &end; [gamma!(&a), gamma!(&c)]);
            assert_eq!(
                crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                expected
            );
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
        let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        assert_eq!(result, expected);
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((result).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            result
        );
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
            let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression();
            assert_eq!(result, expected);
            assert_eq!(
                crate::tensor::SymbolicTensor::infer((result).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression(),
                result
            );
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
