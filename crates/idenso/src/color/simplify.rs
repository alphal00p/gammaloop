use std::{cell::Cell, collections::VecDeque, sync::LazyLock};

use spenso::{
    chain,
    network::{library::symbolic::ETS, parsing::ParseState, tags::SPENSO_TAG as T},
    rep_,
    shadowing::{self, ProjectorExpander, TensorCollectFilter},
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        partial::{PartialStructure, PartialStructureExt},
        representation::{LibraryRep, RepName},
        slot::{DualSlotTo, IsAbstractSlot, ParseableAind, Slot, SlotMatch, SlotMatcher},
    },
    trace, trace_sym,
};
#[cfg(test)]
use symbolica::function;
use symbolica::{
    atom::{AddView, Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
    id::Replacement,
};
use symbolica_utils::PatternReplacement;

use crate::{
    W_, color_f, color_t,
    representations::{ColorAdjoint, ColorFundamental, ColorSextet},
    shorthands::chain::Chain,
    tensor::{SymbolicTensor, inference::TensorInferenceError},
};

use super::{CS, ColorSimplifier, ColorSimplifySettings};
use crate::tensor::simplification::{ReductionStatus, observation::SettledRegions};

mod states;
mod trace;
use trace::LineContext;

static TRACE_TERMINALS: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    [Replacement::new(
        // Empty color trace: Tr_rep(1) -> dim(rep).
        trace!(rep_!(0; W_.d_)).to_pattern(),
        Atom::var(W_.d_),
    )]
});

pub(crate) struct ColorAlgebraSimplifier {
    pub settings: ColorSimplifySettings,
    pub(crate) dummies: ParseState<AbstractIndex>,
    certificate: RowCertificate,
    /// Nesting of the term-local repeat, see [`Self::rewrite_terms_at`].
    repeat_depth: Cell<usize>,
    /// Whether a row of this call received canonical states.
    merged: Cell<bool>,
    /// Canonical forms of the states of this call, by first-occurrence
    /// labelling; see `states.rs`.
    canonical_cache: std::cell::RefCell<std::collections::HashMap<Atom, Atom>>,
}

/// Bound on the nesting of the term-local repeat. Every local rule shortens a
/// line, removes structure constants of a short cycle or a repeated pair, so
/// the bound only guards termination; reaching it leaves the rest of the
/// term to the next pass.
const TERM_REPEAT_DEPTH: usize = 64;

/// Whether the row being rewritten is left a colour fixed point. Only the
/// one-shot decomposition of a whole isolated term qualifies (see trace.rs).
#[derive(Default)]
struct RowCertificate {
    /// Set by the trace kernel after decomposing an isolated line completely.
    terminal: Cell<bool>,
    /// Cleared by any other rewrite and by a term left for a later round.
    fixed: Cell<bool>,
    /// Set by a rule that expands a term into many: a trace decomposition
    /// step (including a long adjoint cycle cut through it), a cross-line
    /// Fierz identity or a sum distribution. Its terms wait for the next pass,
    /// on the merged row, instead of the term-local repeat.
    expansion: Cell<bool>,
    /// Set by a rewrite of a whole term other than the complete
    /// decomposition of an isolated line, and by a long-cycle cut: all of the
    /// row's monomials are given canonical dummy labels so that equal states
    /// merge.
    merge: Cell<bool>,
    /// Set by a complete one-shot decomposition of a whole term, whose
    /// monomials are distinct.
    terminal_term: Cell<bool>,
}

impl RowCertificate {
    /// Apply a notation pass to a certified row; a change revokes the
    /// certificate, since the new notation may enable further rules.
    fn unless_changed(&self, expression: Atom, pass: impl FnOnce(Atom) -> Atom) -> Atom {
        if !self.fixed.get() {
            return pass(expression);
        }
        let result = pass(expression.clone());
        if result != expression {
            self.fixed.set(false);
        }
        result
    }

    /// Apply a substitution of port-free invariants to a certified row. No
    /// rule reads an invariant, so new work needs a new arithmetic context
    /// for colour factors: a collapsed sum or a changed power around them.
    /// Merged or vanishing terms leave fewer candidates and keep the row fixed.
    fn unless_colour_shape_changed(
        &self,
        expression: Atom,
        pass: impl FnOnce(Atom) -> Atom,
    ) -> Atom {
        if !self.fixed.get() {
            return pass(expression);
        }
        let shape = ColourShape::of(expression.as_view());
        let result = pass(expression);
        if ColourShape::of(result.as_view()) != shape {
            self.fixed.set(false);
        }
        result
    }
}

/// The sums and powers that contain colour ports or pending colour work.
#[derive(PartialEq)]
struct ColourShape {
    sums: usize,
    powers: Vec<Atom>,
}

impl ColourShape {
    fn of(expression: AtomView<'_>) -> Self {
        let mut shape = Self {
            sums: 0,
            powers: Vec::new(),
        };
        shape.visit(expression);
        shape
    }

    fn visit(&mut self, expression: AtomView<'_>) -> bool {
        match expression {
            AtomView::Add(sum) => {
                let colour = sum
                    .iter()
                    .fold(false, |colour, term| self.visit(term) | colour);
                self.sums += usize::from(colour);
                colour
            }
            AtomView::Mul(product) => product
                .iter()
                .fold(false, |colour, factor| self.visit(factor) | colour),
            AtomView::Pow(_) => {
                let colour = atom_contains_color_port_or_work(expression);
                if colour {
                    self.powers.push(expression.to_owned());
                }
                colour
            }
            AtomView::Fun(_) => atom_contains_color_port_or_work(expression),
            _ => false,
        }
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Also reports whether the result is certified to be a colour fixed
    /// point, so that the planner need not confirm it with a no-op round.
    ///
    /// Signed graph zeros are not pruned here: every identity is exact, so a
    /// row with an odd automorphism still reduces to zero algebraically, and
    /// the graph-automorphism pass cost more than the work it saved. Index
    /// canonicalization keeps that pass.
    pub(crate) fn simplify_color_parts(
        &self,
        settings: ColorSimplifySettings,
        settled: &mut SettledRegions,
    ) -> Result<(Self, bool), TensorInferenceError> {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ColorKernel,
        );
        let simplifier = ColorAlgebraSimplifier::new(settings, Self::reserved_dummies([self]));
        let certify = settings.one_shot_traces;
        // Foreign colour factors that no identity reaches are not colour
        // work: they stay factors, outside every row and kernel step.
        let foreign = foreign_colour_factors(self.expression.as_view());
        // Cache hits are completed, unchanged inputs and hence fixed points.
        let fixed_point = Cell::new(certify);
        let mut row = |selected: Self, complete: bool| {
            settled.run(selected, complete, |selected| {
                simplifier.certificate.fixed.set(certify);
                let expression =
                    simplifier.step_beside(selected.expression.as_view(), complete, &foreign);
                let expression = simplifier
                    .certificate
                    .unless_changed(expression, restore_explicit_su_n_generator_chains);
                let rewritten = selected.with_identity_result(expression, None)?;
                // A retained sum can include non-colour coefficients. Only
                // colour connections are prerequisites of this identity;
                // the shared contractor preserves the other representations.
                let contraction = crate::tensor::ContractSettings {
                    representations: Some(&[
                        ColorFundamental {}.into(),
                        ColorAdjoint {}.into(),
                        ColorSextet {}.into(),
                    ]),
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                };
                // The prerequisite contraction preceded this kernel and
                // the identities emit no vectors. A row without a metric
                // or identity line, or whose metrics join external ports
                // only, has nothing to contract and needs no network.
                let (mut root, status) = if rewritten.expression.contains_symbol(ETS.metric)
                    && !rewritten.only_external_metric_sources(contraction)
                {
                    let contracted = rewritten.contract_parts(contraction)?;
                    (contracted.root, contracted.status)
                } else if simplifier.certificate.fixed.get() {
                    (rewritten, ReductionStatus::Complete)
                } else {
                    (rewritten, selected.proofs.frontier)
                };
                root.proofs.frontier = status;
                // Only rows with a complete frontier form the value; the
                // others are intermediate states of the factor order. An
                // unchanged row is trivially a fixed point, and so is a row
                // left with group invariants only. Otherwise its
                // contraction acts only on certified decompositions.
                if complete && fixed_point.get() {
                    fixed_point.set(
                        status == ReductionStatus::Complete
                            && (simplifier.certificate.fixed.get()
                                || root.expression == selected.expression
                                || !atom_contains_color_port_or_work(root.expression.as_view())),
                    );
                }
                Ok(root)
            })
        };
        let value = match self.colour_term(&foreign) {
            // One product: its colour factors are a single complete row. The
            // kernel spans and distributes their colour sums itself, so the
            // shared collector's factor-by-factor frontier is not needed.
            // Scalar spectators stay outside the row.
            Some((colour, spectators)) => {
                if spectators.is_empty() {
                    row(self.clone(), true)?
                } else {
                    let root = row(self.with_identity_result(colour, None)?, true)?;
                    let mut value = self.with_identity_result(
                        Atom::mul_many(spectators.into_iter().chain([root.expression.as_view()])),
                        None,
                    )?;
                    value.proofs.frontier = root.proofs.frontier;
                    value
                }
            }
            None => self.collect_with_map(
                crate::tensor::CollectionMode::Factored,
                None,
                |value| atom_contains_color_node(value) && !foreign.contains(&value),
                |selected, complete, _| row(selected, complete),
            )?,
        };
        // The planner multiplies each row's states by the row's number.
        // Distribute it, so that equal canonical states of different rows merge.
        let value = if simplifier.merged.get()
            && let Some(merged) =
                ColorAlgebraSimplifier::merging_row_numbers(value.expression.as_view())
        {
            value.with_identity_result(merged, None)?
        } else {
            value
        };
        let certified = fixed_point.get() && value.proofs.frontier == ReductionStatus::Complete;
        Ok((value, certified))
    }
}

impl SymbolicTensor<PartialStructure> {
    /// A domain that is one product of colour factors, also sums and powers
    /// with colour, and of scalar spectators: the product of its colour
    /// factors with exact numbers, and the other spectators. `None` for a
    /// sum, without colour, or with a spectator that has an index or work
    /// of its own, such as a foreign tensor, also one of the `foreign`
    /// colour factors anywhere in its sums: the shared collector separates
    /// those, and its rows keep them beside the kernel step.
    fn colour_term(&self, foreign: &[AtomView<'_>]) -> Option<(Atom, Vec<AtomView<'_>>)> {
        let expression = self.expression.as_view();
        if !matches!(expression, AtomView::Mul(_) | AtomView::Fun(_)) || !foreign.is_empty() {
            return None;
        }
        let observed = self.reduction_observations();
        let mut colour = Vec::new();
        let mut spectators = Vec::new();
        for factor in multiplicative_factor_views(expression) {
            if atom_contains_color_node(factor)
                || matches!(factor, AtomView::Num(number)
                    if matches!(number.get_coeff_view(), CoefficientView::Natural(..) | CoefficientView::Large(..)))
            {
                colour.push(factor);
            } else if observed.certify_scalar_region(factor).is_some() {
                spectators.push(factor);
            } else {
                return None;
            }
        }
        if !colour
            .iter()
            .any(|factor| atom_contains_color_node(*factor))
        {
            return None;
        }
        let colour = if spectators.is_empty() {
            self.expression.clone()
        } else {
            Atom::mul_many(colour)
        };
        Some((colour, spectators))
    }
}

impl ColorAlgebraSimplifier {
    pub(crate) fn new(settings: ColorSimplifySettings, dummies: ParseState<AbstractIndex>) -> Self {
        Self {
            settings,
            dummies,
            certificate: RowCertificate::default(),
            repeat_depth: Cell::new(0),
            merged: Cell::new(false),
            canonical_cache: Default::default(),
        }
    }

    fn color_adjoint_dummy_like(&self, slot: &Atom) -> Option<Atom> {
        let dimension = color_adjoint_dimension(slot)?;
        Some(ColorAdjoint {}.to_symbolic([dimension, self.dummies.fresh_index().to_atom()]))
    }

    fn color_adjoint_dummy_for_pair(&self, left: &Atom, right: &Atom) -> Option<Atom> {
        let left_dimension = color_adjoint_dimension(left)?;
        let right_dimension = color_adjoint_dimension(right)?;
        (left_dimension == right_dimension).then(|| {
            ColorAdjoint {}.to_symbolic([left_dimension, self.dummies.fresh_index().to_atom()])
        })
    }

    /// Apply the local colour identities. Contraction and cross-domain
    /// scheduling belong to the shared tensor owner.
    pub(crate) fn step(&self, expression: AtomView<'_>, complete: bool) -> Atom {
        self.step_within(expression, complete, &[])
    }

    /// [`Self::step`] on a row beside other factors: the canonical
    /// relabelling of its states avoids their labels.
    fn step_within(
        &self,
        expression: AtomView<'_>,
        complete: bool,
        beside: &[AtomView<'_>],
    ) -> Atom {
        // Without an indexed fundamental slot there is no generator to join:
        // a closed trace names only its representation.
        let collected = if !has_fundamental_slot(expression) {
            expression.to_owned()
        } else if let Some(joined) = joined_fundamental_lines(expression) {
            joined
        } else {
            expression.join_chains(ColorFundamental {}.into())
        };
        self.certificate.merge.set(false);
        let rewritten = self.rewrite_terms(collected.as_view(), complete);
        // A rewritten row with a complete frontier is a whole top-level term:
        // its generated states can take canonical dummy labels and merge. A
        // complete decomposition of an isolated line generates distinct
        // states; relabel only a few of them, so small results stay canonical.
        let limit = if self.certificate.merge.take() {
            states::STATE_LIMIT
        } else {
            states::DISTINCT_STATE_LIMIT
        };
        let rewritten = if complete
            && rewritten != collected
            && let Some(canonical) = self.canonical_states(rewritten.as_view(), limit, beside)
        {
            self.merged.set(true);
            canonical
        } else {
            rewritten
        };
        // Only invariants are substituted; most rows of a round carry none.
        if self.settings.substitute_cof_dimension_invariants
            && [CS.cas, CS.idx, CS.gram]
                .iter()
                .any(|&invariant| rewritten.contains_symbol(invariant))
        {
            self.certificate
                .unless_colour_shape_changed(rewritten, |rewritten| {
                    rewritten.to_cof_dimension_invariants()
                })
        } else {
            rewritten
        }
    }

    /// [`Self::step`] on the colour factors of a row; its `foreign` factors,
    /// which no identity reaches, stay factors of the result as found. A
    /// term of a sum is a row with them. The canonical relabelling avoids
    /// the labels of every foreign factor, in this row or outside it.
    fn step_beside(
        &self,
        expression: AtomView<'_>,
        complete: bool,
        foreign: &[AtomView<'_>],
    ) -> Atom {
        if foreign.is_empty() {
            return self.step(expression, complete);
        }
        let (kept, colour): (Vec<_>, Vec<_>) = multiplicative_factor_views(expression)
            .into_iter()
            .partition(|factor| foreign.contains(factor));
        if kept.is_empty() {
            return self.step_within(expression, complete, foreign);
        }
        let stepped = self.step_within(Atom::mul_many(colour).as_view(), complete, foreign);
        Atom::mul_many(kept.into_iter().chain([stepped.as_view()]))
    }

    fn rewrite_terms(&self, expr: AtomView<'_>, complete: bool) -> Atom {
        self.rewrite_terms_at(expr, complete, true)
    }

    /// `top` marks whole top-level terms of the row, as in [`Self::rewrite_node`].
    fn rewrite_terms_at(&self, expr: AtomView<'_>, complete: bool, top: bool) -> Atom {
        // Terminal trace rules can create sums; product rules such as f*f -> CA*g
        // then need to run on each generated term instead of on the whole Add.
        if let AtomView::Add(add) = expr {
            let terms = add
                .iter()
                .map(|term| self.rewrite_terms_at(term, complete, top))
                .collect::<Vec<_>>();
            if !expr.needs_normalization()
                && terms
                    .iter()
                    .zip(add.iter())
                    .all(|(new, old)| new.as_view() == old)
            {
                return expr.to_owned();
            }
            return Self::sum_rewritten_terms(terms);
        }

        let outer = self.certificate.expansion.take();
        let rewritten = self.rewrite_node(expr, complete, top);
        let expanded = self.certificate.expansion.replace(outer);
        let terminal_term = self.certificate.terminal_term.take();
        if let Some(rewritten) = rewritten {
            if top && complete && !terminal_term {
                self.certificate.merge.set(true);
            }
            if !(top && complete) || expanded {
                return rewritten;
            }
            // FORM's repeat: a local identity on a whole top-level term is
            // followed by the local identities on the terms it generated,
            // instead of waiting for the next planner round.
            let depth = self.repeat_depth.get();
            if depth == TERM_REPEAT_DEPTH {
                self.certificate.fixed.set(false);
                return rewritten;
            }
            self.repeat_depth.set(depth + 1);
            let repeated = self.rewrite_terms_at(rewritten.as_view(), complete, true);
            self.repeat_depth.set(depth);
            return repeated;
        }
        // A term left in place, for instance until a contraction, can still
        // change in a later round.
        self.certificate.fixed.set(false);

        // Try product-level rewrites first so trace*f contractions can fire before the
        // trace terminal expands into symmetric trace and f terms.
        let mut slots = SlotMatcher::default();
        let mut opaque_level = None;
        expr.replace_map(|arg, context, out| {
            if opaque_level.is_some_and(|level| context.function_level > level) {
                return;
            }
            opaque_level = None;
            if !matches!(slots.classify(arg), SlotMatch::Other) {
                opaque_level = Some(context.function_level);
                return;
            }
            // The root already declined above. Descendants still run in the
            // existing top-down order, stopping below each successful rewrite.
            if context.parent_type.is_some()
                && let Some(rewritten) = self.rewrite_node(arg, complete, false)
            {
                **out = rewritten;
            }
        })
    }

    fn sum_rewritten_terms(terms: Vec<Atom>) -> Atom {
        fn exact_coefficients(term: AtomView<'_>) -> bool {
            match term {
                AtomView::Num(number) => matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                ),
                AtomView::Add(sum) => sum.iter().all(exact_coefficients),
                AtomView::Mul(product) => product.iter().all(|factor| {
                    !matches!(factor, AtomView::Num(_)) || exact_coefficients(factor)
                }),
                _ => true,
            }
        }
        // Merging exact terms once avoids repeatedly normalizing the growing
        // prefix. Rounded coefficients retain their original addition order.
        if terms.iter().all(|term| exact_coefficients(term.as_view())) {
            Atom::add_many(terms)
        } else {
            terms.into_iter().fold(Atom::Zero, |sum, term| sum + term)
        }
    }

    /// `top` marks a whole top-level term of the row. With a complete
    /// frontier nothing outside such a term can reach its colour ports.
    fn rewrite_node(&self, arg: AtomView<'_>, complete: bool, top: bool) -> Option<Atom> {
        let whole_term = top && complete;
        let rewritten = self
            .simplify_product(arg, complete, whole_term)
            .or_else(|| self.simplify_chain_node(arg))
            .or_else(|| {
                let line = if whole_term {
                    LineContext::Isolated
                } else {
                    LineContext::Unknown
                };
                self.settings
                    .evaluate_traces
                    .then(|| self.simplify_trace_node(arg, complete, line))
                    .flatten()
            })
            .or_else(|| self.simplify_power(arg));
        // Only a complete decomposition of the whole term keeps the row fixed.
        let terminal = self.certificate.terminal.take() && whole_term;
        if rewritten.is_some() && !terminal {
            self.certificate.fixed.set(false);
        }
        self.certificate
            .terminal_term
            .set(rewritten.is_some() && terminal);
        rewritten
    }

    fn simplify_chain_node(&self, chain: AtomView) -> Option<Atom> {
        let (start, end, factors) = chain_parts(chain)?;
        if factors.is_empty()
            || factors
                .iter()
                .all(|factor| is_chain_identity_factor(factor.as_view()))
        {
            return Some(color_metric(start.clone(), end.clone()));
        }

        if let Some(identity_index) = factors
            .iter()
            .position(|factor| is_chain_identity_factor(factor.as_view()))
        {
            return Some(chain_with_factors(
                start,
                end,
                factors_excluding_indices(&factors, &[identity_index]),
            ));
        }

        if let Some(rewritten) = self.simplify_antisymmetric_chain_projector(&start, &end, &factors)
        {
            return Some(rewritten);
        }

        for i in 0..factors.len().saturating_sub(1) {
            let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
                continue;
            };
            let Some(right) = color_generator_adjoint(factors[i + 1].as_view()) else {
                continue;
            };
            if left != right {
                continue;
            }

            return Some(
                fundamental_casimir(fundamental_chain_dimension(&start, &end)?)
                    * chain_with_removed_range(&start, &end, &factors, i, i + 2),
            );
        }

        if factors.len() > 2
            && let Some(rewritten) = Self::simplify_separated_chain_casimir(&start, &end, &factors)
        {
            return Some(rewritten);
        }

        None
    }

    fn simplify_trace_node(
        &self,
        trace: AtomView,
        complete: bool,
        line: LineContext,
    ) -> Option<Atom> {
        let (rep, factors) = trace_parts(trace)?;

        if matches!(rep.as_view(), AtomView::Fun(f) if f.get_symbol() == CS.adjoint_rep) {
            let mut prefactor = Atom::one();
            let mut degree = Some(0usize);
            let normalized = factors
                .iter()
                .map(|factor| {
                    if let Some((coefficient, factor)) =
                        Self::normalize_adjoint_factor(rep.as_view(), factor.as_view())
                    {
                        prefactor *= coefficient;
                        if let Some(value) = degree.as_mut() {
                            *value += projector_parts(factor.as_view())
                                .map_or(1, |(_, factors)| factors.len());
                        }
                        factor
                    } else {
                        degree = None;
                        factor.clone()
                    }
                })
                .collect::<Vec<_>>();
            if !prefactor.is_one() || normalized != factors {
                return Some(prefactor * trace_with_factors(rep, normalized));
            }
            if let Some(degree) = degree {
                // Every certified raw adjoint matrix is antisymmetric. Trace
                // transposition reverses the ordered blocks; a symmetric block
                // only contributes the parity of its number of matrices.
                let reverse =
                    trace_with_factors(rep.clone(), factors.iter().rev().cloned().collect());
                if reverse.as_view() == trace && degree % 2 == 1 {
                    return Some(Atom::Zero);
                }
                if reverse.as_view() < trace {
                    return Some(if degree % 2 == 0 { reverse } else { -reverse });
                }
            }
        }

        if factors.is_empty()
            || factors
                .iter()
                .all(|factor| is_chain_identity_factor(factor.as_view()))
        {
            // Tr(1) and Tr(identity-line factors) collapse to the traced
            // representation dimension.
            return trace_terminal_dimension(rep.as_view());
        }

        if let Some(identity_index) = factors
            .iter()
            .position(|factor| is_chain_identity_factor(factor.as_view()))
        {
            // An identity line inside a longer trace is neutral:
            // Tr(... 1 ... ) -> Tr(...).
            return Some(trace_with_factors(
                rep,
                factors_excluding_indices(&factors, &[identity_index]),
            ));
        }

        // Raw adjoint matrices are F^a_bc = f^{bca} = i T_A^a_bc.
        // Their odd symmetric traces vanish. Keep the even symmetric traces
        // as invariant tensors; unfolding them would undo the terminal form.
        if let [factor] = factors.as_slice()
            && let Some(args) = color_symmetric_trace_arg_views(rep.as_view(), factor.as_view())
        {
            if symmetric_trace_phase(rep.as_view(), args.len()) == 0 {
                return Some(Atom::Zero);
            }
            if args.len() <= 2 {
                return self.simplify_generator_trace(
                    &rep,
                    &args
                        .into_iter()
                        .map(|slot| slot.to_owned())
                        .collect::<Vec<_>>(),
                    line,
                );
            }
        }

        // Through degree four, every cyclic ordering of a repeated adjoint pair
        // reduces by the adjacent or separated Casimir rules below. Preserve open
        // symmetric invariants and higher-degree projectors without factorial expansion.
        if let [factor] = factors.as_slice()
            && let Some(args) = color_symmetric_trace_arg_views(rep.as_view(), factor.as_view())
            && args.len() <= 4
            && let Some((left, right)) = args.iter().enumerate().find_map(|(i, a)| {
                args[i + 1..]
                    .iter()
                    .position(|b| a == b)
                    .map(|j| (i, i + 1 + j))
            })
        {
            if let AtomView::Fun(representation) = rep.as_view()
                && representation.get_symbol() == CS.adjoint_rep
            {
                // The degree-four adjoint invariant has a direct contracted
                // trace identity (color.h); raw f words are not fundamental T
                // matrices and must not enter the fundamental trace expansion.
                let open = args
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != left && *i != right)
                    .map(|(_, slot)| (*slot).to_owned())
                    .collect::<Vec<_>>();
                return Some(
                    Atom::num((5, 6))
                        * quadratic_casimir(rep.clone()).pow(2)
                        * color_metric(open[0].clone(), open[1].clone()),
                );
            }
            return Some(trace.expand_projectors());
        }

        let adjoint = matches!(rep.as_view(), AtomView::Fun(f) if f.get_symbol() == CS.adjoint_rep);
        let generators = if adjoint {
            factors
                .iter()
                .map(|factor| {
                    adjoint_generator_slot(rep.as_view(), factor.as_view())
                        .map(|slot| slot.to_owned())
                })
                .collect::<Option<Vec<_>>>()
        } else {
            factors
                .iter()
                .map(|factor| color_generator_adjoint(factor.as_view()))
                .collect()
        };
        if let Some(generators) = &generators {
            // Short Casimir and terminal rules require the same compatible
            // adjoint space as the general word recurrence.
            Self::trace_generator_is_adjoint(&rep, generators)?;
        }

        if let Some(rewritten) = self.simplify_antisymmetric_trace_projector(&rep, &factors) {
            return Some(rewritten);
        }

        if factors.len() > 2
            && let Some(rewritten) = Self::simplify_adjacent_trace_casimir(&rep, &factors)
        {
            return Some(rewritten);
        }

        if factors.len() > 3
            && let Some(rewritten) = Self::simplify_separated_trace_casimir(&rep, &factors)
        {
            return Some(rewritten);
        }

        if let Some(generators) = &generators
            && let Some(rewritten) = self.simplify_repeated_generator_trace(&rep, generators)
        {
            return Some(rewritten);
        }

        // Later selected factors can still contract with this ordered trace.
        // Keep its terminal decomposition until the frontier is complete;
        // identity, Casimir and product rewrites remain available on prefixes.
        if !complete {
            return None;
        }
        let Some(generators) = generators else {
            return self.simplify_prefixed_generator_trace(&rep, &factors, line);
        };
        if adjoint {
            return self.simplify_generator_trace(&rep, &generators, line);
        }

        // These terminals decompose an isolated line completely, as the
        // one-shot recursion does for longer ones.
        if generators.len() <= 4 && line == LineContext::Isolated && self.settings.one_shot_traces {
            self.certificate.terminal.set(true);
        }
        match generators.as_slice() {
            // Tr(T^a) -> 0.
            [_] => Some(Atom::Zero),
            // Tr(T^a T^b) -> idx(2,rep) g^{ab}.
            [a, b] => Some(quadratic_index(rep.clone()) * color_metric(a.clone(), b.clone())),
            // Tr(T^a T^b T^c) ->
            //   Tr(sym(T^a,T^b,T^c)) + i/2 idx(2,rep) f^{abc}.
            [a, b, c] => Some(
                color_symmetric_trace(&rep, [a.clone(), b.clone(), c.clone()])
                    + Atom::i() * Atom::num(1) / Atom::num(2)
                        * quadratic_index(rep.clone())
                        * color_f!(a.clone(), b.clone(), c.clone()),
            ),
            [a, b, c, d] => self.simplify_four_generator_trace_terminal(&rep, a, b, c, d),
            _ => self.simplify_generator_trace(&rep, &generators, line),
        }
    }

    fn normalize_adjoint_factor(rep: AtomView<'_>, factor: AtomView<'_>) -> Option<(Atom, Atom)> {
        if let Some((mut coefficient, symbol, factors)) = projector_factor(factor)
            && symbol == *shadowing::SYM
        {
            let normalized = factors
                .iter()
                .map(|factor| {
                    let (slot, phase) = adjoint_generator(rep, factor.as_view())?;
                    coefficient *= phase;
                    Some(color_f!(
                        Atom::var(T.chain_in),
                        Atom::var(T.chain_out),
                        slot
                    ))
                })
                .collect::<Option<Vec<_>>>()?;
            return Some((coefficient, shadowing::sym(normalized)));
        }
        let (slot, coefficient) = adjoint_generator(rep, factor)?;
        Some((
            coefficient,
            color_f!(Atom::var(T.chain_in), Atom::var(T.chain_out), slot),
        ))
    }

    fn simplify_adjacent_trace_casimir(rep: &Atom, factors: &[Atom]) -> Option<Atom> {
        for i in 0..factors.len() {
            let next = (i + 1) % factors.len();
            let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
                continue;
            };
            let Some(right) = color_generator_adjoint(factors[next].as_view()) else {
                continue;
            };
            if left != right {
                continue;
            }

            // Adjacent equal generators inside a fundamental trace:
            // Tr(... T^a T^a ...) -> cas(2,rep) Tr(...).
            return Some(
                quadratic_casimir(rep.clone())
                    * trace_with_factors(
                        rep.clone(),
                        factors_excluding_indices(factors, &[i, next]),
                    ),
            );
        }

        None
    }

    fn simplify_antisymmetric_chain_projector(
        &self,
        start: &Atom,
        end: &Atom,
        factors: &[Atom],
    ) -> Option<Atom> {
        // `antisym(T^a,T^b)` is the normalized commutator: i/2 f^{abx} T^x.
        for (position, factor) in factors.iter().enumerate() {
            let Some((prefactor, args)) = color_antisymmetric_generator_args(factor.as_view())
            else {
                continue;
            };
            let [a, b] = args.as_slice() else {
                continue;
            };
            let x = self.color_adjoint_dummy_for_pair(a, b)?;

            let mut replacement_factors = factors.to_vec();
            replacement_factors[position] = color_t!(x.clone());
            return Some(
                prefactor * Atom::i() / Atom::num(2)
                    * color_f!(a.clone(), b.clone(), x.clone())
                    * chain_with_factors(start.clone(), end.clone(), replacement_factors),
            );
        }

        None
    }

    fn simplify_antisymmetric_trace_projector(&self, rep: &Atom, factors: &[Atom]) -> Option<Atom> {
        if let [factor] = factors {
            let (prefactor, args) = color_antisymmetric_generator_args(factor.as_view())?;
            return match args.as_slice() {
                // Tr(antisym(T^a,T^b)) is the trace of a commutator.
                [_, _] => Some(Atom::Zero),
                // Tr(antisym(T^a,T^b,T^c)) -> i/2 idx(2,rep) f^{abc}.
                [a, b, c] => Some(
                    prefactor
                        * Atom::i()
                        * quadratic_index(rep.clone())
                        * color_f!(a.clone(), b.clone(), c.clone())
                        / Atom::num(2),
                ),
                _ => None,
            };
        }

        for (position, factor) in factors.iter().enumerate() {
            let Some((prefactor, args)) = color_antisymmetric_generator_args(factor.as_view())
            else {
                continue;
            };
            let [a, b] = args.as_slice() else {
                continue;
            };
            let x = self.color_adjoint_dummy_for_pair(a, b)?;

            let mut replacement_factors = factors.to_vec();
            replacement_factors[position] = color_t!(x.clone());
            // In a longer trace, antisym(T^a,T^b) is the normalized
            // commutator: i/2 f^{abx} T^x.
            return Some(
                prefactor * Atom::i() / Atom::num(2)
                    * color_f!(a.clone(), b.clone(), x.clone())
                    * trace_with_factors(rep.clone(), replacement_factors),
            );
        }

        None
    }

    fn simplify_separated_trace_casimir(rep: &Atom, factors: &[Atom]) -> Option<Atom> {
        // A trace has no distinguished starting point. Test the closing
        // connection too, independently of the labels chosen by cyclic ordering.
        for i in 0..factors.len() {
            let middle = (i + 1) % factors.len();
            let right_index = (i + 2) % factors.len();
            let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
                continue;
            };
            let Some(right) = color_generator_adjoint(factors[right_index].as_view()) else {
                continue;
            };
            let Some(middle) = color_generator_adjoint(factors[middle].as_view()) else {
                continue;
            };
            if left != right || color_adjoint_dimension(&left) != color_adjoint_dimension(&middle) {
                continue;
            }

            // Separated equal generators with one generator between them:
            // Tr(... T^a T^b T^a ...) -> (cas(2,rep) - cas(2,adj)/2) Tr(... T^b ...).
            return Some(
                (quadratic_casimir(rep.clone())
                    - adjoint_casimir_for_dimension(color_adjoint_dimension(&left)?)
                        / Atom::num(2))
                    * trace_with_factors(
                        rep.clone(),
                        factors_excluding_indices(factors, &[i, right_index]),
                    ),
            );
        }

        None
    }

    fn simplify_separated_chain_casimir(
        start: &Atom,
        end: &Atom,
        factors: &[Atom],
    ) -> Option<Atom> {
        for i in 0..factors.len().saturating_sub(2) {
            let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
                continue;
            };
            let Some(right) = color_generator_adjoint(factors[i + 2].as_view()) else {
                continue;
            };
            let Some(middle) = color_generator_adjoint(factors[i + 1].as_view()) else {
                continue;
            };
            if left != right || color_adjoint_dimension(&left) != color_adjoint_dimension(&middle) {
                continue;
            }

            let coefficient = fundamental_casimir(fundamental_chain_dimension(start, end)?)
                - adjoint_casimir_for_dimension(color_adjoint_dimension(&left)?) / Atom::num(2);

            return Some(
                coefficient
                    * chain_with_factors(
                        start.clone(),
                        end.clone(),
                        factors_excluding_indices(factors, &[i, i + 2]),
                    ),
            );
        }

        None
    }

    fn simplify_four_generator_trace_terminal(
        &self,
        rep: &Atom,
        a: &Atom,
        b: &Atom,
        c: &Atom,
        d: &Atom,
    ) -> Option<Atom> {
        let x = self.color_adjoint_dummy_like(a)?;

        // Four-generator terminal decomposition:
        // Tr(T^a T^b T^c T^d) is split into a fully symmetric trace, two
        // f * symmetric-trace terms, and two f*f terms.
        Some(
            color_symmetric_trace(rep, [a.clone(), b.clone(), c.clone(), d.clone()])
                + Atom::i() / Atom::num(2)
                    * color_symmetric_trace(rep, [a.clone(), b.clone(), x.clone()])
                    * color_f!(c.clone(), d.clone(), x.clone())
                + Atom::i() / Atom::num(2)
                    * color_symmetric_trace(rep, [c.clone(), d.clone(), x.clone()])
                    * color_f!(a.clone(), b.clone(), x.clone())
                - quadratic_index(rep.clone()) / Atom::num(6)
                    * color_f!(a.clone(), c.clone(), x.clone())
                    * color_f!(b.clone(), d.clone(), x.clone())
                + quadratic_index(rep.clone()) / Atom::num(3)
                    * color_f!(a.clone(), d.clone(), x.clone())
                    * color_f!(b.clone(), c.clone(), x.clone()),
        )
    }

    fn simplify_product(
        &self,
        product: AtomView,
        complete: bool,
        whole_term: bool,
    ) -> Option<Atom> {
        if !matches!(product, AtomView::Mul(_)) {
            return None;
        }
        let product = ProductView::parse(product);

        Self::simplify_symmetric_prefix_structure_product(&product)
            .or_else(|| Self::join_color_chain_product(&product))
            .or_else(|| {
                self.settings
                    .evaluate_traces
                    .then(|| self.simplify_trace_structure_product(&product))
                    .flatten()
            })
            .or_else(|| self.simplify_chain_structure_product(&product))
            .or_else(|| Self::simplify_symmetric_structure_product(&product))
            .or_else(|| Self::simplify_symmetric_structure_pair_product(&product))
            .or_else(|| self.simplify_two_f_loop_product(&product))
            .or_else(|| self.simplify_adjoint_loop_product(&product, whole_term))
            .or_else(|| Self::simplify_symmetric_invariant_product(&product))
            .or_else(|| {
                self.settings
                    .expand_cross_chain_fierz
                    .then(|| Self::simplify_cross_chain_fierz_product(&product, complete))
                    .flatten()
                    .inspect(|_| self.certificate.expansion.set(true))
            })
            .or_else(|| self.distribute_color_sum_factor(&product, complete, whole_term))
            .or_else(|| self.simplify_embedded_color_node(&product, complete, whole_term))
    }

    /// Multiply a colour sum factor into the other factors, so that a product
    /// rule can reach each summand. A rewrite does not revisit its own result
    /// in the same pass, so a summand carrying a nested sum used to wait for
    /// the next colour round, one nesting level per round. In one-shot mode
    /// the sum is distributed only when a rule can span a summand and the
    /// other factors; on a complete frontier the generated terms are then
    /// rewritten in the same pass, including nested sums such a rule needs.
    fn distribute_color_sum_factor(
        &self,
        product: &ProductView,
        complete: bool,
        whole_term: bool,
    ) -> Option<Atom> {
        let selective = self.settings.one_shot_traces;
        let (index, sum) = product.color_sum_factor(selective)?;
        let outside = SlotCounts::of_factors(
            product
                .factors
                .iter()
                .enumerate()
                .filter(|(other, _)| *other != index)
                .map(|(_, factor)| factor.atom),
        );
        let distributed = Self::sum_rewritten_terms(
            sum.iter()
                .map(|term| {
                    let term =
                        product.replacing_one(index, self.without_captured_dummies(term, &outside));
                    if selective && complete {
                        // Distributing a whole top-level term yields whole terms.
                        self.rewrite_terms_at(term.as_view(), complete, whole_term)
                    } else {
                        term
                    }
                })
                .collect(),
        );
        self.certificate.expansion.set(true);
        Some(distributed)
    }

    fn join_color_chain_product(product: &ProductView) -> Option<Atom> {
        for (left_index, left_factor) in product.factors.iter().enumerate() {
            let Some(left_chain) = &left_factor.chain else {
                continue;
            };
            let Some((left_end_dim, left_end_index, true)) = color_fundamental_slot(left_chain.end)
            else {
                continue;
            };

            for (right_index, right_factor) in product.factors.iter().enumerate() {
                if right_index == left_index {
                    continue;
                }
                let Some(right_chain) = &right_factor.chain else {
                    continue;
                };
                let Some((right_start_dim, right_start_index, false)) =
                    color_fundamental_slot(right_chain.start)
                else {
                    continue;
                };
                if left_end_dim != right_start_dim || left_end_index != right_start_index {
                    continue;
                }

                let replacement = chain!(
                    left_chain.start,
                    right_chain.end;
                    left_chain.factors.iter().cloned().chain(right_chain.factors.iter().cloned())
                );
                return Some(product.replacing_pair(left_index, right_index, replacement));
            }
        }

        None
    }

    /// Rewrite one colour node of the product. Lines go first, the one with
    /// the fewest distinct external generators leading, as color.h's cOlTT
    /// picks the trace with the fewest different indices: decomposing it
    /// attaches few new structures to the other lines.
    fn simplify_embedded_color_node(
        &self,
        product: &ProductView,
        complete: bool,
        whole_term: bool,
    ) -> Option<Atom> {
        let mut order = (0..product.len()).collect::<Vec<_>>();
        order.sort_by_cached_key(|&index| product.factors[index].line_cost());
        for index in order {
            let factor = &product.factors[index];
            // The parsed views already tell which node rule can apply.
            let rewritten = if factor.chain.is_some() {
                self.simplify_chain_node(factor.atom)
            } else if factor.trace.is_some() {
                self.settings
                    .evaluate_traces
                    .then(|| {
                        let line = product.line_context(|other| other == index, whole_term);
                        self.simplify_trace_node(factor.atom, complete, line)
                    })
                    .flatten()
            } else {
                self.simplify_power(factor.atom)
            };
            let Some(rewritten) = rewritten else {
                continue;
            };

            return Some(product.replacing_one(index, rewritten));
        }

        None
    }

    /// With an incomplete frontier a chain may still be continued by a
    /// factor outside the product, so only traces count as whole lines.
    fn simplify_cross_chain_fierz_product(product: &ProductView, complete: bool) -> Option<Atom> {
        let lines = product
            .factors
            .iter()
            .enumerate()
            .filter_map(|(index, factor)| {
                let (dimension, factors) = if let Some(chain) = &factor.chain {
                    (
                        fundamental_chain_dimension_view(chain.start, chain.end)?,
                        &chain.factors,
                    )
                } else if let Some(trace) = &factor.trace {
                    let AtomView::Fun(rep) = trace.rep else {
                        return None;
                    };
                    if rep.get_symbol() != CS.fundamental_rep || rep.get_nargs() != 1 {
                        return None;
                    }
                    (rep.iter().next()?.to_owned(), &trace.factors)
                } else {
                    return None;
                };
                Some((index, factor, dimension, factors, factor.generator_slots()?))
            })
            // FORM joins T*T and closes T(i,i) into traces first: a chain
            // continued by another line or a metric is not a whole line.
            .filter(|(index, factor, ..)| {
                factor
                    .chain
                    .as_ref()
                    .is_none_or(|chain| complete && product.chain_is_maximal(*index, chain))
            })
            .collect::<Vec<_>>();

        for (position, (left_index, left, left_dimension, left_factors, left_slots)) in
            lines.iter().enumerate()
        {
            for (right_index, right, right_dimension, right_factors, right_slots) in
                &lines[position + 1..]
            {
                if left_dimension != right_dimension {
                    continue;
                }
                for (left_position, left_generator) in left_factors.iter().enumerate() {
                    let Some(left_adjoint) = color_generator_adjoint(*left_generator) else {
                        continue;
                    };
                    if left_slots
                        .iter()
                        .filter(|slot| **slot == left_adjoint.as_view())
                        .count()
                        != 1
                        || right_slots
                            .iter()
                            .filter(|slot| **slot == left_adjoint.as_view())
                            .count()
                            != 1
                    {
                        continue;
                    }
                    for (right_position, right_generator) in right_factors.iter().enumerate() {
                        if color_generator_adjoint(*right_generator).as_ref() != Some(&left_adjoint)
                        {
                            continue;
                        }
                        let left_before = &left_factors[..left_position];
                        let left_after = &left_factors[left_position + 1..];
                        let right_before = &right_factors[..right_position];
                        let right_after = &right_factors[right_position + 1..];
                        let dimension = left_dimension.to_owned();
                        let rep = fundamental_rep(dimension.clone());
                        let left_trace = trace_with_factors(
                            rep.clone(),
                            left_after
                                .iter()
                                .chain(left_before)
                                .map(|factor| (*factor).to_owned())
                                .collect(),
                        );
                        let right_trace = trace_with_factors(
                            rep.clone(),
                            right_after
                                .iter()
                                .chain(right_before)
                                .map(|factor| (*factor).to_owned())
                                .collect(),
                        );
                        // Cutting a trace at the contracted generator leaves its
                        // remaining word in cyclic order: after, then before.
                        let (crossed, uncrossed) = match (&left.chain, &right.chain) {
                            (Some(left), Some(right)) => (
                                chain_with_factor_view_slices(
                                    left.start,
                                    right.end,
                                    &[left_before, right_after],
                                ) * chain_with_factor_view_slices(
                                    right.start,
                                    left.end,
                                    &[right_before, left_after],
                                ),
                                chain_with_factor_view_slices(
                                    left.start,
                                    left.end,
                                    &[left_before, left_after],
                                ) * chain_with_factor_view_slices(
                                    right.start,
                                    right.end,
                                    &[right_before, right_after],
                                ),
                            ),
                            (Some(left), None) => (
                                chain_with_factor_view_slices(
                                    left.start,
                                    left.end,
                                    &[left_before, right_after, right_before, left_after],
                                ),
                                chain_with_factor_view_slices(
                                    left.start,
                                    left.end,
                                    &[left_before, left_after],
                                ) * right_trace,
                            ),
                            (None, Some(right)) => (
                                chain_with_factor_view_slices(
                                    right.start,
                                    right.end,
                                    &[right_before, left_after, left_before, right_after],
                                ),
                                left_trace
                                    * chain_with_factor_view_slices(
                                        right.start,
                                        right.end,
                                        &[right_before, right_after],
                                    ),
                            ),
                            (None, None) => (
                                trace_with_factors(
                                    rep,
                                    left_after
                                        .iter()
                                        .chain(left_before)
                                        .chain(right_after)
                                        .chain(right_before)
                                        .map(|factor| (*factor).to_owned())
                                        .collect(),
                                ),
                                left_trace * right_trace,
                            ),
                        };
                        let replacement = fundamental_index(dimension.clone())
                            * (crossed - uncrossed / dimension);
                        return Some(product.replacing_pair(
                            *left_index,
                            *right_index,
                            replacement,
                        ));
                    }
                }
            }
        }
        None
    }

    fn simplify_power(&self, power: AtomView) -> Option<Atom> {
        let AtomView::Pow(pow) = power else {
            return None;
        };
        let (base, exponent) = pow.get_base_exp();
        let exponent = positive_integer(exponent)?;
        if exponent >= 2
            && let AtomView::Add(sum) = base
            && let Some(ports) = colour_power_ports(base)
        {
            // X^n = sum_i x_i X^(n-1). The copies contract through the ports
            // of X; each x_i gets fresh internal dummies so that it cannot
            // capture a same-spelled dummy of the remaining copy.
            let remaining = if exponent == 2 {
                base.to_owned()
            } else {
                base.to_owned().pow(Atom::num(exponent - 1))
            };
            return Some(Atom::add_many(sum.iter().map(|term| {
                self.with_fresh_dummies(term, |slot| !ports.contains(slot)) * &remaining
            })));
        }
        if exponent != 2 {
            return None;
        }

        if let Some(invariant) = SymmetricInvariantView::parse(base)
            && invariant.args.len() >= 3
            && invariant.has_distinct_args()
        {
            return Some(
                Atom::num(invariant.phase * invariant.phase)
                    * color_symmetric_product(
                        invariant.args.len(),
                        invariant.rep.to_owned(),
                        invariant.rep.to_owned(),
                    ),
            );
        }

        let args = structure_constant_args(base)?;
        let dimension = color_structure_dimension(&args)?;
        Some(adjoint_casimir_for_dimension(dimension.clone()) * dimension)
    }

    /// Rename the indices of the slots of `term` selected by `rename` to fresh
    /// dummies, one per label. Slots are whole: a dummy of another
    /// representation may reuse a kept slot's label and is still renamed.
    fn with_fresh_dummies(
        &self,
        term: AtomView<'_>,
        rename: impl Fn(&Slot<LibraryRep, AbstractIndex>) -> bool,
    ) -> Atom {
        let mut slots = SlotMatcher::default();
        let mut fresh = std::collections::HashMap::new();
        term.replace_map(|node, _context, out| {
            let Ok(mut slot) = slots.parse::<LibraryRep, AbstractIndex>(node) else {
                return;
            };
            if !rename(&slot) {
                return;
            }
            slot.aind = *fresh
                .entry(slot.aind)
                .or_insert_with(|| self.dummies.fresh_index());
            **out = slot.to_atom();
        })
    }

    /// A summand multiplied into the other factors of a product: its own
    /// dummies that another factor also spells, for instance in two copies
    /// of one sum, are renamed so that the product cannot capture them. So
    /// are the labels of a power, a closed scope, that the summand or
    /// another factor spells elsewhere.
    fn without_captured_dummies(&self, summand: AtomView<'_>, outside: &SlotCounts<'_>) -> Atom {
        let closed = self.with_closed_powers(summand, outside);
        let summand = closed.as_ref().map_or(summand, Atom::as_view);
        let captured = SlotCounts::shared_with(summand, outside).captured_by(outside);
        if captured.is_empty() {
            return summand.to_owned();
        }
        let mut matcher = SlotMatcher::default();
        let folded = |slot: Slot<LibraryRep, AbstractIndex>| {
            if slot.rep_name().is_dual() {
                slot.dual()
            } else {
                slot
            }
        };
        let captured = captured
            .into_iter()
            .filter_map(|node| matcher.parse::<LibraryRep, AbstractIndex>(node).ok())
            .map(folded)
            .collect::<std::collections::HashSet<_>>();
        self.with_fresh_dummies(summand, |slot| captured.contains(&folded(*slot)))
    }

    /// A power of a tensor contracts its copies with each other, so all of
    /// its labels are bound inside it. Give fresh labels to each power of
    /// `expression` that shares one with the rest of `expression`, with
    /// `outside` or with an earlier power, so that no rule opening it can
    /// capture them; `None` when no power needs it.
    fn with_closed_powers(
        &self,
        expression: AtomView<'_>,
        outside: &SlotCounts<'_>,
    ) -> Option<Atom> {
        fn powers<'a>(expression: AtomView<'a>, found: &mut Vec<AtomView<'a>>) {
            match expression {
                AtomView::Add(sum) => sum.iter().for_each(|term| powers(term, found)),
                AtomView::Mul(product) => product.iter().for_each(|factor| powers(factor, found)),
                AtomView::Pow(_) => found.push(expression),
                _ => {}
            }
        }
        let counts = SlotCounts::of(expression);
        if counts.bound.is_empty() {
            return None;
        }
        let mut found = Vec::new();
        powers(expression, &mut found);
        let mut taken = counts
            .counts
            .iter()
            .map(|(key, ..)| *key)
            .chain(outside.counts.iter().map(|(key, ..)| *key))
            .chain(outside.bound.iter().copied())
            .collect::<Vec<_>>();
        let mut renamed = std::collections::HashMap::new();
        for power in found {
            let labels = SlotCounts::of(power).bound;
            if labels.iter().any(|label| taken.contains(label)) {
                renamed
                    .entry(power.to_owned())
                    .or_insert_with(|| self.with_fresh_dummies(power, |_| true));
            } else {
                taken.extend(labels);
            }
        }
        if renamed.is_empty() {
            return None;
        }
        Some(expression.replace_map(|node, _context, out| {
            if matches!(node, AtomView::Pow(_))
                && let Some(fresh) = renamed.get(&node.to_owned())
            {
                **out = fresh.clone();
            }
        }))
    }

    fn simplify_trace_structure_product(&self, product: &ProductView) -> Option<Atom> {
        let mut best: Option<(AbsorbedLine<'_>, bool)> = None;
        for (trace_index, trace_factor) in product.factors.iter().enumerate() {
            let Some(trace) = &trace_factor.trace else {
                continue;
            };
            if trace.factors.len() < 2 {
                continue;
            }
            let Some(slots) = trace_factor.generator_slots() else {
                continue;
            };
            let adjoint = matches!(trace.rep, AtomView::Fun(f) if f.get_symbol() == CS.adjoint_rep);
            let words = trace
                .factors
                .iter()
                .map(|&factor| {
                    if adjoint {
                        adjoint_generator_slot(trace.rep, factor)
                    } else {
                        color_generator_adjoint_view(factor)
                    }
                })
                .collect::<Vec<_>>();
            if let Some((distance, structure, first)) =
                Self::closest_absorbed_pair(product, trace_index, &words, &slots, true)
                && best
                    .as_ref()
                    .is_none_or(|(closest, _)| distance < closest.distance)
            {
                let line = AbsorbedLine {
                    distance,
                    line: trace_index,
                    structure,
                    first,
                    words,
                };
                best = Some((line, adjoint));
            }
        }
        let (
            AbsorbedLine {
                distance,
                line: trace_index,
                structure: f_index,
                first,
                words,
            },
            adjoint,
        ) = best?;
        let trace = product.factors[trace_index].trace.as_ref()?;
        // Read the word from the first leg of the pair: X^a W X^c R.
        let n = trace.factors.len();
        let rotated = (0..n).map(|k| (first + k) % n).collect::<Vec<_>>();
        let replacement = self.absorbed_pair_terms(
            product,
            f_index,
            &[],
            &rotated
                .iter()
                .map(|&k| trace.factors[k])
                .collect::<Vec<_>>(),
            &rotated.iter().map(|&k| words[k]).collect::<Vec<_>>(),
            distance,
            adjoint,
            |factors| trace_with_factors(trace.rep.to_owned(), factors),
        )?;
        Some(product.replacing_pair(trace_index, f_index, replacement))
    }

    fn simplify_chain_structure_product(&self, product: &ProductView) -> Option<Atom> {
        let mut best: Option<AbsorbedLine<'_>> = None;
        for (chain_index, chain_factor) in product.factors.iter().enumerate() {
            let Some(chain) = &chain_factor.chain else {
                continue;
            };
            let Some(slots) = chain_factor.generator_slots() else {
                continue;
            };
            let words = chain
                .factors
                .iter()
                .map(|&factor| color_generator_adjoint_view(factor))
                .collect::<Vec<_>>();
            if let Some((distance, structure, first)) =
                Self::closest_absorbed_pair(product, chain_index, &words, &slots, false)
                && best
                    .as_ref()
                    .is_none_or(|closest| distance < closest.distance)
            {
                best = Some(AbsorbedLine {
                    distance,
                    line: chain_index,
                    structure,
                    first,
                    words,
                });
            }
        }
        let AbsorbedLine {
            distance,
            line: chain_index,
            structure: f_index,
            first,
            words,
        } = best?;
        let chain = product.factors[chain_index].chain.as_ref()?;
        let replacement = self.absorbed_pair_terms(
            product,
            f_index,
            &chain.factors[..first],
            &chain.factors[first..],
            &words[first..],
            distance,
            false,
            |factors| chain_with_factors(chain.start.to_owned(), chain.end.to_owned(), factors),
        )?;
        Some(product.replacing_pair(chain_index, f_index, replacement))
    }

    /// The closest pair of generators of one line that is also a pair of legs
    /// of a structure constant of the product, as (distance, structure index,
    /// first position). `words` holds each factor's generator slot; a trace
    /// is `cyclic` and is read along its shorter arc. The pair and the
    /// generators between them must be unique slots of the line. Beyond
    /// adjacent pairs every factor of the line must be a generator: a
    /// symmetric block belongs to the symmetric-prefix rule.
    fn closest_absorbed_pair(
        product: &ProductView<'_>,
        line_index: usize,
        words: &[Option<AtomView<'_>>],
        slots: &[AtomView<'_>],
        cyclic: bool,
    ) -> Option<(usize, usize, usize)> {
        let n = words.len();
        let unique = |slot: &AtomView<'_>| slots.iter().filter(|other| *other == slot).count() == 1;
        let plain = words.iter().all(Option::is_some);
        let mut best: Option<(usize, usize, usize)> = None;
        for (f_index, f_factor) in product.factors.iter().enumerate() {
            if f_index == line_index {
                continue;
            }
            let Some(structure) = &f_factor.structure else {
                continue;
            };
            if color_structure_dimension(&structure.args).is_none() {
                continue;
            }
            for first in 0..n {
                let Some(a) = words[first] else {
                    continue;
                };
                if !structure.args.contains(&a) || !unique(&a) {
                    continue;
                }
                let longest = if cyclic { n / 2 } else { n - 1 - first };
                for distance in 1..=longest {
                    if best.is_some_and(|(closest, ..)| closest <= distance)
                        || (distance > 1 && !plain)
                    {
                        break;
                    }
                    let Some(c) = words[(first + distance) % n] else {
                        continue;
                    };
                    if c == a
                        || !structure.args.contains(&c)
                        || !unique(&c)
                        || (1..distance)
                            .any(|k| !words[(first + k) % n].is_some_and(|w| unique(&w)))
                    {
                        continue;
                    }
                    best = Some((distance, f_index, first));
                }
            }
        }
        best
    }

    /// Absorb the structure constant f at `f_index` into a line read as
    /// `before`, then X^a W X^c R in `word` (with `slots` its generator
    /// slots), where W has `distance - 1` generators. With
    /// [X^a, X^b] = κ f^{abc} X^c (κ = i for T, -1 for raw adjoint F):
    ///   X^a W X^c R f^{ace} = κ (C_A/2) X^e W R
    ///       + κ Σ_k f^{w_k c y_k} f^{ace} X^a W[w_k -> y_k] R,
    /// moving X^c past W and using X^a X^c f^{ace} = ½ [X^a, X^c] f^{ace}.
    /// The stored f keeps its orientation in the commutator terms; only the
    /// Casimir term reads it, positionally.
    #[allow(clippy::too_many_arguments)]
    fn absorbed_pair_terms(
        &self,
        product: &ProductView<'_>,
        f_index: usize,
        before: &[AtomView<'_>],
        word: &[AtomView<'_>],
        slots: &[Option<AtomView<'_>>],
        distance: usize,
        adjoint: bool,
        line: impl Fn(Vec<Atom>) -> Atom,
    ) -> Option<Atom> {
        let structure = product.factors[f_index].structure.as_ref()?;
        let (a, c) = (slots[0]?, slots[distance]?);
        let (target, orientation) =
            Self::structure_target_for_generator_pair(&structure.args, &a, &c)?;
        let casimir = adjoint_casimir_for_dimension(color_structure_dimension(&structure.args)?);
        let kappa = if adjoint { -Atom::one() } else { Atom::i() };
        let generator = |slot: &Atom| {
            if adjoint {
                color_f!(Atom::var(T.chain_in), Atom::var(T.chain_out), slot)
            } else {
                color_t!(slot)
            }
        };
        let between = &word[1..distance];
        let after = &word[distance + 1..];
        let owned = |factors: &[AtomView<'_>]| {
            factors
                .iter()
                .map(|factor| factor.to_owned())
                .collect::<Vec<_>>()
        };
        let mut terms = vec![
            orientation * &kappa * casimir / Atom::num(2)
                * line(
                    owned(before)
                        .into_iter()
                        .chain(std::iter::once(generator(&target)))
                        .chain(owned(between))
                        .chain(owned(after))
                        .collect(),
                ),
        ];
        for k in 1..distance {
            let w = slots[k]?;
            let y = self.color_adjoint_dummy_like(&c.to_owned())?;
            let mut moved = owned(between);
            moved[k - 1] = generator(&y);
            terms.push(
                &kappa
                    * color_f!(w, c, &y)
                    * product.factors[f_index].atom
                    * line(
                        owned(before)
                            .into_iter()
                            .chain(std::iter::once(word[0].to_owned()))
                            .chain(moved)
                            .chain(owned(after))
                            .collect(),
                    ),
            );
        }
        Some(Atom::add_many(terms))
    }

    fn structure_target_for_generator_pair(
        args: &[AtomView<'_>; 3],
        left: &AtomView<'_>,
        right: &AtomView<'_>,
    ) -> Option<(Atom, Atom)> {
        let left = args.iter().position(|arg| arg == left)?;
        let right = args.iter().position(|arg| arg == right)?;
        if left == right {
            return None;
        }
        let target = (0..3).find(|&i| i != left && i != right)?;
        Some((
            args[target].to_owned(),
            StructureView::orientation(left, right),
        ))
    }

    fn simplify_symmetric_structure_product(product: &ProductView) -> Option<Atom> {
        for (symmetric_index, symmetric_factor) in product.factors.iter().enumerate() {
            let Some(symmetric) = &symmetric_factor.symmetric_invariant else {
                continue;
            };

            for (f_index, f_factor) in product.factors.iter().enumerate() {
                if f_index == symmetric_index {
                    continue;
                }
                let Some(structure) = &f_factor.structure else {
                    continue;
                };
                if color_structure_dimension(&structure.args.map(|arg| arg.to_owned())).is_none() {
                    continue;
                }
                // A repeated slot is bound inside the symmetric trace. It
                // cannot connect to a same-spelled index outside that trace.
                let common_count = structure
                    .args
                    .iter()
                    .filter(|arg| symmetric.args.iter().filter(|slot| *slot == *arg).count() == 1)
                    .count();
                if common_count >= 2 {
                    return Some(Atom::Zero);
                }
                if common_count != 1 || !symmetric.has_distinct_args() {
                    continue;
                }
                let bridge = structure
                    .args
                    .iter()
                    .find(|arg| symmetric.args.contains(arg))
                    .expect("one shared structure-constant leg");
                // Invariance gives sum_i f(x,a_i,b) D(...,b,...) = 0.
                // Contract every a_i with a second symmetric invariant: all
                // terms coincide, so each is zero. This includes
                // f(x,a,d) D_R(a,b,c) D_S(d,b,c,e), independently of the
                // representations and the order used to decompose traces.
                for (other_index, other_factor) in product.factors.iter().enumerate() {
                    if other_index == symmetric_index {
                        continue;
                    }
                    let Some(other) = &other_factor.symmetric_invariant else {
                        continue;
                    };
                    if !other.has_distinct_args() {
                        continue;
                    }
                    if !other.args.contains(bridge)
                        && symmetric
                            .args
                            .iter()
                            .all(|arg| arg == bridge || other.args.contains(arg))
                        && structure.args.iter().any(|arg| other.args.contains(arg))
                    {
                        return Some(Atom::Zero);
                    }
                }
            }
        }

        None
    }

    /// A rank-three symmetric invariant whose two legs meet two structure
    /// constants sharing one index (color.h, simpli):
    /// d_R^{abx} f^{aik} f^{bjk} = (C_A/2) d_R^{ijx}.
    /// Invariance of d_R under the adjoint action gives the coefficient; it
    /// is the same for the Hermitian symmetric trace of any representation.
    /// Raw adjoint words have no odd symmetric invariant.
    fn simplify_symmetric_structure_pair_product(product: &ProductView) -> Option<Atom> {
        let structures = product
            .factors
            .iter()
            .enumerate()
            .filter_map(|(index, factor)| Some((index, factor.structure.as_ref()?)))
            .collect::<Vec<_>>();
        if structures.len() < 2 {
            return None;
        }
        for (symmetric_index, factor) in product.factors.iter().enumerate() {
            let Some(symmetric) = &factor.symmetric_invariant else {
                continue;
            };
            if symmetric.args.len() != 3 || symmetric.phase != 1 || !symmetric.has_distinct_args() {
                continue;
            }
            for (a_leg, b_leg, x_leg) in [(0, 1, 2), (0, 2, 1), (1, 2, 0)] {
                let [a, b, x] = [a_leg, b_leg, x_leg].map(|leg| symmetric.args[leg]);
                // A structure constant on two legs of d_R annihilates it; that
                // rule runs first.
                let on_leg = |structure: &StructureView<'_>, leg: AtomView<'_>| {
                    let position = structure.args.iter().position(|arg| *arg == leg)?;
                    (structure
                        .args
                        .iter()
                        .filter(|arg| symmetric.args.contains(arg))
                        .count()
                        == 1)
                        .then_some(position)
                };
                for &(left_index, left) in &structures {
                    let Some(left_leg) = on_leg(left, a) else {
                        continue;
                    };
                    for &(right_index, right) in &structures {
                        if right_index == left_index {
                            continue;
                        }
                        let Some(right_leg) = on_leg(right, b) else {
                            continue;
                        };
                        let shared = (0..3)
                            .filter(|&position| position != left_leg)
                            .flat_map(|position| {
                                (0..3)
                                    .filter(move |&other| {
                                        other != right_leg
                                            && left.args[position] == right.args[other]
                                    })
                                    .map(move |other| (position, other))
                            })
                            .collect::<Vec<_>>();
                        let [(left_shared, right_shared)] = shared.as_slice() else {
                            continue;
                        };
                        let Some(dimension) =
                            color_structure_dimension(&left.args).filter(|dimension| {
                                color_structure_dimension(&right.args).as_ref() == Some(dimension)
                            })
                        else {
                            continue;
                        };
                        // Read both stored orientations as f(a,i,k), f(b,j,k).
                        let i = 3 - left_leg - left_shared;
                        let j = 3 - right_leg - right_shared;
                        let sign = StructureView::orientation(left_leg, i)
                            * StructureView::orientation(right_leg, j);
                        let ports = [left.args[i], right.args[j], x].map(|port| port.to_owned());
                        let invariant = if let Some(trace) = &factor.trace {
                            color_symmetric_trace(&trace.rep.to_owned(), ports)
                        } else {
                            CS.symmetric_d(symmetric.rep, ports.to_vec())
                        };
                        let mut excluded = vec![false; product.len()];
                        for index in [symmetric_index, left_index, right_index] {
                            excluded[index] = true;
                        }
                        return Some(product.excluding(
                            &excluded,
                            sign * adjoint_casimir_for_dimension(dimension) / Atom::num(2)
                                * invariant,
                        ));
                    }
                }
            }
        }
        None
    }

    fn simplify_two_f_loop_product(&self, product: &ProductView) -> Option<Atom> {
        for (left_index, left_factor) in product.factors.iter().enumerate() {
            let Some(left_structure) = &left_factor.structure else {
                continue;
            };
            for (right_index, right_factor) in
                product.factors.iter().enumerate().skip(left_index + 1)
            {
                let Some(right_structure) = &right_factor.structure else {
                    continue;
                };
                let Some(replacement) = left_structure.contract_loop(right_structure) else {
                    continue;
                };
                return Some(product.replacing_pair(left_index, right_index, replacement));
            }
        }

        None
    }

    fn simplify_adjoint_loop_product(
        &self,
        product: &ProductView,
        whole_term: bool,
    ) -> Option<Atom> {
        let structures = product
            .factors
            .iter()
            .enumerate()
            .filter_map(|(index, factor)| {
                let args = factor.structure.as_ref()?.args.map(|arg| arg.to_owned());
                let dimension = color_structure_dimension(&args)?;
                Some((index, args, dimension))
            })
            .collect::<Vec<_>>();
        if structures.len() < 3 {
            return None;
        }

        // Complete typed slots define the edges. Prefer a shortest cycle:
        // its trace decomposition introduces the fewest new invariant slots.
        // This also handles triangles without enumerating every vertex triple.
        let mut neighbours = vec![Vec::new(); structures.len()];
        for (left, (_, left_args, left_dimension)) in structures.iter().enumerate() {
            for (right, (_, right_args, right_dimension)) in
                structures.iter().enumerate().skip(left + 1)
            {
                if left_dimension == right_dimension
                    && common_structure_positions(left_args, right_args).count() == 1
                {
                    neighbours[left].push(right);
                    neighbours[right].push(left);
                }
            }
        }
        // Among equally short cycles prefer the one with the fewest ports into
        // generator lines: its cut attaches the fewest new legs to them.
        let line_slots = product
            .factors
            .iter()
            .filter_map(ProductFactor::generator_slots)
            .flatten()
            .collect::<Vec<_>>();
        let ports_into_lines = |cycle: &[usize]| {
            cycle
                .iter()
                .flat_map(|&vertex| &structures[vertex].1)
                .filter(|arg| line_slots.contains(&arg.as_view()))
                .count()
        };
        let mut cycle = Vec::new();
        let mut cycle_ports = 0;
        'roots: for root in 0..structures.len() {
            let mut parent = vec![usize::MAX; structures.len()];
            let mut depth = vec![0; structures.len()];
            parent[root] = root;
            let mut pending = VecDeque::from([root]);
            while let Some(vertex) = pending.pop_front() {
                for &next in &neighbours[vertex] {
                    if parent[next] == usize::MAX {
                        parent[next] = vertex;
                        depth[next] = depth[vertex] + 1;
                        pending.push_back(next);
                    } else if next != parent[vertex] && parent[next] != vertex {
                        let (mut left, mut right) = (vec![vertex], vec![next]);
                        let (mut a, mut b) = (vertex, next);
                        while a != b {
                            if depth[a] >= depth[b] {
                                a = parent[a];
                                left.push(a);
                            } else {
                                b = parent[b];
                                right.push(b);
                            }
                        }
                        right.pop();
                        left.extend(right.into_iter().rev());
                        if !cycle.is_empty() && left.len() > cycle.len() {
                            continue;
                        }
                        let ports = ports_into_lines(&left);
                        if cycle.is_empty() || left.len() < cycle.len() || ports < cycle_ports {
                            cycle = left;
                            cycle_ports = ports;
                            if cycle.len() == 3 && cycle_ports == 0 {
                                break 'roots;
                            }
                        }
                    }
                }
            }
        }
        if cycle.is_empty() || (cycle.len() > 3 && !self.settings.evaluate_traces) {
            return None;
        }
        // A cut of a cycle of length >= 5 is an expansion. While an ordered
        // generator line of the same term is still being decomposed, its
        // commutators keep attaching to the cycle's ports; cutting first
        // multiplies every later insertion by the cut. Cut after the line,
        // on the merged row. Shorter cycles stay local, as in FORM's repeat.
        if cycle.len() >= 5 && product.has_pending_line() {
            return None;
        }

        let links = (0..cycle.len())
            .map(|i| {
                common_structure_positions(
                    &structures[cycle[i]].1,
                    &structures[cycle[(i + 1) % cycle.len()]].1,
                )
                .next()
                .expect("adjacent cycle vertices share a port")
            })
            .collect::<Vec<_>>();
        let mut external = Vec::with_capacity(cycle.len());
        let mut prefactor = Atom::one();
        for i in 0..cycle.len() {
            let args = &structures[cycle[i]].1;
            let incoming = links[(i + cycle.len() - 1) % cycle.len()].right;
            let outgoing = links[i].left;
            if incoming == outgoing {
                return None;
            }
            let open = (0..3).find(|p| *p != incoming && *p != outgoing)?;
            external.push(args[open].clone());
            prefactor *= StructureView::orientation(incoming, outgoing);
        }
        if prefactor.is_zero() {
            return None;
        }
        let dimension = &structures[cycle[0]].2;
        let mut excluded = vec![false; product.len()];
        for &vertex in &cycle {
            excluded[structures[vertex].0] = true;
        }
        // Permuted copies of one long cycle in other terms decompose into
        // the same states, even when each cut is a one-shot decomposition.
        if cycle.len() >= 5 {
            self.certificate.merge.set(true);
        }
        let replacement = match external.as_slice() {
            [a, b, c] => {
                adjoint_casimir_for_dimension(dimension.clone()) / Atom::num(2) * color_f!(a, b, c)
            }
            [a, b, c, d] => {
                // Keep the compact color.h quartic identity. Reuse a removed
                // cycle edge for its new dummy; neither endpoint survives.
                let dummy = &structures[cycle[0]].1[links[0].left];
                let symmetric = trace_sym!(adjoint_rep(dimension.clone()); external.iter().map(|port| {
                    color_f!(Atom::var(T.chain_in), Atom::var(T.chain_out), port)
                }));
                symmetric
                    + adjoint_casimir_for_dimension(dimension.clone()) / Atom::num(6)
                        * (color_f!(a, b, dummy) * color_f!(d, c, dummy)
                            + color_f!(a, dummy, d) * color_f!(c, b, dummy))
            }
            // The shortest cycle's ports are its only connections to the
            // other factors; a colour-free remainder leaves the line isolated.
            _ => self.simplify_generator_trace(
                &adjoint_rep(dimension.clone()),
                &external,
                product.line_context(|index| excluded[index], whole_term),
            )?,
        };
        Some(product.excluding(&excluded, prefactor * replacement))
    }

    fn simplify_symmetric_invariant_product(product: &ProductView) -> Option<Atom> {
        for (left_index, left_factor) in product.factors.iter().enumerate() {
            let Some(left) = &left_factor.symmetric_invariant else {
                continue;
            };
            for (right_index, right_factor) in
                product.factors.iter().enumerate().skip(left_index + 1)
            {
                let Some(right) = &right_factor.symmetric_invariant else {
                    continue;
                };
                if left.args.len() < 3
                    || right.args.len() < 3
                    || !left.has_distinct_args()
                    || !right.has_distinct_args()
                {
                    continue;
                };

                let left_args = left
                    .args
                    .iter()
                    .map(|arg| arg.to_owned())
                    .collect::<Vec<_>>();
                let right_args = right
                    .args
                    .iter()
                    .map(|arg| arg.to_owned())
                    .collect::<Vec<_>>();
                let (common, left_open, right_open) =
                    symmetric_common_and_open_args(&left_args, &right_args);
                // An invariant with a single adjoint port is zero: the
                // adjoint of SU(N) has no invariant vector. This also avoids
                // leaving a label-order-dependent d_R(3)*d_S(4) remainder.
                if matches!(
                    (left_open.as_slice(), right_open.as_slice()),
                    ([], [_]) | ([_], [])
                ) {
                    return Some(Atom::Zero);
                }
                if left_args.len() != right_args.len() {
                    continue;
                }
                // Contract equal-rank symmetric traces into the corresponding
                // scalar invariant family, leaving a metric for one open pair.
                if common.len() == left_args.len() {
                    let replacement = Atom::num(left.phase * right.phase)
                        * color_symmetric_product(
                            left_args.len(),
                            left.rep.to_owned(),
                            right.rep.to_owned(),
                        );
                    return Some(product.replacing_pair(left_index, right_index, replacement));
                }
                let ([left_open], [right_open]) = (left_open.as_slice(), right_open.as_slice())
                else {
                    continue;
                };
                let dimension = color_adjoint_dimension(left_open)
                    .or_else(|| color_adjoint_dimension(right_open))
                    .or_else(|| common.iter().find_map(color_adjoint_dimension))?;
                let replacement = Atom::num(left.phase * right.phase)
                    * color_symmetric_product(
                        left_args.len(),
                        left.rep.to_owned(),
                        right.rep.to_owned(),
                    )
                    * color_metric(left_open.clone(), right_open.clone())
                    / dimension;
                return Some(product.replacing_pair(left_index, right_index, replacement));
            }
        }

        None
    }
}

/// The line of a product that absorbs its closest structure constant: the
/// pair's distance, the line and structure factor indices, the position of
/// the pair's first leg, and each line factor's generator slot.
struct AbsorbedLine<'a> {
    distance: usize,
    line: usize,
    structure: usize,
    first: usize,
    words: Vec<Option<AtomView<'a>>>,
}

#[derive(Clone, Debug)]
struct ChainView<'a> {
    start: AtomView<'a>,
    end: AtomView<'a>,
    factors: Vec<AtomView<'a>>,
}

impl<'a> ChainView<'a> {
    fn parse(expr: AtomView<'a>) -> Option<Self> {
        let AtomView::Fun(f) = expr else {
            return None;
        };
        if f.get_symbol() != T.chain {
            return None;
        }

        let args = f.iter().collect::<Vec<_>>();
        let [start, end, factors @ ..] = args.as_slice() else {
            return None;
        };
        Some(Self {
            start: *start,
            end: *end,
            factors: factors.to_vec(),
        })
    }
}

#[derive(Clone, Debug)]
struct TraceView<'a> {
    rep: AtomView<'a>,
    factors: Vec<AtomView<'a>>,
}

impl<'a> TraceView<'a> {
    fn parse(expr: AtomView<'a>) -> Option<Self> {
        let AtomView::Fun(f) = expr else {
            return None;
        };
        let (rep, factors) = shadowing::trace_parts(f)?;
        Some(Self { rep, factors })
    }
}

#[derive(Clone, Debug)]
struct StructureView<'a> {
    args: [AtomView<'a>; 3],
}

impl<'a> StructureView<'a> {
    fn contract_loop(&self, right: &Self) -> Option<Atom> {
        // Most pairs share fewer than two ports. Reject those on borrowed
        // indices before validating dimensions or constructing the Casimir.
        let mut common = common_structure_positions(&self.args, &right.args);
        let first = common.next()?;
        let second = common.next()?;
        let closed = common.next().is_some();
        let dimension = color_structure_dimension(&self.args)?;
        if color_structure_dimension(&right.args).as_ref() != Some(&dimension) {
            return None;
        }
        let prefactor = Self::orientation(first.left, second.left)
            * Self::orientation(first.right, second.right)
            * adjoint_casimir_for_dimension(dimension.clone());
        if closed {
            Some(prefactor * dimension)
        } else {
            let left_open = (0..3).find(|index| *index != first.left && *index != second.left)?;
            let right_open =
                (0..3).find(|index| *index != first.right && *index != second.right)?;
            Some(
                prefactor
                    * color_metric(
                        self.args[left_open].to_owned(),
                        right.args[right_open].to_owned(),
                    ),
            )
        }
    }

    fn orientation(first: usize, second: usize) -> Atom {
        // The third position is fixed. Read the antisymmetric permutation
        // directly; rebuilding f and extracting a coefficient must not depend
        // on the spelling or intrinsic normalization of admitted index labels.
        Atom::num(if second == (first + 1) % 3 { 1 } else { -1 })
    }

    fn parse(expr: AtomView<'a>) -> Option<Self> {
        let AtomView::Fun(f) = expr else {
            return None;
        };
        if f.get_symbol() != CS.f || f.get_nargs() != 3 {
            return None;
        }

        let mut args = f.iter();
        Some(Self {
            args: [
                args.next().unwrap(),
                args.next().unwrap(),
                args.next().unwrap(),
            ],
        })
    }
}

#[derive(Clone, Debug)]
struct SymmetricInvariantView<'a> {
    rep: AtomView<'a>,
    args: Vec<AtomView<'a>>,
    // Raw adjoint symmetric words carry i^rank relative to the Hermitian
    // invariant used by gram. Odd adjoint symmetric invariants are zero.
    phase: i8,
}

impl<'a> SymmetricInvariantView<'a> {
    fn has_distinct_args(&self) -> bool {
        // Repeated slots contract inside this invariant. They are not free
        // ports that a second copy may connect to form a Gram invariant.
        self.args
            .iter()
            .enumerate()
            .all(|(i, arg)| !self.args[..i].contains(arg))
    }

    fn parse(expr: AtomView<'a>) -> Option<Self> {
        if let Some(trace) = TraceView::parse(expr) {
            let [factor] = trace.factors.as_slice() else {
                return None;
            };
            let args = color_symmetric_trace_arg_views(trace.rep, *factor)?;
            ColorAlgebraSimplifier::trace_generator_is_adjoint(&trace.rep.to_owned(), &args)?;
            return Some(Self {
                rep: trace.rep,
                phase: symmetric_trace_phase(trace.rep, args.len()),
                args,
            });
        }

        let AtomView::Fun(f) = expr else {
            return None;
        };
        if f.get_symbol() != CS.d || f.get_nargs() < 4 {
            return None;
        }

        let args = f.iter().collect::<Vec<_>>();
        // Explicit d_R also supports representations without a generator-word
        // kernel. Its contracted index space must still be homogeneous.
        let _ = color_structure_dimension(&args[1..])?;
        Some(Self {
            rep: args[0],
            args: args[1..].to_vec(),
            phase: 1,
        })
    }
}

#[derive(Clone, Debug)]
struct ProductFactor<'a> {
    atom: AtomView<'a>,
    chain: Option<ChainView<'a>>,
    trace: Option<TraceView<'a>>,
    structure: Option<StructureView<'a>>,
    symmetric_invariant: Option<SymmetricInvariantView<'a>>,
}

impl<'a> ProductFactor<'a> {
    fn parse(atom: AtomView<'a>) -> Self {
        let chain = ChainView::parse(atom);
        let trace = TraceView::parse(atom);
        let structure = StructureView::parse(atom);
        let symmetric_invariant = SymmetricInvariantView::parse(atom);
        Self {
            atom,
            chain,
            trace,
            structure,
            symmetric_invariant,
        }
    }

    /// Lines before other factors, ordered by their numbers of external
    /// (once-occurring) and of all generator slots.
    fn line_cost(&self) -> (bool, usize, usize) {
        let Some(slots) = self.generator_slots() else {
            return (true, 0, 0);
        };
        let external = slots
            .iter()
            .filter(|slot| slots.iter().filter(|other| other == slot).count() == 1)
            .count();
        (false, external, slots.len())
    }

    /// Borrow all generator side ports before selecting a connection to another
    /// factor. A repeated port contracts within this line and is not an external
    /// leg, including when its other occurrence lies in a symmetric block.
    fn generator_slots(&self) -> Option<Vec<AtomView<'a>>> {
        let (representation, factors) = if let Some(trace) = &self.trace {
            (trace.rep, &trace.factors)
        } else if let Some(chain) = &self.chain {
            let _ = fundamental_chain_dimension_view(chain.start, chain.end)?;
            (chain.start, &chain.factors)
        } else {
            return None;
        };
        let adjoint =
            matches!(representation, AtomView::Fun(f) if f.get_symbol() == CS.adjoint_rep);
        let mut slots = Vec::new();
        for &factor in factors {
            let slot = if adjoint {
                adjoint_generator_slot(representation, factor)
            } else {
                color_generator_adjoint_view(factor)
            };
            if let Some(slot) = slot {
                slots.push(slot);
            } else {
                slots.extend(color_symmetric_trace_arg_views(representation, factor)?);
            }
        }
        Some(slots)
    }
}

#[derive(Clone, Debug)]
struct ProductView<'a> {
    factors: Vec<ProductFactor<'a>>,
}

impl<'a> ProductView<'a> {
    fn parse(expr: AtomView<'a>) -> Self {
        Self {
            factors: multiplicative_factor_views(expr)
                .into_iter()
                .map(ProductFactor::parse)
                .collect(),
        }
    }

    fn len(&self) -> usize {
        self.factors.len()
    }

    fn replacing_pair(&self, left_index: usize, right_index: usize, replacement: Atom) -> Atom {
        let kept = || {
            self.factors
                .iter()
                .enumerate()
                .filter(move |(index, _)| *index != left_index && *index != right_index)
                .map(|(_, factor)| factor.atom)
        };
        Self::eliminating_emitted_metrics(&replacement, kept())
            .unwrap_or_else(|| Atom::mul_many(std::iter::once(replacement.as_view()).chain(kept())))
    }

    fn replacing_one(&self, target_index: usize, replacement: Atom) -> Atom {
        let kept = self
            .factors
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != target_index)
            .map(|(_, factor)| factor.atom);
        Self::eliminating_emitted_metrics(&replacement, kept).unwrap_or_else(|| {
            Atom::mul_many(self.factors.iter().enumerate().map(|(index, factor)| {
                if index == target_index {
                    replacement.as_view()
                } else {
                    factor.atom
                }
            }))
        })
    }

    /// Multiply the kept factors by an identity's replacement, eliminating
    /// each adjoint metric the replacement emits as a factor as FORM does
    /// with d_: g(x,y) X(..y..) = X(..x..) when y occurs exactly once among
    /// the other factors, uniformly in every summand of a sum, and x at most
    /// once. A metric with both ends open, or whose other end repeats, sits
    /// in a power or in a foreign tensor, or is spelled by a summand-local
    /// dummy, is left to the shared contractor. `None` when the
    /// replacement emits no such metric.
    fn eliminating_emitted_metrics<'b>(
        replacement: &Atom,
        kept: impl Iterator<Item = AtomView<'b>>,
    ) -> Option<Atom> {
        let emitted = multiplicative_factor_views(replacement.as_view());
        if !emitted
            .iter()
            .any(|factor| emitted_adjoint_metric(*factor).is_some())
        {
            return None;
        }
        let mut metrics = Vec::new();
        let mut factors = Vec::new();
        for factor in emitted {
            match emitted_adjoint_metric(factor) {
                Some(ends) => metrics.push(ends),
                None => factors.push(factor.to_owned()),
            }
        }
        factors.extend(kept.map(|factor| factor.to_owned()));
        let mut remaining = Vec::new();
        for [x, y] in metrics {
            let eliminated = [(&y, &x), (&x, &y)].into_iter().any(|(from, to)| {
                let Some(counts) = factors
                    .iter()
                    .map(|factor| colour_slot_occurrences(factor.as_view(), from.as_view()))
                    .collect::<Option<Vec<_>>>()
                else {
                    return false;
                };
                if counts.iter().sum::<usize>() != 1 {
                    return false;
                }
                // The kept end may occur at most once more, uniformly: a
                // summand-local dummy spelled like it would be captured.
                if factors
                    .iter()
                    .map(|factor| colour_slot_occurrences(factor.as_view(), to.as_view()))
                    .sum::<Option<usize>>()
                    .is_none_or(|count| count > 1)
                {
                    return false;
                }
                let position = counts.iter().position(|&count| count == 1).unwrap();
                factors[position] = factors[position].replace_map(|arg, _context, out| {
                    if arg == from.as_view() {
                        **out = to.clone();
                    }
                });
                true
            });
            if !eliminated {
                remaining.push(color_metric(x, y));
            }
        }
        Some(Atom::mul_many(factors.iter().chain(&remaining)))
    }

    /// A chain is maximal when no other chain or metric of the product carries
    /// one of its fundamental endpoints, and its endpoints do not close it.
    fn chain_is_maximal(&self, index: usize, chain: &ChainView<'_>) -> bool {
        let endpoint = |slot: AtomView<'_>| {
            color_fundamental_slot(slot).map(|(dimension, index, _)| (dimension, index))
        };
        let (Some(start), Some(end)) = (endpoint(chain.start), endpoint(chain.end)) else {
            return false;
        };
        if start == end {
            return false;
        }
        self.factors.iter().enumerate().all(|(other, factor)| {
            if other == index {
                return true;
            }
            let ports = if let Some(chain) = &factor.chain {
                vec![chain.start, chain.end]
            } else if let AtomView::Fun(metric) = factor.atom
                && metric.get_symbol() == ETS.metric
            {
                metric.iter().collect()
            } else {
                return true;
            };
            ports
                .into_iter()
                .filter_map(endpoint)
                .all(|port| port != start && port != end)
        })
    }

    /// Whether a generator trace of this product still awaits decomposition:
    /// an ordered word or a symmetric prefix with an ordered suffix. Open
    /// chains are terminal and a lone symmetric block is an invariant.
    fn has_pending_line(&self) -> bool {
        self.factors.iter().any(|factor| {
            factor
                .trace
                .as_ref()
                .is_some_and(|trace| trace.factors.len() >= 2)
                && factor.generator_slots().is_some()
        })
    }

    /// A line is isolated when nothing else in a whole top-level term has a
    /// colour port or pending colour work; no other factor can reach it.
    fn line_context(&self, line: impl Fn(usize) -> bool, whole_term: bool) -> LineContext {
        if whole_term
            && self.factors.iter().enumerate().all(|(index, factor)| {
                line(index) || !atom_contains_color_port_or_work(factor.atom)
            })
        {
            LineContext::Isolated
        } else {
            LineContext::Unknown
        }
    }

    fn excluding(&self, excluded: &[bool], replacement: Atom) -> Atom {
        let kept = || {
            self.factors
                .iter()
                .enumerate()
                .filter(|(index, _)| !excluded[*index])
                .map(|(_, factor)| factor.atom)
        };
        Self::eliminating_emitted_metrics(&replacement, kept())
            .unwrap_or_else(|| Atom::mul_many(kept()) * replacement)
    }

    /// The colour sum factor with the fewest summands that may be
    /// distributed. `selective` replaces the f-only connection test below by
    /// [`Self::rule_can_span_sum`] for sums with colour ports; port-free
    /// invariant sums keep the existing distribution and with it the
    /// established expanded output.
    fn color_sum_factor(&self, selective: bool) -> Option<(usize, AddView<'a>)> {
        let mut sums = self
            .factors
            .iter()
            .enumerate()
            .filter(|(_, factor)| matches!(factor.atom, AtomView::Add(_)))
            .collect::<Vec<_>>();
        sums.sort_by_key(|(_, factor)| match factor.atom {
            AtomView::Add(add) => add.get_nargs(),
            _ => 0,
        });
        sums.into_iter()
            .find_map(|(index, factor)| match factor.atom {
                AtomView::Add(add) if add.iter().any(atom_contains_color_node) => {
                    if selective && add.iter().any(atom_contains_color_port_or_work) {
                        return self.rule_can_span_sum(index).then_some((index, add));
                    }
                    let mut matcher = SlotMatcher::default();
                    let mut outside = Vec::new();
                    let mut only_structures = true;
                    for (other_index, other) in self.factors.iter().enumerate() {
                        if other_index == index {
                            continue;
                        }
                        if let Some(structure) = &other.structure {
                            let parsed = structure
                                .args
                                .iter()
                                .map(|slot| matcher.parse::<LibraryRep, AbstractIndex>(*slot))
                                .collect::<Result<Vec<_>, _>>();
                            let Ok(parsed) = parsed else {
                                only_structures = false;
                                break;
                            };
                            outside.extend(parsed);
                        } else if !OrderedStructure::<LibraryRep, AbstractIndex>::syntactic_structure_from_atom(
                            other.atom,
                            &mut matcher,
                        )
                        .is_ok_and(|structure| structure.canonical().is_scalar())
                        {
                            only_structures = false;
                            break;
                        }
                    }
                    // An f-only region needs at least two connections into
                    // this sum for a loop or a trace/symmetry contraction.
                    // A single bridge cannot authorize distribution. Reuse
                    // admitted homogeneous-sum inference: it inspects only
                    // the first summand and performs no validation walk.
                    // Other tensors, missing interfaces and Fierz/chain
                    // connections retain the existing conservative behavior.
                    if only_structures && outside.is_empty() {
                        return None;
                    }
                    if only_structures
                        && let Ok(structure) = OrderedStructure::<LibraryRep, AbstractIndex>::syntactic_structure_from_atom(
                            factor.atom,
                            &mut matcher,
                        )
                        && outside
                            .iter()
                            .filter(|outer| structure.canonical().external_structure_iter().any(|inner| outer.matches(&inner)))
                            .take(2)
                            .count()
                            < 2
                    {
                        return None;
                    }
                    Some((index, add))
                }
                _ => None,
            })
    }

    /// Every product rule acts on contracted colour ports, so a rule spanning
    /// a summand and other factors uses the sum's external slots. Fundamental
    /// lines can join or Fierz-contract through one such connection. Adjoint
    /// structures, including the adjoint chains and traces written by the
    /// shared contractor, need a loop through the sum and hence two. With
    /// fewer, distribution only spreads the other factors over the summands.
    /// An interface that cannot be inferred keeps the sum distributable.
    fn rule_can_span_sum(&self, sum_index: usize) -> bool {
        let sum = SlotCounts::of(self.factors[sum_index].atom);
        if sum.opaque {
            return true;
        }
        let mut connections = 0;
        for (index, factor) in self.factors.iter().enumerate() {
            if index == sum_index {
                continue;
            }
            let counts = SlotCounts::shared_with(factor.atom, &sum);
            if counts.opaque {
                return true;
            }
            let shared = counts.ports().filter(|key| sum.get(key) == Some(1)).count();
            if shared == 0 {
                continue;
            }
            if atom_contains_fundamental_line(factor.atom) {
                return true;
            }
            connections += shared;
        }
        connections >= 2
    }
}

fn chain_parts(chain: AtomView) -> Option<(Atom, Atom, Vec<Atom>)> {
    let AtomView::Fun(f) = chain else {
        return None;
    };
    if f.get_symbol() != T.chain {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    let [start, end, factors @ ..] = args.as_slice() else {
        return None;
    };
    Some((start.clone(), end.clone(), factors.to_vec()))
}

fn trace_parts(trace: AtomView) -> Option<(Atom, Vec<Atom>)> {
    let AtomView::Fun(f) = trace else {
        return None;
    };
    let (rep, factors) = shadowing::trace_parts(f)?;
    Some((
        rep.to_owned(),
        factors
            .into_iter()
            .map(|factor| factor.to_owned())
            .collect(),
    ))
}

fn color_generator_adjoint(factor: AtomView) -> Option<Atom> {
    color_generator_adjoint_view(factor).map(|arg| arg.to_owned())
}

pub(super) fn color_generator_adjoint_view(factor: AtomView<'_>) -> Option<AtomView<'_>> {
    let AtomView::Fun(f) = factor else {
        return None;
    };
    if f.get_symbol() != CS.t || f.get_nargs() != 3 {
        return None;
    }

    let args = f.iter().collect::<Vec<_>>();
    if !has_chain_endpoints(args[1], args[2]) {
        return None;
    }

    Some(args[0])
}

fn structure_constant_args(factor: AtomView) -> Option<[Atom; 3]> {
    let AtomView::Fun(f) = factor else {
        return None;
    };
    if f.get_symbol() != CS.f || f.get_nargs() != 3 {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    Some([args[0].clone(), args[1].clone(), args[2].clone()])
}

fn symmetric_trace_phase(rep: AtomView, rank: usize) -> i8 {
    if matches!(rep, AtomView::Fun(f) if f.get_symbol() == CS.adjoint_rep) {
        match rank % 4 {
            0 => 1,
            2 => -1,
            _ => 0,
        }
    } else {
        1
    }
}

fn color_symmetric_trace_arg_views<'a>(
    rep: AtomView<'a>,
    projector: AtomView<'a>,
) -> Option<Vec<AtomView<'a>>> {
    let AtomView::Fun(f) = projector else {
        return None;
    };
    if f.get_symbol() != *shadowing::SYM {
        return None;
    }

    let args = if let AtomView::Fun(representation) = rep
        && representation.get_symbol() == CS.adjoint_rep
        && representation.get_nargs() == 1
    {
        f.iter()
            .map(|factor| adjoint_generator_slot(rep, factor))
            .collect::<Option<Vec<_>>>()?
    } else {
        f.iter()
            .map(color_generator_adjoint_view)
            .collect::<Option<Vec<_>>>()?
    };
    if !args.is_empty() {
        let _ = color_structure_dimension(&args)?;
    }
    Some(args)
}

fn adjoint_generator_slot<'a>(rep: AtomView<'a>, factor: AtomView<'a>) -> Option<AtomView<'a>> {
    let (slot, coefficient) = adjoint_generator(rep, factor)?;
    coefficient.is_one().then_some(slot)
}

fn adjoint_generator<'a>(rep: AtomView<'a>, factor: AtomView<'a>) -> Option<(AtomView<'a>, Atom)> {
    let AtomView::Fun(representation) = rep else {
        return None;
    };
    let mut coefficient = Atom::one();
    let mut matrix = None;
    for factor in multiplicative_factor_views(factor) {
        if matches!(factor, AtomView::Num(_)) {
            coefficient *= factor;
        } else if matrix.replace(factor).is_some() {
            return None;
        }
    }
    let AtomView::Fun(structure) = matrix? else {
        return None;
    };
    if representation.get_symbol() != CS.adjoint_rep
        || representation.get_nargs() != 1
        || structure.get_symbol() != CS.f
        || structure.get_nargs() != 3
    {
        return None;
    }
    let args = structure.iter().collect::<Vec<_>>();
    let input = args
        .iter()
        .position(|arg| is_chain_endpoint(*arg, T.chain_in))?;
    let output = args
        .iter()
        .position(|arg| is_chain_endpoint(*arg, T.chain_out))?;
    let slot = args[(0..3).find(|&i| i != input && i != output)?];
    let AtomView::Fun(index) = slot else {
        return None;
    };
    if index.get_symbol() == CS.adjoint_rep
        && index.get_nargs() == 2
        && index.iter().next()? == representation.iter().next()?
    {
        return Some((
            slot,
            coefficient * StructureView::orientation(input, output),
        ));
    }
    None
}

fn color_antisymmetric_generator_args(factor: AtomView) -> Option<(Atom, Vec<Atom>)> {
    let (prefactor, projector_symbol, factors) = projector_factor(factor)?;
    if projector_symbol != *shadowing::ANTISYM {
        return None;
    }

    let args = factors
        .iter()
        .map(|factor| color_generator_adjoint(factor.as_view()))
        .collect::<Option<Vec<_>>>()?;
    Some((prefactor, args))
}

fn projector_factor(factor: AtomView) -> Option<(Atom, Symbol, Vec<Atom>)> {
    if let Some((projector_symbol, factors)) = projector_parts(factor) {
        return Some((Atom::num(1), projector_symbol, factors));
    }

    let AtomView::Mul(product) = factor else {
        return None;
    };

    let mut prefactor = Atom::num(1);
    let mut projector = None;
    for factor in product.iter() {
        if let Some(parts) = projector_parts(factor) {
            if projector.is_some() {
                return None;
            }
            projector = Some(parts);
        } else if matches!(factor, AtomView::Num(_)) {
            prefactor *= factor.to_owned();
        } else {
            return None;
        }
    }

    let (projector_symbol, factors) = projector?;
    Some((prefactor, projector_symbol, factors))
}

fn projector_parts(projector: AtomView) -> Option<(Symbol, Vec<Atom>)> {
    let AtomView::Fun(f) = projector else {
        return None;
    };
    if f.get_symbol() != *shadowing::SYM && f.get_symbol() != *shadowing::ANTISYM {
        return None;
    }

    Some((f.get_symbol(), f.iter().map(|arg| arg.to_owned()).collect()))
}

fn color_fundamental_slot(slot: AtomView) -> Option<(Atom, Atom, bool)> {
    if let Some((dimension, index)) = representation_slot(slot, CS.fundamental_rep) {
        return Some((dimension, index, false));
    }

    let AtomView::Fun(f) = slot else {
        return None;
    };
    if f.get_symbol() != AIND_SYMBOLS.dind || f.get_nargs() != 1 {
        return None;
    }

    representation_slot(f.iter().next()?, CS.fundamental_rep)
        .map(|(dimension, index)| (dimension, index, true))
}

fn color_adjoint_dimension(slot: &impl AtomCore) -> Option<Atom> {
    representation_slot(slot.as_atom_view(), CS.adjoint_rep).map(|(dimension, _)| dimension)
}

fn color_structure_dimension(args: &[impl AtomCore]) -> Option<Atom> {
    let mut dimensions = args.iter().map(color_adjoint_dimension);
    let dimension = dimensions.next()??;
    dimensions
        .all(|candidate| candidate.is_some_and(|candidate| candidate == dimension))
        .then_some(dimension)
}

/// The dimension and index of a `symbol(dimension, index)` slot.
fn representation_slot_view(
    slot: AtomView<'_>,
    symbol: Symbol,
) -> Option<(AtomView<'_>, AtomView<'_>)> {
    let AtomView::Fun(f) = slot else {
        return None;
    };
    if f.get_symbol() != symbol || f.get_nargs() != 2 {
        return None;
    }
    let mut args = f.iter();
    Some((args.next()?, args.next()?))
}

fn representation_slot(slot: AtomView, symbol: Symbol) -> Option<(Atom, Atom)> {
    representation_slot_view(slot, symbol)
        .map(|(dimension, index)| (dimension.to_owned(), index.to_owned()))
}

fn trace_terminal_dimension(rep: AtomView) -> Option<Atom> {
    let trace = trace!(rep.to_owned());
    let simplified = trace.replace_multiple_repeat(TRACE_TERMINALS.as_ref());
    (simplified != trace).then_some(simplified)
}

fn has_chain_endpoints(left: AtomView, right: AtomView) -> bool {
    is_chain_endpoint(left, T.chain_in) && is_chain_endpoint(right, T.chain_out)
}

fn is_chain_identity_factor(factor: AtomView) -> bool {
    // `chain!` represents an empty line as a metric over the placeholder
    // endpoints; collapse it to the physical endpoint metric/dimension.
    let AtomView::Fun(f) = factor else {
        return false;
    };
    if f.get_symbol() != ETS.metric || f.get_nargs() != 2 {
        return false;
    }

    let args = f.iter().collect::<Vec<_>>();
    has_chain_endpoints(args[0], args[1])
}

fn is_chain_endpoint(arg: AtomView, expected: Symbol) -> bool {
    matches!(arg, AtomView::Var(symbol) if symbol.get_symbol() == expected)
}

fn color_metric(left: Atom, right: Atom) -> Atom {
    ETS.metric_literal(left, right)
}

/// The two distinct ends of an adjoint metric factor.
fn emitted_adjoint_metric(factor: AtomView<'_>) -> Option<[Atom; 2]> {
    let AtomView::Fun(metric) = factor else {
        return None;
    };
    if metric.get_symbol() != ETS.metric || metric.get_nargs() != 2 {
        return None;
    }
    let [x, y] = [0, 1].map(|i| metric.iter().nth(i).unwrap().to_owned());
    (x != y
        && color_adjoint_dimension(&x).is_some()
        && color_adjoint_dimension(&x) == color_adjoint_dimension(&y))
    .then_some([x, y])
}

/// Occurrences of an adjoint slot in a colour factor, the same in every
/// summand of a sum. `None` where renaming the occurrence is not a
/// contraction of one copy: in a power, in a foreign function or in summands
/// with different counts.
fn colour_slot_occurrences(factor: AtomView<'_>, slot: AtomView<'_>) -> Option<usize> {
    if factor == slot {
        return Some(1);
    }
    match factor {
        AtomView::Add(sum) => {
            let mut counts = sum.iter().map(|term| colour_slot_occurrences(term, slot));
            let first = counts.next()??;
            counts.all(|count| count == Some(first)).then_some(first)
        }
        AtomView::Mul(product) => product
            .iter()
            .map(|factor| colour_slot_occurrences(factor, slot))
            .sum(),
        AtomView::Fun(_) if representation_slot(factor, CS.adjoint_rep).is_some() => Some(0),
        AtomView::Fun(function)
            if [
                CS.f,
                CS.t,
                CS.d,
                T.chain,
                T.trace,
                ETS.metric,
                *shadowing::SYM,
                *shadowing::ANTISYM,
                *shadowing::CYCLIC,
            ]
            .contains(&function.get_symbol()) =>
        {
            function
                .iter()
                .map(|arg| colour_slot_occurrences(arg, slot))
                .sum()
        }
        AtomView::Fun(_) | AtomView::Pow(_) => (!factor.contains(slot)).then_some(0),
        _ => Some(0),
    }
}

fn quadratic_casimir(rep: Atom) -> Atom {
    CS.cas(Atom::num(2), rep)
}

fn quadratic_index(rep: Atom) -> Atom {
    CS.idx(Atom::num(2), rep)
}

fn fundamental_rep(dimension: Atom) -> Atom {
    ColorFundamental {}.to_symbolic([dimension])
}

fn adjoint_rep(dimension: Atom) -> Atom {
    ColorAdjoint {}.to_symbolic([dimension])
}

fn fundamental_casimir(dimension: Atom) -> Atom {
    quadratic_casimir(fundamental_rep(dimension))
}

fn fundamental_index(dimension: Atom) -> Atom {
    quadratic_index(fundamental_rep(dimension))
}

fn adjoint_casimir_for_dimension(dimension: Atom) -> Atom {
    quadratic_casimir(adjoint_rep(dimension))
}

fn color_symmetric_trace(rep: &Atom, factors: impl IntoIterator<Item = Atom>) -> Atom {
    let sym_factors = factors
        .into_iter()
        .map(|factor| color_t!(factor.clone()))
        .collect::<Vec<_>>();
    trace_sym!(rep.clone(); sym_factors)
}

fn color_symmetric_product(rank: usize, left_rep: Atom, right_rep: Atom) -> Atom {
    CS.gram(Atom::num(rank as i64), left_rep, right_rep)
}

fn fundamental_chain_dimension(start: &Atom, end: &Atom) -> Option<Atom> {
    fundamental_chain_dimension_view(start.as_view(), end.as_view())
}

pub(super) fn fundamental_chain_dimension_view(
    start: AtomView<'_>,
    end: AtomView<'_>,
) -> Option<Atom> {
    let Some((start_dimension, _, false)) = color_fundamental_slot(start) else {
        return None;
    };
    let Some((end_dimension, _, true)) = color_fundamental_slot(end) else {
        return None;
    };
    (start_dimension == end_dimension).then_some(start_dimension)
}

fn is_su_n_fundamental_chain(start: &Atom, end: &Atom) -> bool {
    fundamental_chain_dimension(start, end).is_some_and(|dimension| dimension == Atom::var(CS.nc))
}

fn is_su_n_adjoint_slot(slot: &Atom) -> bool {
    let n = Atom::var(CS.nc);
    let su_n_adjoint_dimension = n.clone().pow(Atom::num(2)) - Atom::num(1);
    color_adjoint_dimension(slot).is_some_and(|dimension| dimension == su_n_adjoint_dimension)
}

fn restore_explicit_su_n_generator_chains(expression: Atom) -> Atom {
    expression.replace_map(|arg, _context, out| {
        let Some((start, end, factors)) = chain_parts(arg) else {
            return;
        };
        let [factor] = factors.as_slice() else {
            return;
        };
        if !is_su_n_fundamental_chain(&start, &end) {
            return;
        }
        let Some(adjoint) = color_generator_adjoint(factor.as_view()) else {
            return;
        };
        if !is_su_n_adjoint_slot(&adjoint) {
            return;
        }

        **out = CS.explicit_t(adjoint, start, end);
    })
}

fn symmetric_common_and_open_args(
    left: &[Atom],
    right: &[Atom],
) -> (Vec<Atom>, Vec<Atom>, Vec<Atom>) {
    let mut right_used = vec![false; right.len()];
    let mut common = Vec::new();
    let mut left_open = Vec::new();

    for left_arg in left {
        if let Some((right_index, _)) = right
            .iter()
            .enumerate()
            .find(|(index, right_arg)| !right_used[*index] && *right_arg == left_arg)
        {
            right_used[right_index] = true;
            common.push(left_arg.clone());
        } else {
            left_open.push(left_arg.clone());
        }
    }

    let right_open = right
        .iter()
        .enumerate()
        .filter(|(index, _)| !right_used[*index])
        .map(|(_, right_arg)| right_arg.clone())
        .collect();

    (common, left_open, right_open)
}

fn chain_with_removed_range(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    from: usize,
    to: usize,
) -> Atom {
    let remaining = factors_excluding_range(factors, from, to);
    if remaining.is_empty() {
        color_metric(start.clone(), end.clone())
    } else {
        chain!(start.clone(), end.clone(); remaining)
    }
}

fn chain_with_factors(start: Atom, end: Atom, factors: Vec<Atom>) -> Atom {
    if factors.is_empty() {
        color_metric(start, end)
    } else {
        chain!(start, end; factors)
    }
}

fn chain_with_factor_view_slices(
    start: AtomView<'_>,
    end: AtomView<'_>,
    slices: &[&[AtomView<'_>]],
) -> Atom {
    let factors = slices
        .iter()
        .flat_map(|slice| slice.iter().map(|factor| factor.to_owned()))
        .collect::<Vec<_>>();
    chain_with_factors(start.to_owned(), end.to_owned(), factors)
}

fn trace_with_factors(rep: Atom, factors: Vec<Atom>) -> Atom {
    if factors.is_empty() {
        trace_terminal_dimension(rep.as_view()).unwrap_or(rep)
    } else {
        trace!(rep; factors)
    }
}

fn factors_excluding_range(factors: &[Atom], from: usize, to: usize) -> Vec<Atom> {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| *index < from || *index >= to)
        .map(|(_, factor)| factor.clone())
        .collect()
}

fn factors_excluding_indices(factors: &[Atom], excluded: &[usize]) -> Vec<Atom> {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| !excluded.contains(index))
        .map(|(_, factor)| factor.clone())
        .collect()
}

fn multiplicative_factor_views(expr: AtomView<'_>) -> Vec<AtomView<'_>> {
    match expr {
        AtomView::Mul(mul) => mul.iter().collect(),
        _ => vec![expr],
    }
}

fn atom_contains_color_node(expr: AtomView<'_>) -> bool {
    color_node_search(expr, false)
}

/// Scalar invariants such as Casimirs have no ports and no rule of their own,
/// so they leave every other colour factor of a term isolated.
fn atom_contains_color_port_or_work(expr: AtomView<'_>) -> bool {
    color_node_search(expr, true)
}

/// Whether the expression holds an object of the colour algebra itself:
/// a generator, a structure constant, a symmetric or scalar invariant, a
/// line or trace, or a metric, such as the kernel emits. A foreign tensor
/// with colour slots holds none.
fn atom_contains_colour_algebra(expr: AtomView<'_>) -> bool {
    let heads = [
        CS.f, CS.d, CS.t, CS.gram, CS.cas, CS.idx, T.chain, T.trace, ETS.metric,
    ];
    let mut found = false;
    expr.visitor(&mut |node| {
        found |= matches!(node, AtomView::Fun(function) if heads.contains(&function.get_symbol()));
        !found
    });
    found
}

/// The factors with colour slots but no colour algebra, such as a foreign
/// tensor V(coad(8,z)), that share no slot, directly or through other
/// factors, with a factor holding colour algebra: in the product, in a term
/// of a sum, or in a term of a colour sum factor whose other factors they do
/// not meet either. No colour identity reaches them, so they are not part of
/// any colour row: they stay factors as found, like scalar spectators.
/// Unreadable slots connect everything.
fn foreign_colour_factors(expression: AtomView<'_>) -> Vec<AtomView<'_>> {
    let mut search = ForeignColourSearch::default();
    search.visit(expression, &[]);
    // A factor spelled alike elsewhere must be unreachable there too.
    search
        .unreached
        .retain(|factor| !search.reached.contains(factor));
    search.unreached
}

#[derive(Default)]
struct ForeignColourSearch<'a> {
    unreached: Vec<AtomView<'a>>,
    reached: Vec<AtomView<'a>>,
}

impl<'a> ForeignColourSearch<'a> {
    /// `outside` holds the factors around `expression`: a factor meeting
    /// them is reached through them.
    fn visit(&mut self, expression: AtomView<'a>, outside: &[AtomView<'a>]) {
        if let AtomView::Add(sum) = expression {
            for summand in sum.iter() {
                self.visit(summand, outside);
            }
            return;
        }
        let factors = multiplicative_factor_views(expression);
        let foreign = factors
            .iter()
            .map(|&factor| {
                atom_contains_color_node(factor) && !atom_contains_colour_algebra(factor)
            })
            .collect::<Vec<_>>();
        if foreign.contains(&true) {
            let slots = factors
                .iter()
                .map(|&factor| SlotCounts::of(factor))
                .collect::<Vec<_>>();
            let around = SlotCounts::of_factors(outside.iter().copied());
            let opaque = around.opaque || slots.iter().any(|counts| counts.opaque);
            let mut reached = (0..factors.len())
                .map(|index| {
                    opaque
                        || (!foreign[index] && atom_contains_color_node(factors[index]))
                        || slots[index].meets(&around)
                })
                .collect::<Vec<_>>();
            let mut frontier = (0..factors.len())
                .filter(|&index| reached[index])
                .collect::<Vec<_>>();
            while let Some(index) = frontier.pop() {
                for other in 0..factors.len() {
                    if !reached[other] && slots[index].meets(&slots[other]) {
                        reached[other] = true;
                        frontier.push(other);
                    }
                }
            }
            for (index, &factor) in factors.iter().enumerate() {
                if foreign[index] {
                    if reached[index] {
                        self.reached.push(factor);
                    } else {
                        self.unreached.push(factor);
                    }
                }
            }
        }
        // The terms of a colour sum are rows of their own.
        for (index, &factor) in factors.iter().enumerate() {
            if matches!(factor, AtomView::Add(_)) && atom_contains_colour_algebra(factor) {
                let around = outside
                    .iter()
                    .copied()
                    .chain(
                        factors
                            .iter()
                            .enumerate()
                            .filter(|&(other, _)| other != index)
                            .map(|(_, &other)| other),
                    )
                    .collect::<Vec<_>>();
                self.visit(factor, &around);
            }
        }
    }
}

/// The chains and closed traces that `join_chains` collects from a product of
/// fundamental line segments, explicit generators t(a, cof(d,i), dind(cof(d,j)))
/// and chains from cof(d,i) to dind(cof(d,j)), found by following their
/// indices instead of matching pairs of factors. Metrics are not line
/// segments there either. `None` when another factor has an indexed
/// fundamental slot, or a slot does not join exactly one end to one start:
/// the shared collector handles those.
fn joined_fundamental_lines(expression: AtomView<'_>) -> Option<Atom> {
    Some(joined_fundamental_line_terms(expression)?.normalize_chains())
}

/// The joins of [`joined_fundamental_lines`] in every term of a sum, before
/// the chain normalization.
fn joined_fundamental_line_terms(expression: AtomView<'_>) -> Option<Atom> {
    if let AtomView::Add(sum) = expression {
        return sum
            .iter()
            .map(|term| {
                if has_fundamental_slot(term) {
                    joined_fundamental_line_terms(term)
                } else {
                    Some(term.to_owned())
                }
            })
            .collect::<Option<Vec<_>>>()
            .map(Atom::add_many);
    }
    let mut segments = Vec::new();
    let mut others = Vec::new();
    for factor in multiplicative_factor_views(expression) {
        if let Some(segment) = fundamental_line_segment(factor) {
            segments.push(segment);
        } else if !has_fundamental_slot(factor)
            || matches!(factor, AtomView::Fun(metric) if metric.get_symbol() == ETS.metric)
        {
            others.push(factor.to_owned());
        } else if matches!(factor, AtomView::Add(_)) {
            // A sum's lines join within its terms only.
            others.push(joined_fundamental_line_terms(factor)?);
        } else {
            return None;
        }
    }
    if segments.is_empty() {
        return Some(Atom::mul_many(&others));
    }
    let mut successor = Vec::with_capacity(segments.len());
    let mut has_predecessor = vec![false; segments.len()];
    for (_, end, _) in &segments {
        let mut next = segments
            .iter()
            .enumerate()
            .filter(|(_, (start, ..))| start == end)
            .map(|(index, _)| index);
        let first = next.next();
        if next.next().is_some() {
            return None;
        }
        if let Some(first) = first {
            if has_predecessor[first] {
                return None;
            }
            has_predecessor[first] = true;
        }
        successor.push(first);
    }
    let mut used = vec![false; segments.len()];
    let mut lines = Vec::new();
    // Open lines start where no segment ends.
    for first in (0..segments.len()).filter(|&index| !has_predecessor[index]) {
        let mut chain = FunctionBuilder::new(T.chain).add_arg(segments[first].0);
        let mut word = Vec::new();
        let mut last = first;
        let mut current = Some(first);
        while let Some(index) = current {
            used[index] = true;
            word.extend(segments[index].2.iter().cloned());
            last = index;
            current = successor[index];
        }
        chain = chain.add_arg(
            FunctionBuilder::new(AIND_SYMBOLS.dind)
                .add_arg(segments[last].1)
                .finish(),
        );
        for factor in word {
            chain = chain.add_arg(factor);
        }
        lines.push(chain.finish());
    }
    // The remaining segments close into traces.
    for first in 0..segments.len() {
        if used[first] {
            continue;
        }
        let mut word = Vec::new();
        let mut current = first;
        while !used[current] {
            used[current] = true;
            word.extend(segments[current].2.iter().cloned());
            current = successor[current]?;
        }
        let (dimension, _) = representation_slot_view(segments[first].0, CS.fundamental_rep)?;
        lines.push(shadowing::trace(
            fundamental_rep(dimension.to_owned()),
            word,
        ));
    }
    Some(Atom::mul_many(others.iter().chain(&lines)))
}

/// A fundamental line segment from cof(d,i) to dind(cof(d,j)): an explicit
/// generator t(a, cof(d,i), dind(cof(d,j))) or a nonempty chain between those
/// slots, as its start, the base slot of its end and its word.
fn fundamental_line_segment(
    factor: AtomView<'_>,
) -> Option<(AtomView<'_>, AtomView<'_>, Vec<Atom>)> {
    let AtomView::Fun(function) = factor else {
        return None;
    };
    let chain = function.get_symbol() == T.chain;
    if !chain && (function.get_symbol() != CS.t || function.get_nargs() != 3) {
        return None;
    }
    let mut args = function.iter();
    let adjoint = if chain {
        None
    } else {
        let adjoint = args.next()?;
        representation_slot_view(adjoint, CS.adjoint_rep)?;
        Some(adjoint)
    };
    let (start, end) = (args.next()?, args.next()?);
    let (start_dimension, _) = representation_slot_view(start, CS.fundamental_rep)?;
    let AtomView::Fun(dual) = end else {
        return None;
    };
    if dual.get_symbol() != AIND_SYMBOLS.dind || dual.get_nargs() != 1 {
        return None;
    }
    let end = dual.iter().next()?;
    let (end_dimension, _) = representation_slot_view(end, CS.fundamental_rep)?;
    if start_dimension != end_dimension {
        return None;
    }
    let word = match adjoint {
        Some(adjoint) => vec![CS.chain_t(adjoint)],
        None => args.map(|factor| factor.to_owned()).collect(),
    };
    (!word.is_empty()).then_some((start, end, word))
}

/// Whether an indexed fundamental slot, `cof(d,i)` possibly dualized, occurs.
fn has_fundamental_slot(expr: AtomView<'_>) -> bool {
    let mut found = false;
    expr.visitor(&mut |node| {
        found |= matches!(node, AtomView::Fun(function)
            if function.get_symbol() == CS.fundamental_rep && function.get_nargs() == 2);
        !found
    });
    found
}

/// Generators or ports of a non-adjoint colour line, which the single-port
/// chain join and Fierz rules can reach.
fn atom_contains_fundamental_line(expr: AtomView<'_>) -> bool {
    let lines = [
        CS.t,
        LibraryRep::from(ColorFundamental {}).symbol(),
        LibraryRep::from(ColorSextet {}).symbol(),
    ];
    let mut found = false;
    expr.visitor(&mut |node| {
        found |= matches!(node, AtomView::Fun(function) if lines.contains(&function.get_symbol()));
        !found
    });
    found
}

fn color_node_search(expr: AtomView<'_>, skip_invariants: bool) -> bool {
    let mut slots = SlotMatcher::default();
    let mut selected = false;
    expr.visitor(&mut |node| {
        if selected || !matches!(slots.classify(node), SlotMatch::Other) {
            return false;
        }
        if let AtomView::Fun(function) = node {
            let symbol = function.get_symbol();
            if skip_invariants && [CS.gram, CS.cas, CS.idx].contains(&symbol) {
                return false;
            }
            selected = TensorCollectFilter::Reps([
                ColorAdjoint {}.into(),
                ColorFundamental {}.into(),
                ColorSextet {}.into(),
            ])
            .matches(node)
                || [CS.f, CS.d, CS.t, CS.gram, CS.cas, CS.idx].contains(&symbol);
        }
        !selected
    });
    selected
}

#[derive(Clone, Copy)]
struct CommonStructurePosition {
    left: usize,
    right: usize,
}

fn common_structure_positions<'a, A: AtomCore>(
    left: &'a [A; 3],
    right: &'a [A; 3],
) -> impl Iterator<Item = CommonStructurePosition> + 'a {
    let mut right_used = [false; 3];
    left.iter().enumerate().filter_map(move |(left, left_arg)| {
        let right = right.iter().enumerate().position(|(right, right_arg)| {
            !right_used[right] && right_arg.as_atom_view() == left_arg.as_atom_view()
        })?;
        right_used[right] = true;
        Some(CommonStructurePosition { left, right })
    })
}

type ColourPorts = std::collections::HashSet<Slot<LibraryRep, AbstractIndex>>;

/// Occurrences of each explicit slot of an expression, a dual folded onto its
/// slot. Summands are alternatives, so a sum counts its largest summand: a
/// slot counted more than twice is then over-used in some term of the
/// expanded expression. A power of a tensor contracts its copies with each
/// other, so it is a closed scope: its labels are `bound`, never counted
/// with the occurrences around it. A malformed slot makes the counts
/// `opaque`. Products hold few distinct slots, so the counts are short lists.
#[derive(Default)]
struct SlotCounts<'a> {
    counts: Vec<(SlotKey<'a>, usize, AtomView<'a>)>,
    bound: Vec<SlotKey<'a>>,
    opaque: bool,
}

/// A slot's representation head, dimension and index, without its variance.
type SlotKey<'a> = (Symbol, AtomView<'a>, AtomView<'a>);

impl<'a> SlotCounts<'a> {
    fn of(expression: AtomView<'a>) -> Self {
        Self::of_factors(std::iter::once(expression))
    }

    fn of_factors(factors: impl IntoIterator<Item = AtomView<'a>>) -> Self {
        let mut matcher = SlotMatcher::default();
        let mut counts = Self::default();
        for factor in factors {
            counts.count(factor, &mut matcher, None);
        }
        counts
    }

    /// The counts of the slots of `expression` that `other` counts too.
    fn shared_with(expression: AtomView<'a>, other: &Self) -> Self {
        let mut counts = Self::default();
        if !other.counts.is_empty() {
            counts.count(expression, &mut SlotMatcher::default(), Some(other));
        }
        counts
    }

    fn get(&self, key: &SlotKey<'a>) -> Option<usize> {
        self.counts
            .iter()
            .find(|(seen, ..)| seen == key)
            .map(|(_, count, _)| *count)
    }

    /// Whether the expression spells this slot, also inside a power.
    fn spells(&self, key: &SlotKey<'a>) -> bool {
        self.get(key).is_some() || self.bound.contains(key)
    }

    /// Combine `count` occurrences of a slot with its present count.
    fn merge(
        &mut self,
        key: SlotKey<'a>,
        count: usize,
        node: AtomView<'a>,
        combine: fn(usize, usize) -> usize,
    ) {
        match self.counts.iter_mut().find(|(seen, ..)| *seen == key) {
            Some((_, present, _)) => *present = combine(*present, count),
            None => self.counts.push((key, count, node)),
        }
    }

    /// Add the occurrences in `expression` to these counts, of the slots
    /// `only` counts if given.
    fn count(&mut self, expression: AtomView<'a>, matcher: &mut SlotMatcher, only: Option<&Self>) {
        match expression {
            AtomView::Add(sum) => {
                let mut largest = Self::default();
                for summand in sum.iter() {
                    let mut counts = Self::default();
                    counts.count(summand, matcher, only);
                    for (key, count, node) in counts.counts {
                        largest.merge(key, count, node, usize::max);
                    }
                    largest.bound.extend(counts.bound);
                    largest.opaque |= counts.opaque;
                }
                for (key, count, node) in largest.counts {
                    self.merge(key, count, node, |present, count| present + count);
                }
                self.bound.extend(largest.bound);
                self.opaque |= largest.opaque;
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    self.count(factor, matcher, only);
                }
            }
            AtomView::Pow(power) => {
                let mut scope = Self::default();
                scope.count(power.get_base_exp().0, matcher, only);
                self.bound.extend(
                    scope
                        .counts
                        .into_iter()
                        .map(|(key, ..)| key)
                        .chain(scope.bound),
                );
                self.opaque |= scope.opaque;
            }
            AtomView::Fun(function) => match matcher.classify(expression) {
                SlotMatch::Explicit(slot) => {
                    let key = (slot.representation().head(), slot.dimension(), slot.index());
                    if only.is_none_or(|only| only.spells(&key)) {
                        self.merge(key, 1, expression, |present, count| present + count);
                    }
                }
                // A compact representation such as a trace's names no slot.
                SlotMatch::Opaque => {
                    self.opaque |= matcher.compact_representation(expression).is_none();
                }
                SlotMatch::Other => {
                    for arg in function.iter() {
                        self.count(arg, matcher, only);
                    }
                }
            },
            _ => {}
        }
    }

    /// The slots occurring once: in a sum, the ports of every summand.
    fn ports(&self) -> impl Iterator<Item = &SlotKey<'a>> {
        self.counts
            .iter()
            .filter(|(_, count, _)| *count == 1)
            .map(|(key, ..)| key)
    }

    #[cfg(test)]
    fn over_used(&self) -> bool {
        self.counts.iter().any(|(_, count, _)| *count > 2)
    }

    /// Whether the two expressions spell a common slot, also inside powers.
    fn meets(&self, other: &Self) -> bool {
        self.counts.iter().any(|(key, ..)| other.spells(key))
            || self.bound.iter().any(|key| other.spells(key))
    }

    /// The slots occurring at least twice here, outside powers, that
    /// `outside` also spells.
    fn captured_by(&self, outside: &Self) -> Vec<AtomView<'a>> {
        self.counts
            .iter()
            .filter(|(key, count, _)| *count >= 2 && outside.spells(key))
            .map(|(.., node)| *node)
            .collect()
    }
}

/// The port slots of a colour sum, with their duals, when it has any. A power
/// of such a sum contracts its copies with each other.
fn colour_power_ports(base: AtomView<'_>) -> Option<ColourPorts> {
    if !atom_contains_color_port_or_work(base) {
        return None;
    }
    let ports = crate::tensor::inference::InterfaceInference::replacement_interface(base)
        .ok()?
        .slots()
        .ok()?;
    (!ports.is_empty()).then(|| {
        ports
            .into_iter()
            .flat_map(|slot| [slot.dual(), slot])
            .collect()
    })
}

fn positive_integer(expr: AtomView) -> Option<i64> {
    let AtomView::Num(number) = expr else {
        return None;
    };
    let CoefficientView::Natural(value, 1, 0, 1) = number.get_coeff_view() else {
        return None;
    };

    (value > 0).then_some(value)
}

#[cfg(test)]
mod reconstruction_tests {
    use super::*;
    use std::sync::{Arc, Mutex};
    use symbolica::{domains::float::Float, parse, parse_lit, symbol};

    #[test]
    fn legacy_symmetric_invariant_reader_requires_a_uniform_adjoint_index_space() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        for representation in [
            fundamental_rep(Atom::num(3)),
            ColorSextet {}.to_symbolic([Atom::num(6)]),
        ] {
            for degree in [4, 6] {
                for uniform in [false, true] {
                    let slots = (0..degree)
                        .map(|index| {
                            ColorAdjoint {}.to_symbolic([
                                Atom::num(if uniform || index != 0 { 8 } else { 3 }),
                                Atom::num(99851 + index),
                            ])
                        })
                        .collect::<Vec<_>>();
                    // This is the legacy internal d_R representation. Public
                    // admission uses a symmetric trace; its representation
                    // argument must not be invented as an additional open port.
                    let invariant = CS.symmetric_d(&representation, slots);
                    assert_eq!(
                        SymmetricInvariantView::parse(invariant.as_view()).is_some(),
                        uniform
                    );
                    let result = simplifier.simplify_power(invariant.pow(2).as_view());
                    let gram = CS.gram(Atom::num(degree), &representation, &representation);
                    if uniform {
                        assert_eq!(result, Some(gram));
                    } else {
                        assert_eq!(result, None);
                    }
                }
            }
        }
    }

    /// An f on two legs of a line at distance d shortens the line by one in
    /// each of its 1 + (d - 1) terms; the closest pair is taken first.
    #[test]
    fn structure_constants_on_distant_legs_shorten_the_line() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let slots = adjoint_slots("distant_legs", ["a0", "a1", "a2", "a3", "a4", "a5", "x"]);
        let (incoming, outgoing) = fundamental_endpoints(61);
        let longest = |result: &Atom| {
            let mut lengths = Vec::new();
            result.as_view().visitor(&mut |node| {
                if let Some((_, factors)) = trace_parts(node) {
                    lengths.push(factors.len());
                } else if let Some((_, _, factors)) = chain_parts(node) {
                    lengths.push(factors.len());
                }
                true
            });
            lengths
        };
        for (first, second, distance) in [(0, 2, 2), (0, 3, 3), (3, 0, 3), (4, 1, 3)] {
            let f = color_f!(&slots[first], &slots[second], &slots[6]);
            for adjoint in [false, true] {
                let rep = if adjoint {
                    adjoint_rep(Atom::num(8))
                } else {
                    fundamental_rep(Atom::num(3))
                };
                let word = slots[..6]
                    .iter()
                    .map(|slot| ColorAlgebraSimplifier::trace_generator_factor(adjoint, slot))
                    .collect::<Vec<_>>();
                let source = trace!(&rep; &word) * &f;
                SymbolicTensor::infer(source.clone()).unwrap();
                let result = simplifier
                    .simplify_trace_structure_product(&ProductView::parse(source.as_view()))
                    .unwrap();
                let lengths = longest(&result);
                assert_eq!(lengths.len(), distance, "{result}");
                assert!(lengths.iter().all(|&length| length == 5), "{result}");
            }
            if first < second {
                let source =
                    chain!(&incoming, &outgoing; slots[..6].iter().map(|slot| color_t!(slot))) * &f;
                SymbolicTensor::infer(source.clone()).unwrap();
                let result = simplifier
                    .simplify_chain_structure_product(&ProductView::parse(source.as_view()))
                    .unwrap();
                let lengths = longest(&result);
                assert_eq!(lengths.len(), distance, "{result}");
                assert!(lengths.iter().all(|&length| length == 5), "{result}");
            }
        }
    }

    #[test]
    fn structure_word_connections_require_unique_ports_in_the_whole_line() {
        crate::test_support::test_initialize();
        let slots = (1..7)
            .map(|i| ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(i)]))
            .collect::<Vec<_>>();
        let incoming = ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(11)]);
        let outgoing =
            spenso::dind!(ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(12)]));
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        for repeated in [false, true] {
            let word = [0, 1, if repeated { 0 } else { 2 }, 3].map(|i| color_t!(&slots[i]));
            let line = trace!(fundamental_rep(Atom::num(3)); &word);
            SymbolicTensor::infer(line.clone()).unwrap();
            let input = &line * color_f!(&slots[0], &slots[1], &slots[4]);
            assert_eq!(
                simplifier
                    .simplify_trace_structure_product(&ProductView::parse(input.as_view()))
                    .is_some(),
                !repeated,
            );
            let line = chain!(&incoming, &outgoing; &word);
            SymbolicTensor::infer(line.clone()).unwrap();
            let input = line * color_f!(&slots[0], &slots[1], &slots[4]);
            assert_eq!(
                simplifier
                    .simplify_chain_structure_product(&ProductView::parse(input.as_view()))
                    .is_some(),
                !repeated,
            );
        }
    }

    #[test]
    fn cross_line_fierz_never_uses_a_port_contracted_inside_either_line() {
        crate::test_support::test_initialize();
        let slots = (21..28)
            .map(|i| ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(i)]))
            .collect::<Vec<_>>();
        let incoming = ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(31)]);
        let outgoing =
            spenso::dind!(ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(32)]));
        for repeated in [false, true] {
            let factors = [0, 1, if repeated { 0 } else { 2 }, 3].map(|i| color_t!(&slots[i]));
            let other =
                trace!(fundamental_rep(Atom::num(3)); [color_t!(&slots[0]), color_t!(&slots[4])]);
            for line in [
                trace!(fundamental_rep(Atom::num(3)); &factors),
                chain!(&incoming, &outgoing; &factors),
            ] {
                // The raw product deliberately repeats an internal spelling;
                // each admitted line still has its own established boundary.
                SymbolicTensor::infer(line.clone()).unwrap();
                let input = line * &other;
                assert_eq!(
                    ColorAlgebraSimplifier::simplify_cross_chain_fierz_product(
                        &ProductView::parse(input.as_view()),
                        true,
                    )
                    .is_some(),
                    !repeated,
                );
            }
        }
    }

    fn adjoint_slots<const N: usize>(scope: &str, labels: [&str; N]) -> [Atom; N] {
        labels.map(|label| {
            ColorAdjoint {}.to_symbolic([
                Atom::num(8),
                Atom::var(symbol!(&format!("{scope}::{label}"))),
            ])
        })
    }

    fn fundamental_endpoints(first: i64) -> (Atom, Atom) {
        (
            ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(first)]),
            spenso::dind!(ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(first + 1)])),
        )
    }

    #[test]
    fn structure_loops_preserve_orientation_and_validate_all_ports() {
        crate::test_support::test_initialize();
        let [a, b, x, y, z] = adjoint_slots("borrowed_structure_loop", ["a", "b", "x", "y", "z"]);
        let permutations = [
            ([0, 1, 2], 1),
            ([1, 2, 0], 1),
            ([2, 0, 1], 1),
            ([1, 0, 2], -1),
            ([0, 2, 1], -1),
            ([2, 1, 0], -1),
        ];
        let casimir = adjoint_casimir_for_dimension(Atom::num(8));
        for (open, result) in [(&y, color_metric(x.clone(), y.clone())), (&x, Atom::num(8))] {
            for (left_order, left_sign) in permutations {
                for (right_order, right_sign) in permutations {
                    let left = StructureView {
                        args: left_order.map(|i| [&a, &b, &x][i].as_view()),
                    };
                    let right = StructureView {
                        args: right_order.map(|i| [&a, &b, open][i].as_view()),
                    };
                    assert_eq!(
                        left.contract_loop(&right),
                        Some(Atom::num(left_sign * right_sign) * &casimir * &result)
                    );
                }
            }
        }
        let left = StructureView {
            args: [&a, &b, &x].map(Atom::as_view),
        };
        let single = StructureView {
            args: [&a, &y, &z].map(Atom::as_view),
        };
        assert!(left.contract_loop(&single).is_none());
        let wrong_dimension = ColorAdjoint {}.to_symbolic([Atom::num(15), Atom::num(701)]);
        let mixed = StructureView {
            args: [&a, &b, &wrong_dimension].map(Atom::as_view),
        };
        assert!(left.contract_loop(&mixed).is_none());
    }

    #[test]
    fn emitted_adjoint_metrics_are_eliminated_against_their_one_partner() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let casimir = adjoint_casimir_for_dimension(Atom::num(8));
        let (incoming, outgoing) = fundamental_endpoints(41);
        // Wildcard-like labels stay literal: renaming compares atoms and
        // builds no pattern from them.
        for labels in [["a", "b", "x", "y"], ["a_", "b_", "x_", "y_"]] {
            let [a, b, x, y] = adjoint_slots("emitted_metric", labels);
            let line = |slot: &Atom| chain!(&incoming, &outgoing; [color_t!(slot)]);
            let source = color_f!(&a, &b, &x) * color_f!(&a, &b, &y) * line(&y);
            SymbolicTensor::infer(source.clone()).unwrap();
            let result = simplifier.step(source.as_view(), true);
            assert_eq!(result, &casimir * line(&x), "{labels:?}");
            // With both ends open the metric is the result.
            let open = simplifier.step(
                (color_f!(&a, &b, &x) * color_f!(&a, &b, &y)).as_view(),
                true,
            );
            assert_eq!(
                open,
                &casimir * color_metric(x.clone(), y.clone()),
                "{labels:?}"
            );
            // A partner inside a power stands for two copies, and one inside
            // a foreign tensor is not the kernel's; the contractor keeps both.
            let [c, z] = adjoint_slots("emitted_metric_partner", ["c", "z"]);
            let foreign = function!(spenso::tensor_symbol!("emitted_metric_partner::V"), &y);
            for partner in [color_f!(&y, &c, &z).pow(Atom::num(2)), foreign] {
                let replacement = &casimir * color_metric(x.clone(), y.clone());
                let kept = ProductView::eliminating_emitted_metrics(
                    &replacement,
                    std::iter::once(partner.as_view()),
                );
                assert_eq!(kept, Some(replacement * &partner), "{labels:?}");
            }
        }
    }

    #[test]
    fn symmetric_invariance_annihilates_a_bridge_between_traces() {
        use crate::tensor::{AlgebraContraction, AlgebraSettings};
        crate::test_support::test_initialize();
        let rep = fundamental_rep(Atom::num(3));
        let adjoint = adjoint_rep(Atom::num(8));
        let [a, b, c, d, e, x, y] =
            adjoint_slots("invariant_bridge", ["a", "b", "c", "d", "e", "x", "y"]);
        let left = color_symmetric_trace(&rep, [a.clone(), b.clone(), c.clone()]);
        let right = shadowing::trace_sym(
            &adjoint,
            [d.clone(), b.clone(), c.clone(), e.clone()]
                .map(|slot| ColorAlgebraSimplifier::trace_generator_factor(true, &slot)),
        );
        let settings = AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        let source = color_f!(&x, &a, &d) * &left * &right;
        assert!(!source.is_zero());
        let reduced = SymbolicTensor::infer(source)
            .unwrap()
            .simplify_algebra(&settings)
            .unwrap();
        assert!(reduced.expression().is_zero());

        // A missing connection does not sum the action on every leg of D_R.
        let open = shadowing::trace_sym(
            &adjoint,
            [d.clone(), b.clone(), y, e]
                .map(|slot| ColorAlgebraSimplifier::trace_generator_factor(true, &slot)),
        );
        let source = color_f!(x, a, d) * left * open;
        assert!(
            ColorAlgebraSimplifier::simplify_symmetric_structure_product(&ProductView::parse(
                source.as_view()
            ))
            .is_none()
        );
        assert!(
            !SymbolicTensor::infer(source)
                .unwrap()
                .simplify_algebra(&settings)
                .unwrap()
                .expression()
                .is_zero()
        );
    }

    #[test]
    fn symmetric_invariants_cannot_leave_one_adjoint_port() {
        use crate::tensor::{AlgebraContraction, AlgebraSettings};
        crate::test_support::test_initialize();
        let [a, b, c, d, e] = adjoint_slots("invariant_vector", ["a", "b", "c", "d", "e"]);
        let rep = fundamental_rep(Atom::num(3));
        let adjoint = adjoint_rep(Atom::num(8));
        let left = CS.symmetric_d(&rep, vec![a.clone(), b.clone(), c.clone()]);
        let right = CS.symmetric_d(&adjoint, vec![a.clone(), b.clone(), c.clone(), d.clone()]);
        let settings = AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        let source = &left * right;
        assert!(!source.is_zero());
        assert!(
            SymbolicTensor::infer(source)
                .unwrap()
                .simplify_algebra(&settings)
                .unwrap()
                .expression()
                .is_zero()
        );
        let open = CS.symmetric_d(&adjoint, vec![a, b, d, e]);
        let source = left * open;
        assert!(
            ColorAlgebraSimplifier::simplify_symmetric_invariant_product(&ProductView::parse(
                source.as_view()
            ))
            .is_none()
        );
    }

    #[test]
    fn quark_box_zero_is_independent_of_index_order() {
        use crate::tensor::{AlgebraContraction, AlgebraSettings};
        crate::test_support::test_initialize();
        let rep = fundamental_rep(Atom::num(3));
        let ports = adjoint_slots(
            "quark_box_zero",
            ["a", "b", "c", "d", "e", "h", "k", "x", "y"],
        );
        let settings = AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        // FK3426/FK3430 after propagator metrics: a four-generator quark
        // loop joined to four three-gluon vertices, with two open color ports.
        for reverse in [false, true] {
            for shift in 0..ports.len() {
                let order = |i: usize| {
                    if reverse {
                        (ports.len() - 1 - i + shift) % ports.len()
                    } else {
                        (i + shift) % ports.len()
                    }
                };
                let [a, b, c, d, e, h, k, x, y] =
                    std::array::from_fn::<_, 9, _>(|i| &ports[order(i)]);
                let source = trace!(&rep, color_t!(a), color_t!(c), color_t!(b), color_t!(d))
                    * color_f!(x, d, e)
                    * color_f!(y, e, h)
                    * color_f!(a, b, k)
                    * color_f!(c, h, k);
                let reduced = SymbolicTensor::infer(source)
                    .unwrap()
                    .simplify_algebra(&settings)
                    .unwrap();
                assert!(
                    reduced.expression().is_zero(),
                    "reverse={reverse}, shift={shift}: {}",
                    reduced.expression()
                );
            }
        }
    }

    /// A cut of a cycle of length five waits while a generator trace of the
    /// same term is still ordered; a terminal symmetric invariant does not
    /// hold it back.
    #[test]
    fn long_cycle_cuts_wait_for_a_pending_line() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let internal = adjoint_slots("pending_line", ["x0", "x1", "x2", "x3", "x4"]);
        let external = adjoint_slots("pending_line", ["e0", "e1", "e2", "e3", "e4"]);
        let [p, q] = adjoint_slots("pending_line", ["p", "q"]);
        let cycle = Atom::mul_many(
            (0..5).map(|i| color_f!(&internal[i], &internal[(i + 1) % 5], &external[i])),
        );
        let rep = fundamental_rep(Atom::num(3));
        let ordered = trace!(&rep, color_t!(&external[0]), color_t!(&p), color_t!(&q));
        let symmetric = color_symmetric_trace(&rep, [external[0].clone(), p, q]);
        for (line, waits) in [(ordered, true), (symmetric, false)] {
            let product = &cycle * line;
            SymbolicTensor::infer(product.clone()).unwrap();
            let cut = simplifier
                .simplify_adjoint_loop_product(&ProductView::parse(product.as_view()), true);
            assert_eq!(cut.is_none(), waits, "{product}");
        }
    }

    /// X^2 = sum_i x_i X renames the internal dummies of each x_i. A dummy of
    /// another representation spelled like a colour port is internal too.
    #[test]
    fn colour_power_renames_dummies_that_share_a_port_label() {
        use crate::IndexTooling;
        use crate::tensor::{AlgebraContraction, AlgebraSettings};
        crate::test_support::test_initialize();
        let adjoint = ColorAdjoint {}.new_rep(symbol!("power_ports::Na"));
        let label = |name: &str| Atom::var(symbol!(&format!("power_ports::{name}")));
        let [p, q] = ["p", "q"].map(|name| adjoint.to_symbolic([label(name)]));
        let rep = ColorFundamental {}
            .new_rep(symbol!("power_ports::Nc"))
            .to_symbolic([]);
        let [v, w] = ["V", "W"].map(|name| spenso::tensor_symbol!(&format!("power_ports::{name}")));
        let euclidean = |name: &str| {
            spenso::structure::representation::Euclidean {}
                .new_rep(4)
                .to_symbolic([label(name)])
        };
        let dot = |name: &str| function!(v, euclidean(name)) * function!(w, euclidean(name));
        let sum = trace!(&rep, color_t!(&p), color_t!(&q)) * dot("p")
            + color_metric(p.clone(), q.clone()) * dot("r");
        let result = SymbolicTensor::infer(sum.pow(Atom::num(2)))
            .unwrap()
            .simplify_algebra(&AlgebraSettings {
                color: Some(ColorSimplifySettings::default()),
                contract: AlgebraContraction::None,
                ..Default::default()
            })
            .unwrap();
        // Capturing euc(4,p) would leave it on three ports.
        let readmitted = SymbolicTensor::infer(result.expression().clone()).unwrap();
        let fresh = ParseState::<AbstractIndex>::default();
        let [d0, d1] = [fresh.fresh_index(), fresh.fresh_index()].map(Atom::from);
        let paired = |index: &Atom| {
            let slot = spenso::structure::representation::Euclidean {}
                .new_rep(4)
                .to_symbolic([index.clone()]);
            function!(v, &slot) * function!(w, slot)
        };
        let expected = (quadratic_index(rep.clone()) + Atom::one()).pow(Atom::num(2))
            * Atom::var(symbol!("power_ports::Na"))
            * paired(&d0)
            * paired(&d1);
        let difference = (readmitted.expression().clone() - expected)
            .canonize(AbstractIndex::Dummy)
            .unwrap()
            .expand();
        assert!(difference.is_zero(), "{}", result.expression());
    }

    /// Distributing a sum into a product renames the summand's own dummies
    /// that another factor spells too, as in two copies of one sum; the
    /// canonical relabelling declines such a product instead.
    #[test]
    fn distributed_summands_never_capture_outer_labels() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let [a, b, c] = adjoint_slots("capture", ["a", "b", "c"]);
        let rep = fundamental_rep(Atom::num(3));
        let copy = trace!(&rep, color_t!(&a), color_t!(&c))
            * trace!(&rep, color_t!(&b), color_t!(&c))
            + trace!(&rep, color_t!(&a), color_t!(&c)) * color_metric(b.clone(), c.clone());
        let other = trace!(&rep, color_t!(&a), color_t!(&c))
            * trace!(&rep, color_t!(&b), color_t!(&c))
            + Atom::num(2)
                * trace!(&rep, color_t!(&a), color_t!(&c))
                * color_metric(b.clone(), c.clone());
        let outside = SlotCounts::of(other.as_view());
        let AtomView::Add(sum) = copy.as_view() else {
            panic!("{copy}")
        };
        for summand in sum.iter() {
            assert!(SlotCounts::of((summand.to_owned() * &other).as_view()).over_used());
            let renamed = simplifier.without_captured_dummies(summand, &outside);
            assert!(!renamed.contains(c.as_view()), "{renamed}");
            let term = &renamed * &other;
            assert!(!SlotCounts::of(term.as_view()).over_used(), "{term}");
        }
        let shadowed = &copy * trace!(&rep, color_t!(&c), color_t!(&b));
        assert!(SlotCounts::of(shadowed.as_view()).over_used());
        assert!(states::colour_states_for_test(shadowed.as_view()).is_none());
    }

    /// Foreign colour factors are the tensors with colour slots and no
    /// colour algebra that no colour algebra reaches through shared slots:
    /// in the product, in each term of a sum and in each term of a colour sum
    /// factor beside its outer factors. A metric is colour algebra.
    #[test]
    fn foreign_colour_factors_are_unreached_by_colour_algebra() {
        crate::test_support::test_initialize();
        let [a, b, c, d, y, z] = adjoint_slots("foreign", ["a", "b", "c", "d", "y", "z"]);
        let [v, u, m] = [
            spenso::tensor_symbol!("idenso::foreign_test::V"),
            spenso::tensor_symbol!("idenso::foreign_test::U"),
            spenso::tensor_symbol!("idenso::foreign_test::M"),
        ];
        let x = Atom::var(symbol!("idenso::foreign_test::x"));
        let loop_ = color_f!(&a, &b, &c) * color_f!(&a, &b, &c);
        let colour_sum = &x * color_f!(&a, &b, &c) + color_f!(&a, &c, &b);
        let cases = [
            (function!(v, &z) * &loop_, vec![function!(v, &z)]),
            (
                function!(v, &c) * color_f!(&a, &b, &c) * color_f!(&a, &b, &d) * function!(u, &d),
                vec![],
            ),
            (
                function!(v, &z) * function!(u, &z) * &loop_,
                vec![function!(v, &z), function!(u, &z)],
            ),
            (function!(m, &y, &y) * &loop_, vec![function!(m, &y, &y)]),
            (
                function!(v, &z) * color_f!(&a, &b, &c) * &colour_sum + function!(u, &z) * &loop_,
                vec![function!(v, &z), function!(u, &z)],
            ),
            (
                color_f!(&a, &b, &c)
                    * (function!(v, &z) * color_f!(&a, &b, &c)
                        + function!(u, &z) * color_f!(&a, &c, &b)),
                vec![function!(v, &z), function!(u, &z)],
            ),
            (
                color_f!(&a, &b, &c)
                    * (function!(v, &c) * color_f!(&a, &b, &y)
                        + function!(u, &y) * color_f!(&a, &b, &c)),
                vec![function!(u, &y)],
            ),
            (
                function!(v, &z)
                    * ETS.metric_literal(&z, &y)
                    * color_f!(&a, &b, &y)
                    * color_f!(&a, &b, &d),
                vec![],
            ),
            (&x * &loop_, vec![]),
        ];
        for (expression, expected) in cases {
            let mut found = foreign_colour_factors(expression.as_view())
                .into_iter()
                .map(|factor| factor.to_owned())
                .collect::<Vec<_>>();
            found.sort();
            let mut expected = expected;
            expected.sort();
            assert_eq!(found, expected, "{expression}");
        }
    }

    /// A power of a tensor is a closed scope. A summand holding one beside a
    /// port of the same spelling keeps the port; the power's labels are
    /// renamed when another factor or the rest of the summand spells them,
    /// and the power names no port of the summand.
    #[test]
    fn distributed_powers_are_closed_scopes() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let [c, e, g, x] = adjoint_slots("power_scope", ["c", "e", "g", "x"]);
        let rep = fundamental_rep(Atom::num(3));
        let power = trace!(&rep, color_t!(&c), color_t!(&x)).pow(Atom::num(2));
        let port = trace!(&rep, color_t!(&c), color_t!(&e), color_t!(&g));
        let summand = &power * &port;
        let other = color_f!(&c, &g, &e);
        let outside = SlotCounts::of(other.as_view());
        let counts = SlotCounts::of(summand.as_view());
        assert_eq!(counts.ports().count(), 3, "{summand}");
        assert!(counts.bound.len() == 2 && outside.spells(&counts.bound[0]));
        let renamed = simplifier.without_captured_dummies(summand.as_view(), &outside);
        let AtomView::Mul(product) = renamed.as_view() else {
            panic!("{renamed}")
        };
        let fresh_power = product
            .iter()
            .find(|factor| matches!(factor, AtomView::Pow(_)))
            .unwrap();
        assert!(!fresh_power.contains(c.as_view()) && !fresh_power.contains(x.as_view()));
        assert!(
            product.iter().any(|factor| factor == port.as_view()),
            "{renamed}"
        );
        // A power that shares nothing keeps its labels.
        let lone = &power * trace!(&rep, color_t!(&e), color_t!(&g));
        let unrelated = color_f!(&e, &g, &adjoint_slots("power_scope", ["z"])[0]);
        assert_eq!(
            simplifier
                .without_captured_dummies(lone.as_view(), &SlotCounts::of(unrelated.as_view())),
            lone
        );
    }

    /// Following the indices of explicit generators and fundamental chains
    /// collects the same lines as the pattern-based `join_chains`: closed
    /// loops, open lines, metrics left alone, lines joined within a sum's
    /// terms only, and terms of a sum.
    #[test]
    fn fundamental_line_join_matches_the_shared_collector() {
        crate::test_support::test_initialize();
        let [a1, a2, a3, a4, x] = adjoint_slots("line_join", ["a1", "a2", "a3", "a4", "x"]);
        let slot = |index: i64| ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(index)]);
        let t = |adjoint: &Atom, start: i64, end: i64| {
            CS.explicit_t(adjoint.clone(), slot(start), spenso::dind!(slot(end)))
        };
        let chain = chain!(slot(1), spenso::dind!(slot(2)); [CS.chain_t(a1.clone())]);
        let metric = ETS.metric_literal(slot(2), spenso::dind!(slot(3)));
        let loop_term = t(&a1, 1, 2) * t(&a2, 2, 1);
        let cases = [
            loop_term.clone(),
            t(&a1, 1, 2) * t(&a2, 2, 3) * t(&a3, 3, 1) * color_f!(&a1, &a2, &x),
            &chain * t(&a2, 2, 3),
            &chain * t(&a2, 2, 1),
            t(&a1, 1, 2) * &metric,
            (t(&a1, 1, 2) * color_f!(&a1, &a3, &x) + t(&a2, 1, 2) * color_f!(&a2, &a3, &x))
                * t(&a4, 2, 1),
            loop_term + t(&a3, 1, 2) * t(&a4, 2, 1),
            t(&a1, 4, 4),
        ];
        for expression in cases {
            assert_eq!(
                joined_fundamental_lines(expression.as_view()),
                Some(expression.join_chains(ColorFundamental {}.into())),
                "{expression}"
            );
        }
        // A foreign matrix with fundamental slots is left to the shared collector.
        let foreign = function!(
            spenso::tensor_symbol!("line_join::M"),
            slot(2),
            spenso::dind!(slot(1))
        );
        assert!(joined_fundamental_lines((t(&a1, 1, 2) * foreign).as_view()).is_none());
    }

    /// Of two traces sharing generators, the one with fewer distinct external
    /// generators is rewritten first (color.h's cOlTT), whatever the factor
    /// order of the product.
    #[test]
    fn the_line_with_fewer_external_generators_goes_first() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let slots = adjoint_slots(
            "line_order",
            ["a0", "a1", "a2", "a3", "a4", "a5", "z1", "z2"],
        );
        let rep = fundamental_rep(Atom::num(3));
        let word = |indices: [usize; 6]| trace!(&rep; indices.map(|index| color_t!(&slots[index])));
        // Six external generators against two: a0, a3 and two contracted pairs.
        let open = word([0, 1, 2, 3, 4, 5]);
        let closed = word([0, 3, 6, 6, 7, 7]);
        let product = &open * &closed;
        let view = ProductView::parse(product.as_view());
        // The product order alone would rewrite the open trace.
        assert_eq!(view.factors[0].atom, open.as_view());
        let rewritten = simplifier
            .simplify_embedded_color_node(&view, true, true)
            .unwrap();
        assert!(rewritten.contains(open.as_view()), "{rewritten}");
        assert!(!rewritten.contains(closed.as_view()), "{rewritten}");
    }

    /// color.h's two-block identity on a trace whose pairs are two apart.
    #[test]
    fn repeated_blocks_reverse_into_two_terms() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let [a, b, c, d] = adjoint_slots("block", ["a", "b", "c", "d"]);
        let rep = fundamental_rep(Atom::num(3));
        let generators = [&a, &b, &c, &a, &b, &d].map(|slot| slot.clone());
        let result = simplifier
            .simplify_repeated_generator_trace(&rep, &generators)
            .unwrap();
        let reversed = trace!(&rep; [&a, &b, &c, &b, &a, &d].map(|slot| color_t!(slot)));
        let bridged = trace!(&rep; [&a, &c, &a, &d].map(|slot| color_t!(slot)));
        assert_eq!(
            result,
            reversed - adjoint_casimir_for_dimension(Atom::num(8)) / Atom::num(2) * bridged
        );
    }

    #[test]
    fn term_repeat_finishes_local_rules_and_its_bound_leaves_the_row_unfixed() {
        crate::test_support::test_initialize();
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let [a, b, c, d, x, y] = adjoint_slots("term_repeat", ["a", "b", "c", "d", "x", "y"]);
        // Two bubbles in series: each two-f loop leaves a metric that the
        // next identity needs eliminated.
        let source = color_f!(&a, &b, &x)
            * color_f!(&a, &b, &y)
            * color_f!(&x, &c, &d)
            * color_f!(&y, &c, &d);
        SymbolicTensor::infer(source.clone()).unwrap();
        let casimir = adjoint_casimir_for_dimension(Atom::num(8));
        assert_eq!(
            simplifier.rewrite_terms(source.as_view(), true),
            casimir.clone().pow(Atom::num(2)) * Atom::num(8)
        );
        simplifier.repeat_depth.set(TERM_REPEAT_DEPTH);
        simplifier.certificate.fixed.set(true);
        let bounded = simplifier.rewrite_terms(source.as_view(), true);
        assert!(bounded.contains_symbol(CS.f));
        assert!(!simplifier.certificate.fixed.get());
    }

    #[test]
    fn cross_line_fierz_waits_for_maximal_chains() {
        crate::test_support::test_initialize();
        let [a, b, c, e] = adjoint_slots("maximal_fierz", ["a", "b", "c", "e"]);
        let rep = fundamental_rep(Atom::num(3));
        let other = trace!(rep; [color_t!(&a), color_t!(&e)]);
        let (start, end) = fundamental_endpoints(51);
        let line = chain!(&start, &end; [color_t!(&a), color_t!(&b)]);
        let continued = |slot: AtomView<'_>| {
            let (_, index, _) = color_fundamental_slot(slot).unwrap();
            ColorFundamental {}.to_symbolic([Atom::num(3), index])
        };
        let next_end =
            spenso::dind!(ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(59)]));
        let next = chain!(continued(end.as_view()), &next_end; [color_t!(&c)]);
        let metric = color_metric(continued(end.as_view()), next_end.clone());
        let fierz = |source: Atom| {
            SymbolicTensor::infer(source.clone()).unwrap();
            ColorAlgebraSimplifier::simplify_cross_chain_fierz_product(
                &ProductView::parse(source.as_view()),
                true,
            )
        };
        assert!(fierz(&line * &other).is_some());
        // A chain continued by another chain or a metric is joined first.
        assert!(fierz(&line * &next * &other).is_none());
        assert!(fierz(&line * &metric * &other).is_none());
        // Before the frontier is complete a continuation may still arrive.
        assert!(
            ColorAlgebraSimplifier::simplify_cross_chain_fierz_product(
                &ProductView::parse((&line * &other).as_view()),
                false,
            )
            .is_none()
        );
    }

    // Preserve the prior reconstruction schedule as an independent oracle for
    // rounding and user normalizers, which exact polynomial equality cannot test.
    fn prefix_rewrite(simplifier: &ColorAlgebraSimplifier, expression: AtomView<'_>) -> Atom {
        if let AtomView::Add(sum) = expression {
            return sum.iter().fold(Atom::Zero, |result, term| {
                result + prefix_rewrite(simplifier, term)
            });
        }
        if let Some(result) = simplifier.rewrite_node(expression, true, true) {
            return result;
        }
        expression.to_owned().replace_map(|node, _context, out| {
            if let Some(result) = simplifier.rewrite_node(node, true, false) {
                **out = result;
            }
        })
    }

    #[test]
    fn trace_casimir_pairs_are_cyclic_and_require_generator_middle() {
        let reps = crate::test_support::test_initialize();
        let [a, b, c, d] = [1, 2, 3, 4].map(|i| color_t!(reps.coad_da.to_symbolic([Atom::num(i)])));
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let mut adjacent = vec![a.clone(), b.clone(), c.clone(), a.clone()];
        let expected = quadratic_casimir(rep.clone())
            * trace_with_factors(rep.clone(), vec![b.clone(), c.clone()]);
        for _ in 0..adjacent.len() {
            assert_eq!(
                ColorAlgebraSimplifier::simplify_adjacent_trace_casimir(&rep, &adjacent),
                Some(expected.clone()),
            );
            adjacent.rotate_left(1);
        }
        let mut separated = vec![b.clone(), a.clone(), c.clone(), d.clone(), a.clone()];
        let expected = (quadratic_casimir(rep.clone())
            - adjoint_casimir_for_dimension(reps.coad_da.dim.to_symbolic()) / Atom::num(2))
            * trace_with_factors(rep.clone(), vec![b, c.clone(), d]);
        for _ in 0..separated.len() {
            assert_eq!(
                ColorAlgebraSimplifier::simplify_separated_trace_casimir(&rep, &separated),
                Some(expected.clone()),
            );
            separated.rotate_left(1);
        }
        let foreign = parse!("M(in,out)", default_namespace = "spenso");
        let factors = vec![a.clone(), foreign, a, c];
        // Keep the chain test linear and the trace long enough that the equal
        // endpoints have no alternative one-generator path around the cycle.
        let start = reps.cof_nc.to_symbolic([Atom::num(8)]);
        let end = reps.cof_nc.dual().to_symbolic([Atom::num(9)]);
        assert!(
            ColorAlgebraSimplifier::simplify_separated_chain_casimir(&start, &end, &factors)
                .is_none()
        );
        let mut factors = factors;
        factors.push(color_t!(reps.coad_da.to_symbolic([Atom::num(5)])));
        assert!(ColorAlgebraSimplifier::simplify_separated_trace_casimir(&rep, &factors).is_none());
    }

    #[test]
    fn incomplete_frontiers_keep_traces_for_product_rules() {
        crate::test_support::test_initialize();
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        let left = parse!(
            "trace(cof(Nc), cyclic(t(coad(A,b),in,out), t(coad(A,a),in,out), t(coad(A,c),in,out)))",
            default_namespace = "spenso"
        );
        let right = parse!(
            "trace(cof(Nc), cyclic(t(coad(A,d),in,out), t(coad(A,a),in,out), t(coad(A,e),in,out)))",
            default_namespace = "spenso"
        );
        assert_eq!(simplifier.step(left.as_view(), false), left);
        let terminal = simplifier.step(left.as_view(), true);
        assert_ne!(terminal, left);
        assert_eq!(simplifier.step(terminal.as_view(), true), terminal);

        let product = &left * &right;
        let prefix = simplifier.step(product.as_view(), false);
        assert_ne!(prefix, product);
        // Fierz is a product identity and must not be disabled with the terminal.
        assert_eq!(prefix, simplifier.step(product.as_view(), true));
    }

    #[test]
    fn rewritten_color_sum_preserves_callbacks_and_factored_spectators() {
        crate::test_support::test_initialize();
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        let calls = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&calls);
        let hook = symbol!(
            "color_reconstruction_hook",
            norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
        );
        let terminal = parse_lit!(trace(cof(Nc)), default_namespace = "spenso");
        let spectator = parse_lit!((x + y) ^ 7);
        let input = Atom::add_many((0..12).map(|i| function!(hook, i, &terminal) * &spectator));
        calls.lock().unwrap().clear();
        let expected = prefix_rewrite(&simplifier, input.as_view());
        let transcript = calls.lock().unwrap().clone();
        assert!(!transcript.is_empty());
        calls.lock().unwrap().clear();
        let actual = simplifier.rewrite_terms(input.as_view(), true);
        assert_eq!(actual, expected);
        assert_eq!(*calls.lock().unwrap(), transcript);
        calls.lock().unwrap().clear();
        assert_eq!(simplifier.rewrite_terms(actual.as_view(), true), actual);
        assert!(calls.lock().unwrap().is_empty());
        assert_eq!(
            actual,
            Atom::add_many((0..12).map(|i| function!(hook, i, Atom::var(CS.nc)) * &spectator))
        );
    }

    #[test]
    fn unchanged_color_sum_still_normalizes_unchecked_input() {
        crate::test_support::test_initialize();
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        let mut input = Atom::new();
        input.to_add();
        assert!(input.as_view().needs_normalization());
        assert_eq!(simplifier.rewrite_terms(input.as_view(), true), Atom::Zero);
    }

    #[test]
    fn rewritten_color_sums_keep_rounded_prefix_order() {
        crate::test_support::test_initialize();
        let x = parse_lit!(x);
        let rounded = |value| Atom::num(Float::parse(value, Some(11)).unwrap());
        for terms in [
            vec![rounded("1e10"), rounded("1"), rounded("-1e10")],
            vec![
                rounded("1e10") * &x,
                rounded("1") * &x,
                rounded("-1e10") * &x,
            ],
            vec![
                rounded("0.1") * &x + parse_lit!(y),
                rounded("0.2") * &x - parse_lit!(y),
                rounded("0.3") * &x,
            ],
        ] {
            let expected = terms.iter().fold(Atom::Zero, |sum, term| sum + term);
            assert_eq!(ColorAlgebraSimplifier::sum_rewritten_terms(terms), expected);
        }
        assert_eq!(
            ColorAlgebraSimplifier::sum_rewritten_terms(Vec::new()),
            Atom::Zero
        );
    }

    #[test]
    fn color_products_match_ordered_replacement_with_rounded_coefficients() {
        crate::test_support::test_initialize();
        let rounded = Atom::num(Float::parse("0.1", Some(11)).unwrap());
        let spectator = parse_lit!((x + y) ^ 7);
        for source in [
            parse_lit!(a * b * c) * &spectator,
            &rounded * parse_lit!(a * b * c) * &spectator,
        ] {
            let product = ProductView::parse(source.as_view());
            for replacement in [
                Atom::Zero,
                parse_lit!(a * c ^ 2 / 3),
                &rounded * parse_lit!(a * c ^ 2),
                parse_lit!(a + b),
            ] {
                for left in 0..product.len() {
                    let expected = product.factors.iter().enumerate().fold(
                        Atom::num(1),
                        |value, (index, factor)| {
                            value
                                * if index == left {
                                    replacement.as_view()
                                } else {
                                    factor.atom
                                }
                        },
                    );
                    assert_eq!(product.replacing_one(left, replacement.clone()), expected);
                    for right in left + 1..product.len() {
                        let expected = product
                            .factors
                            .iter()
                            .enumerate()
                            .filter(|(index, _)| *index != left && *index != right)
                            .fold(replacement.clone(), |value, (_, factor)| {
                                value * factor.atom
                            });
                        assert_eq!(
                            product.replacing_pair(left, right, replacement.clone()),
                            expected
                        );
                    }
                }
            }
            for mask in 0..1usize << product.len() {
                let excluded = (0..product.len())
                    .map(|i| mask & (1 << i) != 0)
                    .collect::<Vec<_>>();
                let expected = product
                    .factors
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| !excluded[*i])
                    .fold(Atom::num(1), |value, (_, factor)| value * factor.atom);
                assert_eq!(product.excluding(&excluded, Atom::num(1)), expected);
            }
        }
    }
}
