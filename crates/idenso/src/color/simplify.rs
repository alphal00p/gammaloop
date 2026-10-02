use std::{collections::VecDeque, sync::LazyLock};

use spenso::{
    chain,
    network::{
        library::symbolic::ETS,
        parsing::{ParseSettings, ParseState, ShorthandParsing},
        tags::SPENSO_TAG as T,
    },
    rep_,
    shadowing::{self, ProjectorExpander, TensorCollectFilter},
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        partial::PartialStructure,
        representation::{LibraryRep, RepName},
        slot::{DualSlotTo, ParseableAind, SlotMatch, SlotMatcher},
    },
    trace, trace_sym,
};
#[cfg(test)]
use symbolica::function;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    id::Replacement,
};
use symbolica_utils::PatternReplacement;

use crate::{
    W_, color_f, color_t,
    representations::{ColorAdjoint, ColorFundamental, ColorSextet},
    shorthands::chain::Chain,
    tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor, inference::TensorInferenceError},
};

use super::{CS, ColorSimplifier, ColorSimplifySettings};

mod trace;

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
}

impl SymbolicTensor<PartialStructure> {
    pub(crate) fn simplify_color_parts(
        &self,
        settings: ColorSimplifySettings,
        settled: &mut crate::tensor::simplification::observation::SettledRegions,
    ) -> Result<Self, TensorInferenceError> {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::ColorKernel,
        );
        let simplifier = ColorAlgebraSimplifier {
            settings,
            dummies: Self::reserved_dummies([self]),
        };
        self.collect_with_map(
            crate::tensor::CollectionMode::Factored,
            None,
            atom_contains_color_node,
            |selected, complete, _| {
                settled.run(selected, complete, |selected| {
                    // Recognize signed graph zeros before a loop identity opens
                    // them into sums and hides their antisymmetric automorphisms.
                    let mut network = selected
                        .expression
                        .as_view()
                        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                            // Signed-zero detection only needs the visible topology.
                            // Opening a symmetric trace here would enumerate every
                            // projector permutation before the colour kernel runs.
                            shorthand_parsing: ShorthandParsing::Opaque,
                            ..Default::default()
                        })
                        .map_err(|error| TensorInferenceError::Invalid(error.to_string()))?;
                    let expression = if network.remove_antisymmetric_zero_terms() {
                        network
                            .simple_execute::<()>()
                            .map_err(|error| TensorInferenceError::Invalid(error.to_string()))?
                    } else {
                        selected.expression.clone()
                    };
                    let expression = simplifier.step(expression.as_view(), complete);
                    let expression = restore_explicit_su_n_generator_chains(expression);
                    let rewritten = selected.with_identity_result(expression, None)?;
                    // A retained sum can include non-colour coefficients. Only
                    // colour connections are prerequisites of this identity;
                    // the shared contractor preserves the other representations.
                    let mut contracted =
                        rewritten.contract_parts(crate::tensor::ContractSettings {
                            representations: Some(&[
                                ColorFundamental {}.into(),
                                ColorAdjoint {}.into(),
                                ColorSextet {}.into(),
                            ]),
                            collect_chains: false,
                            collect_traces: false,
                            ..Default::default()
                        })?;
                    contracted.root.proofs.frontier = contracted.status;
                    Ok(contracted.root)
                })
            },
        )
    }
}

impl ColorAlgebraSimplifier {
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
        let collected = expression.join_chains(ColorFundamental {}.into());
        let rewritten = self.rewrite_terms(collected.as_view(), complete);
        if self.settings.substitute_cof_dimension_invariants {
            rewritten.to_cof_dimension_invariants()
        } else {
            rewritten
        }
    }

    fn rewrite_terms(&self, expr: AtomView<'_>, complete: bool) -> Atom {
        // Terminal trace rules can create sums; product rules such as f*f -> CA*g
        // then need to run on each generated term instead of on the whole Add.
        if let AtomView::Add(add) = expr {
            let terms = add
                .iter()
                .map(|term| self.rewrite_terms(term, complete))
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

        if let Some(rewritten) = self.rewrite_node(expr, complete) {
            return rewritten;
        }

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
                && let Some(rewritten) = self.rewrite_node(arg, complete)
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

    fn rewrite_node(&self, arg: AtomView<'_>, complete: bool) -> Option<Atom> {
        self.simplify_product(arg, complete)
            .or_else(|| self.simplify_chain_node(arg))
            .or_else(|| {
                self.settings
                    .evaluate_traces
                    .then(|| self.simplify_trace_node(arg, complete))
                    .flatten()
            })
            .or_else(|| self.simplify_power(arg))
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

    fn simplify_trace_node(&self, trace: AtomView, complete: bool) -> Option<Atom> {
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
            return self.simplify_prefixed_generator_trace(&rep, &factors);
        };
        if adjoint {
            return self.simplify_generator_trace(&rep, &generators);
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
            _ => self.simplify_generator_trace(&rep, &generators),
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

    fn simplify_product(&self, product: AtomView, complete: bool) -> Option<Atom> {
        let product = ProductView::parse(product);
        if product.len() < 2 {
            return None;
        }

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
            .or_else(|| self.simplify_two_f_loop_product(&product))
            .or_else(|| self.simplify_adjoint_loop_product(&product))
            .or_else(|| Self::simplify_symmetric_invariant_product(&product))
            .or_else(|| {
                self.settings
                    .expand_cross_chain_fierz
                    .then(|| Self::simplify_cross_chain_fierz_product(&product))
                    .flatten()
            })
            .or_else(|| product.distribute_color_sum_factor())
            .or_else(|| self.simplify_embedded_color_node(&product, complete))
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

    fn simplify_embedded_color_node(&self, product: &ProductView, complete: bool) -> Option<Atom> {
        for (index, factor) in product.factors.iter().enumerate() {
            let rewritten = self
                .simplify_chain_node(factor.atom)
                .or_else(|| {
                    self.settings
                        .evaluate_traces
                        .then(|| self.simplify_trace_node(factor.atom, complete))
                        .flatten()
                })
                .or_else(|| self.simplify_power(factor.atom));
            let Some(rewritten) = rewritten else {
                continue;
            };

            return Some(product.replacing_one(index, rewritten));
        }

        None
    }

    fn simplify_cross_chain_fierz_product(product: &ProductView) -> Option<Atom> {
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
        if positive_integer(exponent)? != 2 {
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

    fn simplify_trace_structure_product(&self, product: &ProductView) -> Option<Atom> {
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
            let generator = |factor| {
                if adjoint {
                    adjoint_generator_slot(trace.rep, factor)
                } else {
                    color_generator_adjoint_view(factor)
                }
            };

            for (f_index, f_factor) in product.factors.iter().enumerate() {
                if f_index == trace_index {
                    continue;
                }
                let Some(structure) = &f_factor.structure else {
                    continue;
                };
                let Some(structure_dimension) =
                    color_structure_dimension(&structure.args.map(|arg| arg.to_owned()))
                else {
                    continue;
                };

                for pair_index in 0..trace.factors.len() {
                    let first = trace.factors[pair_index];
                    let second = trace.factors[(pair_index + 1) % trace.factors.len()];
                    let Some(a) = generator(first) else {
                        continue;
                    };
                    let Some(b) = generator(second) else {
                        continue;
                    };
                    if [a, b].iter().any(|candidate| {
                        slots.iter().filter(|slot| *slot == candidate).count() != 1
                    }) {
                        continue;
                    }
                    let Some((target, structure_prefactor)) =
                        Self::structure_target_for_generator_pair(&structure.args, &a, &b)
                    else {
                        continue;
                    };

                    let rest = (0..trace.factors.len() - 2).map(|offset| {
                        trace.factors[(pair_index + 2 + offset) % trace.factors.len()].to_owned()
                    });
                    let target = if adjoint {
                        color_f!(Atom::var(T.chain_in), Atom::var(T.chain_out), target)
                    } else {
                        color_t!(target)
                    };
                    let replacement_factors =
                        std::iter::once(target).chain(rest).collect::<Vec<_>>();

                    // Fundamental [T,T] = i f T; raw adjoint [F,F] = -f F.
                    let replacement = structure_prefactor
                        * if adjoint { -Atom::one() } else { Atom::i() }
                        * adjoint_casimir_for_dimension(structure_dimension.clone())
                        / Atom::num(2)
                        * trace_with_factors(trace.rep.to_owned(), replacement_factors);
                    return Some(product.replacing_pair(trace_index, f_index, replacement));
                }
            }
        }

        None
    }

    fn simplify_chain_structure_product(&self, product: &ProductView) -> Option<Atom> {
        for (chain_index, chain_factor) in product.factors.iter().enumerate() {
            let Some(chain) = &chain_factor.chain else {
                continue;
            };
            let Some(slots) = chain_factor.generator_slots() else {
                continue;
            };

            for pair_index in 0..chain.factors.len().saturating_sub(1) {
                let Some(left) = color_generator_adjoint_view(chain.factors[pair_index]) else {
                    continue;
                };
                let Some(right) = color_generator_adjoint_view(chain.factors[pair_index + 1])
                else {
                    continue;
                };
                if [left, right]
                    .iter()
                    .any(|candidate| slots.iter().filter(|slot| *slot == candidate).count() != 1)
                {
                    continue;
                }

                for (f_index, f_factor) in product.factors.iter().enumerate() {
                    if f_index == chain_index {
                        continue;
                    }
                    let Some(structure) = &f_factor.structure else {
                        continue;
                    };
                    let Some(structure_dimension) =
                        color_structure_dimension(&structure.args.map(|arg| arg.to_owned()))
                    else {
                        continue;
                    };
                    let Some((target, structure_prefactor)) =
                        Self::structure_target_for_generator_pair(&structure.args, &left, &right)
                    else {
                        continue;
                    };

                    let coefficient = structure_prefactor
                        * Atom::i()
                        * adjoint_casimir_for_dimension(structure_dimension)
                        / Atom::num(2);
                    let chain_factors = chain
                        .factors
                        .iter()
                        .map(|factor| factor.to_owned())
                        .collect::<Vec<_>>();
                    let replacement = coefficient
                        * chain_replacing_factor_pair(
                            &chain.start.to_owned(),
                            &chain.end.to_owned(),
                            &chain_factors,
                            pair_index,
                            color_t!(target.clone()),
                        );
                    return Some(product.replacing_pair(chain_index, f_index, replacement));
                }
            }
        }

        None
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
            }
        }

        None
    }

    fn simplify_two_f_loop_product(&self, product: &ProductView) -> Option<Atom> {
        for (left_index, left_factor) in product.factors.iter().enumerate() {
            let Some(left_structure) = &left_factor.structure else {
                continue;
            };
            let left = left_structure.args.map(|arg| arg.to_owned());
            let Some(dimension) = color_structure_dimension(&left) else {
                continue;
            };

            for (right_index, right_factor) in
                product.factors.iter().enumerate().skip(left_index + 1)
            {
                let Some(right_structure) = &right_factor.structure else {
                    continue;
                };
                let right = right_structure.args.map(|arg| arg.to_owned());
                if color_structure_dimension(&right).as_ref() != Some(&dimension) {
                    continue;
                }
                let Some(replacement) = two_structure_loop_contraction(
                    &left,
                    &right,
                    adjoint_casimir_for_dimension(dimension.clone()),
                ) else {
                    continue;
                };
                return Some(product.replacing_pair(left_index, right_index, replacement));
            }
        }

        None
    }

    fn simplify_adjoint_loop_product(&self, product: &ProductView) -> Option<Atom> {
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
                    && common_structure_positions(left_args, right_args).len() == 1
                {
                    neighbours[left].push(right);
                    neighbours[right].push(left);
                }
            }
        }
        let mut cycle = Vec::new();
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
                        if cycle.is_empty() || left.len() < cycle.len() {
                            cycle = left;
                            if cycle.len() == 3 {
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

        let links = (0..cycle.len())
            .map(|i| {
                common_structure_positions(
                    &structures[cycle[i]].1,
                    &structures[cycle[(i + 1) % cycle.len()]].1,
                )[0]
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
            _ => self.simplify_generator_trace(&adjoint_rep(dimension.clone()), &external)?,
        };
        let mut excluded = vec![false; product.len()];
        for vertex in cycle {
            excluded[structures[vertex].0] = true;
        }
        Some(product.excluding(&excluded) * prefactor * replacement)
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
                if left.args.len() != right.args.len()
                    || left.args.len() < 3
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
                // Contract equal-rank symmetric traces into the corresponding
                // scalar invariant family, leaving a metric for one open pair.
                let (common, left_open, right_open) =
                    symmetric_common_and_open_args(&left_args, &right_args);
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

        let args = f.iter().collect::<Vec<_>>();
        Some(Self {
            args: [args[0], args[1], args[2]],
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
        Atom::mul_many(
            std::iter::once(replacement.as_view()).chain(
                self.factors
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| *index != left_index && *index != right_index)
                    .map(|(_, factor)| factor.atom),
            ),
        )
    }

    fn replacing_one(&self, target_index: usize, replacement: Atom) -> Atom {
        Atom::mul_many(self.factors.iter().enumerate().map(|(index, factor)| {
            if index == target_index {
                replacement.as_view()
            } else {
                factor.atom
            }
        }))
    }

    fn excluding(&self, excluded: &[bool]) -> Atom {
        Atom::mul_many(
            self.factors
                .iter()
                .enumerate()
                .filter(|(index, _)| !excluded[*index])
                .map(|(_, factor)| factor.atom),
        )
    }

    fn distribute_color_sum_factor(&self) -> Option<Atom> {
        let (sum_index, sum) = self
            .factors
            .iter()
            .enumerate()
            .find_map(|(index, factor)| match factor.atom {
                AtomView::Add(add) if add.iter().any(atom_contains_color_node) => {
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
            })?;

        Some(ColorAlgebraSimplifier::sum_rewritten_terms(
            sum.iter()
                .map(|term| self.replacing_one(sum_index, term.to_owned()))
                .collect(),
        ))
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

fn representation_slot(slot: AtomView, symbol: Symbol) -> Option<(Atom, Atom)> {
    let AtomView::Fun(f) = slot else {
        return None;
    };
    if f.get_symbol() != symbol || f.get_nargs() != 2 {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    Some((args[0].clone(), args[1].clone()))
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

fn chain_replacing_factor_pair(
    start: &Atom,
    end: &Atom,
    factors: &[Atom],
    pair_index: usize,
    replacement: Atom,
) -> Atom {
    let mut remaining = factors.to_vec();
    remaining.splice(pair_index..pair_index + 2, [replacement]);
    chain_with_factors(start.clone(), end.clone(), remaining)
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
    let mut slots = SlotMatcher::default();
    let mut selected = false;
    expr.visitor(&mut |node| {
        if selected || !matches!(slots.classify(node), SlotMatch::Other) {
            return false;
        }
        if let AtomView::Fun(function) = node {
            let symbol = function.get_symbol();
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

fn two_structure_loop_contraction(
    left: &[Atom; 3],
    right: &[Atom; 3],
    adjoint_casimir: Atom,
) -> Option<Atom> {
    let common = common_structure_positions(left, right);
    match common.as_slice() {
        [first, second] => {
            let left_open = (0..3).find(|index| *index != first.left && *index != second.left)?;
            let right_open =
                (0..3).find(|index| *index != first.right && *index != second.right)?;
            let left_prefactor = StructureView::orientation(first.left, second.left);
            let right_prefactor = StructureView::orientation(first.right, second.right);
            Some(
                left_prefactor
                    * right_prefactor
                    * adjoint_casimir
                    * color_metric(left[left_open].clone(), right[right_open].clone()),
            )
        }
        [first, second, _] => Some(
            StructureView::orientation(first.left, second.left)
                * StructureView::orientation(first.right, second.right)
                * adjoint_casimir
                * color_structure_dimension(left)?,
        ),
        _ => None,
    }
}

#[derive(Clone, Copy)]
struct CommonStructurePosition {
    left: usize,
    right: usize,
}

fn common_structure_positions(left: &[Atom; 3], right: &[Atom; 3]) -> Vec<CommonStructurePosition> {
    let mut right_used = [false; 3];
    let mut common = Vec::new();

    for (left_index, left_arg) in left.iter().enumerate() {
        if let Some((right_index, _)) = right
            .iter()
            .enumerate()
            .find(|(right_index, right_arg)| !right_used[*right_index] && *right_arg == left_arg)
        {
            right_used[right_index] = true;
            common.push(CommonStructurePosition {
                left: left_index,
                right: right_index,
            });
        }
    }

    common
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
        let simplifier = ColorAlgebraSimplifier {
            settings: Default::default(),
            dummies: ParseState::default(),
        };
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

    #[test]
    fn structure_word_connections_require_unique_ports_in_the_whole_line() {
        crate::test_support::test_initialize();
        let slots = (1..7)
            .map(|i| ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(i)]))
            .collect::<Vec<_>>();
        let incoming = ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(11)]);
        let outgoing =
            spenso::dind!(ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(12)]));
        let simplifier = ColorAlgebraSimplifier {
            settings: Default::default(),
            dummies: ParseState::default(),
        };
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
                        &ProductView::parse(input.as_view())
                    )
                    .is_some(),
                    !repeated,
                );
            }
        }
    }

    // Preserve the prior reconstruction schedule as an independent oracle for
    // rounding and user normalizers, which exact polynomial equality cannot test.
    fn prefix_rewrite(simplifier: &ColorAlgebraSimplifier, expression: AtomView<'_>) -> Atom {
        if let AtomView::Add(sum) = expression {
            return sum.iter().fold(Atom::Zero, |result, term| {
                result + prefix_rewrite(simplifier, term)
            });
        }
        if let Some(result) = simplifier.rewrite_node(expression, true) {
            return result;
        }
        expression.to_owned().replace_map(|node, _context, out| {
            if let Some(result) = simplifier.rewrite_node(node, true) {
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
        let simplifier = ColorAlgebraSimplifier {
            settings: ColorSimplifySettings::default(),
            dummies: ParseState::default(),
        };
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
        let simplifier = ColorAlgebraSimplifier {
            settings: ColorSimplifySettings::default(),
            dummies: ParseState::default(),
        };
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
        let simplifier = ColorAlgebraSimplifier {
            settings: ColorSimplifySettings::default(),
            dummies: ParseState::default(),
        };
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
                assert_eq!(product.excluding(&excluded), expected);
            }
        }
    }
}
