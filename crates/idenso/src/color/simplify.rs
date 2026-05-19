use std::sync::LazyLock;

use itertools::Itertools;
use spenso::{
    chain,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    rep_,
    shadowing::{self, TensorCollectExt},
    structure::{abstract_index::AIND_SYMBOLS, representation::RepName},
    tensors::parametric::atomcore::PatternReplacement,
    trace, trace_sym,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    function,
    id::{Context, Replacement},
    symbol,
    utils::Settable,
};

use crate::{
    W_, color_f, color_t,
    representations::{ColorAdjoint, ColorFundamental},
    shorthands::{chain::Chain, metric::MetricSimplifier},
};

use super::{CS, ColorSimplifier, ColorSimplifySettings};
use crate::rep_symbols::RS;

static TRACE_TERMINALS: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    [Replacement::new(
        // Empty color trace: Tr_rep(1) -> dim(rep).
        trace!(rep_!(0; W_.d_)).to_pattern(),
        Atom::var(W_.d_),
    )]
});

/// Applies SU(N) simplification rules to raw color tensors and chain/trace form.
pub fn color_simplify_impl(expression: AtomView<'_>) -> Atom {
    color_simplify_with_impl(expression, ColorSimplifySettings::default())
}

/// Applies SU(N) simplification rules with explicit chain/trace settings.
pub fn color_simplify_with_impl(expression: AtomView<'_>, settings: ColorSimplifySettings) -> Atom {
    ColorAlgebraSimplifier {
        settings,
        fundamental_dimension: find_fundamental_dimension(expression),
    }
    .run(expression)
}

struct ColorAlgebraSimplifier {
    settings: ColorSimplifySettings,
    fundamental_dimension: Option<Atom>,
}

impl ColorAlgebraSimplifier {
    fn run(&self, expression: AtomView<'_>) -> Atom {
        let mut current = expression.to_owned().collect_metrics().simplify_metrics();

        loop {
            let next = self.apply_once(current.as_view());
            if next == current {
                return restore_explicit_default_generator_chains(next)
                    .collect_metrics()
                    .simplify_metrics();
            }
            current = next;
        }
    }

    fn apply_once(&self, expression: AtomView<'_>) -> Atom {
        let collected = self.collect_lines(expression);
        self.rewrite_terms(collected.as_view())
            .collect_color()
            .collect_metrics()
            .simplify_metrics()
    }

    fn rewrite_terms(&self, expr: AtomView<'_>) -> Atom {
        // Terminal trace rules can create sums; product rules such as f*f -> CA*g
        // then need to run on each generated term instead of on the whole Add.
        if let AtomView::Add(add) = expr {
            return add
                .iter()
                .map(|term| self.rewrite_terms(term))
                .fold(Atom::Zero, |sum, term| sum + term);
        }

        if let Some(rewritten) = self.rewrite_node(expr) {
            return rewritten;
        }

        // Try product-level rewrites first so trace*f contractions can fire before the
        // trace terminal expands into symmetric trace and f terms.
        expr.to_owned().replace_map(|arg, _context, out| {
            if let Some(rewritten) = self.rewrite_node(arg) {
                **out = rewritten;
            }
        })
    }

    fn collect_lines(&self, expr: AtomView<'_>) -> Atom {
        let rep = ColorFundamental {}.into();
        expr.to_owned()
            .chainify(rep)
            .collect_chains(ColorFundamental {}.into())
            .replace_map(&Self::close_trace_chain)
            .replace_map(&Self::collapse_identity_chain)
    }

    fn close_trace_chain(arg: AtomView, _context: &Context, out: &mut Settable<'_, Atom>) {
        let Some((start, end, factors)) = chain_parts(arg) else {
            return;
        };
        let Some((dimension, start_index, false)) = color_fundamental_slot(start.as_view()) else {
            return;
        };
        let Some((end_dimension, end_index, true)) = color_fundamental_slot(end.as_view()) else {
            return;
        };
        if dimension != end_dimension || start_index != end_index {
            return;
        }

        **out = trace!(ColorFundamental {}.to_symbolic([dimension]); factors);
    }

    fn collapse_identity_chain(arg: AtomView, _context: &Context, out: &mut Settable<'_, Atom>) {
        let Some((start, end, factors)) = chain_parts(arg) else {
            return;
        };
        if factors
            .iter()
            .all(|factor| is_chain_identity_factor(factor.as_view()))
        {
            **out = color_metric(start, end);
        }
    }

    fn rewrite_node(&self, arg: AtomView<'_>) -> Option<Atom> {
        self.simplify_product(arg)
            .or_else(|| Self::simplify_chain_node(arg))
            .or_else(|| {
                self.settings
                    .evaluate_traces
                    .then(|| Self::simplify_trace_node(arg))
                    .flatten()
            })
            .or_else(|| self.simplify_power(arg))
    }

    fn simplify_chain_node(chain: AtomView) -> Option<Atom> {
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

        if let Some(rewritten) =
            Self::simplify_antisymmetric_chain_projector(&start, &end, &factors)
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
                Atom::var(CS.cf) * chain_with_removed_range(&start, &end, &factors, i, i + 2),
            );
        }

        if factors.len() > 2
            && let Some(rewritten) = Self::simplify_separated_chain_casimir(&start, &end, &factors)
        {
            return Some(rewritten);
        }

        None
    }

    fn simplify_trace_node(trace: AtomView) -> Option<Atom> {
        let (rep, factors) = trace_parts(trace)?;

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

        if let Some(rewritten) = Self::simplify_antisymmetric_trace_projector(&rep, &factors) {
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

        let generators = factors
            .iter()
            .map(|factor| color_generator_adjoint(factor.as_view()))
            .collect::<Option<Vec<_>>>()?;

        match generators.as_slice() {
            // Tr(T^a) -> 0.
            [_] => Some(Atom::Zero),
            // Tr(T^a T^b) -> TR g^{ab}.
            [a, b] => Some(Atom::var(CS.tr) * color_metric(a.clone(), b.clone())),
            // Tr(T^a T^b T^c) ->
            //   Tr(sym(T^a,T^b,T^c)) + i/2 TR f^{abc}.
            [a, b, c] => Some(
                color_symmetric_trace(&rep, [a.clone(), b.clone(), c.clone()])
                    + Atom::i() * Atom::num(1) / Atom::num(2)
                        * Atom::var(CS.tr)
                        * color_f!(a.clone(), b.clone(), c.clone()),
            ),
            [a, b, c, d] => Self::simplify_four_generator_trace_terminal(&rep, a, b, c, d),
            _ => None,
        }
    }

    fn simplify_adjacent_trace_casimir(rep: &Atom, factors: &[Atom]) -> Option<Atom> {
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

            // Adjacent equal generators inside a fundamental trace:
            // Tr(... T^a T^a ...) -> CF Tr(...).
            return Some(
                Atom::var(CS.cf)
                    * trace_with_factors(rep.clone(), factors_excluding_range(factors, i, i + 2)),
            );
        }

        None
    }

    fn simplify_antisymmetric_chain_projector(
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
            let x = color_adjoint_dummy_for_pair(a, b)?;

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

    fn simplify_antisymmetric_trace_projector(rep: &Atom, factors: &[Atom]) -> Option<Atom> {
        if let [factor] = factors {
            let (prefactor, args) = color_antisymmetric_generator_args(factor.as_view())?;
            return match args.as_slice() {
                // Tr(antisym(T^a,T^b)) is the trace of a commutator.
                [_, _] => Some(Atom::Zero),
                // Tr(antisym(T^a,T^b,T^c)) -> i/2 TR f^{abc}.
                [a, b, c] => Some(
                    prefactor
                        * Atom::i()
                        * Atom::var(CS.tr)
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
            let x = color_adjoint_dummy_for_pair(a, b)?;

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
        for i in 0..factors.len().saturating_sub(2) {
            let Some(left) = color_generator_adjoint(factors[i].as_view()) else {
                continue;
            };
            let Some(right) = color_generator_adjoint(factors[i + 2].as_view()) else {
                continue;
            };
            if left != right {
                continue;
            }

            // Separated equal generators with one generator between them:
            // Tr(... T^a T^b T^a ...) -> (CF - CA/2) Tr(... T^b ...).
            return Some(
                (Atom::var(CS.cf) - Atom::var(CS.ca) / Atom::num(2))
                    * trace_with_factors(
                        rep.clone(),
                        factors_excluding_indices(factors, &[i, i + 2]),
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
            if left != right {
                continue;
            }

            let coefficient =
                if is_default_fundamental_chain(start, end) && is_default_adjoint_slot(&left) {
                    -Atom::var(CS.tr) / Atom::var(CS.nc)
                } else {
                    Atom::var(CS.cf) - Atom::var(CS.ca) / Atom::num(2)
                };

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
        rep: &Atom,
        a: &Atom,
        b: &Atom,
        c: &Atom,
        d: &Atom,
    ) -> Option<Atom> {
        let x = color_adjoint_dummy_like(a)?;

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
                - Atom::var(CS.tr) / Atom::num(6)
                    * color_f!(a.clone(), c.clone(), x.clone())
                    * color_f!(b.clone(), d.clone(), x.clone())
                + Atom::var(CS.tr) / Atom::num(3)
                    * color_f!(a.clone(), d.clone(), x.clone())
                    * color_f!(b.clone(), c.clone(), x.clone()),
        )
    }

    fn simplify_product(&self, product: AtomView) -> Option<Atom> {
        let product = ProductView::parse(product);
        if product.len() < 2 {
            return None;
        }

        product
            .distribute_color_sum_factor()
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
            .or_else(|| self.simplify_three_f_loop_product(&product))
            .or_else(|| Self::simplify_symmetric_invariant_product(&product))
            .or_else(|| self.simplify_embedded_color_node(&product))
            .or_else(|| {
                self.settings
                    .expand_cross_chain_fierz
                    .then(|| Self::simplify_cross_chain_fierz_product(&product))
                    .flatten()
            })
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

    fn simplify_embedded_color_node(&self, product: &ProductView) -> Option<Atom> {
        for (index, factor) in product.factors.iter().enumerate() {
            let rewritten = Self::simplify_chain_node(factor.atom)
                .or_else(|| {
                    self.settings
                        .evaluate_traces
                        .then(|| Self::simplify_trace_node(factor.atom))
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
        for (left_index, left_factor) in product.factors.iter().enumerate() {
            let Some(left_chain) = &left_factor.chain else {
                continue;
            };
            let Some(left_dimension) =
                fundamental_chain_dimension_view(left_chain.start, left_chain.end)
            else {
                continue;
            };

            for (right_index, right_factor) in
                product.factors.iter().enumerate().skip(left_index + 1)
            {
                let Some(right_chain) = &right_factor.chain else {
                    continue;
                };
                let Some(right_dimension) =
                    fundamental_chain_dimension_view(right_chain.start, right_chain.end)
                else {
                    continue;
                };
                if left_dimension != right_dimension {
                    continue;
                }

                for (left_generator_index, left_generator) in left_chain.factors.iter().enumerate()
                {
                    let Some(left_adjoint) = color_generator_adjoint(*left_generator) else {
                        continue;
                    };
                    for (right_generator_index, right_generator) in
                        right_chain.factors.iter().enumerate()
                    {
                        let Some(right_adjoint) = color_generator_adjoint(*right_generator) else {
                            continue;
                        };
                        if left_adjoint != right_adjoint {
                            continue;
                        }

                        let left_before = &left_chain.factors[..left_generator_index];
                        let left_after = &left_chain.factors[left_generator_index + 1..];
                        let right_before = &right_chain.factors[..right_generator_index];
                        let right_after = &right_chain.factors[right_generator_index + 1..];

                        let crossed_left = chain_with_factor_view_slices(
                            left_chain.start,
                            right_chain.end,
                            &[left_before, right_after],
                        );
                        let crossed_right = chain_with_factor_view_slices(
                            right_chain.start,
                            left_chain.end,
                            &[right_before, left_after],
                        );
                        let uncrossed_left = chain_with_factor_view_slices(
                            left_chain.start,
                            left_chain.end,
                            &[left_before, left_after],
                        );
                        let uncrossed_right = chain_with_factor_view_slices(
                            right_chain.start,
                            right_chain.end,
                            &[right_before, right_after],
                        );

                        let replacement = Atom::var(CS.tr)
                            * (crossed_left * crossed_right
                                - uncrossed_left * uncrossed_right / left_dimension);

                        return Some(product.replacing_pair(left_index, right_index, replacement));
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

        if let Some(invariant) = color_symmetric_invariant(base)
            && invariant.args.len() >= 3
        {
            return Some(color_symmetric_product(
                invariant.args.len(),
                invariant.rep.clone(),
                invariant.rep,
            ));
        }

        let args = structure_constant_args(base)?;
        let dimension = color_structure_dimension(&args)?;
        Some(self.adjoint_casimir() * dimension)
    }

    fn adjoint_casimir(&self) -> Atom {
        match &self.fundamental_dimension {
            Some(dimension) if dimension != &Atom::var(CS.nc) => {
                Atom::num(2) * Atom::var(CS.tr) * dimension
            }
            _ => Atom::var(CS.ca),
        }
    }

    fn simplify_trace_structure_product(&self, product: &ProductView) -> Option<Atom> {
        for (trace_index, trace_factor) in product.factors.iter().enumerate() {
            let Some(trace) = &trace_factor.trace else {
                continue;
            };
            let [first, second, rest @ ..] = trace.factors.as_slice() else {
                continue;
            };
            let Some(a) = color_generator_adjoint_view(*first) else {
                continue;
            };
            let Some(b) = color_generator_adjoint_view(*second) else {
                continue;
            };

            for (f_index, f_factor) in product.factors.iter().enumerate() {
                if f_index == trace_index {
                    continue;
                }
                let Some(structure) = &f_factor.structure else {
                    continue;
                };
                let Some((target, structure_prefactor)) =
                    Self::structure_target_for_generator_pair(&structure.args, &a, &b)
                else {
                    continue;
                };

                // f^{abx} Tr(T^a T^b rest) -> i CA/2 Tr(T^x rest).
                let replacement = structure_prefactor * Atom::i() * self.adjoint_casimir()
                    / Atom::num(2)
                    * trace!(
                        trace.rep.to_owned();
                        std::iter::once(color_t!(target))
                            .chain(rest.iter().map(|factor| factor.to_owned()))
                    );
                return Some(product.replacing_pair(trace_index, f_index, replacement));
            }
        }

        None
    }

    fn simplify_chain_structure_product(&self, product: &ProductView) -> Option<Atom> {
        for (chain_index, chain_factor) in product.factors.iter().enumerate() {
            let Some(chain) = &chain_factor.chain else {
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

                for (f_index, f_factor) in product.factors.iter().enumerate() {
                    if f_index == chain_index {
                        continue;
                    }
                    let Some(structure) = &f_factor.structure else {
                        continue;
                    };
                    let Some((target, structure_prefactor)) =
                        Self::structure_target_for_generator_pair(&structure.args, &left, &right)
                    else {
                        continue;
                    };

                    let coefficient =
                        structure_prefactor * Atom::i() * self.adjoint_casimir() / Atom::num(2);
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
        let actual = color_f!(args[0].to_owned(), args[1].to_owned(), args[2].to_owned());

        args.iter().find_map(|target| {
            let oriented = color_f!(target.to_owned(), left.to_owned(), right.to_owned());
            let prefactor = oriented.coefficient(actual.as_view());
            (!prefactor.is_zero()).then_some((target.to_owned(), prefactor))
        })
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
                let common_count = symmetric
                    .args
                    .iter()
                    .filter(|arg| structure.args.iter().any(|candidate| candidate == *arg))
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

            for (right_index, right_factor) in
                product.factors.iter().enumerate().skip(left_index + 1)
            {
                let Some(right_structure) = &right_factor.structure else {
                    continue;
                };
                let right = right_structure.args.map(|arg| arg.to_owned());
                let Some(replacement) =
                    two_structure_loop_contraction(&left, &right, self.adjoint_casimir())
                else {
                    continue;
                };
                return Some(product.replacing_pair(left_index, right_index, replacement));
            }
        }

        None
    }

    fn simplify_three_f_loop_product(&self, product: &ProductView) -> Option<Atom> {
        for indices in (0..product.factors.len()).combinations(3) {
            let Some(f_args) = indices
                .iter()
                .map(|index| {
                    product.factors[*index]
                        .structure
                        .as_ref()
                        .map(|structure| structure.args.map(|arg| arg.to_owned()))
                })
                .collect::<Option<Vec<_>>>()
            else {
                continue;
            };

            let all_args = f_args.iter().flatten().cloned().collect::<Vec<_>>();
            let externals = all_args
                .iter()
                .filter(|arg| {
                    all_args
                        .iter()
                        .filter(|candidate| *candidate == *arg)
                        .count()
                        == 1
                })
                .cloned()
                .collect::<Vec<_>>();
            if externals.len() != 3 {
                continue;
            }

            let [a, b, c] = externals.as_slice() else {
                continue;
            };
            let replacement =
                self.adjoint_casimir() / Atom::num(2) * color_f!(a.clone(), b.clone(), c.clone());
            let mut excluded = vec![false; product.factors.len()];
            for index in indices {
                excluded[index] = true;
            }
            return Some(product.excluding(&excluded) * replacement);
        }

        None
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
                if left.args.len() != right.args.len() || left.args.len() < 3 {
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
                    let replacement = color_symmetric_product(
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
                let replacement = color_symmetric_product(
                    left_args.len(),
                    left.rep.to_owned(),
                    right.rep.to_owned(),
                ) * color_metric(left_open.clone(), right_open.clone())
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
}

impl<'a> SymmetricInvariantView<'a> {
    fn parse(expr: AtomView<'a>) -> Option<Self> {
        if let Some(trace) = TraceView::parse(expr) {
            let [factor] = trace.factors.as_slice() else {
                return None;
            };
            let args = color_symmetric_trace_arg_views(*factor)?;
            return Some(Self {
                rep: trace.rep,
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
        Some(Self {
            rep: args[0],
            args: args[1..].to_vec(),
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
        self.factors
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != left_index && *index != right_index)
            .fold(replacement, |product, (_, factor)| {
                product * factor.atom.to_owned()
            })
    }

    fn replacing_one(&self, target_index: usize, replacement: Atom) -> Atom {
        self.factors
            .iter()
            .enumerate()
            .fold(Atom::num(1), |product, (index, factor)| {
                if index == target_index {
                    product * replacement.clone()
                } else {
                    product * factor.atom.to_owned()
                }
            })
    }

    fn excluding(&self, excluded: &[bool]) -> Atom {
        self.factors
            .iter()
            .enumerate()
            .filter(|(index, _)| !excluded[*index])
            .fold(Atom::num(1), |product, (_, factor)| {
                product * factor.atom.to_owned()
            })
    }

    fn distribute_color_sum_factor(&self) -> Option<Atom> {
        let (sum_index, sum) = self
            .factors
            .iter()
            .enumerate()
            .find_map(|(index, factor)| match factor.atom {
                AtomView::Add(add) if add.iter().any(atom_contains_color_node) => {
                    Some((index, add))
                }
                _ => None,
            })?;

        Some(
            sum.iter()
                .map(|term| self.replacing_one(sum_index, term.to_owned()))
                .fold(Atom::Zero, |sum, term| sum + term),
        )
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

fn color_generator_adjoint_view(factor: AtomView<'_>) -> Option<AtomView<'_>> {
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

#[derive(Clone, Debug)]
struct ColorSymmetricInvariant {
    rep: Atom,
    args: Vec<Atom>,
}

fn color_symmetric_invariant(invariant: AtomView) -> Option<ColorSymmetricInvariant> {
    if let Some((rep, factors)) = trace_parts(invariant) {
        let [factor] = factors.as_slice() else {
            return None;
        };
        let args = color_symmetric_trace_args(factor.as_view())?;
        return Some(ColorSymmetricInvariant { rep, args });
    }

    let AtomView::Fun(f) = invariant else {
        return None;
    };
    if f.get_symbol() != CS.d || f.get_nargs() < 4 {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    Some(ColorSymmetricInvariant {
        rep: args[0].clone(),
        args: args[1..].to_vec(),
    })
}

fn color_symmetric_trace_args(projector: AtomView) -> Option<Vec<Atom>> {
    color_symmetric_trace_arg_views(projector).map(|args| {
        args.into_iter()
            .map(|arg| arg.to_owned())
            .collect::<Vec<_>>()
    })
}

fn color_symmetric_trace_arg_views(projector: AtomView<'_>) -> Option<Vec<AtomView<'_>>> {
    let AtomView::Fun(f) = projector else {
        return None;
    };
    if f.get_symbol() != *shadowing::SYM {
        return None;
    }

    f.iter()
        .map(color_generator_adjoint_view)
        .collect::<Option<Vec<_>>>()
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

fn color_adjoint_dimension(slot: &Atom) -> Option<Atom> {
    representation_slot(slot.as_view(), CS.adjoint_rep).map(|(dimension, _)| dimension)
}

fn default_adjoint_dimension() -> Atom {
    Atom::var(CS.nc).pow(Atom::num(2)) - Atom::num(1)
}

fn color_structure_dimension(args: &[Atom; 3]) -> Option<Atom> {
    let dimensions = args
        .iter()
        .filter_map(color_adjoint_dimension)
        .dedup()
        .collect::<Vec<_>>();
    match dimensions.as_slice() {
        [] => Some(default_adjoint_dimension()),
        [dimension] => Some(dimension.clone()),
        _ => None,
    }
}

fn color_adjoint_dummy_like(slot: &Atom) -> Option<Atom> {
    let dimension = color_adjoint_dimension(slot)?;
    Some(ColorAdjoint {}.to_symbolic([dimension, Atom::var(CS.trace_dummy)]))
}

fn color_adjoint_dummy_for_pair(left: &Atom, right: &Atom) -> Option<Atom> {
    let left_dimension = color_adjoint_dimension(left)?;
    let right_dimension = color_adjoint_dimension(right)?;
    (left_dimension == right_dimension)
        .then(|| ColorAdjoint {}.to_symbolic([left_dimension, Atom::var(CS.trace_dummy)]))
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

fn find_fundamental_dimension(expression: AtomView<'_>) -> Option<Atom> {
    let generator_pattern =
        color_t!([CS.adj_, RS.a_], [CS.nc_, RS.i_], [CS.nc_, RS.j_]).to_pattern();
    let mut fundamental_dimension = None;

    for matched in expression.pattern_match(&generator_pattern, None, None) {
        let candidate = matched[&CS.nc_].clone();
        if let Some(existing) = &fundamental_dimension {
            if existing != &candidate {
                panic!("Mismatched Nc values in expression")
            }
        } else {
            fundamental_dimension = Some(candidate);
        }
    }

    fundamental_dimension
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
    function!(ETS.metric, left, right)
}

fn color_symmetric_trace(rep: &Atom, factors: impl IntoIterator<Item = Atom>) -> Atom {
    let sym_factors = factors
        .into_iter()
        .map(|factor| color_t!(factor.clone()))
        .collect::<Vec<_>>();
    trace_sym!(rep.clone(); sym_factors)
}

fn color_symmetric_product(rank: usize, left_rep: Atom, right_rep: Atom) -> Atom {
    function!(color_symmetric_product_symbol(rank), left_rep, right_rep)
}

fn color_symmetric_product_symbol(rank: usize) -> Symbol {
    if rank == 3 {
        CS.d33
    } else {
        symbol!(format!("spenso::d{rank}{rank}"))
    }
}

fn fundamental_chain_dimension(start: &Atom, end: &Atom) -> Option<Atom> {
    fundamental_chain_dimension_view(start.as_view(), end.as_view())
}

fn fundamental_chain_dimension_view(start: AtomView<'_>, end: AtomView<'_>) -> Option<Atom> {
    let Some((start_dimension, _, false)) = color_fundamental_slot(start) else {
        return None;
    };
    let Some((end_dimension, _, true)) = color_fundamental_slot(end) else {
        return None;
    };
    (start_dimension == end_dimension).then_some(start_dimension)
}

fn is_default_fundamental_chain(start: &Atom, end: &Atom) -> bool {
    fundamental_chain_dimension(start, end).is_some_and(|dimension| dimension == Atom::var(CS.nc))
}

fn is_default_adjoint_slot(slot: &Atom) -> bool {
    color_adjoint_dimension(slot).is_some_and(|dimension| dimension == default_adjoint_dimension())
}

fn restore_explicit_default_generator_chains(expression: Atom) -> Atom {
    expression.replace_map(|arg, _context, out| {
        let Some((start, end, factors)) = chain_parts(arg) else {
            return;
        };
        let [factor] = factors.as_slice() else {
            return;
        };
        if !is_default_fundamental_chain(&start, &end) {
            return;
        }
        let Some(adjoint) = color_generator_adjoint(factor.as_view()) else {
            return;
        };
        if !is_default_adjoint_slot(&adjoint) {
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
    match expr {
        AtomView::Fun(f) => {
            let symbol = f.get_symbol();
            symbol == CS.f
                || symbol == CS.d
                || symbol == CS.t
                || symbol == T.chain
                || shadowing::trace_parts(f).is_some()
                || f.iter().any(atom_contains_color_node)
        }
        AtomView::Add(add) => add.iter().any(atom_contains_color_node),
        AtomView::Mul(mul) => mul.iter().any(atom_contains_color_node),
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            atom_contains_color_node(base) || atom_contains_color_node(exponent)
        }
        AtomView::Num(_) | AtomView::Var(_) => false,
    }
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
            let common = [left[first.left].clone(), left[second.left].clone()];
            let left_actual = color_f!(left[0].clone(), left[1].clone(), left[2].clone());
            let right_actual = color_f!(right[0].clone(), right[1].clone(), right[2].clone());
            let left_oriented = color_f!(
                left[left_open].clone(),
                common[0].clone(),
                common[1].clone()
            );
            let right_oriented = color_f!(
                right[right_open].clone(),
                common[0].clone(),
                common[1].clone()
            );

            let left_prefactor = left_oriented.coefficient(left_actual.as_view());
            let right_prefactor = right_oriented.coefficient(right_actual.as_view());
            if left_prefactor.is_zero() || right_prefactor.is_zero() {
                return None;
            }
            Some(
                left_prefactor
                    * right_prefactor
                    * adjoint_casimir
                    * color_metric(left[left_open].clone(), right[right_open].clone()),
            )
        }
        [_, _, _] => Some(adjoint_casimir * color_structure_dimension(left)?),
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
