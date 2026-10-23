//! General trace decomposition in the existing colour kernel.
//!
//! A normalized symmetric prefix is multiplied by one ordered generator at a
//! time using equations (30)--(34) of hep-ph/9802376. No projector is expanded
//! into its factorial number of permutations.
//!
//! A line that may still connect to other factors returns after each insertion
//! and leaves the unprocessed suffix in existing trace/projector notation, so
//! the shared planner can apply product rules and contractions before resuming.
//! An isolated line has no such partner: every insertion only creates
//! commutator dummies used later in the same recursion. With
//! [`ColorSimplifySettings::one_shot_traces`] it is decomposed completely in
//! one call.

use symbolica::domains::{integer::Integer, rational::Rational};

use super::*;

/// How far the remaining colour factors of a term can reach a trace line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum LineContext {
    /// Other factors may still connect to the line's ports.
    Unknown,
    /// The line is the only colour factor of a complete top-level term.
    Isolated,
}

/// Fixed data of one symmetric-prefix recursion.
struct TraceRecursion<'a> {
    rep: &'a Atom,
    adjoint: bool,
    /// Bernoulli coefficients up to the largest prefix that can occur.
    coefficients: &'a [Rational],
    /// Continue past the first insertion instead of returning to the planner.
    one_shot: bool,
}

impl ColorAlgebraSimplifier {
    /// Decompose one compatible generator word, preserving raw adjoint phases.
    ///
    /// Fundamental words contain Hermitian T matrices. An adjoint word contains
    /// F^a_bc = f^{bca} = i T_A^a_bc, so [F^a,F^b] = -f^{abc} F^c.
    /// Its emitted symmetric traces also contain raw F matrices: at even degree
    /// n they equal i^n times the normalized Hermitian invariant d_A.
    pub(super) fn simplify_generator_trace(
        &self,
        rep: &Atom,
        generators: &[Atom],
        line: LineContext,
    ) -> Option<Atom> {
        let adjoint = Self::trace_generator_is_adjoint(rep, generators)?;
        let Some(first) = generators.first() else {
            return trace_terminal_dimension(rep.as_view());
        };
        if adjoint
            && let [i1, i2, i3, i4, i5] = generators
            && generators
                .iter()
                .enumerate()
                .all(|(i, slot)| !generators[..i].contains(slot))
        {
            return self.adjoint_five_trace(rep, [i1, i2, i3, i4, i5], line);
        }
        self.trace_line(rep, adjoint, vec![first.clone()], &generators[1..], line)
    }

    /// The five-generator raw adjoint trace, color.h's cOlff(i1,...,i5)
    /// table: four f d_A^{(4)} terms and five f^3 terms, against fifteen
    /// from the symmetric-prefix recursion. d_A^{(4)} is the raw symmetric
    /// trace (i^4 = 1), as in the four-cycle identity.
    fn adjoint_five_trace(
        &self,
        rep: &Atom,
        [i1, i2, i3, i4, i5]: [&Atom; 5],
        line: LineContext,
    ) -> Option<Atom> {
        let k3 = self.color_adjoint_dummy_like(i1)?;
        let k4 = self.color_adjoint_dummy_like(i1)?;
        let quartic = |slots: [&Atom; 4]| trace_sym!(rep.clone(); slots.map(|slot| Self::trace_generator_factor(true, slot)));
        let casimir = quadratic_casimir(rep.clone());
        let result = (quartic([i2, i3, i4, &k3]) * color_f!(&k3, i1, i5)
            - quartic([i1, i4, i5, &k3]) * color_f!(&k3, i2, i3)
            - quartic([i1, i3, i5, &k3]) * color_f!(&k3, i2, i4)
            - quartic([i1, i2, i5, &k3]) * color_f!(&k3, i3, i4))
            / Atom::num(2)
            + casimir / Atom::num(12)
                * (color_f!(i1, &k3, &k4) * color_f!(i2, i3, &k3) * color_f!(i4, i5, &k4)
                    + color_f!(i1, &k3, &k4) * color_f!(i2, i5, &k4) * color_f!(i3, i4, &k3)
                    + color_f!(i1, i3, &k3) * color_f!(i2, i4, &k4) * color_f!(i5, &k3, &k4)
                    - Atom::num(2)
                        * color_f!(i1, i5, &k3)
                        * color_f!(i2, i3, &k4)
                        * color_f!(i4, &k3, &k4)
                    + color_f!(i1, i5, &k3) * color_f!(i2, i4, &k4) * color_f!(i3, &k3, &k4));
        if line == LineContext::Isolated && self.settings.one_shot_traces {
            self.certificate.terminal.set(true);
        }
        self.certificate.expansion.set(true);
        Some(result)
    }

    /// Continue a trace with one symmetric block and an ordered suffix.
    /// Cyclic canonicalization may move the block, so recover its position
    /// instead of treating the first printed factor as a distinguished start.
    pub(super) fn simplify_prefixed_generator_trace(
        &self,
        rep: &Atom,
        factors: &[Atom],
        line: LineContext,
    ) -> Option<Atom> {
        if factors.len() < 2 {
            return None;
        }
        let adjoint = Self::trace_generator_is_adjoint(rep, &[] as &[Atom])?;
        let (position, prefix) = factors.iter().enumerate().find_map(|(position, factor)| {
            let (projector, prefix) = projector_parts(factor.as_view())?;
            (projector == *shadowing::SYM).then_some((position, prefix))
        })?;
        if prefix.is_empty() {
            return None;
        }
        let generator_slot = |factor: &Atom| {
            if adjoint {
                adjoint_generator_slot(rep.as_view(), factor.as_view()).map(|slot| slot.to_owned())
            } else {
                color_generator_adjoint(factor.as_view())
            }
        };
        let prefix = prefix
            .iter()
            .map(generator_slot)
            .collect::<Option<Vec<_>>>()?;
        let suffix = (1..factors.len())
            .map(|offset| generator_slot(&factors[(position + offset) % factors.len()]))
            .collect::<Option<Vec<_>>>()?;
        let generators = prefix.iter().chain(&suffix).cloned().collect::<Vec<_>>();
        Self::trace_generator_is_adjoint(rep, &generators)?;
        self.trace_line(rep, adjoint, prefix, &suffix, line)
    }

    /// Decompose a symmetric prefix and ordered suffix. Repeated ports
    /// contract inside the line, and their open symmetric invariants still
    /// have rules of their own, so only distinct ports qualify for one shot.
    fn trace_line(
        &self,
        rep: &Atom,
        adjoint: bool,
        prefix: Vec<Atom>,
        suffix: &[Atom],
        line: LineContext,
    ) -> Option<Atom> {
        let one_shot = line == LineContext::Isolated && self.settings.one_shot_traces && {
            let mut ports = prefix.iter().chain(suffix).collect::<Vec<_>>();
            ports.sort();
            ports.windows(2).all(|pair| pair[0] != pair[1])
        };
        // The prefix grows by one generator per insertion.
        let coefficients = Self::trace_bernoulli_coefficients(prefix.len() + suffix.len());
        let recursion = TraceRecursion {
            rep,
            adjoint,
            coefficients: &coefficients,
            one_shot,
        };
        let result = self.trace_symmetric_prefix(&recursion, prefix, suffix)?;
        if one_shot {
            self.certificate.terminal.set(true);
        }
        self.certificate.expansion.set(true);
        Some(result)
    }

    /// A structure constant annihilates two legs of the same symmetric block,
    /// even when that block is followed by an unevaluated ordered matrix word.
    pub(super) fn simplify_symmetric_prefix_structure_product(
        product: &ProductView<'_>,
    ) -> Option<Atom> {
        for (trace_index, factor) in product.factors.iter().enumerate() {
            let Some(trace) = &factor.trace else {
                continue;
            };
            // Whole symmetric traces use the existing invariant-product rule.
            if trace.factors.len() < 2 {
                continue;
            }
            // The product factor owns the complete line's slot inventory;
            // equal spellings cannot capture a pair internal to that line.
            let Some(generators) = factor.generator_slots() else {
                continue;
            };
            for &block in &trace.factors {
                let Some(slots) = color_symmetric_trace_arg_views(trace.rep, block) else {
                    continue;
                };
                for (index, other) in product.factors.iter().enumerate() {
                    if index == trace_index {
                        continue;
                    }
                    let Some(structure) = &other.structure else {
                        continue;
                    };
                    if color_structure_dimension(&structure.args.map(|arg| arg.to_owned()))
                        .is_none()
                    {
                        continue;
                    }
                    // Count distinct external legs of f. A slot repeated in
                    // either the prefix or suffix belongs to the trace's scope.
                    let common = structure
                        .args
                        .iter()
                        .filter(|slot| {
                            slots.contains(slot)
                                && generators
                                    .iter()
                                    .filter(|candidate| *candidate == *slot)
                                    .count()
                                    == 1
                        })
                        .count();
                    if common >= 2 {
                        return Some(Atom::Zero);
                    }
                }
            }
        }
        None
    }

    /// Eliminate a repeated generator pair before decomposing an ordered trace.
    /// The identity follows from sum_a T^a X T^a = C_R X -
    /// 1/2 sum_a [T^a,[T^a,X]]. Each term removes two generators.
    pub(super) fn simplify_repeated_generator_trace(
        &self,
        rep: &Atom,
        generators: &[Atom],
    ) -> Option<Atom> {
        if generators.len() <= 2 {
            return None;
        }
        let adjoint = Self::trace_generator_is_adjoint(rep, generators)?;
        // Fundamental adjacent and one-generator gaps retain the established
        // cheap rules. For longer gaps choose the shorter arc of the cycle:
        // the number of commutator terms is binomial(gap, 2), not factorial.
        let first_gap = if adjoint { 0 } else { 2 };
        let (gap, start) = (first_gap..=(generators.len() - 2) / 2).find_map(|gap| {
            (0..generators.len())
                .find(|&start| {
                    generators[start] == generators[(start + gap + 1) % generators.len()]
                })
                .map(|start| (gap, start))
        })?;
        if gap >= 2
            && let Some(block) = Self::simplify_repeated_block_trace(rep, adjoint, generators)
        {
            return Some(block);
        }
        let contracted = &generators[start];
        // The pair was scoped inside the trace. Allocate the replacement
        // connection freshly when moving it outside that notation, rather than
        // capturing a same-spelled index in an unrelated occurrence.
        let bridge = self.color_adjoint_dummy_like(contracted)?;
        let word = (1..generators.len())
            .filter(|&offset| offset != gap + 1)
            .map(|offset| generators[(start + offset) % generators.len()].clone())
            .collect::<Vec<_>>();
        let adjoint_casimir = adjoint_casimir_for_dimension(color_adjoint_dimension(contracted)?);
        let sign = if adjoint { -Atom::one() } else { Atom::one() };
        let mut terms = vec![
            &sign
                * (quadratic_casimir(rep.clone())
                    - Atom::num(Integer::from(gap)) * adjoint_casimir / Atom::num(2))
                * Self::trace_generator_word(rep, adjoint, &word),
        ];
        for left in 0..gap {
            for right in left + 1..gap {
                let a = self.color_adjoint_dummy_like(&word[left])?;
                let b = self.color_adjoint_dummy_like(&word[right])?;
                let mut replaced = word.clone();
                replaced[left] = a.clone();
                replaced[right] = b.clone();
                terms.push(
                    &sign
                        * color_f!(&bridge, &word[left], a)
                        * color_f!(&bridge, &word[right], b)
                        * Self::trace_generator_word(rep, adjoint, &replaced),
                );
            }
        }
        Some(Atom::add_many(terms))
    }

    /// color.h's two-block identity: a repeated ordered pair a b becomes
    /// reversed, and the commutator closes the pair into one bridge,
    ///   Tr(X^a X^b W X^a X^b R) = Tr(X^a X^b W X^b X^a R) + κ² (C_A/2) Tr(X^a W X^a R),
    /// with [X^a, X^b] = κ f^{abk} X^k (κ² = -1 for T, +1 for raw F), since
    /// f^{abk} X^a X^b = (κ/2) C_A X^k. Two terms against binomial(gap, 2) + 1
    /// from the generic pair rule; the reversed pair is closer.
    fn simplify_repeated_block_trace(
        rep: &Atom,
        adjoint: bool,
        generators: &[Atom],
    ) -> Option<Atom> {
        let n = generators.len();
        let twice = |slot: &Atom| generators.iter().filter(|other| *other == slot).count() == 2;
        let (first, second) = (0..n).find_map(|first| {
            let (a, b) = (&generators[first], &generators[(first + 1) % n]);
            if a == b || !twice(a) || !twice(b) {
                return None;
            }
            (2..n - 1)
                .map(|offset| (first + offset) % n)
                .find(|&second| generators[second] == *a && generators[(second + 1) % n] == *b)
                .map(|second| (first, second))
        })?;
        let a = &generators[first];
        let b = &generators[(first + 1) % n];
        let between = (first + 2..)
            .map(|position| position % n)
            .take_while(|&position| position != second)
            .map(|position| generators[position].clone())
            .collect::<Vec<_>>();
        let rest = (second + 2..)
            .map(|position| position % n)
            .take_while(|&position| position != first)
            .map(|position| generators[position].clone())
            .collect::<Vec<_>>();
        let word = |pair: [&Atom; 2], closing: [&Atom; 2]| {
            pair.into_iter()
                .cloned()
                .chain(between.iter().cloned())
                .chain(closing.into_iter().cloned())
                .chain(rest.iter().cloned())
                .collect::<Vec<_>>()
        };
        let reversed = word([a, b], [b, a]);
        let bridged = std::iter::once(a.clone())
            .chain(between.iter().cloned())
            .chain(std::iter::once(a.clone()))
            .chain(rest.iter().cloned())
            .collect::<Vec<_>>();
        let kappa_squared = if adjoint { Atom::one() } else { -Atom::one() };
        let casimir = adjoint_casimir_for_dimension(color_adjoint_dimension(a)?);
        Some(
            Self::trace_generator_word(rep, adjoint, &reversed)
                + kappa_squared * casimir / Atom::num(2)
                    * Self::trace_generator_word(rep, adjoint, &bridged),
        )
    }

    fn trace_generator_word(rep: &Atom, adjoint: bool, generators: &[Atom]) -> Atom {
        trace_with_factors(
            rep.clone(),
            generators
                .iter()
                .map(|slot| Self::trace_generator_factor(adjoint, slot))
                .collect(),
        )
    }

    pub(super) fn trace_generator_factor(adjoint: bool, slot: &Atom) -> Atom {
        if adjoint {
            color_f!(Atom::var(T.chain_in), Atom::var(T.chain_out), slot)
        } else {
            color_t!(slot)
        }
    }

    fn trace_prefixed_word(
        rep: &Atom,
        adjoint: bool,
        mut prefix: Vec<Atom>,
        suffix: &[Atom],
    ) -> Atom {
        if suffix.len() <= 1 {
            prefix.extend_from_slice(suffix);
            return Self::trace_symmetric_terminal(rep, adjoint, prefix);
        }
        let symmetric = shadowing::sym(
            prefix
                .iter()
                .map(|slot| Self::trace_generator_factor(adjoint, slot)),
        );
        trace_with_factors(
            rep.clone(),
            std::iter::once(symmetric)
                .chain(
                    suffix
                        .iter()
                        .map(|slot| Self::trace_generator_factor(adjoint, slot)),
                )
                .collect(),
        )
    }

    pub(super) fn trace_generator_is_adjoint(
        rep: &Atom,
        generators: &[impl AtomCore],
    ) -> Option<bool> {
        let AtomView::Fun(representation) = rep.as_view() else {
            return None;
        };
        if representation.get_nargs() != 1 {
            return None;
        }
        let adjoint = representation.get_symbol() == CS.adjoint_rep;
        if !adjoint && representation.get_symbol() != CS.fundamental_rep {
            return None;
        }
        if !generators.is_empty() {
            let dimension = color_structure_dimension(generators)?;
            if adjoint && representation.iter().next()? != dimension.as_view() {
                return None;
            }
        }
        Some(adjoint)
    }

    /// Coefficients of x/(1-exp(-x)): (-1)^j B_j/j!, including +1/2 at j=1.
    /// Symbolica's Bernoulli helper is private; use its existing arbitrary-size
    /// rational arithmetic for the triangular recurrence and factorial here.
    fn trace_bernoulli_coefficients(degree: usize) -> Vec<Rational> {
        let mut table = Vec::<Rational>::with_capacity(degree + 1);
        let mut coefficients = Vec::with_capacity(degree + 1);
        let mut factorial = Rational::one();
        for n in 0..=degree {
            table.push(Rational::one() / Rational::from(Integer::from(n + 1)));
            for j in (1..=n).rev() {
                table[j - 1] = (&table[j - 1] - &table[j]) * Rational::from(Integer::from(j));
            }
            if n != 0 {
                factorial *= Rational::from(Integer::from(n));
            }
            coefficients.push(&table[0] / &factorial);
        }
        coefficients
    }

    fn trace_symmetric_prefix(
        &self,
        recursion: &TraceRecursion<'_>,
        mut prefix: Vec<Atom>,
        suffix: &[Atom],
    ) -> Option<Atom> {
        let TraceRecursion {
            rep,
            adjoint,
            coefficients,
            ..
        } = *recursion;
        // Cyclicity makes Tr(sym(A) B) fully symmetric, so the last generator
        // needs no commutator expansion. In particular this avoids allocating
        // dummy indices in terms which would vanish by trace cyclicity.
        if suffix.len() <= 1 {
            prefix.extend_from_slice(suffix);
            return Some(Self::trace_symmetric_terminal(rep, adjoint, prefix));
        }
        prefix.sort();
        if let ([a, b], [c, d]) = (prefix.as_slice(), suffix) {
            // Averaging the ordered four-generator identity over a,b gives
            // color.h's two-plus-two terminal rule. It avoids intermediate
            // rank-two traces/metrics and combines the two cubic invariants.
            // For raw F=i*T_A, the fourth-order phase is positive: the f*f
            // coefficient is +CA/12, while the odd symmetric trace vanishes.
            let x = self.color_adjoint_dummy_like(a)?;
            let index = if adjoint {
                quadratic_casimir(rep.clone())
            } else {
                quadratic_index(rep.clone())
            };
            let mut result = Self::trace_symmetric_terminal(
                rep,
                adjoint,
                vec![a.clone(), b.clone(), c.clone(), d.clone()],
            ) + index / Atom::num(12)
                * (color_f!(a, c, &x) * color_f!(b, d, &x)
                    + color_f!(a, d, &x) * color_f!(b, c, &x));
            if !adjoint {
                result += Atom::i() / Atom::num(2)
                    * color_f!(c, d, &x)
                    * color_symmetric_trace(rep, [a.clone(), b.clone(), x]);
            }
            return Some(result);
        }
        let (next, remaining) = suffix.split_first()?;
        let mut terms = Vec::new();
        for (depth, coefficient) in coefficients.iter().enumerate().take(prefix.len() + 1) {
            if coefficient.is_zero() {
                continue;
            }
            // With at most one remaining generator, cyclicity makes every
            // branch fully symmetric. Its degree is known before allocating
            // the commutator tree; odd raw-adjoint invariants vanish exactly.
            if adjoint
                && remaining.len() <= 1
                && (prefix.len() + 1 - depth + remaining.len()) % 2 == 1
            {
                continue;
            }
            let term = self.trace_prefix_commutators(recursion, &prefix, next, depth, remaining)?;
            if !term.is_zero() {
                terms.push(Atom::num(coefficient.clone()) * term);
            }
        }
        Some(Atom::add_many(terms))
    }

    fn trace_prefix_commutators(
        &self,
        recursion: &TraceRecursion<'_>,
        prefix: &[Atom],
        next: &Atom,
        depth: usize,
        suffix: &[Atom],
    ) -> Option<Atom> {
        let adjoint = recursion.adjoint;
        if depth == 0 {
            let mut symmetric = prefix.to_vec();
            symmetric.push(next.clone());
            // Nothing between two insertions of an isolated line is
            // contractible: each dummy is used only later in this recursion.
            return if recursion.one_shot {
                self.trace_symmetric_prefix(recursion, symmetric, suffix)
            } else {
                Some(Self::trace_prefixed_word(
                    recursion.rep,
                    adjoint,
                    symmetric,
                    suffix,
                ))
            };
        }
        let mut terms = Vec::new();
        for (position, generator) in prefix.iter().enumerate() {
            if generator == next || (position != 0 && prefix[position - 1] == *generator) {
                continue;
            }
            // Equal entries in a symmetric prefix give the same term. Account
            // for their multiplicity without repeating the whole recursion.
            let multiplicity = prefix[position..]
                .iter()
                .take_while(|candidate| *candidate == generator)
                .count();
            let contracted = self.color_adjoint_dummy_like(generator)?;
            let mut remaining = prefix.to_vec();
            let _ = remaining.remove(position);
            let nested = self.trace_prefix_commutators(
                recursion,
                &remaining,
                &contracted,
                depth - 1,
                suffix,
            )?;
            if !nested.is_zero() {
                let commutator = if adjoint { -Atom::one() } else { Atom::i() };
                // A rank-two terminal leaves the new dummy in a metric only:
                // f^{g n c} g^{c y} = f^{g n y}, with nothing else to contract.
                let (nested, contracted) =
                    Self::metric_partner(&nested, &contracted).unwrap_or((nested, contracted));
                terms.push(
                    Atom::num(Integer::from(multiplicity))
                        * commutator
                        * color_f!(generator, next, contracted)
                        * nested,
                );
            }
        }
        Some(Atom::add_many(terms))
    }

    /// `term` as `rest * g(dummy, partner)`, with no other occurrence of
    /// `dummy`: the rank-two terminal of a commutator branch.
    fn metric_partner(term: &Atom, dummy: &Atom) -> Option<(Atom, Atom)> {
        let factors = multiplicative_factor_views(term.as_view());
        let (position, partner) = factors.iter().enumerate().find_map(|(position, factor)| {
            let [x, y] = emitted_adjoint_metric(*factor)?;
            if x == *dummy {
                Some((position, y))
            } else if y == *dummy {
                Some((position, x))
            } else {
                None
            }
        })?;
        let rest = factors
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != position)
            .map(|(_, factor)| *factor)
            .collect::<Vec<_>>();
        if rest.iter().any(|factor| factor.contains(dummy.as_view())) {
            return None;
        }
        Some((Atom::mul_many(rest), partner))
    }

    fn trace_symmetric_terminal(rep: &Atom, adjoint: bool, generators: Vec<Atom>) -> Atom {
        match generators.as_slice() {
            [] => trace!(rep.clone()),
            [_] => Atom::Zero,
            [left, right] => {
                let index = if adjoint {
                    -quadratic_casimir(rep.clone())
                } else {
                    quadratic_index(rep.clone())
                };
                index * color_metric(left.clone(), right.clone())
            }
            _ if adjoint && generators.len() % 2 == 1 => Atom::Zero,
            _ if adjoint => trace_sym!(rep.clone(); generators.iter().map(|slot| {
                Self::trace_generator_factor(true, slot)
            })),
            _ => color_symmetric_trace(rep, generators),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A rank-two terminal of a commutator branch contracts its metric into
    /// the branch's structure constant: a one-shot decomposition leaves no
    /// metric with a dummy end, so the planner has nothing to contract.
    #[test]
    fn rank_two_terminals_fold_into_structure_constants() {
        crate::test_support::test_initialize();
        let label = spenso::index_symbol!("idenso::terminal_fold::edge");
        for (adjoint, length) in [(false, 5), (false, 6), (true, 5), (true, 6)] {
            let slots = (0..length)
                .map(|i| ColorAdjoint {}.to_symbolic([Atom::num(8), function!(label, i, 1)]))
                .collect::<Vec<_>>();
            let rep = if adjoint {
                ColorAdjoint {}.to_symbolic([Atom::num(8)])
            } else {
                ColorFundamental {}.to_symbolic([Atom::num(3)])
            };
            let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
            let recursion = TraceRecursion {
                rep: &rep,
                adjoint,
                coefficients: &ColorAlgebraSimplifier::trace_bernoulli_coefficients(length),
                one_shot: true,
            };
            let result = simplifier
                .trace_symmetric_prefix(&recursion, slots[..1].to_vec(), &slots[1..])
                .unwrap();
            let mut dummy_metrics = 0;
            result.visitor(&mut |node| {
                if let Some([x, y]) = emitted_adjoint_metric(node)
                    && (!slots.contains(&x) || !slots.contains(&y))
                {
                    dummy_metrics += 1;
                }
                true
            });
            assert_eq!(
                dummy_metrics, 0,
                "adjoint={adjoint} length={length}: {result}"
            );
        }
    }

    #[test]
    fn exact_bernoulli_coefficients_have_no_machine_integer_factorial_limit() {
        let coefficients = ColorAlgebraSimplifier::trace_bernoulli_coefficients(20);
        for (degree, value) in [
            (0, "1"),
            (1, "1/2"),
            (2, "1/12"),
            (4, "-1/720"),
            (6, "1/30240"),
            (12, "-691/1307674368000"),
            (20, "-174611/802857662698291200000"),
        ] {
            assert_eq!(
                Atom::num(coefficients[degree].clone()),
                Atom::parse(value, "idenso", Default::default()).unwrap()
            );
        }
        for degree in (3..20).step_by(2) {
            assert!(coefficients[degree].is_zero());
        }
    }

    #[test]
    fn cubic_generator_traces_preserve_fundamental_and_raw_adjoint_conventions() {
        let reps = crate::test_support::test_initialize();
        let slots = [1, 2, 3].map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        let fundamental = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let adjoint = reps.coad_da.to_symbolic([] as [Atom; 0]);
        let mut reversed = slots.clone();
        reversed.swap(0, 1);
        for word in [slots, reversed] {
            let f = color_f!(&word[0], &word[1], &word[2]);
            for (rep, expected) in [
                (
                    &fundamental,
                    color_symmetric_trace(&fundamental, word.clone())
                        + Atom::i() / Atom::num(2) * quadratic_index(fundamental.clone()) * &f,
                ),
                (
                    &adjoint,
                    quadratic_casimir(adjoint.clone()) / Atom::num(2) * &f,
                ),
            ] {
                let result = simplifier
                    .simplify_generator_trace(rep, &word, LineContext::Unknown)
                    .unwrap();
                let contracted = SymbolicTensor::infer(result)
                    .unwrap()
                    .contract(crate::tensor::ContractSettings {
                        collect_chains: false,
                        collect_traces: false,
                        ..Default::default()
                    })
                    .unwrap();
                assert!(
                    (contracted.expression - expected)
                        .collect_factors()
                        .factor()
                        .is_zero()
                );
            }
        }
    }

    #[test]
    fn antisymmetric_structure_annihilates_only_shared_symmetric_prefix_legs() {
        crate::test_support::test_initialize();
        let label = spenso::index_symbol!("idenso::symmetric_prefix_zero::edge");
        let slots = (0..5)
            .map(|index| {
                let index = function!(label, Atom::num(index), Atom::num(0));
                let admitted = AbstractIndex::from_view(index.as_view()).unwrap();
                assert_eq!(admitted.to_atom(), index);
                ColorAdjoint {}.to_symbolic([Atom::num(8), index])
            })
            .collect::<Vec<_>>();
        for (rep, adjoint) in [
            (fundamental_rep(Atom::num(3)), false),
            (adjoint_rep(Atom::num(8)), true),
        ] {
            let block = shadowing::sym(
                slots[..2]
                    .iter()
                    .map(|slot| ColorAlgebraSimplifier::trace_generator_factor(adjoint, slot)),
            );
            let trace = trace_with_factors(
                rep,
                vec![
                    block,
                    ColorAlgebraSimplifier::trace_generator_factor(adjoint, &slots[2]),
                    ColorAlgebraSimplifier::trace_generator_factor(adjoint, &slots[3]),
                ],
            );
            for (second, annihilates) in [(1, true), (2, false)] {
                let source = &trace * color_f!(&slots[0], &slots[second], &slots[4]);
                SymbolicTensor::infer(source.clone()).unwrap();
                let product = ProductView::parse(source.as_view());
                assert_eq!(
                    ColorAlgebraSimplifier::simplify_symmetric_prefix_structure_product(&product),
                    annihilates.then_some(Atom::Zero),
                );
            }
            for (prefix, suffix) in [(vec![0, 0, 1], vec![2]), (vec![0, 1], vec![0, 2])] {
                let block = shadowing::sym(prefix.iter().map(|&index| {
                    ColorAlgebraSimplifier::trace_generator_factor(adjoint, &slots[index])
                }));
                let scoped = trace_with_factors(
                    if adjoint {
                        adjoint_rep(Atom::num(8))
                    } else {
                        fundamental_rep(Atom::num(3))
                    },
                    std::iter::once(block)
                        .chain(suffix.iter().map(|&index| {
                            ColorAlgebraSimplifier::trace_generator_factor(adjoint, &slots[index])
                        }))
                        .collect(),
                );
                SymbolicTensor::infer(scoped.clone()).unwrap();
                let source = scoped * color_f!(&slots[0], &slots[1], &slots[4]);
                assert!(
                    ColorAlgebraSimplifier::simplify_symmetric_prefix_structure_product(
                        &ProductView::parse(source.as_view())
                    )
                    .is_none(),
                    "an internal trace pair is not a boundary connection: {source}",
                );
            }
            // A different adjoint space with the same index spellings does not
            // contract with the prefix's slots.
            let foreign = slots[..2]
                .iter()
                .map(|slot| {
                    let (_, index) = representation_slot(slot.as_view(), CS.adjoint_rep).unwrap();
                    ColorAdjoint {}.to_symbolic([Atom::num(7), index])
                })
                .collect::<Vec<_>>();
            let source = trace
                * color_f!(
                    &foreign[0],
                    &foreign[1],
                    ColorAdjoint {}.to_symbolic([Atom::num(7), Atom::num(99)])
                );
            SymbolicTensor::infer(source.clone()).unwrap();
            assert!(
                ColorAlgebraSimplifier::simplify_symmetric_prefix_structure_product(
                    &ProductView::parse(source.as_view())
                )
                .is_none()
            );
        }
    }

    #[test]
    fn generator_trace_leaves_exact_prefix_continuations() {
        let reps = crate::test_support::test_initialize();
        let slots = (0..8)
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]))
            .collect::<Vec<_>>();
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let longest_prefix = |result: &Atom| {
            let mut longest = 0;
            result.as_view().visitor(&mut |node| {
                if let AtomView::Fun(function) = node
                    && function.get_symbol() == *shadowing::SYM
                {
                    longest = longest.max(function.get_nargs());
                }
                true
            });
            longest
        };
        let result = simplifier
            .simplify_generator_trace(&rep, &slots, LineContext::Unknown)
            .unwrap();
        let AtomView::Add(sum) = result.as_view() else {
            panic!("one insertion must leave its symmetric and commutator branches");
        };
        assert_eq!(sum.get_nargs(), 2);
        assert_eq!(longest_prefix(&result), 2);

        // An isolated line has no partner to wait for; one shot leaves no
        // ordered continuation, only complete symmetric invariants.
        let isolated = simplifier
            .simplify_generator_trace(&rep, &slots, LineContext::Isolated)
            .unwrap();
        assert_eq!(ordered_trace_words(&isolated), 0);
        assert_eq!(longest_prefix(&isolated), slots.len());
        // Without one-shot traces, an isolated line also resumes per insertion.
        let incremental = ColorAlgebraSimplifier::new(
            ColorSimplifySettings::default().without_one_shot_traces(),
            ParseState::default(),
        )
        .simplify_generator_trace(&rep, &slots, LineContext::Isolated)
        .unwrap();
        assert_eq!(longest_prefix(&incremental), 2);
    }

    fn ordered_trace_words(result: &Atom) -> usize {
        let mut words = 0;
        result.as_view().visitor(&mut |node| {
            if let AtomView::Fun(function) = node
                && let Some((_, factors)) = shadowing::trace_parts(function)
            {
                words += usize::from(factors.len() > 1);
            }
            true
        });
        words
    }

    /// Fully distributed term count, without expanding the nested result.
    fn flat_terms(view: AtomView<'_>) -> usize {
        match view {
            AtomView::Add(sum) => sum.iter().map(flat_terms).sum(),
            AtomView::Mul(product) => product.iter().map(flat_terms).product(),
            _ => 1,
        }
    }

    #[test]
    fn one_shot_decomposition_matches_the_resumed_continuations() {
        let reps = crate::test_support::test_initialize();
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        // Leaves of the bare recursion for ranks 5, 6 and 7. The planner
        // merges no terms between insertions, so both routes give these.
        // Fundamental counts follow the two-plus-two terminal of 164274b4b.
        for (adjoint, census) in [(false, [28, 179, 1432]), (true, [15, 109, 798])] {
            let rep = if adjoint {
                reps.coad_da.to_symbolic([] as [Atom; 0])
            } else {
                reps.cof_nc.to_symbolic([] as [Atom; 0])
            };
            for (rank, expected) in (5..).zip(census) {
                let slots = (0..rank)
                    .map(|index| reps.coad_da.to_symbolic([Atom::num(99_870 + index)]))
                    .collect::<Vec<_>>();
                let coefficients = ColorAlgebraSimplifier::trace_bernoulli_coefficients(rank);
                let one_shot = simplifier
                    .trace_symmetric_prefix(
                        &TraceRecursion {
                            rep: &rep,
                            adjoint,
                            coefficients: &coefficients,
                            one_shot: true,
                        },
                        vec![slots[0].clone()],
                        &slots[1..],
                    )
                    .unwrap();
                assert_eq!(ordered_trace_words(&one_shot), 0);
                assert_eq!(
                    flat_terms(one_shot.as_view()),
                    expected,
                    "adjoint={adjoint} rank={rank}"
                );

                // Resume each continuation as the planner does, without the
                // contractions it applies between rounds. The five-generator
                // adjoint trace takes color.h's table instead: nine terms.
                let mut resumed = simplifier
                    .simplify_generator_trace(&rep, &slots, LineContext::Unknown)
                    .unwrap();
                if adjoint && rank == 5 {
                    assert_eq!(ordered_trace_words(&resumed), 0);
                    assert_eq!(flat_terms(resumed.as_view()), 9);
                    continue;
                }
                while ordered_trace_words(&resumed) != 0 {
                    resumed = resumed.replace_map(|node, _context, out| {
                        if let Some(result) =
                            simplifier.simplify_trace_node(node, true, LineContext::Unknown)
                        {
                            **out = result;
                        }
                    });
                }
                assert_eq!(
                    flat_terms(resumed.as_view()),
                    expected,
                    "adjoint={adjoint} rank={rank}"
                );
            }
        }
    }

    #[test]
    fn two_plus_two_terminal_matches_symmetrized_ordered_quartic_identity() {
        use crate::IndexTooling;

        let reps = crate::test_support::test_initialize();
        let label = spenso::index_symbol!("quartic_prefix::edge");
        let slots = [0, 1, 2, 3].map(|axis| {
            reps.coad_da
                .to_symbolic([function!(label, Atom::num(axis), Atom::num(1))])
        });
        let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        let coefficients = ColorAlgebraSimplifier::trace_bernoulli_coefficients(2);
        let result = simplifier
            .trace_symmetric_prefix(
                &TraceRecursion {
                    rep: &rep,
                    adjoint: false,
                    coefficients: &coefficients,
                    one_shot: false,
                },
                slots[..2].to_vec(),
                &slots[2..],
            )
            .unwrap();
        // The established ordered-word identity supplies an independent
        // symbolic derivation; only its first two matrices are symmetrized.
        let first = simplifier
            .simplify_four_generator_trace_terminal(
                &rep, &slots[0], &slots[1], &slots[2], &slots[3],
            )
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        let second = simplifier
            .simplify_four_generator_trace_terminal(
                &rep, &slots[1], &slots[0], &slots[2], &slots[3],
            )
            .unwrap()
            .canonize(AbstractIndex::Dummy)
            .unwrap();
        assert!(
            (result.canonize(AbstractIndex::Dummy).unwrap() - (first + second) / Atom::num(2))
                .collect_factors()
                .factor()
                .is_zero()
        );
        assert!(!result.contains_symbol(ETS.metric));
        assert!(result.contains_symbol(CS.idx));
        assert!(result.contains_symbol(T.trace));
    }

    #[test]
    fn prefixed_trace_terminal_uses_cyclicity_and_adjoint_parity() {
        let reps = crate::test_support::test_initialize();
        let slots = [1, 2, 3].map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        for (rep, adjoint) in [
            (reps.cof_nc.to_symbolic([] as [Atom; 0]), false),
            (reps.coad_da.to_symbolic([] as [Atom; 0]), true),
        ] {
            let mut factors = vec![
                shadowing::sym(
                    slots[..2]
                        .iter()
                        .map(|slot| ColorAlgebraSimplifier::trace_generator_factor(adjoint, slot)),
                ),
                ColorAlgebraSimplifier::trace_generator_factor(adjoint, &slots[2]),
            ];
            let expected = if adjoint {
                Atom::Zero
            } else {
                color_symmetric_trace(&rep, slots.clone())
            };
            for _ in 0..factors.len() {
                assert_eq!(
                    simplifier.simplify_prefixed_generator_trace(
                        &rep,
                        &factors,
                        LineContext::Unknown
                    ),
                    Some(expected.clone()),
                );
                factors.rotate_left(1);
            }
        }
    }

    #[test]
    fn odd_adjoint_terminal_branches_allocate_no_commutator_dummies() {
        let reps = crate::test_support::test_initialize();
        let slots = (0..5)
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]))
            .collect::<Vec<_>>();
        // Isolate the second Bernoulli term: three prefix generators and two
        // suffix generators become a symmetric trace of degree three. The
        // whole raw-adjoint branch is zero before any dummy can be allocated.
        let coefficients = [
            Rational::zero(),
            Rational::zero(),
            Rational::one() / Rational::from(Integer::from(12)),
        ];
        for (rep, adjoint) in [
            (reps.coad_da.to_symbolic([] as [Atom; 0]), true),
            (reps.cof_nc.to_symbolic([] as [Atom; 0]), false),
        ] {
            let simplifier = ColorAlgebraSimplifier::new(
                ColorSimplifySettings::default(),
                ParseState::default(),
            );
            // The existing diagnostic view includes the operation-local
            // reservation set; unlike the global dummy counter this is stable
            // when unrelated tests allocate their own fresh indices.
            let before = format!("{:?}", simplifier.dummies);
            let recursion = TraceRecursion {
                rep: &rep,
                adjoint,
                coefficients: &coefficients,
                one_shot: false,
            };
            let result = simplifier
                .trace_symmetric_prefix(&recursion, slots[..3].to_vec(), &slots[3..])
                .unwrap();
            assert_eq!(result.is_zero(), adjoint);
            assert_eq!(format!("{:?}", simplifier.dummies) == before, adjoint);
        }
    }

    #[test]
    fn repeated_generator_pairs_match_independent_closed_su2_components() {
        crate::test_support::test_initialize();
        let slots =
            [1, 2, 3, 4].map(|index| ColorAdjoint {}.to_symbolic([Atom::num(3), Atom::num(index)]));
        let word = slots.iter().chain(&slots).cloned().collect::<Vec<_>>();
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        // For T=Pauli/2, Str(Ta Tb Tc Td)=S_abcd/24; for the real
        // Levi-Civita matrices, Str(Fa Fb Fc Fd)=2*S_abcd/3, where
        // S_abcd=delta_ab*delta_cd+delta_ac*delta_bd+delta_ad*delta_bc.
        // Contract these independent component formulas over all 3^4 slots.
        let mut gram_numerator = 0;
        for a in 0..3 {
            for b in 0..3 {
                for c in 0..3 {
                    for d in 0..3 {
                        let symmetric = i64::from(a == b && c == d)
                            + i64::from(a == c && b == d)
                            + i64::from(a == d && b == c);
                        gram_numerator += symmetric * symmetric;
                    }
                }
            }
        }
        let mixed_quartic_value = Atom::num((gram_numerator, 36));
        assert_eq!(mixed_quartic_value, Atom::num((5, 4)));
        // Independent finite sums over all 3^4 Pauli/2 or real Levi-Civita
        // matrix words give -39/128 and 18, respectively. No simplifier was
        // used to obtain either reference value.
        for (rep, expected) in [
            (fundamental_rep(Atom::num(2)), Atom::num((-39, 128))),
            (adjoint_rep(Atom::num(3)), Atom::num(18)),
        ] {
            let shortened = simplifier
                .simplify_repeated_generator_trace(&rep, &word)
                .unwrap();
            let reduced = SymbolicTensor::infer(shortened)
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                    ..Default::default()
                })
                .unwrap();
            // The general result retains the mixed quartic invariant. Its
            // independent SU(2) value is 5/4: sum the 81 component products of
            // normalized symmetric Pauli/2 and Levi-Civita traces. This is an
            // oracle specialization, not another algebra-kernel rewrite.
            let mixed_quartic = CS.gram(
                Atom::num(4),
                adjoint_rep(Atom::num(3)),
                fundamental_rep(Atom::num(2)),
            );
            let evaluated = reduced
                .expression
                .replace(mixed_quartic.to_pattern())
                .with(mixed_quartic_value.clone());
            assert!(
                (&evaluated - &expected)
                    .collect_factors()
                    .factor()
                    .is_zero(),
                "repeated-pair trace in {rep}: actual={}, expected={expected}, status={:?}",
                reduced.expression,
                reduced.reduction_status(),
            );
        }
    }

    #[test]
    fn generator_trace_requires_compatible_adjoint_spaces() {
        crate::test_support::test_initialize();
        let simplifier =
            ColorAlgebraSimplifier::new(ColorSimplifySettings::default(), ParseState::default());
        let slots = [3, 8]
            .map(|dimension| ColorAdjoint {}.to_symbolic([Atom::num(dimension), Atom::num(1)]));
        assert!(
            simplifier
                .simplify_generator_trace(
                    &fundamental_rep(Atom::num(3)),
                    &slots,
                    LineContext::Unknown
                )
                .is_none()
        );
        assert!(
            simplifier
                .simplify_generator_trace(
                    &adjoint_rep(Atom::num(8)),
                    &slots[..1],
                    LineContext::Unknown
                )
                .is_none()
        );
    }
}
