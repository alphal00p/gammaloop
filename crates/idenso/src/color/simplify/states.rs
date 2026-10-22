//! Canonical generated states.
//!
//! Identities applied to the terms of a row generate colour monomials that
//! coincide up to the names of their dummy indices: cutting a long adjoint
//! cycle whose ports meet a symmetric invariant, or permuted copies of one
//! cycle, generates mostly such duplicates. FORM merges them at its next sort
//! after renumbering dummies. Here each generated monomial receives canonical
//! dummy labels from one fixed pool, and the planner's sum merges equal states
//! and cancels opposite ones. Every contracted adjoint index of a generated
//! monomial is a dummy, whether the kernel or the input named it; open
//! labels keep their names. A one-shot decomposition of an isolated line
//! generates distinct monomials and is not relabelled; nor is a term whose
//! relabelling would neither rename nor merge anything.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use symbolica::{
    domains::{float::Complex, rational::Rational},
    id::AliasedAtom,
};

use super::*;

/// Monomials one generated term may be distributed into for merging.
pub(super) const STATE_LIMIT: usize = 1 << 16;
/// Monomials a complete decomposition of an isolated line may be
/// distributed into: its states are distinct, so relabelling only buys
/// canonical labels.
pub(super) const DISTINCT_STATE_LIMIT: usize = 64;

impl ColorAlgebraSimplifier {
    /// Distribute the colour sums of a generated whole term and give every
    /// resulting monomial canonical dummy labels; `None` when that would
    /// neither rename nor merge anything. Factors without colour ports, such
    /// as spectator sums, invariant coefficients and powers, are not
    /// distributed. The labels of the factors `beside` the term stay
    /// unused, like those of its common factors.
    pub(super) fn canonical_states(
        &self,
        term: AtomView<'_>,
        limit: usize,
        beside: &[AtomView<'_>],
    ) -> Option<Atom> {
        let (distributed, common) = distributed_factors(term)?;
        let core;
        let distributed = if common.is_empty() {
            term
        } else {
            core = Atom::mul_many(distributed);
            core.as_view()
        };
        // Common factors and those beside the term keep their labels: the
        // pool avoids them.
        let outside = common.iter().chain(beside).copied().collect::<Vec<_>>();
        let blocked = MonomialSlots::labels(&outside);
        let blocked = blocked.iter().map(Atom::as_view).collect::<Vec<_>>();
        let states = colour_states(distributed, limit)?;
        let count = states.len();
        let mut renamed = false;
        let mut canonical = Vec::with_capacity(count);
        for factors in states {
            let state = Atom::mul_many(factors);
            let relabelled = self.canonical_state(&state, &blocked)?;
            renamed |= relabelled != state;
            canonical.push(relabelled);
        }
        let merged = Atom::add_many(canonical);
        let merged_count = match merged.as_view() {
            AtomView::Add(sum) => sum.get_nargs(),
            _ => 1,
        };
        (renamed || merged_count < count)
            .then(|| Atom::mul_many(common.into_iter().chain([merged.as_view()])))
    }

    /// Complete the changed contribution before reconciling untouched arithmetic.
    /// The source's colour-free additive contribution is an exact value, not a
    /// set of coefficient spellings. Subtracting it in the protected arithmetic
    /// accounts for row weights and cancellations before it can be reopened.
    /// Nonnumeric coefficient products, function arguments, arithmetic power
    /// payloads and unreached components remain opaque.
    pub(super) fn merging_row_numbers(
        expression: AtomView<'_>,
        source: AtomView<'_>,
        unreached: &[AtomView<'_>],
    ) -> Option<Atom> {
        if expression == source || !matches!(expression, AtomView::Add(_) | AtomView::Mul(_)) {
            return None;
        }
        // Reuse the original product boundary before inspecting its changed
        // child. These factors are immutable spectators, even when indexed.
        let source_factors = multiplicative_factor_views(source);
        let factors = multiplicative_factor_views(expression);
        let common = source_factors
            .iter()
            .copied()
            .filter(|factor| {
                !matches!(factor, AtomView::Num(_))
                    && !atom_contains_color_node(*factor)
                    && factors.contains(factor)
            })
            .collect::<Vec<_>>();
        if !common.is_empty() {
            let child = Atom::mul_many(
                factors
                    .into_iter()
                    .filter(|factor| !common.contains(factor)),
            );
            let before = Atom::mul_many(
                source_factors
                    .into_iter()
                    .filter(|factor| !common.contains(factor)),
            );
            return Self::merging_row_numbers(child.as_view(), before.as_view(), unreached)
                .map(|merged| Atom::mul_many(common.into_iter().chain([merged.as_view()])));
        }
        let mut exact = true;
        source.visitor(&mut |node| {
            if let AtomView::Num(number) = node {
                exact &= matches!(
                    number.get_coeff_view(),
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                );
            }
            matches!(node, AtomView::Add(_) | AtomView::Mul(_))
        });
        if !exact {
            return None;
        }
        // This projection follows the original Add/Mul tree, without distributing
        // a product or entering a function/power. In x*(C+M), for example, the
        // untouched contribution is x*M; in C*M it is zero. Thus equal M's in
        // different additive contexts never share an ownership flag.
        let untouched = source.replace_map(|node, _, out| {
            if unreached.contains(&node) || !atom_contains_color_node(node) {
                **out = node.to_owned();
            } else if !matches!(node, AtomView::Add(_) | AtomView::Mul(_)) {
                **out = Atom::Zero;
            }
        });
        let mut has_sum = false;
        let mut arithmetic = HashSet::new();
        for value in [expression, untouched.as_view()] {
            value.visitor(&mut |node| match node {
                AtomView::Add(_) => {
                    arithmetic.insert(node.to_owned());
                    has_sum = true;
                    true
                }
                AtomView::Mul(_) => {
                    arithmetic.insert(node.to_owned());
                    true
                }
                AtomView::Pow(power) => {
                    // Arithmetic payloads need a whole-power mask. Powers of
                    // leaves are already opaque to both arithmetic operations;
                    // keeping them literal also lets x*x and x^2 share a key.
                    let (base, exponent) = power.get_base_exp();
                    if [base, exponent].into_iter().any(|part| {
                        matches!(part, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
                    }) {
                        arithmetic.insert(node.to_owned());
                    }
                    false
                }
                AtomView::Num(number) => {
                    exact &= matches!(
                        number.get_coeff_view(),
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    );
                    false
                }
                _ => false,
            });
        }
        if !exact || !has_sum {
            return None;
        }
        let mut retained = BTreeSet::new();
        source.visitor(&mut |node| {
            if matches!(node, AtomView::Add(_) | AtomView::Mul(_)) {
                if arithmetic.contains(node.get_data())
                    && (unreached.contains(&node) || !atom_contains_color_node(node))
                {
                    retained.insert(node.to_owned());
                }
                true
            } else {
                false
            }
        });
        retained.extend(
            arithmetic
                .into_iter()
                .filter(|node| matches!(node.as_view(), AtomView::Pow(_))),
        );
        let mut used = expression.get_all_symbols(true);
        used.extend(source.get_all_symbols(true));
        let mut aliases = HashMap::new();
        let mut linear_aliases = HashSet::new();
        let mut serial = 0;
        for original in retained {
            let alias = loop {
                let symbol = symbolica::symbol!(format!("idenso::color::row_coefficient_{serial}"));
                serial += 1;
                if used.insert(symbol) {
                    break Atom::var(symbol);
                }
            };
            let linear_sum = matches!(original.as_view(), AtomView::Add(_))
                || matches!(original.as_view(), AtomView::Mul(product)
                    if product.get_nargs() == 2 && product.get_coefficient().is_some());
            if linear_sum
                && !unreached.contains(&original.as_view())
                && !atom_contains_color_node(original.as_view())
            {
                linear_aliases.insert(alias.clone());
            }
            aliases.insert(original, alias);
        }
        let protect = |value: AtomView<'_>, generated_root: bool| {
            value.replace_map(|node, _, out| {
                if matches!(node, AtomView::Fun(_)) {
                    **out = node.to_owned();
                } else if (!generated_root || node != value || unreached.contains(&node))
                    && let Some(alias) = aliases.get::<[u8]>(node.get_data())
                {
                    **out = alias.clone();
                }
            })
        };
        // A newly generated root is never an old coefficient occurrence, even
        // if a vanished row had exactly that coefficient.
        let mut protected = AliasedAtom::from(protect(expression, true));
        // Mask original occurrences before evaluating the projection. A newly
        // formed x*M must not acquire the mask of an unrelated source product.
        let mut contexts = BTreeSet::new();
        let untouched = source.replace_map(|node, _, out| {
            if unreached.contains(&node) || !atom_contains_color_node(node) {
                **out = protect(node, false);
            } else if !matches!(node, AtomView::Add(_) | AtomView::Mul(_)) {
                **out = Atom::Zero;
            } else if let AtomView::Mul(product) = node
                && product.iter().any(|factor| {
                    matches!(factor, AtomView::Add(_)) && atom_contains_color_node(factor)
                })
            {
                contexts.extend(
                    product
                        .iter()
                        .filter(|factor| {
                            !matches!(factor, AtomView::Num(_))
                                && !atom_contains_color_node(*factor)
                        })
                        .map(|factor| protect(factor, false)),
                );
            }
        });
        let mut powered_contexts = BTreeSet::new();
        untouched.visitor(&mut |node| {
            if let AtomView::Pow(power) = node {
                if contexts.contains(power.get_base().get_data()) {
                    powered_contexts.insert(node.to_owned());
                }
                false
            } else {
                matches!(node, AtomView::Add(_) | AtomView::Mul(_))
            }
        });
        contexts.extend(powered_contexts);
        // Restore the source's factored sum contexts before entering their
        // children. collect_factors only extracts factors shared by *all*
        // summands; here 4+x*(4+M)-x*M must first group the two x terms.
        // coefficient() extracts literal factors and never expands a product.
        let collect_contexts = |value: &Atom| {
            value.replace_map(|node, _, out| {
                if matches!(node, AtomView::Fun(_) | AtomView::Pow(_)) {
                    **out = node.to_owned();
                    return;
                }
                let AtomView::Add(sum) = node else { return };
                for factor in &contexts {
                    let (inside, outside): (Vec<_>, Vec<_>) = sum.iter().partition(|term| {
                        multiplicative_factor_views(*term).contains(&factor.as_view())
                    });
                    if inside.len() > 1 {
                        let grouped = factor * node.coefficient(factor);
                        **out = Atom::add_many(
                            outside
                                .into_iter()
                                .map(|term| term.to_owned())
                                .chain([grouped]),
                        );
                        return;
                    }
                }
            })
        };
        for (original, alias) in &aliases {
            // Descendant powers get whole aliases before descent is stopped;
            // function arguments never see temporary symbols.
            let body = original.replace_map(|node, _, out| {
                if matches!(node, AtomView::Fun(_)) {
                    **out = node.to_owned();
                } else if node != original.as_view()
                    && let Some(alias) = aliases.get::<[u8]>(node.get_data())
                {
                    **out = alias.clone();
                } else if matches!(node, AtomView::Pow(_)) {
                    **out = node.to_owned();
                }
            });
            protected.register_alias(alias.clone(), body);
        }
        let linear = protected
            .get_aliases()
            .iter()
            .filter(|(alias, _)| linear_aliases.contains(*alias))
            .map(|(alias, body)| (alias.clone(), body.clone()))
            .collect::<HashMap<_, _>>();
        // A comparison form opens only sums and rational scaling. Every other
        // factor is a single key, so this cannot form a Cartesian product.
        let coefficients = |value: AtomView<'_>| -> Option<BTreeMap<Atom, Rational>> {
            let mut terms = BTreeMap::new();
            let mut pending = vec![(
                value.to_owned(),
                Complex::new(Rational::one(), Rational::zero()),
            )];
            while let Some((mut body, mut number)) = pending.pop() {
                if let AtomView::Mul(product) = body.as_view()
                    && let Some(coefficient) = product.get_coefficient()
                {
                    number *= Complex::<Rational>::try_from(coefficient).ok()?;
                    body = body.as_view() / coefficient;
                }
                if let Some(sum) = linear.get(&body) {
                    pending.push((sum.clone(), number));
                } else if let AtomView::Add(sum) = body.as_view() {
                    pending.extend(sum.iter().map(|term| (term.to_owned(), number.clone())));
                } else {
                    if matches!(body.as_view(), AtomView::Num(_)) {
                        number *= Complex::<Rational>::try_from(body.as_view()).ok()?;
                        body = Atom::one();
                    }
                    // Real and imaginary parts are independent exact rational
                    // coordinates; the imaginary unit remains part of the key.
                    if !number.re.is_zero() {
                        *terms.entry(body.clone()).or_insert_with(Rational::zero) += number.re;
                    }
                    if !number.im.is_zero() {
                        *terms.entry(Atom::i() * body).or_insert_with(Rational::zero) += number.im;
                    }
                }
            }
            terms.retain(|_, number| !number.is_zero());
            Some(terms)
        };
        let colour_free = !atom_contains_color_node(expression);
        let mut current = protected.get_root() - &untouched;
        loop {
            // Cancel the untouched projection while its scopes are still
            // masked. In particular x*(4+M)-x*M becomes 4*x before M can open.
            let factored = if colour_free {
                collect_contexts(&current).collect_factors()
            } else {
                current.clone()
            };
            let next = factored
                .replace_map(|node, _, out| {
                    if matches!(node, AtomView::Fun(_) | AtomView::Pow(_)) {
                        **out = node.to_owned();
                        return;
                    }
                    if let AtomView::Mul(product) = node
                        && product.iter().any(|factor| {
                            coefficients(factor).is_some_and(|terms| terms.is_empty())
                        })
                    {
                        **out = Atom::Zero;
                        return;
                    }
                    let AtomView::Add(sum) = node else { return };
                    let mut terms = BTreeMap::<Atom, Rational>::new();
                    let mut cancelled = false;
                    for term in sum {
                        let Some(contribution) = coefficients(term) else {
                            return;
                        };
                        cancelled |= contribution.is_empty();
                        for (body, number) in contribution {
                            let previous = terms.entry(body).or_insert_with(Rational::zero);
                            cancelled |= !previous.is_zero()
                                && previous.is_negative() != number.is_negative();
                            *previous += number;
                        }
                    }
                    // A completed child remains factored unless this very sum has
                    // an opposite contribution. A cancellation elsewhere cannot
                    // open the M in x*(4+M), nor an internally reducible spectator.
                    if cancelled {
                        **out = Atom::add_many(
                            terms
                                .into_iter()
                                .map(|(body, number)| Atom::num(number) * body),
                        );
                    }
                })
                .replace_map(|node, _, out| {
                    // Numeric distribution completes coefficient arithmetic,
                    // not the kernel's factored colour decomposition. In an
                    // isolated trace, pushing a Bernoulli coefficient into a
                    // nested colour sum changes its certified arithmetic shape
                    // without merging any states. Actual cancellations above
                    // still change that shape and revoke the certificate.
                    if matches!(node, AtomView::Fun(_) | AtomView::Pow(_)) {
                        **out = node.to_owned();
                    } else if colour_free || !atom_contains_color_port_or_work(node) {
                        **out = node.expand_num();
                    }
                });
            if next == current {
                break;
            }
            current = next;
        }
        // Only surviving changed coefficients may consume an untouched scope.
        // Compare each occurrence separately; an untouched internally zero
        // expression is retained when nothing from the changed subtotal uses it.
        let mut active = coefficients(current.as_view())?;
        if active.is_empty() {
            current = Atom::Zero;
        }
        let occurrences = match untouched.as_view() {
            AtomView::Add(sum) => sum.iter().map(|term| term.to_owned()).collect::<Vec<_>>(),
            _ => vec![untouched],
        };
        let mut retained = Vec::new();
        let mut consumed = false;
        for occurrence in occurrences {
            let original = coefficients(occurrence.as_view())?;
            let mut remainder = original.clone();
            let mut touched = false;
            for (body, number) in &mut remainder {
                let Some(changed) = active.get_mut(body) else {
                    continue;
                };
                if changed.is_zero() || changed.is_negative() == number.is_negative() {
                    continue;
                }
                let amount = changed.abs().min(number.abs());
                let amount = if number.is_negative() {
                    -amount
                } else {
                    amount
                };
                *changed += &amount;
                *number -= amount;
                touched = true;
            }
            if !touched {
                retained.push(occurrence);
                continue;
            }
            consumed = true;
            // A partial cancellation of a rational multiple still preserves
            // the original scope when its remainder has one uniform scale.
            let (body, number) = original.first_key_value().unwrap();
            let scale = &remainder[body] / number;
            if original
                .iter()
                .all(|(body, number)| remainder[body] == &scale * number)
            {
                retained.push(Atom::num(scale) * occurrence);
            } else {
                retained.push(Atom::add_many(
                    remainder
                        .into_iter()
                        .map(|(body, number)| Atom::num(number) * body),
                ));
            }
        }
        if consumed {
            current = Atom::add_many(
                active
                    .into_iter()
                    .map(|(body, number)| Atom::num(number) * body),
            );
        }
        retained.push(current);
        let mut result = Atom::add_many(retained);
        if colour_free {
            loop {
                let next = collect_contexts(&result).collect_factors().expand_num();
                if next == result {
                    break;
                }
                result = next;
            }
        }
        let merged = protected.map_root(|_| result).into_inner();
        (merged.as_view() != expression).then_some(merged)
    }
}

/// The colour monomials of `term` as lists of factors, or `None` beyond
/// `limit` or with inexact coefficients, whose addition order is kept.
fn colour_states(term: AtomView<'_>, limit: usize) -> Option<Vec<Vec<AtomView<'_>>>> {
    match term {
        AtomView::Num(number)
            if !matches!(
                number.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            ) =>
        {
            None
        }
        AtomView::Add(sum) if atom_contains_color_port_or_work(term) => {
            let mut states = Vec::new();
            for summand in sum.iter() {
                states.extend(colour_states(summand, limit)?);
                if states.len() > limit {
                    return None;
                }
            }
            Some(states)
        }
        AtomView::Mul(product) => {
            let mut states = vec![Vec::new()];
            for factor in product.iter() {
                let factor_states = colour_states(factor, limit)?;
                if states.len() * factor_states.len() > limit {
                    return None;
                }
                states = states
                    .iter()
                    .flat_map(|state| {
                        factor_states.iter().map(move |factors| {
                            state.iter().chain(factors).copied().collect::<Vec<_>>()
                        })
                    })
                    .collect();
            }
            Some(states)
        }
        _ => Some(vec![vec![term]]),
    }
}

/// Contract the adjoint metric factors of a monomial against their one
/// partner, as the planner would after the relabelling; `None` without one.
fn eliminating_metrics(monomial: &Atom) -> Option<Atom> {
    let (metrics, rest): (Vec<_>, Vec<_>) = multiplicative_factor_views(monomial.as_view())
        .into_iter()
        .partition(|factor| emitted_adjoint_metric(*factor).is_some());
    if metrics.is_empty() {
        return None;
    }
    ProductView::eliminating_emitted_metrics(&Atom::mul_many(metrics), rest.into_iter())
}

/// One walk over a monomial: its adjoint slots with their numbers of
/// occurrences, the counter labels outside them, and whether it has an
/// ordered generator word, a power of a tensor or a slot used more than
/// twice. A power contracts its copies with each other: it is a closed
/// scope whose labels are not counted with the monomial's. An over-used
/// slot is a summand-local dummy that the distribution of its sum
/// captured, because another factor spells it too.
#[derive(Default)]
struct MonomialSlots<'a> {
    adjoint: Vec<(AtomView<'a>, usize)>,
    blocked: Vec<AtomView<'a>>,
    ordered: bool,
    closed_power: bool,
    over_used: bool,
}

impl<'a> MonomialSlots<'a> {
    /// The adjoint labels and counter labels of some factors, also inside
    /// their powers, which a relabelling beside them must not reuse.
    fn labels(factors: &[AtomView<'_>]) -> Vec<Atom> {
        let mut matcher = SlotMatcher::default();
        let mut labels = Vec::new();
        for factor in factors {
            factor.visitor(&mut |node| {
                if let Some((_, index)) = representation_slot_view(node, CS.adjoint_rep) {
                    labels.push(index.to_owned());
                    return false;
                }
                match matcher.classify(node) {
                    SlotMatch::Explicit(slot) => {
                        if is_counter_dummy(slot.index()) {
                            labels.push(slot.index().to_owned());
                        }
                        false
                    }
                    SlotMatch::Opaque => false,
                    SlotMatch::Other => {
                        if is_counter_dummy(node) {
                            labels.push(node.to_owned());
                        }
                        true
                    }
                }
            });
        }
        labels
    }

    fn of(monomial: AtomView<'a>) -> Self {
        let mut slots = Self::default();
        let mut other = Vec::new();
        slots.visit(monomial, &mut SlotMatcher::default(), &mut other);
        slots.over_used = slots.adjoint.iter().any(|(_, count)| *count > 2)
            || other.iter().any(|(_, count)| *count > 2);
        slots
    }

    fn visit(
        &mut self,
        node: AtomView<'a>,
        matcher: &mut SlotMatcher,
        other: &mut Vec<(SlotKey<'a>, usize)>,
    ) {
        if representation_slot_view(node, CS.adjoint_rep).is_some() {
            match self.adjoint.iter_mut().find(|(seen, _)| *seen == node) {
                Some((_, count)) => *count += 1,
                None => self.adjoint.push((node, 1)),
            }
            return;
        }
        match node {
            AtomView::Fun(function) => {
                if let Some((_, factors)) = shadowing::trace_parts(function) {
                    self.ordered |= factors.len() > 1;
                }
                if let SlotMatch::Explicit(slot) = matcher.classify(node) {
                    let key = (slot.representation().head(), slot.dimension(), slot.index());
                    match other.iter_mut().find(|(seen, _)| *seen == key) {
                        Some((_, count)) => *count += 1,
                        None => other.push((key, 1)),
                    }
                    if is_counter_dummy(slot.index()) {
                        self.blocked.push(slot.index());
                    }
                    return;
                }
                for arg in function.iter() {
                    self.visit(arg, matcher, other);
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    self.visit(factor, matcher, other);
                }
            }
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    self.visit(term, matcher, other);
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let mut scope = Self::default();
                let mut scope_other = Vec::new();
                scope.visit(base, matcher, &mut scope_other);
                scope.visit(exponent, matcher, &mut scope_other);
                self.closed_power |= !scope.adjoint.is_empty() || !scope_other.is_empty();
                self.ordered |= scope.ordered;
                self.blocked.extend(scope.blocked);
            }
            _ => {
                if is_counter_dummy(node) {
                    self.blocked.push(node);
                }
            }
        }
    }
}

/// The factors of a term that the distribution of its colour sums
/// multiplies out, and the common colour factors outside their connected
/// colour component, such as an independent closed trace, which stay
/// factored. (Foreign colour tensors that no identity reaches never enter a
/// row.) `None` when the term has colour sums in more than one connected
/// component: distributing those would multiply out independent sums.
fn distributed_factors(term: AtomView<'_>) -> Option<(Vec<AtomView<'_>>, Vec<AtomView<'_>>)> {
    let AtomView::Mul(product) = term else {
        return Some((vec![term], Vec::new()));
    };
    let mut matcher = SlotMatcher::default();
    let all = product.iter().collect::<Vec<_>>();
    let coloured = all
        .iter()
        .filter(|factor| atom_contains_color_port_or_work(**factor))
        .copied()
        .collect::<Vec<_>>();
    let factors = coloured
        .iter()
        .map(|&factor| {
            let mut slots = Vec::new();
            factor.visitor(&mut |node| match matcher.classify(node) {
                SlotMatch::Explicit(_) => {
                    if let Ok(slot) = matcher.parse::<LibraryRep, AbstractIndex>(node) {
                        slots.push(if slot.rep_name().is_dual() {
                            slot.dual()
                        } else {
                            slot
                        });
                    }
                    false
                }
                SlotMatch::Opaque => false,
                SlotMatch::Other => true,
            });
            (matches!(factor, AtomView::Add(_)), slots)
        })
        .collect::<Vec<_>>();
    if factors.iter().filter(|(sum, _)| *sum).count() == 0 {
        return Some((all, Vec::new()));
    }
    // Union the factors that share a slot.
    let mut component = (0..factors.len()).collect::<Vec<_>>();
    fn root(component: &mut [usize], mut index: usize) -> usize {
        while component[index] != index {
            component[index] = component[component[index]];
            index = component[index];
        }
        index
    }
    for left in 0..factors.len() {
        for right in left + 1..factors.len() {
            if factors[left]
                .1
                .iter()
                .any(|slot| factors[right].1.contains(slot))
            {
                let (a, b) = (root(&mut component, left), root(&mut component, right));
                component[a] = b;
            }
        }
    }
    let mut sum_components = factors
        .iter()
        .enumerate()
        .filter(|(_, (sum, _))| *sum)
        .map(|(index, _)| root(&mut component, index))
        .collect::<Vec<_>>();
    sum_components.sort_unstable();
    sum_components.dedup();
    let [distributed] = sum_components[..] else {
        return None;
    };
    let common = (0..coloured.len())
        .filter(|&index| root(&mut component, index) != distributed)
        .map(|index| coloured[index])
        .collect::<Vec<_>>();
    let distributed = all
        .into_iter()
        .filter(|factor| !common.contains(factor))
        .collect();
    Some((distributed, common))
}

/// A label allocated by a dummy counter, `d_n`, which the pool must avoid
/// where it is not renamed.
fn is_counter_dummy(index: AtomView<'_>) -> bool {
    matches!(index, AtomView::Var(variable)
        if variable.get_symbol().has_tag(spenso::structure::abstract_index::DUMMY_INDEX_TAG))
}

impl ColorAlgebraSimplifier {
    /// Rename the contracted adjoint indices of one monomial to the pool d_0,
    /// d_1, ... and canonize them; open labels keep their names. Adjoint
    /// metrics are contracted first. A monomial with an ordered generator word
    /// is still being decomposed and keeps its labels: the kernel orients such
    /// words by comparing atoms.
    ///
    /// The labelling by first occurrence is cheap and often agrees between
    /// equal states: it keys a cache of canonical forms for the call, so the
    /// graph canonization runs about once per distinct labelling.
    ///
    /// `None` when a slot of the monomial is used more than twice.
    fn canonical_state(&self, monomial: &Atom, blocked: &[AtomView<'_>]) -> Option<Atom> {
        let slots = MonomialSlots::of(monomial.as_view());
        if slots.over_used {
            return None;
        }
        if slots.ordered || slots.closed_power {
            return Some(monomial.clone());
        }
        let eliminated = eliminating_metrics(monomial);
        let (monomial, mut slots) = match &eliminated {
            Some(eliminated) => (eliminated, MonomialSlots::of(eliminated.as_view())),
            None => (monomial, slots),
        };
        slots.blocked.extend_from_slice(blocked);
        let Some((relabelled, pool)) = labelled_by_first_occurrence(monomial, slots) else {
            return Some(monomial.clone());
        };
        // One dummy has only one labelling.
        if pool.len() == 1 {
            return Some(relabelled);
        }
        let hit = self.canonical_cache.borrow().get(&relabelled).cloned();
        if let Some(canonical) = hit {
            return Some(canonical);
        }
        let canonical = match relabelled.canonize_tensors(pool) {
            Ok(canonical) => canonical.canonical_form,
            Err(_) => relabelled.clone(),
        };
        self.canonical_cache
            .borrow_mut()
            .insert(relabelled, canonical.clone());
        Some(canonical)
    }
}

/// Rename the contracted adjoint slots to d_0, d_1, ... in order of first
/// occurrence, skipping dummy labels other slots use. Also returns
/// the renamed slots with their groups: dummies of different adjoint spaces
/// never exchange labels. `None` without such a dummy.
fn labelled_by_first_occurrence(
    monomial: &Atom,
    slots: MonomialSlots<'_>,
) -> Option<(Atom, Vec<(Atom, Atom)>)> {
    let MonomialSlots {
        adjoint,
        mut blocked,
        ..
    } = slots;
    let (dummies, kept): (Vec<_>, Vec<_>) = adjoint.into_iter().partition(|(_, count)| *count == 2);
    if dummies.is_empty() {
        return None;
    }
    blocked.extend(
        kept.iter()
            .filter_map(|(slot, _)| representation_slot_view(*slot, CS.adjoint_rep))
            .map(|(_, index)| index),
    );
    let mut labels = (0..)
        .map(|n| AbstractIndex::Dummy(n).to_atom())
        .filter(|label| blocked.iter().all(|blocked| *blocked != label.as_view()));
    let renamed = dummies
        .iter()
        .map(|(slot, _)| {
            let (dimension, _) = representation_slot_view(*slot, CS.adjoint_rep).unwrap();
            let target = FunctionBuilder::new(CS.adjoint_rep)
                .add_arg(dimension)
                .add_arg(labels.next().unwrap())
                .finish();
            let group = FunctionBuilder::new(CS.adjoint_rep)
                .add_arg(dimension)
                .finish();
            (*slot, target, group)
        })
        .collect::<Vec<_>>();
    let relabelled = monomial.replace_map(|node, _context, out| {
        if let Some((_, target, _)) = renamed.iter().find(|(source, ..)| node == *source) {
            **out = target.clone();
        }
    });
    let pool = renamed
        .into_iter()
        .map(|(_, target, group)| (target, group))
        .collect();
    Some((relabelled, pool))
}

/// The states of a term, as the relabelling distributes them.
#[cfg(test)]
pub(super) fn colour_states_for_test(term: AtomView<'_>) -> Option<usize> {
    let states = colour_states(term, STATE_LIMIT)?;
    let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
    states
        .into_iter()
        .map(|factors| simplifier.canonical_state(&Atom::mul_many(factors), &[]))
        .collect::<Option<Vec<_>>>()
        .map(|states| states.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn adjoint(label: &str) -> Atom {
        ColorAdjoint {}.to_symbolic([
            Atom::num(8),
            Atom::var(symbolica::symbol!(&format!("idenso::states_test::{label}"))),
        ])
    }

    fn kernel_dummy(index: usize) -> Atom {
        ColorAdjoint {}.to_symbolic([
            Atom::num(8),
            AbstractIndex::Dummy(900_000 + index).to_atom(),
        ])
    }

    /// A colour factor outside the connected component of the distributed
    /// sum, such as a foreign tensor with a free adjoint index, stays a
    /// common factor; the pool avoids its labels.
    #[test]
    fn common_colour_factors_stay_outside_the_states() {
        crate::test_support::test_initialize();
        let [r1, r2, r3, r4] = ["r1", "r2", "r3", "r4"].map(adjoint);
        let [c0, c1] = [0, 1].map(kernel_dummy);
        let vector =
            |slot: &Atom| function!(spenso::tensor_symbol!("idenso::states_test::V"), slot);
        let blocked =
            ColorAdjoint {}.to_symbolic([Atom::num(8), AbstractIndex::Dummy(0).to_atom()]);
        let foreign = vector(&adjoint("z")) * vector(&blocked) * vector(&blocked);
        let sum = color_f!(&r1, &r2, &c0) * color_f!(&c0, &r3, &r4)
            + color_f!(&r1, &r3, &c1) * color_f!(&c1, &r2, &r4);
        let term = &foreign * &sum;
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let states = simplifier
            .canonical_states(term.as_view(), STATE_LIMIT, &[])
            .expect("the sum's dummies are renamed");
        let AtomView::Mul(product) = states.as_view() else {
            panic!("{states}")
        };
        assert!(
            product.iter().any(|factor| factor == foreign.as_view()
                || multiplicative_factor_views(foreign.as_view()).contains(&factor)),
            "{states}"
        );
        let AtomView::Add(merged) = product
            .iter()
            .find(|factor| matches!(factor, AtomView::Add(_)))
            .unwrap()
        else {
            unreachable!()
        };
        for state in merged.iter() {
            assert!(!state.contains(blocked.as_view()), "{state}");
            assert!(!state.contains(adjoint("z").as_view()), "{state}");
        }
    }

    /// Two labellings of the same structure-constant tree take one form.
    #[test]
    fn relabelled_trees_share_one_canonical_state() {
        crate::test_support::test_initialize();
        let [r1, r2, r3, r4, r5, r6] = ["r1", "r2", "r3", "r4", "r5", "r6"].map(adjoint);
        let [c0, c1, c2] = [0, 1, 2].map(kernel_dummy);
        let first = color_f!(&r1, &r2, &c0)
            * color_f!(&c0, &r3, &c2)
            * color_f!(&r4, &r5, &c1)
            * color_f!(&c2, &c1, &r6);
        let second = color_f!(&r1, &r2, &c2)
            * color_f!(&c0, &r4, &r5)
            * color_f!(&c0, &c1, &r6)
            * color_f!(&r3, &c2, &c1);
        let simplifier = ColorAlgebraSimplifier::new(Default::default(), ParseState::default());
        let (first, second) = (
            simplifier.canonical_state(&first, &[]).unwrap(),
            simplifier.canonical_state(&second, &[]).unwrap(),
        );
        assert_eq!(first, second, "{first}\n{second}");
    }
}
