//! Syntactic index occurrences, before structure inference removes contractions.

use crate::structure::slot::{SlotMatch, SlotMatcher};
use symbolica::atom::AtomView;

/// Store common symbol IDs separately from arbitrary borrowed index payloads.
struct Indices<'a> {
    symbols: Vec<u32>,
    others: Vec<AtomView<'a>>,
}

#[derive(Default)]
pub(super) struct RepeatedIndices {
    slots: SlotMatcher,
}

impl RepeatedIndices {
    // Decide once at the first hit whether to finish observation. Continuing
    // uses the same traversal, with index bookkeeping disabled for the suffix.
    pub(super) fn contains<'a>(
        expr: AtomView<'a>,
        mut observe: impl FnMut(AtomView<'a>, &SlotMatch<'a>),
        on_first_repeat: impl FnOnce() -> bool,
    ) -> bool {
        let mut decision = None;
        let mut on_first_repeat = Some(on_first_repeat);
        let mut observe_remainder =
            || *decision.get_or_insert_with(|| on_first_repeat.take().unwrap()());
        let mut scanner = Self::default();
        let mut seen = Indices {
            symbols: Vec::with_capacity(16),
            others: Vec::new(),
        };
        // Top-level alternatives do not need their possible-index sets merged.
        let repeated = match expr {
            AtomView::Add(sum) => {
                observe(expr, &SlotMatch::Other);
                let mut terms = sum.iter();
                let repeated = terms.any(|term| {
                    seen.symbols.clear();
                    seen.others.clear();
                    scanner.visit::<true>(term, &mut seen, &mut observe, &mut observe_remainder)
                });
                if repeated && observe_remainder() {
                    scanner.visit_children::<false>(
                        terms,
                        &mut seen,
                        &mut observe,
                        &mut observe_remainder,
                    );
                }
                repeated
            }
            _ => scanner.visit::<true>(expr, &mut seen, &mut observe, &mut observe_remainder),
        };
        // A root power can discover its first repetition in an explicit-slot
        // exponent, with no ancestor left to request a suffix. Resolve that
        // terminal decision too; already handled hits reuse the memoized value.
        if repeated {
            observe_remainder();
        }
        repeated
    }

    fn visit<'a, const TRACK_INDICES: bool>(
        &mut self,
        expr: AtomView<'a>,
        seen: &mut Indices<'a>,
        observe: &mut impl FnMut(AtomView<'a>, &SlotMatch<'a>),
        observe_remainder: &mut impl FnMut() -> bool,
    ) -> bool {
        let classification = match expr {
            AtomView::Fun(_) => self.slots.classify(expr),
            _ => SlotMatch::Other,
        };
        observe(expr, &classification);
        match expr {
            AtomView::Num(_) | AtomView::Var(_) => false,
            AtomView::Fun(fun) => {
                match classification {
                    SlotMatch::Explicit(_) if !TRACK_INDICES => false,
                    SlotMatch::Explicit(slot) => match slot.index() {
                        AtomView::Var(var) => {
                            let index = var.get_symbol_id();
                            if seen.symbols.contains(&index) {
                                return true;
                            }
                            seen.symbols.push(index);
                            false
                        }
                        index => {
                            if seen.others.contains(&index) {
                                return true;
                            }
                            seen.others.push(index);
                            false
                        }
                    },
                    // Dimensions and index payloads are opaque, including malformed slots.
                    SlotMatch::Opaque => false,
                    SlotMatch::Other => self.visit_children::<TRACK_INDICES>(
                        fun.iter(),
                        seen,
                        observe,
                        observe_remainder,
                    ),
                }
            }
            AtomView::Mul(product) => self.visit_children::<TRACK_INDICES>(
                product.iter(),
                seen,
                observe,
                observe_remainder,
            ),
            AtomView::Add(sum) if !TRACK_INDICES => {
                self.visit_children::<false>(sum.iter(), seen, observe, observe_remainder)
            }
            AtomView::Add(sum) => {
                let outer_symbols = seen.symbols.len();
                let outer_others = seen.others.len();
                let mut possible_symbols = Vec::new();
                let mut possible_others = Vec::new();
                let mut terms = sum.iter();
                while let Some(term) = terms.next() {
                    if self.visit::<TRACK_INDICES>(term, seen, observe, observe_remainder) {
                        if observe_remainder() {
                            self.visit_children::<false>(terms, seen, observe, observe_remainder);
                        }
                        return true;
                    }
                    for &index in &seen.symbols[outer_symbols..] {
                        if !possible_symbols.contains(&index) {
                            possible_symbols.push(index);
                        }
                    }
                    for &index in &seen.others[outer_others..] {
                        if !possible_others.contains(&index) {
                            possible_others.push(index);
                        }
                    }
                    seen.symbols.truncate(outer_symbols);
                    seen.others.truncate(outer_others);
                }
                // Later product factors can meet an index from any summand.
                seen.symbols.extend(possible_symbols);
                seen.others.extend(possible_others);
                false
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                if !TRACK_INDICES {
                    self.visit::<false>(base, seen, observe, observe_remainder);
                    self.visit::<false>(exponent, seen, observe, observe_remainder);
                    return false;
                }
                let before = seen.symbols.len() + seen.others.len();
                if self.visit::<TRACK_INDICES>(base, seen, observe, observe_remainder) {
                    if observe_remainder() {
                        self.visit::<false>(exponent, seen, observe, observe_remainder);
                    }
                    return true;
                }
                if seen.symbols.len() + seen.others.len() > before
                    && match i64::try_from(exponent) {
                        Ok(n) => n.unsigned_abs() > 1,
                        Err(_) => matches!(exponent, AtomView::Num(number)
                            if number.get_coeff_view().is_integer()),
                    }
                {
                    if observe_remainder() {
                        self.visit::<false>(exponent, seen, observe, observe_remainder);
                    }
                    return true;
                }
                self.visit::<TRACK_INDICES>(exponent, seen, observe, observe_remainder)
            }
        }
    }

    fn visit_children<'a, const TRACK_INDICES: bool>(
        &mut self,
        mut children: impl Iterator<Item = AtomView<'a>>,
        seen: &mut Indices<'a>,
        observe: &mut impl FnMut(AtomView<'a>, &SlotMatch<'a>),
        observe_remainder: &mut impl FnMut() -> bool,
    ) -> bool {
        let repeated = children
            .any(|child| self.visit::<TRACK_INDICES>(child, seen, observe, observe_remainder));
        if repeated && observe_remainder() {
            self.visit_children::<false>(children, seen, observe, observe_remainder);
        }
        repeated
    }
}

#[cfg(test)]
mod tests {
    use crate::network::parsing::AtomStructureExt;
    use crate::structure::slot::{SlotMatch, SlotMatcher};
    use symbolica::{
        atom::{Atom, AtomCore},
        parser::ParseSettings,
    };
    fn parse(source: &str) -> Atom {
        Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap()
    }
    #[test]
    fn repeated_explicit_indices_follow_product_and_sum_scopes() {
        crate::structure::representation::initialize();
        let cases = [
            ("g(mink(4,a),mink(4,b))", false),
            ("x*g(mink(4,a),mink(4,b))", false),
            ("g(mink(4,a),mink(4,b))*g(mink(4,c),mink(4,d))", false),
            ("g(mink(4,a),mink(4,b))+g(mink(4,a),mink(4,c))", false),
            (
                "g(mink(4,a),mink(4,b))*epsilon(mink(4,c),mink(4,d),mink(4,e),mink(4,f))",
                false,
            ),
            // Metric traces normalize during construction, before the scan.
            ("g(mink(4,a),mink(4,a))", false),
            ("Unknown(mink(4,a),mink(4,a))", true),
            ("g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))", true),
            (
                "g(mink(4,a),mink(4,b))*epsilon(mink(4,b),mink(4,c),mink(4,d),mink(4,e))",
                true,
            ),
            ("g(mink(4,a),mink(4,b))^2", true),
            ("f(g(mink(4,a),mink(4,a)))", false),
            ("f(Unknown(mink(4,a),mink(4,a)))", true),
            (
                "g(mink(4,a),mink(4,b))*(g(mink(4,b),mink(4,c))+g(mink(4,b),mink(4,d)))",
                true,
            ),
            ("g(mink(4,a),P(mink(4)))", false),
            ("g(P(mink(4)),Q(mink(4)))", false),
            (
                "epsilon(mink(4,a),mink(4,b))*epsilon(mink(4,c),mink(4,d))",
                false,
            ),
            ("g(mink(4,a),mink(4,b))*g(mink(5,a),mink(5,b))", true),
            ("Unknown(mink(4,a))*Another(mink(4,b))", false),
            ("Unknown(mink(4,a))*Another(mink(4,a))", true),
            ("f(Unknown(mink(4,a)))", false),
            ("Unknown(mink(4,a))*(A(mink(4,a))+B(mink(4,b)))", true),
            ("Unknown(mink(4,c))*(A(mink(4,a))+B(mink(4,b)))", false),
            (
                "(A(mink(4,a))+B(mink(4,b)))*(C(mink(4,a))+D(mink(4,c)))",
                true,
            ),
            (
                "(A(mink(4,a))+B(mink(4,b)))*(C(mink(4,c))+D(mink(4,d)))",
                false,
            ),
            ("f(A(mink(4,a))+B(mink(4,a)))", false),
            ("T(mink(4,a))+U(mink(4,b),mink(4,b))", true),
            (
                "T(mink(4,a))*(U(mink(4,b))+V(mink(4,c)))*W(mink(4,c))",
                true,
            ),
            (
                "T(mink(4,a))*(U(mink(4,b))+V(mink(4,c)))*W(mink(4,d))",
                false,
            ),
            ("f(mink(4,1))*g(mink(4,1),mink(4,2))", true),
            ("f(1)*h(1)*P(mink(4,a))", false),
            ("T(mink(4,a))^(-2)", true),
            ("T(mink(4,a))^n", false),
            ("T(mink(4,a))^mink(4,a)", true),
            ("scope(T(mink(4,a))^mink(4,a),tail_marker)", true),
            ("(A(mink(4,a))+B(mink(4,b)))^2", true),
            ("T(mink(4,f(a)))*U(mink(4,f(a)))", true),
            ("T(mink(4,f(a)))*U(mink(4,a))", false),
            ("T(mink(f(mink(4,a))))*U(mink(4,a))", false),
            ("T(mink(4,f(mink(4,a))))*U(mink(4,a))", false),
            ("T(mink(4,a,b))*U(mink(4,a))", false),
            ("uind(mink(4,a))*U(mink(4,a))", true),
            ("dind(mink(4,a))*U(mink(4,a))", true),
            ("dind(mink(4,a),extra)*U(mink(4,a))", false),
            ("dind(dind(mink(4,a),extra))*U(mink(4,a))", false),
            ("T(mink(4,a))^18446744073709551616", true),
            ("T(mink(4,a))^(-18446744073709551616)", true),
        ];
        for (source, expected) in cases {
            let expr = parse(source);
            assert_eq!(expr.has_repeated_explicit_indices(), expected, "{source}");
            let mut visited = Vec::new();
            assert_eq!(
                expr.has_repeated_explicit_indices_with_observer(
                    |node, _| visited.push(node),
                    || true
                ),
                expected,
                "observer: {source}"
            );
            assert_eq!(visited.first().copied(), Some(expr.as_view()), "{source}");
            let mut prefix = Vec::new();
            assert_eq!(
                expr.has_repeated_explicit_indices_with_observer(
                    |node, _| prefix.push(node),
                    || false
                ),
                expected,
                "early observer: {source}"
            );
            assert!(visited.starts_with(&prefix), "{source}");
            if !expected {
                assert_eq!(prefix, visited, "no-hit observation must finish: {source}");
            }
            let mut resumed = Vec::new();
            let mut decisions = 0;
            let repeated = expr.has_repeated_explicit_indices_with_observer(
                |node, _| resumed.push(node),
                || {
                    decisions += 1;
                    true
                },
            );
            assert_eq!(repeated, expected, "resumed observer: {source}");
            assert_eq!(decisions, usize::from(expected), "decision count: {source}");
            assert_eq!(resumed, visited, "resumed visit order: {source}");
        }
    }

    #[test]
    fn first_hit_decision_runs_between_the_observed_prefix_and_suffix() {
        use std::cell::RefCell;
        crate::structure::representation::initialize();
        let terminal_hit = parse("T(mink(4,a))^mink(4,a)");
        let symbolica::atom::AtomView::Pow(power) = terminal_hit.as_view() else {
            panic!("the root exponent fixture must retain its power");
        };
        assert!(matches!(
            SlotMatcher::default().classify(power.get_base_exp().1),
            SlotMatch::Explicit(_)
        ));
        for source in [
            "scope(mink(4,a),mink(4,a),tail_marker)",
            "scope(T(mink(4,a))*(U(mink(4,a))+V(mink(4,b))),tail_marker)",
            "T(mink(4,a),mink(4,a))+U(mink(4,b),mink(4,b))",
            "scope(T(mink(4,a))^2,tail_marker)",
            "T(mink(4,a),mink(4,a))^tail_marker",
            "T(mink(4,a))^mink(4,a)",
            "scope(T(mink(4,a))^mink(4,a),tail_marker)",
            "x^T(mink(4,a),mink(4,a))",
        ] {
            let expression = parse(source);
            let mut prefix = Vec::new();
            let mut full = Vec::new();
            assert!(expression.has_repeated_explicit_indices_with_observer(
                |node, _| prefix.push(node),
                || false
            ));
            assert!(
                expression.has_repeated_explicit_indices_with_observer(
                    |node, _| full.push(node),
                    || true
                )
            );
            for resume in [false, true] {
                let events = RefCell::new(Vec::new());
                // Capturing an owned token makes this a genuine FnOnce; a
                // propagated hit must reuse the decision instead of calling it.
                let token = source.to_owned();
                assert!(expression.has_repeated_explicit_indices_with_observer(
                    |node, _| events.borrow_mut().push(Some(node)),
                    || {
                        drop(token);
                        events.borrow_mut().push(None);
                        resume
                    },
                ));
                let events = events.into_inner();
                assert_eq!(events.iter().filter(|event| event.is_none()).count(), 1);
                assert_eq!(
                    events.iter().position(|event| event.is_none()),
                    Some(prefix.len()),
                    "{source}"
                );
                let observed: Vec<_> = events.into_iter().flatten().collect();
                assert_eq!(&observed, if resume { &full } else { &prefix }, "{source}");
            }
        }
    }

    #[test]
    fn index_observer_reports_slots_without_entering_their_payloads() {
        crate::structure::representation::initialize();
        for source in [
            "T(mink(dim(hidden_dimension),index(hidden_index)))",
            "T(mink(4,index(hidden_index)))+U(mink(4,index(hidden_index)))",
            "T(dind(cof(3,index(hidden_index))))",
            "T(mink(dim(hidden_dimension)))",
            "T(mink(4,a,hidden_payload(mink(4,b))))",
            "T(dind(dind(mink(4,a),extra)))",
        ] {
            let expr = parse(source);
            let mut expected_nodes = vec![expr.as_view()];
            let terms: Vec<_> = match expr.as_view() {
                symbolica::atom::AtomView::Add(sum) => sum.iter().collect(),
                _ => vec![expr.as_view()],
            };
            let mut slots = Vec::new();
            for term in terms {
                if term != expr.as_view() {
                    expected_nodes.push(term);
                }
                let symbolica::atom::AtomView::Fun(tensor) = term else {
                    panic!("expected tensor in {source}");
                };
                expected_nodes.extend(tensor.iter());
            }
            let mut visited = Vec::new();
            assert!(!expr.has_repeated_explicit_indices_with_observer(
                |node, classification| {
                    visited.push(node);
                    if !matches!(classification, SlotMatch::Other) {
                        slots.push(node);
                    }
                },
                || true
            ));
            assert_eq!(visited, expected_nodes, "{source}");
            assert_eq!(slots.len(), expr.nterms(), "{source}");
        }
    }

    #[test]
    fn index_observer_finishes_after_the_boolean_predicate_can_stop() {
        crate::structure::representation::initialize();
        // Function arguments retain their order, so the marker lies after the
        // repeated pair regardless of commutative term or factor ordering.
        let expr = parse("index_observer_scope(mink(4,a),mink(4,a),tail_marker)");
        let mut prefix = Vec::new();
        assert!(
            expr.has_repeated_explicit_indices_with_observer(|node, _| prefix.push(node), || false)
        );
        assert_eq!(prefix.len(), 3);
        let mut visited = Vec::new();
        let mut indices = Vec::new();
        assert!(expr.has_repeated_explicit_indices_with_observer(
            |node, classification| {
                visited.push(node);
                if let SlotMatch::Explicit(slot) = classification {
                    indices.push(slot.index());
                }
            },
            || true
        ));
        assert_eq!(visited.len(), 4);
        assert_eq!(visited[..prefix.len()], prefix);
        let marker = parse("tail_marker");
        assert_eq!(visited.last().copied(), Some(marker.as_view()));
        let expected = parse("a");
        assert_eq!(indices, vec![expected.as_view(); 2]);
    }

    #[test]
    fn completed_index_observations_cover_remaining_scopes_and_prune_payloads() {
        crate::structure::representation::initialize();
        for source in [
            "scope(U(mink(4,a))*V(mink(4,a)),tail_marker)",
            "scope(T(mink(4,a))*(U(mink(4,a))+V(mink(4,b))),tail_marker)",
            "U(mink(4,a),mink(4,a))+V(mink(4,b))",
            "scope(U(mink(4,a))^2,tail_marker)",
            "U(mink(4,a),mink(4,a))^tail_marker",
            "x^U(mink(4,a),mink(4,a))",
            "scope(mink(4,a),mink(4,a),T(mink(4,opaque(hidden))),tail_marker)",
            "scope(mink(4,a),mink(4,a),T(mink(4,a,opaque(hidden))),tail_marker)",
        ] {
            let expr = parse(source);
            let mut expected = Vec::new();
            let mut slots = SlotMatcher::default();
            expr.visitor(&mut |node| {
                expected.push(node);
                matches!(slots.classify(node), SlotMatch::Other)
            });
            let mut observed = Vec::new();
            assert!(expr.has_repeated_explicit_indices_with_observer(
                |node, _| {
                    observed.push(node);
                },
                || true
            ));
            assert_eq!(observed, expected, "{source}");
            let mut prefix = Vec::new();
            assert!(expr.has_repeated_explicit_indices_with_observer(
                |node, _| prefix.push(node),
                || false
            ));
            assert!(observed.starts_with(&prefix), "{source}");
        }
    }
}
