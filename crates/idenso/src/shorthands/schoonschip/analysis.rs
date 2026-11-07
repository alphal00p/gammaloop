//! Share the syntactic checks that decide whether normalization has any work.

use spenso::{
    network::{library::symbolic::ETS, parsing::AtomStructureExt, tags::SPENSO_TAG},
    structure::slot::{SlotMatch, SlotMatcher},
};
use symbolica::atom::{AtomView, Symbol};

use crate::tensor::inference::InterfaceInference;

#[cfg(test)]
std::thread_local! {
    pub(super) static SCAN_COUNTS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    pub(crate) static CANDIDATE_NODE_VISITS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Call-local observations of one expression. Intrinsic dot rewrites may retain
/// the contraction candidate: they remove indices rather than introduce them.
/// Compound representation metadata leaves head observations incomplete; admitted
/// index payloads are opaque identities and do not contain tensor candidates.
#[derive(Clone, Debug)]
pub(crate) struct SimplificationCandidates<const N: usize> {
    pub(crate) repeated_indices: bool,
    pub(crate) brackets: bool,
    pub(crate) dots: bool,
    pub(crate) symbols: [bool; N],
    pub(crate) counts: [usize; N],
    pub(crate) complete: bool,
    // Separate from opaque metadata: an early observer may leave normal syntax
    // unvisited after proving only the repeated-index predicate.
    pub(crate) traversal_complete: bool,
    pub(crate) intrinsic: bool,
}

impl<const N: usize> SimplificationCandidates<N> {
    pub(crate) fn scan(
        expression: AtomView<'_>,
        symbols: [Symbol; N],
        on_first_repeat: impl FnOnce() -> bool,
    ) -> Self {
        Self::scan_observing(expression, symbols, on_first_repeat, |_, _| true)
    }

    /// Extend the same traversal with regional or dependency observations.
    /// Consumers must not start a second walk to recover these facts.
    pub(crate) fn scan_observing(
        expression: AtomView<'_>,
        symbols: [Symbol; N],
        on_first_repeat: impl FnOnce() -> bool,
        mut observe: impl for<'a> FnMut(AtomView<'a>, &SlotMatch<'a>) -> bool,
    ) -> Self {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::CandidateScan,
        );
        #[cfg(test)]
        SCAN_COUNTS.with(|count| count.set(count.get() + 1));
        let symbols = symbols.map(|symbol| symbol.get_id());
        let bracket = SPENSO_TAG.bracket.get_id();
        let rank_one = &SPENSO_TAG.rank1;
        let mut rank_one_heads = [None; 16];
        let mut vector_slots = None;
        let mut candidates = Self {
            repeated_indices: false,
            brackets: false,
            dots: false,
            symbols: [false; N],
            counts: [0; N],
            complete: true,
            traversal_complete: false,
            intrinsic: true,
        };
        let mut traversal_complete = true;
        let mut skipped = false;
        candidates.repeated_indices = expression.has_repeated_explicit_indices_with_observer(
            |node, slot| {
                #[cfg(test)]
                CANDIDATE_NODE_VISITS.with(|count| count.set(count.get() + 1));
                if !observe(node, slot) {
                    skipped = true;
                    candidates.allow_all();
                    return false;
                }
                match node {
                    AtomView::Fun(function) => {
                        let id = function.get_symbol_id();
                        candidates.observe_symbol(id, bracket, &symbols);
                        let entry =
                            &mut rank_one_heads[(id.wrapping_mul(0x9e37_79b9) >> 28) as usize];
                        let (tagged, intrinsic) = match *entry {
                            Some((cached, tagged, intrinsic)) if cached == id => {
                                (tagged, intrinsic)
                            }
                            _ => {
                                let tagged = function.get_symbol().has_tag(rank_one);
                                let intrinsic = InterfaceInference::intrinsic_normalization_head(
                                    function.get_symbol(),
                                );
                                *entry = Some((id, tagged, intrinsic));
                                (tagged, intrinsic)
                            }
                        };
                        candidates.intrinsic &= intrinsic;
                        if tagged && !candidates.dots {
                            let slots = vector_slots.get_or_insert_with(SlotMatcher::default);
                            candidates.dots =
                                !slots.vector_argument(function).is_some_and(|argument| {
                                    slots
                                        .compact_representation(argument)
                                        .is_some_and(|representation| representation.is_base())
                                        || matches!(slots.concrete_component(argument), Some(Ok(_)))
                                });
                        }
                        match slot {
                            SlotMatch::Explicit(slot) => {
                                // A variance wrapper can hide the representation
                                // head from the observer's outer function node.
                                if slot.representation().wrapper().is_some() {
                                    candidates.intrinsic &=
                                        InterfaceInference::intrinsic_normalization_head(
                                            slot.representation().head(),
                                        );
                                    candidates.observe_symbol(
                                        slot.representation().head().get_id(),
                                        bracket,
                                        &symbols,
                                    );
                                }
                                // An explicit slot owns its index as one identity. Neither
                                // its spelling nor heads inside it authorize tensor work.
                                if !matches!(slot.dimension(), AtomView::Var(_) | AtomView::Num(_))
                                {
                                    candidates.allow_all();
                                }
                            }
                            SlotMatch::Opaque => {
                                if function.iter().any(|argument| {
                                    !matches!(argument, AtomView::Var(_) | AtomView::Num(_))
                                }) {
                                    candidates.allow_all();
                                }
                            }
                            SlotMatch::Other => {}
                        }
                    }
                    AtomView::Var(variable) => {
                        candidates.observe_symbol(variable.get_symbol_id(), bracket, &symbols);
                    }
                    AtomView::Pow(power) => {
                        // Rank-one bases are observed as functions below. Only a
                        // metric base adds a power identity of its own; scalar
                        // powers still expose any eligible work in their children.
                        if !candidates.dots
                            && let AtomView::Fun(base) = power.get_base()
                        {
                            candidates.dots = base.get_symbol_id() == ETS.metric.get_id();
                        }
                    }
                    _ => {}
                }
                true
            },
            || {
                traversal_complete = on_first_repeat();
                traversal_complete
            },
        );
        candidates.traversal_complete = traversal_complete && !skipped;
        candidates
    }

    pub(crate) fn normalized(&self) -> bool {
        self.traversal_complete && !self.repeated_indices && !self.brackets && !self.dots
    }

    fn observe_symbol(&mut self, id: u32, bracket: u32, symbols: &[u32; N]) {
        if id == bracket {
            self.brackets = true;
        }
        if symbols.contains(&id) {
            for ((seen, count), &symbol) in
                self.symbols.iter_mut().zip(&mut self.counts).zip(symbols)
            {
                if id == symbol {
                    *seen = true;
                    *count += 1;
                }
            }
        }
    }

    fn allow_all(&mut self) {
        self.brackets = true;
        self.dots = true;
        self.complete = false;
        self.intrinsic = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shorthands::{bracket::BracketNormalizer, schoonschip::Schoonschip};
    use symbolica::{
        atom::{Atom, AtomCore},
        parser::ParseSettings,
    };

    #[test]
    fn early_observation_distinguishes_unvisited_tails_from_index_payloads() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        use symbolica::atom::FunctionBuilder;

        crate::test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let callback = spenso::tensor_symbol!(
            "analysis_early_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let vector = SPENSO_TAG.rank_one_tensor_symbol("spenso::analysis_early_vector");
        let slot = spenso::mink!(4, 83719);
        let tail = FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(FunctionBuilder::new(vector).add_arg(&slot).finish().pow(2))
            .finish();
        // Noncommutative function arguments put both tails strictly after the
        // repeated pair, independent of Add/Mul canonical ordering.
        let source = FunctionBuilder::new(symbolica::symbol!("analysis_early_scope"))
            .add_arg(&slot)
            .add_arg(&slot)
            .add_arg(FunctionBuilder::new(callback).add_arg(1).finish())
            .add_arg(tail)
            .finish();
        calls.store(0, Ordering::Relaxed);
        let early = SimplificationCandidates::scan(source.as_view(), [callback], || false);
        assert!(early.repeated_indices && early.complete && early.intrinsic);
        assert!(!early.traversal_complete && !early.normalized());
        assert!(!early.brackets && !early.dots && !early.symbols[0]);
        let full = SimplificationCandidates::scan(source.as_view(), [callback], || true);
        assert!(full.repeated_indices && full.traversal_complete && full.complete);
        assert!(full.brackets && full.dots && full.symbols[0]);
        assert!(!full.intrinsic);
        let mut decisions = 0;
        let resumed = SimplificationCandidates::scan(source.as_view(), [callback], || {
            decisions += 1;
            true
        });
        assert_eq!(decisions, 1);
        assert!(resumed.traversal_complete && resumed.complete && resumed.repeated_indices);
        assert_eq!(resumed.brackets, full.brackets);
        assert_eq!(resumed.dots, full.dots);
        assert_eq!(resumed.intrinsic, full.intrinsic);
        assert_eq!(resumed.symbols, full.symbols);
        assert_eq!(calls.load(Ordering::Relaxed), 0);

        let opaque =
            Atom::parse("T(mink(4,hidden(a)))", "spenso", ParseSettings::symbolica()).unwrap();
        let early = SimplificationCandidates::scan(opaque.as_view(), [], || false);
        assert!(!early.repeated_indices && early.traversal_complete);
        assert!(
            early.complete && early.intrinsic,
            "an index payload does not hide tensor syntax"
        );
    }

    #[test]
    fn intrinsic_observations_exclude_callbacks_and_hidden_payloads() {
        crate::test_support::test_initialize();
        let callback = spenso::tensor_symbol!("analysis_callback", norm = |_, _| {});
        let plain = Atom::parse(
            "g(mink(4,a),mink(4,b))*T(mink(4,a))",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        assert!(SimplificationCandidates::scan(plain.as_view(), [], || true).intrinsic);
        // The observer must continue after finding a repeated index, including
        // into scalar function metadata where normalization is still observable.
        let callback = symbolica::function!(callback, plain.clone());
        for expression in [callback.clone(), &plain * &callback] {
            assert!(!SimplificationCandidates::scan(expression.as_view(), [], || true).intrinsic);
        }
        for source in ["T(mink(f(x),a))", "T(dind(cof(f(x),a)))"] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            assert!(!SimplificationCandidates::scan(expression.as_view(), [], || true).intrinsic);
        }
    }

    #[test]
    fn scalar_powers_do_not_request_dot_normalization() {
        crate::test_support::test_initialize();
        let vector = SPENSO_TAG.rank_one_tensor_symbol("spenso::analysis_power_vector");
        assert!(vector.has_tag(&SPENSO_TAG.rank1));
        let square = Atom::parse(
            "analysis_power_vector(mink(4,a))^2",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        let AtomView::Pow(power) = square.as_view() else {
            panic!("expected a vector square");
        };
        let AtomView::Fun(base) = power.get_base() else {
            panic!("expected a vector function");
        };
        assert_eq!(base.get_symbol(), vector);
        let expected = Atom::parse(
            "g(analysis_power_vector(mink(4)),analysis_power_vector(mink(4)))",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        assert_ne!(square, expected);
        assert_eq!(square.normalize_dots(), expected);
        for source in ["(x+y)^8", "x^n", "unknown(x)^3", "(x+y)^(-2)"] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            let candidates = SimplificationCandidates::scan(expression.as_view(), [], || true);
            assert!(!candidates.dots, "{source}");
            assert_eq!(expression.normalize_dots(), expression, "{source}");
        }
        // Continue through bases, exponents, and the remainder after a repeated
        // index. Index payloads themselves do not expose tensor syntax.
        for source in [
            "g(mink(4,a),mink(4,b))^2",
            "g(mink(4,a),mink(4,b))^(-2)",
            "g(mink(4,a),mink(4,b))^(2/3)",
            "analysis_power_vector(mink(4,a))^2",
            "analysis_power_vector(mink(4,a))^(-3)",
            "(analysis_power_vector(mink(4,a))^2+x)^3",
            "x^(g(mink(4,a),mink(4,b))^2)",
            "scope(mink(4,a),mink(4,a),analysis_power_vector(mink(4,b))^2)",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            assert!(
                SimplificationCandidates::scan(expression.as_view(), [], || true).dots,
                "{source}"
            );
        }
    }

    #[test]
    fn completed_observations_preserve_heads_and_normalization_boundaries() {
        crate::test_support::test_initialize();
        let vector = SPENSO_TAG.rank_one_tensor_symbol("spenso::analysis_observed_vector");
        assert!(vector.has_tag(&SPENSO_TAG.rank1));
        let square = Atom::parse(
            "analysis_observed_vector(mink(4,a))^2",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        assert_ne!(square.normalize_dots(), square);
        let heads = [
            SPENSO_TAG.chain,
            SPENSO_TAG.trace,
            SPENSO_TAG.bracket,
            *crate::epsilon::EPSILON_SYMBOL,
        ];
        for source in [
            "g(mink(D,a),mink(D,b))",
            "g(mink(D,a),mink(D,b))*(g(mink(D,c),mink(D,d))+g(mink(D,c),mink(D,e)))",
            "T(dind(lor(4,a)))",
            "T(dind(mink(4)))",
            "T(mink(spenso::trace(4),a))",
            "g(analysis_observed_vector(mink(4)),analysis_observed_vector(mink(4)))",
            "analysis_observed_vector(mink(4,a))^2",
            "(x+y)^6*g(mink(D,a),mink(D,b))",
            "bracket(analysis_observed_vector(mink(4,a))^2)",
            "g(mink(4,a),mink(4,b))*T(mink(4,b))*later(spenso::trace)",
            "scope(mink(4,a),mink(4,a),bracket(analysis_observed_vector(mink(4,b))^2))",
            "scope(mink(4,a),mink(4,a),analysis_observed_vector(q(mink(4))))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            let candidates = SimplificationCandidates::scan(expression.as_view(), heads, || true);
            assert_eq!(
                candidates.repeated_indices,
                expression.has_repeated_explicit_indices(),
                "{source}"
            );
            if candidates.complete {
                assert_eq!(
                    candidates.symbols,
                    heads.map(|head| expression.contains_symbol(head)),
                    "{source}"
                );
            }
            if !candidates.dots {
                assert_eq!(expression.normalize_dots(), expression, "{source}");
            }
            if !candidates.brackets {
                assert_eq!(
                    BracketNormalizer::normalize(expression.as_view()),
                    expression,
                    "{source}"
                );
            }
        }
        // A wrapped representation is still visible to the old full head scan.
        use spenso::structure::representation::LibraryRep;
        let representation: LibraryRep = crate::representations::Bispinor {}.into();
        let head = representation.symbol();
        let expression =
            Atom::parse("T(dind(bis(4,a)))", "spenso", ParseSettings::symbolica()).unwrap();
        let candidates = SimplificationCandidates::scan(expression.as_view(), [head], || true);
        assert!(!candidates.complete || candidates.symbols == [expression.contains_symbol(head)]);
    }

    #[test]
    fn index_payloads_do_not_supply_candidates_or_callback_boundaries() {
        crate::test_support::test_initialize();
        let callback = symbolica::symbol!("analysis_index_callback", norm = |_, _| {});
        let heads = [SPENSO_TAG.trace, *crate::epsilon::EPSILON_SYMBOL, callback];
        let label = symbolica::function!(callback, heads[0], heads[1]);
        // Representation dimensions also name scalar values, not identity heads.
        let dimensional = symbolica::function!(
            spenso::tensor_symbol!("analysis_dimension_T"),
            spenso::mink!(heads[0], Atom::var(heads[1]))
        );
        let observed = SimplificationCandidates::scan(dimensional.as_view(), heads, || true);
        assert!(observed.complete && observed.intrinsic && observed.normalized());
        assert_eq!(observed.counts, [0; 3]);
        // The syntactic observer also keeps arbitrary written payloads opaque;
        // structure admission separately decides whether these are valid indices.
        for source in [
            "T(mink(spenso::chain,a))",
            "T(mink(4,spenso::trace))",
            "T(mink(4,a,spenso::bracket))",
            "T(mink(4,f(spenso::epsilon)))",
            "T(mink(4,opaque(analysis_power_vector(mink(4,a))^2)))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            let observed = SimplificationCandidates::scan(expression.as_view(), heads, || true);
            assert!(
                observed.complete && observed.intrinsic && observed.normalized(),
                "{source}"
            );
            assert_eq!(observed.counts, [0; 3], "{source}");
            assert_eq!(expression.normalize_dots(), expression, "{source}");
        }
        for index in [Atom::var(heads[0]), label] {
            for port in [
                spenso::mink!(4, &index),
                symbolica::function!(
                    spenso::structure::abstract_index::AIND_SYMBOLS.dind,
                    crate::cof!(3, &index)
                ),
            ] {
                let expression =
                    symbolica::function!(spenso::tensor_symbol!("analysis_index_T"), &port);
                let observed = SimplificationCandidates::scan(expression.as_view(), heads, || true);
                assert!(observed.complete && observed.traversal_complete && observed.intrinsic);
                assert!(observed.normalized());
                assert_eq!(observed.counts, [0; 3]);
            }
        }
    }

    #[test]
    fn compact_vectors_do_not_request_dot_normalization() {
        crate::test_support::test_initialize();
        let vector = SPENSO_TAG.rank_one_tensor_symbol("spenso::analysis_compact_vector");
        let callback = spenso::vector_symbol!(
            "spenso::analysis_compact_callback",
            norm = |_value, _out| {}
        );
        let heads = [vector, callback, SPENSO_TAG.trace];
        for source in [
            "analysis_compact_vector(mink(4))",
            "g(analysis_compact_vector(mink(4)),analysis_compact_callback(mink(4)))",
            "analysis_compact_vector(mink(4))^2",
            "(analysis_compact_vector(mink(4))+analysis_compact_callback(mink(4)))^3",
            "analysis_compact_vector(spenso::trace,mink(4))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            let candidates = SimplificationCandidates::scan(expression.as_view(), heads, || true);
            assert!(candidates.complete, "{source}");
            assert!(!candidates.dots, "{source}");
            assert_eq!(
                candidates.symbols,
                heads.map(|head| expression.contains_symbol(head)),
                "{source}"
            );
            assert_eq!(expression.normalize_dots(), expression, "{source}");
        }
    }

    #[test]
    fn compact_admission_preserves_nested_and_opaque_dot_candidates() {
        crate::test_support::test_initialize();
        let _ = SPENSO_TAG.rank_one_tensor_symbol("spenso::analysis_nested_vector");
        for source in [
            "analysis_nested_vector(analysis_nested_vector(mink(4)))",
            "analysis_nested_vector(mink(4,a))^2",
            "x^(analysis_nested_vector(mink(4,a))^2)",
            "(x+analysis_nested_vector(mink(4,a))^2)^3",
            "g(analysis_nested_vector(mink(4)),analysis_nested_vector(mink(4)))^2",
            "analysis_nested_vector(mink(4,a,b))",
            "analysis_nested_vector(dind(lor(4)))",
            "analysis_nested_vector(mink(D+1))",
            "analysis_nested_vector(meta(analysis_nested_vector(mink(4,a))^2),mink(4))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            assert!(
                SimplificationCandidates::scan(expression.as_view(), [], || true).dots,
                "{source}"
            );
        }
        let expression = Atom::parse(
            "analysis_nested_vector(mink(opaque(spenso::trace)))",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        let candidates =
            SimplificationCandidates::scan(expression.as_view(), [SPENSO_TAG.trace], || true);
        assert!(!candidates.complete);
        assert!(candidates.dots);
    }

    #[test]
    fn concrete_components_do_not_request_dot_normalization() {
        crate::test_support::test_initialize();
        SPENSO_TAG.rank_one_tensor_symbol("spenso::analysis_component_vector");
        let callback = spenso::vector_symbol!(
            "spenso::analysis_component_callback",
            norm = |_value, _out| {}
        );
        for source in [
            "analysis_component_vector(0,cind(0))",
            "analysis_component_vector(1,find(3))^2",
            "n*(1/(analysis_component_vector(0,cind(0))+x)+1/(analysis_component_vector(0,cind(0))+y))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            let candidates = SimplificationCandidates::scan(expression.as_view(), [], || true);
            assert!(candidates.complete && candidates.intrinsic, "{source}");
            assert!(!candidates.repeated_indices && !candidates.dots, "{source}");
            assert_eq!(expression.normalize_dots(), expression, "{source}");
        }
        for source in [
            "analysis_component_vector(cind())",
            "analysis_component_vector(cind(-1))",
            "analysis_component_vector(cind(0,1))",
            "analysis_component_vector(find(x))",
            "analysis_component_vector(meta(analysis_component_vector(mink(4,a))^2),cind(0))",
            "analysis_component_vector(analysis_component_vector(mink(4)))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            assert!(
                SimplificationCandidates::scan(expression.as_view(), [], || true).dots,
                "{source}"
            );
        }
        let expression = Atom::parse(
            "analysis_component_callback(cind(0))",
            "spenso",
            ParseSettings::symbolica(),
        )
        .unwrap();
        let candidates = SimplificationCandidates::scan(expression.as_view(), [callback], || true);
        assert!(!candidates.intrinsic);
        assert_eq!(candidates.symbols, [true]);
    }
}
