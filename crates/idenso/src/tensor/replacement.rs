//! Whole tensor-leaf substitutions with locally checked interfaces.

use super::composition::ExplicitIndexOccurrences;
use super::inference::{InterfaceInference, TensorInferenceError};
use super::{SymbolicTensor, composition};
use spenso::structure::partial::PartialStructure;
use spenso::structure::slot::SlotMatcher;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use symbolica::atom::{Atom, AtomOrView, AtomView, Symbol};
use symbolica::id::{
    AtomMatchIterator, Condition, Match, MatchSettings, Pattern, PatternRestriction, ReplaceWith,
    WrappedMatchStack,
};
use symbolica::state::Workspace;
use symbolica::transformer::TransformerError;

type Result<T> = std::result::Result<T, TensorInferenceError>;

/// A reusable whole-leaf tensor replacement.
///
/// Patterns, conditions and wildcard closure belong to the rule. A literal RHS
/// can also carry its intrinsic interface proof. Wildcards which bind whole
/// slots, argument slices or function heads still require the existing checks
/// after substitution; a rule does not certify those unknown interfaces.
/// Matching state and RHS callback caches are local to each application.
#[derive(Clone)]
pub struct TensorRule<'rhs> {
    pattern: Pattern,
    rhs: ReplaceWith<'rhs>,
    conditions: Condition<PatternRestriction>,
    settings: MatchSettings,
    rhs_cache_size: usize,
    fixed_head: Option<u32>,
    literal_rhs: Option<Arc<Signature>>,
}

impl std::fmt::Debug for TensorRule<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Pattern Display reconstructs an Atom; Debug preserves the stored
        // pattern and never invokes normalization or replacement callbacks.
        write!(
            f,
            "TensorRule({:?} -> {:?}, rhs_cache_size={})",
            self.pattern, self.rhs, self.rhs_cache_size,
        )
    }
}

impl<'rhs> TensorRule<'rhs> {
    pub fn new(
        pattern: Pattern,
        rhs: ReplaceWith<'rhs>,
        conditions: Option<Condition<PatternRestriction>>,
        rhs_cache_size: usize,
    ) -> Result<Self> {
        if let ReplaceWith::Pattern(rhs) = &rhs {
            let mut available = HashSet::new();
            let mut required = HashSet::new();
            Self::wildcards(&pattern, &mut available);
            Self::wildcards(rhs.borrow(), &mut required);
            if let Some(wildcard) = required.difference(&available).min_by_key(|s| s.get_id()) {
                return Err(TensorInferenceError::invalid(format!(
                    "replacement wildcard `{wildcard}` does not occur in the tensor pattern"
                )));
            }
        }
        let literal_rhs = match &rhs {
            ReplaceWith::Pattern(rhs) => match rhs.borrow() {
                Pattern::Literal(atom)
                    if !atom.as_view().needs_normalization()
                        && InterfaceInference::default()
                            .rewrites_preserve_leaf_interfaces(atom.as_view()) =>
                {
                    // Invalid or callback-sensitive RHS expressions still fail
                    // only when a match uses them, as in ordinary application.
                    Signature::observe(atom.as_view(), &mut SlotMatcher::default())
                        .ok()
                        .map(Arc::new)
                }
                _ => None,
            },
            ReplaceWith::Map(_) => None,
        };
        let fixed_head = match &pattern {
            Pattern::Fn(head, _) if head.get_wildcard_level() == 0 => Some(head.get_id()),
            _ => None,
        };
        Ok(Self {
            pattern,
            rhs,
            conditions: conditions.unwrap_or(Condition::True),
            settings: MatchSettings::new()
                .partial(false)
                .min_level(0)
                .max_level(0),
            rhs_cache_size,
            fixed_head,
            literal_rhs,
        })
    }

    // Symbolica's public pattern visitor does not descend into alternatives or
    // count wildcard function heads. This validates closure only; matching and
    // substitution remain entirely owned by Symbolica.
    fn wildcards(pattern: &Pattern, result: &mut HashSet<Symbol>) {
        match pattern {
            Pattern::Wildcard(symbol, _) => {
                result.insert(*symbol);
            }
            Pattern::Fn(head, arguments) => {
                if head.get_wildcard_level() > 0 {
                    result.insert(*head);
                }
                for argument in arguments {
                    Self::wildcards(argument, result);
                }
            }
            Pattern::Pow(arguments) => {
                for argument in arguments.iter() {
                    Self::wildcards(argument, result);
                }
            }
            Pattern::Mul(arguments) | Pattern::Add(arguments) | Pattern::Alternative(arguments) => {
                for argument in arguments {
                    Self::wildcards(argument, result);
                }
            }
            Pattern::Transformer(transformer) => {
                if let Some(input) = &transformer.0 {
                    Self::wildcards(input, result);
                }
            }
            Pattern::Literal(_) => {}
        }
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Apply a reusable rule to whole tensor leaves underneath sums and products.
    ///
    /// Tensor arguments and scalar metadata are opaque. Each reduced RHS must
    /// retain the matched leaf's external ports without adding encoded indices
    /// or increasing their multiplicities. Zero retains the typed interface.
    /// Sources with unresolved ports, user normalization hooks, or unsupported
    /// tensor powers are rejected. Raw RHS expressions must also have normalized
    /// intrinsic leaves and valid encoded interfaces.
    /// Conditions run for every candidate match; RHS callbacks run only on cache
    /// misses. Set the rule's `rhs_cache_size` to zero for side-effecting callbacks.
    /// As with the other certified algebra operations, `self` must already have
    /// a valid established interface and explicit-index multiplicities.
    ///
    /// # Errors
    ///
    /// Rejects unsupported source syntax, invalid source metadata, malformed rule bindings,
    /// or a replacement that changes the matched interface or introduces index uses.
    ///
    /// # Example
    ///
    /// ```
    /// use std::slice;
    /// use idenso::tensor::{SymbolicTensor, TensorRule};
    /// use idenso::tensor::inference::TensorInferenceError;
    /// use spenso::structure::abstract_index::AbstractIndex;
    /// use spenso::structure::representation::{Euclidean, RepName};
    /// use spenso::structure::slot::IsAbstractSlot;
    /// use spenso::tensor_symbol;
    /// use symbolica::atom::AtomCore;
    /// use symbolica::function;
    ///
    /// let slot = Euclidean {}
    ///     .new_rep(3)
    ///     .to_lib()
    ///     .slot::<AbstractIndex, _>(AbstractIndex::Normal(7));
    /// let from = function!(tensor_symbol!("docs::rule_from"), slot.to_atom());
    /// let to = function!(tensor_symbol!("docs::rule_to"), slot.to_atom());
    /// let source = SymbolicTensor::infer(from.clone())?.finish_parts()?;
    /// let rule = TensorRule::new(from.to_pattern(), to.clone().into(), None, 0)?;
    /// let replaced = source.replace(&rule)?;
    /// assert_eq!(replaced.expression(), &to);
    /// assert!(source.replace_rules(slice::from_ref(&rule))?.same_value(&replaced));
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn replace(&self, rule: &TensorRule<'_>) -> Result<Self> {
        self.replace_rules(std::slice::from_ref(rule))
    }

    /// Apply an ordered rule set in one traversal of the original tensor leaves.
    /// The first matching rule owns a leaf, including identity replacements;
    /// emitted RHS leaves are not matched again in the same stage. Conditions
    /// and bounded RHS caches remain local to each rule and this application.
    /// See the [single-rule example](Self::replace), including ordered rule application.
    ///
    /// # Errors
    ///
    /// Returns the same admission and per-binding errors as [`Self::replace`].
    pub fn replace_rules(&self, rules: &[TensorRule<'_>]) -> Result<Self> {
        if rules.is_empty() {
            return Ok(self.clone());
        }
        if !self.rewrites_preserve_interface() {
            return Err(TensorInferenceError::invalid(
                "tensor-safe replacement requires normalized intrinsic tensor leaves with explicit ports and supported powers",
            ));
        }
        self.explicit_index_occurrences()?;
        let mut replacement = TensorReplacement {
            rules: rules.iter().map(RuleApplication::new).collect(),
            source_signatures: HashMap::new(),
            slots: SlotMatcher::default(),
        };
        match replacement.apply(self.expression.as_view())? {
            AtomOrView::View(_) => Ok(self.clone()),
            expression => {
                let mut result = Self::from_validated_parts(
                    expression.into_owned(),
                    self.structure.clone(),
                    Arc::clone(
                        self.proofs
                            .occurrences
                            .get()
                            .expect("published index history"),
                    ),
                );
                // Every binding preserves ports and bounds its explicit indices
                // by the source. Retain that conservative metadata even when
                // normalization cancels leaves or coalesces them into powers.
                result.proofs = self.proofs.algebra_result();
                Ok(result)
            }
        }
    }
}

// A signature describes encoded syntax as well as its external interface.
// Free ports alone cannot detect an internal pair colliding with the context.
struct Signature {
    interface: PartialStructure,
    occurrences: ExplicitIndexOccurrences,
}

impl Signature {
    fn observe(value: AtomView<'_>, slots: &mut SlotMatcher) -> Result<Self> {
        #[cfg(test)]
        tests::SIGNATURE_OBSERVATIONS.with(|count| count.set(count.get() + 1));
        let (interface, occurrences) = InterfaceInference::replacement_observation(value, slots)?;
        Ok(Self {
            interface,
            occurrences,
        })
    }

    fn accepts(&self, rhs: &Self, zero: bool) -> bool {
        (zero || InterfaceInference::additive_interfaces_match(&self.interface, &rhs.interface))
            && rhs.occurrences.is_bounded_by(&self.occurrences)
    }
}

type Bindings<'source> = Vec<(Symbol, Match<'source>)>;

struct RuleApplication<'source, 'rule, 'rhs> {
    rule: &'rule TensorRule<'rhs>,
    matcher: Option<AtomMatchIterator<'source, 'rule>>,
    match_stack: WrappedMatchStack<'source, 'rule>,
    rhs_cache: HashMap<Bindings<'source>, (Atom, Arc<Signature>)>,
}

impl<'source, 'rule, 'rhs> RuleApplication<'source, 'rule, 'rhs> {
    fn new(rule: &'rule TensorRule<'rhs>) -> Self {
        Self {
            rule,
            matcher: None,
            match_stack: WrappedMatchStack::new(&rule.conditions, &rule.settings),
            rhs_cache: HashMap::new(),
        }
    }
}

struct TensorReplacement<'source, 'rule, 'rhs> {
    rules: Vec<RuleApplication<'source, 'rule, 'rhs>>,
    source_signatures: HashMap<AtomView<'source>, Signature>,
    slots: SlotMatcher,
}

impl<'source> TensorReplacement<'source, '_, '_> {
    fn apply(&mut self, value: AtomView<'source>) -> Result<AtomOrView<'source>> {
        match value {
            AtomView::Add(sum) => {
                let terms = sum
                    .iter()
                    .map(|term| self.apply(term))
                    .collect::<Result<Vec<_>>>()?;
                if terms.iter().all(|term| matches!(term, AtomOrView::View(_))) {
                    return Ok(value.into());
                }
                Ok(Atom::add_many(terms.iter().map(AtomOrView::as_view)).into())
            }
            AtomView::Mul(product) => {
                let factors = product
                    .iter()
                    .map(|factor| self.apply(factor))
                    .collect::<Result<Vec<_>>>()?;
                if factors
                    .iter()
                    .all(|factor| matches!(factor, AtomOrView::View(_)))
                {
                    return Ok(value.into());
                }
                Ok(Atom::mul_many(factors.iter().map(AtomOrView::as_view)).into())
            }
            AtomView::Fun(function)
                if self.rules.iter().any(|application| {
                    application
                        .rule
                        .fixed_head
                        .is_none_or(|head| head == function.get_symbol_id())
                }) && composition::is_tensor_leaf_head(function.get_symbol())
                    && !function.get_symbol().is_scalar() =>
            {
                self.replace_leaf(value)
            }
            // A tensor inside metadata, a bracket, or a power is not a whole
            // leaf of this operation's transparent sum/product context.
            _ => Ok(value.into()),
        }
    }

    fn replace_leaf(&mut self, value: AtomView<'source>) -> Result<AtomOrView<'source>> {
        let mut matched = None;
        for application in &mut self.rules {
            if let AtomView::Fun(function) = value
                && application
                    .rule
                    .fixed_head
                    .is_some_and(|head| head != function.get_symbol_id())
            {
                continue;
            }
            application.match_stack.reset();
            // Compile lazily, retaining no-match behavior for unsupported
            // patterns whose stage never reaches an eligible tensor leaf.
            let matcher = application
                .matcher
                .get_or_insert_with(|| AtomMatchIterator::new(&application.rule.pattern));
            matcher.set_new_target(value, &application.match_stack);
            if let Some(used_flags) = matcher.next(&mut application.match_stack) {
                debug_assert!(used_flags.iter().all(|used| *used));
                matched = Some(application);
                break;
            }
        }
        let Some(application) = matched else {
            return Ok(value.into());
        };
        let found = application.match_stack.get_match_stack();
        // This bounded borrowed-key cache saves repeated source observation;
        // its limit is independent of the user-selected RHS callback cache.
        let cache_source = self.source_signatures.len() < 256;
        let uncached;
        let source = match self.source_signatures.entry(value) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let signature = Signature::observe(value, &mut self.slots)?;
                if cache_source {
                    entry.insert(signature)
                } else {
                    uncached = signature;
                    &uncached
                }
            }
        };
        let key = found.get_matches();
        if let Some((rhs, signature)) = application.rhs_cache.get(key) {
            if !source.accepts(signature, rhs.is_zero()) {
                return Err(TensorInferenceError::invalid(
                    "tensor replacement changes the matched interface or introduces explicit indices",
                ));
            }
            return Ok(if rhs.as_view() == value {
                value.into()
            } else {
                rhs.clone().into()
            });
        }
        let rhs = match &application.rule.rhs {
            ReplaceWith::Pattern(pattern) => Workspace::get_local().with(|workspace| {
                let mut result = Atom::new();
                pattern
                    .replace_wildcards_with_matches_impl(workspace, &mut result, found, false, None)
                    .map_err(|error| match error {
                        TransformerError::ValueError(message) => {
                            TensorInferenceError::invalid(message)
                        }
                        TransformerError::Interrupt => TensorInferenceError::Interrupted,
                    })?;
                Ok::<_, TensorInferenceError>(result)
            })?,
            ReplaceWith::Map(map) => map(found),
        };
        let signature = if let Some(signature) = &application.rule.literal_rhs {
            Arc::clone(signature)
        } else {
            if rhs.as_view().needs_normalization()
                || !InterfaceInference::normalization_is_intrinsic(rhs.as_view())
            {
                return Err(TensorInferenceError::invalid(
                    "tensor replacement RHS requires normalized intrinsic leaves",
                ));
            }
            Arc::new(Signature::observe(rhs.as_view(), &mut self.slots)?)
        };
        if !source.accepts(&signature, rhs.is_zero()) {
            return Err(TensorInferenceError::invalid(
                "tensor replacement changes the matched interface or introduces explicit indices",
            ));
        }
        if application.rhs_cache.len() < application.rule.rhs_cache_size {
            application
                .rhs_cache
                .insert(key.to_vec(), (rhs.clone(), signature));
        }
        Ok(if rhs.as_view() == value {
            value.into()
        } else {
            rhs.into()
        })
    }
}

#[cfg(test)]
mod tests {

    use super::super::composition::TensorCompositionError;
    use super::super::inference::tests::{
        INFERENCE_CALLS, OCCURRENCE_FUNCTION_VISITS, SCOPE_VALIDATIONS,
    };
    use super::*;
    use crate::test_support;
    use spenso::network::tags::SPENSO_TAG;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::dimension::Dimension;
    use spenso::structure::partial::{PartialIndex, PartialStructureExt};
    use spenso::structure::representation::{ExtendibleReps, RepName};
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::structure::{HasStructure, OrderedStructure, TensorStructure};
    use spenso::tensor_symbol;
    use std::cell::Cell;
    use std::slice;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};
    use symbolica::atom::{AtomCore, FunctionBuilder};
    use symbolica::id::ConditionResult;
    use symbolica::symbol;

    thread_local! {
        pub(super) static SIGNATURE_OBSERVATIONS: Cell<usize> = const { Cell::new(0) };
    }

    fn leaf(head: Symbol, args: &[Atom]) -> Atom {
        FunctionBuilder::new(head).add_args(args).finish()
    }

    fn slot(index: usize) -> Atom {
        ExtendibleReps::MINKOWSKI
            .new_rep(Dimension::Concrete(4))
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    }

    fn typed(expression: Atom) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor::<PartialStructure>::infer(expression)
            .unwrap()
            .finish_parts()
            .unwrap()
    }

    #[test]
    fn rule_stage_matches_original_leaves_once_in_declared_priority() {
        let f = spenso::tensor_symbol!("rule_stage_f");
        let g = spenso::tensor_symbol!("rule_stage_g");
        let h = spenso::tensor_symbol!("rule_stage_h");
        let a = slot(81871);
        let f = leaf(f, std::slice::from_ref(&a));
        let g = leaf(g, std::slice::from_ref(&a));
        let h = leaf(h, &[a]);
        let rules = [
            TensorRule::new(f.to_pattern(), g.clone().into(), None, 10).unwrap(),
            TensorRule::new(g.to_pattern(), h.clone().into(), None, 10).unwrap(),
        ];
        let source = typed(&f + &g);
        let result = source.replace_rules(&rules).unwrap();
        assert_eq!(result.expression, &g + &h);
        assert_eq!(result.structure, source.structure);
        // A matching identity rule still prevents lower-priority application.
        let rules = [
            TensorRule::new(f.to_pattern(), f.clone().into(), None, 0).unwrap(),
            TensorRule::new(f.to_pattern(), h.into(), None, 0).unwrap(),
        ];
        assert_eq!(
            typed(f.clone()).replace_rules(&rules).unwrap().expression,
            f
        );
    }

    #[test]
    fn rule_stage_conditions_and_callback_caches_are_per_rule() {
        let f = spenso::tensor_symbol!("rule_stage_callback_f");
        let g = spenso::tensor_symbol!("rule_stage_callback_g");
        let h = spenso::tensor_symbol!("rule_stage_callback_h");
        let wildcard = symbolica::symbol!("rule_stage_callback_arg_");
        let a = slot(81881);
        let first = leaf(f, std::slice::from_ref(&a));
        let second = leaf(g, std::slice::from_ref(&a));
        let rhs = leaf(h, &[a]);
        let x = Atom::var(symbolica::symbol!("rule_stage_callback_x"));
        let y = Atom::var(symbolica::symbol!("rule_stage_callback_y"));
        let input = typed(&x * (&first + &second) + &y * (&first + &second));
        let calls = Arc::new(Mutex::new(Vec::new()));
        let rules = [f, g].map(|head| {
            let calls = Arc::clone(&calls);
            let rhs = rhs.clone();
            TensorRule::new(
                leaf(head, &[Atom::var(wildcard)]).to_pattern(),
                ReplaceWith::Map(Box::new(move |_| {
                    calls.lock().unwrap().push(head);
                    rhs.clone()
                })),
                None,
                1,
            )
            .unwrap()
        });
        let result = input.replace_rules(&rules).unwrap();
        assert_eq!(result.expression, 2 * &x * &rhs + 2 * &y * rhs);
        let calls = calls.lock().unwrap();
        assert_eq!(calls.len(), 2);
        assert!(calls.contains(&f) && calls.contains(&g));
    }

    #[test]
    fn tensor_rule_debug_does_not_normalize_patterns_or_call_rhs() {
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = calls.clone();
        let head = symbolica::symbol!("rule_debug_callback"; Scalar; norm=move|_,_| {
            seen.fetch_add(1, Ordering::Relaxed);
        });
        let seen = calls.clone();
        let rule = TensorRule::new(
            Pattern::Fn(head, vec![Pattern::Literal(Atom::num(3))]),
            ReplaceWith::Map(Box::new(move |_| {
                seen.fetch_add(1, Ordering::Relaxed);
                Atom::Zero
            })),
            None,
            0,
        )
        .unwrap();
        let first = format!("{rule:?}");
        assert_eq!(first, format!("{rule:?}"));
        assert!(first.starts_with("TensorRule(Fn("));
        assert!(first.ends_with(" -> Map, rhs_cache_size=0)"));
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn tensor_rule_reuses_literal_certification_without_caching_target_compatibility() {
        let f = spenso::tensor_symbol!("rule_literal_f");
        let g = spenso::tensor_symbol!("rule_literal_g");
        let a = slot(81901);
        let source = typed(leaf(f, std::slice::from_ref(&a)));
        let rhs = leaf(g, std::slice::from_ref(&a));
        SIGNATURE_OBSERVATIONS.with(|count| count.set(0));
        let rule =
            TensorRule::new(source.expression.to_pattern(), rhs.clone().into(), None, 0).unwrap();
        assert!(rule.literal_rhs.is_some());
        assert_eq!(SIGNATURE_OBSERVATIONS.with(Cell::get), 1);
        for _ in 0..2 {
            let result = source.replace(&rule).unwrap();
            assert_eq!(result.expression, rhs);
            assert_eq!(result.structure, source.structure);
        }
        assert_eq!(SIGNATURE_OBSERVATIONS.with(Cell::get), 3);

        let extra = slot(81903);
        let wide = leaf(f, &[a, extra]);
        let pattern = Pattern::Alternative(vec![source.expression.to_pattern(), wide.to_pattern()]);
        let rule = TensorRule::new(pattern, rhs.into(), None, 100).unwrap();
        // Both literal alternatives bind nothing, so the second use hits the
        // RHS cache. Its additional port must still cause rejection.
        let source = typed(source.expression * wide);
        assert!(source.replace(&rule).is_err());
    }

    #[test]
    fn tensor_rule_checks_closure_in_alternatives_and_function_heads() {
        let f = spenso::tensor_symbol!("rule_closure_f");
        let g = spenso::tensor_symbol!("rule_closure_g");
        let h = spenso::tensor_symbol!("rule_closure_h");
        let argument = symbolica::symbol!("rule_closure_argument_");
        let head = symbolica::symbol!("rule_closure_head_");
        let missing = symbolica::symbol!("rule_closure_missing_");
        let call = |head| Pattern::Fn(head, vec![Pattern::Wildcard(argument, false)]);
        let alternatives = Pattern::Alternative(vec![call(f), call(g)]);
        let rule = TensorRule::new(alternatives.clone(), call(h).into(), None, 10).unwrap();
        assert!(rule.literal_rhs.is_none());
        let a = slot(81911);
        let source = typed(leaf(f, std::slice::from_ref(&a)) + leaf(g, std::slice::from_ref(&a)));
        assert_eq!(source.replace(&rule).unwrap().expression, 2 * leaf(h, &[a]));
        assert!(TensorRule::new(alternatives, call(missing).into(), None, 10).is_err());
        assert!(TensorRule::new(call(head), call(head).into(), None, 10).is_ok());
        assert!(TensorRule::new(call(head), call(missing).into(), None, 10).is_err());
    }

    #[test]
    fn tensor_rule_keeps_callback_caches_local_to_each_application() {
        let f = spenso::tensor_symbol!("rule_callback_f");
        let g = spenso::tensor_symbol!("rule_callback_g");
        let argument = symbolica::symbol!("rule_callback_argument_");
        let a = slot(81921);
        let input = leaf(f, std::slice::from_ref(&a));
        let output = leaf(g, &[a]);
        let x = Atom::var(symbolica::symbol!("rule_callback_x"));
        let y = Atom::var(symbolica::symbol!("rule_callback_y"));
        let source = typed(&x * &input + &y * &input);
        let mut schedules = Vec::new();
        for (cache_size, expected_calls) in [(0, 4), (1, 2)] {
            let calls = Arc::new(AtomicUsize::new(0));
            let conditions = Arc::new(AtomicUsize::new(0));
            let seen = Arc::clone(&calls);
            let rhs = output.clone();
            let seen_condition = Arc::clone(&conditions);
            let rule = TensorRule::new(
                leaf(f, &[Atom::var(argument)]).to_pattern(),
                ReplaceWith::Map(Box::new(move |_| {
                    seen.fetch_add(1, Ordering::Relaxed);
                    rhs.clone()
                })),
                Some(Condition::match_stack(move |_| {
                    seen_condition.fetch_add(1, Ordering::Relaxed);
                    ConditionResult::True
                })),
                cache_size,
            )
            .unwrap();
            assert_eq!(calls.load(Ordering::Relaxed), 0);
            assert_eq!(conditions.load(Ordering::Relaxed), 0);
            for _ in 0..2 {
                assert_eq!(
                    source.replace(&rule).unwrap().expression,
                    &x * &output + &y * &output
                );
            }
            assert_eq!(calls.load(Ordering::Relaxed), expected_calls);
            schedules.push(conditions.load(Ordering::Relaxed));
        }
        assert!(schedules[0] >= 4);
        assert_eq!(schedules[0], schedules[1]);
    }

    #[test]
    fn tensor_replacement_preserves_bound_ports_logical_order_and_typed_zero() {
        let t = spenso::tensor_symbol!("safe_replace_bound_t");
        let u = spenso::tensor_symbol!("safe_replace_bound_u");
        let p = spenso::vector_symbol!("safe_replace_bound_p");
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let compact = leaf(p, &[m.to_symbolic([])]);
        let metadata = leaf(
            symbolica::symbol!("safe_replace_routing"; Scalar),
            std::slice::from_ref(&compact),
        );
        let a = slot(81101);
        let b = slot(81103);
        for args in [
            vec![metadata.clone(), compact.clone(), a.clone(), b.clone()],
            vec![metadata.clone(), a.clone(), compact.clone(), b.clone()],
            vec![metadata.clone(), a.clone(), b.clone(), compact.clone()],
        ] {
            let expression = leaf(t, &args);
            let mut source = typed(expression.clone());
            source.structure = PartialStructure::from_logical_slots(
                source.structure.logical_slots().into_iter().rev(),
            );
            let pattern = expression.to_pattern();
            let result = source
                .replace(
                    &TensorRule::new(pattern.clone(), leaf(u, &args).into(), None, 10).unwrap(),
                )
                .unwrap();
            assert_eq!(result.expression, leaf(u, &args));
            assert_eq!(result.structure, source.structure);
            let zero = source
                .replace(&TensorRule::new(pattern.clone(), Atom::Zero.into(), None, 10).unwrap())
                .unwrap();
            assert!(zero.expression.is_zero());
            assert_eq!(zero.structure, source.structure);
        }
    }

    #[test]
    fn tensor_replacement_checks_each_branch_and_forbids_new_internal_indices() {
        let f = spenso::tensor_symbol!("safe_replace_branch_f");
        let g = spenso::tensor_symbol!("safe_replace_branch_g");
        let h = spenso::tensor_symbol!("safe_replace_branch_h");
        let a = slot(81201);
        let b = slot(81203);
        let fa = leaf(f, std::slice::from_ref(&a));
        let ga = leaf(g, std::slice::from_ref(&a));
        let ha = leaf(h, std::slice::from_ref(&a));
        let source = typed(fa.clone());
        let pattern = fa.to_pattern();
        let same = &ga + &ha;
        assert_eq!(
            source
                .replace(&TensorRule::new(pattern.clone(), same.clone().into(), None, 0).unwrap())
                .unwrap()
                .expression,
            same
        );
        assert!(
            source
                .replace(
                    &TensorRule::new(
                        pattern.clone(),
                        (ga + leaf(h, std::slice::from_ref(&b))).into(),
                        None,
                        0
                    )
                    .unwrap()
                )
                .is_err()
        );
        let diagonal = leaf(f, &[a.clone(), a]);
        let scalar = typed(diagonal.clone());
        assert!(scalar.is_scalar());
        assert!(
            scalar
                .replace(
                    &TensorRule::new(
                        diagonal.to_pattern(),
                        (leaf(g, std::slice::from_ref(&b)) * leaf(h, &[b])).into(),
                        None,
                        0
                    )
                    .unwrap()
                )
                .is_err()
        );
        assert_eq!(
            scalar
                .replace(
                    &TensorRule::new(diagonal.to_pattern(), Atom::one().into(), None, 0).unwrap()
                )
                .unwrap()
                .expression,
            Atom::one()
        );
    }

    #[test]
    fn tensor_replacement_cancellation_retains_source_dummy_history_without_output_observation() {
        test_support::test_initialize();
        let f = tensor_symbol!("rule_cancellation_f");
        let g = tensor_symbol!("rule_cancellation_g");
        let index = AbstractIndex::Normal(81321);
        let port = ExtendibleReps::MINKOWSKI
            .new_rep(Dimension::Concrete(4))
            .slot(index);
        let a = port.to_atom();
        let diagonal = leaf(f, &[a.clone(), a.clone()]);
        let x = Atom::var(symbol!("rule_cancellation_x"));
        let admitted = typed(&diagonal * &x);
        // Publication must carry history before replacement can erase its payload.
        let history = Arc::clone(admitted.proofs.occurrences.get().unwrap());
        let source =
            SymbolicTensor::from_validated_parts(admitted.expression, admitted.structure, history);
        assert!(source.proofs.occurrences.get().is_some());
        let collision = typed(leaf(g, &[a]));
        let rule = TensorRule::new(diagonal.to_pattern(), x.pow(-1).into(), None, 0).unwrap();
        assert!(rule.literal_rhs.is_some());
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
        INFERENCE_CALLS.with(|count| count.set(0));
        SIGNATURE_OBSERVATIONS.with(|count| count.set(0));

        let result = source.replace(&rule).unwrap();
        assert_eq!(result.expression, Atom::one());
        assert_eq!(result.structure, source.structure);
        // Only the matched leaf requires a signature; the normalized output
        // inherits admission and shares the completed source metadata.
        assert_eq!(SIGNATURE_OBSERVATIONS.with(Cell::get), 1);
        assert_eq!(SCOPE_VALIDATIONS.with(Cell::get), 1);
        assert_eq!(INFERENCE_CALLS.with(Cell::get), 0);
        let visits = OCCURRENCE_FUNCTION_VISITS.with(Cell::get);
        let cloned = result.clone();
        for value in [&result, &cloned] {
            assert!(Arc::ptr_eq(
                source.proofs.occurrences.get().unwrap(),
                value.proofs.occurrences.get().unwrap(),
            ));
            assert_eq!(value.explicit_index_occurrences().unwrap().count(port), 2);
            assert!(matches!(
                value.multiply(&collision),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
        assert_eq!(SCOPE_VALIDATIONS.with(Cell::get), 1);
        assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(Cell::get), visits);
        assert_eq!(INFERENCE_CALLS.with(Cell::get), 0);
    }

    #[test]
    fn tensor_replacement_carries_coalesced_powers_and_rejects_source_powers_before_rhs() {
        let f = tensor_symbol!("safe_replace_power_f");
        let g = tensor_symbol!("safe_replace_power_g");
        let a = slot(81301);
        let fa = leaf(f, slice::from_ref(&a));
        let ga = leaf(g, &[a]);
        let source = typed(&fa * &ga);
        let rule = TensorRule::new(fa.to_pattern(), ga.clone().into(), None, 0).unwrap();
        assert!(rule.literal_rhs.is_some());
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        INFERENCE_CALLS.with(|count| count.set(0));
        let result = source.replace(&rule).unwrap();
        assert_eq!(result.expression, ga.clone().pow(2));
        assert!(result.is_scalar());
        assert!(Arc::ptr_eq(
            source.proofs.occurrences.get().unwrap(),
            result.proofs.occurrences.get().unwrap(),
        ));
        assert_eq!(SCOPE_VALIDATIONS.with(Cell::get), 1);
        assert_eq!(INFERENCE_CALLS.with(Cell::get), 0);
        assert!(result.proofs.rewrite.get().is_none());
        assert!(result.proofs.algebra.get().is_none());
        assert!(!result.rewrites_preserve_interface());
        assert!(
            SymbolicTensor::<PartialStructure>::new(&fa * ga.pow(2), source.structure).is_err()
        );
    }

    #[test]
    fn tensor_replacement_checks_powered_rhs_without_certifying_leaf_rewrites() {
        crate::test_support::test_initialize();
        let head = spenso::tensor_symbol!("rule_powered_rhs_vertex");
        let q = spenso::tensor_symbol!("rule_powered_rhs_q");
        let a = slot(81351);
        let vertex = leaf(head, &[a.clone(), a.clone()]);
        let qa = leaf(q, std::slice::from_ref(&a));
        let source = typed(vertex.clone());
        let square = qa.clone().pow(2);
        let result = source
            .replace(&TensorRule::new(vertex.to_pattern(), square.clone().into(), None, 8).unwrap())
            .unwrap();
        assert_eq!(result.expression(), &square);
        assert!(result.is_scalar());
        assert!(!result.rewrites_preserve_interface());
        for invalid in [qa.pow(3), leaf(q, &[slot(81352)]).pow(2)] {
            assert!(
                source
                    .replace(
                        &TensorRule::new(vertex.to_pattern(), invalid.into(), None, 8).unwrap()
                    )
                    .is_err()
            );
        }
        let open = leaf(head, &[a]);
        assert!(
            typed(open.clone())
                .replace(&TensorRule::new(open.to_pattern(), square.into(), None, 8).unwrap())
                .is_err()
        );
    }

    #[test]
    fn tensor_replacement_cache_keeps_conditions_and_identity_flags() {
        let f = spenso::tensor_symbol!("safe_replace_cache_f");
        let fa = leaf(f, &[slot(81401)]);
        let x = Atom::var(symbolica::symbol!("safe_replace_cache_x"));
        let y = Atom::var(symbolica::symbol!("safe_replace_cache_y"));
        let args = symbolica::symbol!("safe_replace_cache_args__");
        let first = leaf(f, &[Atom::num(1), Atom::num(2), slot(81401)]);
        let second = leaf(f, &[Atom::num(3), Atom::num(4), slot(81401)]);
        for (mut source, pattern, expected_calls) in [
            (typed(&x * &fa + &y * &fa), fa.to_pattern(), [2, 1, 1]),
            (
                typed(&x * &first + &y * &first + &x * &second + &y * &second),
                leaf(f, &[Atom::var(args)]).to_pattern(),
                // Capacity one retains one repeated binding and leaves the
                // other uncached. Both keys contain Match::Multiple arguments.
                [4, 3, 2],
            ),
        ] {
            source.is_metric = true;
            source.is_composite = false;
            let mut condition_counts = Vec::new();
            for (cache, expected_calls) in [0, 1, 10].into_iter().zip(expected_calls) {
                let calls = Arc::new(AtomicUsize::new(0));
                let conditions = Arc::new(AtomicUsize::new(0));
                let observed = calls.clone();
                let rhs_atom = fa.clone();
                let rhs = ReplaceWith::Map(Box::new(move |matches| {
                    observed.fetch_add(1, Ordering::Relaxed);
                    if let Some(binding) = matches.get(args) {
                        let Match::Multiple(_, arguments) = binding else {
                            panic!("expected a ranged argument binding");
                        };
                        FunctionBuilder::new(f)
                            .add_args(arguments.iter().copied())
                            .finish()
                    } else {
                        rhs_atom.clone()
                    }
                }));
                let observed = conditions.clone();
                let condition = Condition::match_stack(move |_| {
                    observed.fetch_add(1, Ordering::Relaxed);
                    ConditionResult::True
                });
                let result = source
                    .replace(
                        &TensorRule::new(
                            pattern.clone(),
                            rhs.clone(),
                            Some(condition.clone()),
                            cache,
                        )
                        .unwrap(),
                    )
                    .unwrap();
                assert_eq!(result.expression, source.expression);
                assert_eq!(result.structure, source.structure);
                assert_eq!((result.is_metric, result.is_composite), (true, false));
                assert_eq!(calls.load(Ordering::Relaxed), expected_calls);
                condition_counts.push(conditions.load(Ordering::Relaxed));
            }
            assert!(condition_counts[0] >= 2);
            assert!(
                condition_counts
                    .iter()
                    .all(|&count| count == condition_counts[0])
            );
        }
    }

    #[test]
    fn tensor_replacement_fixed_head_keeps_match_and_rhs_callback_schedule() {
        use spenso::network::library::symbolic::ETS;
        use symbolica::id::MatchStack;

        let f = spenso::tensor_symbol!("safe_replace_fixed_f");
        let g = spenso::tensor_symbol!("safe_replace_fixed_g");
        let other = spenso::tensor_symbol!("safe_replace_fixed_other");
        let p = spenso::vector_symbol!("safe_replace_fixed_p");
        let q = spenso::vector_symbol!("safe_replace_fixed_q");
        let a = slot(81411);
        let fa = leaf(f, std::slice::from_ref(&a));
        let wildcard = symbolica::symbol!("safe_replace_fixed_slot_");
        let pattern = leaf(f, &[Atom::var(wildcard)]).to_pattern();
        assert!(matches!(pattern, Pattern::Fn(_, _)));
        let compact = ExtendibleReps::MINKOWSKI
            .new_rep(Dimension::Concrete(4))
            .to_symbolic([]);
        let dot = ETS.metric(
            leaf(p, std::slice::from_ref(&compact)),
            leaf(q, std::slice::from_ref(&compact)),
        );
        let x = Atom::var(symbolica::symbol!("safe_replace_fixed_x"));
        let y = Atom::var(symbolica::symbol!("safe_replace_fixed_y"));
        let source = typed((&x * &fa + &y * &fa + leaf(other, &[a])) * dot);
        for cache in [0, 10] {
            let mut outcomes = Vec::new();
            for certified in [false, true] {
                let events = Arc::new(Mutex::new(Vec::new()));
                let seen = events.clone();
                let condition = Condition::match_stack(move |matches| {
                    seen.lock().unwrap().push(format!("condition:{matches}"));
                    ConditionResult::True
                });
                let seen = events.clone();
                let rhs = move |matches: &MatchStack<'_>| {
                    seen.lock().unwrap().push(format!("rhs:{matches}"));
                    leaf(g, &[matches.get_atom(wildcard).unwrap().to_owned()])
                };
                let expression = if certified {
                    source
                        .replace(
                            &TensorRule::new(
                                pattern.clone(),
                                ReplaceWith::Map(Box::new(rhs)),
                                Some(condition.clone()),
                                cache,
                            )
                            .unwrap(),
                        )
                        .unwrap()
                        .expression
                } else {
                    source
                        .expression
                        .replace(&pattern)
                        .when(&condition)
                        .rhs_cache_size(cache)
                        .with_map(rhs)
                };
                let transcript = events.lock().unwrap().clone();
                assert_eq!(
                    transcript
                        .iter()
                        .filter(|event| event.starts_with("rhs:"))
                        .count(),
                    if cache == 0 { 2 } else { 1 }
                );
                outcomes.push((expression, transcript));
            }
            assert_eq!(outcomes[0], outcomes[1], "cache size {cache}");
        }
    }

    #[test]
    fn tensor_replacement_wildcard_heads_keep_tag_constraints() {
        let p = spenso::vector_symbol!("safe_replace_head_p");
        let q = spenso::vector_symbol!("safe_replace_head_q");
        let t = spenso::tensor_symbol!("safe_replace_head_t");
        let r = spenso::vector_symbol!("safe_replace_head_r");
        let head = symbolica::symbol!("safe_replace_head_");
        let argument = symbolica::symbol!("safe_replace_head_argument_");
        let a = slot(81421);
        let ta = leaf(t, std::slice::from_ref(&a));
        let source =
            typed(leaf(p, std::slice::from_ref(&a)) + leaf(q, std::slice::from_ref(&a)) + &ta);
        let pattern = Pattern::Fn(head, vec![Pattern::Wildcard(argument, false)]);
        let condition = head.filter_tag(SPENSO_TAG.rank1.clone());
        let rhs = leaf(r, &[Atom::var(argument)]).to_pattern();
        let expected = 2 * leaf(r, &[a]) + ta;
        assert_eq!(
            source
                .expression
                .replace(&pattern)
                .when(&condition)
                .with(&rhs),
            expected
        );
        assert_eq!(
            source
                .replace(
                    &TensorRule::new(
                        pattern.clone(),
                        ReplaceWith::Pattern((&rhs).into()),
                        Some(condition.clone()),
                        10
                    )
                    .unwrap()
                )
                .unwrap()
                .expression,
            expected
        );
    }

    #[test]
    fn tensor_replacement_non_function_patterns_keep_optional_matching() {
        let f = spenso::tensor_symbol!("safe_replace_optional_f");
        let g = spenso::tensor_symbol!("safe_replace_optional_g");
        let a = slot(81431);
        let fa = leaf(f, std::slice::from_ref(&a));
        let ga = leaf(g, &[a]);
        let source = typed(fa.clone());
        let x = Atom::var(symbolica::symbol!("safe_replace_optional_x"));
        let y = Atom::var(symbolica::symbol!("safe_replace_optional_y"));
        let repeated = typed(&x * &fa + &y * &fa);
        let argument = symbolica::symbol!("safe_replace_optional_argument_");
        let optional = symbolica::symbol!("safe_replace_optional_");
        let function = Pattern::Fn(f, vec![Pattern::Wildcard(argument, false)]);
        for pattern in [
            Pattern::Literal(fa),
            Pattern::Mul(vec![function.clone(), Pattern::Wildcard(optional, true)]),
            Pattern::Add(vec![function.clone(), Pattern::Wildcard(optional, true)]),
            Pattern::Pow(Box::new([function, Pattern::Wildcard(optional, true)])),
        ] {
            let expected = source
                .expression
                .replace(&pattern)
                .partial(false)
                .with(ga.clone());
            assert_eq!(expected, ga);
            assert_eq!(
                source
                    .replace(
                        &TensorRule::new(pattern.clone(), ga.clone().into(), None, 10).unwrap()
                    )
                    .unwrap()
                    .expression,
                expected
            );
            assert_eq!(
                repeated
                    .replace(&TensorRule::new(pattern.clone(), ga.clone().into(), None, 0).unwrap())
                    .unwrap()
                    .expression,
                &x * &ga + &y * &ga
            );
        }
    }

    #[test]
    fn tensor_replacement_reuses_matcher_without_retaining_failed_bindings() {
        use symbolica::id::MatchStack;

        let f = spenso::tensor_symbol!("safe_reuse_f");
        let g = spenso::tensor_symbol!("safe_reuse_g");
        let x = symbolica::symbol!("safe_reuse_x_");
        let y = symbolica::symbol!("safe_reuse_y_");
        let index = symbolica::symbol!("safe_reuse_index_");
        let function = |x, y| {
            Pattern::Fn(
                f,
                vec![
                    Pattern::Wildcard(x, false),
                    Pattern::Wildcard(y, false),
                    Pattern::Wildcard(index, false),
                ],
            )
        };
        let pattern = Pattern::Alternative(vec![function(x, y), function(y, x)]);
        let a = slot(81441);
        let first = leaf(f, &[Atom::num(1), Atom::num(2), a.clone()]);
        let rejected = leaf(f, &[Atom::num(2), Atom::num(2), a.clone()]);
        let last = leaf(f, &[Atom::num(3), Atom::num(4), a.clone()]);
        let short = leaf(f, &[Atom::num(5), a.clone()]);
        let long = leaf(f, &[Atom::num(6), Atom::num(7), Atom::num(8), a.clone()]);
        let scalar = Atom::var(symbolica::symbol!("safe_reuse_scalar"));
        let source = typed(&first + &scalar * &first + &rejected + last + &short + &long);
        let expected_first = leaf(g, &[Atom::num(2), Atom::num(1), a.clone()]);
        let expected = &expected_first
            + &scalar * &expected_first
            + rejected
            + leaf(g, &[Atom::num(4), Atom::num(3), a])
            + short
            + long;
        for cache in [0, 10] {
            let mut outcomes = Vec::new();
            for certified in [false, true] {
                let events = Arc::new(Mutex::new(Vec::new()));
                let completed = Arc::new(Mutex::new(HashMap::<String, usize>::new()));
                let seen = events.clone();
                let condition = Condition::match_stack(move |matches| {
                    let key = format!("{matches}");
                    seen.lock().unwrap().push(format!("condition:{key}"));
                    let (Some(left), Some(right), Some(_)) = (
                        matches.get_atom(x),
                        matches.get_atom(y),
                        matches.get_atom(index),
                    ) else {
                        return ConditionResult::Inconclusive;
                    };
                    let mut completed = completed.lock().unwrap();
                    let count = completed.entry(key).or_default();
                    *count += 1;
                    // The last insertion remains inconclusive, so acceptance
                    // requires next()'s final condition check and backtracking.
                    if *count % 2 == 1 {
                        ConditionResult::Inconclusive
                    } else if i64::try_from(left).unwrap() > i64::try_from(right).unwrap() {
                        ConditionResult::True
                    } else {
                        ConditionResult::False
                    }
                });
                let seen = events.clone();
                let rhs = move |matches: &MatchStack<'_>| {
                    seen.lock().unwrap().push(format!("rhs:{matches}"));
                    FunctionBuilder::new(g)
                        .add_arg(matches.get_atom(x).unwrap())
                        .add_arg(matches.get_atom(y).unwrap())
                        .add_arg(matches.get_atom(index).unwrap())
                        .finish()
                };
                let result = if certified {
                    source
                        .replace(
                            &TensorRule::new(
                                pattern.clone(),
                                ReplaceWith::Map(Box::new(rhs)),
                                Some(condition.clone()),
                                cache,
                            )
                            .unwrap(),
                        )
                        .unwrap()
                        .expression
                } else {
                    source
                        .expression
                        .replace(&pattern)
                        .when(&condition)
                        .rhs_cache_size(cache)
                        .with_map(rhs)
                };
                assert_eq!(result, expected);
                let transcript = events.lock().unwrap().clone();
                assert_eq!(
                    transcript
                        .iter()
                        .filter(|event| event.starts_with("rhs:"))
                        .count(),
                    if cache == 0 { 3 } else { 2 }
                );
                outcomes.push((result, transcript));
            }
            assert_eq!(outcomes[0], outcomes[1], "cache size {cache}");
        }
    }

    #[test]
    fn tensor_replacement_reuses_symmetric_argument_search() {
        use spenso::network::library::symbolic::ETS;

        // User symmetric heads are outside the source proof; the intrinsic
        // metric exercises the same unordered argument matcher on valid tensors.
        let f = ETS.metric;
        let g = spenso::tensor_symbol!("safe_reuse_symmetric_result");
        let x = symbolica::symbol!("safe_reuse_symmetric_x_");
        let [a, b, c, d, e] = [81451, 81452, 81453, 81454, 81455].map(slot);
        let missed = leaf(f, &[d, e]);
        let source =
            typed(leaf(f, &[a.clone(), b.clone()]) * leaf(f, &[a.clone(), c.clone()]) * &missed);
        let pattern = Pattern::Fn(
            f,
            vec![Pattern::Wildcard(x, false), Pattern::Literal(a.clone())],
        );
        let rhs = leaf(g, &[Atom::var(x), a.clone()]).to_pattern();
        let expected = leaf(g, &[b, a.clone()]) * leaf(g, &[c, a]) * missed;
        assert_eq!(source.expression.replace(&pattern).with(&rhs), expected);
        assert_eq!(
            source
                .replace(
                    &TensorRule::new(
                        pattern.clone(),
                        ReplaceWith::Pattern((&rhs).into()),
                        None,
                        0
                    )
                    .unwrap()
                )
                .unwrap()
                .expression,
            expected
        );
    }

    #[test]
    fn tensor_replacement_constructs_matcher_only_for_eligible_leaves() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        let f = spenso::tensor_symbol!("safe_reuse_lazy_f");
        let scalar = symbolica::symbol!("safe_reuse_lazy_scalar"; Scalar);
        let fa = leaf(f, &[slot(81461)]);
        let source = typed(fa.clone());
        let pattern = Pattern::Transformer(Box::new((None, Vec::new())));
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = calls.clone();
        let rhs = ReplaceWith::Map(Box::new(move |_| {
            seen.fetch_add(1, Ordering::Relaxed);
            Atom::Zero
        }));
        for source in [
            SymbolicTensor::from_produced_parts(
                Atom::one(),
                PartialStructure::from_logical_slots([]),
            )
            .unwrap(),
            SymbolicTensor::from_produced_parts(
                leaf(scalar, std::slice::from_ref(&fa)),
                PartialStructure::from_logical_slots([]),
            )
            .unwrap(),
            SymbolicTensor::new(Atom::Zero, source.structure.clone()).unwrap(),
        ] {
            let result = source
                .replace(&TensorRule::new(pattern.clone(), rhs.clone(), None, 0).unwrap())
                .unwrap();
            assert_eq!(result.expression, source.expression);
            assert_eq!(result.structure, source.structure);
        }
        // The existing matcher rejects transformers on the LHS once reached.
        assert!(
            catch_unwind(AssertUnwindSafe(|| {
                source.replace(&TensorRule::new(pattern.clone(), rhs.clone(), None, 0).unwrap())
            }))
            .is_err()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn tensor_replacement_cached_rhs_is_checked_against_each_actual_target() {
        let p = spenso::tensor_symbol!("safe_replace_alternative_p");
        let q = spenso::tensor_symbol!("safe_replace_alternative_q");
        let r = spenso::tensor_symbol!("safe_replace_alternative_r");
        let a = slot(81501);
        let b = slot(81503);
        let wildcard = Atom::var(symbolica::symbol!("safe_replace_alternative_slot_"));
        let pattern = leaf(p, std::slice::from_ref(&wildcard))
            .to_pattern()
            .alternative(leaf(q, &[wildcard, b.clone()]).to_pattern())
            .unwrap();
        let source = typed(leaf(p, std::slice::from_ref(&a)) * leaf(q, &[a.clone(), b]));
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = calls.clone();
        let rhs_atom = leaf(r, &[a]);
        let rhs = ReplaceWith::Map(Box::new(move |_| {
            seen.fetch_add(1, Ordering::Relaxed);
            rhs_atom.clone()
        }));
        assert!(
            source
                .replace(&TensorRule::new(pattern.clone(), rhs.clone(), None, 10).unwrap())
                .is_err()
        );
        // Both patterns produce the same bindings, but the second target owns
        // an additional port. Caching the first compatibility result is invalid.
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn tensor_replacement_rejects_callback_sources_before_rhs_and_callback_rank_loss() {
        let f = spenso::tensor_symbol!("safe_replace_callback_f");
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = calls.clone();
        let callback = symbolica::symbol!("safe_replace_hidden_callback"; Scalar; norm=move|_,_| { seen.fetch_add(1,Ordering::Relaxed); });
        let a = slot(81601);
        let raw = leaf(f, &[leaf(callback, &[Atom::num(3)]), a.clone()]);
        let source = typed(raw.clone());
        calls.store(0, Ordering::Relaxed);
        let rhs_calls = Arc::new(AtomicUsize::new(0));
        let seen = rhs_calls.clone();
        let rhs = ReplaceWith::Map(Box::new(move |_| {
            seen.fetch_add(1, Ordering::Relaxed);
            Atom::Zero
        }));
        let absent = spenso::tensor_symbol!("safe_replace_callback_absent");
        let absent_pattern = leaf(
            absent,
            &[Atom::var(symbolica::symbol!(
                "safe_replace_callback_argument_"
            ))],
        )
        .to_pattern();
        assert!(matches!(absent_pattern, Pattern::Fn(_, _)));
        for pattern in [raw.to_pattern(), absent_pattern] {
            assert!(
                source
                    .replace(&TensorRule::new(pattern.clone(), rhs.clone(), None, 0).unwrap())
                    .is_err()
            );
        }
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        assert_eq!(rhs_calls.load(Ordering::Relaxed), 0);
        let ordinary = leaf(f, &[a]);
        let source = typed(ordinary.clone());
        assert!(
            source
                .replace(
                    &TensorRule::new(ordinary.to_pattern(), Atom::one().into(), None, 0).unwrap()
                )
                .is_err()
        );
        let placeholder = Atom::var(SPENSO_TAG.chain_in);
        assert!(
            source
                .replace(
                    &TensorRule::new(ordinary.to_pattern(), placeholder.into(), None, 0).unwrap()
                )
                .is_err()
        );
    }

    #[test]
    fn tensor_replacement_rejects_bare_tagged_names_even_in_scalar_rhs() {
        let f = spenso::tensor_symbol!("safe_replace_bare_f");
        let u = spenso::tensor_symbol!("safe_replace_bare_u");
        let a = slot(81651);
        let source = typed(leaf(f, &[a.clone(), a]));
        let pattern = source.expression.to_pattern();
        for rhs in [Atom::var(u), Atom::var(u) + Atom::one()] {
            let error = source
                .replace(&TensorRule::new(pattern.clone(), rhs.into(), None, 0).unwrap())
                .unwrap_err();
            assert!(error.to_string().contains("must be called"), "{error}");
        }
    }

    #[test]
    fn tensor_replacement_keeps_scalar_metadata_opaque_and_declines_unresolved_ports() {
        let f = spenso::tensor_symbol!("safe_replace_opaque_f");
        let fa = leaf(f, &[slot(81701)]);
        let scalar_head = symbolica::symbol!("safe_replace_opaque_scalar"; Scalar);
        let source = SymbolicTensor::from_produced_parts(
            leaf(scalar_head, std::slice::from_ref(&fa)),
            PartialStructure::from_logical_slots([]),
        )
        .unwrap();
        assert_eq!(
            source
                .replace(&TensorRule::new(fa.to_pattern(), Atom::Zero.into(), None, 0).unwrap())
                .unwrap()
                .expression,
            source.expression
        );
        let m = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let unresolved = SymbolicTensor::new(
            leaf(f, &[m.to_symbolic([])]),
            PartialStructure::from_logical_slots([m.slot(PartialIndex::open(0))]),
        )
        .unwrap();
        assert!(
            unresolved
                .replace(
                    &TensorRule::new(
                        unresolved.expression.to_pattern(),
                        Atom::Zero.into(),
                        None,
                        0
                    )
                    .unwrap()
                )
                .is_err()
        );
    }

    #[test]
    fn tensor_replacement_matches_normalized_rounded_container_reconstruction() {
        let f = spenso::tensor_symbol!("safe_replace_rounded_f");
        let g = spenso::tensor_symbol!("safe_replace_rounded_g");
        let a = slot(81801);
        let fa = leaf(f, std::slice::from_ref(&a));
        let ga = leaf(g, &[a]);
        let rounded = |s| Atom::num(symbolica::domains::float::Float::parse(s, Some(53)).unwrap());
        let x = Atom::var(symbolica::symbol!("safe_replace_rounded_x"));
        let y = Atom::var(symbolica::symbol!("safe_replace_rounded_y"));
        let rhs = rounded("0.2") * &ga;
        for expression in [
            rounded("0.3") * (rounded("0.1") * &x) * &fa,
            rounded("0.3") * &fa + rounded("0.1") * &ga,
            (&x + &y) * &fa * &ga,
        ] {
            let source = typed(expression.clone());
            let expected = expression.replace(fa.to_pattern()).with(rhs.clone());
            let actual = source
                .replace(&TensorRule::new(fa.to_pattern(), rhs.clone().into(), None, 0).unwrap())
                .unwrap();
            assert_eq!(actual.expression, expected);
        }
    }

    #[test]
    fn carried_rule_proof_survives_identity_and_is_cleared_by_mutable_structure() {
        let f = tensor_symbol!("carried_rule_source");
        let g = tensor_symbol!("carried_rule_result");
        let a = slot(81941);
        let inferred =
            SymbolicTensor::<PartialStructure>::infer(leaf(f, slice::from_ref(&a))).unwrap();
        assert_eq!(inferred.proofs.validated.get(), Some(&true));
        let (expression, structure) = inferred.into_parts();
        // Detached storage must establish its own facts even when its parts
        // came from a previously checked value.
        let source = SymbolicTensor::from_declared_parts(expression, structure).unwrap();
        let rule = TensorRule::new(
            source.expression.to_pattern(),
            leaf(g, &[a]).into(),
            None,
            0,
        )
        .unwrap();
        SCOPE_VALIDATIONS.with(|count| count.set(0));
        assert!(source.rewrites_preserve_interface());
        let established = SCOPE_VALIDATIONS.with(Cell::get);
        assert_eq!(established, 0);
        assert!(source.rewrites_preserve_interface());
        let same = source
            .with_checked_expression(source.expression.clone())
            .unwrap();
        assert!(same.rewrites_preserve_interface());
        assert_eq!(SCOPE_VALIDATIONS.with(Cell::get), established);

        let changed = source.replace(&rule).unwrap();
        assert_eq!(changed.proofs.validated.get(), Some(&true));
        assert!(changed.proofs.rewrite.get().is_none());
        assert!(changed.proofs.algebra.get().is_none());
        assert!(changed.proofs.intrinsic.get().is_none());
        assert_eq!(
            changed.proofs.normalized_product.get(),
            source.proofs.normalized_product.get()
        );
        let after_binding_checks = SCOPE_VALIDATIONS.with(Cell::get);
        assert!(changed.rewrites_preserve_interface());
        assert_eq!(SCOPE_VALIDATIONS.with(Cell::get), after_binding_checks);
        // Mutable structural access exists on the ordered network carrier; the
        // partial-interface carrier exposes only checked transformations.
        let mut ordered = SymbolicTensor {
            structure: changed.structure.canonical().clone(),
            expression: changed.expression,
            is_metric: changed.is_metric,
            is_composite: changed.is_composite,
            proofs: changed.proofs,
        };
        let mapped = ordered.clone().map_same_structure(|value| value);
        assert!(mapped.proofs.rewrite.get().is_none());
        // Even obtaining mutable access without an edit clears established facts.
        let _ = ordered.mut_structure();
        assert!(ordered.proofs.rewrite.get().is_none());
        assert!(ordered.proofs.validated.get().is_none());
        assert!(ordered.proofs.algebra.get().is_none());
        *ordered.mut_structure() = OrderedStructure::new(Vec::new()).into_canonical();
        let (expression, structure) = ordered.into_parts();
        assert!(
            SymbolicTensor::from_declared_parts(
                expression,
                PartialStructure::from_logical_slots(structure.external_structure()),
            )
            .is_err()
        );

        let (expression, _) = source.into_parts();
        assert!(
            SymbolicTensor::from_declared_parts(
                expression,
                PartialStructure::from_logical_slots([]),
            )
            .is_err()
        );
    }
}
