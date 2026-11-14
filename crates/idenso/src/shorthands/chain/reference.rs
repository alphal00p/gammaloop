//! The chainify implementation before its matcher pre-checks, kept verbatim
//! as the reference that every test's chainify call is compared against.

use std::cell::{Cell, RefCell};

use spenso::{
    chain,
    network::{library::symbolic::ETS, parsing::ParseState, tags::SPENSO_TAG as T},
    structure::{
        abstract_index::AbstractIndex,
        representation::LibraryRep,
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol};

use super::{chainified_slot_function, chainify_endpoints, chainify_rule};
use crate::W_;

mod tests;

/// What a chainify call runs on this thread.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mode {
    /// The pre-checked implementation, compared with the reference.
    Checked,
    /// Only the pre-checked implementation, for timing.
    PreChecked,
    /// Only the reference implementation, for timing.
    Reference,
}

/// Matcher work of the pre-checked chainify: `[skipped, run]` for the
/// endpoint rule at a function and for an inner product's selection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct MatcherCounts {
    pub(crate) rule: [usize; 2],
    pub(crate) selection: [usize; 2],
}

thread_local! {
    static MODE: Cell<Mode> = const { Cell::new(Mode::Checked) };
    static INSIDE: Cell<bool> = const { Cell::new(false) };
    static CROSS_CHECKS: Cell<usize> = const { Cell::new(0) };
    static COUNTS: Cell<MatcherCounts> = const {
        Cell::new(MatcherCounts { rule: [0; 2], selection: [0; 2] })
    };
    static RECORDED: RefCell<Option<Vec<(Atom, LibraryRep)>>> = const { RefCell::new(None) };
}

pub(super) fn count(update: impl FnOnce(&mut MatcherCounts)) {
    COUNTS.with(|counts| {
        let mut current = counts.get();
        update(&mut current);
        counts.set(current);
    });
}

/// Matcher work since the previous call on this thread.
pub(crate) fn take_counts() -> MatcherCounts {
    COUNTS.with(|counts| counts.replace(MatcherCounts::default()))
}

/// Number of chainify results compared with the reference on this thread.
pub(crate) fn cross_checks() -> usize {
    CROSS_CHECKS.with(Cell::get)
}

/// Run `body` with chainify in `mode` on this thread.
pub(crate) fn with_mode<R>(mode: Mode, body: impl FnOnce() -> R) -> R {
    let previous = MODE.with(|current| current.replace(mode));
    let result = body();
    MODE.with(|current| current.set(previous));
    result
}

/// Run `body` and return every chainify input it reaches on this thread.
pub(crate) fn record<R>(body: impl FnOnce() -> R) -> (R, Vec<(Atom, LibraryRep)>) {
    let previous = RECORDED.with(|recorded| recorded.replace(Some(Vec::new())));
    let result = body();
    let inputs = RECORDED
        .with(|recorded| recorded.replace(previous))
        .expect("recording was enabled");
    (result, inputs)
}

/// Called on entry to chainify. `None` lets the pre-checked body run; it
/// always does when re-entered from here. Otherwise this returns the result
/// for this thread's mode, comparing both implementations byte for byte
/// when checked.
pub(super) fn intercept(expression: AtomView<'_>, representation: LibraryRep) -> Option<Atom> {
    if INSIDE.with(Cell::get) {
        return None;
    }
    RECORDED.with(|recorded| {
        if let Some(inputs) = recorded.borrow_mut().as_mut() {
            inputs.push((expression.to_owned(), representation));
        }
    });
    let pre_checked = || {
        INSIDE.with(|inside| inside.set(true));
        let chained = super::Chain::chainify(&expression, representation);
        INSIDE.with(|inside| inside.set(false));
        chained
    };
    match MODE.with(Cell::get) {
        Mode::PreChecked => Some(pre_checked()),
        Mode::Reference => Some(chainify(&expression, representation)),
        Mode::Checked => {
            let chained = pre_checked();
            let expected = chainify(&expression, representation);
            assert!(
                chained.as_view().get_data() == expected.as_view().get_data(),
                "chainify differs from its reference for {representation:?} on\n{expression}\n\
                 pre-checked: {chained}\nreference:   {expected}"
            );
            CROSS_CHECKS.with(|checks| checks.set(checks.get() + 1));
            Some(chained)
        }
    }
}

pub(super) fn chainify(this: &AtomView<'_>, representation: LibraryRep) -> Atom {
    // Both endpoints must use this representation. Terminal metric/epsilon
    // expressions often have no such slots; avoid building a replacement
    // and invoking its matcher at every unrelated function in that case.
    if !this.contains_symbol(representation.symbol()) {
        return this.to_owned();
    }

    let in_index = Atom::var(W_.i_);
    let out_index = Atom::var(W_.j_);
    let endpoints = chainify_endpoints(representation);
    let replacement = chainify_rule(representation);
    let opaque_heads = [T.chain, T.trace];
    let existing_chain = chain!(&in_index, &out_index, W_.a___).to_pattern();
    // An unchanged pruning assignment in replace_map rebuilds ancestors
    // and invokes their normalizers. Keep borrowed children until this
    // existing chain rule actually changes a node.
    fn rewrite<'a>(
        atom: AtomView<'a>,
        slots: &mut SlotMatcher,
        opaque_heads: &[Symbol],
        apply: &impl Fn(AtomView<'_>) -> Atom,
        open: &impl Fn(symbolica::atom::representation::FunView<'_>) -> Option<Atom>,
    ) -> AtomOrView<'a> {
        if !matches!(slots.classify(atom), SlotMatch::Other) {
            return atom.into();
        }
        let children = match atom {
            AtomView::Fun(function) => {
                if let Some(opened) = open(function) {
                    return rewrite(opened.as_view(), slots, opaque_heads, apply, open)
                        .into_owned()
                        .into();
                }
                if opaque_heads.contains(&function.get_symbol()) {
                    return atom.into();
                }
                let changed = apply(atom);
                if changed.as_view() != atom {
                    return changed.into();
                }
                function.iter().collect::<Vec<_>>()
            }
            AtomView::Add(sum) => sum.iter().collect(),
            AtomView::Mul(product) => product.iter().collect(),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                vec![base, exponent]
            }
            _ => return atom.into(),
        };
        let children = children
            .into_iter()
            .map(|child| rewrite(child, slots, opaque_heads, apply, open))
            .collect::<Vec<_>>();
        if children
            .iter()
            .all(|child| matches!(child, AtomOrView::View(_)))
        {
            return atom.into();
        }
        match atom {
            AtomView::Fun(function) => {
                let mut output = FunctionBuilder::new(function.get_symbol());
                for child in children {
                    output = output.add_arg(child);
                }
                output.finish().into()
            }
            AtomView::Add(_) => Atom::add_many(children).into(),
            AtomView::Mul(_) => Atom::mul_many(children).into(),
            AtomView::Pow(_) => children[0].as_view().pow(children[1].as_view()).into(),
            _ => unreachable!("only compound expressions have children"),
        }
    }
    let state = std::cell::OnceCell::new();
    rewrite(
        *this,
        &mut SlotMatcher::default(),
        &opaque_heads,
        &|atom| match chainified_slot_function(atom, representation) {
            Some(chained) => chained.unwrap_or_else(|| atom.to_owned()),
            None => atom.replace_multiple([&replacement]),
        },
        &|function| {
            if function.get_symbol() != T.dot && function.get_symbol() != ETS.metric {
                return None;
            }
            let (left, right, _, _) = spenso::structure::slot::SlotMatcher::default()
                .compact_inner_product_parts::<AbstractIndex>(function)?;
            // Select with the existing endpoint rule, including its opaque-head
            // and placeholder restrictions. Existing chains already own their
            // endpoints.
            // Either operand may carry the selected endpoints, as in gamma·p.
            let selected = [left, right].into_iter().any(|operand| {
                operand
                    .pattern_match(&replacement.pat, replacement.conditions.as_ref(), None)
                    .next()
                    .is_some()
                    || operand
                        .pattern_match(&existing_chain, Some(&endpoints), None)
                        .next()
                        .is_some()
            });
            if !selected {
                return None;
            }
            let state = state.get_or_init(|| {
                let state = ParseState::<AbstractIndex>::default();
                state.reserve_indices(*this);
                state
            });
            state.materialize_inner_product(function)
        },
    )
    .into_owned()
}
