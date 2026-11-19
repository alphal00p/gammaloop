//! Frontier-level closing of quark loops (opt-in,
//! [`crate::tensor::AlgebraSettings::gamma_closing`]; inert when unset).
//!
//! Without closing, the GAMMA pass evaluates a quark loop as an OPEN trace:
//! every Lorentz port that leaves the loop (a gluon line) stays a free index,
//! the trace expands into metrics and momenta on those ports, and the
//! factorized frontier contracts that expansion with the vertex sums and
//! propagators of the gluon network. Closing reverses the order:
//!
//! 1. the GAMMA arm writes a single-cycle quark loop whose environment holds a
//!    tensor sum as ONE trace word (the gamma cycle in matrix order, its slash
//!    vectors -- routed sums included -- as letters, the gluon ports open),
//!    and does not evaluate it while a sum still reaches its ports;
//! 2. the factorized frontier never selects a word: the word rides along with
//!    its bound ports while the environment merges into states;
//! 3. at emission the word of every state, closed by the state's bound
//!    vectors, is evaluated once per distinct binding by the gamma kernel
//!    (routed sums as opaque stand-in vectors) and folded into the state's
//!    polynomial; the closed states are summed as one polynomial, grouped by
//!    their values; the stand-ins are spelled as their sums once per variable.
//!
//! Every step is an identity; values equal the unclosed path's. A word still
//! open at emission (a port on the root) is emitted as a word and evaluated by
//! the next GAMMA round. Symbolica's cooperative resource budget is checked
//! before every word evaluation; an exhausted budget aborts as anywhere in the
//! frontier, and the thread's closing switch is restored on unwinding.

use std::collections::HashMap;

use spenso::structure::partial::PartialStructure;
use spenso::structure::representation::{LibraryRep, Minkowski};
use spenso::structure::slot::{SlotMatch, SlotMatcher};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use crate::representations::Bispinor;
use crate::tensor::SymbolicTensor;
use crate::tensor::inference::TensorInferenceError;

/// Frontier-level closing of quark loops ([`crate::tensor::AlgebraSettings::gamma_closing`]).
/// No parameters today; constructed with `GammaClosing::default()`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct GammaClosing {}

fn budget_checkpoint() {
    if let Err(error) = symbolica::utils::resources::checkpoint() {
        panic!("resource budget exhausted while writing a Dirac word: {error}");
    }
}

struct Symbols {
    minkowski: Symbol,
    bispinor: Symbol,
}

fn symbols() -> Symbols {
    Symbols {
        minkowski: LibraryRep::from(Minkowski {}).symbol(),
        bispinor: LibraryRep::from(Bispinor {}).symbol(),
    }
}

/// Free slots of an expression (multiset with the contracted pairs of a
/// product removed): Minkowski slots, and whether a bispinor slot occurs.
#[derive(Default)]
struct Slots {
    minkowski: Vec<Atom>,
    bispinor: bool,
    /// An explicit slot of another representation (colour, ...).
    other: bool,
}

fn free_slots(value: AtomView<'_>, matcher: &mut SlotMatcher, symbols: &Symbols) -> Slots {
    let mut out = Slots::default();
    match value {
        AtomView::Add(add) => {
            // every term carries the sum's interface: read the first one, and
            // a bispinor slot anywhere marks the factor as Dirac
            let mut first = true;
            for term in add.iter() {
                let slots = free_slots(term, matcher, symbols);
                out.bispinor |= slots.bispinor;
                out.other |= slots.other;
                if first {
                    out.minkowski = slots.minkowski;
                    first = false;
                }
            }
        }
        AtomView::Mul(mul) => {
            let mut all = Vec::new();
            for factor in mul.iter() {
                let slots = free_slots(factor, matcher, symbols);
                out.bispinor |= slots.bispinor;
                out.other |= slots.other;
                all.extend(slots.minkowski);
            }
            out.minkowski = once_only(all);
        }
        AtomView::Pow(pow) => {
            let (base, _) = pow.get_base_exp();
            let slots = free_slots(base, matcher, symbols);
            out.bispinor = slots.bispinor;
            out.other = slots.other;
            // a tensor power is never a closing candidate: keep its slots so
            // that the term is seen as connected, and decline below
            out.minkowski = slots.minkowski;
        }
        // any bispinor representation (an explicit slot, or the stripped one
        // of a collected trace / chain) marks a Dirac factor
        AtomView::Fun(function) if function.get_symbol() == symbols.bispinor => {
            out.bispinor = true;
        }
        AtomView::Fun(function) => match matcher.classify(value) {
            SlotMatch::Explicit(slot) => {
                let head = slot.representation().head();
                if head == symbols.minkowski {
                    out.minkowski.push(value.to_owned());
                } else if head == symbols.bispinor {
                    out.bispinor = true;
                } else {
                    out.other = true;
                }
            }
            SlotMatch::Opaque => {}
            SlotMatch::Other => {
                let mut all = Vec::new();
                for argument in function.iter() {
                    let slots = free_slots(argument, matcher, symbols);
                    out.bispinor |= slots.bispinor;
                    out.other |= slots.other;
                    all.extend(slots.minkowski);
                }
                out.minkowski = once_only(all);
            }
        },
        _ => {}
    }
    out
}

/// The slots of a multiset that occur exactly once.
fn once_only(mut slots: Vec<Atom>) -> Vec<Atom> {
    if slots.len() < 2 {
        return slots;
    }
    slots.sort();
    let mut out = Vec::with_capacity(slots.len());
    let mut i = 0;
    while i < slots.len() {
        let mut j = i + 1;
        while j < slots.len() && slots[j] == slots[i] {
            j += 1;
        }
        if j - i == 1 {
            out.push(slots[i].clone());
        }
        i = j;
    }
    out
}

/// The Dirac region of a product as ONE cycle of `gamma(bis_in, bis_out,
/// letter)` factors in matrix order, its slash vectors substituted into the
/// letters (stripped Lorentz representation), ports kept as slots: the
/// bispinor representation (stripped) and the letters; `None` unless every
/// bispinor factor is such a gamma and they form a single cycle.
fn word_letters(
    factors: &[AtomView<'_>],
    slots: &[Slots],
    slash: &[bool],
    symbols: &Symbols,
) -> Option<(Atom, Vec<Atom>)> {
    let gamma_symbol = super::AGS.gamma;
    let mut entries: Vec<(Atom, Atom, Atom)> = Vec::new();
    for (i, factor) in factors.iter().enumerate() {
        if !slots[i].bispinor {
            continue;
        }
        let AtomView::Fun(function) = factor else {
            return None;
        };
        if function.get_symbol() != gamma_symbol || function.get_nargs() != 3 {
            return None;
        }
        let arguments: Vec<AtomView<'_>> = function.iter().collect();
        for spin in &arguments[..2] {
            let AtomView::Fun(slot) = spin else {
                return None;
            };
            if slot.get_symbol() != symbols.bispinor || slot.get_nargs() != 2 {
                return None;
            }
        }
        entries.push((
            arguments[0].to_owned(),
            arguments[1].to_owned(),
            arguments[2].to_owned(),
        ));
    }
    if entries.len() < 2 {
        return None;
    }
    let by_input: HashMap<Atom, usize> = entries
        .iter()
        .enumerate()
        .map(|(i, entry)| (entry.0.clone(), i))
        .collect();
    let mut order = vec![0usize];
    let mut current = 0usize;
    loop {
        let next = *by_input.get(&entries[current].1)?;
        if next == 0 {
            break;
        }
        if order.len() >= entries.len() {
            return None;
        }
        order.push(next);
        current = next;
    }
    if order.len() != entries.len() {
        return None;
    }
    let representation = stripped_slot(&entries[0].0)?;
    let mut letters = Vec::with_capacity(order.len());
    for &i in &order {
        let letter = entries[i].2.clone();
        let vector = (0..factors.len())
            .find(|&j| slash[j] && slots[j].minkowski[0] == letter)
            .map(|j| factors[j].to_owned());
        match vector {
            Some(vector) => {
                letters.push(vector.replace(letter.clone()).with(stripped_slot(&letter)?))
            }
            None => letters.push(letter),
        }
    }
    Some((representation, letters))
}

/// `rep(dim, index)` -> `rep(dim)`.
fn stripped_slot(slot: &Atom) -> Option<Atom> {
    let AtomView::Fun(function) = slot.as_view() else {
        return None;
    };
    let dimension = function.iter().next()?;
    Some(
        symbolica::atom::FunctionBuilder::new(function.get_symbol())
            .add_arg(dimension)
            .finish(),
    )
}

// ---- the closing switch -------------------------------------------------------

thread_local! {
    /// Whether frontier closing is active on this thread.
    static FRONTIER: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether frontier closing is active on this thread (set for the duration of
/// a `simplify_algebra` whose settings ask for it; read by the frontier).
#[inline]
pub(crate) fn frontier_active() -> bool {
    FRONTIER.with(std::cell::Cell::get)
}

/// Run `operation` with frontier closing active (restored afterwards, also on
/// unwinding).
pub(crate) fn with_frontier<R>(operation: impl FnOnce() -> R) -> R {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            FRONTIER.with(|cell| cell.set(self.0));
        }
    }
    let _restore = Restore(FRONTIER.with(|cell| cell.replace(true)));
    operation()
}

/// A Dirac trace word: `trace(rep, gamma(in, out, letter), ...)` whose every
/// factor is a gamma (or a sum of scaled gammas: a slash of a routed sum).
/// Colour traces, gamma5 words and other traces are never words.
pub(crate) fn is_word(value: AtomView<'_>) -> bool {
    let AtomView::Fun(function) = value else {
        return false;
    };
    let Some((_, factors)) = spenso::shadowing::trace_parts(function) else {
        return false;
    };
    let gamma_symbol = super::AGS.gamma;
    let is_gamma =
        |view: AtomView<'_>| matches!(view, AtomView::Fun(f) if f.get_symbol() == gamma_symbol);
    !factors.is_empty()
        && factors.iter().all(|factor| match factor {
            AtomView::Fun(_) => is_gamma(*factor),
            AtomView::Add(sum) => sum.iter().all(|term| match term {
                AtomView::Mul(product) => product.iter().filter(|f| is_gamma(*f)).count() == 1,
                other => is_gamma(other),
            }),
            _ => false,
        })
}

/// What the GAMMA arm does in frontier mode.
pub(crate) enum FrontierStep {
    /// A Dirac region was written as a word: the new expression.
    Written(Box<SymbolicTensor<PartialStructure>>),
    /// A word still meets a tensor sum on its own ports: leave it to the
    /// frontier (only while STRUCTURE is still to run).
    Pending,
    /// Nothing for frontier closing: evaluate as usual.
    Nothing,
}

/// The GAMMA arm's step in frontier mode: [`FrontierStep::Written`] when it
/// rewrote at least one Dirac region into a word; [`FrontierStep::Pending`]
/// when a word already written still meets a tensor sum through its own
/// Lorentz ports (the arm then runs the structural contraction, which closes
/// it, instead of evaluating it); [`FrontierStep::Nothing`] lets the arm
/// evaluate as usual (a word whose ports reach no sum included).
pub(crate) fn frontier_step(
    tensor: &SymbolicTensor<PartialStructure>,
) -> Result<FrontierStep, TensorInferenceError> {
    let symbols = symbols();
    let mut matcher = SlotMatcher::default();
    let expression = tensor.expression.as_view();
    let terms: Vec<AtomView<'_>> = match expression {
        AtomView::Add(add) => add.iter().collect(),
        _ => vec![expression],
    };
    let mut changed = false;
    let mut pending = false;
    let mut out = Atom::Zero;
    for term in terms {
        budget_checkpoint();
        let AtomView::Mul(mul) = term else {
            out += term.to_owned();
            continue;
        };
        let factors: Vec<AtomView<'_>> = mul.iter().collect();
        let slots: Vec<Slots> = factors
            .iter()
            .map(|factor| free_slots(*factor, &mut matcher, &symbols))
            .collect();
        // a tensor sum (a vertex sum, a projector) somewhere in the product
        let tensor_sum = |skip: &[bool]| {
            (0..factors.len()).any(|i| {
                !skip[i]
                    && !slots[i].bispinor
                    && matches!(factors[i], AtomView::Add(_))
                    && !slots[i].minkowski.is_empty()
            })
        };
        if let Some(word) = factors.iter().position(|f| is_word(*f)) {
            // already a word: wait for the frontier only while a tensor sum is
            // connected to the word's own ports (through Lorentz slots)
            let mut reached: std::collections::HashSet<Atom> =
                slots[word].minkowski.iter().cloned().collect();
            let mut frontier: Vec<Atom> = reached.iter().cloned().collect();
            let mut seen = vec![false; factors.len()];
            seen[word] = true;
            let mut meets_sum = false;
            while let Some(slot) = frontier.pop() {
                for (i, s) in slots.iter().enumerate() {
                    if seen[i] || s.bispinor || !s.minkowski.contains(&slot) {
                        continue;
                    }
                    seen[i] = true;
                    meets_sum |= matches!(factors[i], AtomView::Add(_));
                    for other in &s.minkowski {
                        if reached.insert(other.clone()) {
                            frontier.push(other.clone());
                        }
                    }
                }
            }
            pending |= meets_sum;
            out += term.to_owned();
            continue;
        }
        let dirac: Vec<usize> = (0..factors.len()).filter(|&i| slots[i].bispinor).collect();
        let ports = once_only(
            dirac
                .iter()
                .flat_map(|&i| slots[i].minkowski.iter().cloned())
                .collect(),
        );
        let mut slash = vec![false; factors.len()];
        for (i, s) in slots.iter().enumerate() {
            // a slash: ONE Lorentz slot, no other explicit slot (a colour
            // index stays outside the word)
            if !s.bispinor
                && !s.other
                && s.minkowski.len() == 1
                && matches!(factors[i], AtomView::Fun(_) | AtomView::Add(_))
                && ports.contains(&s.minkowski[0])
            {
                slash[i] = true;
            }
        }
        // the word only pays when the environment holds a tensor sum
        if dirac.is_empty() || !tensor_sum(&slash) {
            out += term.to_owned();
            continue;
        }
        let Some((representation, letters)) = word_letters(&factors, &slots, &slash, &symbols)
        else {
            out += term.to_owned();
            continue;
        };
        let word =
            spenso::trace!(&representation; letters.iter().map(|letter| crate::gamma!(letter)));
        let mut rest = Atom::num(1);
        for (i, factor) in factors.iter().enumerate() {
            if !slots[i].bispinor && !slash[i] {
                rest = &rest * &factor.to_owned();
            }
        }
        #[cfg(test)]
        tests::WORDS.with(|count| count.set(count.get() + 1));
        out += &rest * &word;
        changed = true;
        pending = true;
    }
    if changed {
        return Ok(FrontierStep::Written(Box::new(
            tensor.with_identity_result(out, None)?,
        )));
    }
    Ok(if pending {
        FrontierStep::Pending
    } else {
        FrontierStep::Nothing
    })
}

/// Evaluate a CLOSED word (every Lorentz slot paired or slashed) with the
/// gamma kernel. A sum of momenta among the letters (a routed slash) enters
/// as an opaque stand-in vector, named after the sum (the same sum gets the
/// same stand-in everywhere), and is NOT substituted back: the value and the
/// stand-ins with their sums are returned, so that a caller can keep the
/// value compact over the stand-ins and substitute once at the end (each
/// stand-in then distributes inside its dots only). `None` when the word is
/// not a closed gamma word.
pub(crate) fn evaluate_closed_word(word: AtomView<'_>) -> Option<(Atom, Vec<(Atom, Atom)>)> {
    let symbols = symbols();
    let AtomView::Fun(function) = word else {
        return None;
    };
    if !is_word(word) {
        return None;
    }
    // trace(rep, cyclic(gamma(..), ..)): the representation and the letters
    let (representation, arguments) = spenso::shadowing::trace_parts(function)?;
    let representation = representation.to_owned();
    let gamma_symbol = super::AGS.gamma;
    let mut letters: Vec<Atom> = Vec::new();
    let mut explicit: Vec<Atom> = Vec::new();
    let mut vector_representation: Option<Atom> = None;
    // a gamma letter's vector; a slash of a SUM of momenta is written as a sum
    // of gammas (gamma is linear in its vector): read it back as one letter
    let letter_of = |argument: AtomView<'_>| -> Option<Atom> {
        let vector = |gamma: AtomView<'_>| -> Option<Atom> {
            let AtomView::Fun(gamma) = gamma else {
                return None;
            };
            (gamma.get_symbol() == gamma_symbol && gamma.get_nargs() == 3)
                .then(|| gamma.iter().nth(2).map(|v| v.to_owned()))
                .flatten()
        };
        match argument {
            AtomView::Fun(_) => vector(argument),
            AtomView::Add(sum) => {
                let mut total = Atom::Zero;
                for term in sum.iter() {
                    let (coefficient, gamma) = match term {
                        AtomView::Mul(product) => {
                            let mut coefficient = Atom::num(1);
                            let mut gamma = None;
                            for factor in product.iter() {
                                if matches!(factor, AtomView::Fun(f) if f.get_symbol() == gamma_symbol)
                                {
                                    if gamma.is_some() {
                                        return None;
                                    }
                                    gamma = Some(factor);
                                } else {
                                    coefficient = &coefficient * &factor.to_owned();
                                }
                            }
                            (coefficient, gamma?)
                        }
                        _ => (Atom::num(1), term),
                    };
                    total += &coefficient * &vector(gamma)?;
                }
                Some(total)
            }
            _ => None,
        }
    };
    for argument in arguments {
        let letter_atom = letter_of(argument)?;
        let letter = letter_atom.as_view();
        if let AtomView::Fun(slot) = letter
            && slot.get_symbol() == symbols.minkowski
        {
            if slot.get_nargs() == 2 {
                explicit.push(letter.to_owned());
                vector_representation =
                    vector_representation.or_else(|| stripped_slot(&letter.to_owned()));
            }
        } else if let AtomView::Fun(vector) = letter
            && let Some(last) = vector.iter().last()
            && matches!(last, AtomView::Fun(r) if r.get_symbol() == symbols.minkowski && r.get_nargs() == 1)
        {
            vector_representation = vector_representation.or_else(|| Some(last.to_owned()));
        }
        letters.push(letter.to_owned());
    }
    // closed: every explicit slot occurs exactly twice
    if !once_only(explicit).is_empty() {
        return None;
    }
    if vector_representation.is_none() {
        // only summed letters: read the stripped representation inside them
        for letter in &letters {
            letter.as_view().visitor(&mut |node| {
                if vector_representation.is_none()
                    && let AtomView::Fun(r) = node
                    && r.get_symbol() == symbols.minkowski
                    && r.get_nargs() == 1
                {
                    vector_representation = Some(node.to_owned());
                }
                vector_representation.is_none()
            });
        }
    }
    let mut opaque: Vec<(Atom, Atom)> = Vec::new();
    let word_letters: Vec<Atom> = letters
        .iter()
        .map(|letter| {
            if !matches!(letter.as_view(), AtomView::Add(_)) {
                return Some(letter.clone());
            }
            if let Some((stand_in, _)) = opaque.iter().find(|(_, sum)| sum == letter) {
                return Some(stand_in.clone());
            }
            let symbol = spenso::network::tags::SPENSO_TAG.rank_one_tensor_symbol(&format!(
                "idenso::closing::v{:016x}",
                stand_in_key(letter)
            ));
            let stand_in = symbolica::atom::FunctionBuilder::new(symbol)
                .add_arg(vector_representation.as_ref()?)
                .finish();
            opaque.push((stand_in.clone(), letter.clone()));
            Some(stand_in)
        })
        .collect::<Option<_>>()?;
    let settings = super::GammaSimplifySettings::canonical();
    let simplifier = super::DiracSimplifier::new(&settings);
    let closed =
        spenso::trace!(&representation; word_letters.iter().map(|letter| crate::gamma!(letter)));
    let (value, _) = simplifier.evaluate_terminal_trace::<false>(closed.as_view())?;
    #[cfg(test)]
    tests::CLOSED.with(|count| count.set(count.get() + 1));
    Some((value, opaque))
}

/// A deterministic key of a sum of momenta (its canonical spelling).
fn stand_in_key(sum: &Atom) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    sum.to_canonical_string().hash(&mut hasher);
    hasher.finish()
}

/// Put the sums back for their stand-ins (after the value was merged).
pub(crate) fn substitute_stand_ins(value: Atom, stand_ins: &[(Atom, Atom)]) -> Atom {
    if stand_ins.is_empty() {
        return value;
    }
    let replacements: Vec<symbolica::id::Replacement> = stand_ins
        .iter()
        .map(|(stand_in, sum)| symbolica::id::Replacement::new(stand_in.to_pattern(), sum.clone()))
        .collect();
    value.replace_multiple(&replacements)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use spenso::g;
    use spenso::network::tags::SPENSO_TAG as T;
    use spenso::trace;
    use symbolica::atom::{Atom, AtomCore, FunctionBuilder};

    use super::GammaClosing;
    use crate::dirac::GammaSimplifySettings;
    use crate::gamma;
    use crate::tensor::{AlgebraContraction, AlgebraSettings, SymbolicTensor};
    use crate::test_support::test_initialize;

    thread_local! {
        /// Words written by the GAMMA arm and closed words evaluated, on this thread.
        pub(super) static WORDS: Cell<u64> = const { Cell::new(0) };
        pub(super) static CLOSED: Cell<u64> = const { Cell::new(0) };
    }

    fn vector(name: &str, slot: &Atom) -> Atom {
        FunctionBuilder::new(T.rank_one_tensor_symbol(&format!("idenso::closing_test::{name}")))
            .add_arg(slot)
            .finish()
    }

    /// gamma(b0,b1,s0) gamma(b1,b2,s6) gamma(b2,b3,s1) gamma(b3,b0,s7) slash(s6) p2(s7)
    fn explicit_cycle_with(r: &crate::test_support::TestReps, s: &[Atom], slash: &Atom) -> Atom {
        explicit_cycle(r, s) / vector("p1", &s[6]) * slash
    }

    /// gamma(b0,b1,s0) gamma(b1,b2,s6) gamma(b2,b3,s1) gamma(b3,b0,s7) p1(s6) p2(s7)
    fn explicit_cycle(r: &crate::test_support::TestReps, s: &[Atom]) -> Atom {
        let b: Vec<_> = (0..4)
            .map(|index| r.bis4.pattern(Atom::num((20 + index) as i64)))
            .collect();
        let gamma = |a: &Atom, c: &Atom, mu: &Atom| {
            FunctionBuilder::new(crate::dirac::AGS.gamma)
                .add_arg(a)
                .add_arg(c)
                .add_arg(mu)
                .finish()
        };
        gamma(&b[0], &b[1], &s[0])
            * gamma(&b[1], &b[2], &s[6])
            * gamma(&b[2], &b[3], &s[1])
            * gamma(&b[3], &b[0], &s[7])
            * vector("p1", &s[6])
            * vector("p2", &s[7])
    }

    /// The expanded dot polynomial and the two certificates gammalooprs reads
    /// (closed scalar, terminal polynomial) of the simplified tensor.
    fn simplified_full(input: &Atom, closing: Option<GammaClosing>) -> (Atom, bool, bool) {
        let settings = AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            contract: AlgebraContraction::Fully,
            gamma_closing: closing,
            ..Default::default()
        };
        let tensor = SymbolicTensor::infer(input.clone())
            .unwrap()
            .simplify_algebra(&settings)
            .unwrap();
        let certificates = (
            tensor.is_certified_closed_scalar(),
            tensor.is_certified_terminal_polynomial(),
        );
        let value = tensor.to_dots().unwrap().into_expression().expand();
        (value, certificates.0, certificates.1)
    }

    fn simplified(input: &Atom, closing: Option<GammaClosing>) -> Atom {
        simplified_full(input, closing).0
    }

    /// The inputs: open words with a metric/vector sum, three-gluon vertex
    /// sums (one and two), an explicit gamma cycle, a routed-sum slash letter,
    /// an open word next to an unrelated interface sum, and a word whose own
    /// ports reach no sum.
    fn inputs() -> Vec<Atom> {
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let s: Vec<_> = (0..8)
            .map(|index| r.mink_d.pattern(Atom::num(index)))
            .collect();
        let d = r.mink_d.to_symbolic([]);
        let slash = |name: &str| vector(name, &d);
        let x = Atom::var(symbolica::symbol!("idenso::closing_test::x"));
        let p1 = slash("p1");
        let p2 = slash("p2");
        let open = trace!(&spin; [gamma!(&s[0]), gamma!(&p1), gamma!(&s[1]), gamma!(&p2)]);
        let open3 = trace!(&spin; [gamma!(&s[0]), gamma!(&p1), gamma!(&s[1]), gamma!(&p2), gamma!(&s[2]), gamma!(&p1)]);
        let vertex = |a: &Atom, b: &Atom, c: &Atom| {
            // a three-gluon-like vertex: g(a,b)(k1-k2)(c) + cyclic
            g!(a, b) * (vector("k1", c) - vector("k2", c))
                + g!(b, c) * (vector("k2", a) - vector("k3", a))
                + g!(c, a) * (vector("k3", b) - vector("k1", b))
        };
        vec![
            &open * (g!(&s[0], &s[1]) * &x + vector("p3", &s[0]) * vector("p4", &s[1])),
            &open3 * vertex(&s[0], &s[1], &s[2]),
            &open3
                * vertex(&s[0], &s[1], &s[5])
                * (g!(&s[5], &s[2]) * &x + vector("p4", &s[5]) * vector("p3", &s[2])),
            &open * vertex(&s[0], &s[4], &s[5]) * vertex(&s[1], &s[4], &s[5]),
            // an explicit gamma cycle (the planner's raw form): two vertex
            // letters, two slashes as separate vectors, closed by a vertex sum
            explicit_cycle(r, &s) * vertex(&s[0], &s[1], &s[5]) * vector("p4", &s[5]),
            // a routed-sum slash letter: (p3 + p4) on the cycle's second slot
            explicit_cycle_with(r, &s, &(vector("p3", &s[6]) + vector("p4", &s[6])))
                * vertex(&s[0], &s[1], &s[5])
                * vector("p2", &s[5]),
            // an open word next to an unrelated interface sum (s3 stays open)
            &open
                * (g!(&s[0], &s[1]) * &x + vector("p3", &s[0]) * vector("p4", &s[1]))
                * (vector("k1", &s[3]) + vector("k2", &s[3])),
            // a word whose own ports reach no sum, a sum elsewhere in the term:
            // written, never pending, evaluated by the ordinary gamma pass
            &open * g!(&s[0], &s[1]) * (vector("k1", &s[3]) + vector("k2", &s[3])),
        ]
    }

    /// Frontier closing gives the unclosed path's polynomial and both of its
    /// certificates (closed scalar, terminal polynomial), writes and closes
    /// words, and leaves no word under
    /// `AlgebraContraction::None` or `Minimal`.
    #[test]
    fn closing_equals_the_unclosed_path() {
        let inputs = inputs();
        WORDS.with(|count| count.set(0));
        CLOSED.with(|count| count.set(0));
        for input in &inputs {
            let (plain, closed_scalar, terminal) = simplified_full(input, None);
            let (closed, closed_closed_scalar, closed_terminal) =
                simplified_full(input, Some(GammaClosing::default()));
            assert_eq!((&closed - &plain).expand(), Atom::Zero, "input {input}");
            assert_eq!(
                (closed_closed_scalar, closed_terminal),
                (closed_scalar, terminal),
                "certificates, input {input}"
            );
            // STRUCTURE never enabled (None), or minimal: no word may stay
            // unevaluated
            for contract in [AlgebraContraction::None, AlgebraContraction::Minimal] {
                let run = |closing: Option<GammaClosing>| {
                    let settings = AlgebraSettings {
                        gamma: Some(GammaSimplifySettings::canonical()),
                        epsilon: true,
                        contract: contract.clone(),
                        gamma_closing: closing,
                        ..Default::default()
                    };
                    SymbolicTensor::infer(input.clone())
                        .unwrap()
                        .simplify_algebra(&settings)
                        .unwrap()
                        .into_expression()
                };
                let plain_run = run(None);
                let closed_run = run(Some(GammaClosing::default()));
                assert!(
                    !closed_run.as_view().to_string().contains("trace("),
                    "contract {contract:?} leaves a word: {closed_run}"
                );
                let difference =
                    (&simplified(&closed_run, None) - &simplified(&plain_run, None)).expand();
                assert_eq!(
                    difference,
                    Atom::Zero,
                    "contract {contract:?}, input {input}"
                );
            }
        }
        assert!(WORDS.with(Cell::get) > 0, "no word written");
        assert!(CLOSED.with(Cell::get) > 0, "no word closed");
    }

    /// A forced exhaustion of the cooperative resource budget aborts the
    /// closing contraction with the budget's panic; the thread's closing switch
    /// is restored, and the same input then gives the unclosed path's value.
    #[test]
    fn an_exhausted_budget_aborts_cleanly() {
        use symbolica::utils::resources::{ResourceBudget, ResourceLimits};
        let input = inputs().swap_remove(3);
        let budget = ResourceBudget::new(ResourceLimits {
            max_elapsed: Some(std::time::Duration::ZERO),
            ..Default::default()
        });
        let aborted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            budget.run(|| simplified(&input, Some(GammaClosing::default())))
        }));
        let payload = aborted.expect_err("an exhausted budget must abort");
        let message = payload
            .downcast_ref::<String>()
            .cloned()
            .or_else(|| payload.downcast_ref::<&str>().map(|text| text.to_string()))
            .unwrap_or_default();
        assert!(message.contains("resource budget exhausted"), "{message}");
        assert!(!super::frontier_active(), "closing switch left on");
        let closed = simplified(&input, Some(GammaClosing::default()));
        assert_eq!((&closed - &simplified(&input, None)).expand(), Atom::Zero);
    }
}
