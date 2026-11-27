//! Preserve tensor binding scopes while scalar expansion distributes terms.
use super::ExplicitIndexOccurrences;
use crate::tensor::SymbolicTensor;
use crate::tensor::inference::InterfaceInference;
use crate::tensor::network::RenamedLeafPorts;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::PartialStructure;
use spenso::structure::representation::LibraryRep;
use spenso::structure::slot::{IsAbstractSlot, Slot, SlotMatcher};
use std::collections::{HashMap, HashSet};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::coefficient::CoefficientView;

/// One dummy renaming, keyed by the base orientation of the renamed slots.
pub(in crate::tensor) type ScopeRenaming = HashMap<Slot<LibraryRep, AbstractIndex>, AbstractIndex>;

/// Fresh dummy names for scope separation. Each renamed leaf occurrence
/// carries the raw ports which its syntax cannot reproduce, so a producer can
/// publish them for exactly the leaves that the renamings rebuilt.
pub(in crate::tensor) struct FreshDummies<'a, F> {
    next: F,
    leaves: RenamedLeafPorts<'a>,
}

impl<'a, F: FnMut() -> AbstractIndex> FreshDummies<'a, F> {
    pub(in crate::tensor) fn new(next: F, leaves: RenamedLeafPorts<'a>) -> Self {
        Self { next, leaves }
    }

    /// The facts carried by the renamed leaves so far.
    pub(in crate::tensor) fn into_leaf_ports(self) -> RenamedLeafPorts<'a> {
        self.leaves
    }
}

impl ExplicitIndexOccurrences {
    /// Keep each closed scope which an expansion distributes in its own factors.
    ///
    /// A sum binds its internal pairs, so the same name may also occur in a
    /// neighbouring factor or in another term of that sum. Expansion moves
    /// those scopes into one product, where equal names would contract across
    /// the former boundary. Only scopes which the selected expansion actually
    /// distributes are rewritten, following the backend's rules: with `var`,
    /// subexpressions without it stay untouched, and only integer powers whose
    /// magnitude fits in `u32` are expanded. A colliding closed pair of a
    /// distributed sum is renamed with fresh `dummies`, which carry the raw
    /// ports of each renamed leaf. A distributed power of a closed scalar sum
    /// keeps the terms' names disjoint and brackets each closed term, so the
    /// multinomial expansion keeps every pair inside its own power scope.
    /// Returns `None` when nothing can collide.
    pub(in crate::tensor) fn separate_expansion_scopes(
        value: AtomView<'_>,
        var: Option<AtomView<'_>>,
        via_poly: bool,
        dummies: &mut FreshDummies<'_, impl FnMut() -> AbstractIndex>,
    ) -> Option<Atom> {
        let mut slots = SlotMatcher::default();
        if via_poly {
            Self::separate_polynomial_scopes(value, var, &mut slots, dummies)
        } else {
            Self::separate_scopes(value, var, var, &mut slots, dummies)
        }
    }

    /// Polynomial expansion keeps expanded terms and terms without `var`, and
    /// converts every other term of the outer sum as a whole.
    fn separate_polynomial_scopes(
        value: AtomView<'_>,
        var: Option<AtomView<'_>>,
        slots: &mut SlotMatcher,
        dummies: &mut FreshDummies<'_, impl FnMut() -> AbstractIndex>,
    ) -> Option<Atom> {
        if let AtomView::Add(sum) = value {
            return Self::separate_terms(sum.iter(), |term| {
                Self::separate_polynomial_scopes(term, var, slots, dummies)
            });
        }
        if AtomCore::is_expanded(&value, var) || !Self::expands_in(value, var) {
            return None;
        }
        // Coefficients of a selected variable may still be multiplied out.
        Self::separate_scopes(value, None, var, slots, dummies)
    }

    fn separate_terms<'a>(
        terms: impl Iterator<Item = AtomView<'a>>,
        separate: impl FnMut(AtomView<'a>) -> Option<Atom>,
    ) -> Option<Atom> {
        let terms = terms.collect::<Vec<_>>();
        let separated = terms.iter().copied().map(separate).collect::<Vec<_>>();
        if separated.iter().all(Option::is_none) {
            return None;
        }
        Some(Atom::add_many(terms.into_iter().zip(separated).map(
            |(term, separated)| separated.unwrap_or_else(|| term.to_owned()),
        )))
    }

    /// Whether expansion in `var` visits `value`, as the backend decides.
    fn expands_in(value: AtomView<'_>, var: Option<AtomView<'_>>) -> bool {
        match var {
            None => true,
            Some(AtomView::Var(variable)) => {
                AtomCore::contains_symbol(&value, variable.get_symbol())
            }
            Some(var) => AtomCore::contains(&value, var),
        }
    }

    /// The integer power which the backend expands, with its sign.
    fn expanded_power(exponent: AtomView<'_>) -> Option<(bool, u32)> {
        if let AtomView::Num(number) = exponent
            && let CoefficientView::Natural(power, 1, 0, 1) = number.get_coeff_view()
        {
            return u32::try_from(power.unsigned_abs())
                .ok()
                .map(|magnitude| (power < 0, magnitude));
        }
        None
    }

    /// Whether a product factor is distributed over the other factors.
    fn distributes(factor: AtomView<'_>, var: Option<AtomView<'_>>) -> bool {
        match factor {
            AtomView::Add(_) => Self::expands_in(factor, var),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                matches!(base, AtomView::Add(_))
                    && Self::expands_in(base, var)
                    && Self::expanded_power(exponent)
                        .is_some_and(|(negative, power)| !negative && power > 1)
            }
            _ => false,
        }
    }

    fn separate_scopes(
        value: AtomView<'_>,
        var: Option<AtomView<'_>>,
        preserve: Option<AtomView<'_>>,
        slots: &mut SlotMatcher,
        dummies: &mut FreshDummies<'_, impl FnMut() -> AbstractIndex>,
    ) -> Option<Atom> {
        if !Self::expands_in(value, var) {
            return None;
        }
        match value {
            AtomView::Add(sum) => Self::separate_terms(sum.iter(), |term| {
                Self::separate_scopes(term, var, preserve, slots, dummies)
            }),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let separated_base = Self::separate_scopes(base, var, preserve, slots, dummies);
                let separated_exponent =
                    Self::separate_scopes(exponent, var, preserve, slots, dummies);
                let base_view = separated_base.as_ref().map_or(base, Atom::as_view);
                let exponent_view = separated_exponent.as_ref().map_or(exponent, Atom::as_view);
                // An admitted square has exactly two copies. Freshen only
                // its closed internal pairs in those two copies;
                // closed scalar powers keep the compact multinomial path below.
                if matches!(base_view, AtomView::Add(_))
                    && Self::expanded_power(exponent_view) == Some((false, 2))
                {
                    let closed = Self::closed_slots(base_view, slots);
                    if !closed.is_empty() && closed != Self::written_slots(base_view, slots) {
                        let left = Self::renamed(base_view, &closed, slots, dummies);
                        let right = Self::renamed(base_view, &closed, slots, dummies);
                        return Some(Atom::mul_many([left, right]));
                    }
                }
                let raised = match (base_view, Self::expanded_power(exponent_view)) {
                    (AtomView::Add(sum), Some((_, power))) if power > 1 => {
                        Self::separate_raised_terms(
                            base_view,
                            sum.iter(),
                            var,
                            preserve,
                            slots,
                            dummies,
                        )
                    }
                    _ => None,
                };
                if raised.is_none() && separated_base.is_none() && separated_exponent.is_none() {
                    return None;
                }
                Some(
                    raised
                        .as_ref()
                        .map_or(base_view, Atom::as_view)
                        .pow(exponent_view),
                )
            }
            AtomView::Mul(product) => {
                let mut changed = false;
                let mut factors = Vec::new();
                for factor in product.iter() {
                    match Self::separate_scopes(factor, var, preserve, slots, dummies) {
                        Some(separated) => {
                            changed = true;
                            factors.push(separated);
                        }
                        None => factors.push(factor.to_owned()),
                    }
                }
                if factors
                    .iter()
                    .any(|factor| Self::distributes(factor.as_view(), var))
                {
                    let mut written = factors
                        .iter()
                        .map(|factor| Self::written_slots(factor.as_view(), slots))
                        .collect::<Vec<_>>();
                    for position in 0..factors.len() {
                        if !Self::distributes(factors[position].as_view(), var) {
                            continue;
                        }
                        let mut shared = written[position]
                            .iter()
                            .filter(|slot| {
                                written
                                    .iter()
                                    .enumerate()
                                    .any(|(other, names)| other != position && names.contains(slot))
                            })
                            .copied()
                            .collect::<HashSet<_>>();
                        if shared.is_empty() {
                            continue;
                        }
                        let closed = Self::closed_slots(factors[position].as_view(), slots);
                        shared.retain(|slot| closed.contains(slot));
                        if shared.is_empty() {
                            continue;
                        }
                        factors[position] =
                            Self::renamed(factors[position].as_view(), &shared, slots, dummies);
                        written[position] = Self::written_slots(factors[position].as_view(), slots);
                        changed = true;
                    }
                }
                changed.then(|| Atom::mul_many(&factors))
            }
            _ => None,
        }
    }

    /// Prepare a closed scalar sum for a multinomial power. The backend raises
    /// each term containing `var` separately, together with a sole remaining
    /// term; several remaining terms stay one closed sum. Raised terms receive
    /// disjoint dummy names, and each indexed product is bracketed so that its
    /// powers remain closed scalars instead of distributing over its factors.
    fn separate_raised_terms<'a>(
        base: AtomView<'a>,
        terms: impl Iterator<Item = AtomView<'a>>,
        var: Option<AtomView<'_>>,
        preserve: Option<AtomView<'_>>,
        slots: &mut SlotMatcher,
        dummies: &mut FreshDummies<'_, impl FnMut() -> AbstractIndex>,
    ) -> Option<Atom> {
        let closed = Self::closed_slots(base, slots);
        if closed.is_empty() || closed != Self::written_slots(base, slots) {
            return None;
        }
        let (mut raised, rest): (Vec<_>, Vec<_>) =
            terms.partition(|term| Self::expands_in(*term, var));
        let kept = if let [term] = rest.as_slice() {
            raised.push(*term);
            Vec::new()
        } else {
            rest
        };
        let mut changed = false;
        let mut used = HashSet::new();
        let mut terms = kept.iter().map(|term| term.to_owned()).collect::<Vec<_>>();
        for term in raised {
            let written = Self::written_slots(term, slots);
            if written.is_empty() {
                terms.push(term.to_owned());
                continue;
            }
            let shared = written.intersection(&used).copied().collect::<HashSet<_>>();
            let renamed =
                (!shared.is_empty()).then(|| Self::renamed(term, &shared, slots, dummies));
            let term = renamed.as_ref().map_or(term, Atom::as_view);
            used.extend(Self::written_slots(term, slots));
            let scalar = Self::raised_scalar_expression(term, preserve, slots, dummies);
            changed |= renamed.is_some() || scalar.as_view() != term;
            terms.push(scalar);
        }
        changed.then(|| Atom::add_many(&terms))
    }

    /// Bracket indexed factors together while keeping the selected scalar
    /// coefficient visible to polynomial conversion. Already closed selected
    /// sums remain outside; each of their indexed terms is bracketed, since
    /// the backend raises such a sum with its term rather than as a power base.
    fn raised_scalar_expression(
        term: AtomView<'_>,
        preserve: Option<AtomView<'_>>,
        slots: &mut SlotMatcher,
        dummies: &mut FreshDummies<'_, impl FnMut() -> AbstractIndex>,
    ) -> Atom {
        if let (AtomView::Mul(product), Some(selected)) = (term, preserve) {
            let mut coefficients = Vec::new();
            let mut indexed = Vec::new();
            for factor in product.iter() {
                let written = Self::written_slots(factor, slots);
                if Self::expands_in(factor, Some(selected))
                    && (written.is_empty()
                        || Self::closed_slots(factor, slots) == written
                        || (matches!(factor, AtomView::Fun(_))
                            && dummies.leaves.is_closed_leaf(factor)))
                {
                    // This closed sum is raised with its term, not as the base
                    // of a power, so the backend does not keep its terms
                    // without the selector together: bracket each indexed term.
                    let prepared = if let AtomView::Add(sum) = factor {
                        Self::separate_raised_terms(
                            factor,
                            sum.iter(),
                            None,
                            preserve,
                            slots,
                            dummies,
                        )
                    } else if written.is_empty() {
                        None
                    } else {
                        Some(SymbolicTensor::<PartialStructure>::scalar_power_expression(
                            factor.to_owned(),
                            Some(true),
                        ))
                    };
                    coefficients.push(prepared.unwrap_or_else(|| factor.to_owned()));
                } else {
                    indexed.push(factor.to_owned());
                }
            }
            if !coefficients.is_empty() {
                if !indexed.is_empty() {
                    coefficients.push(SymbolicTensor::<PartialStructure>::scalar_power_expression(
                        Atom::mul_many(indexed),
                        Some(true),
                    ));
                }
                return Atom::mul_many(coefficients);
            }
        }
        SymbolicTensor::<PartialStructure>::scalar_power_expression(term.to_owned(), Some(true))
    }

    /// Explicit identities written anywhere in `value`, keyed by base orientation.
    /// Only arguments which inference reads as indices count: those of tagged
    /// tensors and composite heads. Any other function, including an
    /// unregistered one with index-shaped arguments, is an opaque scalar.
    fn written_slots(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> HashSet<Slot<LibraryRep, AbstractIndex>> {
        let mut written = HashSet::new();
        Self::insert_written_slots(value, slots, &mut written);
        written
    }

    fn insert_written_slots(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        written: &mut HashSet<Slot<LibraryRep, AbstractIndex>>,
    ) {
        if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(value) {
            written.insert(slot.rep().base().slot(slot.aind()));
            return;
        }
        match value {
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                for argument in function.iter() {
                    if InterfaceInference::explicit_index_argument(symbol, argument, slots) {
                        Self::insert_written_slots(argument, slots, written);
                    }
                }
            }
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    Self::insert_written_slots(term, slots, written);
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    Self::insert_written_slots(factor, slots, written);
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::insert_written_slots(base, slots, written);
                Self::insert_written_slots(exponent, slots, written);
            }
            _ => {}
        }
    }

    /// Written identities which `value` binds internally: none of them is
    /// visible to a surrounding product.
    fn closed_slots(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> HashSet<Slot<LibraryRep, AbstractIndex>> {
        let mut closed = Self::written_slots(value, slots);
        if !closed.is_empty() {
            let exported = Self::from_atom(value, slots);
            closed.retain(|slot| exported.count(*slot) == 0);
        }
        closed
    }

    fn renamed(
        value: AtomView<'_>,
        closed: &HashSet<Slot<LibraryRep, AbstractIndex>>,
        slots: &mut SlotMatcher,
        dummies: &mut FreshDummies<'_, impl FnMut() -> AbstractIndex>,
    ) -> Atom {
        let renaming = closed
            .iter()
            .map(|slot| (*slot, (dummies.next)()))
            .collect::<ScopeRenaming>();
        Self::rename_index_arguments(
            value,
            &renaming,
            slots,
            &mut dummies.leaves,
            InterfaceInference::explicit_index_argument,
        )
        .unwrap_or_else(|| value.to_owned())
    }

    /// Rename explicit indices and carry proofs when rebuilding a leaf.
    /// `indexed` selects which function arguments participate in the rewrite.
    pub(in crate::tensor) fn rename_index_arguments(
        value: AtomView<'_>,
        renaming: &ScopeRenaming,
        slots: &mut SlotMatcher,
        leaves: &mut RenamedLeafPorts<'_>,
        indexed: fn(Symbol, AtomView<'_>, &mut SlotMatcher) -> bool,
    ) -> Option<Atom> {
        if value.get_symbol().is_some_and(|symbol| symbol.is_scalar()) {
            return None;
        }
        if let Ok(mut slot) = slots.parse::<LibraryRep, AbstractIndex>(value) {
            slot.aind = *renaming.get(&slot.rep().base().slot(slot.aind()))?;
            return Some(slot.to_atom());
        }
        let rename = |child| Self::rename_index_arguments(child, renaming, slots, leaves, indexed);
        match value {
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                let arguments = Self::renamed_children(function.iter(), |argument| {
                    if indexed(symbol, argument, slots) {
                        Self::rename_index_arguments(argument, renaming, slots, leaves, indexed)
                    } else {
                        None
                    }
                })?;
                let rebuilt = FunctionBuilder::new(symbol).add_args(&arguments).finish();
                leaves.observe(function, &arguments, rebuilt.as_view(), |port| {
                    let mut slot = slots.parse::<LibraryRep, AbstractIndex>(port).ok()?;
                    slot.aind = *renaming.get(&slot.rep().base().slot(slot.aind()))?;
                    Some(slot.to_atom())
                });
                Some(rebuilt)
            }
            AtomView::Add(sum) => Self::renamed_children(sum.iter(), rename).map(Atom::add_many),
            AtomView::Mul(product) => {
                Self::renamed_children(product.iter(), rename).map(Atom::mul_many)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::renamed_children([base, exponent].into_iter(), rename)
                    .map(|parts| parts[0].pow(&parts[1]))
            }
            _ => None,
        }
    }

    /// Only allocate a rebuilt argument list after the first change.
    fn renamed_children<'a>(
        children: impl Iterator<Item = AtomView<'a>> + Clone,
        mut rename: impl FnMut(AtomView<'a>) -> Option<Atom>,
    ) -> Option<Vec<Atom>> {
        let mut output = None;
        for (position, child) in children.clone().enumerate() {
            match rename(child) {
                Some(renamed) => output
                    .get_or_insert_with(|| {
                        children
                            .clone()
                            .take(position)
                            .map(|node| node.to_owned())
                            .collect::<Vec<_>>()
                    })
                    .push(renamed),
                None if let Some(output) = &mut output => output.push(child.to_owned()),
                None => {}
            }
        }
        output
    }
}
