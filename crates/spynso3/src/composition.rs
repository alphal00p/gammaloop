use std::collections::{HashMap, HashSet};

use spenso::{
    network::{
        parsing::{
            AtomStructureExt, ChainNestingError, StrictTensorFilter, StructureInferenceMode,
        },
        tags::SPENSO_TAG,
    },
    shadowing,
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialSlot, PartialStructure, PartialStructureExt},
        representation::{LibraryRep, RepName, Representation},
        slot::{DummyAind, IsAbstractSlot, ParseableAind, Slot},
    },
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use thiserror::Error;

/// An atom paired with its ordered, tensor-aware external interface.
#[derive(Clone, Debug)]
pub struct StructuredAtom {
    pub atom: Atom,
    pub interface: PartialStructure,
}

impl StructuredAtom {
    pub fn new(atom: Atom, interface: PartialStructure) -> Self {
        Self {
            atom,
            interface: interface.canonicalize_open_ports(),
        }
    }

    pub fn rank(&self) -> usize {
        self.interface.canonical().order()
    }

    pub fn is_scalar(&self) -> bool {
        self.rank() == 0
    }

    /// Build a presentation-only atom whose surviving tensor ports follow the
    /// public logical interface rather than `OrderedStructure`'s canonical
    /// storage order.
    pub(crate) fn presentation_atom(&self) -> Atom {
        let slots = self.interface.logical_slots();
        let mut state = PortRewriteState::new(slots.len());
        reorder_presentation_ports(self.atom.as_view(), &slots, &mut state, false)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PortPair {
    pub left: usize,
    pub right: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MatrixChannel {
    pub input: usize,
    pub output: usize,
}

#[derive(Debug, Error)]
pub enum TensorCompositionError {
    #[error("port position {position} is outside an interface of rank {rank}")]
    InvalidPort { position: usize, rank: usize },
    #[error("ports {left} and {right} carry incompatible representations")]
    IncompatiblePorts { left: usize, right: usize },
    #[error("explicit ports {left} and {right} do not carry the same abstract index")]
    UnequalExplicitIndices { left: usize, right: usize },
    #[error("tensor multiplication is ambiguous; compatible port pairs are {0:?}")]
    Ambiguous(Vec<PortPair>),
    #[error("the selected tensor does not have an unambiguous two-ended channel")]
    NoMatrixChannel,
    #[error("a two-ended channel requires distinct ports, got ({input}, {output})")]
    DegenerateChannel { input: usize, output: usize },
    #[error("ports ({input}, {output}) do not form an input-to-output propagation channel")]
    InvalidChannelOrientation { input: usize, output: usize },
    #[error("interface port {position} could not be located consistently in the tensor atom")]
    MissingInterfacePort { position: usize },
    #[error("explicit index `{index}` occurs on {occurrences} compatible tensor ports")]
    InvalidExplicitMultiplicity {
        index: AbstractIndex,
        occurrences: usize,
    },
    #[error(transparent)]
    InvalidChainNesting(#[from] ChainNestingError),
    #[error(
        "composition would nest an existing chain or trace; only the primary channel of a root chain can be extended"
    )]
    NestedChainLike,
}

fn representation(slot: &PartialSlot) -> Representation<LibraryRep> {
    slot.rep()
}

fn representations_match(left: &PartialSlot, right: &PartialSlot) -> bool {
    representation(left).matches(&representation(right))
}

fn automatically_contractible(left: &PartialSlot, right: &PartialSlot) -> bool {
    if !representations_match(left, right) {
        return false;
    }

    match (left.aind, right.aind) {
        (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) => left == right,
        _ => true,
    }
}

pub fn compatible_pairs(left: &PartialStructure, right: &PartialStructure) -> Vec<PortPair> {
    let left = left.logical_slots();
    let right = right.logical_slots();
    left.iter()
        .enumerate()
        .flat_map(|(left_position, left_slot)| {
            right
                .iter()
                .enumerate()
                .filter(move |(_, right_slot)| automatically_contractible(left_slot, right_slot))
                .map(move |(right_position, _)| PortPair {
                    left: left_position,
                    right: right_position,
                })
        })
        .collect()
}

fn is_structured_scalar(value: AtomView<'_>) -> bool {
    if matches!(value, AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.dot) {
        return true;
    }

    value
        .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>(
            StructureInferenceMode::Fast,
        )
        .is_ok_and(|structure| structure.canonical().is_scalar())
}

fn port_matches(value: AtomView<'_>, expected: PartialSlot) -> bool {
    match expected.aind {
        PartialIndex::Explicit(index) => Slot::<LibraryRep, AbstractIndex>::try_from(value)
            .is_ok_and(|slot| slot.rep() == expected.rep() && slot.aind() == index),
        PartialIndex::Open(_) => Representation::<LibraryRep>::try_from(value)
            .is_ok_and(|representation| representation == expected.rep()),
    }
}

fn chain_endpoints_match(value: AtomView<'_>, input: PartialSlot, output: PartialSlot) -> bool {
    let AtomView::Fun(function) = value else {
        return false;
    };
    if function.get_symbol() != SPENSO_TAG.chain {
        return false;
    }
    let mut arguments = function.iter();
    matches!(arguments.next(), Some(value) if port_matches(value, input))
        && matches!(arguments.next(), Some(value) if port_matches(value, output))
}

fn has_transparent_chain_channel(
    value: AtomView<'_>,
    input: PartialSlot,
    output: PartialSlot,
) -> bool {
    match value {
        AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain => {
            chain_endpoints_match(value, input, output)
        }
        AtomView::Fun(function) if function.get_symbol().has_tag(&SPENSO_TAG.broadcast) => {
            matches!(function.iter().collect::<Vec<_>>().as_slice(), [argument] if has_transparent_chain_channel(*argument, input, output))
        }
        AtomView::Fun(function)
            if function.get_symbol() == *shadowing::SYM
                || function.get_symbol() == *shadowing::ANTISYM
                || function.get_symbol() == *shadowing::CYCLIC =>
        {
            let tensor_arguments = function
                .iter()
                .filter(|argument| {
                    argument.is_tensorial(StrictTensorFilter::Tagged)
                        || has_transparent_chain_channel(*argument, input, output)
                })
                .collect::<Vec<_>>();
            matches!(tensor_arguments.as_slice(), [argument] if has_transparent_chain_channel(*argument, input, output))
        }
        AtomView::Mul(product) => {
            let tensor_factors = product
                .iter()
                .filter(|factor| {
                    (factor.is_tensorial(StrictTensorFilter::Tagged)
                        || matches!(
                            factor,
                            AtomView::Fun(function)
                                if function.get_symbol() == *shadowing::SYM
                                    || function.get_symbol() == *shadowing::ANTISYM
                                    || function.get_symbol() == *shadowing::CYCLIC
                        ))
                        && !is_structured_scalar(*factor)
                })
                .collect::<Vec<_>>();
            matches!(tensor_factors.as_slice(), [factor] if has_transparent_chain_channel(*factor, input, output))
        }
        AtomView::Add(sum) => {
            let terms = sum.iter().collect::<Vec<_>>();
            !terms.is_empty()
                && terms
                    .into_iter()
                    .all(|term| has_transparent_chain_channel(term, input, output))
        }
        _ => false,
    }
}

fn contains_chain_like(value: AtomView<'_>) -> bool {
    match value {
        AtomView::Add(sum) => sum.iter().any(contains_chain_like),
        AtomView::Mul(product) => product.iter().any(contains_chain_like),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            contains_chain_like(base) || contains_chain_like(exponent)
        }
        AtomView::Fun(function) => {
            function.get_symbol() == SPENSO_TAG.chain
                || function.get_symbol() == SPENSO_TAG.trace
                || function.iter().any(contains_chain_like)
        }
        _ => false,
    }
}

fn root_chain_channel_is_live(value: &StructuredAtom, channel: MatrixChannel) -> bool {
    if channel
        != (MatrixChannel {
            input: 0,
            output: 1,
        })
    {
        return false;
    }
    let slots = value.interface.logical_slots();
    slots.len() >= 2 && chain_endpoints_match(value.atom.as_view(), slots[0], slots[1])
}

pub fn matrix_channel(value: &StructuredAtom) -> Option<MatrixChannel> {
    let slots = value.interface.logical_slots();
    if slots.len() >= 2 && has_transparent_chain_channel(value.atom.as_view(), slots[0], slots[1]) {
        return Some(MatrixChannel {
            input: 0,
            output: 1,
        });
    }

    let pairs = (0..slots.len())
        .flat_map(|left| ((left + 1)..slots.len()).map(move |right| (left, right)))
        .filter(|&(left, right)| representations_match(&slots[left], &slots[right]))
        .collect::<Vec<_>>();
    let [(first, second)] = pairs.as_slice() else {
        return None;
    };

    let first_rep = representation(&slots[*first]);
    let second_rep = representation(&slots[*second]);
    if first_rep.rep.is_self_dual() || (first_rep.rep.is_base() && second_rep.rep.is_dual()) {
        Some(MatrixChannel {
            input: *first,
            output: *second,
        })
    } else {
        Some(MatrixChannel {
            input: *second,
            output: *first,
        })
    }
}

fn validate_position(
    slots: &[PartialSlot],
    position: usize,
) -> Result<PartialSlot, TensorCompositionError> {
    slots
        .get(position)
        .copied()
        .ok_or(TensorCompositionError::InvalidPort {
            position,
            rank: slots.len(),
        })
}

fn shared_index(
    left: PartialSlot,
    right: PartialSlot,
    positions: PortPair,
) -> Result<Option<AbstractIndex>, TensorCompositionError> {
    if !representations_match(&left, &right) {
        return Err(TensorCompositionError::IncompatiblePorts {
            left: positions.left,
            right: positions.right,
        });
    }

    match (left.aind, right.aind) {
        (PartialIndex::Explicit(left), PartialIndex::Explicit(right)) if left == right => {
            Ok(Some(left))
        }
        (PartialIndex::Explicit(_), PartialIndex::Explicit(_)) => {
            Err(TensorCompositionError::UnequalExplicitIndices {
                left: positions.left,
                right: positions.right,
            })
        }
        (PartialIndex::Explicit(index), PartialIndex::Open(_))
        | (PartialIndex::Open(_), PartialIndex::Explicit(index)) => Ok(Some(index)),
        (PartialIndex::Open(_), PartialIndex::Open(_)) => Ok(None),
    }
}

/// Allocate a dummy whose serialized index symbol is absent from the supplied
/// atoms and interfaces.
///
/// `AbstractIndex::Dummy(n)` is stored in Symbolica atoms as the symbol `d_n`,
/// which parses back as an `AbstractIndex::Symbol`. Comparing the round-tripped
/// value prevents a generated contraction index from aliasing a user-supplied
/// explicit `d_n` index.
pub(crate) fn fresh_dummy_index<'a>(
    atoms: impl IntoIterator<Item = &'a Atom>,
    interfaces: impl IntoIterator<Item = &'a PartialStructure>,
) -> AbstractIndex {
    let mut occupied = HashSet::new();
    for atom in atoms {
        let _ = atom.replace_map(|value, _, _| {
            if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
                occupied.insert(slot.aind());
            }
        });
    }
    for interface in interfaces {
        occupied.extend(interface.logical_slots().into_iter().filter_map(|slot| {
            if let PartialIndex::Explicit(index) = slot.aind {
                let atom = index.to_atom();
                AbstractIndex::try_from(atom.as_view()).ok()
            } else {
                None
            }
        }));
    }

    loop {
        let atom = AbstractIndex::new_dummy().to_atom();
        let marker = AbstractIndex::try_from(atom.as_view())
            .expect("fresh dummy atoms must remain valid abstract indices");
        if !occupied.contains(&marker) {
            return marker;
        }
    }
}

fn is_composite_head(symbol: Symbol) -> bool {
    symbol == SPENSO_TAG.bracket
        || symbol == SPENSO_TAG.chain
        || symbol == SPENSO_TAG.trace
        || symbol.has_tag(&SPENSO_TAG.broadcast)
        || symbol == *shadowing::SYM
        || symbol == *shadowing::ANTISYM
        || symbol == *shadowing::CYCLIC
}

fn is_tensor_leaf_head(symbol: Symbol) -> bool {
    symbol.has_tag(&SPENSO_TAG.tensor) && !is_composite_head(symbol)
}

fn direct_structural_port(value: AtomView<'_>) -> bool {
    Slot::<LibraryRep, AbstractIndex>::try_from(value).is_ok()
        || Representation::<LibraryRep>::try_from(value).is_ok()
}

fn explicit_occurrences(value: AtomView<'_>, target: Slot<LibraryRep, AbstractIndex>) -> usize {
    if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
        return usize::from(
            slot.aind() == target.aind()
                && (slot.rep() == target.rep() || slot.rep().matches(&target.rep())),
        );
    }

    match value {
        AtomView::Add(sum) => sum
            .iter()
            .map(|term| explicit_occurrences(term, target))
            .max()
            .unwrap_or_default(),
        AtomView::Mul(product) => product.iter().fold(0, |occurrences, factor| {
            occurrences
                .saturating_add(explicit_occurrences(factor, target))
                .min(3)
        }),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            explicit_occurrences(base, target).max(explicit_occurrences(exponent, target))
        }
        AtomView::Fun(function) => function.iter().fold(0, |occurrences, argument| {
            let structured = direct_structural_port(argument)
                || argument.is_tensorial(StrictTensorFilter::Tagged)
                || matches!(
                    argument,
                    AtomView::Fun(nested)
                        if nested.get_symbol() == *shadowing::SYM
                            || nested.get_symbol() == *shadowing::ANTISYM
                            || nested.get_symbol() == *shadowing::CYCLIC
                );
            if (is_tensor_leaf_head(function.get_symbol())
                || is_composite_head(function.get_symbol()))
                && structured
            {
                occurrences
                    .saturating_add(explicit_occurrences(argument, target))
                    .min(3)
            } else {
                occurrences
            }
        }),
        _ => 0,
    }
}

/// Reject Einstein indices that occur more than twice in any additive branch.
///
/// The walk follows tensor-bearing arguments only, so representation-valued
/// scalar metadata is not mistaken for a structural port. Additive branches
/// are checked independently because each summand owns the same public
/// interface rather than contributing another occurrence of it.
pub(crate) fn validate_explicit_index_occurrences(
    atom: &Atom,
) -> Result<(), TensorCompositionError> {
    let mut candidates = HashSet::new();
    let _ = atom.replace_map(|value, _, _| {
        if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
            candidates.insert(slot);
        }
    });
    for target in candidates {
        let occurrences = explicit_occurrences(atom.as_view(), target);
        if occurrences > 2 {
            return Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: target.aind(),
                occurrences,
            });
        }
    }
    Ok(())
}

#[derive(Clone)]
struct PortRewriteState {
    claimed: Vec<bool>,
    applied: Vec<bool>,
}

impl PortRewriteState {
    fn new(rank: usize) -> Self {
        Self {
            claimed: vec![false; rank],
            applied: vec![false; rank],
        }
    }
}

fn matching_interface_position(
    value: AtomView<'_>,
    slots: &[PartialSlot],
    claimed: &[bool],
) -> Option<usize> {
    if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(value) {
        return slots.iter().enumerate().position(|(position, expected)| {
            !claimed[position]
                && expected.rep() == slot.rep()
                && matches!(expected.aind, PartialIndex::Explicit(index) if index == slot.aind())
        });
    }

    let representation = Representation::<LibraryRep>::try_from(value).ok()?;
    slots.iter().enumerate().position(|(position, expected)| {
        !claimed[position]
            && expected.rep() == representation
            && matches!(expected.aind, PartialIndex::Open(_))
    })
}

fn reorder_presentation_ports(
    value: AtomView<'_>,
    slots: &[PartialSlot],
    state: &mut PortRewriteState,
    preserve_additive_order: bool,
) -> Atom {
    match value {
        AtomView::Add(add) => {
            let initial = state.clone();
            let mut common = None::<PortRewriteState>;
            let sum = add.iter().fold(Atom::Zero, |sum, term| {
                let mut local = initial.clone();
                let term = reorder_presentation_ports(term, slots, &mut local, true);
                if let Some(common) = &mut common {
                    for (common, local) in common.claimed.iter_mut().zip(&local.claimed) {
                        *common &= local;
                    }
                } else {
                    common = Some(local);
                }
                sum + term
            });
            *state = common.unwrap_or(initial);
            sum
        }
        AtomView::Mul(mul) => mul.iter().fold(Atom::num(1), |product, factor| {
            let tensorial = factor.is_tensorial(StrictTensorFilter::Tagged)
                || matches!(
                    factor,
                    AtomView::Fun(fun)
                        if fun.get_symbol() == *shadowing::SYM
                            || fun.get_symbol() == *shadowing::ANTISYM
                            || fun.get_symbol() == *shadowing::CYCLIC
                );
            product
                * if tensorial {
                    reorder_presentation_ports(factor, slots, state, preserve_additive_order)
                } else {
                    factor.to_owned()
                }
        }),
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            reorder_presentation_ports(base, slots, state, preserve_additive_order)
                .pow(exponent.to_owned())
        }
        // A compact dot has no public ports of its own.
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot => value.to_owned(),
        AtomView::Fun(fun) if is_tensor_leaf_head(fun.get_symbol()) => {
            let mut arguments = fun
                .iter()
                .map(|argument| argument.to_owned())
                .collect::<Vec<_>>();
            let matched = arguments
                .iter()
                .enumerate()
                .filter_map(|(argument_position, argument)| {
                    let interface_position =
                        matching_interface_position(argument.as_view(), slots, &state.claimed)?;
                    state.claimed[interface_position] = true;
                    Some((argument_position, interface_position, argument.clone()))
                })
                .collect::<Vec<_>>();
            // Within a sum, each term's explicit argument order records its map to the
            // shared interface; normalizing it would turn A(i,j) + A(j,i) into A(i,j) + A(i,j).
            if !preserve_additive_order {
                let mut ordered = matched.clone();
                ordered.sort_by_key(|(_, interface_position, _)| *interface_position);
                for ((argument_position, _, _), (_, _, argument)) in
                    matched.into_iter().zip(ordered)
                {
                    arguments[argument_position] = argument;
                }
            }
            FunctionBuilder::new(fun.get_symbol())
                .add_args(arguments)
                .finish()
        }
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace => {
            let mut arguments = fun.iter();
            let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
            if let Some(representation) = arguments.next() {
                rebuilt = rebuilt.add_arg(representation);
            }
            for argument in arguments {
                rebuilt = rebuilt.add_arg(reorder_presentation_ports(
                    argument,
                    slots,
                    state,
                    preserve_additive_order,
                ));
            }
            rebuilt.finish()
        }
        AtomView::Fun(fun) => {
            let symbol = fun.get_symbol();
            let mut rebuilt = FunctionBuilder::new(symbol);
            for argument in fun.iter() {
                if direct_structural_port(argument) {
                    if let Some(position) =
                        matching_interface_position(argument, slots, &state.claimed)
                    {
                        state.claimed[position] = true;
                    }
                    rebuilt = rebuilt.add_arg(argument);
                } else if is_composite_head(symbol)
                    && (argument.is_tensorial(StrictTensorFilter::Tagged)
                        || matches!(
                            argument,
                            AtomView::Fun(nested)
                                if nested.get_symbol() == *shadowing::SYM
                                    || nested.get_symbol() == *shadowing::ANTISYM
                                    || nested.get_symbol() == *shadowing::CYCLIC
                        ))
                {
                    rebuilt = rebuilt.add_arg(reorder_presentation_ports(
                        argument,
                        slots,
                        state,
                        preserve_additive_order,
                    ));
                } else {
                    rebuilt = rebuilt.add_arg(argument);
                }
            }
            rebuilt.finish()
        }
        _ => value.to_owned(),
    }
}

fn rewrite_ports(
    value: AtomView<'_>,
    slots: &[PartialSlot],
    replacements: &HashMap<usize, Atom>,
    state: &mut PortRewriteState,
) -> Atom {
    match value {
        AtomView::Add(add) => {
            let initial = state.clone();
            let mut common = None::<PortRewriteState>;
            let sum = add.iter().fold(Atom::Zero, |sum, term| {
                let mut local = initial.clone();
                let term = rewrite_ports(term, slots, replacements, &mut local);
                if let Some(common) = &mut common {
                    for (common, local) in common.claimed.iter_mut().zip(&local.claimed) {
                        *common &= local;
                    }
                    for (common, local) in common.applied.iter_mut().zip(&local.applied) {
                        *common &= local;
                    }
                } else {
                    common = Some(local);
                }
                sum + term
            });
            *state = common.unwrap_or(initial);
            sum
        }
        AtomView::Mul(mul) => mul.iter().fold(Atom::num(1), |product, factor| {
            let tensorial = factor.is_tensorial(StrictTensorFilter::Tagged)
                || matches!(
                    factor,
                    AtomView::Fun(fun)
                        if fun.get_symbol() == *shadowing::SYM
                            || fun.get_symbol() == *shadowing::ANTISYM
                            || fun.get_symbol() == *shadowing::CYCLIC
                );
            product
                * if tensorial {
                    rewrite_ports(factor, slots, replacements, state)
                } else {
                    factor.to_owned()
                }
        }),
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            rewrite_ports(base, slots, replacements, state).pow(exponent.to_owned())
        }
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace => {
            let mut args = fun.iter();
            let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
            if let Some(rep) = args.next() {
                rebuilt = rebuilt.add_arg(rep);
            }
            for arg in args {
                rebuilt = rebuilt.add_arg(rewrite_ports(arg, slots, replacements, state));
            }
            rebuilt.finish()
        }
        // A dot's rank-one channel is internal and therefore does not consume
        // positions in the dot expression's external interface.
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot => value.to_owned(),
        AtomView::Fun(fun) => {
            let tensor_leaf = is_tensor_leaf_head(fun.get_symbol());
            let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
            for arg in fun.iter() {
                if direct_structural_port(arg) {
                    if let Some(position) = matching_interface_position(arg, slots, &state.claimed)
                    {
                        state.claimed[position] = true;
                        if let Some(replacement) = replacements.get(&position) {
                            state.applied[position] = true;
                            rebuilt = rebuilt.add_arg(replacement);
                        } else {
                            rebuilt = rebuilt.add_arg(arg);
                        }
                    } else {
                        // Contracted dummy slots and representation-valued metadata
                        // are not part of the surviving public interface.
                        rebuilt = rebuilt.add_arg(arg);
                    }
                } else if tensor_leaf {
                    // Scalar metadata belongs to this tensor leaf. Nested
                    // representations inside it are not public tensor ports.
                    rebuilt = rebuilt.add_arg(arg);
                } else if is_composite_head(fun.get_symbol())
                    && !arg.is_tensorial(StrictTensorFilter::Tagged)
                    && !matches!(
                        arg,
                        AtomView::Fun(nested)
                            if nested.get_symbol() == *shadowing::SYM
                                || nested.get_symbol() == *shadowing::ANTISYM
                                || nested.get_symbol() == *shadowing::CYCLIC
                    )
                {
                    rebuilt = rebuilt.add_arg(arg);
                } else {
                    rebuilt = rebuilt.add_arg(rewrite_ports(arg, slots, replacements, state));
                }
            }
            rebuilt.finish()
        }
        _ => value.to_owned(),
    }
}

fn rewrite_interface_ports(
    value: &StructuredAtom,
    replacements: &HashMap<usize, Atom>,
) -> Result<Atom, TensorCompositionError> {
    if value.atom.as_view().is_zero() {
        return Ok(Atom::Zero);
    }

    let slots = value.interface.logical_slots();
    let mut state = PortRewriteState::new(slots.len());
    let atom = rewrite_ports(value.atom.as_view(), &slots, replacements, &mut state);
    if let Some(position) = replacements
        .keys()
        .copied()
        .find(|position| !state.applied[*position])
    {
        return Err(TensorCompositionError::MissingInterfacePort { position });
    }
    Ok(atom)
}

fn collect_interface_positions(
    value: AtomView<'_>,
    slots: &[PartialSlot],
    claimed: &mut [bool],
    positions: &mut Vec<usize>,
) {
    if direct_structural_port(value) {
        if let Some(position) = matching_interface_position(value, slots, claimed) {
            claimed[position] = true;
            positions.push(position);
        }
        return;
    }

    match value {
        // Every summand has the same public interface. Inspecting the first one
        // avoids collecting the same external ports repeatedly.
        AtomView::Add(add) => {
            if let Some(term) = add.iter().next() {
                collect_interface_positions(term, slots, claimed, positions);
            }
        }
        AtomView::Mul(mul) => {
            for factor in mul.iter() {
                let tensorial = factor.is_tensorial(StrictTensorFilter::Tagged)
                    || matches!(
                        factor,
                        AtomView::Fun(fun)
                            if fun.get_symbol() == *shadowing::SYM
                                || fun.get_symbol() == *shadowing::ANTISYM
                                || fun.get_symbol() == *shadowing::CYCLIC
                    );
                if tensorial {
                    collect_interface_positions(factor, slots, claimed, positions);
                }
            }
        }
        AtomView::Pow(pow) => {
            let (base, _) = pow.get_base_exp();
            collect_interface_positions(base, slots, claimed, positions);
        }
        // A compact dot has no public ports of its own.
        AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot => {}
        AtomView::Fun(fun) => {
            let symbol = fun.get_symbol();
            let tensor_leaf = is_tensor_leaf_head(symbol);
            let skip = usize::from(symbol == SPENSO_TAG.trace);
            for argument in fun.iter().skip(skip) {
                let tensor_argument = !tensor_leaf
                    && is_composite_head(symbol)
                    && (argument.is_tensorial(StrictTensorFilter::Tagged)
                        || matches!(
                            argument,
                            AtomView::Fun(nested)
                                if nested.get_symbol() == *shadowing::SYM
                                    || nested.get_symbol() == *shadowing::ANTISYM
                                    || nested.get_symbol() == *shadowing::CYCLIC
                        ));
                if direct_structural_port(argument) || tensor_argument {
                    collect_interface_positions(argument, slots, claimed, positions);
                }
            }
        }
        _ => {}
    }
}

/// Align partial-interface metadata with the canonical factor order stored in
/// an atom. Cyclic trace normalization may rotate factors, so retaining the
/// pre-normalization interface order would make `as_tensor(to_expression())`
/// observe a different logical interface.
fn interface_in_atom_order(
    atom: &Atom,
    interface: &PartialStructure,
) -> Result<PartialStructure, TensorCompositionError> {
    let slots = interface.logical_slots();
    let mut claimed = vec![false; slots.len()];
    let mut positions = Vec::with_capacity(slots.len());
    collect_interface_positions(atom.as_view(), &slots, &mut claimed, &mut positions);
    if let Some(position) = claimed.iter().position(|claimed| !claimed) {
        return Err(TensorCompositionError::MissingInterfacePort { position });
    }
    Ok(PartialStructure::from_logical_slots(
        positions.into_iter().map(|position| slots[position]),
    ))
}

/// Replace selected unresolved ports by explicit indices in every summand.
///
/// Positions are interpreted in the public logical interface, never in the
/// canonical storage order of `OrderedStructure`.
pub(crate) fn materialize_interface_ports(
    value: &StructuredAtom,
    replacements: &HashMap<usize, AbstractIndex>,
) -> Result<Atom, TensorCompositionError> {
    let slots = value.interface.logical_slots();
    let replacements = replacements
        .iter()
        .map(|(&position, &index)| {
            let slot = validate_position(&slots, position)?;
            let atom = match slot.aind {
                PartialIndex::Explicit(current) if current == index => port_atom(slot),
                PartialIndex::Explicit(_) => {
                    return Err(TensorCompositionError::UnequalExplicitIndices {
                        left: position,
                        right: position,
                    });
                }
                PartialIndex::Open(_) => slot.rep().slot::<AbstractIndex, _>(index).to_atom(),
            };
            Ok((position, atom))
        })
        .collect::<Result<HashMap<_, _>, TensorCompositionError>>()?;
    let atom = rewrite_interface_ports(value, &replacements)?;
    validate_explicit_index_occurrences(&atom)?;
    Ok(atom)
}

/// Replace selected public ports by explicit graph indices.
///
/// Unlike `materialize_interface_ports`, this also permits changing an
/// already-explicit index. It is used to keep a tensor network's execution
/// graph collision-free while its semantic interface remains unchanged.
pub(crate) fn reindex_interface_ports(
    value: &StructuredAtom,
    replacements: &HashMap<usize, AbstractIndex>,
) -> Result<StructuredAtom, TensorCompositionError> {
    let mut logical = value.interface.logical_slots();
    let replacements = replacements
        .iter()
        .map(|(&position, &index)| {
            let slot = validate_position(&logical, position)?;
            logical[position].set_aind(PartialIndex::Explicit(index));
            Ok((
                position,
                slot.rep().slot::<AbstractIndex, _>(index).to_atom(),
            ))
        })
        .collect::<Result<HashMap<_, _>, TensorCompositionError>>()?;
    let atom = rewrite_interface_ports(value, &replacements)?;
    validate_explicit_index_occurrences(&atom)?;
    Ok(StructuredAtom::new(
        atom,
        PartialStructure::from_logical_slots(logical),
    ))
}

fn without_positions(interface: &PartialStructure, positions: &[usize]) -> PartialStructure {
    PartialStructure::from_logical_slots(
        interface
            .logical_slots()
            .into_iter()
            .enumerate()
            .filter(|(position, _)| !positions.contains(position))
            .map(|(_, slot)| slot),
    )
}

/// Canonicalize a root chain whose two explicit endpoints have become one contraction.
pub(crate) fn normalize_closed_root_chain(
    value: StructuredAtom,
) -> Result<StructuredAtom, TensorCompositionError> {
    let AtomView::Fun(function) = value.atom.as_view() else {
        return Ok(value);
    };
    if function.get_symbol() != SPENSO_TAG.chain {
        return Ok(value);
    }
    let arguments = function.iter().collect::<Vec<_>>();
    let [start, end, factors @ ..] = arguments.as_slice() else {
        return Ok(value);
    };
    let (Ok(start), Ok(end)) = (
        Slot::<LibraryRep, AbstractIndex>::try_from(*start),
        Slot::<LibraryRep, AbstractIndex>::try_from(*end),
    ) else {
        return Ok(value);
    };
    let start_rep = start.rep();
    let end_rep = end.rep();
    if start.aind() != end.aind()
        || !start_rep.matches(&end_rep)
        || (!start_rep.rep.is_self_dual() && !(start_rep.rep.is_base() && end_rep.rep.is_dual()))
    {
        return Ok(value);
    }

    let channel = MatrixChannel {
        input: 0,
        output: 1,
    };
    let interface = if root_chain_channel_is_live(&value, channel) {
        without_positions(&value.interface, &[channel.input, channel.output])
    } else {
        value.interface.clone()
    };
    let atom = shadowing::trace(
        start_rep.to_symbolic([]),
        factors.iter().map(|factor| factor.to_owned()),
    );
    let interface = interface_in_atom_order(&atom, &interface)?;
    Ok(StructuredAtom::new(atom, interface))
}

fn concatenate_interfaces(
    left: impl IntoIterator<Item = PartialSlot>,
    right: impl IntoIterator<Item = PartialSlot>,
) -> PartialStructure {
    PartialStructure::from_logical_slots(left.into_iter().chain(right))
}

pub fn outer(left: &StructuredAtom, right: &StructuredAtom) -> StructuredAtom {
    let atom = if left.atom.as_view().is_zero() || right.atom.as_view().is_zero() {
        Atom::Zero
    } else {
        FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(&left.atom)
            .add_arg(&right.atom)
            .finish()
    };
    StructuredAtom::new(
        atom,
        concatenate_interfaces(
            left.interface.logical_slots(),
            right.interface.logical_slots(),
        ),
    )
}

pub fn contract(
    left: &StructuredAtom,
    right: &StructuredAtom,
    positions: PortPair,
) -> Result<StructuredAtom, TensorCompositionError> {
    let left_slots = left.interface.logical_slots();
    let right_slots = right.interface.logical_slots();
    let left_slot = validate_position(&left_slots, positions.left)?;
    let right_slot = validate_position(&right_slots, positions.right)?;
    let shared = shared_index(left_slot, right_slot, positions)?;
    let spectator_pairs = left_slots
        .iter()
        .enumerate()
        .filter(|(position, _)| *position != positions.left)
        .flat_map(|(left_position, left_slot)| {
            right_slots
                .iter()
                .enumerate()
                .filter(move |(position, _)| *position != positions.right)
                .filter_map(move |(right_position, right_slot)| {
                    matches!(
                        (left_slot.aind, right_slot.aind),
                        (PartialIndex::Explicit(left), PartialIndex::Explicit(right))
                            if left == right && representations_match(left_slot, right_slot)
                    )
                    .then_some(PortPair {
                        left: left_position,
                        right: right_position,
                    })
                })
        })
        .collect::<Vec<_>>();
    if spectator_pairs.iter().enumerate().any(|(position, pair)| {
        spectator_pairs[position + 1..]
            .iter()
            .any(|candidate| pair.left == candidate.left || pair.right == candidate.right)
    }) {
        return Err(TensorCompositionError::Ambiguous(spectator_pairs));
    }
    let contracted_left = spectator_pairs
        .iter()
        .map(|pair| pair.left)
        .chain([positions.left])
        .collect::<HashSet<_>>();
    let contracted_right = spectator_pairs
        .iter()
        .map(|pair| pair.right)
        .chain([positions.right])
        .collect::<HashSet<_>>();
    let interface = PartialStructure::from_logical_slots(
        left_slots
            .into_iter()
            .enumerate()
            .filter(|(position, _)| !contracted_left.contains(position))
            .map(|(_, slot)| slot)
            .chain(
                right_slots
                    .into_iter()
                    .enumerate()
                    .filter(|(position, _)| !contracted_right.contains(position))
                    .map(|(_, slot)| slot),
            ),
    );
    if left.atom.as_view().is_zero() || right.atom.as_view().is_zero() {
        return Ok(StructuredAtom::new(Atom::Zero, interface));
    }

    let index = shared.unwrap_or_else(|| {
        fresh_dummy_index(
            [&left.atom, &right.atom],
            [&left.interface, &right.interface],
        )
    });
    let materialized_left =
        materialize_interface_ports(left, &HashMap::from([(positions.left, index)]))?;
    let materialized_right =
        materialize_interface_ports(right, &HashMap::from([(positions.right, index)]))?;

    let atom = if left.rank() == 1 && right.rank() == 1 {
        // Preserve the established compact dot spelling; the parser performs
        // the same shared-dummy materialization represented above.
        let compact_left = rewrite_interface_ports(
            left,
            &HashMap::from([(positions.left, left_slot.rep().to_symbolic([]))]),
        )?;
        let compact_right = rewrite_interface_ports(
            right,
            &HashMap::from([(positions.right, right_slot.rep().to_symbolic([]))]),
        )?;
        FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(compact_left)
            .add_arg(compact_right)
            .finish()
    } else {
        FunctionBuilder::new(SPENSO_TAG.bracket)
            .add_arg(materialized_left)
            .add_arg(materialized_right)
            .finish()
    };
    validate_explicit_index_occurrences(&atom)?;
    Ok(StructuredAtom::new(atom, interface))
}

pub(crate) fn port_atom(slot: PartialSlot) -> Atom {
    match slot.aind {
        PartialIndex::Explicit(index) => slot.rep().slot::<AbstractIndex, _>(index).to_atom(),
        PartialIndex::Open(_) => slot.rep().to_symbolic([]),
    }
}

fn channel_factor(
    value: &StructuredAtom,
    channel: MatrixChannel,
) -> Result<Atom, TensorCompositionError> {
    let replacements = HashMap::from([
        (channel.input, Atom::var(SPENSO_TAG.chain_in)),
        (channel.output, Atom::var(SPENSO_TAG.chain_out)),
    ]);
    rewrite_interface_ports(value, &replacements)
}

pub(crate) fn chain_factors(
    value: &StructuredAtom,
    channel: MatrixChannel,
) -> Result<Vec<Atom>, TensorCompositionError> {
    value.atom.validate_chain_like_nesting()?;
    if root_chain_channel_is_live(value, channel) {
        let AtomView::Fun(function) = value.atom.as_view() else {
            unreachable!("a live root chain channel requires a chain function")
        };
        return Ok(function
            .iter()
            .skip(2)
            .map(|factor| factor.to_owned())
            .collect());
    }
    if contains_chain_like(value.atom.as_view()) {
        return Err(TensorCompositionError::NestedChainLike);
    }
    Ok(vec![channel_factor(value, channel)?])
}

pub fn compose(
    left: &StructuredAtom,
    right: &StructuredAtom,
    left_channel: MatrixChannel,
    right_channel: MatrixChannel,
) -> Result<StructuredAtom, TensorCompositionError> {
    if left_channel.input == left_channel.output {
        return Err(TensorCompositionError::DegenerateChannel {
            input: left_channel.input,
            output: left_channel.output,
        });
    }
    if right_channel.input == right_channel.output {
        return Err(TensorCompositionError::DegenerateChannel {
            input: right_channel.input,
            output: right_channel.output,
        });
    }
    let left_slots = left.interface.logical_slots();
    let right_slots = right.interface.logical_slots();
    let left_input = validate_position(&left_slots, left_channel.input)?;
    let left_output = validate_position(&left_slots, left_channel.output)?;
    let right_input = validate_position(&right_slots, right_channel.input)?;
    let right_output = validate_position(&right_slots, right_channel.output)?;
    if !representations_match(&left_input, &left_output) {
        return Err(TensorCompositionError::IncompatiblePorts {
            left: left_channel.input,
            right: left_channel.output,
        });
    }
    let left_input_rep = representation(&left_input);
    let left_output_rep = representation(&left_output);
    if !left_input_rep.rep.is_self_dual()
        && !(left_input_rep.rep.is_base() && left_output_rep.rep.is_dual())
    {
        return Err(TensorCompositionError::InvalidChannelOrientation {
            input: left_channel.input,
            output: left_channel.output,
        });
    }
    if !representations_match(&right_input, &right_output) {
        return Err(TensorCompositionError::IncompatiblePorts {
            left: right_channel.input,
            right: right_channel.output,
        });
    }
    let right_input_rep = representation(&right_input);
    let right_output_rep = representation(&right_output);
    if !right_input_rep.rep.is_self_dual()
        && !(right_input_rep.rep.is_base() && right_output_rep.rep.is_dual())
    {
        return Err(TensorCompositionError::InvalidChannelOrientation {
            input: right_channel.input,
            output: right_channel.output,
        });
    }
    let shared = shared_index(
        left_output,
        right_input,
        PortPair {
            left: left_channel.output,
            right: right_channel.input,
        },
    )?;
    if let Some(index) = shared {
        let target = left_output.rep().slot(index);
        let occurrences = explicit_occurrences(left.atom.as_view(), target)
            .saturating_add(explicit_occurrences(right.atom.as_view(), target))
            .saturating_add(usize::from(matches!(
                left_output.aind,
                PartialIndex::Open(_)
            )))
            .saturating_add(usize::from(matches!(
                right_input.aind,
                PartialIndex::Open(_)
            )));
        if occurrences > 2 {
            return Err(TensorCompositionError::InvalidExplicitMultiplicity { index, occurrences });
        }
    }

    let spectator_pairs = left_slots
        .iter()
        .enumerate()
        .filter(|(position, _)| *position != left_channel.input && *position != left_channel.output)
        .flat_map(|(left_position, left_slot)| {
            right_slots
                .iter()
                .enumerate()
                .filter(move |(position, _)| {
                    *position != right_channel.input && *position != right_channel.output
                })
                .filter_map(move |(right_position, right_slot)| {
                    matches!(
                        (left_slot.aind, right_slot.aind),
                        (PartialIndex::Explicit(left), PartialIndex::Explicit(right))
                            if left == right && representations_match(left_slot, right_slot)
                    )
                    .then_some(PortPair {
                        left: left_position,
                        right: right_position,
                    })
                })
        })
        .collect::<Vec<_>>();
    if spectator_pairs.iter().enumerate().any(|(position, pair)| {
        spectator_pairs[position + 1..]
            .iter()
            .any(|candidate| pair.left == candidate.left || pair.right == candidate.right)
    }) {
        return Err(TensorCompositionError::Ambiguous(spectator_pairs));
    }
    let contracted_left = spectator_pairs
        .iter()
        .map(|pair| pair.left)
        .collect::<HashSet<_>>();
    let contracted_right = spectator_pairs
        .iter()
        .map(|pair| pair.right)
        .collect::<HashSet<_>>();

    let interface = PartialStructure::from_logical_slots(
        [left_input, right_output]
            .into_iter()
            .chain(
                left_slots
                    .into_iter()
                    .enumerate()
                    .filter(|(position, _)| {
                        *position != left_channel.input
                            && *position != left_channel.output
                            && !contracted_left.contains(position)
                    })
                    .map(|(_, slot)| slot),
            )
            .chain(
                right_slots
                    .into_iter()
                    .enumerate()
                    .filter(|(position, _)| {
                        *position != right_channel.input
                            && *position != right_channel.output
                            && !contracted_right.contains(position)
                    })
                    .map(|(_, slot)| slot),
            ),
    );
    if left.atom.as_view().is_zero() || right.atom.as_view().is_zero() {
        return Ok(StructuredAtom::new(Atom::Zero, interface));
    }

    let factors = chain_factors(left, left_channel)?
        .into_iter()
        .chain(chain_factors(right, right_channel)?);
    let atom = SPENSO_TAG.chain(port_atom(left_input), port_atom(right_output), factors);
    validate_explicit_index_occurrences(&atom)?;
    normalize_closed_root_chain(StructuredAtom::new(atom, interface))
}

pub fn multiply(
    left: &StructuredAtom,
    right: &StructuredAtom,
) -> Result<StructuredAtom, TensorCompositionError> {
    if left.is_scalar() || right.is_scalar() {
        let interface = if left.is_scalar() {
            right.interface.clone()
        } else {
            left.interface.clone()
        };
        return Ok(StructuredAtom::new(
            left.atom.as_ref() * right.atom.as_ref(),
            interface,
        ));
    }

    if let (Some(left_channel), Some(right_channel)) = (matrix_channel(left), matrix_channel(right))
    {
        let left_slots = left.interface.logical_slots();
        let right_slots = right.interface.logical_slots();
        if automatically_contractible(
            &left_slots[left_channel.output],
            &right_slots[right_channel.input],
        ) {
            return compose(left, right, left_channel, right_channel);
        }
    }

    match compatible_pairs(&left.interface, &right.interface).as_slice() {
        [] => Ok(outer(left, right)),
        [positions] => contract(left, right, *positions),
        candidates => Err(TensorCompositionError::Ambiguous(candidates.to_vec())),
    }
}

pub fn trace(
    value: &StructuredAtom,
    channel: MatrixChannel,
) -> Result<StructuredAtom, TensorCompositionError> {
    value.atom.validate_chain_like_nesting()?;
    if channel.input == channel.output {
        return Err(TensorCompositionError::DegenerateChannel {
            input: channel.input,
            output: channel.output,
        });
    }
    let slots = value.interface.logical_slots();
    let input = validate_position(&slots, channel.input)?;
    let output = validate_position(&slots, channel.output)?;
    if !representations_match(&input, &output) {
        return Err(TensorCompositionError::IncompatiblePorts {
            left: channel.input,
            right: channel.output,
        });
    }
    let input_rep = representation(&input);
    let output_rep = representation(&output);
    if !input_rep.rep.is_self_dual() && !(input_rep.rep.is_base() && output_rep.rep.is_dual()) {
        return Err(TensorCompositionError::InvalidChannelOrientation {
            input: channel.input,
            output: channel.output,
        });
    }
    let shared = shared_index(
        input,
        output,
        PortPair {
            left: channel.input,
            right: channel.output,
        },
    )?;
    let interface = without_positions(&value.interface, &[channel.input, channel.output]);
    if value.atom.as_view().is_zero() {
        return Ok(StructuredAtom::new(Atom::Zero, interface));
    }

    if !root_chain_channel_is_live(value, channel)
        && matches!(value.atom.as_view(), AtomView::Fun(function)
            if function.get_symbol() == SPENSO_TAG.chain
                || function.get_symbol() == SPENSO_TAG.trace)
    {
        // The root shorthand already owns the global `in`/`out` placeholders.
        // Close a spectator pair in place instead of adding a second consumer.
        let index = shared.unwrap_or_else(|| fresh_dummy_index([&value.atom], [&value.interface]));
        let atom = materialize_interface_ports(
            value,
            &HashMap::from([(channel.input, index), (channel.output, index)]),
        )?;
        let interface = interface_in_atom_order(&atom, &interface)?;
        return Ok(StructuredAtom::new(atom, interface));
    }
    if contains_chain_like(value.atom.as_view()) && !root_chain_channel_is_live(value, channel) {
        return Err(TensorCompositionError::NestedChainLike);
    }

    let factors = chain_factors(value, channel)?;
    let atom = shadowing::trace(input.rep().to_symbolic([]), factors);
    let interface = interface_in_atom_order(&atom, &interface)?;
    Ok(StructuredAtom::new(atom, interface))
}

pub fn trace_unique(value: &StructuredAtom) -> Result<StructuredAtom, TensorCompositionError> {
    let channel = matrix_channel(value).ok_or(TensorCompositionError::NoMatrixChannel)?;
    trace(value, channel)
}

#[cfg(test)]
mod tests {
    use super::*;
    use idenso::{
        color::CS,
        dirac::AGS,
        representations::{Bispinor, ColorAdjoint, ColorFundamental},
    };
    use spenso::network::parsing::{NetworkParse, ParseSettings};
    use spenso::structure::{
        dimension::Dimension,
        partial::PartialIndex,
        representation::{ExtendibleReps, Minkowski, RepName},
    };
    use symbolica::atom::Symbol;

    fn rep() -> Representation<LibraryRep> {
        ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(4))
    }

    fn tensor(name: &str, reps: &[Representation<LibraryRep>]) -> StructuredAtom {
        let symbol: Symbol = SPENSO_TAG.tensor_symbol(name);
        let ports = reps
            .iter()
            .enumerate()
            .map(|(position, rep)| rep.slot(PartialIndex::open(position)))
            .collect::<Vec<_>>();
        partial_tensor(symbol, &ports, &ports)
    }

    fn partial_tensor(
        symbol: Symbol,
        atom_ports: &[PartialSlot],
        logical_ports: &[PartialSlot],
    ) -> StructuredAtom {
        let atom = atom_ports
            .iter()
            .fold(FunctionBuilder::new(symbol), |builder, &slot| {
                builder.add_arg(port_atom(slot))
            })
            .finish();
        StructuredAtom::new(
            atom,
            PartialStructure::from_logical_slots(logical_ports.iter().copied()),
        )
    }

    fn factor_slots(atom: &Atom, symbol: Symbol) -> Vec<Slot<LibraryRep, AbstractIndex>> {
        let AtomView::Fun(product) = atom.as_view() else {
            panic!("expected an ordered product")
        };
        product
            .iter()
            .find_map(|factor| {
                let AtomView::Fun(fun) = factor else {
                    return None;
                };
                (fun.get_symbol() == symbol).then(|| {
                    fun.iter()
                        .filter_map(|argument| {
                            Slot::<LibraryRep, AbstractIndex>::try_from(argument).ok()
                        })
                        .collect()
                })
            })
            .unwrap()
    }

    #[test]
    fn rank_one_product_uses_dot() {
        let left = tensor("left", &[rep()]);
        let right = tensor("right", &[rep()]);
        let result = multiply(&left, &right).unwrap();

        assert!(result.is_scalar());
        assert!(
            matches!(result.atom.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot)
        );
    }

    #[test]
    fn explicit_rank_one_contraction_keeps_compact_dot_syntax() {
        let index = AbstractIndex::Normal(7);
        let port = rep().slot(PartialIndex::Explicit(index));
        let left_symbol = SPENSO_TAG.tensor_symbol("explicit_dot_left");
        let right_symbol = SPENSO_TAG.tensor_symbol("explicit_dot_right");
        let left = partial_tensor(left_symbol, &[port], &[port]);
        let right = partial_tensor(right_symbol, &[port], &[port]);
        let result = contract(&left, &right, PortPair { left: 0, right: 0 }).unwrap();

        let AtomView::Fun(dot) = result.atom.as_view() else {
            panic!("expected a dot")
        };
        assert_eq!(dot.get_symbol(), SPENSO_TAG.dot);
        for operand in dot.iter() {
            let AtomView::Fun(vector) = operand else {
                panic!("expected a compact vector")
            };
            let argument = vector.iter().next().unwrap();
            assert!(Representation::<LibraryRep>::try_from(argument).is_ok());
            assert!(Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err());
        }

        let network = result
            .atom
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.state.is_scalar());
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn dual_rank_one_dot_materializes_with_opposite_orientations() {
        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let dual = base.dual();
        let left = tensor("dual_dot_left", &[base]);
        let right = tensor("dual_dot_right", &[dual]);
        let result = multiply(&left, &right).unwrap();

        assert!(
            matches!(result.atom.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.dot)
        );
        let network = result
            .atom
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(
            network.state.is_scalar(),
            "state {:?}, dangling {:?}",
            network.state,
            network.graph.dangling_indices()
        );
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn matrix_products_flatten_chains() {
        let first = tensor("first", &[rep(), rep()]);
        let second = tensor("second", &[rep(), rep()]);
        let third = tensor("third", &[rep(), rep()]);
        let chain = multiply(&first, &second).unwrap();
        let chain = multiply(&chain, &third).unwrap();

        let AtomView::Fun(fun) = chain.atom.as_view() else {
            panic!("expected a chain")
        };
        assert_eq!(fun.get_symbol(), SPENSO_TAG.chain);
        assert_eq!(fun.get_nargs(), 5);
        assert_eq!(chain.rank(), 2);
    }

    #[test]
    fn transparent_chain_wrappers_expose_but_cannot_extend_the_chain_channel() {
        let first = tensor("wrapped_chain_first", &[rep(), rep()]);
        let second = tensor("wrapped_chain_second", &[rep(), rep()]);
        let third = tensor("wrapped_chain_third", &[rep(), rep()]);
        let chain = multiply(&first, &second).unwrap();
        let scalar = Atom::var(symbolica::symbol!("wrapped_chain_scalar"));
        let broadcast = FunctionBuilder::new(spenso::broadcast_symbol!("wrapped_chain_broadcast"))
            .add_arg(&chain.atom)
            .finish();

        for atom in [scalar * chain.atom.as_ref(), broadcast] {
            let wrapped = StructuredAtom::new(atom, chain.interface.clone());
            assert_eq!(
                matrix_channel(&wrapped),
                Some(MatrixChannel {
                    input: 0,
                    output: 1,
                })
            );

            assert!(matches!(
                multiply(&wrapped, &third),
                Err(TensorCompositionError::NestedChainLike)
            ));
            assert!(matches!(
                trace(
                    &wrapped,
                    MatrixChannel {
                        input: 0,
                        output: 1,
                    },
                ),
                Err(TensorCompositionError::NestedChainLike)
            ));
        }
    }

    #[test]
    fn contracted_chain_endpoints_are_not_treated_as_the_public_channel() {
        let channel = rep();
        let endpoint = channel
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(61))
            .to_atom();
        let spectators = [
            ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4)),
            ColorAdjoint {}.new_rep(8).cast(),
            Bispinor {}.new_rep(4).cast(),
        ];
        let factor = spectators.iter().fold(
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol("closed_chain_spectators"))
                .add_arg(Atom::var(SPENSO_TAG.chain_in))
                .add_arg(Atom::var(SPENSO_TAG.chain_out)),
            |builder, representation| builder.add_arg(representation.to_symbolic([])),
        );
        let value = StructuredAtom::new(
            SPENSO_TAG.chain(&endpoint, &endpoint, [factor.finish()]),
            PartialStructure::from_logical_slots(spectators.into_iter().enumerate().map(
                |(position, representation)| representation.slot(PartialIndex::open(position)),
            )),
        );

        assert_eq!(matrix_channel(&value), None);
    }

    #[test]
    fn matrix_composition_contracts_equal_explicit_spectators() {
        let channel = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(47);
        let left_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let right_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_spectator_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_spectator_right"),
            &right_ports,
            &right_ports,
        );

        let result = multiply(&left, &right).unwrap();

        assert_eq!(result.rank(), 2);
        assert!(
            matches!(result.atom.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.chain)
        );
        let materialized = materialize_interface_ports(
            &result,
            &HashMap::from([
                (0, AbstractIndex::Normal(61)),
                (1, AbstractIndex::Normal(67)),
            ]),
        )
        .unwrap();
        let network = materialized
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn explicit_contraction_contracts_equal_explicit_spectators() {
        let channel = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(51);
        let left_ports = [
            channel.slot(PartialIndex::open(0)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let right_ports = [
            channel.slot(PartialIndex::open(0)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_contract_spectator_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("explicit_contract_spectator_right"),
            &right_ports,
            &right_ports,
        );

        let result = contract(&left, &right, PortPair { left: 0, right: 0 }).unwrap();

        assert!(result.is_scalar());
        let network = result
            .atom
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn matrix_composition_rejects_ambiguous_explicit_spectators() {
        let channel = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let index = AbstractIndex::Normal(53);
        let left_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let right_ports = [
            channel.slot(PartialIndex::open(0)),
            channel.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::Explicit(index)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("ambiguous_spectator_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("ambiguous_spectator_right"),
            &right_ports,
            &right_ports,
        );

        assert!(matches!(
            compose(
                &left,
                &right,
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
            ),
            Err(TensorCompositionError::Ambiguous(candidates)) if candidates.len() == 2
        ));
    }

    #[test]
    fn multiple_compatible_spectators_are_ambiguous() {
        let left = tensor("left", &[rep(), rep(), rep()]);
        let right = tensor("right", &[rep()]);
        let error = multiply(&left, &right).unwrap_err();

        assert!(
            matches!(error, TensorCompositionError::Ambiguous(candidates) if candidates.len() == 3)
        );
    }

    #[test]
    fn outer_materialization_uses_logical_operand_order() {
        let q = tensor("z_outer_q", &[rep()]);
        let p = tensor("a_outer_p", &[rep()]);
        let product = outer(&q, &p);
        let indexed =
            materialize_interface_ports(&product, &HashMap::from([(0, AbstractIndex::Normal(17))]))
                .unwrap();
        let q_symbol = SPENSO_TAG.tensor_symbol("z_outer_q");
        let p_symbol = SPENSO_TAG.tensor_symbol("a_outer_p");
        let AtomView::Fun(product) = indexed.as_view() else {
            panic!("expected an ordered product")
        };
        assert_eq!(product.get_symbol(), SPENSO_TAG.bracket);
        let mut q_is_indexed = false;
        let mut p_is_open = false;
        for factor in product.iter() {
            let AtomView::Fun(fun) = factor else {
                continue;
            };
            let argument = fun.iter().next().unwrap();
            if fun.get_symbol() == q_symbol {
                q_is_indexed = Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_ok();
            } else if fun.get_symbol() == p_symbol {
                p_is_open = Representation::<LibraryRep>::try_from(argument).is_ok()
                    && Slot::<LibraryRep, AbstractIndex>::try_from(argument).is_err();
            }
        }
        assert!(q_is_indexed && p_is_open);
    }

    #[test]
    fn materialization_matches_logical_ports_in_canonical_atom_order() {
        let euc = rep();
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let euc_port = euc.slot(PartialIndex::open(0));
        let mink_port = mink.slot(PartialIndex::open(1));
        let value = partial_tensor(
            SPENSO_TAG.tensor_symbol("canonical_port_order"),
            &[euc_port, mink_port],
            &[mink_port, euc_port],
        );
        let indexed =
            materialize_interface_ports(&value, &HashMap::from([(0, AbstractIndex::Normal(19))]))
                .unwrap();

        let AtomView::Fun(tensor) = indexed.as_view() else {
            panic!("expected a tensor")
        };
        let args = tensor.iter().collect::<Vec<_>>();
        assert_eq!(
            Representation::<LibraryRep>::try_from(args[0]).unwrap(),
            euc
        );
        assert_eq!(
            Slot::<LibraryRep, AbstractIndex>::try_from(args[1]).unwrap(),
            mink.slot(AbstractIndex::Normal(19))
        );
    }

    #[test]
    fn materialization_does_not_rewrite_tensor_scalar_metadata() {
        let representation = rep();
        let port = representation.slot(PartialIndex::open(0));
        let metadata = FunctionBuilder::new(symbolica::symbol!("composition_scalar_metadata"))
            .add_arg(representation.to_symbolic([]))
            .finish();
        let atom = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("tensor_with_scalar_metadata"))
            .add_arg(metadata)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let value = StructuredAtom::new(atom, PartialStructure::from_logical_slots([port]));

        let indexed =
            materialize_interface_ports(&value, &HashMap::from([(0, AbstractIndex::Normal(21))]))
                .unwrap();
        let AtomView::Fun(tensor) = indexed.as_view() else {
            panic!("expected a tensor")
        };
        let arguments = tensor.iter().collect::<Vec<_>>();
        let AtomView::Fun(metadata) = arguments[0] else {
            panic!("expected scalar metadata")
        };

        assert!(Representation::<LibraryRep>::try_from(metadata.iter().next().unwrap()).is_ok());
        assert_eq!(
            Slot::<LibraryRep, AbstractIndex>::try_from(arguments[1]).unwrap(),
            representation.slot(AbstractIndex::Normal(21))
        );
    }

    #[test]
    fn materialization_does_not_rewrite_scalar_factor_metadata() {
        let representation = rep();
        let port = representation.slot(PartialIndex::open(0));
        let metadata_symbol = symbolica::symbol!("composition_external_scalar_metadata");
        let tensor_symbol = SPENSO_TAG.tensor_symbol("tensor_beside_scalar_metadata");
        let metadata = FunctionBuilder::new(metadata_symbol)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let tensor = FunctionBuilder::new(tensor_symbol)
            .add_arg(representation.to_symbolic([]))
            .finish();
        let value = StructuredAtom::new(
            metadata * tensor,
            PartialStructure::from_logical_slots([port]),
        );

        let indexed =
            materialize_interface_ports(&value, &HashMap::from([(0, AbstractIndex::Normal(22))]))
                .unwrap();
        let AtomView::Mul(product) = indexed.as_view() else {
            panic!("expected a product")
        };
        let arguments = product
            .iter()
            .filter_map(|factor| {
                let AtomView::Fun(function) = factor else {
                    return None;
                };
                Some((function.get_symbol(), function.iter().next().unwrap()))
            })
            .collect::<HashMap<_, _>>();

        assert!(Representation::<LibraryRep>::try_from(arguments[&metadata_symbol]).is_ok());
        assert_eq!(
            Slot::<LibraryRep, AbstractIndex>::try_from(arguments[&tensor_symbol]).unwrap(),
            representation.slot(AbstractIndex::Normal(22))
        );
    }

    #[test]
    fn structured_zero_preserves_interface_metadata_through_tensor_operations() {
        let representation = rep();
        let vector_port = representation.slot(PartialIndex::open(0));
        let zero_vector = StructuredAtom::new(
            Atom::Zero,
            PartialStructure::from_logical_slots([vector_port]),
        );
        let vector = tensor("structured_zero_vector_partner", &[representation]);

        assert!(
            materialize_interface_ports(
                &zero_vector,
                &HashMap::from([(0, AbstractIndex::Normal(24))]),
            )
            .unwrap()
            .as_view()
            .is_zero()
        );
        let dotted = contract(&zero_vector, &vector, PortPair { left: 0, right: 0 }).unwrap();
        assert!(dotted.is_scalar());
        assert!(dotted.atom.as_view().is_zero());

        let outer_product = outer(&zero_vector, &vector);
        assert_eq!(outer_product.rank(), 2);
        assert!(outer_product.atom.as_view().is_zero());

        let matrix_ports = [
            representation.slot(PartialIndex::open(0)),
            representation.slot(PartialIndex::open(1)),
        ];
        let zero_matrix = StructuredAtom::new(
            Atom::Zero,
            PartialStructure::from_logical_slots(matrix_ports),
        );
        let matrix = tensor("structured_zero_matrix_partner", &[representation; 2]);
        let composed = compose(
            &zero_matrix,
            &matrix,
            MatrixChannel {
                input: 0,
                output: 1,
            },
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        assert_eq!(composed.rank(), 2);
        assert!(composed.atom.as_view().is_zero());

        let traced = trace(
            &zero_matrix,
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        assert!(traced.is_scalar());
        assert!(traced.atom.as_view().is_zero());
    }

    #[test]
    fn surviving_ports_skip_internal_contraction_dummies() {
        let left_symbol = SPENSO_TAG.tensor_symbol("dummy_skip_left");
        let right_symbol = SPENSO_TAG.tensor_symbol("dummy_skip_right");
        let ports = [
            rep().slot(PartialIndex::open(0)),
            rep().slot(PartialIndex::open(1)),
        ];
        let left = partial_tensor(left_symbol, &ports, &ports);
        let right = partial_tensor(right_symbol, &ports, &ports);
        let contracted = contract(&left, &right, PortPair { left: 0, right: 0 }).unwrap();
        let indexed = materialize_interface_ports(
            &contracted,
            &HashMap::from([(0, AbstractIndex::Normal(23))]),
        )
        .unwrap();

        let left_slots = factor_slots(&indexed, left_symbol);
        let right_slots = factor_slots(&indexed, right_symbol);
        assert_eq!(left_slots.len(), 2);
        assert_eq!(right_slots.len(), 1);
        assert_eq!(left_slots[0].aind(), right_slots[0].aind());
        assert_ne!(left_slots[0].aind(), AbstractIndex::Normal(23));
        assert_eq!(left_slots[1].aind(), AbstractIndex::Normal(23));
    }

    #[test]
    fn open_contractions_allocate_fresh_dummies() {
        let ports = [
            rep().slot(PartialIndex::open(0)),
            rep().slot(PartialIndex::open(1)),
        ];
        let left_symbol = SPENSO_TAG.tensor_symbol("fresh_dummy_left");
        let right_symbol = SPENSO_TAG.tensor_symbol("fresh_dummy_right");
        let left = partial_tensor(left_symbol, &ports, &ports);
        let right = partial_tensor(right_symbol, &ports, &ports);
        let first = contract(&left, &right, PortPair { left: 0, right: 0 }).unwrap();
        let second = contract(&left, &right, PortPair { left: 0, right: 0 }).unwrap();

        let first_dummy = factor_slots(&first.atom, left_symbol)[0].aind();
        let first_partner = factor_slots(&first.atom, right_symbol)[0].aind();
        let second_dummy = factor_slots(&second.atom, left_symbol)[0].aind();
        let second_partner = factor_slots(&second.atom, right_symbol)[0].aind();
        assert_eq!(first_dummy, first_partner);
        assert_eq!(second_dummy, second_partner);
        assert_ne!(first_dummy, second_dummy);
    }

    #[test]
    fn open_contraction_dummy_does_not_alias_an_explicit_dummy_symbol() {
        let next_dummy = match AbstractIndex::new_dummy() {
            AbstractIndex::Dummy(index) => AbstractIndex::Dummy(index + 1),
            _ => unreachable!(),
        };
        let explicit = AbstractIndex::try_from(next_dummy.to_atom().as_view()).unwrap();
        let open = rep().slot(PartialIndex::open(0));
        let spectator = rep().slot(PartialIndex::Explicit(explicit));
        let left_symbol = SPENSO_TAG.tensor_symbol("dummy_collision_left");
        let right_symbol = SPENSO_TAG.tensor_symbol("dummy_collision_right");
        let left = partial_tensor(left_symbol, &[open, spectator], &[open, spectator]);
        let right = partial_tensor(right_symbol, &[open], &[open]);

        let contracted = contract(&left, &right, PortPair { left: 0, right: 0 }).unwrap();
        let left_slots = factor_slots(&contracted.atom, left_symbol);
        let right_slots = factor_slots(&contracted.atom, right_symbol);

        assert_eq!(left_slots[0].aind(), right_slots[0].aind());
        assert_ne!(left_slots[0].aind(), explicit);
        assert_eq!(left_slots[1].aind(), explicit);
    }

    #[test]
    fn explicit_indices_must_agree_before_contraction() {
        let left_port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(3)));
        let right_port = rep().slot(PartialIndex::Explicit(AbstractIndex::Normal(4)));
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("unequal_explicit_left"),
            &[left_port],
            &[left_port],
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("unequal_explicit_right"),
            &[right_port],
            &[right_port],
        );

        assert!(matches!(
            contract(&left, &right, PortPair { left: 0, right: 0 }),
            Err(TensorCompositionError::UnequalExplicitIndices { left: 0, right: 0 })
        ));
    }

    #[test]
    fn index_materialization_rejects_a_third_explicit_occurrence() {
        let representation = rep();
        let index = AbstractIndex::Normal(27);
        let open = representation.slot(PartialIndex::open(0));
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let value = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_index_occurrence"),
            &[open, explicit, explicit],
            &[open],
        );

        assert!(matches!(
            materialize_interface_ports(&value, &HashMap::from([(0, index)])),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: actual,
                occurrences: 3
            }) if actual == index
        ));
    }

    #[test]
    fn explicit_contract_rejects_a_third_explicit_occurrence() {
        let representation = rep();
        let index = AbstractIndex::Normal(31);
        let open = representation.slot(PartialIndex::open(0));
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_contract_occurrence_left"),
            &[open, explicit],
            &[open, explicit],
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_contract_occurrence_right"),
            &[explicit],
            &[explicit],
        );

        assert!(matches!(
            contract(&left, &right, PortPair { left: 0, right: 0 }),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: actual,
                occurrences: 3
            }) if actual == index
        ));
    }

    #[test]
    fn explicit_compose_rejects_a_third_explicit_occurrence() {
        let representation = rep();
        let index = AbstractIndex::Normal(33);
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let left_ports = [
            representation.slot(PartialIndex::open(0)),
            representation.slot(PartialIndex::open(1)),
            explicit,
        ];
        let right_ports = [explicit, representation.slot(PartialIndex::open(0))];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_compose_occurrence_left"),
            &left_ports,
            &left_ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("third_compose_occurrence_right"),
            &right_ports,
            &right_ports,
        );

        assert!(matches!(
            compose(
                &left,
                &right,
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
            ),
            Err(TensorCompositionError::InvalidExplicitMultiplicity {
                index: actual,
                occurrences: 3
            }) if actual == index
        ));
    }

    #[test]
    fn explicit_occurrences_are_counted_per_additive_branch() {
        let representation = rep();
        let index = AbstractIndex::Normal(35);
        let explicit = representation.slot(PartialIndex::Explicit(index));
        let terms = ["first", "second", "third"].map(|name| {
            partial_tensor(SPENSO_TAG.tensor_symbol(name), &[explicit], &[explicit]).atom
        });
        let sum = terms.into_iter().fold(Atom::Zero, |sum, term| sum + term);

        validate_explicit_index_occurrences(&sum).unwrap();
    }

    #[test]
    fn materialization_reuses_the_selected_index_in_every_summand() {
        let port = rep().slot(PartialIndex::open(0));
        let first = partial_tensor(
            SPENSO_TAG.tensor_symbol("sum_materialization_first"),
            &[port],
            &[port],
        );
        let second = partial_tensor(
            SPENSO_TAG.tensor_symbol("sum_materialization_second"),
            &[port],
            &[port],
        );
        let sum = StructuredAtom::new(
            first.atom + second.atom,
            PartialStructure::from_logical_slots([port]),
        );
        let materialized =
            materialize_interface_ports(&sum, &HashMap::from([(0, AbstractIndex::Normal(29))]))
                .unwrap();

        let AtomView::Add(sum) = materialized.as_view() else {
            panic!("expected a sum")
        };
        assert_eq!(sum.get_nargs(), 2);
        for term in sum.iter() {
            let AtomView::Fun(tensor) = term else {
                panic!("expected a tensor summand")
            };
            let slot =
                Slot::<LibraryRep, AbstractIndex>::try_from(tensor.iter().next().unwrap()).unwrap();
            assert_eq!(slot.aind(), AbstractIndex::Normal(29));
        }
    }

    #[test]
    fn matrix_roles_follow_self_dual_order_and_dual_orientation() {
        let self_dual = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let self_dual_ports = [
            spectator.slot(PartialIndex::open(0)),
            self_dual.slot(PartialIndex::open(1)),
            self_dual.slot(PartialIndex::open(2)),
        ];
        let self_dual_tensor = partial_tensor(
            SPENSO_TAG.tensor_symbol("self_dual_channel"),
            &self_dual_ports,
            &self_dual_ports,
        );
        assert_eq!(
            matrix_channel(&self_dual_tensor),
            Some(MatrixChannel {
                input: 1,
                output: 2
            })
        );

        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let dual = base.dual();
        let dual_ports = [
            dual.slot(PartialIndex::open(0)),
            spectator.slot(PartialIndex::open(1)),
            base.slot(PartialIndex::open(2)),
        ];
        let dual_tensor = partial_tensor(
            SPENSO_TAG.tensor_symbol("dual_channel"),
            &dual_ports,
            &dual_ports,
        );
        assert_eq!(
            matrix_channel(&dual_tensor),
            Some(MatrixChannel {
                input: 2,
                output: 0
            })
        );
    }

    #[test]
    fn gamma_and_color_chain_factors_keep_registered_argument_order() {
        let mink: Representation<LibraryRep> = Minkowski {}.new_rep(4).cast();
        let bis: Representation<LibraryRep> = Bispinor {}.new_rep(4).cast();
        let gamma_logical = [
            mink.slot(PartialIndex::open(0)),
            bis.slot(PartialIndex::open(1)),
            bis.slot(PartialIndex::open(2)),
        ];
        let gamma = partial_tensor(
            AGS.gamma,
            &[gamma_logical[1], gamma_logical[2], gamma_logical[0]],
            &gamma_logical,
        );
        let gamma_factor = chain_factors(&gamma, matrix_channel(&gamma).unwrap())
            .unwrap()
            .pop()
            .unwrap();
        let AtomView::Fun(gamma_factor) = gamma_factor.as_view() else {
            panic!("expected a gamma factor")
        };
        let gamma_args = gamma_factor.iter().collect::<Vec<_>>();
        assert!(matches!(gamma_args[0], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_in));
        assert!(
            matches!(gamma_args[1], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_out)
        );
        assert_eq!(
            Representation::<LibraryRep>::try_from(gamma_args[2]).unwrap(),
            mink
        );

        let adjoint: Representation<LibraryRep> = ColorAdjoint {}.new_rep(8).cast();
        let fundamental: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let antifundamental = fundamental.dual();
        let color_ports = [
            adjoint.slot(PartialIndex::open(0)),
            fundamental.slot(PartialIndex::open(1)),
            antifundamental.slot(PartialIndex::open(2)),
        ];
        let color = partial_tensor(CS.t, &color_ports, &color_ports);
        let color_factor = chain_factors(&color, matrix_channel(&color).unwrap())
            .unwrap()
            .pop()
            .unwrap();
        let AtomView::Fun(color_factor) = color_factor.as_view() else {
            panic!("expected a color-generator factor")
        };
        let color_args = color_factor.iter().collect::<Vec<_>>();
        assert_eq!(
            Representation::<LibraryRep>::try_from(color_args[0]).unwrap(),
            adjoint
        );
        assert!(matches!(color_args[1], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_in));
        assert!(
            matches!(color_args[2], AtomView::Var(v) if v.get_symbol() == SPENSO_TAG.chain_out)
        );
    }

    #[test]
    fn tracing_a_chain_spectator_keeps_its_original_endpoints() {
        let channel_rep = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let ports = [
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::open(2)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("spectator_trace_left"),
            &ports,
            &ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("spectator_trace_right"),
            &ports,
            &ports,
        );
        let chain = compose(
            &left,
            &right,
            MatrixChannel {
                input: 0,
                output: 1,
            },
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        let traced = trace(
            &chain,
            MatrixChannel {
                input: 2,
                output: 3,
            },
        )
        .unwrap();

        assert_eq!(traced.rank(), 2);
        assert!(
            matches!(traced.atom.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.chain)
        );
        traced.atom.validate_chain_like_nesting().unwrap();
        let materialized = materialize_interface_ports(
            &traced,
            &HashMap::from([
                (0, AbstractIndex::Normal(31)),
                (1, AbstractIndex::Normal(37)),
            ]),
        )
        .unwrap();
        let network = materialized
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(network.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn tracing_chain_endpoints_then_spectators_keeps_one_root_trace() {
        let channel_rep = rep();
        let spectator = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let ports = [
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
            spectator.slot(PartialIndex::open(2)),
        ];
        let left = partial_tensor(
            SPENSO_TAG.tensor_symbol("double_trace_left"),
            &ports,
            &ports,
        );
        let right = partial_tensor(
            SPENSO_TAG.tensor_symbol("double_trace_right"),
            &ports,
            &ports,
        );
        let chain = compose(
            &left,
            &right,
            MatrixChannel {
                input: 0,
                output: 1,
            },
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();

        let primary = trace(
            &chain,
            MatrixChannel {
                input: 0,
                output: 1,
            },
        )
        .unwrap();
        let traced = trace_unique(&primary).unwrap();

        assert!(traced.is_scalar());
        assert!(
            matches!(traced.atom.as_view(), AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.trace)
        );
        traced.atom.validate_chain_like_nesting().unwrap();
        let network = traced
            .atom
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.state.is_scalar());
        assert!(network.graph.dangling_indices().is_empty());
    }

    #[test]
    fn trace_interface_follows_canonical_cyclic_factor_rotation() {
        let channel_rep: Representation<LibraryRep> = Bispinor {}.new_rep(4).cast();
        let spectator: Representation<LibraryRep> = Minkowski {}.new_rep(4).cast();
        let left_index = AbstractIndex::Normal(43);
        let right_index = AbstractIndex::Normal(41);
        let left_ports = [
            spectator.slot(PartialIndex::Explicit(left_index)),
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
        ];
        let right_ports = [
            spectator.slot(PartialIndex::Explicit(right_index)),
            channel_rep.slot(PartialIndex::open(0)),
            channel_rep.slot(PartialIndex::open(1)),
        ];
        let left = partial_tensor(
            AGS.gamma,
            &[left_ports[1], left_ports[2], left_ports[0]],
            &left_ports,
        );
        let right = partial_tensor(
            AGS.gamma,
            &[right_ports[1], right_ports[2], right_ports[0]],
            &right_ports,
        );

        let traced = trace_unique(
            &compose(
                &left,
                &right,
                MatrixChannel {
                    input: 1,
                    output: 2,
                },
                MatrixChannel {
                    input: 1,
                    output: 2,
                },
            )
            .unwrap(),
        )
        .unwrap();

        assert_eq!(
            traced
                .interface
                .logical_slots()
                .into_iter()
                .map(|slot| slot.aind)
                .collect::<Vec<_>>(),
            vec![
                PartialIndex::Explicit(right_index),
                PartialIndex::Explicit(left_index),
            ]
        );
    }

    #[test]
    fn trace_closes_chain_in_canonical_form() {
        let first = tensor("first", &[rep(), rep()]);
        let second = tensor("second", &[rep(), rep()]);
        let chain = multiply(&first, &second).unwrap();
        let traced = trace_unique(&chain).unwrap();

        assert!(traced.is_scalar());
        assert!(
            matches!(traced.atom.as_view(), AtomView::Fun(fun) if fun.get_symbol() == SPENSO_TAG.trace)
        );
        assert!(
            shadowing::trace_parts(match traced.atom.as_view() {
                AtomView::Fun(fun) => fun,
                _ => unreachable!(),
            })
            .is_some()
        );
    }

    #[test]
    fn explicit_channels_require_distinct_endpoints() {
        let matrix = tensor("degenerate_channel", &[rep(), rep()]);
        let channel = MatrixChannel {
            input: 0,
            output: 0,
        };

        assert!(matches!(
            trace(&matrix, channel),
            Err(TensorCompositionError::DegenerateChannel { .. })
        ));
        assert!(matches!(
            compose(
                &matrix,
                &matrix,
                channel,
                MatrixChannel {
                    input: 0,
                    output: 1,
                },
            ),
            Err(TensorCompositionError::DegenerateChannel { .. })
        ));
    }

    #[test]
    fn explicit_channels_validate_representation_and_dual_orientation() {
        let minkowski = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let incompatible_left = tensor("incompatible_channel_left", &[rep(), minkowski]);
        let incompatible_right = tensor("incompatible_channel_right", &[minkowski, rep()]);
        let channel = MatrixChannel {
            input: 0,
            output: 1,
        };
        assert!(matches!(
            compose(&incompatible_left, &incompatible_right, channel, channel),
            Err(TensorCompositionError::IncompatiblePorts { .. })
        ));

        let base: Representation<LibraryRep> = ColorFundamental {}.new_rep(3).cast();
        let reversed = tensor("reversed_dual_channel", &[base.dual(), base]);
        assert!(matches!(
            trace(&reversed, channel),
            Err(TensorCompositionError::InvalidChannelOrientation {
                input: 0,
                output: 1
            })
        ));
    }
}
