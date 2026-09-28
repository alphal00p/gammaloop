use std::{cmp::Reverse, time::Instant};

use linnet::half_edge::subgraph::{SubSetLike, subset::SubSet};

use spenso::{
    algebra::ScalarMul,
    contraction::{Contract, ContractionError},
    network::{
        ContractionStrategy, ProductContraction, TensorNetworkError,
        graph::{NetworkGraph, NetworkLeaf, NetworkOperation},
        library::{DummyKey, DummyLibrary},
        parsing::StructureFromAtom,
        store::NetworkStore,
    },
    shadowing::TensorCollectExt,
    structure::{
        ApplyPendingIndexPermutation, Canonicalized, OrderedStructure, SlotIndex,
        StructureContract, TensorStructure,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
};

use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    id::Pattern,
};

use crate::tensor::{SymbolicTensor, inference::InterfaceInference};

use super::{
    Schoonschip, SchoonschipSettings,
    api::NetworkSchoonschip,
    utils::{
        TRACE_SCHOONSCHIP, disable_direct_sum_contractions, distribute_smallest_expanded_sum_side,
        expression_size, is_sum, multiplicative_factors, product_excluding,
        trace_contraction_ordering, trace_direct_sum_term_expressions, trace_direct_sum_terms,
        trace_finish_contracts, trace_sum_contractions,
    },
};

pub struct Schoonschipify<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool>;
fn expression_order_metric_name<const METRIC: u8>() -> &'static str {
    match METRIC {
        ORDER_MIN_LARGEST_OPERAND_BYTES => "min_largest_operand_bytes",
        ORDER_MIN_PRODUCT_TERMS => "min_product_terms",
        ORDER_MIN_PRODUCT_BYTES => "min_product_bytes",
        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES => {
            "smallest_degree_min_largest_operand_bytes"
        }
        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS => "smallest_degree_min_product_terms",
        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES => "smallest_degree_min_product_bytes",
        _ => "unknown",
    }
}
fn tensor_slot_pos<Aind: AbsInd + ParseableAind>(
    tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    slot: &LibrarySlot<Aind>,
) -> Option<SlotIndex> {
    let slot_atom = slot.to_atom();
    (0..tensor.structure.order())
        .map(SlotIndex::from)
        .find(|&pos| {
            tensor
                .structure
                .get_slot(pos)
                .is_some_and(|s| s.to_atom() == slot_atom)
        })
}

fn parse_tensor_factor<Aind: AbsInd + DummyAind + ParseableAind>(
    factor: &Atom,
) -> Option<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
    SymbolicTensor::parse(factor.as_view())
        .ok()
        .map(Canonicalized::into_canonical)
}

fn direct_contract_factor_replacement<Aind: AbsInd + ParseableAind>(
    factor: &Atom,
    factor_tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    factor_slot: &LibrarySlot<Aind>,
    target_slot: &LibrarySlot<Aind>,
) -> Option<(Atom, Atom)> {
    let factor_pos = tensor_slot_pos(factor_tensor, factor_slot)?;

    if factor_tensor.is_metric && factor_tensor.structure.order() == 2 {
        let free_slot = metric_free_slot(factor_tensor, factor_pos)?;
        return Some((target_slot.to_atom(), free_slot.to_atom()));
    }

    if factor_tensor.structure.order() == 1 {
        let factor_slot = factor_tensor.structure.get_slot(factor_pos)?;
        let stripped = factor_slot.rep().base().to_symbolic([]);
        let contracted_expr = factor.replace(factor_slot.to_atom()).with(stripped);
        return Some((target_slot.to_atom(), contracted_expr));
    }

    None
}

fn apply_replacements_transitively(expr: &Atom, replacements: &[(Atom, Atom)]) -> Atom {
    let mut result = expr.clone();
    for _ in 0..=replacements.len() {
        let mut next = result.clone();
        for (from, to) in replacements {
            next = next.replace(from.clone()).with(to.clone());
        }
        if next == result {
            break;
        }
        result = next;
    }
    result
}

fn expression_contains_atom(expr: &Atom, atom: &Atom) -> bool {
    match atom.to_pattern() {
        Pattern::Literal(literal) => expr.contains(&literal),
        pattern => expr.replace(pattern).match_iter().next().is_some(),
    }
}

fn residual_contract_slots<Aind: AbsInd + ParseableAind>(
    expr: &Atom,
    slot_pairs: &[(LibrarySlot<Aind>, LibrarySlot<Aind>)],
) -> Vec<Atom> {
    let mut residual = Vec::new();
    for slot in slot_pairs
        .iter()
        .flat_map(|(sum_slot, target_slot)| [sum_slot.to_atom(), target_slot.to_atom()])
    {
        if !residual.contains(&slot) {
            residual.push(slot);
        }
    }
    residual.retain(|slot| expression_contains_atom(expr, slot));
    residual
}

fn product_from_factors<'a>(factors: impl IntoIterator<Item = &'a Atom>) -> Atom {
    factors
        .into_iter()
        .fold(Atom::num(1), |product, factor| product * factor)
}

const LOCAL_BOUNDARY_MAX_EXPANDED_TERMS: usize = 512;
const LOCAL_BOUNDARY_MAX_PRODUCT_BYTES: usize = 100_000;
const LOCAL_BOUNDARY_MAX_TARGET_BYTES: usize = 500_000;

fn locally_expand_residual_product_boundary(expr: &Atom, residual_slots: &[Atom]) -> Option<Atom> {
    let trace = trace_direct_sum_terms();
    match expr.as_view() {
        AtomView::Add(add) => {
            let mut changed = false;
            let mut changed_terms = 0usize;
            let mut sum = Atom::Zero;
            for term in add.iter() {
                let term = term.to_owned();
                let term_has_residual = residual_slots
                    .iter()
                    .any(|slot| expression_contains_atom(&term, slot));
                if let Some(cleaned) =
                    locally_expand_residual_product_boundary(&term, residual_slots)
                {
                    changed = true;
                    changed_terms += 1;
                    sum += cleaned;
                } else if term_has_residual {
                    if trace {
                        eprintln!(
                            "local_boundary add abort term_has_residual=true bytes={}",
                            term.as_view().get_byte_size()
                        );
                    }
                    return None;
                } else {
                    sum += term;
                }
            }
            if trace {
                eprintln!(
                    "local_boundary add terms={} changed_terms={} bytes={}",
                    add.iter().count(),
                    changed_terms,
                    expr.as_view().get_byte_size()
                );
            }
            changed.then_some(sum)
        }
        AtomView::Mul(_) => {
            let factors = multiplicative_factors(expr.as_view());
            for slot in residual_slots {
                let containing_indices: Vec<_> = factors
                    .iter()
                    .enumerate()
                    .filter_map(|(index, factor)| {
                        expression_contains_atom(factor, slot).then_some(index)
                    })
                    .collect();
                if containing_indices.len() < 2 {
                    continue;
                }

                let estimated_expanded_terms = containing_indices
                    .iter()
                    .try_fold(1usize, |terms, index| {
                        terms.checked_mul(factors[*index].nterms())
                    })
                    .unwrap_or(usize::MAX);
                if estimated_expanded_terms > LOCAL_BOUNDARY_MAX_EXPANDED_TERMS {
                    if trace {
                        eprintln!(
                            "local_boundary skip slot={} factors={} containing={} estimated_expanded_terms={} max_terms={}",
                            slot,
                            factors.len(),
                            containing_indices.len(),
                            estimated_expanded_terms,
                            LOCAL_BOUNDARY_MAX_EXPANDED_TERMS
                        );
                    }
                    return None;
                }

                let estimated_product_bytes: usize = containing_indices
                    .iter()
                    .map(|index| factors[*index].as_view().get_byte_size())
                    .sum();
                if estimated_product_bytes > LOCAL_BOUNDARY_MAX_PRODUCT_BYTES {
                    if trace {
                        eprintln!(
                            "local_boundary skip_product_size slot={} estimated_product_bytes={} max_bytes={}",
                            slot, estimated_product_bytes, LOCAL_BOUNDARY_MAX_PRODUCT_BYTES
                        );
                    }
                    return None;
                }

                let local_product =
                    product_from_factors(containing_indices.iter().map(|index| &factors[*index]));
                let expanded_local_product = local_product.expand();
                if trace {
                    eprintln!(
                        "local_boundary mul slot={} factors={} containing={} local_terms={} local_bytes={} expanded_terms={} expanded_bytes={}",
                        slot,
                        factors.len(),
                        containing_indices.len(),
                        local_product.nterms(),
                        local_product.as_view().get_byte_size(),
                        expanded_local_product.nterms(),
                        expanded_local_product.as_view().get_byte_size()
                    );
                }
                if expanded_local_product.nterms() > LOCAL_BOUNDARY_MAX_EXPANDED_TERMS {
                    if trace {
                        eprintln!(
                            "local_boundary skip_expanded slot={} expanded_terms={} max_terms={}",
                            slot,
                            expanded_local_product.nterms(),
                            LOCAL_BOUNDARY_MAX_EXPANDED_TERMS
                        );
                    }
                    return None;
                }
                if expanded_local_product == local_product {
                    continue;
                }

                let cleaned_local_product = expanded_local_product.schoonschip();
                if trace {
                    eprintln!(
                        "local_boundary cleaned_local terms={} bytes={} residual_slots={:?}",
                        cleaned_local_product.nterms(),
                        cleaned_local_product.as_view().get_byte_size(),
                        residual_slots
                            .iter()
                            .filter(|slot| expression_contains_atom(&cleaned_local_product, slot))
                            .map(ToString::to_string)
                            .collect::<Vec<_>>()
                    );
                }
                if expression_contains_atom(&cleaned_local_product, slot) {
                    return None;
                }
                let mut inserted_cleaned_local_product = false;
                let mut recombined = Atom::num(1);
                for (index, factor) in factors.iter().enumerate() {
                    if containing_indices.contains(&index) {
                        if !inserted_cleaned_local_product {
                            recombined *= &cleaned_local_product;
                            inserted_cleaned_local_product = true;
                        }
                    } else {
                        recombined *= factor;
                    }
                }

                return Some(recombined.schoonschip());
            }

            None
        }
        _ => None,
    }
}

fn cleanup_residual_target_boundaries<Aind: AbsInd + ParseableAind>(
    cleaned_target: &Atom,
    slot_pairs: &[(LibrarySlot<Aind>, LibrarySlot<Aind>)],
    initial_residual_slots: &[Atom],
) -> Option<Atom> {
    let target_bytes = cleaned_target.as_view().get_byte_size();
    if target_bytes > LOCAL_BOUNDARY_MAX_TARGET_BYTES {
        if trace_direct_sum_terms() {
            eprintln!(
                "local_boundary skip_target_size target_bytes={} max_bytes={}",
                target_bytes, LOCAL_BOUNDARY_MAX_TARGET_BYTES
            );
        }
        return None;
    }

    let mut cleaned = cleaned_target.clone();
    let mut residual_slots = initial_residual_slots.to_vec();
    for _ in 0..=slot_pairs.len() {
        if residual_slots.is_empty() {
            break;
        }
        let Some(next) = locally_expand_residual_product_boundary(&cleaned, &residual_slots) else {
            break;
        };
        if next == cleaned {
            break;
        }
        cleaned = next;
        residual_slots = residual_contract_slots(&cleaned, slot_pairs);
    }

    (cleaned != *cleaned_target).then_some(cleaned)
}

fn direct_contract_expanded_sum_side<Aind: AbsInd + DummyAind + ParseableAind + 'static>(
    expanded_sum_side: &Atom,
    target_expr: &Atom,
    slot_pairs: &[(LibrarySlot<Aind>, LibrarySlot<Aind>)],
) -> Option<Atom> {
    let trace_terms = trace_direct_sum_terms();
    let terms: Vec<_> = match expanded_sum_side.as_view() {
        AtomView::Add(add) => add.iter().map(|term| term.to_owned()).collect(),
        _ => vec![expanded_sum_side.clone()],
    };
    let input_term_count = terms.len();

    let mut sum = Atom::Zero;
    let mut residual_term_count = 0usize;
    let trace_term_expressions = trace_direct_sum_term_expressions();
    for (term_index, term) in terms.into_iter().enumerate() {
        let term_start = trace_terms.then(Instant::now);
        let factors = multiplicative_factors(term.as_view());
        let parsed_factors: Vec<_> = factors.iter().map(parse_tensor_factor::<Aind>).collect();
        let mut consumed = vec![false; factors.len()];
        let mut consumed_pairs = vec![false; slot_pairs.len()];
        let mut replacements = Vec::new();
        if trace_terms {
            eprintln!(
                "direct_sum_term start term_index={} factors={} term_bytes={} target_terms={} target_bytes={}",
                term_index,
                factors.len(),
                term.as_view().get_byte_size(),
                target_expr.nterms(),
                target_expr.as_view().get_byte_size()
            );
        }

        for (index, (factor, factor_tensor)) in
            factors.iter().zip(parsed_factors.iter()).enumerate()
        {
            let Some(factor_tensor) = factor_tensor else {
                continue;
            };

            let matched_pairs: Vec<_> = slot_pairs
                .iter()
                .enumerate()
                .filter(|(pair_index, (sum_slot, _))| {
                    !consumed_pairs[*pair_index]
                        && tensor_slot_pos(factor_tensor, sum_slot).is_some()
                })
                .collect();

            match matched_pairs.as_slice() {
                [] => {}
                [(_, (sum_slot, target_slot))] => {
                    replacements.push(direct_contract_factor_replacement(
                        factor,
                        factor_tensor,
                        sum_slot,
                        target_slot,
                    )?);
                    consumed[index] = true;
                    consumed_pairs[matched_pairs[0].0] = true;
                }
                [
                    (first_pair_index, (_, first_target_slot)),
                    (second_pair_index, (_, second_target_slot)),
                ] if factor_tensor.is_metric && factor_tensor.structure.order() == 2 => {
                    replacements.push((first_target_slot.to_atom(), second_target_slot.to_atom()));
                    consumed[index] = true;
                    consumed_pairs[*first_pair_index] = true;
                    consumed_pairs[*second_pair_index] = true;
                }
                _ => return None,
            }
        }

        if consumed_pairs.iter().any(|consumed| !consumed) {
            return None;
        }

        let remaining = product_excluding(&factors, &consumed);
        let target = apply_replacements_transitively(target_expr, &replacements);
        let cleanup_start = trace_terms.then(Instant::now);
        let cleaned_target = target.schoonschip();
        if let Some(start) = cleanup_start {
            eprintln!(
                "direct_sum_term target_cleaned term_index={} target_terms={} target_bytes={} elapsed={:.3?}",
                term_index,
                cleaned_target.nterms(),
                cleaned_target.as_view().get_byte_size(),
                start.elapsed()
            );
        }
        let compact_reconstructed = &remaining * &cleaned_target;
        let compact_cleanup_start = trace_terms.then(Instant::now);
        let compact_cleaned = compact_reconstructed.schoonschip();
        if let Some(start) = compact_cleanup_start {
            eprintln!(
                "direct_sum_term compact_cleaned term_index={} terms={} bytes={} elapsed={:.3?}",
                term_index,
                compact_cleaned.nterms(),
                compact_cleaned.as_view().get_byte_size(),
                start.elapsed()
            );
        }
        let residual_start = trace_terms.then(Instant::now);
        let compact_residual_slots = residual_contract_slots(&compact_cleaned, slot_pairs);
        if let Some(start) = residual_start {
            eprintln!(
                "direct_sum_term compact_residual_scan term_index={} residual_slots={:?} elapsed={:.3?}",
                term_index,
                compact_residual_slots
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>(),
                start.elapsed()
            );
        }
        let (reconstructed, cleaned) = if compact_residual_slots.is_empty() {
            (compact_reconstructed, compact_cleaned)
        } else {
            let fallback_distribute = || {
                if trace_terms {
                    eprintln!("direct_sum_term fallback_distribute term_index={term_index}");
                }
                let reconstructed = distribute_smallest_expanded_sum_side(&remaining, &target);
                let fallback_cleanup_start = trace_terms.then(Instant::now);
                let pattern_cleaned = reconstructed.schoonschip();
                let cleaned = if residual_contract_slots(&pattern_cleaned, slot_pairs).is_empty() {
                    pattern_cleaned
                } else {
                    recursive_schoonschip::<true, false, Aind>(&reconstructed).ok()?
                };
                if let Some(start) = fallback_cleanup_start {
                    eprintln!(
                        "direct_sum_term fallback_cleaned term_index={} terms={} bytes={} elapsed={:.3?}",
                        term_index,
                        cleaned.nterms(),
                        cleaned.as_view().get_byte_size(),
                        start.elapsed()
                    );
                }
                Some((reconstructed, cleaned))
            };

            if let Some(locally_cleaned_target) = cleanup_residual_target_boundaries(
                &cleaned_target,
                slot_pairs,
                &compact_residual_slots,
            ) {
                let local_cleanup_start = trace_terms.then(Instant::now);
                let local_reconstructed = &remaining * &locally_cleaned_target;
                let local_cleaned = local_reconstructed.schoonschip();
                let local_residual_slots = residual_contract_slots(&local_cleaned, slot_pairs);
                if let Some(start) = local_cleanup_start {
                    eprintln!(
                        "direct_sum_term local_boundary_cleaned term_index={} terms={} bytes={} residual_slots={:?} elapsed={:.3?}",
                        term_index,
                        local_cleaned.nterms(),
                        local_cleaned.as_view().get_byte_size(),
                        local_residual_slots
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>(),
                        start.elapsed()
                    );
                }

                if local_residual_slots.is_empty() {
                    (local_reconstructed, local_cleaned)
                } else {
                    fallback_distribute()?
                }
            } else {
                fallback_distribute()?
            }
        };
        if let Some(start) = term_start {
            eprintln!(
                "direct_sum_term done term_index={} cleaned_terms={} cleaned_bytes={} elapsed={:.3?}",
                term_index,
                cleaned.nterms(),
                cleaned.as_view().get_byte_size(),
                start.elapsed()
            );
        }
        if trace_terms {
            let residual_slots = residual_contract_slots(&cleaned, slot_pairs);
            if !residual_slots.is_empty() {
                residual_term_count += 1;
                eprintln!(
                    "direct_sum_term residual term_index={} factors={} consumed={:?} consumed_pairs={:?} term_bytes={} remaining_bytes={} target_bytes={} cleaned_bytes={} residual_slots={:?}",
                    term_index,
                    factors.len(),
                    consumed,
                    consumed_pairs,
                    term.as_view().get_byte_size(),
                    remaining.as_view().get_byte_size(),
                    target.as_view().get_byte_size(),
                    cleaned.as_view().get_byte_size(),
                    residual_slots
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                );
                if trace_term_expressions {
                    eprintln!(
                        "direct_sum_term residual_expr term_index={}\nterm={}\nremaining={}\ntarget={}\nreconstructed={}\ncleaned={}",
                        term_index,
                        term.as_view().to_plain_string(),
                        remaining.as_view().to_plain_string(),
                        target.as_view().to_plain_string(),
                        reconstructed.as_view().to_plain_string(),
                        cleaned.as_view().to_plain_string()
                    );
                }
            }
        }
        sum += cleaned;
    }

    if trace_terms {
        let residual_slots = residual_contract_slots(&sum, slot_pairs);
        eprintln!(
            "direct_sum_terms summary terms={} residual_term_count={} sum_terms={} residual_slots={:?}",
            input_term_count,
            residual_term_count,
            sum.nterms(),
            residual_slots
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
        );
    }

    Some(sum)
}

fn contracted_slot_pairs<Aind: AbsInd + ParseableAind>(
    left: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    right: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    left_positions: &SubSet<SlotIndex>,
    right_positions: &SubSet<SlotIndex>,
) -> Option<Vec<(LibrarySlot<Aind>, LibrarySlot<Aind>)>> {
    if left_positions.n_included() != right_positions.n_included() {
        return None;
    }

    let right_slots: Vec<_> = right_positions
        .included_iter()
        .map(|pos| right.structure.get_slot(pos))
        .collect::<Option<_>>()?;
    let mut used_right_slots = vec![false; right_slots.len()];
    let mut pairs = Vec::new();

    for left_pos in left_positions.included_iter() {
        let left_slot = left.structure.get_slot(left_pos)?;
        let left_slot_atom = left_slot.to_atom();
        let right_index = right_slots
            .iter()
            .enumerate()
            .find_map(|(index, right_slot)| {
                (!used_right_slots[index] && right_slot.to_atom() == left_slot_atom)
                    .then_some(index)
            })?;
        used_right_slots[right_index] = true;
        pairs.push((left_slot, right_slots[right_index]));
    }

    Some(pairs)
}

fn structure_contains_slot<Aind: AbsInd + ParseableAind>(
    structure: &OrderedStructure<LibraryRep, Aind>,
    slot_atom: &Atom,
) -> bool {
    (0..structure.order()).any(|pos| {
        structure
            .get_slot(SlotIndex::from(pos))
            .is_some_and(|slot| slot.to_atom() == *slot_atom)
    })
}

fn removed_slots_still_in_expression<Aind: AbsInd + ParseableAind>(
    left: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    right: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    left_positions: &SubSet<SlotIndex>,
    right_positions: &SubSet<SlotIndex>,
    result_structure: &OrderedStructure<LibraryRep, Aind>,
    result: &Atom,
) -> Vec<String> {
    let mut residual = Vec::new();
    for slot in left_positions
        .included_iter()
        .filter_map(|pos| left.structure.get_slot(pos))
        .chain(
            right_positions
                .included_iter()
                .filter_map(|pos| right.structure.get_slot(pos)),
        )
    {
        let slot_atom = slot.to_atom();
        if !structure_contains_slot(result_structure, &slot_atom)
            && result
                .replace(slot_atom.clone())
                .match_iter()
                .next()
                .is_some()
        {
            let slot_name = slot_atom.to_string();
            if !residual.contains(&slot_name) {
                residual.push(slot_name);
            }
        }
    }
    residual
}

fn direct_contract_smallest_expanded_sum_side<
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>(
    left: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    right: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    left_positions: &SubSet<SlotIndex>,
    right_positions: &SubSet<SlotIndex>,
    left_expr: &Atom,
    right_expr: &Atom,
) -> Option<Atom> {
    let slot_pairs = contracted_slot_pairs(left, right, left_positions, right_positions)?;

    if expression_size(left_expr) <= expression_size(right_expr) {
        direct_contract_expanded_sum_side::<Aind>(
            &left_expr.collect_tensors(),
            right_expr,
            &slot_pairs,
        )
    } else {
        let reversed_slot_pairs: Vec<_> = slot_pairs
            .into_iter()
            .map(|(left_slot, right_slot)| (right_slot, left_slot))
            .collect();
        direct_contract_expanded_sum_side::<Aind>(
            &right_expr.collect_tensors(),
            left_expr,
            &reversed_slot_pairs,
        )
    }
}

pub(super) struct SchoonschipSmallestDegree<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
>;

pub(super) struct SchoonschipLargestDegree<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
>;

pub(super) struct SchoonschipExpressionOrder<
    const METRIC: u8,
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
>;

pub(super) const ORDER_MIN_LARGEST_OPERAND_BYTES: u8 = 0;
pub(super) const ORDER_MIN_PRODUCT_TERMS: u8 = 1;
pub(super) const ORDER_MIN_PRODUCT_BYTES: u8 = 2;
pub(super) const ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES: u8 = 3;
pub(super) const ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS: u8 = 4;
pub(super) const ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES: u8 = 5;

impl<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool>
    SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>
{
    fn simplify_scalar_tensors<Aind: AbsInd + DummyAind + ParseableAind + 'static>(
        executor: &mut NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
    ) -> Result<(), ContractionError> {
        if !RECURSE {
            return Ok(());
        }

        let settings = if DEPTH_FIRST {
            SchoonschipSettings::depth_first(Some(1))
        } else {
            SchoonschipSettings::breadth_first(Some(1))
        }
        .into_single_pass();

        for tensor in &mut executor.tensors {
            if tensor.structure.is_scalar() && tensor.is_composite {
                let shortcut =
                    if NetworkSchoonschip::scalar_requires_network(tensor.expression.as_view()) {
                        false
                    } else {
                        tensor.expression.normalize_dots() == tensor.expression
                    };
                tensor.invalidate_proofs();
                tensor.expression = if shortcut {
                    if EXPANDSUMS {
                        tensor.expression.expand_num()
                    } else {
                        tensor.expression.clone()
                    }
                } else {
                    tensor
                        .expression
                        .schoonschip_with_net::<EXPANDSUMS, Aind>(&settings)
                        .map_err(|error| ContractionError::Other(eyre::Report::new(error)))?
                };
                tensor.is_composite = false;
                tensor.is_metric = false;
            }
        }
        Ok(())
    }
}

impl<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>
    ContractionStrategy<
        NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
        DummyLibrary<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
        DummyKey,
        symbolica::atom::Symbol,
        Aind,
    > for SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>
where
    SymbolicTensor<OrderedStructure<LibraryRep, Aind>>: ScalarMul<Atom, Output = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
        + ApplyPendingIndexPermutation<Output = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
{
    fn contract(
        executor: &mut NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
        graph: &NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>,
        operation: &NetworkOperation<symbolica::atom::Symbol>,
        lib: &DummyLibrary<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
    ) -> Result<NetworkLeaf<DummyKey, Aind>, TensorNetworkError<DummyKey, symbolica::atom::Symbol>>
    {
        let trace = trace_contraction_ordering();
        let start = trace.then(Instant::now);
        Self::simplify_scalar_tensors(executor)?;
        if let Some(start) = start {
            eprintln!(
                "smallest_degree phase=pre_simplify_scalars elapsed={:.3?}",
                start.elapsed()
            );
        }
        let start = trace.then(Instant::now);
        let mut product = ProductContraction::from_operation(graph, operation)?;
        let scalar_didsmth = product.contract_scalars(executor, graph, lib)?;
        if let Some(start) = start {
            eprintln!(
                "smallest_degree phase=pre_contract_scalars changed={scalar_didsmth} elapsed={:.3?}",
                start.elapsed()
            );
        }

        let mut step = 0usize;
        loop {
            let start = trace.then(Instant::now);
            let smth = product.contract_one_by_degree::<
                 false,
                 false,
                 _,
                 _,
                 _,
                 _,
                 Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                 (),
                 _,
                 _,
             >(executor, graph, lib)?;
            if let Some(start) = start {
                eprintln!(
                    "smallest_degree phase=tensor_step step={step} changed={smth} elapsed={:.3?}",
                    start.elapsed()
                );
            }
            if !smth {
                break;
            }
            product.contract_scalars(executor, graph, lib)?;
            step += 1;
        }

        let start = trace.then(Instant::now);
        Self::simplify_scalar_tensors(executor)?;
        if let Some(start) = start {
            eprintln!(
                "smallest_degree phase=post_simplify_scalars elapsed={:.3?}",
                start.elapsed()
            );
        }
        let start = trace.then(Instant::now);
        let scalar_didsmth = product.contract_scalars(executor, graph, lib)?;
        if let Some(start) = start {
            eprintln!(
                "smallest_degree phase=post_contract_scalars changed={scalar_didsmth} elapsed={:.3?}",
                start.elapsed()
            );
        }

        product.finish(executor, graph, lib)
    }
}

fn expression_order_score<const METRIC: u8>(
    degree: u32,
    left_size: (usize, usize),
    right_size: (usize, usize),
) -> (u8, u128, u128, Reverse<u32>, u128) {
    let (left_bytes, left_terms) = (left_size.0 as u128, left_size.1 as u128);
    let (right_bytes, right_terms) = (right_size.0 as u128, right_size.1 as u128);
    let non_internal_penalty = u8::from(degree == 0);
    let max_operand_bytes = left_bytes.max(right_bytes);
    let sum_operand_bytes = left_bytes + right_bytes;
    let product_terms = left_terms * right_terms;
    let product_bytes = left_bytes * right_bytes;

    match METRIC {
        ORDER_MIN_LARGEST_OPERAND_BYTES => (
            non_internal_penalty,
            max_operand_bytes,
            sum_operand_bytes,
            Reverse(degree),
            product_terms,
        ),
        ORDER_MIN_PRODUCT_TERMS => (
            non_internal_penalty,
            product_terms,
            max_operand_bytes,
            Reverse(degree),
            sum_operand_bytes,
        ),
        ORDER_MIN_PRODUCT_BYTES => (
            non_internal_penalty,
            product_bytes,
            product_terms,
            Reverse(degree),
            max_operand_bytes,
        ),
        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES => (
            non_internal_penalty,
            degree as u128,
            max_operand_bytes,
            Reverse(degree),
            sum_operand_bytes,
        ),
        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS => (
            non_internal_penalty,
            degree as u128,
            product_terms,
            Reverse(degree),
            max_operand_bytes,
        ),
        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES => (
            non_internal_penalty,
            degree as u128,
            product_bytes,
            Reverse(degree),
            product_terms,
        ),
        _ => unreachable!("unknown symbolic contraction ordering metric"),
    }
}

impl<
    const METRIC: u8,
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>
    ContractionStrategy<
        NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
        DummyLibrary<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
        DummyKey,
        symbolica::atom::Symbol,
        Aind,
    > for SchoonschipExpressionOrder<METRIC, EXPANDSUMS, RECURSE, DEPTH_FIRST>
where
    SymbolicTensor<OrderedStructure<LibraryRep, Aind>>: ScalarMul<Atom, Output = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
        + ApplyPendingIndexPermutation<Output = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
{
    fn contract(
        executor: &mut NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
        graph: &NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>,
        operation: &NetworkOperation<symbolica::atom::Symbol>,
        lib: &DummyLibrary<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
    ) -> Result<NetworkLeaf<DummyKey, Aind>, TensorNetworkError<DummyKey, symbolica::atom::Symbol>>
    {
        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        )?;
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars(executor, graph, lib)?;

        let trace_ordering = trace_contraction_ordering();
        let mut step = 0usize;

        loop {
            product.materialize_libraries(executor, graph, lib)?;
            let edge_to_contract = product.best_tensor_pair_by::<_, _, false>(
                executor,
                |_, _, degree, left, right| {
                    expression_order_score::<METRIC>(
                        degree,
                        expression_size(&left.expression),
                        expression_size(&right.expression),
                    )
                },
            );

            let Some((left, right, degree, score)) = edge_to_contract else {
                break;
            };

            if trace_ordering {
                let describe = |operand| {
                    let index = product.local_tensor_index(operand).unwrap();
                    let tensor = &executor.tensors[index];
                    let (bytes, terms) = expression_size(&tensor.expression);
                    format!(
                        "tensor#{index} degree_expr_terms={terms} bytes={bytes} structure={}",
                        tensor.structure
                    )
                };
                eprintln!(
                    "order_contract metric={} step={} degree={} score={:?} left={} right={}",
                    expression_order_metric_name::<METRIC>(),
                    step,
                    degree,
                    score,
                    describe(left),
                    describe(right),
                );
            }

            product
                 .contract_pair::<
                     _,
                     _,
                     _,
                     _,
                     Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                     (),
                     _,
                     _,
                 >(left, right, executor, graph, lib)?;
            product.contract_scalars(executor, graph, lib)?;
            step += 1;
        }

        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        )?;
        product.contract_scalars(executor, graph, lib)?;

        product.finish(executor, graph, lib)
    }
}

impl<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>
    ContractionStrategy<
        NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
        DummyLibrary<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
        DummyKey,
        symbolica::atom::Symbol,
        Aind,
    > for SchoonschipLargestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>
where
    SymbolicTensor<OrderedStructure<LibraryRep, Aind>>: ScalarMul<Atom, Output = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>
        + ApplyPendingIndexPermutation<Output = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
{
    fn contract(
        executor: &mut NetworkStore<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, Atom>,
        graph: &NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>,
        operation: &NetworkOperation<symbolica::atom::Symbol>,
        lib: &DummyLibrary<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>,
    ) -> Result<NetworkLeaf<DummyKey, Aind>, TensorNetworkError<DummyKey, symbolica::atom::Symbol>>
    {
        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        )?;
        let mut product = ProductContraction::from_operation(graph, operation)?;
        product.contract_scalars(executor, graph, lib)?;

        while product.contract_one_by_degree::<
             false,
             true,
             _,
             _,
             _,
             _,
             Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
             (),
             _,
             _,
         >(executor, graph, lib)?
         {
             product.contract_scalars(executor, graph, lib)?;
         }

        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        )?;
        product.contract_scalars(executor, graph, lib)?;

        product.finish(executor, graph, lib)
    }
}

fn recursive_schoonschip_settings<const DEPTH_FIRST: bool>() -> SchoonschipSettings {
    if DEPTH_FIRST {
        SchoonschipSettings::depth_first(Some(1))
    } else {
        SchoonschipSettings::breadth_first(Some(1))
    }
}

fn recursive_schoonschip<
    const EXPANDSUMS: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>(
    expr: &Atom,
) -> Result<Atom, ContractionError> {
    expr.schoonschip_with_net::<EXPANDSUMS, Aind>(
        &recursive_schoonschip_settings::<DEPTH_FIRST>().into_single_pass(),
    )
    .map_err(|error| ContractionError::Other(eyre::Report::new(error)))
}

fn finish_contract<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>(
    mut result: SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    recurse_result: bool,
    validate_result_interface: bool,
) -> Result<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>, ContractionError> {
    if validate_result_interface {
        // A replacement callback can erase a free port before recursion builds
        // its next graph. Never let that graph consume stale merge metadata.
        result
            .validate_rewritten_interface(&result.expression)
            .map_err(eyre::Report::new)?;
    }
    let trace = trace_finish_contracts();
    if recurse_result {
        let start = trace.then(Instant::now);
        result.invalidate_proofs();
        result.expression =
            recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&result.expression)?;
        if let Some(start) = start {
            eprintln!(
                "finish_contract recursive terms={} bytes={} elapsed={:.3?}",
                result.expression.nterms(),
                result.expression.as_view().get_byte_size(),
                start.elapsed()
            );
        }
    }
    let start = trace.then(Instant::now);
    result.invalidate_proofs();
    result.expression = result.expression.normalize_dots();
    if validate_result_interface {
        result
            .validate_rewritten_interface(&result.expression)
            .map_err(eyre::Report::new)?;
    }
    if let Some(start) = start {
        eprintln!(
            "finish_contract normalize terms={} bytes={} elapsed={:.3?}",
            result.expression.nterms(),
            result.expression.as_view().get_byte_size(),
            start.elapsed()
        );
    }
    Ok(result)
}

fn single_contracted_pos(positions: &SubSet<SlotIndex>) -> Option<SlotIndex> {
    // `merge` marks the slots from one side that participate in the current
    // contraction. The metric shortcut is only unambiguous for one contracted
    // metric leg and one contracted tensor leg.
    (positions.n_included() == 1)
        .then(|| positions.included_iter().next())
        .flatten()
}

fn metric_free_slot<Aind: AbsInd>(
    metric: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    contracted_pos: SlotIndex,
) -> Option<LibrarySlot<Aind>> {
    if metric.structure.order() != 2 {
        return None;
    }

    (0..metric.structure.order())
        .map(SlotIndex::from)
        .find(|&pos| pos != contracted_pos)
        .and_then(|pos| metric.structure.get_slot(pos))
}

fn contract_metric_into_tensor<Aind: AbsInd + ParseableAind>(
    metric: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    metric_positions: &SubSet<SlotIndex>,
    tensor_positions: &SubSet<SlotIndex>,
    tensor_expr: &Atom,
    structure: OrderedStructure<LibraryRep, Aind>,
) -> Option<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
    if metric.is_composite || !metric.is_metric {
        return None;
    }

    // A rank-two metric acts as an index relabeling operator:
    // g(i, j) * T(..., i, ...) -> T(..., j, ...).
    // The actual contracted leg comes from the network merge information; the
    // remaining metric leg is the slot that should be propagated into `tensor`.
    let metric_pos = single_contracted_pos(metric_positions)?;
    let tensor_pos = single_contracted_pos(tensor_positions)?;
    let free_metric_slot = metric_free_slot(metric, metric_pos)?;
    let contracted_tensor_slot = tensor.structure.get_slot(tensor_pos)?;

    Some(SymbolicTensor {
        proofs: Default::default(),
        structure,
        is_composite: tensor.is_composite,
        is_metric: tensor.is_metric,
        expression: tensor_expr
            .replace(contracted_tensor_slot.to_atom())
            .with(free_metric_slot.to_atom()),
    })
}

fn contract_rank_one_into_tensor<Aind: AbsInd + ParseableAind>(
    rank_one: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    rank_one_positions: &SubSet<SlotIndex>,
    tensor_positions: &SubSet<SlotIndex>,
    rank_one_expr: &Atom,
    tensor_expr: &Atom,
    structure: OrderedStructure<LibraryRep, Aind>,
) -> Option<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>> {
    if rank_one.is_composite || rank_one.structure.order() != 1 {
        return None;
    }

    // The rank-one shortcut is only a contraction when the merge information
    // says the vector leg is actually consumed. Without this guard an
    // unconnected vector would be dropped when it is merely multiplied by
    // another tensor.
    let rank_one_pos = single_contracted_pos(rank_one_positions)?;
    let tensor_pos = single_contracted_pos(tensor_positions)?;
    let rank_one_slot = rank_one.structure.get_slot(rank_one_pos)?;
    let contracted_tensor_slot = tensor.structure.get_slot(tensor_pos)?;
    let stripped = rank_one_slot.rep().base().to_symbolic([]);
    let contracted_expr = rank_one_expr
        .replace(rank_one_slot.to_atom())
        .with(stripped);

    Some(SymbolicTensor {
        proofs: Default::default(),
        structure,
        is_composite: true,
        is_metric: tensor.is_metric,
        expression: tensor_expr
            .replace(contracted_tensor_slot.to_atom())
            .with(contracted_expr)
            .normalize_dots(),
    })
}

impl<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>
    Contract<
        SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
        Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
    > for SymbolicTensor<OrderedStructure<LibraryRep, Aind>>
{
    type LCM = SymbolicTensor<OrderedStructure<LibraryRep, Aind>>;
    fn contract(
        &self,
        other: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    ) -> Result<Self::LCM, ContractionError> {
        if TRACE_SCHOONSCHIP {
            println!(
                "Contracting  {} {}rank {} with rank {} {} {}: \n{}\nwith\n{}\n gives:",
                if self.is_composite { "composite " } else { "" },
                if self.is_metric { "metric " } else { "" },
                self.structure.order(),
                if other.is_composite { "composite " } else { "" },
                if other.is_metric { "metric " } else { "" },
                other.structure.order(),
                self.expression,
                other.expression
            );
        }

        let validate_result_interface =
            !InterfaceInference::normalization_is_intrinsic(self.expression.as_view())
                || !InterfaceInference::normalization_is_intrinsic(other.expression.as_view());
        let (sexpr, oexpr) = if RECURSE && DEPTH_FIRST {
            (
                recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&self.expression)?,
                recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&other.expression)?,
            )
        } else {
            (self.expression.clone(), other.expression.clone())
        };
        if RECURSE && DEPTH_FIRST && validate_result_interface {
            self.validate_rewritten_interface(&sexpr)
                .map_err(eyre::Report::new)?;
            other
                .validate_rewritten_interface(&oexpr)
                .map_err(eyre::Report::new)?;
        }

        let (structure, pos_self, pos_other, _) = self.structure.merge(&other.structure)?;
        let trace_finish = trace_finish_contracts();
        if trace_finish {
            eprintln!(
                "contract_enter self_terms={} other_terms={} self_bytes={} other_bytes={} pos_self={} pos_other={} self_metric={} other_metric={} self_composite={} other_composite={}",
                sexpr.nterms(),
                oexpr.nterms(),
                sexpr.as_view().get_byte_size(),
                oexpr.as_view().get_byte_size(),
                pos_self.n_included(),
                pos_other.n_included(),
                self.is_metric,
                other.is_metric,
                self.is_composite,
                other.is_composite
            );
        }

        if self.structure.is_scalar() || other.structure.is_scalar() {
            let (sexpr, oexpr) = if RECURSE && !DEPTH_FIRST {
                (
                    recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&self.expression)?,
                    recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&other.expression)?,
                )
            } else {
                (sexpr, oexpr)
            };
            if RECURSE && !DEPTH_FIRST && validate_result_interface {
                self.validate_rewritten_interface(&sexpr)
                    .map_err(eyre::Report::new)?;
                other
                    .validate_rewritten_interface(&oexpr)
                    .map_err(eyre::Report::new)?;
            }

            return finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
                SymbolicTensor {
                    proofs: Default::default(),
                    structure,
                    is_composite: true,
                    is_metric: false,
                    expression: &sexpr * &oexpr,
                },
                false,
                validate_result_interface,
            );
        }

        // A typed zero keeps its merged free ports, but must not be inserted
        // as a vector argument by the shortcuts below. Scalar operands retain
        // their existing recursive callback schedule above.
        if sexpr.as_view().is_zero() || oexpr.as_view().is_zero() {
            return Ok(Self {
                proofs: Default::default(),
                structure,
                is_composite: true,
                is_metric: false,
                expression: Atom::Zero,
            });
        }

        // Metrics get first chance after scalar handling because they preserve
        // the non-metric tensor expression shape and only rewrite one slot.
        if trace_finish {
            eprintln!("contract_try metric_self_into_other");
        }
        if let Some(result) = contract_metric_into_tensor(
            self,
            other,
            &pos_self,
            &pos_other,
            &oexpr,
            structure.clone(),
        ) {
            return finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
                result,
                RECURSE && !DEPTH_FIRST,
                validate_result_interface,
            );
        }

        if trace_finish {
            eprintln!("contract_try metric_other_into_self");
        }
        if let Some(result) = contract_metric_into_tensor(
            other,
            self,
            &pos_other,
            &pos_self,
            &sexpr,
            structure.clone(),
        ) {
            return finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
                result,
                RECURSE && !DEPTH_FIRST,
                validate_result_interface,
            );
        }

        // Rank-one tensors are the original Schoonschip contraction shortcut:
        // their contracted slot is stripped from the vector and inserted into
        // the tensor slot selected by the merge.
        if trace_finish {
            eprintln!("contract_try rank_one_self_into_other");
        }
        if let Some(result) = contract_rank_one_into_tensor(
            self,
            other,
            &pos_self,
            &pos_other,
            &sexpr,
            &oexpr,
            structure.clone(),
        ) {
            return finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
                result,
                RECURSE && !DEPTH_FIRST,
                validate_result_interface,
            );
        }

        if trace_finish {
            eprintln!("contract_try rank_one_other_into_self");
        }
        if let Some(result) = contract_rank_one_into_tensor(
            other,
            self,
            &pos_other,
            &pos_self,
            &oexpr,
            &sexpr,
            structure.clone(),
        ) {
            return finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
                result,
                RECURSE && !DEPTH_FIRST,
                validate_result_interface,
            );
        }

        if trace_finish {
            eprintln!("contract_try generic_product");
        }
        let expression = if EXPANDSUMS
            && pos_self.n_included() > 0
            && is_sum(&sexpr)
            && is_sum(&oexpr)
        {
            // Only distribute genuine sum-by-sum contractions. One-sided
            // sums such as p(mu) * sum(mu) or g(mu,nu) * sum(mu) should be
            // handled by the network-informed slot replacement above, not
            // forced through expansion.
            let trace = trace_sum_contractions();
            let start = trace.then(Instant::now);
            let direct = if disable_direct_sum_contractions() {
                None
            } else {
                direct_contract_smallest_expanded_sum_side::<Aind>(
                    other, self, &pos_other, &pos_self, &oexpr, &sexpr,
                )
            };

            if let Some(result) = direct {
                if let Some(start) = start {
                    let residual_removed = removed_slots_still_in_expression(
                        self, other, &pos_self, &pos_other, &structure, &result,
                    );
                    eprintln!(
                        "sum_contract direct slots={} left_terms={} right_terms={} left_bytes={} right_bytes={} out_terms={} out_bytes={} residual_removed_slots={:?} elapsed={:.3?}",
                        pos_self.n_included(),
                        oexpr.nterms(),
                        sexpr.nterms(),
                        oexpr.as_view().get_byte_size(),
                        sexpr.as_view().get_byte_size(),
                        result.nterms(),
                        result.as_view().get_byte_size(),
                        residual_removed,
                        start.elapsed()
                    );
                }
                result
            } else {
                let fallback_start = trace.then(Instant::now);
                let settings = recursive_schoonschip_settings::<DEPTH_FIRST>().into_single_pass();
                let result = distribute_smallest_expanded_sum_side(&oexpr, &sexpr)
                    .schoonschip_with_net::<false, Aind>(&settings)
                    .map_err(|error| ContractionError::Other(eyre::Report::new(error)))?;
                if let Some(start) = fallback_start {
                    let residual_removed = removed_slots_still_in_expression(
                        self, other, &pos_self, &pos_other, &structure, &result,
                    );
                    eprintln!(
                        "sum_contract fallback slots={} left_terms={} right_terms={} left_bytes={} right_bytes={} out_terms={} out_bytes={} residual_removed_slots={:?} elapsed={:.3?}",
                        pos_self.n_included(),
                        oexpr.nterms(),
                        sexpr.nterms(),
                        oexpr.as_view().get_byte_size(),
                        sexpr.as_view().get_byte_size(),
                        result.nterms(),
                        result.as_view().get_byte_size(),
                        residual_removed,
                        start.elapsed()
                    );
                }
                result
            }
        } else {
            &oexpr * &sexpr
        };

        finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
            Self {
                proofs: Default::default(),
                structure,
                is_composite: true,
                is_metric: false,
                expression,
            },
            RECURSE && !DEPTH_FIRST,
            validate_result_interface,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::network::library::symbolic::ETS;
    use std::sync::{Arc, Mutex};
    use symbolica::atom::FunctionBuilder;

    fn callback_contraction<const RECURSE: bool, const DEPTH_FIRST: bool>(
        metric: &SymbolicTensor,
        tensor: &SymbolicTensor,
    ) -> Result<SymbolicTensor, ContractionError> {
        <SymbolicTensor as Contract<
            SymbolicTensor,
            Schoonschipify<false, RECURSE, DEPTH_FIRST>,
        >>::contract(metric, tensor)
    }

    #[test]
    fn network_rejects_callback_rank_loss_before_reusing_merged_ports() {
        crate::test_support::test_initialize();
        let parse = |source| {
            Atom::parse(
                source,
                "network_callback_rank_loss",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap()
        };
        let a = parse("spenso::mink(4,a)");
        let b = parse("spenso::mink(4,b)");
        let events = Arc::new(Mutex::new(Vec::new()));
        let seen = Arc::clone(&events);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "network_callback_rank_loss_T",
            norm = move |node, out| {
                if let AtomView::Fun(function) = node
                    && function.iter().last() == Some(target.as_view())
                {
                    seen.lock().unwrap().push(node.to_owned());
                    **out = Atom::one();
                }
            }
        );
        let metric_atom = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let tensor_atom = FunctionBuilder::new(head).add_arg(&a).finish();
        let metric = SymbolicTensor::parse(metric_atom.as_view())
            .unwrap()
            .into_canonical();
        let tensor = SymbolicTensor::parse(tensor_atom.as_view())
            .unwrap()
            .into_canonical();
        for contract in [
            callback_contraction::<false, false>,
            callback_contraction::<true, false>,
            callback_contraction::<true, true>,
        ] {
            events.lock().unwrap().clear();
            let error = contract(&metric, &tensor).unwrap_err();
            assert!(error.to_string().contains("compatible tensor interface"));
            assert_eq!(events.lock().unwrap().len(), 1);
        }
    }

    #[test]
    fn network_callback_validation_preserves_zero_ports_and_does_not_replay() {
        crate::test_support::test_initialize();
        let parse = |source| {
            Atom::parse(
                source,
                "network_callback_zero",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap()
        };
        let a = parse("spenso::mink(4,a)");
        let b = parse("spenso::mink(4,b)");
        let events = Arc::new(Mutex::new(Vec::new()));
        let seen = Arc::clone(&events);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "network_callback_zero_T",
            norm = move |node, out| {
                if let AtomView::Fun(function) = node
                    && function.iter().last() == Some(target.as_view())
                {
                    seen.lock().unwrap().push(node.to_owned());
                    **out = Atom::Zero;
                }
            }
        );
        let metric_atom = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let tensor_atom = FunctionBuilder::new(head).add_arg(&a).finish();
        let vector_atom = FunctionBuilder::new(spenso::vector_symbol!("network_callback_zero_Q"))
            .add_arg(&b)
            .finish();
        let metric = SymbolicTensor::parse(metric_atom.as_view())
            .unwrap()
            .into_canonical();
        let tensor = SymbolicTensor::parse(tensor_atom.as_view())
            .unwrap()
            .into_canonical();
        let vector = SymbolicTensor::parse(vector_atom.as_view())
            .unwrap()
            .into_canonical();
        let planned = metric.structure.merge(&tensor.structure).unwrap().0;
        for contract in [
            callback_contraction::<false, false>,
            callback_contraction::<true, false>,
            callback_contraction::<true, true>,
        ] {
            events.lock().unwrap().clear();
            let result = contract(&metric, &tensor).unwrap();
            assert_eq!(result.expression, Atom::Zero);
            assert_eq!(result.structure, planned);
            assert_eq!(events.lock().unwrap().len(), 1);
            let closed = contract(&result, &vector).unwrap();
            assert_eq!(closed.expression, Atom::Zero);
            assert!(closed.structure.is_scalar());
            assert_eq!(events.lock().unwrap().len(), 1);
        }
    }

    #[test]
    fn network_callback_validation_observes_indices_and_retains_planned_order() {
        crate::test_support::test_initialize();
        let parse = |source| {
            Atom::parse(
                source,
                "network_callback_order",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap()
        };
        let a = parse("spenso::mink(4,a)");
        let b = parse("spenso::mink(4,b)");
        let c = parse("spenso::mink(4,c)");
        let events = Arc::new(Mutex::new(Vec::new()));
        let seen = Arc::clone(&events);
        let head = spenso::tensor_symbol!(
            "network_callback_order_T",
            norm = move |node, _out| {
                seen.lock().unwrap().push(node.to_owned());
            }
        );
        let metric_atom = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let tensor_atom = FunctionBuilder::new(head).add_arg(&c).add_arg(&a).finish();
        let metric = SymbolicTensor::parse(metric_atom.as_view())
            .unwrap()
            .into_canonical();
        let tensor = SymbolicTensor::parse(tensor_atom.as_view())
            .unwrap()
            .into_canonical();
        let planned = metric.structure.merge(&tensor.structure).unwrap().0;
        events.lock().unwrap().clear();
        let expected = tensor_atom
            .replace(a.to_pattern())
            .with(b.to_pattern())
            .normalize_dots();
        let expected_events = std::mem::take(&mut *events.lock().unwrap());
        let result = callback_contraction::<false, false>(&metric, &tensor).unwrap();
        assert_eq!(result.expression, expected);
        assert_eq!(result.structure, planned);
        assert_eq!(*events.lock().unwrap(), expected_events);
        result
            .validate_rewritten_interface(&result.expression)
            .unwrap();
        assert_eq!(*events.lock().unwrap(), expected_events);
        let wrong_index = result
            .expression
            .replace(b.to_pattern())
            .with(a.to_pattern());
        assert!(result.validate_rewritten_interface(&wrong_index).is_err());
    }

    #[test]
    fn observed_network_interface_accepts_generic_index_storage() {
        use spenso::structure::abstract_index::{AbstractIndex, AbstractIndexError};

        #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        struct WrappedIndex(AbstractIndex);
        impl std::fmt::Display for WrappedIndex {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                std::fmt::Display::fmt(&self.0, f)
            }
        }
        impl AbsInd for WrappedIndex {}
        impl ParseableAind for WrappedIndex {
            type Error = AbstractIndexError;
            fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
                AbstractIndex::from_view(view).map(Self)
            }
            fn to_atom(&self) -> Atom {
                self.0.to_atom()
            }
        }

        crate::test_support::test_initialize();
        let head = spenso::tensor_symbol!("network_generic_index_T");
        let slot = Atom::parse(
            "mink(4,network_generic_index_a)",
            "spenso",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let expression = FunctionBuilder::new(head).add_arg(&slot).finish();
        let port = LibrarySlot::<WrappedIndex>::try_from(slot.as_view()).unwrap();
        let value = SymbolicTensor {
            proofs: Default::default(),
            structure: OrderedStructure::new(vec![port]).into_canonical(),
            expression,
            is_composite: false,
            is_metric: false,
        };
        value
            .validate_rewritten_interface(&value.expression)
            .unwrap();
        assert!(value.validate_rewritten_interface(&Atom::one()).is_err());
        value.validate_rewritten_interface(&Atom::Zero).unwrap();
        assert_eq!(value.structure.external_structure(), vec![port]);
    }

    #[test]
    fn observed_network_interface_keeps_distinct_open_owners() {
        use spenso::structure::abstract_index::AbstractIndex;

        crate::test_support::test_initialize();
        let slot = Atom::parse(
            "mink(4,network_open_owner_a)",
            "spenso",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let rep = LibrarySlot::<AbstractIndex>::try_from(slot.as_view())
            .unwrap()
            .rep();
        let first = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 101,
            axis: 0,
        });
        let second = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 102,
            axis: 0,
        });
        let events = Arc::new(Mutex::new(Vec::new()));
        let seen = Arc::clone(&events);
        let head = spenso::tensor_symbol!(
            "network_distinct_open_owners_T",
            norm = move |node, _out| {
                seen.lock().unwrap().push(node.to_owned());
            }
        );
        let expression = FunctionBuilder::new(head)
            .add_arg(first.to_atom())
            .add_arg(second.to_atom())
            .finish();
        let value = SymbolicTensor {
            proofs: Default::default(),
            structure: OrderedStructure::new(vec![first, second]).into_canonical(),
            expression,
            is_composite: false,
            is_metric: false,
        };
        let construction_events = events.lock().unwrap().clone();
        value
            .validate_rewritten_interface(&value.expression)
            .unwrap();
        assert_eq!(*events.lock().unwrap(), construction_events);
        assert_eq!(value.structure.order(), 2);
        assert_ne!(first.aind(), second.aind());
        let planned = value.structure.clone();
        for port in [first, second] {
            let dropped = FunctionBuilder::new(head).add_arg(port.to_atom()).finish();
            let before_validation = events.lock().unwrap().clone();
            assert!(value.validate_rewritten_interface(&dropped).is_err());
            assert_eq!(*events.lock().unwrap(), before_validation);
            assert_eq!(value.structure, planned);
        }
    }

    #[test]
    fn network_callback_rejects_changed_encoded_open_owner() {
        use spenso::structure::{
            abstract_index::AbstractIndex,
            dimension::Dimension,
            representation::{ExtendibleReps, RepName},
        };

        crate::test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let first = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 301,
            axis: 0,
        });
        let second = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 302,
            axis: 0,
        });
        let changed = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 303,
            axis: 0,
        });
        let contracted = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(41));
        let target = first.to_atom();
        let replacement = FunctionBuilder::new(spenso::tensor_symbol!("network_changed_owner_U"))
            .add_arg(second.to_atom())
            .add_arg(changed.to_atom())
            .finish();
        let events = Arc::new(Mutex::new(Vec::new()));
        let seen = Arc::clone(&events);
        let head = spenso::tensor_symbol!(
            "network_changed_owner_T",
            norm = move |node, out| {
                if let AtomView::Fun(function) = node
                    && function.iter().last() == Some(target.as_view())
                {
                    seen.lock().unwrap().push(node.to_owned());
                    **out = replacement.clone();
                }
            }
        );
        let metric_atom = FunctionBuilder::new(ETS.metric)
            .add_arg(contracted.to_atom())
            .add_arg(first.to_atom())
            .finish();
        let tensor_atom = FunctionBuilder::new(head)
            .add_arg(second.to_atom())
            .add_arg(contracted.to_atom())
            .finish();
        let metric = SymbolicTensor::parse(metric_atom.as_view())
            .unwrap()
            .into_canonical();
        let tensor = SymbolicTensor::parse(tensor_atom.as_view())
            .unwrap()
            .into_canonical();
        for contract in [
            callback_contraction::<false, false>,
            callback_contraction::<true, false>,
            callback_contraction::<true, true>,
        ] {
            events.lock().unwrap().clear();
            let error = contract(&metric, &tensor).unwrap_err();
            assert!(error.to_string().contains("compatible tensor interface"));
            assert_eq!(events.lock().unwrap().len(), 1);
        }
    }

    #[test]
    fn network_callback_preserves_both_mixed_representation_argument_orders() {
        use spenso::structure::{
            ScalarTensor,
            abstract_index::AbstractIndex,
            dimension::Dimension,
            representation::{ExtendibleReps, RepName},
        };

        crate::test_support::test_initialize();
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let first = mink.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 401,
            axis: 0,
        });
        let second = euc.slot::<AbstractIndex, _>(AbstractIndex::Open {
            owner: 402,
            axis: 0,
        });
        let events = Arc::new(Mutex::new(Vec::new()));
        let seen = Arc::clone(&events);
        let head = spenso::tensor_symbol!(
            "network_mixed_owner_order_T",
            norm = move |node, _out| {
                seen.lock().unwrap().push(node.to_owned());
            }
        );
        let scalar = SymbolicTensor::new_scalar(Atom::num(2));
        for ports in [[first, second], [second, first]] {
            let expression = FunctionBuilder::new(head)
                .add_args(ports.map(|slot| slot.to_atom()))
                .finish();
            let value = SymbolicTensor::parse(expression.as_view())
                .unwrap()
                .into_canonical();
            let expected = (&scalar.expression * &expression).normalize_dots();
            let before = events.lock().unwrap().clone();
            let result = callback_contraction::<false, false>(&scalar, &value).unwrap();
            assert_eq!(result.expression, expected);
            assert_eq!(result.structure, value.structure);
            assert_eq!(*events.lock().unwrap(), before);
        }
    }

    #[test]
    fn residual_indices_preserve_literal_and_pattern_scope() {
        crate::test_support::test_initialize();
        let parse = |source| {
            Atom::parse(
                source,
                "spenso",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap()
        };
        for (source, index, expected) in [
            ("T(mink(4,a))", "mink(4,a)", true),
            ("T(mink(4,b))", "mink(4,a)", false),
            ("T(mink(4,meta(mink(4,a))))", "mink(4,a)", true),
            ("(x+y)^mink(4,a)", "mink(4,a)", true),
            ("T(mink(4,a,extra))", "mink(4,a)", false),
            ("T(mink(4,b))", "mink(4,a_)", true),
            ("T(mink(D,b))", "mink(4,a_)", false),
            ("T(mink(4,b,c))", "mink(4,a__)", true),
        ] {
            assert_eq!(
                expression_contains_atom(&parse(source), &parse(index)),
                expected,
                "{source}: {index}"
            );
        }
    }
}
