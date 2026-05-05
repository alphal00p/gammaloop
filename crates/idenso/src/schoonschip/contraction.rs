use std::{cmp::Reverse, time::Instant};

use linnet::half_edge::subgraph::{SubSetLike, subset::SubSet};

use spenso::{
    algebra::ScalarMul,
    contraction::{Contract, ContractionError},
    network::{
        ContractScalars, ContractionStrategy, SingleSmallestDegree, TensorNetworkError,
        contract::SingleLargestDegree,
        graph::{NetworkGraph, NetworkLeaf, NetworkNode, NetworkOp},
        library::{DummyKey, DummyLibrary},
        parsing::Parse,
        store::NetworkStore,
    },
    structure::{
        OrderedStructure, SlotIndex, StructureContract, TensorStructure,
        permuted::PermuteTensor,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
};

use symbolica::atom::{Atom, AtomCore, AtomView};

use crate::tensor::SymbolicTensor;

use super::{
    Schoonschip, SchoonschipSettings,
    utils::{
        TRACE_SCHOONSCHIP, disable_direct_sum_contractions, distribute_smallest_expanded_sum_side,
        expression_size, is_sum, multiplicative_factors, product_excluding,
        trace_contraction_ordering, trace_direct_sum_term_expressions, trace_direct_sum_terms,
        trace_sum_contractions,
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
    tensor: &SymbolicTensor<Aind>,
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

fn parse_tensor_factor<Aind: AbsInd + ParseableAind>(
    factor: &Atom,
) -> Option<SymbolicTensor<Aind>> {
    SymbolicTensor::parse(factor.as_view())
        .ok()
        .map(|parsed| parsed.structure)
}

fn direct_contract_factor_into_expr<Aind: AbsInd + ParseableAind>(
    factor: &Atom,
    factor_tensor: &SymbolicTensor<Aind>,
    factor_slot: &LibrarySlot<Aind>,
    target_expr: &Atom,
    target_slot: &LibrarySlot<Aind>,
) -> Option<Atom> {
    let factor_pos = tensor_slot_pos(factor_tensor, factor_slot)?;

    if factor_tensor.is_metric && factor_tensor.structure.order() == 2 {
        let free_slot = metric_free_slot(factor_tensor, factor_pos)?;
        return Some(
            target_expr
                .replace(target_slot.to_atom())
                .with(free_slot.to_atom()),
        );
    }

    if factor_tensor.structure.order() == 1 {
        let factor_slot = factor_tensor.structure.get_slot(factor_pos)?;
        let stripped = factor_slot.rep().base().to_symbolic([]);
        let contracted_expr = factor.replace(factor_slot.to_atom()).with(stripped);
        return Some(
            target_expr
                .replace(target_slot.to_atom())
                .with(contracted_expr),
        );
    }

    None
}

fn direct_contract_expanded_sum_side<Aind: AbsInd + ParseableAind>(
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
        let factors = multiplicative_factors(term.as_view());
        let parsed_factors: Vec<_> = factors.iter().map(parse_tensor_factor::<Aind>).collect();
        let mut consumed = vec![false; factors.len()];
        let mut consumed_pairs = vec![false; slot_pairs.len()];
        let mut target = target_expr.clone();

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
                    target = direct_contract_factor_into_expr(
                        factor,
                        factor_tensor,
                        sum_slot,
                        &target,
                        target_slot,
                    )?;
                    consumed[index] = true;
                    consumed_pairs[matched_pairs[0].0] = true;
                }
                [
                    (first_pair_index, (_, first_target_slot)),
                    (second_pair_index, (_, second_target_slot)),
                ] if factor_tensor.is_metric && factor_tensor.structure.order() == 2 => {
                    target = target
                        .replace(first_target_slot.to_atom())
                        .with(second_target_slot.to_atom());
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
        let reconstructed = &remaining * &target;
        let cleaned = reconstructed.schoonschip();
        if trace_terms {
            let residual_slots: Vec<_> = slot_pairs
                .iter()
                .flat_map(|(sum_slot, target_slot)| [sum_slot.to_atom(), target_slot.to_atom()])
                .filter(|slot| cleaned.replace(slot.clone()).match_iter().next().is_some())
                .map(|slot| slot.to_string())
                .collect();
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
        let residual_slots: Vec<_> = slot_pairs
            .iter()
            .flat_map(|(sum_slot, target_slot)| [sum_slot.to_atom(), target_slot.to_atom()])
            .filter(|slot| sum.replace(slot.clone()).match_iter().next().is_some())
            .map(|slot| slot.to_string())
            .collect();
        eprintln!(
            "direct_sum_terms summary terms={} residual_term_count={} sum_terms={} residual_slots={:?}",
            input_term_count,
            residual_term_count,
            sum.nterms(),
            residual_slots
        );
    }

    Some(sum)
}

fn contracted_slot_pairs<Aind: AbsInd + ParseableAind>(
    left: &SymbolicTensor<Aind>,
    right: &SymbolicTensor<Aind>,
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
    left: &SymbolicTensor<Aind>,
    right: &SymbolicTensor<Aind>,
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

fn direct_contract_smallest_expanded_sum_side<Aind: AbsInd + ParseableAind>(
    left: &SymbolicTensor<Aind>,
    right: &SymbolicTensor<Aind>,
    left_positions: &SubSet<SlotIndex>,
    right_positions: &SubSet<SlotIndex>,
    left_expr: &Atom,
    right_expr: &Atom,
) -> Option<Atom> {
    let slot_pairs = contracted_slot_pairs(left, right, left_positions, right_positions)?;

    if expression_size(left_expr) <= expression_size(right_expr) {
        direct_contract_expanded_sum_side::<Aind>(&left_expr.expand(), right_expr, &slot_pairs)
    } else {
        let reversed_slot_pairs: Vec<_> = slot_pairs
            .into_iter()
            .map(|(left_slot, right_slot)| (right_slot, left_slot))
            .collect();
        direct_contract_expanded_sum_side::<Aind>(
            &right_expr.expand(),
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
        executor: &mut NetworkStore<SymbolicTensor<Aind>, Atom>,
    ) {
        if !RECURSE {
            return;
        }

        let settings = if DEPTH_FIRST {
            SchoonschipSettings::depth_first(Some(1)).without_parse_inner_products()
        } else {
            SchoonschipSettings::breadth_first(Some(1)).without_parse_inner_products()
        };

        for tensor in &mut executor.tensors {
            if tensor.structure.is_scalar() && tensor.is_composite {
                tensor.expression = tensor
                    .expression
                    .schoonschip_with_net::<EXPANDSUMS, true, Aind>(&settings);
                tensor.is_composite = false;
                tensor.is_metric = false;
            }
        }
    }
}

impl<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>
    ContractionStrategy<
        NetworkStore<SymbolicTensor<Aind>, Atom>,
        DummyLibrary<SymbolicTensor<Aind>>,
        DummyKey,
        symbolica::atom::Symbol,
        Aind,
    > for SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>
where
    SymbolicTensor<Aind>: ScalarMul<Atom, Output = SymbolicTensor<Aind>>
        + PermuteTensor<Permuted = SymbolicTensor<Aind>>,
{
    fn contract(
        executor: &mut NetworkStore<SymbolicTensor<Aind>, Atom>,
        graph: NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>,
        lib: &DummyLibrary<SymbolicTensor<Aind>>,
    ) -> Result<
        (NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>, bool),
        TensorNetworkError<DummyKey, symbolica::atom::Symbol>,
    > {
        Self::simplify_scalar_tensors(executor);
        let (mut graph, mut didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        while {
            let (newgraph, smth) = SingleSmallestDegree::<
                false,
                Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
            >::contract(executor, graph, lib)?;
            graph = newgraph;
            smth
        } {
            didsmth = true;
        }

        Self::simplify_scalar_tensors(executor);
        let (graph, scalar_didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        Ok((graph, didsmth || scalar_didsmth))
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
        NetworkStore<SymbolicTensor<Aind>, Atom>,
        DummyLibrary<SymbolicTensor<Aind>>,
        DummyKey,
        symbolica::atom::Symbol,
        Aind,
    > for SchoonschipExpressionOrder<METRIC, EXPANDSUMS, RECURSE, DEPTH_FIRST>
where
    SymbolicTensor<Aind>: ScalarMul<Atom, Output = SymbolicTensor<Aind>>
        + PermuteTensor<Permuted = SymbolicTensor<Aind>>,
{
    fn contract(
        executor: &mut NetworkStore<SymbolicTensor<Aind>, Atom>,
        graph: NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>,
        lib: &DummyLibrary<SymbolicTensor<Aind>>,
    ) -> Result<
        (NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>, bool),
        TensorNetworkError<DummyKey, symbolica::atom::Symbol>,
    > {
        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        );
        let (mut graph, mut didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        let trace_ordering = trace_contraction_ordering();
        let mut step = 0usize;

        loop {
            graph.sync_order();

            let node_size = |_nid, node: &NetworkNode<DummyKey, symbolica::atom::Symbol>| match node
            {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    Some(expression_size(&executor.tensors[*index].expression))
                }
                NetworkNode::Leaf(NetworkLeaf::LibraryKey(_)) => None,
                _ => None,
            };

            let mut last_tensor = None;
            let edge_to_contract = graph
                .graph
                .iter_nodes()
                .filter(|(_, _, node)| node.is_tensor())
                .filter_map(|(nid1, adjacencies, n1)| {
                    let mut degree = 0;
                    let mut first = None;
                    for hedge in adjacencies {
                        if graph.graph[[&hedge]].is_slot() && graph.graph.inv(hedge) != hedge {
                            first = Some(hedge);
                            degree += 1;
                        }
                    }

                    let nid2 = if degree == 0 {
                        if let Some(last_tensor) = last_tensor {
                            last_tensor
                        } else {
                            last_tensor = Some(nid1);
                            return None;
                        }
                    } else {
                        graph.graph.involved_node_id(first?)?
                    };

                    let n2 = &graph.graph[nid2];
                    last_tensor = Some(nid1);

                    Some((
                        expression_order_score::<METRIC>(
                            degree,
                            node_size(nid1, n1)?,
                            node_size(nid2, n2)?,
                        ),
                        degree,
                        nid1,
                        n1,
                        nid2,
                        n2,
                    ))
                })
                .min_by_key(|(score, _, _, _, _, _)| *score);

            let Some((score, degree, nid1, n1, nid2, n2)) = edge_to_contract else {
                break;
            };

            if trace_ordering {
                let describe = |node: &NetworkNode<DummyKey, symbolica::atom::Symbol>| match node {
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                        let tensor = &executor.tensors[*index];
                        let (bytes, terms) = expression_size(&tensor.expression);
                        format!(
                            "tensor#{index} degree_expr_terms={terms} bytes={bytes} structure={}",
                            tensor.structure
                        )
                    }
                    other => format!("{other:?}"),
                };
                eprintln!(
                    "order_contract metric={} step={} degree={} score={:?} left={} right={}",
                    expression_order_metric_name::<METRIC>(),
                    step,
                    degree,
                    score,
                    describe(n1),
                    describe(n2),
                );
            }

            let new_node = match (n1, n2) {
                (NetworkNode::Leaf(_), NetworkNode::Op(NetworkOp::Product))
                | (NetworkNode::Op(NetworkOp::Product), NetworkNode::Leaf(_)) => {
                    return Err(TensorNetworkError::SlotEdgeToProdNode);
                }
                (NetworkNode::Leaf(l1), NetworkNode::Leaf(l2)) => match (l1, l2) {
                    (NetworkLeaf::Scalar(_), _) | (_, NetworkLeaf::Scalar(_)) => {
                        return Err(TensorNetworkError::SlotEdgeToScalarNode);
                    }
                    (NetworkLeaf::LocalTensor(l1), NetworkLeaf::LocalTensor(l2)) => {
                        let contracted = <SymbolicTensor<Aind> as Contract<
                            SymbolicTensor<Aind>,
                            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                        >>::contract(
                            &executor.tensors[*l1], &executor.tensors[*l2]
                        )?;
                        let pos = executor.tensors.len();
                        executor.tensors.push(contracted);
                        NetworkLeaf::LocalTensor(pos)
                    }
                    (NetworkLeaf::LibraryKey(_), _) | (_, NetworkLeaf::LibraryKey(_)) => {
                        return Err(TensorNetworkError::FailedContractMsg(
                            "expression-order contraction does not support library tensors"
                                .to_owned(),
                        ));
                    }
                },
                (a, b) => {
                    return Err(TensorNetworkError::CannotContractEdgeBetween(
                        a.clone(),
                        b.clone(),
                    ));
                }
            };

            graph.identify_nodes_without_self_edges_merge_heads(
                &[nid1, nid2],
                NetworkNode::Leaf(new_node),
            );
            didsmth = true;
            step += 1;
        }

        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        );
        let (graph, scalar_didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        Ok((graph, didsmth || scalar_didsmth))
    }
}

impl<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>
    ContractionStrategy<
        NetworkStore<SymbolicTensor<Aind>, Atom>,
        DummyLibrary<SymbolicTensor<Aind>>,
        DummyKey,
        symbolica::atom::Symbol,
        Aind,
    > for SchoonschipLargestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>
where
    SymbolicTensor<Aind>: ScalarMul<Atom, Output = SymbolicTensor<Aind>>
        + PermuteTensor<Permuted = SymbolicTensor<Aind>>,
{
    fn contract(
        executor: &mut NetworkStore<SymbolicTensor<Aind>, Atom>,
        graph: NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>,
        lib: &DummyLibrary<SymbolicTensor<Aind>>,
    ) -> Result<
        (NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>, bool),
        TensorNetworkError<DummyKey, symbolica::atom::Symbol>,
    > {
        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        );
        let (mut graph, mut didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        while {
            let (newgraph, smth) = SingleLargestDegree::<
                false,
                Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
            >::contract(executor, graph, lib)?;
            graph = newgraph;
            smth
        } {
            didsmth = true;
        }

        SchoonschipSmallestDegree::<EXPANDSUMS, RECURSE, DEPTH_FIRST>::simplify_scalar_tensors(
            executor,
        );
        let (graph, scalar_didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        Ok((graph, didsmth || scalar_didsmth))
    }
}

fn recursive_schoonschip_settings<const DEPTH_FIRST: bool>() -> SchoonschipSettings {
    if DEPTH_FIRST {
        SchoonschipSettings::depth_first(Some(1)).without_parse_inner_products()
    } else {
        SchoonschipSettings::breadth_first(Some(1)).without_parse_inner_products()
    }
}

fn recursive_schoonschip<
    const EXPANDSUMS: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>(
    expr: &Atom,
) -> Atom {
    expr.schoonschip_with_net::<EXPANDSUMS, true, Aind>(&recursive_schoonschip_settings::<
        DEPTH_FIRST,
    >())
}

fn finish_contract<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>(
    mut result: SymbolicTensor<Aind>,
    recurse_result: bool,
) -> Result<SymbolicTensor<Aind>, ContractionError> {
    if recurse_result {
        result.expression =
            recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&result.expression);
    }
    result.expression = result.expression.normalize_dots();
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
    metric: &SymbolicTensor<Aind>,
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
    metric: &SymbolicTensor<Aind>,
    tensor: &SymbolicTensor<Aind>,
    metric_positions: &SubSet<SlotIndex>,
    tensor_positions: &SubSet<SlotIndex>,
    tensor_expr: &Atom,
    structure: OrderedStructure<LibraryRep, Aind>,
) -> Option<SymbolicTensor<Aind>> {
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
        structure,
        is_composite: tensor.is_composite,
        is_metric: tensor.is_metric,
        expression: tensor_expr
            .replace(contracted_tensor_slot.to_atom())
            .with(free_metric_slot.to_atom()),
    })
}

fn contract_rank_one_into_tensor<Aind: AbsInd + ParseableAind>(
    rank_one: &SymbolicTensor<Aind>,
    tensor: &SymbolicTensor<Aind>,
    rank_one_positions: &SubSet<SlotIndex>,
    tensor_positions: &SubSet<SlotIndex>,
    rank_one_expr: &Atom,
    tensor_expr: &Atom,
    structure: OrderedStructure<LibraryRep, Aind>,
) -> Option<SymbolicTensor<Aind>> {
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
> Contract<SymbolicTensor<Aind>, Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>>
    for SymbolicTensor<Aind>
{
    type LCM = SymbolicTensor<Aind>;
    fn contract(&self, other: &SymbolicTensor<Aind>) -> Result<Self::LCM, ContractionError> {
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

        let (sexpr, oexpr) = if RECURSE && DEPTH_FIRST {
            (
                recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&self.expression),
                recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&other.expression),
            )
        } else {
            (self.expression.clone(), other.expression.clone())
        };

        let (structure, pos_self, pos_other, _) = self.structure.merge(&other.structure)?;

        if self.structure.is_scalar() || other.structure.is_scalar() {
            let (sexpr, oexpr) = if RECURSE && !DEPTH_FIRST {
                (
                    recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&self.expression),
                    recursive_schoonschip::<EXPANDSUMS, DEPTH_FIRST, Aind>(&other.expression),
                )
            } else {
                (sexpr, oexpr)
            };

            return finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
                SymbolicTensor {
                    structure,
                    is_composite: true,
                    is_metric: false,
                    expression: &sexpr * &oexpr,
                },
                false,
            );
        }

        // Metrics get first chance after scalar handling because they preserve
        // the non-metric tensor expression shape and only rewrite one slot.
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
            );
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
            );
        }

        // Rank-one tensors are the original Schoonschip contraction shortcut:
        // their contracted slot is stripped from the vector and inserted into
        // the tensor slot selected by the merge.
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
            );
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
            );
        }

        let expression = &oexpr * &sexpr;
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
                direct_contract_smallest_expanded_sum_side(
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
                let result = distribute_smallest_expanded_sum_side(&oexpr, &sexpr)
                    .schoonschip_with_net::<false, false, Aind>(&recursive_schoonschip_settings::<
                        DEPTH_FIRST,
                    >());
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
            expression
        };

        finish_contract::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(
            Self {
                structure,
                is_composite: true,
                is_metric: false,
                expression,
            },
            RECURSE && !DEPTH_FIRST,
        )
    }
}
