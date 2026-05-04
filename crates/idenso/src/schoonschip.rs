use std::{cmp::Reverse, time::Instant};

use linnet::half_edge::subgraph::{SubSetLike, subset::SubSet};
use spenso::{
    algebra::ScalarMul,
    contraction::{Contract, ContractionError},
    network::{
        ContractScalars, ContractionStrategy, ExecutionResult, Sequential, SingleSmallestDegree,
        TensorNetworkError, TensorOrScalarOrKey,
        contract::SingleLargestDegree,
        graph::{NetworkGraph, NetworkLeaf, NetworkNode, NetworkOp},
        library::{DummyKey, DummyLibrary, function_lib::Wrap, symbolic::ETS},
        parsing::{Parse, ParseSettings},
        store::NetworkStore,
        tags::SPENSO_TAG as T,
    },
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        OrderedStructure, SlotIndex, StructureContract, TensorStructure,
        permuted::PermuteTensor,
        representation::{LibraryRep, LibrarySlot},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::MatchStack,
};

use crate::{
    W_,
    metric::not_slot,
    tensor::{SymbolicNetParse, SymbolicTensor},
};

const TRACE_SCHOONSCHIP: bool = false;

pub struct Schoonschipify<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool>;

fn is_sum(expr: &Atom) -> bool {
    matches!(expr.as_view(), AtomView::Add(_))
}

fn expression_size(expr: &Atom) -> (usize, usize) {
    (expr.as_view().get_byte_size(), expr.nterms())
}

fn trace_sum_contractions() -> bool {
    std::env::var_os("IDENSO_TRACE_SUM_CONTRACTIONS").is_some()
}

fn disable_direct_sum_contractions() -> bool {
    std::env::var_os("IDENSO_DISABLE_DIRECT_SUM_CONTRACTIONS").is_some()
}

fn trace_contraction_ordering() -> bool {
    std::env::var_os("IDENSO_TRACE_CONTRACTION_ORDERING").is_some()
}

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

fn distribute_expanded_left_sum(expanded_left: &Atom, right: &Atom) -> Atom {
    expanded_left.terms().fold(Atom::Zero, |sum, term| {
        let term = term.to_owned();
        sum + &term * right
    })
}

fn distribute_expanded_right_sum(left: &Atom, expanded_right: &Atom) -> Atom {
    expanded_right.terms().fold(Atom::Zero, |sum, term| {
        let term = term.to_owned();
        sum + left * &term
    })
}

fn distribute_smallest_expanded_sum_side(left: &Atom, right: &Atom) -> Atom {
    match (is_sum(left), is_sum(right)) {
        (true, true) => {
            if expression_size(left) <= expression_size(right) {
                distribute_expanded_left_sum(&left.expand(), right)
            } else {
                distribute_expanded_right_sum(left, &right.expand())
            }
        }
        (true, false) => distribute_expanded_left_sum(&left.expand(), right),
        (false, true) => distribute_expanded_right_sum(left, &right.expand()),
        (false, false) => left * right,
    }
}

fn multiplicative_factors(expr: AtomView<'_>) -> Vec<Atom> {
    match expr {
        AtomView::Mul(mul) => mul.iter().map(|factor| factor.to_owned()).collect(),
        _ => vec![expr.to_owned()],
    }
}

fn product_excluding(factors: &[Atom], excluded: &[bool]) -> Atom {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| !excluded[*index])
        .fold(Atom::num(1), |product, (_, factor)| product * factor)
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
    let terms: Vec<_> = match expanded_sum_side.as_view() {
        AtomView::Add(add) => add.iter().map(|term| term.to_owned()).collect(),
        _ => vec![expanded_sum_side.clone()],
    };

    let mut sum = Atom::Zero;
    for term in terms {
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

        sum += product_excluding(&factors, &consumed) * target.normalize_dots();
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

fn positive_even_power(exp: AtomView<'_>) -> bool {
    matches!(i64::try_from(exp), Ok(exp) if exp > 0 && exp % 2 == 0)
}

fn positive_odd_power(exp: AtomView<'_>) -> bool {
    matches!(i64::try_from(exp), Ok(exp) if exp > 0 && exp % 2 == 1)
}

fn matched_power(matches: &MatchStack<'_>, power: symbolica::atom::Symbol) -> i64 {
    i64::try_from(&matches.get(power).unwrap().to_atom()).unwrap()
}

fn matched_pattern(pattern: &Atom, matches: &MatchStack<'_>) -> Atom {
    pattern.to_pattern().replace_wildcards_with_matches(matches)
}

fn pow_if_needed(base: Atom, exponent: i64) -> Atom {
    match exponent {
        0 => Atom::num(1),
        1 => base,
        _ => base.pow(Atom::num(exponent)),
    }
}

fn even_power_replacement(base: Atom, exponent: i64) -> Atom {
    pow_if_needed(base, exponent / 2)
}

fn odd_power_replacement(base: Atom, square: Atom, exponent: i64) -> Atom {
    let half_power = exponent / 2;
    if half_power == 0 {
        base
    } else {
        pow_if_needed(square, half_power) * base
    }
}

struct SchoonschipSmallestDegree<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
>;

struct SchoonschipLargestDegree<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
>;

struct SchoonschipExpressionOrder<
    const METRIC: u8,
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
>;

const ORDER_MIN_LARGEST_OPERAND_BYTES: u8 = 0;
const ORDER_MIN_PRODUCT_TERMS: u8 = 1;
const ORDER_MIN_PRODUCT_BYTES: u8 = 2;
const ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES: u8 = 3;
const ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS: u8 = 4;
const ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES: u8 = 5;

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

pub struct SchoonschipSettings {
    depth_limit: Option<usize>,
    mode: SchoonschipMode,
    parse_inner_products: bool,
    expand_contracted_sums: bool,
    contraction_order: SchoonschipContractionOrder,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SchoonschipMode {
    SinglePass,
    Recursive(SchoonschipTraversal),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SchoonschipTraversal {
    DepthFirst,
    BreadthFirst,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SchoonschipContractionOrder {
    #[default]
    SmallestDegree,
    LargestDegree,
    MinLargestOperandBytes,
    MinProductTerms,
    MinProductBytes,
    SmallestDegreeMinLargestOperandBytes,
    SmallestDegreeMinProductTerms,
    SmallestDegreeMinProductBytes,
}

impl Default for SchoonschipSettings {
    fn default() -> Self {
        Self::partial()
    }
}

impl SchoonschipSettings {
    pub fn new(depth_limit: Option<usize>) -> Self {
        Self::breadth_first(depth_limit)
    }

    pub fn depth_first(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst),
            parse_inner_products: true,
            expand_contracted_sums: false,
            contraction_order: SchoonschipContractionOrder::default(),
        }
    }

    pub fn breadth_first(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst),
            parse_inner_products: true,
            expand_contracted_sums: false,
            contraction_order: SchoonschipContractionOrder::default(),
        }
    }

    pub fn single_pass(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::SinglePass,
            parse_inner_products: true,
            expand_contracted_sums: false,
            contraction_order: SchoonschipContractionOrder::default(),
        }
    }

    pub fn partial() -> Self {
        Self::new(Some(1))
    }

    pub fn with_depth(depth_limit: usize) -> Self {
        Self::new(Some(depth_limit))
    }

    pub fn breadth_first_with_depth(depth_limit: usize) -> Self {
        Self::breadth_first(Some(depth_limit))
    }

    pub fn full() -> Self {
        Self::single_pass(None)
    }

    pub fn with_expanded_contracted_sums(mut self) -> Self {
        self.expand_contracted_sums = true;
        self
    }

    pub fn with_contraction_order(mut self, order: SchoonschipContractionOrder) -> Self {
        self.contraction_order = order;
        self
    }

    fn without_parse_inner_products(mut self) -> Self {
        self.parse_inner_products = false;
        self
    }
}

pub trait Schoonschip {
    fn schoonschip(&self) -> Atom;

    fn normalize_dots(&self) -> Atom;
    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom;

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom;
}

impl Schoonschip for Atom {
    fn schoonschip(&self) -> Atom {
        self.as_view().schoonschip()
    }

    fn normalize_dots(&self) -> Atom {
        self.as_view().normalize_dots()
    }

    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_once_with_net::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(settings)
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_with_net::<EXPANDSUMS, DEEPEST, Aind>(settings)
    }
}

impl Schoonschip for AtomView<'_> {
    fn normalize_dots(&self) -> Atom {
        let index_cond = T.index_fiter(W_.i_);
        let other_index_cond = T.index_fiter(W_.j_);
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
        let stripped = T.rep_::<0, _>([W_.d_]);
        let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
        let self_dual_j = T.self_dual_::<0, _>([W_.d_, W_.j_]);
        let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
        let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
        let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);
        let dualizable_dual_j = T.dualizable_dual_::<0, _>([W_.d_, W_.j_]);
        fn rank_one_with_slot(slot: &Atom) -> Atom {
            T.rank1_::<0, _>([Atom::var(W_.c___), slot.clone()])
        }
        let self_dual_vector = rank_one_with_slot(&self_dual);
        let self_dual_vector_stripped = rank_one_with_slot(&self_dual_stripped);
        let dualizable_vector = rank_one_with_slot(&dualizable);
        let dualizable_dual_vector = rank_one_with_slot(&dualizable_dual);
        let dualizable_vector_stripped = rank_one_with_slot(&dualizable_stripped);

        let self_dual_square = ETS.metric(&self_dual_vector_stripped, &self_dual_vector_stripped);
        let self_dual_metric = function!(ETS.metric, &self_dual, &self_dual_j);
        let dualizable_metric = function!(ETS.metric, &dualizable, &dualizable_dual_j);

        // p(..,q(...,rep(d)))-> g(p(..,rep(d)),q(..,rep(d)))
        self.replace(T.rank1_::<0, _>([
            Atom::var(W_.c___),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
        ]))
        .with(ETS.metric(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &stripped]),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
        ))
        // g(rep(d,nu),p(..,rep(d)))->p(..,rep(d,nu))
        .replace(function!(
            ETS.metric,
            &self_dual,
            &self_dual_vector_stripped
        ))
        .when(index_cond.clone())
        .repeat()
        .with(self_dual_vector.clone())
        .replace(function!(
            ETS.metric,
            &dualizable,
            &dualizable_vector_stripped
        ))
        .when(index_cond.clone())
        .repeat()
        .with(dualizable_vector.clone())
        .replace(function!(
            ETS.metric,
            &dualizable_dual,
            &dualizable_vector_stripped
        ))
        .when(index_cond.clone())
        .repeat()
        .with(dualizable_dual_vector.clone())
        // g(rep(d,nu),p(..))->p(..,rep(d,nu))
        .replace(function!(
            ETS.metric,
            &self_dual,
            T.rank1_::<0, _>([Atom::var(W_.c___)])
        ))
        .when(index_cond.clone())
        .repeat()
        .with(self_dual_vector.clone())
        // Powers of a schoonschipped vector are normalized by parity:
        // p(mu)^(2n) -> g(p,p)^n and
        // p(mu)^(2n+1) -> g(p,p)^n p(mu).
        .replace(self_dual_vector.clone().pow(Atom::var(W_.n_)))
        .when(index_cond.clone() & W_.n_.filter_single(positive_even_power))
        .with_map({
            let self_dual_square = self_dual_square.clone();
            move |matches| {
                even_power_replacement(
                    matched_pattern(&self_dual_square, matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        .replace(self_dual_vector.clone().pow(Atom::var(W_.n_)))
        .when(index_cond.clone() & W_.n_.filter_single(positive_odd_power))
        .with_map({
            let self_dual_square = self_dual_square.clone();
            let self_dual_vector = self_dual_vector.clone();
            move |matches| {
                odd_power_replacement(
                    matched_pattern(&self_dual_vector, matches),
                    matched_pattern(&self_dual_square, matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        // Metric powers follow the same parity split:
        // g(mu,nu)^(2n) -> d^n and
        // g(mu,nu)^(2n+1) -> d^n g(mu,nu).
        .replace(self_dual_metric.clone().pow(Atom::var(W_.n_)))
        .when(
            index_cond.clone()
                & other_index_cond.clone()
                & W_.n_.filter_single(positive_even_power),
        )
        .with_map(move |matches| {
            even_power_replacement(
                matched_pattern(&Atom::var(W_.d_), matches),
                matched_power(matches, W_.n_),
            )
        })
        .replace(self_dual_metric.clone().pow(Atom::var(W_.n_)))
        .when(
            index_cond.clone() & other_index_cond.clone() & W_.n_.filter_single(positive_odd_power),
        )
        .with_map({
            let self_dual_metric = self_dual_metric.clone();
            move |matches| {
                odd_power_replacement(
                    matched_pattern(&self_dual_metric, matches),
                    matched_pattern(&Atom::var(W_.d_), matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        .replace(dualizable_metric.clone().pow(Atom::var(W_.n_)))
        .when(
            index_cond.clone()
                & other_index_cond.clone()
                & W_.n_.filter_single(positive_even_power),
        )
        .with_map(move |matches| {
            even_power_replacement(
                matched_pattern(&Atom::var(W_.d_), matches),
                matched_power(matches, W_.n_),
            )
        })
        .replace(dualizable_metric.clone().pow(Atom::var(W_.n_)))
        .when(index_cond.clone() & other_index_cond & W_.n_.filter_single(positive_odd_power))
        .with_map({
            let dualizable_metric = dualizable_metric.clone();
            move |matches| {
                odd_power_replacement(
                    matched_pattern(&dualizable_metric, matches),
                    matched_pattern(&Atom::var(W_.d_), matches),
                    matched_power(matches, W_.n_),
                )
            }
        })
        // Plain metric traces are the remaining non-power case.
        .replace(function!(ETS.metric, &self_dual, &self_dual))
        .when(&index_cond)
        .with(Atom::var(W_.d_))
        .replace(function!(ETS.metric, &dualizable, &dualizable_dual))
        .when(&index_cond)
        .with(Atom::var(W_.d_))
    }

    fn schoonschip(&self) -> Atom {
        let index_cond = T.index_fiter(W_.i_);
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
        let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
        let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
        let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
        let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

        // `schoonschip()` is the broad bare-symbolic pass. It may use product
        // patterns over plain functions; the network path deliberately calls
        // `normalize_dots()` instead so these broad patterns do not pre-empt
        // network-backed contractions.
        self.replace(function!(W_.f_, W_.a___, &self_dual) * function!(W_.g_, W_.b___, &self_dual))
            .when(index_cond.clone() & not_slot(W_.a___) & not_slot(W_.b___))
            .with(ETS.metric(
                function!(W_.f_, W_.a___, &self_dual_stripped),
                function!(W_.g_, W_.b___, &self_dual_stripped),
            ))
            .replace(function!(W_.f_, W_.a___, &self_dual).pow(Atom::num(2)))
            .when(index_cond.clone() & not_slot(W_.a___))
            .with(ETS.metric(
                function!(W_.f_, W_.a___, &self_dual_stripped),
                function!(W_.f_, W_.a___, &self_dual_stripped),
            ))
            .replace(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                    * T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual]),
            )
            .when(&index_cond)
            .with(ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual_stripped]),
            ))
            .replace(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                    .pow(Atom::num(2)),
            )
            .when(&index_cond)
            .with(ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
            ))
            .replace(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable])
                    * T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_dual]),
            )
            .when(&index_cond)
            .with(ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_stripped]),
            ))
            // Bare product contraction: vector(rep(d,i)) * T(..,rep(d,i),..)
            // becomes T(..,vector(rep(d)),..). Network contraction has a more
            // structured version of this rule and does not use this pass.
            .replace(
                function!(W_.a_, W_.a___, &self_dual, W_.b___)
                    * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual]),
            )
            .when(&index_cond)
            .repeat()
            .with(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                W_.b___
            ))
            .replace(
                function!(W_.a_, W_.a___, &dualizable, W_.b___)
                    * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]),
            )
            .when(&index_cond)
            .repeat()
            .with(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                W_.b___
            ))
            .replace(
                function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)
                    * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]),
            )
            .repeat()
            .with(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped]),
                W_.b___
            ))
            .normalize_dots()
    }

    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        let normalized = self.normalize_dots();
        let mut net = normalized
            .as_view()
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                depth_limit: settings.depth_limit,
                take_first_term_from_sum: false,
                parse_inner_products: settings.parse_inner_products,
                parse_composite_scalars_as_tensors: RECURSE,
                ..Default::default()
            })
            .unwrap();
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        match settings.contraction_order {
            SchoonschipContractionOrder::SmallestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::LargestDegree => net
                .execute::<
                    Sequential,
                    SchoonschipLargestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::MinLargestOperandBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_LARGEST_OPERAND_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::MinProductTerms => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_PRODUCT_TERMS,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::MinProductBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_MIN_PRODUCT_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_LARGEST_OPERAND_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::SmallestDegreeMinProductTerms => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_TERMS,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
            SchoonschipContractionOrder::SmallestDegreeMinProductBytes => net
                .execute::<
                    Sequential,
                    SchoonschipExpressionOrder<
                        ORDER_SMALLEST_DEGREE_MIN_PRODUCT_BYTES,
                        EXPANDSUMS,
                        RECURSE,
                        DEPTH_FIRST,
                    >,
                    _,
                    _,
                    _,
                >(&lib, &Wrap {})
                .unwrap(),
        };

        match net.result().unwrap() {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(a) => match a {
                TensorOrScalarOrKey::Key { .. } => panic!("unexpected library key result"),
                TensorOrScalarOrKey::Scalar(s) => s.clone(),
                TensorOrScalarOrKey::Tensor { tensor, .. } => tensor.expression.clone(),
            },
        }
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        if let AtomView::Add(add) = self {
            return add
                .iter()
                .map(|term| term.schoonschip_with_net::<EXPANDSUMS, DEEPEST, Aind>(settings))
                .fold(Atom::Zero, |sum, term| sum + term);
        }

        let new = if settings.expand_contracted_sums {
            match (settings.mode, DEEPEST) {
                (SchoonschipMode::SinglePass, _) | (_, false) => {
                    self.schoonschip_once_with_net::<true, false, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst), true) => {
                    self.schoonschip_once_with_net::<true, true, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst), true) => {
                    self.schoonschip_once_with_net::<true, true, false, Aind>(settings)
                }
            }
        } else {
            match (settings.mode, DEEPEST) {
                (SchoonschipMode::SinglePass, _) | (_, false) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, false, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst), true) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, true, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst), true) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, true, false, Aind>(settings)
                }
            }
        };

        if TRACE_SCHOONSCHIP {
            println!(
                "New: {}",
                new.printer(SpensoPrintSettings::compact().nice_symbolica())
            );
        }

        new.normalize_dots()
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_snapshot;
    use spenso::{
        shadowing::symbolica_utils::AtomCoreExt,
        structure::{
            abstract_index::AbstractIndex,
            representation::{Minkowski, RepName, Representation},
            slot::IsAbstractSlot,
        },
    };
    use symbolica::symbol;

    use super::*;

    #[test]
    fn simple_dot() {
        let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
        let p = T.rank_one_tensor_symbol("P");
        let q = T.rank_one_tensor_symbol("Q");

        let p1 = function!(p, 1, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p2 = function!(p, 2, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p1_2 = function!(p, 1, mink.slot::<AbstractIndex, _>(2).to_atom());
        let p2_2 = function!(p, 2, mink.slot::<AbstractIndex, _>(2).to_atom());
        let p1_stripped = function!(p, 1, mink.to_symbolic([]));
        let p2_stripped = function!(p, 2, mink.to_symbolic([]));

        let q2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(1).to_atom()
        );
        let q2_2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(2).to_atom()
        );

        let q3 = function!(q, 3, mink.slot::<AbstractIndex, _>(1).to_atom());
        let q3_2 = function!(q, 3, mink.slot::<AbstractIndex, _>(2).to_atom());

        let result = (&p1 * &q2)
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
        assert_snapshot!(result.to_bare_ordered_string(),@"g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let result = (&p1 * &p1).schoonschip();
        assert_snapshot!(result.to_bare_ordered_string(), @"g(P(1,mink(D)),P(1,mink(D)))");

        let result = function!(p, 1, &p2_stripped).normalize_dots();
        assert_snapshot!(result.to_bare_ordered_string(), @"g(P(1,mink(D)),P(2,mink(D)))");

        let result = ETS
            .metric(mink.slot::<AbstractIndex, _>(1).to_atom(), &p1_stripped)
            .normalize_dots();
        assert_snapshot!(result.to_bare_ordered_string(), @"P(1,mink(D,1))");

        let result = p1.clone().pow(Atom::num(4)).normalize_dots();
        assert_snapshot!(result.to_bare_ordered_string(), @"(g(P(1,mink(D)),P(1,mink(D))))^2");

        let result = p1.clone().pow(Atom::num(3)).normalize_dots();
        assert_snapshot!(result.to_bare_ordered_string(), @"P(1,mink(D,1))*g(P(1,mink(D)),P(1,mink(D)))");

        let metric = ETS.metric(
            mink.slot::<AbstractIndex, _>(1).to_atom(),
            mink.slot::<AbstractIndex, _>(2).to_atom(),
        );
        let result = metric.clone().pow(Atom::num(4)).normalize_dots();
        assert_snapshot!(result.to_bare_ordered_string(), @"D^2");

        let result = metric.pow(Atom::num(3)).normalize_dots();
        assert_snapshot!(result.to_bare_ordered_string(), @"D*g(mink(D,1),mink(D,2))");

        let result = (&p1 * (&q2 + &p2 * &q3_2 * &q2_2))
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::partial());
        assert_snapshot!(result.to_bare_ordered_string(),@"g(P(1,mink(D)),P(2,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let result = (&p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2)))
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
        assert_snapshot!(result.to_bare_ordered_string(),@"(g(P(2,mink(D)),Q(2,bla,mink(D)))+g(Q(2,bla,mink(D)),Q(3,mink(D))))*g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let expr = (p1 + q3 * p1_2 * q2_2) * (q2 + p2);

        let result =
            expr.schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
        assert_snapshot!(result.to_bare_ordered_string(),@"(P(1,mink(D,1))+Q(3,mink(D,1))*g(P(1,mink(D)),Q(2,bla,mink(D))))*(P(2,mink(D,1))+Q(2,bla,mink(D,1)))");

        let result = expr.schoonschip_with_net::<false, true, AbstractIndex>(
            &SchoonschipSettings::full().with_expanded_contracted_sums(),
        );
        let result = result.to_bare_ordered_string();
        assert_snapshot!(result, @"g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(P(2,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))");
        assert!(!result.contains("mink(D,1)"));

        let result = expr.schoonschip_with_net::<false, true, AbstractIndex>(
            &SchoonschipSettings::partial().with_expanded_contracted_sums(),
        );
        let result = result.to_bare_ordered_string();
        assert_snapshot!(result, @"g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(P(2,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))");
        assert!(!result.contains("mink(D,1)"));
    }

    #[test]
    fn benchmark_modes_output() {
        let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
        let p = T.rank_one_tensor_symbol("P");
        let q = T.rank_one_tensor_symbol("Q");

        let p1 = function!(p, 1, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p2 = function!(p, 2, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p2_2 = function!(p, 2, mink.slot::<AbstractIndex, _>(2).to_atom());

        let q2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(1).to_atom()
        );
        let q2_2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(2).to_atom()
        );
        let q3_2 = function!(q, 3, mink.slot::<AbstractIndex, _>(2).to_atom());

        let expr = &p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2));

        let single_pass_depth_one = expr
            .clone()
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::single_pass(
                Some(1),
            ))
            .to_bare_ordered_string();
        let depth_first_depth_one = expr
            .clone()
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::partial())
            .to_bare_ordered_string();
        let breadth_first_depth_one = expr
            .clone()
            .schoonschip_with_net::<false, true, AbstractIndex>(
                &SchoonschipSettings::breadth_first(Some(1)),
            )
            .to_bare_ordered_string();
        let full_top = expr
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full())
            .to_bare_ordered_string();

        assert_ne!(single_pass_depth_one, full_top);
        assert_eq!(depth_first_depth_one, full_top);
        assert_eq!(breadth_first_depth_one, full_top);
    }
}
