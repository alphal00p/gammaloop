//! Checked term iteration shared by collection and explicit materialization.
//!
//! This is the former ComponentSum input tape. Leaves are borrowed payloads or
//! existing graph occurrences; it owns no tensor topology and builds no Atoms.

use crate::{
    network::graph::{NetworkGraph, NetworkNode, NetworkOp},
    structure::slot::AbsInd,
};
use ahash::AHashMap;
use linnet::{
    half_edge::{NodeIndex, tree::SimpleTraversalTree},
    tree::child_vec::ChildVecStore,
};
use std::{fmt::Debug, hash::Hash};
use symbolica::domains::{
    Ring,
    rational::{Q, Rational},
};

#[derive(Clone, Copy)]
enum InputFactor {
    Leaf(usize),
    Number(usize),
    Add(usize, usize),
    Mul(usize, usize),
    Pow(usize, u16),
}

/// An admitted leaf of the existing graph, or its exact numeric coefficient.
/// Equal keys must denote the same immutable payload, admission class, and
/// numeric coefficient throughout the operation.
pub enum TermLeaf<L> {
    Value(L),
    Number(L, Rational),
    /// A previously compiled subtree of this same tape. This permits literal
    /// alias bodies to reuse the existing arithmetic owner without Atom expansion.
    Reference(usize, (usize, usize)),
}

/// One operation's bounded arithmetic tape, with caller-owned leaf semantics.
/// Discard the tape after any refused operation: partial scratch state may remain.
pub struct TermTape<L> {
    nodes: Vec<InputFactor>,
    children: Vec<usize>,
    positions: AHashMap<L, usize>,
    leaves: Vec<L>,
    numbers: Vec<Rational>,
    powers: Vec<u16>,
    active: Vec<usize>,
    numeric_depth: usize,
    terms: AHashMap<Vec<(usize, u16)>, (usize, Rational)>,
    bytes: usize,
    distributed: bool,
}

impl<L> Default for TermTape<L> {
    fn default() -> Self {
        Self {
            nodes: Vec::new(),
            children: Vec::new(),
            positions: AHashMap::new(),
            leaves: Vec::new(),
            numbers: Vec::new(),
            powers: Vec::new(),
            active: Vec::new(),
            numeric_depth: 0,
            terms: AHashMap::new(),
            bytes: 0,
            distributed: false,
        }
    }
}

impl<L: Copy + Eq + Hash> TermTape<L> {
    pub const MAX_GENERATED_TERMS: usize = 1_000_000;
    pub const MAX_GENERATED_FACTOR_BYTES: usize = 64 * 1024 * 1024;
    pub const MAX_EXPANSION_DEPTH: usize = 256;

    /// Compile an existing parsed subtree directly into this arithmetic tape.
    ///
    /// The caller owns leaf admission and can stop at an entire scalar/foreign
    /// subtree. Unselected Product/Sum/positive Power operations use the graph's
    /// already established child order; function payloads are never opened by
    /// this iterator. The traversal is supplied by the same graph owner, so no
    /// Atom is rebuilt or reparsed to discover arithmetic structure.
    pub fn compile_graph<K: Debug, F: Debug, Aind: AbsInd>(
        &mut self,
        graph: &NetworkGraph<K, F, Aind>,
        tree: &SimpleTraversalTree<ChildVecStore<()>>,
        root: NodeIndex,
        leaf: &mut impl FnMut(&mut Self, NodeIndex, &NetworkNode<K, F, Aind>) -> Option<TermLeaf<L>>,
    ) -> Option<(usize, (usize, usize))> {
        self.compile_graph_node(graph, tree, root, leaf, 0)
    }

    fn compile_graph_node<K: Debug, F: Debug, Aind: AbsInd>(
        &mut self,
        graph: &NetworkGraph<K, F, Aind>,
        tree: &SimpleTraversalTree<ChildVecStore<()>>,
        node: NodeIndex,
        leaf: &mut impl FnMut(&mut Self, NodeIndex, &NetworkNode<K, F, Aind>) -> Option<TermLeaf<L>>,
        depth: usize,
    ) -> Option<(usize, (usize, usize))> {
        if depth > Self::MAX_EXPANSION_DEPTH {
            return None;
        }
        let data = &graph.graph[node];
        if let Some(value) = leaf(self, node, data) {
            return match value {
                TermLeaf::Value(value) => self.leaf(value),
                TermLeaf::Number(value, coefficient) => self.number(value, coefficient),
                TermLeaf::Reference(node, size) => {
                    (node < self.nodes.len()).then_some((node, Self::bounded_size(size.0, size.1)?))
                }
            };
        }
        match data {
            NetworkNode::Op(NetworkOp::Sum | NetworkOp::Product) => {
                let children = tree
                    .iter_children(node, &graph.graph)
                    .map(|child| self.compile_graph_node(graph, tree, child, leaf, depth + 1))
                    .collect::<Option<Vec<_>>>()?;
                self.group(
                    children.into_iter(),
                    matches!(data, NetworkNode::Op(NetworkOp::Sum)),
                )
            }
            NetworkNode::Op(NetworkOp::Power(power)) if *power > 0 => {
                let mut children = tree.iter_children(node, &graph.graph);
                let child = children.next()?;
                if children.next().is_some() {
                    return None;
                }
                let base = self.compile_graph_node(graph, tree, child, leaf, depth + 1)?;
                self.power(base, u16::try_from(*power).ok()?)
            }
            // Nonlinear/broadcast scopes are admitted as an opaque leaf by the
            // domain owner, after any requested inner operation has completed.
            // Numeric negation and nonpositive powers likewise need an explicit
            // leaf decision rather than silently changing their grammar.
            _ => None,
        }
    }

    pub fn leaves(&self) -> &[L] {
        &self.leaves
    }
    pub fn distributed(&self) -> bool {
        self.distributed
    }
    pub fn clear_terms(&mut self) {
        self.terms.clear();
    }

    /// Preserve the first generated occurrence order after exact coalescence.
    pub fn take_terms(&mut self) -> Vec<(Vec<(usize, u16)>, Rational)> {
        let mut terms = std::mem::take(&mut self.terms)
            .into_iter()
            .filter(|(_, (_, coefficient))| !coefficient.is_zero())
            .collect::<Vec<_>>();
        terms.sort_unstable_by_key(|(_, (position, _))| *position);
        terms
            .into_iter()
            .map(|(factors, (_, coefficient))| (factors, coefficient))
            .collect()
    }

    pub fn bounded_size(terms: usize, factors: usize) -> Option<(usize, usize)> {
        (terms <= Self::MAX_GENERATED_TERMS
            && factors.checked_mul(std::mem::size_of::<(L, u16)>())?
                <= Self::MAX_GENERATED_FACTOR_BYTES)
            .then_some((terms, factors))
    }
    pub fn product_size(left: (usize, usize), right: (usize, usize)) -> Option<(usize, usize)> {
        Self::bounded_size(
            left.0.checked_mul(right.0)?,
            left.1
                .checked_mul(right.0)?
                .checked_add(right.1.checked_mul(left.0)?)?,
        )
    }
    /// Account for leaf-recognition scratch under the same operation budget.
    pub fn reserve_bytes(&mut self, bytes: usize) -> Option<()> {
        let total = self.bytes.checked_add(bytes)?;
        if total > Self::MAX_GENERATED_FACTOR_BYTES {
            return None;
        }
        self.bytes = total;
        Some(())
    }
    fn node(&mut self, factor: InputFactor) -> Option<usize> {
        self.reserve_bytes(4 * std::mem::size_of::<InputFactor>())?;
        let position = self.nodes.len();
        self.nodes.push(factor);
        Some(position)
    }
    pub fn known_leaf(&self, value: &L) -> Option<usize> {
        self.positions.get(value).copied()
    }
    pub fn leaf(&mut self, value: L) -> Option<(usize, (usize, usize))> {
        if let Some(node) = self.known_leaf(&value) {
            return Some((node, (1, 1)));
        }
        self.reserve_bytes(
            4 * std::mem::size_of::<(L, usize)>()
                + std::mem::size_of::<(L, Rational, u16, usize)>(),
        )?;
        let position = self.leaves.len();
        self.leaves.push(value);
        self.powers.push(0);
        let node = self.node(InputFactor::Leaf(position))?;
        self.positions.insert(value, node);
        Some((node, (1, 1)))
    }
    pub fn number(&mut self, value: L, coefficient: Rational) -> Option<(usize, (usize, usize))> {
        if let Some(node) = self.known_leaf(&value) {
            return Some((node, (1, 1)));
        }
        self.reserve_bytes(
            4 * std::mem::size_of::<(L, usize)>()
                + std::mem::size_of::<(L, Rational, u16, usize)>(),
        )?;
        let position = self.numbers.len();
        self.numbers.push(coefficient);
        let node = self.node(InputFactor::Number(position))?;
        self.positions.insert(value, node);
        Some((node, (1, 1)))
    }
    pub fn group(
        &mut self,
        children: impl ExactSizeIterator<Item = (usize, (usize, usize))>,
        sum: bool,
    ) -> Option<(usize, (usize, usize))> {
        self.reserve_bytes(
            children
                .len()
                .checked_mul(2 * std::mem::size_of::<usize>())?,
        )?;
        let mut nodes = Vec::with_capacity(children.len());
        let mut size: (usize, usize) = if sum { (0, 0) } else { (1, 0) };
        for (node, next) in children {
            size = if sum {
                Self::bounded_size(size.0.checked_add(next.0)?, size.1.checked_add(next.1)?)?
            } else {
                Self::product_size(size, next)?
            };
            nodes.push(node);
        }
        let start = self.children.len();
        self.children.extend(nodes);
        let end = self.children.len();
        let node = self.node(if sum {
            InputFactor::Add(start, end)
        } else {
            InputFactor::Mul(start, end)
        })?;
        Some((node, size))
    }
    pub fn power(
        &mut self,
        base: (usize, (usize, usize)),
        exponent: u16,
    ) -> Option<(usize, (usize, usize))> {
        if exponent == 0 {
            return None;
        }
        let (base, mut factor) = base;
        let mut remaining = exponent;
        let mut size = (1, 0);
        while remaining != 0 {
            if remaining % 2 == 1 {
                size = Self::product_size(size, factor)?;
            }
            remaining /= 2;
            if remaining != 0 {
                factor = Self::product_size(factor, factor)?;
            }
        }
        Some((self.node(InputFactor::Pow(base, exponent))?, size))
    }

    pub fn distribute(
        &mut self,
        pending: &mut Vec<(usize, u16)>,
        coefficient: &Rational,
        depth: usize,
    ) -> Option<()> {
        if depth > Self::MAX_EXPANSION_DEPTH {
            return None;
        }
        let Some((node, exponent)) = pending.pop() else {
            let mut key: Vec<_> = self
                .active
                .iter()
                .map(|&factor| (factor, self.powers[factor]))
                .collect();
            key.sort_unstable_by_key(|&(factor, _)| factor);
            if let Some((_, existing)) = self.terms.get_mut(key.as_slice()) {
                *existing += coefficient;
            } else {
                self.reserve_bytes(
                    key.capacity()
                        .checked_mul(std::mem::size_of::<(usize, u16)>())?
                        .checked_add(
                            4 * std::mem::size_of::<(Vec<(usize, u16)>, (usize, Rational))>(),
                        )?,
                )?;
                let order = self.terms.len();
                self.terms.insert(key, (order, coefficient.clone()));
            }
            return Some(());
        };
        let length = pending.len();
        match self.nodes[node] {
            InputFactor::Add(start, end) => {
                // Common root spectators stay factored. Branch-local scalar
                // coefficients are distributed with this indexed core.
                self.distributed |= !pending.is_empty()
                    || !self.active.is_empty()
                    || self.numeric_depth != 0
                    || exponent > 1;
                if exponent > 1 {
                    pending.push((node, exponent - 1));
                }
                for position in start..end {
                    pending.push((self.children[position], 1));
                    self.distribute(pending, coefficient, depth + 1)?;
                    pending.pop();
                }
                pending.truncate(length);
            }
            InputFactor::Mul(start, end) => {
                for position in start..end {
                    pending.push((self.children[position], exponent));
                }
                self.distribute(pending, coefficient, depth + 1)?;
                pending.truncate(length);
            }
            InputFactor::Pow(base, power) => {
                pending.push((base, power.checked_mul(exponent)?));
                self.distribute(pending, coefficient, depth + 1)?;
                pending.pop();
            }
            InputFactor::Number(position) => {
                let mut next = coefficient.clone();
                next *= Q.pow(&self.numbers[position], u64::from(exponent));
                self.numeric_depth += 1;
                let result = self.distribute(pending, &next, depth + 1);
                self.numeric_depth -= 1;
                result?;
            }
            InputFactor::Leaf(factor) => {
                let previous = self.powers[factor];
                self.powers[factor] = previous.checked_add(exponent)?;
                if previous == 0 {
                    self.active.push(factor);
                }
                self.distribute(pending, coefficient, depth + 1)?;
                self.powers[factor] = previous;
                if previous == 0 {
                    self.active.pop();
                }
            }
        }
        pending.push((node, exponent));
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        network::graph::{NAdd, NMul, NetworkLeaf, ScalarRef},
        structure::abstract_index::AbstractIndex,
    };
    use symbolica::{
        atom::{Atom, AtomCore, Symbol},
        symbol,
    };

    type Graph = NetworkGraph<i8, Symbol, AbstractIndex>;

    #[test]
    fn parsed_arithmetic_coalesces_equal_leaves_before_emission() {
        // (x + 2) (x - y); source graph occurrences of x share one tape leaf.
        let left = Graph::scalar(0).n_add([Graph::scalar(2)]);
        let minus_y = Graph::scalar(3).n_mul([Graph::scalar(1)]);
        let right = Graph::scalar(0).n_add([minus_y]);
        let graph = left.n_mul([right]);
        let tree = graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = graph.graph.node_id(graph.head());
        let mut tape = TermTape::default();
        let (node, size) = tape
            .compile_graph(&graph, &tree, root, &mut |_, _, node| match node {
                NetworkNode::Leaf(NetworkLeaf::Scalar(ScalarRef::Store(index))) => {
                    Some(match *index {
                        2 => TermLeaf::Number(*index, Rational::from(2)),
                        3 => TermLeaf::Number(*index, Rational::from(-1)),
                        _ => TermLeaf::Value(*index),
                    })
                }
                _ => None,
            })
            .unwrap();
        assert_eq!(size.0, 4);
        tape.distribute(&mut vec![(node, 1)], &Rational::one(), 0)
            .unwrap();
        let x = Atom::var(symbol!("term_tape_test::x"));
        let y = Atom::var(symbol!("term_tape_test::y"));
        let factored = graph
            .to_expression_at(&tree, root, &mut |_, value| match value {
                NetworkNode::Leaf(NetworkLeaf::Scalar(ScalarRef::Store(index))) => {
                    Ok(Some(match index {
                        0 => x.clone(),
                        1 => y.clone(),
                        2 => Atom::num(2),
                        3 => Atom::num(-1),
                        _ => unreachable!(),
                    }))
                }
                NetworkNode::Op(_) => Ok(None),
                _ => unreachable!(),
            })
            .unwrap();
        assert_eq!(factored, (&x + Atom::num(2)) * (&x - &y));
        let atoms = [&x, &y];
        let rows = tape.take_terms();
        assert_eq!(rows.len(), 4);
        let actual = Atom::add_many(rows.into_iter().map(|(factors, coefficient)| {
            Atom::num(coefficient)
                * Atom::mul_many(
                    factors
                        .into_iter()
                        .map(|(factor, power)| atoms[tape.leaves()[factor]].pow(Atom::num(power))),
                )
        }));
        assert_eq!(actual, ((&x + Atom::num(2)) * (&x - &y)).expand());
    }

    #[test]
    fn literal_definition_uses_the_same_tape_without_materializing_its_sum() {
        // Two literal uses share one already compiled (x+y) template. Its
        // alternatives participate in the surrounding product on this tape.
        let body = Graph::scalar(0).n_add([Graph::scalar(1)]);
        let tree = body.expr_tree().cast::<ChildVecStore<()>>();
        let root = body.graph.node_id(body.head());
        let mut tape = TermTape::default();
        let (definition, size) = tape
            .compile_graph(&body, &tree, root, &mut |_, _, node| match node {
                NetworkNode::Leaf(NetworkLeaf::Scalar(ScalarRef::Store(index))) => {
                    Some(TermLeaf::Value(*index))
                }
                _ => None,
            })
            .unwrap();
        let uses = Graph::scalar(2).n_mul([Graph::scalar(2)]);
        let tree = uses.expr_tree().cast::<ChildVecStore<()>>();
        let root = uses.graph.node_id(uses.head());
        let (product, _) = tape
            .compile_graph(&uses, &tree, root, &mut |_, _, node| {
                matches!(node, NetworkNode::Leaf(_))
                    .then_some(TermLeaf::Reference(definition, size))
            })
            .unwrap();
        tape.distribute(&mut vec![(product, 1)], &Rational::one(), 0)
            .unwrap();
        let rows = tape.take_terms();
        assert_eq!(tape.leaves().len(), 2);
        assert_eq!(rows.len(), 3);
        assert_eq!(
            rows.iter().map(|(_, value)| value).sum::<Rational>(),
            Rational::from(4)
        );
        assert!(
            rows.iter()
                .any(|(powers, value)| powers.len() == 2 && *value == 2)
        );

        let invalid = Graph::scalar(0);
        let tree = invalid.expr_tree().cast::<ChildVecStore<()>>();
        let root = invalid.graph.node_id(invalid.head());
        assert!(
            TermTape::<usize>::default()
                .compile_graph(&invalid, &tree, root, &mut |_, _, _| {
                    Some(TermLeaf::Reference(0, (1, 1)))
                })
                .is_none()
        );
    }

    #[test]
    fn opaque_graph_scope_does_not_observe_its_function_payload() {
        let hidden = Graph::scalar(0)
            .n_add([Graph::scalar(1)])
            .function(symbol!("term_tape_test::foreign"));
        let graph = hidden.n_mul([Graph::scalar(2)]);
        let tree = graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = graph.graph.node_id(graph.head());
        let mut observed = Vec::new();
        let mut tape = TermTape::default();
        let (node, _) = tape
            .compile_graph(&graph, &tree, root, &mut |_, id, data| {
                observed.push(id);
                match data {
                    NetworkNode::Op(NetworkOp::Function(_)) => Some(TermLeaf::Value(3usize)),
                    NetworkNode::Leaf(NetworkLeaf::Scalar(ScalarRef::Store(index))) => {
                        Some(TermLeaf::Value(*index))
                    }
                    _ => None,
                }
            })
            .unwrap();
        assert_eq!(
            observed.len(),
            3,
            "root, opaque function, and visible scalar only"
        );
        tape.distribute(&mut vec![(node, 1)], &Rational::one(), 0)
            .unwrap();
        assert_eq!(tape.take_terms().len(), 1);
        // Follow the graph's established product traversal order.
        assert_eq!(tape.leaves(), &[2, 3]);
    }

    #[test]
    fn graph_powers_coalesce_and_refuse_excessive_distribution() {
        let sum = Graph::scalar(0).n_add([Graph::scalar(1)]);
        let graph = sum.clone().pow(2);
        let tree = graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = graph.graph.node_id(graph.head());
        let leaf = |_: &mut TermTape<usize>,
                    _: NodeIndex,
                    node: &NetworkNode<i8, Symbol, AbstractIndex>| match node {
            NetworkNode::Leaf(NetworkLeaf::Scalar(ScalarRef::Store(index))) => {
                Some(TermLeaf::Value(*index))
            }
            _ => None,
        };
        let mut tape = TermTape::default();
        let (node, size) = tape
            .compile_graph(&graph, &tree, root, &mut { leaf })
            .unwrap();
        assert_eq!(size, (4, 8));
        tape.distribute(&mut vec![(node, 1)], &Rational::one(), 0)
            .unwrap();
        let rows = tape.take_terms();
        assert_eq!(rows.len(), 3);
        assert_eq!(
            rows.iter().map(|(_, c)| c).sum::<Rational>(),
            Rational::from(4)
        );
        assert!(rows.iter().any(|(powers, c)| powers.len() == 2 && *c == 2));
        assert!(tape.distributed());
        tape.clear_terms();
        assert!(tape.take_terms().is_empty());

        for exponent in [-1, 32] {
            let graph = sum.clone().pow(exponent);
            let tree = graph.expr_tree().cast::<ChildVecStore<()>>();
            let root = graph.graph.node_id(graph.head());
            let mut tape = TermTape::default();
            assert!(
                tape.compile_graph(&graph, &tree, root, &mut { leaf })
                    .is_none()
            );
            assert!(tape.take_terms().is_empty());
        }
    }
}
