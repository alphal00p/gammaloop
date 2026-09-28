//! State merging on the collector's existing port graph.
//!
//! Remaining factors are opaque terminals of the same component reducer used
//! for a completed monomial. Only the selected factor contributes alternatives;
//! substitutions remain port bindings until a distinct output variable is emitted.

use super::*;
use crate::tensor::SymbolicNet;
use linnet::{
    half_edge::{
        NodeIndex,
        subgraph::{Inclusion, SuBitGraph},
        tree::SimpleTraversalTree,
    },
    tree::child_vec::ChildVecStore,
};
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        library::DummyLibrary,
        parsing::{ParseSettings, ShorthandParsing},
        store::TensorScalarStore,
    },
    structure::{
        OrderedStructure, TensorStructure, abstract_index::AbstractIndex, partial::PartialIndex,
        slot::IsAbstractSlot,
    },
};
use symbolica::{atom::AtomOrView, coefficient::CoefficientView};

// A state contains occurrence-local port bindings, not rebuilt expressions.
// Variables and the incidence scratch still belong to ComponentSum.
#[derive(Clone, PartialEq, Eq, Hash)]
struct Remaining<'a> {
    factors: Vec<(usize, Vec<Argument<'a>>)>,
    residual: Vec<(usize, u16)>,
}

type FactorTerm = (Vec<(usize, u16)>, Rational);

/// Partial output must distinguish deferred work from an exhausted frontier.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum ContractionStatus {
    Complete,
    Deferred,
    Capped,
}

/// Exact output may retain pending contractions when a frontier budget is met.
/// Algebraic validity and completion are separate facts.
pub(crate) struct FactorizedContraction<Root = Atom> {
    pub(crate) root: Root,
    pub(crate) aliases: Vec<(Atom, Atom)>,
    pub(crate) status: ContractionStatus,
    pub(crate) literal_relabellings: Vec<(Atom, Atom)>,
}

impl SlotContraction {
    /// Plan an already admitted interface; the typed caller owns source validation.
    pub(crate) fn contract_factorized(
        &self,
        source: AtomView<'_>,
        order: Option<&[usize]>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        if source.needs_normalization()
            || self.metric.get_evaluation_info().is_some()
            || !InterfaceInference::normalization_is_intrinsic(source)
        {
            return None;
        }
        // The existing partial parser is the topology and layout owner. A
        // product exposes its immediate factors. Arithmetic is parsed once;
        // the term tape opens only the selected factor, retaining other
        // subtrees as occurrence-local graph leaves.
        let settings = ParseSettings {
            precontract_scalars: false,
            depth_limit: None,
            shorthand_parsing: ShorthandParsing::Opaque,
            parse_composite_scalars_as_tensors: true,
            ..ParseSettings::default()
        };
        type Tensor<'src> = SymbolicTensor<OrderedStructure, AtomOrView<'src>>;
        let network = SymbolicNet::<AbstractIndex, AtomOrView<'_>>::try_from_view::<
            OrderedStructure,
            _,
        >(source, &DummyLibrary::<Tensor<'_>>::new(), &settings)
        .ok()?;
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        self.contract_graph(&network, &tree, root, order, Some(source), rank_one)
    }

    fn power_scope(
        &self,
        network: &SymbolicNet<AbstractIndex, AtomOrView<'_>>,
        tree: &SimpleTraversalTree<ChildVecStore<()>>,
        node: NodeIndex,
    ) -> Option<(NodeIndex, i8, bool)> {
        let NetworkNode::Op(NetworkOp::Power(power)) = network.graph.graph[node] else {
            return None;
        };
        let mut children = tree.iter_children(node, &network.graph.graph);
        let child = children.next()?;
        if children.next().is_some() {
            return None;
        }
        if power <= 0 {
            return Some((child, power, true));
        }
        let members = tree
            .iter_preorder_tree_nodes(&network.graph.graph, child)
            .collect::<std::collections::BTreeSet<_>>();
        let mut internal = false;
        for &member in &members {
            for hedge in network.graph.graph.iter_crown(member) {
                if !network.graph.graph[[&hedge]].is_slot() {
                    continue;
                }
                if network.graph.graph.inv(hedge) == hedge
                    || !members
                        .contains(&network.graph.graph.node_id(network.graph.graph.inv(hedge)))
                {
                    return None;
                }
                internal = true;
            }
        }
        // Already compact scalar coefficients keep their existing tape keys,
        // so equal coefficients still merge across alpha-equivalent terms.
        (internal || matches!(network.graph.graph[child], NetworkNode::Op(_)))
            .then_some((child, power, internal))
    }

    fn emit_graph(
        &self,
        network: &SymbolicNet<AbstractIndex, AtomOrView<'_>>,
        tree: &SimpleTraversalTree<ChildVecStore<()>>,
        root: NodeIndex,
        scoped: &AHashMap<NodeIndex, FactorizedContraction>,
    ) -> Option<Atom> {
        network
            .graph
            .to_expression_at(tree, root, &mut |node, value| {
                if let Some(result) = scoped.get(&node) {
                    return Ok(Some(result.root.clone()));
                }
                match value {
                    NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => Ok(Some(
                        network.store.tensors[*index]
                            .expression
                            .as_view()
                            .to_owned(),
                    )),
                    NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => Ok(Some(
                        network.store.get_scalar_ref(*index).as_view().to_owned(),
                    )),
                    NetworkNode::Op(_) => Ok(None),
                    _ => Err(spenso::network::TensorNetworkError::Other(eyre::eyre!(
                        "symbolic contraction requires a stored occurrence"
                    ))),
                }
            })
            .ok()
    }

    fn contract_graph(
        &self,
        network: &SymbolicNet<AbstractIndex, AtomOrView<'_>>,
        tree: &SimpleTraversalTree<ChildVecStore<()>>,
        root: NodeIndex,
        order: Option<&[usize]>,
        original: Option<AtomView<'_>>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        // A scalar base owns its internal dummy pairs independently of the
        // exponent. Reduce it once before applying the power: distributing
        // first would identify pairs from different copies. Open positive
        // powers retain the existing cross-copy slot contraction.
        if let Some((child, power, internal)) = self.power_scope(network, tree, root) {
            if !internal {
                return Some(FactorizedContraction {
                    root: match original {
                        Some(source) => source.to_owned(),
                        None => self.emit_graph(network, tree, root, &AHashMap::new())?,
                    },
                    aliases: Vec::new(),
                    status: ContractionStatus::Complete,
                    literal_relabellings: Vec::new(),
                });
            }
            let mut result = self.contract_graph(network, tree, child, None, None, rank_one)?;
            result.root = result.root.pow(power);
            return Some(result);
        }
        let nodes = if matches!(
            network.graph.graph[root],
            NetworkNode::Op(NetworkOp::Product)
        ) {
            // Product input heads were connected from the end by the existing
            // parser. Public factor ordinals follow the normalized source.
            let mut children = tree
                .iter_children(root, &network.graph.graph)
                .collect::<Vec<_>>();
            children.reverse();
            children
        } else {
            vec![root]
        };
        let traversal = tree
            .iter_preorder_tree_nodes(&network.graph.graph, root)
            .collect::<Vec<_>>();
        let mut owned = AHashMap::new();
        for &node in traversal.iter().rev() {
            let contains = match &network.graph.graph[node] {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    matches!(network.store.tensors[*index].expression.as_view(), AtomView::Fun(function)
                        if function.get_symbol() == self.metric || function.get_symbol().has_tag(&self.tags.rank1))
                        || network.graph.graph.iter_crown(node).any(|hedge| {
                            network.graph.graph[[&hedge]].is_slot()
                                && network.graph.graph.inv(hedge) != hedge
                                && network.graph.graph.node_id(network.graph.graph.inv(hedge))
                                    == node
                        })
                }
                NetworkNode::Op(NetworkOp::Sum | NetworkOp::Product | NetworkOp::Power(_)) => {
                    let children = tree
                        .iter_children(node, &network.graph.graph)
                        .collect::<Vec<_>>();
                    children.iter().any(|child| owned[child])
                        || (matches!(
                            network.graph.graph[node],
                            NetworkNode::Op(NetworkOp::Product)
                        ) && network
                            .graph
                            .slot_components(tree, &children)
                            .iter()
                            .any(|part| part.len() > 1))
                }
                _ => false,
            };
            owned.insert(node, contains);
        }
        let mut scoped = AHashMap::new();
        let mut pending = vec![root];
        while let Some(node) = pending.pop() {
            if !owned[&node] {
                continue;
            }
            if self.power_scope(network, tree, node).is_some() {
                scoped.insert(
                    node,
                    self.contract_graph(network, tree, node, None, None, rank_one)?,
                );
            } else {
                pending.extend(tree.iter_children(node, &network.graph.graph));
            }
        }
        let mut selected = nodes
            .iter()
            .enumerate()
            .filter_map(|(i, node)| owned[node].then_some(i))
            .collect::<Vec<_>>();
        if let Some(order) = order {
            let mut seen = vec![false; nodes.len()];
            for &position in order {
                if *seen.get(position)? {
                    return None;
                }
                seen[position] = true;
            }
            if seen.iter().any(|entry| !entry) {
                return None;
            }
            selected
                .sort_by_key(|position| order.iter().position(|entry| entry == position).unwrap());
        }
        // Top-level factors and maximal unselected composite subtrees share
        // one occurrence table. Leaves keep their existing borrowed payload.
        let mut opaque_nodes = nodes.clone();
        let mut opaque_positions = nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(i, n)| (n, i))
            .collect::<AHashMap<_, _>>();
        let mut pending = selected.iter().map(|i| nodes[*i]).collect::<Vec<_>>();
        while let Some(node) = pending.pop() {
            if !owned[&node] || scoped.contains_key(&node) {
                if matches!(network.graph.graph[node], NetworkNode::Op(_))
                    && !opaque_positions.contains_key(&node)
                {
                    opaque_positions.insert(node, opaque_nodes.len());
                    opaque_nodes.push(node);
                }
            } else {
                pending.extend(tree.iter_children(node, &network.graph.graph));
            }
        }
        let external: SuBitGraph = network.graph.graph.external_filter();
        let mut port_atoms = Vec::new();
        let mut interfaces = Vec::new();
        for &node in &opaque_nodes {
            let members = tree
                .iter_preorder_tree_nodes(&network.graph.graph, node)
                .collect::<std::collections::BTreeSet<_>>();
            let mut ports = Vec::new();
            for &member in &members {
                for hedge in network.graph.graph.iter_crown(member) {
                    let NetworkEdge::Slot(slot) = network.graph.graph[[&hedge]] else {
                        continue;
                    };
                    if !external.includes(&hedge)
                        && members
                            .contains(&network.graph.graph.node_id(network.graph.graph.inv(hedge)))
                    {
                        continue;
                    }
                    // A sewn directed edge stores one shared label. The leaf's
                    // existing storage axis retains this endpoint's variance.
                    let slot = match &network.graph.graph[member] {
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                            *network.store.tensors[*index]
                                .structure
                                .external_structure()
                                .get(usize::from(network.graph.slot_order[hedge.0]))?
                        }
                        _ => slot,
                    };
                    ports.push(slot);
                }
            }
            // These are fully explicit internal factor interfaces. Port
            // bindings use the matching literal slot, not a fresh Atom scan.
            ports.sort();
            port_atoms.push(ports.iter().map(|slot| slot.to_atom()).collect::<Vec<_>>());
            interfaces.push(PartialStructure::from_logical_slots(
                ports
                    .into_iter()
                    .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
            ));
        }
        let emit = |node: usize| self.emit_graph(network, tree, NodeIndex(node), &scoped);
        let positions = nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(i, n)| (n, i))
            .collect::<AHashMap<_, _>>();
        let components = network.graph.slot_components(tree, &nodes);
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(self, Intake::Contraction, &mut slots);
        sum.rank_one = rank_one;
        sum.factor_emission = Some(&emit);
        sum.factor_roots.resize(opaque_nodes.len(), None);
        let mut remaining = Vec::new();
        for (position, &node) in opaque_nodes.iter().enumerate() {
            let ports = port_atoms[position]
                .iter()
                .map(Atom::as_view)
                .collect::<Vec<_>>();
            remaining.push((
                position,
                ports
                    .iter()
                    .copied()
                    .map(Argument::Original)
                    .collect::<Vec<_>>(),
            ));
            sum.opaque_factors
                .push((node.0, ports, interfaces[position].clone()));
        }
        // Complete admission precedes distribution/cancellation. Foreign and
        // scalar subtrees remain graph occurrences until an output is emitted.
        for &position in &selected {
            let (tape, _) = sum.input.compile_graph(
                &network.graph,
                tree,
                nodes[position],
                &mut |_, node, value| {
                    if let Some(result) = scoped.get(&node)
                        && let Some(&position) = opaque_positions.get(&node)
                        && port_atoms[position].is_empty()
                    {
                        return Some(TermLeaf::Value(InputLeaf::Scalar(result.root.as_view())));
                    }
                    let literal = match value {
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                            Some(network.store.tensors[*index].expression.as_view())
                        }
                        NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                            Some(network.store.get_scalar_ref(*index).as_view())
                        }
                        _ => None,
                    };
                    if let Some(value) = literal {
                        if let AtomView::Num(number) = value {
                            if let Ok(coefficient) = Rational::try_from(value) {
                                Some(TermLeaf::Number(InputLeaf::Atom(value), coefficient))
                            } else if matches!(
                                number.get_coeff_view(),
                                CoefficientView::Natural(..) | CoefficientView::Large(..)
                            ) {
                                // The tape's numeric field is real rational. Exact
                                // Gaussian coefficients remain opaque scalar leaves.
                                Some(TermLeaf::Value(InputLeaf::Scalar(value)))
                            } else {
                                None
                            }
                        } else {
                            Some(TermLeaf::Value(InputLeaf::Atom(value)))
                        }
                    } else if !owned[&node] || scoped.contains_key(&node) {
                        opaque_positions
                            .get(&node)
                            .map(|&position| TermLeaf::Value(InputLeaf::Subtree(position)))
                    } else {
                        None
                    }
                },
            )?;
            sum.factor_roots[position] = Some(tape);
        }
        // Scalar spectators and disconnected foreign sums never enter
        // another component's states. Connectivity comes from the same graph
        // that owns the parsed occurrences, not a second map of Atom labels.
        let mut roots = Vec::new();
        let mut definitions = scoped
            .values()
            .flat_map(|value| value.aliases.iter().cloned())
            .collect::<Vec<_>>();
        let mut complete = scoped
            .values()
            .all(|value| value.status == ContractionStatus::Complete);
        for component in components {
            // Alpha representatives carry literal dummy labels. They may merge
            // alternatives of this component, but cannot be reused in another
            // disconnected component whose contractions have independent labels.
            sum.alpha_terms.clear();
            let component = component
                .into_iter()
                .map(|node| positions[&node])
                .collect::<Vec<_>>();
            let selected = selected
                .iter()
                .copied()
                .filter(|position| component.contains(position))
                .collect::<Vec<_>>();
            if selected.is_empty() {
                roots.push(Atom::mul_many(
                    component
                        .iter()
                        .map(|&position| sum.emit_factor(position, &remaining[position].1))
                        .collect::<Option<Vec<_>>>()?,
                ));
                continue;
            }
            let result = sum.contract_states(
                Remaining {
                    factors: component
                        .into_iter()
                        .map(|position| remaining[position].clone())
                        .collect(),
                    residual: Vec::new(),
                },
                &selected,
            )?;
            roots.push(result.root);
            definitions.extend(result.aliases);
            complete &= result.status == ContractionStatus::Complete;
        }
        // Admission and the component reducer still run, but an already
        // contracted expression must not gain another layer of weight aliases.
        // `contracted` accumulates across all states; resetting incidence scratch
        // deliberately does not reset it. Scoped powers retain their own results.
        if complete
            && scoped.is_empty()
            && !sum.contracted
            && sum
                .literal_relabellings
                .borrow()
                .iter()
                .all(|(source, target)| source == target)
        {
            return Some(FactorizedContraction {
                // Graph product emission may sort anonymous bracket factors.
                // The existing no-work proof permits retaining exact source syntax.
                root: match original {
                    Some(source) => source.to_owned(),
                    None => self.emit_graph(network, tree, root, &scoped)?,
                },
                aliases: Vec::new(),
                status: ContractionStatus::Complete,
                literal_relabellings: Vec::new(),
            });
        }
        Some(FactorizedContraction {
            root: Atom::mul_many(roots),
            aliases: definitions,
            status: if complete {
                ContractionStatus::Complete
            } else {
                ContractionStatus::Capped
            },
            literal_relabellings: scoped
                .values()
                .flat_map(|value| value.literal_relabellings.iter().cloned())
                .chain(std::mem::take(&mut *sum.literal_relabellings.borrow_mut()))
                .collect(),
        })
    }
}

impl<'a> ComponentSum<'a, '_> {
    fn opaque_factor(&mut self, position: usize, arguments: &[Argument<'a>]) -> Option<()> {
        let tensor = self.tensors.len();
        self.tensors
            .push((TensorSource::Factor(position), arguments.to_vec()));
        for (port, &argument) in arguments.iter().enumerate() {
            if let Argument::Original(slot) = argument {
                let (space, index) = self.resolve_endpoint(slot)?;
                let node = self.endpoint(slot, space, index)?;
                self.nodes[node]
                    .terminals
                    .push(Terminal::Tensor(tensor, port));
            }
        }
        Some(())
    }

    fn residual_factor(&mut self, variable: usize, exponent: u16) -> Option<()> {
        match self.variables[variable].clone() {
            Variable::Vector(head, slot) => {
                let (space, index) = self.resolve_endpoint(slot)?;
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            Variable::Metric([first, second]) => {
                self.metric(
                    Argument::Original(first),
                    Argument::Original(second),
                    exponent,
                )?;
            }
            Variable::Tensor(source, arguments) => {
                if exponent != 1 {
                    return None;
                }
                let tensor = self.tensors.len();
                for (position, &argument) in arguments.iter().enumerate() {
                    if let Argument::Original(slot) = argument
                        && matches!(self.slots.classify(slot), SlotMatch::Explicit(_))
                    {
                        let (space, index) = self.resolve_endpoint(slot)?;
                        let node = self.endpoint(slot, space, index)?;
                        self.nodes[node]
                            .terminals
                            .push(Terminal::Tensor(tensor, position));
                    }
                }
                self.tensors.push((source, arguments));
            }
            value => self.variable(value, exponent),
        }
        Some(())
    }

    fn reset_components(&mut self, coefficient: &Rational) {
        self.nodes.clear();
        self.occurrences.clear();
        self.factors.clear();
        self.tensors.clear();
        self.metrics = 0;
        self.coefficient = coefficient.clone();
    }

    fn factor_terms(&mut self, position: usize) -> Option<Vec<FactorTerm>> {
        let node = self.factor_roots.get(position).copied().flatten()?;
        self.input.clear_terms();
        self.input
            .distribute(&mut vec![(node, 1)], &Rational::one(), 0)?;
        Some(self.input.take_terms())
    }

    fn contract_states(
        &mut self,
        remaining: Remaining<'a>,
        order: &[usize],
    ) -> Option<FactorizedContraction> {
        let mut states = vec![(remaining, Atom::num(1))];
        let mut definitions = Vec::new();
        let mut emitted: AHashMap<usize, Atom> = AHashMap::new();
        let mut weights = AHashMap::<Atom, Atom>::new();
        let mut generated = 0usize;
        let mut definition_bytes = 0usize;
        let mut complete = true;
        for &selected in order {
            let terms = self.factor_terms(selected)?;
            // Predict growth of the whole frontier, not just the local
            // template's term iterator. This includes graph-state storage and
            // existing definitions, not a hard heap bound on Atom buffers.
            // Flush exact remaining factors before an exponential next level;
            // explicit materialization may expand them.
            let state_bytes = states
                .iter()
                .map(|(remaining, _)| {
                    std::mem::size_of::<Remaining<'a>>()
                        + remaining
                            .factors
                            .iter()
                            .map(|(_, arguments)| {
                                std::mem::size_of::<(usize, Vec<Argument<'a>>)>()
                                    + arguments.len() * std::mem::size_of::<Argument<'a>>()
                            })
                            .sum::<usize>()
                        + remaining.residual.len() * std::mem::size_of::<(usize, u16)>()
                })
                .max()
                .unwrap_or(0);
            let count = states.len().saturating_mul(terms.len());
            let predicted = count.saturating_mul(state_bytes.saturating_mul(4).saturating_add(256));
            const MAX_FRONTIER_BYTES: usize = if cfg!(test) {
                64 * 1024
            } else {
                64 * 1024 * 1024
            };
            if generated.saturating_add(count) > Self::MAX_GENERATED_TERMS
                || predicted.saturating_add(definition_bytes) > MAX_FRONTIER_BYTES
            {
                complete = false;
                break;
            }
            generated += count;
            let mut positions: AHashMap<Remaining<'a>, usize> = AHashMap::new();
            let mut next: Vec<(Remaining<'a>, Vec<Atom>)> = Vec::new();
            for (remaining, weight) in states {
                let source = remaining
                    .factors
                    .iter()
                    .find(|(position, _)| *position == selected)?;
                self.overrides = self.opaque_factors[selected]
                    .1
                    .iter()
                    .copied()
                    .zip(source.1.iter().copied())
                    .collect();
                for (factors, coefficient) in &terms {
                    self.reset_components(coefficient);
                    for &(factor, exponent) in factors {
                        let leaf = self.input.leaves()[factor];
                        match leaf {
                            InputLeaf::Scalar(value) => {
                                self.variable(Variable::Scalar(value), exponent);
                            }
                            InputLeaf::Atom(value) if self.scalar_spectator(value) => {
                                self.variable(Variable::Scalar(value), exponent);
                            }
                            InputLeaf::Atom(value) => self.factor(value, exponent)?,
                            InputLeaf::Subtree(position) => {
                                let ports = &self.opaque_factors[position].1;
                                if ports.is_empty() {
                                    self.variable(Variable::Subtree(position), exponent);
                                } else {
                                    if exponent != 1 {
                                        return None;
                                    }
                                    let arguments = ports
                                        .iter()
                                        .map(|&slot| self.argument(slot))
                                        .collect::<Vec<_>>();
                                    self.opaque_factor(position, &arguments)?;
                                }
                            }
                        }
                    }
                    for &(variable, exponent) in &remaining.residual {
                        self.residual_factor(variable, exponent)?;
                    }
                    for (position, arguments) in &remaining.factors {
                        if *position != selected {
                            self.opaque_factor(*position, arguments)?;
                        }
                    }
                    self.reduce_components()?;
                    let mut key = Remaining {
                        factors: Vec::new(),
                        residual: Vec::new(),
                    };
                    let mut coefficient = vec![Atom::num(self.coefficient.clone()), weight.clone()];
                    for &(variable, exponent) in &self.factors {
                        match &self.variables[variable] {
                            Variable::Tensor(TensorSource::Factor(position), arguments) => {
                                key.factors.push((*position, arguments.clone()));
                            }
                            value @ (Variable::Scalar(_)
                            | Variable::Subtree(_)
                            | Variable::Dot(_, _)) => {
                                let atom = if let Some(value) = emitted.get(&variable) {
                                    value.clone()
                                } else {
                                    let atom = self.emit_variable(value)?;
                                    emitted.insert(variable, atom.clone());
                                    atom
                                };
                                coefficient.push(atom.pow(exponent));
                            }
                            _ => key.residual.push((variable, exponent)),
                        }
                    }
                    key.factors.sort_unstable_by_key(|(position, _)| *position);
                    key.residual.sort_unstable();
                    let coefficient = Atom::mul_many(coefficient);
                    if coefficient.is_zero() {
                        continue;
                    }
                    if let Some(&position) = positions.get(&key) {
                        next[position].1.push(coefficient);
                    } else {
                        positions.insert(key.clone(), next.len());
                        next.push((key, vec![coefficient]));
                    }
                }
            }
            self.overrides.clear();
            states = Vec::with_capacity(next.len());
            for (key, coefficients) in next {
                let body = Atom::add_many(coefficients);
                if body.is_zero() {
                    continue;
                }
                // Different remaining states often carry the same coefficient.
                // Keep one literal definition for that value, and avoid aliases
                // whose body is just a number or an existing weight handle.
                if matches!(body.as_view(), AtomView::Num(_)) {
                    states.push((key, body));
                    continue;
                }
                if let Some(handle) = weights.get(&body) {
                    states.push((key, handle.clone()));
                    continue;
                }
                let body =
                    SymbolicTensor::checked_parts(body, PartialStructure::from_logical_slots([]))
                        .ok()?;
                let handle = body.alias_handle().ok()?;
                definition_bytes = definition_bytes
                    .saturating_add(body.expression.as_view().get_byte_size())
                    .saturating_add(handle.expression.as_view().get_byte_size());
                weights.insert(body.expression.clone(), handle.expression.clone());
                weights.insert(handle.expression.clone(), handle.expression.clone());
                definitions.push((handle.expression.clone(), body.expression));
                states.push((key, handle.expression));
            }
        }
        let mut terms = Vec::with_capacity(states.len());
        for (remaining, weight) in states {
            let mut factors = vec![weight];
            for (position, arguments) in remaining.factors {
                factors.push(self.emit_factor(position, &arguments)?);
            }
            for (variable, exponent) in remaining.residual {
                let atom = if let Some(value) = emitted.get(&variable) {
                    value.clone()
                } else {
                    let atom = self.emit_variable(&self.variables[variable])?;
                    emitted.insert(variable, atom.clone());
                    atom
                };
                factors.push(atom.pow(exponent));
            }
            terms.push(Atom::mul_many(factors));
        }
        Some(FactorizedContraction {
            root: Atom::add_many(terms),
            aliases: definitions,
            status: if complete {
                ContractionStatus::Complete
            } else {
                ContractionStatus::Capped
            },
            literal_relabellings: Vec::new(),
        })
    }

    pub(super) fn argument(&self, value: AtomView<'a>) -> Argument<'a> {
        self.overrides
            .iter()
            .find_map(|&(source, target)| (source == value).then_some(target))
            .unwrap_or(Argument::Original(value))
    }

    fn emit_argument(&self, argument: Argument<'a>) -> Atom {
        match argument {
            Argument::Original(value) => value.to_owned(),
            Argument::Vector(space, head) => {
                let (representation, dimension) = self.spaces[space];
                self.emit_vector(head, representation.to_symbolic([dimension]).as_view())
            }
        }
    }

    pub(super) fn emit_factor(&self, position: usize, arguments: &[Argument<'a>]) -> Option<Atom> {
        let (node, ports, interface) = &self.opaque_factors[position];
        let source = (self.factor_emission?)(*node)?;
        let replacements = ports
            .iter()
            .zip(arguments)
            .enumerate()
            .filter(|(_, (port, argument))| **argument != Argument::Original(**port))
            .map(|(position, (_, &argument))| (position, self.emit_argument(argument)))
            .collect::<std::collections::HashMap<_, _>>();
        if replacements.is_empty() {
            return Some(source);
        }
        let source = SymbolicTensor::from_normalized_parts(source, interface.clone());
        let expression =
            crate::tensor::composition::rewrite_interface_ports(&source, &replacements, &mut |source, result| {
                if matches!(source, AtomView::Fun(function) if function.get_symbol() == spenso::tensor_symbol!("idenso::tensor_alias")) {
                    self.literal_relabellings.borrow_mut().push((source.to_owned(), result.clone()));
                }
            }).ok()?;
        Some(crate::shorthands::schoonschip::DotNormalizer::run(
            expression.as_view(),
        ))
    }

    pub(super) fn metric(
        &mut self,
        first: Argument<'a>,
        second: Argument<'a>,
        exponent: u16,
    ) -> Option<()> {
        match (first, second) {
            (Argument::Vector(space, a), Argument::Vector(other, b)) => {
                if space != other {
                    return None;
                }
                self.dot(space, a, b, exponent);
            }
            (Argument::Vector(space, head), Argument::Original(slot))
            | (Argument::Original(slot), Argument::Vector(space, head)) => {
                let (other, index) = self.resolve_endpoint(slot)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    self.dot(space, head, head, exponent / 2);
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let node = self.endpoint(slot, space, index)?;
                    self.nodes[node].terminals.push(Terminal::Vector(head));
                }
            }
            (Argument::Original(first), Argument::Original(second)) => {
                let (space, first_index) = self.resolve_endpoint(first)?;
                let (other, second_index) = self.resolve_endpoint(second)?;
                if space != other {
                    return None;
                }
                if exponent > 1 {
                    if first_index == second_index {
                        return None;
                    }
                    self.factor(self.spaces[space].1, exponent / 2)?;
                    self.contracted = true;
                }
                if exponent % 2 == 1 {
                    let a = self.endpoint(first, space, first_index)?;
                    let b = self.endpoint(second, space, second_index)?;
                    let a = self.root(a);
                    let b = self.root(b);
                    self.metrics += 1;
                    self.nodes[a].metrics += 1;
                    self.nodes[b].parent = a;
                }
            }
        }
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::{input, setup};
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use ahash::AHashSet;
    use symbolica::atom::AliasedAtom;

    fn resolve(contracted: FactorizedContraction) -> Atom {
        let mut result = AliasedAtom::from(contracted.root);
        for (handle, body) in contracted.aliases {
            result.register_alias(handle, body);
        }
        result.into_inner()
    }

    #[test]
    fn negative_power_contracts_its_existing_graph_scope() {
        let contractor = setup();
        for (source, expected) in [
            (
                "spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-1",
                "spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^-1",
            ),
            (
                "x+spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-2",
                "x+spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^-2",
            ),
        ] {
            let source = input(source);
            let result = contractor
                .contract_factorized(source.as_view(), None, true)
                .unwrap();
            assert_eq!(resolve(result), input(expected), "{source}");
        }
    }

    #[test]
    fn factorized_states_retain_exact_complex_coefficients() {
        let contractor = setup();
        let source = input("1i*p(mink(4,a))*q(mink(4,a))+2i*r(mink(4,b))*s(mink(4,b))");
        let result = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(result.status == ContractionStatus::Complete);
        assert_eq!(
            resolve(result),
            input("1i*g(p(mink(4)),q(mink(4)))+2i*g(r(mink(4)),s(mink(4)))")
        );
        let inexact = input("1.25*p(mink(4,a))*q(mink(4,a))+2.5*r(mink(4,b))*s(mink(4,b))");
        assert!(
            contractor
                .contract_factorized(inexact.as_view(), None, true)
                .is_none()
        );
    }

    #[test]
    fn factorized_states_use_the_component_reducer_in_both_orders() {
        let contractor = setup();
        for source in [
            "g(mink(4,a),mink(4,b))*(p(mink(4,a))+q(mink(4,a)))*r(mink(4,b))",
            "(g(mink(4,a),mink(4,b))*p(mink(4,c))+g(mink(4,a),mink(4,c))*q(mink(4,b)))*(p(mink(4,a))*q(mink(4,b))*r(mink(4,c))+q(mink(4,a))*r(mink(4,b))*s(mink(4,c)))",
            "(g(mink(4,a),mink(4,b))+p(mink(4,a))*q(mink(4,b)))*(g(mink(4,a),mink(4,b))+r(mink(4,a))*s(mink(4,b)))",
            "(g(mink(4,a),mink(4,b))+t(mink(4,a),mink(4,b)))*p(mink(4,a))*q(mink(4,b))",
        ] {
            let source = input(source);
            let expected = source.expand().schoonschip().expand();
            let count = match source.as_view() {
                AtomView::Mul(product) => product.iter().len(),
                _ => 1,
            };
            let reversed = (0..count).rev().collect::<Vec<_>>();
            for order in [None, Some(reversed.as_slice())] {
                let contracted = contractor
                    .contract_factorized(source.as_view(), order, true)
                    .unwrap();
                let mut bodies = AHashSet::new();
                let handles = contracted
                    .aliases
                    .iter()
                    .map(|(handle, _)| handle)
                    .collect::<AHashSet<_>>();
                for (_, body) in &contracted.aliases {
                    assert!(bodies.insert(body));
                    assert!(!matches!(body.as_view(), AtomView::Num(_)));
                    assert!(!handles.contains(body));
                }
                assert_eq!(resolve(contracted).expand(), expected, "{source}");
            }
        }
    }

    #[test]
    fn factorized_states_keep_foreign_and_disconnected_sums() {
        let contractor = setup();
        let spectator = input("(1+x)*(t(mink(4,c))+routing(z)*u(mink(4,c)))");
        let first = input("(p(mink(4,a))+q(mink(4,a)))*r(mink(4,a))");
        let second = input("(p(mink(4,b))+s(mink(4,b)))*q(mink(4,b))");
        let source = Atom::mul_many([&spectator, &first, &second]);
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        let result = resolve(contracted);
        let expected = Atom::mul_many([
            spectator,
            first.expand().schoonschip(),
            second.expand().schoonschip(),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn factorized_port_bindings_remain_occurrence_local() {
        let contractor = setup();
        let source =
            input("(p(mink(4,a))+q(mink(4,a)))*(p(mink(4,b))+q(mink(4,b)))*t(mink(4,a),mink(4,b))");
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert_eq!(
            resolve(contracted).expand(),
            source.expand().schoonschip().expand()
        );
        assert!(
            contractor
                .contract_factorized(source.as_view(), Some(&[0, 0, 1]), true)
                .is_none()
        );
    }

    #[test]
    fn factorized_frontier_budget_retains_an_exact_unexpanded_remainder() {
        let contractor = setup();
        let slots = (0..9).map(|i| format!("mink(4,a{i})")).collect::<Vec<_>>();
        let source = input(&format!(
            "t({})*{}",
            slots.join(","),
            slots
                .iter()
                .map(|slot| format!("(x*p({slot})+y*q({slot}))"))
                .collect::<Vec<_>>()
                .join("*")
        ));
        let contracted = contractor
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(!contracted.aliases.is_empty());
        assert!(contracted.status == ContractionStatus::Capped);
        use spenso::network::parsing::AtomStructureExt;
        assert!(
            contracted.root.has_repeated_explicit_indices(),
            "the bounded frontier must keep its unprocessed contractions"
        );
        let result = resolve(contracted);
        assert_eq!(
            result.expand().schoonschip().expand(),
            source.expand().schoonschip().expand()
        );
    }
    #[test]
    fn scalar_power_without_graph_incidence_keeps_its_original_source() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbolica::symbol!("no_work_power_x"));
        let y = Atom::var(symbolica::symbol!("no_work_power_y"));
        let z = Atom::var(symbolica::symbol!("no_work_power_z"));
        let source = spenso::bracket!(x + y, z).pow(3);
        let result = SlotContraction::new()
            .contract_factorized(source.as_view(), None, true)
            .unwrap();
        assert!(result.status == ContractionStatus::Complete);
        assert!(result.aliases.is_empty());
        assert_eq!(result.root, source);
    }
}
