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
#[derive(Clone)]
pub(crate) struct FactorizedContraction<Root = Atom> {
    pub(crate) root: Root,
    pub(crate) aliases: Vec<(Atom, Atom)>,
    pub(crate) status: ContractionStatus,
    pub(crate) literal_relabellings: Vec<(Atom, Atom)>,
}

// A parse belongs to one arithmetic scope. Opening an occurrence is cached on
// that scope, so borrowed payloads and its port atoms outlive every DP state.
// This is parser scratch, not another symbolic tensor representation.
struct ContractionScope<'src> {
    source: AtomView<'src>,
    network: SymbolicNet<AbstractIndex, AtomOrView<'src>>,
    tree: SimpleTraversalTree<ChildVecStore<()>>,
    root: NodeIndex,
    ports: AHashMap<NodeIndex, (Vec<Atom>, PartialStructure)>,
    opened: AHashMap<NodeIndex, std::cell::OnceCell<Option<Box<Self>>>>,
    powers: AHashMap<NodeIndex, std::cell::OnceCell<Option<FactorizedContraction>>>,
    work: AHashMap<NodeIndex, std::cell::OnceCell<bool>>,
}

impl<'src> ContractionScope<'src> {
    fn parse(source: AtomView<'src>) -> Option<Self> {
        let settings = ParseSettings {
            precontract_scalars: true,
            depth_limit: Some(1),
            shorthand_parsing: ShorthandParsing::Opaque,
            parse_composite_scalars_as_tensors: true,
            ..ParseSettings::default()
        };
        type Tensor<'src> = SymbolicTensor<OrderedStructure, AtomOrView<'src>>;
        let network = SymbolicNet::<AbstractIndex, AtomOrView<'src>>::try_from_view::<
            OrderedStructure,
            _,
        >(source, &DummyLibrary::<Tensor<'src>>::new(), &settings)
        .ok()?;
        let tree = network.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = network.graph.graph.node_id(network.graph.head());
        let nodes = tree
            .iter_preorder_tree_nodes(&network.graph.graph, root)
            .collect::<Vec<_>>();
        let external: SuBitGraph = network.graph.graph.external_filter();
        let mut ports = AHashMap::new();
        for &node in &nodes {
            let members = tree
                .iter_preorder_tree_nodes(&network.graph.graph, node)
                .collect::<std::collections::BTreeSet<_>>();
            let mut boundary = Vec::new();
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
                    // Sewing shares one edge label; the occurrence's storage
                    // axis still owns this endpoint's variance.
                    let slot = match &network.graph.graph[member] {
                        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                            *network.store.tensors[*index]
                                .structure
                                .external_structure()
                                .get(usize::from(network.graph.slot_order[hedge.0]))?
                        }
                        _ => slot,
                    };
                    boundary.push(slot);
                }
            }
            boundary.sort();
            ports.insert(
                node,
                (
                    boundary.iter().map(|slot| slot.to_atom()).collect(),
                    PartialStructure::from_logical_slots(
                        boundary
                            .into_iter()
                            .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
                    ),
                ),
            );
        }
        Some(Self {
            source,
            network,
            tree,
            root,
            ports,
            opened: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            powers: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
            work: nodes
                .iter()
                .map(|&node| (node, Default::default()))
                .collect(),
        })
    }

    fn payload(&self, node: NodeIndex) -> Option<&AtomOrView<'src>> {
        match self.network.graph.graph[node] {
            NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                Some(&self.network.store.tensors[index].expression)
            }
            NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                Some(self.network.store.get_scalar_ref(index))
            }
            _ => None,
        }
    }

    fn literal(&self, node: NodeIndex) -> Option<AtomView<'_>> {
        self.payload(node)
            .map(AtomOrView::as_view)
            .or_else(|| (node == self.root).then_some(self.source))
    }

    fn open(&self, node: NodeIndex) -> Option<&Self> {
        self.opened[&node]
            .get_or_init(|| {
                // Composite arithmetic leaves retain their original view. The
                // parser's owned scalar coefficients contain no tensor interior.
                let AtomOrView::View(value) = self.payload(node)? else {
                    return None;
                };
                let opened = Self::parse(*value)?;
                if matches!(
                    opened.network.graph.graph[opened.root],
                    NetworkNode::Leaf(_)
                ) && opened.literal(opened.root)? == *value
                {
                    // Unsupported powers and opaque scalar syntax do not expose
                    // arithmetic. Decline instead of reopening the same leaf.
                    return None;
                }
                Some(Box::new(opened))
            })
            .as_deref()
    }

    fn contains_work(&self, node: NodeIndex, contractor: &SlotContraction) -> bool {
        *self.work[&node].get_or_init(|| {
            let Some(value) = self.literal(node) else {
                return false;
            };
            let observed = SimplificationCandidates::scan(value, [contractor.metric], || true);
            observed.symbols[0] || observed.dots || observed.repeated_indices
        })
    }

    fn factors(&self, root: NodeIndex) -> Vec<NodeIndex> {
        if matches!(
            self.network.graph.graph[root],
            NetworkNode::Op(NetworkOp::Product)
        ) {
            let mut nodes = self
                .tree
                .iter_children(root, &self.network.graph.graph)
                .collect::<Vec<_>>();
            nodes.reverse();
            nodes
        } else {
            vec![root]
        }
    }
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
        let scope = ContractionScope::parse(source)?;
        self.contract_graph(&scope, scope.root, order, rank_one)
    }

    fn power_scope(
        &self,
        scope: &ContractionScope<'_>,
        node: NodeIndex,
    ) -> Option<(NodeIndex, i8, bool)> {
        let NetworkNode::Op(NetworkOp::Power(power)) = scope.network.graph.graph[node] else {
            return None;
        };
        let mut children = scope.tree.iter_children(node, &scope.network.graph.graph);
        let child = children.next()?;
        if children.next().is_some() || !scope.ports[&child].0.is_empty() {
            return None;
        }
        let source = scope.literal(child)?;
        use spenso::network::parsing::AtomStructureExt;
        Some((child, power, source.has_repeated_explicit_indices()))
    }

    fn scoped_power<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        rank_one: bool,
    ) -> Option<&'a FactorizedContraction> {
        let (child, power, internal) = self.power_scope(scope, node)?;
        scope.powers[&node]
            .get_or_init(|| {
                if !internal {
                    // A transparent bracket can remain after local dot pairing
                    // has removed its last explicit index. Lower just that
                    // selected wrapper, retaining ordinary scalar powers verbatim.
                    if matches!(scope.literal(child)?, AtomView::Fun(function)
                        if function.get_symbol() == self.tags.bracket)
                        && scope.contains_work(child, self)
                    {
                        let opened = scope.open(child)?;
                        let mut result =
                            self.contract_graph(opened, opened.root, None, rank_one)?;
                        result.root = result.root.pow(power);
                        return Some(result);
                    }
                    return Some(FactorizedContraction {
                        root: scope.literal(node)?.to_owned(),
                        aliases: Vec::new(),
                        status: ContractionStatus::Complete,
                        literal_relabellings: Vec::new(),
                    });
                }
                // Contract a closed base once in its own scope before applying the
                // exponent. Copies must never share their internal dummy pairs.
                let mut result = self.contract_graph(scope, child, None, rank_one)?;
                result.root = result.root.pow(power);
                Some(result)
            })
            .as_ref()
    }

    fn factor_order(
        &self,
        scope: &ContractionScope<'_>,
        nodes: &[NodeIndex],
        mut selected: Vec<usize>,
        order: Option<&[usize]>,
    ) -> Option<Vec<usize>> {
        if let Some(order) = order {
            let factors = match scope.source {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                value => vec![value],
            };
            let mut positions = order.to_vec();
            positions.sort_unstable();
            if positions != (0..factors.len()).collect::<Vec<_>>() {
                return None;
            }
            let ranks = selected
                .iter()
                .map(|&position| {
                    let literal = scope.literal(nodes[position])?;
                    let ordinal = factors.iter().position(|&factor| factor == literal)?;
                    Some((position, order.iter().position(|&entry| entry == ordinal)?))
                })
                .collect::<Option<AHashMap<_, _>>>()?;
            selected.sort_by_key(|position| ranks[position]);
            return Some(selected);
        }
        // Minimize the live boundary of the accumulated factor set. Incidence
        // comes entirely from the depth-one graph; scoring never opens a leaf.
        let graph = &scope.network.graph.graph;
        let edges = nodes
            .iter()
            .map(|&node| {
                scope
                    .tree
                    .iter_preorder_tree_nodes(graph, node)
                    .flat_map(|member| graph.iter_crown(member))
                    .filter(|hedge| scope.network.graph.graph[[hedge]].is_slot())
                    .filter(|&hedge| {
                        graph.inv(hedge) == hedge || graph.node_id(graph.inv(hedge)) != node
                    })
                    .map(|hedge| hedge.0.min(graph.inv(hedge).0))
                    .collect::<std::collections::BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        let mut boundary = std::collections::BTreeSet::new();
        let mut planned = Vec::with_capacity(selected.len());
        while !selected.is_empty() {
            let best = (0..selected.len()).min_by_key(|&i| {
                let position = selected[i];
                (
                    boundary.symmetric_difference(&edges[position]).count(),
                    std::cmp::Reverse(boundary.intersection(&edges[position]).count()),
                    position,
                )
            })?;
            let position = selected.remove(best);
            boundary = boundary
                .symmetric_difference(&edges[position])
                .copied()
                .collect();
            planned.push(position);
        }
        Some(planned)
    }

    fn compile_scope<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        tape: &mut TermTape<InputLeaf<'a>>,
        sum: &mut ComponentSum<'a, '_>,
        positions: &mut AHashMap<(usize, NodeIndex), usize>,
        scoped: &mut AHashMap<(usize, NodeIndex), &'a FactorizedContraction>,
    ) -> Option<(usize, (usize, usize))> {
        tape.compile_graph(
            &scope.network.graph,
            &scope.tree,
            node,
            &mut |tape, node, value| {
                if let Some((child, power, _)) = self.power_scope(scope, node) {
                    let result = self.scoped_power(scope, node, sum.rank_one)?;
                    let unchanged = scope.literal(node).map_or_else(
                        || {
                            scope
                                .literal(child)
                                .is_some_and(|base| result.root == base.pow(power))
                        },
                        |source| result.root.as_view() == source,
                    );
                    // Visiting an unchanged scalar power is bookkeeping, not
                    // a contraction. It must not disable exact source reuse.
                    if !unchanged
                        || result.status != ContractionStatus::Complete
                        || !result.aliases.is_empty()
                        || result
                            .literal_relabellings
                            .iter()
                            .any(|(source, target)| source != target)
                    {
                        scoped.insert((scope as *const _ as usize, node), result);
                    }
                    return Some(TermLeaf::Value(InputLeaf::Scalar(result.root.as_view())));
                }
                if !matches!(value, NetworkNode::Leaf(_)) {
                    return None;
                }
                let literal = scope.literal(node)?;
                if matches!(
                    literal,
                    AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_)
                ) || matches!(literal, AtomView::Fun(function) if function.get_symbol() == self.tags.bracket) {
                    if scope.contains_work(node, self) {
                        let opened = scope.open(node)?;
                        let (node, size) =
                            self.compile_scope(opened, opened.root, tape, sum, positions, scoped)?;
                        Some(TermLeaf::Reference(node, size))
                    } else {
                        let position = sum.register_factor(scope, node, positions)?;
                        Some(TermLeaf::Value(InputLeaf::Subtree(position)))
                    }
                } else if let AtomView::Num(number) = literal {
                    if let Ok(coefficient) = Rational::try_from(literal) {
                        Some(TermLeaf::Number(InputLeaf::Atom(literal), coefficient))
                    } else if matches!(
                        number.get_coeff_view(),
                        CoefficientView::Natural(..) | CoefficientView::Large(..)
                    ) {
                        Some(TermLeaf::Value(InputLeaf::Scalar(literal)))
                    } else {
                        None
                    }
                } else {
                    Some(TermLeaf::Value(InputLeaf::Atom(literal)))
                }
            },
        )
    }

    fn contract_graph<'a>(
        &self,
        scope: &'a ContractionScope<'_>,
        root: NodeIndex,
        order: Option<&[usize]>,
        rank_one: bool,
    ) -> Option<FactorizedContraction> {
        if self.power_scope(scope, root).is_some() {
            return self.scoped_power(scope, root, rank_one).cloned();
        }
        let nodes = scope.factors(root);
        let selected = nodes
            .iter()
            .enumerate()
            .filter_map(|(i, &node)| scope.contains_work(node, self).then_some(i))
            .collect();
        let selected = self.factor_order(scope, &nodes, selected, order)?;
        let components = scope.network.graph.slot_components(&scope.tree, &nodes);
        let node_positions = nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(i, n)| (n, i))
            .collect::<AHashMap<_, _>>();
        let mut slots = SlotMatcher::default();
        let mut sum = ComponentSum::new(self, Intake::Contraction, &mut slots);
        sum.rank_one = rank_one;
        let mut positions = AHashMap::new();
        let mut remaining = Vec::new();
        for &node in &nodes {
            let position = sum.register_factor(scope, node, &mut positions)?;
            remaining.push((
                position,
                sum.opaque_factors[position]
                    .1
                    .iter()
                    .copied()
                    .map(Argument::Original)
                    .collect::<Vec<_>>(),
            ));
        }
        let mut scoped = AHashMap::new();
        let mut compile = |sum: &mut ComponentSum<'a, '_>, position: usize| {
            let mut tape = std::mem::take(&mut sum.input);
            let result = self.compile_scope(
                scope,
                nodes[position],
                &mut tape,
                sum,
                &mut positions,
                &mut scoped,
            );
            sum.input = tape;
            result.map(|(node, _)| node)
        };
        let mut roots = Vec::new();
        let mut definitions = Vec::new();
        let mut status = ContractionStatus::Complete;
        for component in components {
            sum.alpha_terms.clear();
            let component = component
                .into_iter()
                .map(|node| node_positions[&node])
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
                &mut compile,
            )?;
            roots.push(result.root);
            definitions.extend(result.aliases);
            status = status.max(result.status);
        }
        for result in scoped.values() {
            definitions.extend(result.aliases.iter().cloned());
            status = status.max(result.status);
        }
        if status == ContractionStatus::Complete
            && scoped.is_empty()
            && !sum.contracted
            && sum
                .literal_relabellings
                .borrow()
                .iter()
                .all(|(source, target)| source == target)
        {
            let source = scope.literal(root)?;
            // The bracket normalizer owns whether a product can be exposed:
            // closed scalar wrappers disappear, while AUTO operand order stays.
            let root = if matches!(source, AtomView::Fun(function)
                if function.get_symbol() == self.tags.bracket)
            {
                crate::shorthands::bracket::BracketNormalizer::normalize(source)
            } else {
                source.to_owned()
            };
            return Some(FactorizedContraction {
                root,
                aliases: Vec::new(),
                status,
                literal_relabellings: Vec::new(),
            });
        }
        Some(FactorizedContraction {
            root: Atom::mul_many(roots),
            aliases: definitions,
            status,
            literal_relabellings: scoped
                .values()
                .flat_map(|value| value.literal_relabellings.iter().cloned())
                .chain(std::mem::take(&mut *sum.literal_relabellings.borrow_mut()))
                .collect(),
        })
    }
}

impl<'a> ComponentSum<'a, '_> {
    fn register_factor(
        &mut self,
        scope: &'a ContractionScope<'_>,
        node: NodeIndex,
        positions: &mut AHashMap<(usize, NodeIndex), usize>,
    ) -> Option<usize> {
        let key = (scope as *const _ as usize, node);
        if let Some(&position) = positions.get(&key) {
            return Some(position);
        }
        let position = self.opaque_factors.len();
        let (ports, interface) = &scope.ports[&node];
        self.opaque_factors.push((
            scope.literal(node)?,
            ports.iter().map(Atom::as_view).collect(),
            interface.clone(),
        ));
        self.factor_roots.push(None);
        positions.insert(key, position);
        Some(position)
    }

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
        compile: &mut impl FnMut(&mut Self, usize) -> Option<usize>,
    ) -> Option<FactorizedContraction> {
        let mut states = vec![(remaining, Atom::num(1))];
        let mut definitions = Vec::new();
        let mut emitted: AHashMap<usize, Atom> = AHashMap::new();
        let mut weights = AHashMap::<Atom, Atom>::new();
        let mut generated = 0usize;
        let mut definition_bytes = 0usize;
        let mut complete = true;
        for &selected in order {
            if self.factor_roots[selected].is_none() {
                self.factor_roots[selected] = Some(compile(self, selected)?);
            }
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
        let (source, ports, interface) = &self.opaque_factors[position];
        let source = source.to_owned();
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

    #[test]
    fn scalar_rational_sums_keep_their_factorization() {
        let contractor = setup();
        for scalar in ["z", "q(0,spenso::cind(0))", "g(p(mink(4)),q(mink(4)))"] {
            let source = input(&format!(
                "n*(1/({scalar}+x)+1/({scalar}+y))+m*(1/({scalar}-x)+1/({scalar}-y))"
            ));
            for rank_one in [false, true] {
                let result = contractor
                    .contract_factorized(source.as_view(), None, rank_one)
                    .unwrap();
                assert!(result.status == ContractionStatus::Complete);
                assert!(result.aliases.is_empty());
                // An exact Atom comparison protects the factored scalar input;
                // expanding both sides would hide this regression.
                assert_eq!(result.root, source, "{scalar}, rank_one={rank_one}");

                let tensor = SymbolicTensor::infer(source.clone()).unwrap();
                let settings = crate::tensor::ContractionSettings::default();
                let settings = if rank_one {
                    settings
                } else {
                    settings.without_rank_one_tensors()
                };
                for settings in [settings, settings.with_order(&[0])] {
                    let result = tensor.contract(settings).unwrap();
                    assert_eq!(result.root(), tensor);
                    assert!(result.expression.get_aliases().is_empty());
                    assert_eq!(result.resolved().unwrap().expression, source);
                }
            }
        }
    }
}

#[cfg(test)]
mod lazy_tests {
    use super::*;
    use spenso::network::tags::SPENSO_TAG;

    #[test]
    fn partial_arithmetic_scopes_have_only_one_level() {
        let _contractor = super::super::tests::setup();
        for (source, nodes) in [
            ("p(spenso::mink(4,a))*q(spenso::mink(4,a))+x*y", 3),
            ("(p(spenso::mink(4,a))*q(spenso::mink(4,a))+x)^2", 2),
            (
                "spenso::bracket(p(spenso::mink(4,a))*q(spenso::mink(4,a)))^-1",
                2,
            ),
        ] {
            let source = super::super::tests::input(source);
            let scope = ContractionScope::parse(source.as_view()).unwrap();
            assert_eq!(scope.ports.len(), nodes, "{source}");
            assert_eq!(scope.network.graph.graph.n_nodes(), nodes, "{source}");
            assert!(scope.opened.values().all(|entry| entry.get().is_none()));
        }
    }

    #[test]
    fn fixed_component_scalar_sums_stay_closed() {
        let contractor = super::super::tests::setup();
        let source =
            super::super::tests::input("n*(1/(q(0,spenso::cind(0))+x)+1/(q(0,spenso::cind(0))+y))");
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let result = contractor
            .contract_graph(&scope, scope.root, None, true)
            .unwrap();
        assert!(result.status == ContractionStatus::Complete);
        assert_eq!(result.root, source);
        assert!(scope.opened.values().all(|entry| entry.get().is_none()));
    }

    #[test]
    fn graph_order_starts_with_a_small_boundary_and_accepts_an_override() {
        let contractor = super::super::tests::setup();
        let source = super::super::tests::input(
            "spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))*q(spenso::mink(4,b))",
        );
        let scope = ContractionScope::parse(source.as_view()).unwrap();
        let nodes = scope.factors(scope.root);
        let selected = (0..nodes.len()).collect::<Vec<_>>();
        let order = contractor
            .factor_order(&scope, &nodes, selected.clone(), None)
            .unwrap();
        assert_eq!(scope.ports[&nodes[order[0]]].0.len(), 1);
        let explicit = selected.iter().rev().copied().collect::<Vec<_>>();
        assert_eq!(
            contractor
                .factor_order(&scope, &nodes, selected, Some(&explicit))
                .unwrap(),
            explicit
        );
        let expected =
            super::super::tests::input("spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))");
        for order in [None, Some(explicit.as_slice())] {
            let result = contractor
                .contract_graph(&scope, scope.root, order, true)
                .unwrap();
            assert!(result.status == ContractionStatus::Complete);
            let mut resolved = symbolica::atom::AliasedAtom::from(result.root);
            for (handle, body) in result.aliases {
                resolved.register_alias(handle, body);
            }
            assert_eq!(resolved.into_inner(), expected);
        }
    }

    #[test]
    fn unselected_thousand_term_factors_stay_closed() {
        crate::test_support::test_initialize();
        SPENSO_TAG.rank_one_tensor_symbol("lazy_contraction::p");
        SPENSO_TAG.tensor_symbol("lazy_contraction::t");
        SPENSO_TAG.tensor_symbol("lazy_contraction::u");
        let parse = |source: &str| {
            Atom::parse(
                source,
                "lazy_contraction",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap()
        };
        for foreign in [false, true] {
            let first = Atom::add_many((0..1000).map(|i| {
                parse(&if foreign {
                    format!("t({i},spenso::mink(4,c))")
                } else {
                    format!("x({i})")
                })
            }));
            let second = Atom::add_many((0..1000).map(|i| {
                parse(&if foreign {
                    format!("u({i},spenso::mink(4,d))")
                } else {
                    format!("y({i})")
                })
            }));
            let core = parse("spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))");
            let source = Atom::mul_many([core, first.clone(), second.clone()]);
            let scope = ContractionScope::parse(source.as_view()).unwrap();
            assert_eq!(scope.ports.len(), 5, "one product and four factor leaves");
            assert_eq!(scope.network.graph.graph.n_nodes(), 5);
            let protected = scope
                .factors(scope.root)
                .into_iter()
                .filter(|&node| {
                    let value = scope.literal(node).unwrap();
                    value == first.as_view() || value == second.as_view()
                })
                .collect::<Vec<_>>();
            assert_eq!(protected.len(), 2);
            let result = SlotContraction::new()
                .contract_graph(&scope, scope.root, None, true)
                .unwrap();
            assert!(result.status == ContractionStatus::Complete);
            assert!(result.aliases.is_empty());
            assert_eq!(
                result.root,
                Atom::mul_many([parse("p(spenso::mink(4,b))"), first, second])
            );
            for node in protected {
                assert!(
                    scope.opened[&node].get().is_none(),
                    "unselected factor was opened"
                );
            }
            assert!(scope.opened.values().all(|opened| opened.get().is_none()));
        }
    }
}
