//! Signed tensor canonicalization on Spenso's symbolic network.
//!
//! Tensor syntax, operation grouping, and contraction topology belong to the
//! network parser. Graphica detects negative automorphisms of antisymmetric
//! slots, and Symbolica owns dummy labeling and expression rebuilding. Typed
//! callers carry source leaf facts through those known dummy substitutions.
use super::composition;
use super::inference::{InterfaceInference, LeafInference, TensorInferenceError};
use super::{SymbolicNet, SymbolicNetExt, SymbolicTensor, TensorProofs};
use crate::{CanonicalizationError, IndexTooling, NetworkToolingError};
use ahash::{AHashMap, AHashSet, AHasher};
use linnet::half_edge::NodeIndex;
use linnet::half_edge::involution::{Hedge, HedgePair};
use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph};
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::permutation::Permutation;
use linnet::tree::child_pointer::ParentChildStore;
use linnet::tree::child_vec::ChildVecStore;
use spenso::network::TensorNetworkError;
use spenso::network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp};
use spenso::network::library::DummyLibrary;
use spenso::network::parsing::{LeafInterfaceCache, ParseSettings, StrictTensorFilter};
use spenso::network::store::{NetworkStoreAccess, TensorScalarStore};
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, LibrarySlot, RepName, Representation};
use spenso::structure::slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind, Slot};
use spenso::structure::{HasName, OrderedStructure, TensorStructure};
use std::collections::{BTreeMap, BTreeSet};
use std::hash::{Hash, Hasher};
use std::sync::Arc;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::graph::{Graph, HiddenData};
use symbolica::id::{AliasedAtom, Pattern, Replacement};
use symbolica::symbol;

impl SymbolicTensor<PartialStructure> {
    /// Canonically label internal indices while retaining the ordered external interface.
    ///
    /// The shared network canonicalizer owns signed zeros and dummy labels. The
    /// tensor finisher preserves typed zeros. Debug builds check the user-normalizer
    /// free-index preservation contract.
    pub fn canonize(&self) -> Result<Self, CanonicalizationError> {
        self.canonize_with(&mut AbstractIndex::Dummy)
    }

    fn canonize_with<Aind: AbsInd + ParseableAind + DummyAind>(
        &self,
        new_dummy: &mut dyn FnMut(usize) -> Aind,
    ) -> Result<Self, CanonicalizationError> {
        if self.expression.is_zero() {
            return Ok(self.clone());
        }
        let (proofs, leaves) = self.producer_context();
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let network = SymbolicNet::<Aind>::try_from_admitted_view_with_context::<
            OrderedStructure<LibraryRep, Aind>,
            _,
        >(
            self.expression.as_view(),
            &DummyLibrary::<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>::new(),
            &settings,
            leaves,
            Some(&proofs),
        )
        .map_err(|error| NetworkToolingError::Parse {
            reason: error.to_string(),
        })?;
        let (expression, proofs) =
            Self::canonize_parsed(network, new_dummy, Some((self.structure.clone(), proofs)))?;
        let proofs = proofs.expect("admitted canonicalization retains its producer facts");
        Ok(self
            .borrowed_with_proofs(&proofs)
            .with_materialized_expression("canonize", expression)?)
    }

    fn canonical_stage(
        expression: Atom,
        structure: PartialStructure,
        proofs: &TensorProofs,
    ) -> Self {
        let mut stage = Self::from_validated_parts(
            expression,
            structure,
            Arc::clone(
                proofs
                    .occurrences
                    .get()
                    .expect("canonical source carries index history"),
            ),
        );
        stage.proofs.inherit_interfaces([proofs]);
        stage
    }

    pub(crate) fn canonize_parsed<Aind: AbsInd + ParseableAind + DummyAind>(
        mut net: SymbolicNet<Aind>,
        mut new_dummy: &mut dyn FnMut(usize) -> Aind,
        mut admission: Option<(PartialStructure, TensorProofs)>,
    ) -> Result<(Atom, Option<TensorProofs>), CanonicalizationError> {
        net.remove_antisymmetric_zero_terms();

        // Keep duality wrappers in the expression graph. The canonicalizer
        // follows the same base index inside dind(...), while the wrapper
        // preserves which occurrence is covariant. Erasing that distinction
        // can collapse opposite color metrics into an invalid tensor square.
        let mut dummies = BTreeMap::new();
        for (p, _, d) in net.graph.graph.iter_edges() {
            if p.is_paired()
                && let NetworkEdge::Slot(s) = d.data
            {
                let slot = if s.rep_name().is_dual() { s.dual() } else { *s };
                dummies.insert(slot.to_atom(), slot.rep());
            }
        }

        let external = net.graph.dangling_indices();
        // Reconstruct the admitted operation graph, without running tensor algebra
        // merely to hand its syntax to Symbolica's canonicalizer.
        let tree = net.graph.expr_tree().cast::<ChildVecStore<()>>();
        let root = net.graph.graph.node_id(net.graph.head());
        let expr = net
            .graph
            .to_expression_at(&tree, root, &mut |_, node| match node {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
                    Ok(Some(net.store.tensors[*index].expression.clone()))
                }
                NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => {
                    Ok(Some(net.store.get_scalar_ref(*index).clone()))
                }
                NetworkNode::Op(_) => Ok(None),
                _ => Err(TensorNetworkError::Other(eyre::eyre!(
                    "canonicalization requires stored symbolic leaves"
                ))),
            })
            .map_err(|error| CanonicalizationError::Prepare {
                reason: error.to_string(),
            })?;
        if let Some((interface, proofs)) = &mut admission {
            proofs.publish_interfaces(
                expr.as_view(),
                interface.clone(),
                net.store.tensors.iter().map(|tensor| &tensor.proofs),
            );
        }
        // Symbolica rejects tensor-containing negative powers even when their
        // base is a closed scalar network. Canonicalize those independent
        // contractions recursively, then keep them opaque during outer labeling.
        let mut aliases = BTreeMap::new();
        let mut used = expr.get_all_symbols(true);
        let mut serial = 0;
        let mut failure = None;
        let mut alias_proofs = Vec::new();
        let mut select_alias = |atom: AtomView<'_>| {
            let AtomView::Pow(power) = atom else {
                return None;
            };
            let (base, exponent) = power.get_base_exp();
            if i64::try_from(exponent).is_ok_and(|n| n > 0) {
                return None;
            }
            let result = if let Some((_, proofs)) = &admission {
                proofs
                    .admitted_interface(base)
                    .and_then(|interface| proofs.region_with_interface(base, interface))
                    .ok_or_else(|| CanonicalizationError::Prepare {
                        reason: "canonical scalar scope lacks admitted metadata".into(),
                    })
                    .and_then(|source| {
                        source.canonize_with(&mut new_dummy).map(|canonical| {
                            alias_proofs.push(canonical.proofs);
                            canonical.expression
                        })
                    })
            } else {
                base.canonize::<Aind>(&mut new_dummy as &mut dyn FnMut(usize) -> Aind)
            };
            let canonical = match result {
                Ok(canonical) => canonical.pow(exponent),
                Err(error) => {
                    failure = Some(error);
                    return None;
                }
            };
            let alias = aliases.entry(canonical).or_insert_with(|| {
                loop {
                    let name = symbol!(&format!("idenso::canonical_scalar_{serial}"));
                    serial += 1;
                    if used.insert(name) {
                        break Atom::var(name);
                    }
                }
            });
            Some(alias.clone())
        };
        let (expr, mapped) = if let Some((interface, proofs)) = &admission {
            let stage = Self::canonical_stage(expr, interface.clone(), proofs);
            let (expression, mapped) = stage.map_admitted_expression(
                |atom| Ok(select_alias(atom)),
                |interface| Ok(interface.clone()),
            )?;
            (expression, Some(mapped))
        } else {
            (
                expr.replace_map(|atom, _, out| {
                    if let Some(alias) = select_alias(atom) {
                        **out = alias;
                    }
                }),
                None,
            )
        };
        if let Some(error) = failure {
            return Err(error);
        }
        if let Some((interface, proofs)) = &mut admission {
            proofs.publish_interfaces(expr.as_view(), interface.clone(), mapped.iter());
        }
        let mut templates = admission
            .as_ref()
            .map(|(_, proofs)| CanonicalLeafTemplates::new(proofs, dummies.keys()))
            .transpose()?;
        // Spenso supplies validated incidence and external slots. Symbolica
        // canonizes the reconstructed expression, including dummy pairs inside
        // nested sums, without distributing any products of sums.
        let canonical =
            expr.canonize_tensors(dummies)
                .map_err(|error| CanonicalizationError::Indices {
                    reason: error.to_string(),
                })?;
        if let Some((interface, proofs)) = &mut admission {
            proofs.publish_interfaces(canonical.canonical_form.as_view(), interface.clone(), []);
            templates
                .as_ref()
                .expect("admitted canonicalization has leaf templates")
                .retain_output(&canonical.canonical_form, interface, proofs)?;
        }

        // Sum branches also have paired copies of genuinely external slots.
        // Preserve those labels after Symbolica has validated their incidence.
        let mut reserved = BTreeSet::new();
        for slot in external {
            let slot = if slot.rep_name().is_dual() {
                slot.dual()
            } else {
                slot
            };
            let slot = slot.to_atom();
            reserved.insert(slot);
        }

        let mut reps = vec![];
        let mut typed_replacements = Vec::new();
        let mut next_dummy = 0;

        // Only relabel names that Symbolica validated as contracted indices;
        // keep their replacements distinct from the genuine external slots.
        // Canceled terms can leave unused names in the validated index pool.
        for (d, r) in canonical
            .dummy_indices
            .into_iter()
            .filter(|(index, _)| canonical.canonical_form.contains(index))
        {
            let target = loop {
                let candidate = r.slot::<Aind, Aind>(new_dummy(next_dummy)).to_atom();
                next_dummy += 1;
                if reserved.insert(candidate.clone()) {
                    break candidate;
                }
            };
            if admission.is_some() {
                let from = Slot::<LibraryRep, AbstractIndex>::try_from(d.as_view())
                    .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
                let to = Slot::<LibraryRep, AbstractIndex>::try_from(target.as_view())
                    .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
                templates
                    .as_mut()
                    .expect("admitted canonicalization has leaf templates")
                    .add_slot(to);
                typed_replacements.push((from, to));
                typed_replacements.push((from.dual(), to.dual()));
            }
            reps.push(Replacement::new(
                Pattern::Literal(d),
                Pattern::Literal(target),
            ));
        }

        if let Some((interface, mut proofs)) = admission {
            proofs.publish_interfaces(canonical.canonical_form.as_view(), interface.clone(), []);
            let atoms = typed_replacements
                .iter()
                .map(|(from, to)| (from.to_atom(), to.to_atom()))
                .collect::<Vec<_>>();
            let stage = Self::canonical_stage(canonical.canonical_form, interface.clone(), &proofs);
            let (expression, mapped) = stage.map_admitted_expression(
                |atom| {
                    Ok(atoms
                        .iter()
                        .find(|(from, _)| from.as_view() == atom)
                        .map(|(_, to)| to.clone()))
                },
                |structure| {
                    Ok(PartialStructure::from_logical_slots(
                        structure.logical_slots().into_iter().map(|slot| {
                            let PartialIndex::Explicit(index) = slot.aind else {
                                return slot;
                            };
                            typed_replacements
                                .iter()
                                .find(|(from, _)| *from == slot.rep().slot(index))
                                .map_or(slot, |(_, to)| {
                                    to.rep().slot(PartialIndex::Explicit(to.aind()))
                                })
                        }),
                    ))
                },
            )?;
            proofs.publish_interfaces(expression.as_view(), interface.clone(), [&mapped]);
            templates
                .as_ref()
                .expect("admitted canonicalization has leaf templates")
                .retain_output(&expression, &interface, &mut proofs)?;
            let stage = Self::canonical_stage(expression, interface.clone(), &proofs);
            let (expression, mapped) = stage.map_admitted_expression(
                |atom| {
                    Ok(aliases
                        .iter()
                        .find(|(_, alias)| alias.as_view() == atom)
                        .map(|(original, _)| original.clone()))
                },
                |structure| Ok(structure.clone()),
            )?;
            proofs.publish_interfaces(
                expression.as_view(),
                interface.clone(),
                alias_proofs.iter().chain([&mapped]),
            );
            templates
                .as_ref()
                .expect("admitted canonicalization has leaf templates")
                .retain_output(&expression, &interface, &mut proofs)?;
            Ok((expression, Some(proofs)))
        } else {
            let mut scalar_definitions = AliasedAtom::from(Atom::Zero);
            for (original, alias) in aliases {
                scalar_definitions.register_alias(alias, original);
            }
            Ok((
                scalar_definitions
                    .map_root(|_| canonical.canonical_form.replace_multiple(&reps))
                    .into_inner(),
                None,
            ))
        }
    }
}

/// Source leaves may only change the dummy slots supplied to the canonicalizer.
/// Keep the immutable source cache alive while publishing a separate output cache;
/// matching never rebuilds a function or discovers a port from output syntax.
struct CanonicalLeafTemplates {
    source: LeafInterfaceCache,
    pool: AHashMap<Atom, Representation<LibraryRep>>,
    user_normalizers: bool,
}

impl CanonicalLeafTemplates {
    fn new<'a>(
        proofs: &TensorProofs,
        dummies: impl IntoIterator<Item = &'a Atom>,
    ) -> Result<Self, TensorInferenceError> {
        let source = proofs.leaf_interfaces.read_cache().unwrap_or_default();
        let user_normalizers = source
            .lock()
            .expect("leaf interface cache is not poisoned")
            .keys()
            .any(|atom| {
                atom.get_all_symbols(true)
                    .into_iter()
                    .any(|head| !InterfaceInference::intrinsic_normalization_head(head))
            });
        let mut result = Self {
            source,
            pool: AHashMap::new(),
            user_normalizers,
        };
        for dummy in dummies {
            let slot = Slot::<LibraryRep, AbstractIndex>::try_from(dummy.as_view())
                .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
            result.add_slot(slot);
        }
        Ok(result)
    }

    fn add_slot(&mut self, slot: Slot<LibraryRep, AbstractIndex>) {
        for slot in [slot, slot.dual()] {
            self.pool.insert(slot.to_atom(), slot.rep());
        }
    }

    fn children(value: AtomView<'_>) -> Vec<AtomView<'_>> {
        match value {
            AtomView::Fun(function) => function.iter().collect(),
            AtomView::Add(sum) => sum.iter().collect(),
            AtomView::Mul(product) => product.iter().collect(),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                vec![base, exponent]
            }
            _ => Vec::new(),
        }
    }

    fn fingerprint(&self, value: AtomView<'_>) -> u64 {
        fn hash_node<'a>(
            templates: &CanonicalLeafTemplates,
            value: AtomView<'a>,
            labels: &mut AHashMap<&'a [u8], usize>,
            hash: &mut AHasher,
        ) {
            if let Some(rep) = templates.pool.get(value.get_data()) {
                0_u8.hash(hash);
                rep.hash(hash);
                let next = labels.len();
                labels.entry(value.get_data()).or_insert(next).hash(hash);
                return;
            }
            match value {
                AtomView::Fun(function) => {
                    1_u8.hash(hash);
                    function.get_symbol().hash(hash);
                }
                AtomView::Add(_) => 2_u8.hash(hash),
                AtomView::Mul(_) => 3_u8.hash(hash),
                AtomView::Pow(_) => 4_u8.hash(hash),
                _ => {
                    5_u8.hash(hash);
                    value.get_data().hash(hash);
                }
            }
            let children = CanonicalLeafTemplates::children(value);
            children.len().hash(hash);
            for child in children {
                hash_node(templates, child, labels, hash);
            }
        }
        let mut hash = AHasher::default();
        hash_node(self, value, &mut AHashMap::new(), &mut hash);
        hash.finish()
    }

    fn bind<'a, 'b>(
        &self,
        source: AtomView<'a>,
        target: AtomView<'b>,
        mapping: &mut AHashMap<&'a [u8], AtomView<'b>>,
        used: &mut AHashSet<&'b [u8]>,
    ) -> bool {
        if let Some(rep) = self.pool.get(source.get_data()) {
            if self.pool.get(target.get_data()) != Some(rep) {
                return false;
            }
            if let Some(previous) = mapping.get(source.get_data()) {
                return *previous == target;
            }
            if !used.insert(target.get_data()) {
                return false;
            }
            mapping.insert(source.get_data(), target);
            return true;
        }
        let same_kind = match (source, target) {
            (AtomView::Fun(left), AtomView::Fun(right)) => left.get_symbol() == right.get_symbol(),
            (AtomView::Add(_), AtomView::Add(_))
            | (AtomView::Mul(_), AtomView::Mul(_))
            | (AtomView::Pow(_), AtomView::Pow(_)) => true,
            _ => return source == target,
        };
        let left = Self::children(source);
        let right = Self::children(target);
        same_kind
            && left.len() == right.len()
            && left
                .into_iter()
                .zip(right)
                .all(|(left, right)| self.bind(left, right, mapping, used))
    }

    fn retain_output(
        &self,
        expression: &Atom,
        interface: &PartialStructure,
        proofs: &mut TensorProofs,
    ) -> Result<(), TensorInferenceError> {
        if self.user_normalizers {
            // The kernel can invoke user code while rebuilding an indexed leaf.
            // Observe that explicit user transformation once, without replaying it.
            let observed = SymbolicTensor::from_observed_parts(
                expression.clone(),
                interface.clone(),
                LeafInference::ObserveEncoded,
            );
            #[cfg(debug_assertions)]
            assert!(
                observed.is_ok(),
                "normalization callback contract violated in canonize: {observed:?}"
            );
            proofs.inherit_interfaces([&observed?.proofs]);
            return Ok(());
        }
        self.retain_matching(expression, proofs, true)
    }

    fn retain_matching(
        &self,
        expression: &Atom,
        proofs: &TensorProofs,
        require_all: bool,
    ) -> Result<(), TensorInferenceError> {
        let source = self
            .source
            .lock()
            .expect("leaf interface cache is not poisoned");
        let mut candidates = AHashMap::<_, Vec<_>>::new();
        for (atom, (tensorial, ports)) in source.iter() {
            if *tensorial && matches!(atom.as_view(), AtomView::Fun(_)) {
                candidates
                    .entry(self.fingerprint(atom.as_view()))
                    .or_default()
                    .push((atom, ports));
            }
        }
        let target = proofs.parser_leaf_interfaces();
        let mut pending = vec![expression.as_view()];
        while let Some(value) = pending.pop() {
            let AtomView::Fun(function) = value else {
                pending.extend(Self::children(value));
                continue;
            };
            if target
                .lock()
                .expect("leaf interface cache is not poisoned")
                .get(value.get_data())
                .is_some_and(|(tensorial, _)| *tensorial)
            {
                continue;
            }
            let mut matched = false;
            for (template, ports) in candidates
                .get(&self.fingerprint(value))
                .into_iter()
                .flatten()
            {
                let mut mapping = AHashMap::new();
                if self.bind(
                    template.as_view(),
                    value,
                    &mut mapping,
                    &mut AHashSet::new(),
                ) {
                    proofs.retain_leaf_ports(
                        value,
                        ports
                            .iter()
                            .map(|port| {
                                mapping
                                    .get(port.as_view().get_data())
                                    .map_or_else(|| port.clone(), |value| (*value).to_owned())
                            })
                            .collect(),
                    )?;
                    matched = true;
                    break;
                }
            }
            if matched {
                continue;
            }
            if require_all && composition::is_tensor_leaf_head(function.get_symbol()) {
                return Err(TensorInferenceError::invalid(
                    "canonicalizer output has no retained source leaf template",
                ));
            }
            pending.extend(function.iter());
        }
        Ok(())
    }
}

impl TensorProofs {
    /// Carry constructor roles through the dummy permutations of a local
    /// identity. The identity supplies the slots it renamed and separately
    /// owns any builtin constructors it creates.
    pub(crate) fn retain_permuted_leaf_ports<'a>(
        &self,
        expression: &Atom,
        slots: impl IntoIterator<Item = &'a Atom>,
        output: &TensorProofs,
    ) -> Result<(), TensorInferenceError> {
        CanonicalLeafTemplates::new(self, slots)?.retain_matching(expression, output, false)
    }
}

const CYCLE_SYMMETRY_MARKER: usize = usize::MAX - 2;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TensorColor {
    head: Symbol,
    arguments: Vec<Atom>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TensorGraphNode<Aind> {
    Product,
    Tensor(TensorColor),
    /// A nonlinear function, sum, or unresolved leaf, identified per occurrence.
    Opaque(NodeIndex),
    Slot(SlotColor<Aind>),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SlotColor<Aind> {
    Internal(Representation<LibraryRep>),
    External(LibrarySlot<Aind>),
}

type TensorGraph<Aind> =
    Graph<TensorGraphNode<Aind>, HiddenData<(usize, Option<Representation<LibraryRep>>), usize>>;

struct SlotCopies<Aind> {
    slot: LibrarySlot<Aind>,
    power_path: Vec<NodeIndex>,
    vertices: Vec<usize>,
}

type SlotsByHedge<Aind> = BTreeMap<Hedge, SlotCopies<Aind>>;

#[cfg(test)]
thread_local! {
    static PROJECTED_SLOT_VISITS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Prune signed zeros on the already parsed network, without distributing sums.
///
/// Each mutation removes an operation/subtree, or replaces a nonzero tensor by
/// zero. Rebuild only the graph traversal after compacting node/hedge storage;
/// tensor syntax and interfaces are never parsed again. Nonlinear functions
/// and nonpositive powers retain their opaque boundary.
pub(super) fn remove_antisymmetric_zero_terms<Aind: AbsInd + ParseableAind>(
    network: &mut SymbolicNet<Aind>,
) -> bool {
    if !contains_antisymmetric_tensor(network) {
        return false;
    }
    let mut changed = false;
    loop {
        let tree: SimpleTraversalTree<ParentChildStore<()>> = network.graph.expr_tree().cast();
        let root = network.graph.graph.node_id(network.graph.head());
        if !prune_one(network, &tree, root) {
            return changed;
        }
        changed = true;
    }
}

fn is_zero<Aind: AbsInd>(network: &SymbolicNet<Aind>, node: NodeIndex) -> bool {
    match &network.graph.graph[node] {
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => {
            network.store.tensors[*index].expression.is_zero()
        }
        NetworkNode::Leaf(NetworkLeaf::Scalar(index)) => network.store.scalar_ref(*index).is_zero(),
        _ => false,
    }
}

fn prune_one<Aind: AbsInd + ParseableAind>(
    network: &mut SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    node: NodeIndex,
) -> bool {
    if matches!(
        network.graph.graph[node],
        NetworkNode::Op(NetworkOp::Function(_) | NetworkOp::Power(..=0))
    ) || is_zero(network, node)
    {
        return false;
    }
    let children = tree
        .iter_children(node, network.graph.graph.as_ref())
        .collect::<Vec<_>>();
    for child in &children {
        if prune_one(network, tree, *child) {
            return true;
        }
    }
    let zeros = children
        .iter()
        .filter(|child| is_zero(network, **child))
        .count();
    match network.graph.graph[node] {
        NetworkNode::Op(NetworkOp::Sum) if zeros != 0 => {
            if zeros == children.len() {
                replace_with_zero(network, tree, node);
            } else {
                // A zero child has already become a leaf. Remove both halves
                // of its head and Sum-input slots, leaving the Sum output and
                // every surviving arm's incidence unchanged.
                let mut removed: SuBitGraph = network.graph.graph.empty_subgraph();
                for child in children.iter().filter(|child| is_zero(network, **child)) {
                    for hedge in network.graph.graph.iter_crown(*child) {
                        removed.add(hedge);
                        removed.add(network.graph.graph.inv(hedge));
                    }
                }
                network.graph.delete(&removed);
            }
            return true;
        }
        NetworkNode::Op(NetworkOp::Sum) if children.len() == 1 => {
            return remove_singleton_sum(network, tree, node, children[0]);
        }
        NetworkNode::Op(NetworkOp::Product | NetworkOp::Neg | NetworkOp::Power(1..))
            if zeros != 0 =>
        {
            replace_with_zero(network, tree, node);
            return true;
        }
        _ => {}
    }
    if has_odd_automorphism(network, tree, node) {
        replace_with_zero(network, tree, node);
        return true;
    }
    false
}

fn replace_with_zero<Aind: AbsInd>(
    network: &mut SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    root: NodeIndex,
) {
    let nodes = tree
        .iter_preorder_tree_nodes(network.graph.graph.as_ref(), root)
        .collect::<Vec<_>>();
    let inside = nodes
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let mut boundary = Vec::new();
    for node in &nodes {
        for hedge in network.graph.graph.iter_crown(*node) {
            let NetworkEdge::Slot(edge_slot) = network.graph.graph[[&hedge]] else {
                continue;
            };
            let other = network.graph.graph.inv(hedge);
            if other != hedge && inside.contains(&network.graph.graph.node_id(other)) {
                continue;
            }
            // Sewing a dual pair can retain the opposite endpoint's descriptor.
            // Recover this endpoint from its already established storage axis.
            let slot = match &network.graph.graph[*node] {
                NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => network.store.tensors[*index]
                    .structure
                    .external_structure_iter()
                    .nth(usize::from(network.graph.slot_order[hedge.0]))
                    .expect("a parsed tensor endpoint has a storage axis"),
                NetworkNode::Op(NetworkOp::Sum) => {
                    // A Sum output forwards each arm's identical interface.
                    // Its input seam still carries that descriptor even when
                    // an outer product retained the dual output descriptor.
                    network
                        .graph
                        .graph
                        .iter_crown(*node)
                        .find_map(|input| {
                            let NetworkEdge::Slot(slot) = network.graph.graph[[&input]] else {
                                return None;
                            };
                            let other = network.graph.graph.inv(input);
                            (other != input
                                && inside.contains(&network.graph.graph.node_id(other))
                                && (slot == edge_slot || slot.matches(&edge_slot)))
                            .then_some(slot)
                        })
                        .expect("a parsed Sum output has a matching input port")
                }
                _ => edge_slot,
            };
            boundary.push((hedge, slot));
        }
    }
    let layout = OrderedStructure::new(boundary.iter().map(|(_, slot)| *slot).collect());
    let slots = layout.canonical().external_structure();
    for (position, slot) in slots.iter().enumerate() {
        let found = boundary
            .iter()
            .position(|(_, current)| current == slot)
            .expect("the canonical zero interface permutes its boundary ports");
        let (hedge, _) = boundary.remove(found);
        network.graph.slot_order[hedge.0] =
            u8::try_from(position).expect("parsed slot position fits u8");
    }
    let index = network.store.tensors.len();
    network.graph.identify_nodes_without_self_edges(
        &nodes,
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)),
    );
    // Identification intentionally preserves pre-existing self loops. A zero
    // replacement has already consumed those internal contractions, so remove
    // their paired ports too. Recover its occurrence after node compaction.
    let node = network
        .graph
        .graph
        .iter_nodes()
        .find_map(|(node, _, data)| {
            matches!(data, NetworkNode::Leaf(NetworkLeaf::LocalTensor(found)) if *found == index)
                .then_some(node)
        })
        .expect("the replacement occurrence survives graph compaction");
    let mut internal: SuBitGraph = network.graph.graph.empty_subgraph();
    for hedge in network.graph.graph.iter_crown(node) {
        let other = network.graph.graph.inv(hedge);
        if other != hedge && network.graph.graph.node_id(other) == node {
            internal.add(hedge);
        }
    }
    network.graph.delete(&internal);
    // Keep the proven boundary until the existing execution owner consumes it;
    // deleting an open zero's ports here would corrupt the enclosing Sum.
    network.store.tensors.push(SymbolicTensor {
        work: Default::default(),
        proofs: Default::default(),
        expression: Atom::zero(),
        structure: layout.into_canonical(),
        is_metric: false,
        is_composite: true,
    });
}

fn remove_singleton_sum<Aind: AbsInd>(
    network: &mut SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    node: NodeIndex,
    child: NodeIndex,
) -> bool {
    let descendants = tree
        .iter_preorder_tree_nodes(network.graph.graph.as_ref(), child)
        .collect::<std::collections::BTreeSet<_>>();
    let crown = network.graph.graph.iter_crown(node).collect::<Vec<_>>();
    let mut seams = Vec::new();
    for input in &crown {
        let inside = network.graph.graph.inv(*input);
        if inside == *input || !descendants.contains(&network.graph.graph.node_id(inside)) {
            continue;
        }
        let input_data = network.graph.graph[[input]];
        let Some(output) = crown.iter().copied().find(|output| {
            let outside = network.graph.graph.inv(*output);
            let same_port = match (network.graph.graph[[output]], input_data) {
                (NetworkEdge::Head, NetworkEdge::Head) => true,
                (NetworkEdge::Slot(output), NetworkEdge::Slot(input)) => {
                    output == input || output.matches(&input)
                }
                _ => false,
            };
            same_port
                && (outside == *output
                    || !descendants.contains(&network.graph.graph.node_id(outside)))
        }) else {
            return false;
        };
        seams.push((output, *input));
    }
    if seams.len() * 2 != crown.len() {
        return false;
    }
    // Bypass the existing Sum seams, retaining the outside edge's descriptor
    // and orientation. No tensor storage axis or source logical witness moves.
    for (output, input) in seams {
        let outside = network.graph.graph.inv(output);
        let inside = network.graph.graph.inv(input);
        let input_data = network
            .graph
            .graph
            .get_edge_data_full(input)
            .map(|data| *data);
        network
            .graph
            .graph
            .split_edge(input, input_data)
            .expect("Sum input is paired");
        if outside != output {
            let output_data = network
                .graph
                .graph
                .get_edge_data_full(output)
                .map(|data| *data);
            network
                .graph
                .graph
                .split_edge(output, output_data)
                .expect("Sum output is paired");
            network
                .graph
                .graph
                .connect_identities(outside, inside, |flow, data, _, _| (flow, data));
        }
    }
    let mut removed: SuBitGraph = network.graph.graph.empty_subgraph();
    for hedge in crown {
        removed.add(hedge);
    }
    network.graph.delete(&removed);
    true
}

fn contains_antisymmetric_tensor<Aind: AbsInd>(network: &SymbolicNet<Aind>) -> bool {
    network.graph.graph.iter_nodes().any(|(_, _, node)| {
        let NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) = node else {
            return false;
        };
        network.store.tensors[*index]
            .name()
            .is_some_and(|head| head.is_antisymmetric())
    })
}

fn has_odd_automorphism<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    root: NodeIndex,
) -> bool {
    let Some(graph) = project_network(network, tree, root) else {
        return false;
    };
    // An odd stabilizer generator proves that the monomial equals its negative.
    let canonical = graph.canonize();
    canonical
        .orbit_generators
        .iter()
        .any(|generator| generator_is_odd(&canonical.graph, generator))
}

fn project_network<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    root: NodeIndex,
) -> Option<TensorGraph<Aind>> {
    let mut graph = Graph::new();
    let mut slot_copies = BTreeMap::new();
    project_expression(network, tree, root, &mut graph, &mut slot_copies)?;
    // An antisymmetric tensor elsewhere in the network cannot make this
    // projected region vanish. Avoid wiring an unrelated region's slots.
    if !graph.nodes().iter().any(|node| {
        matches!(&node.data,
        TensorGraphNode::Tensor(color) if color.head.is_antisymmetric())
    }) {
        return None;
    }
    connect_slots(network, &mut graph, &slot_copies)?;
    Some(graph)
}

fn project_expression<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    node: NodeIndex,
    graph: &mut TensorGraph<Aind>,
    slot_copies: &mut SlotsByHedge<Aind>,
) -> Option<usize> {
    match &network.graph.graph[node] {
        NetworkNode::Op(NetworkOp::Product) => {
            let header = graph.add_node(TensorGraphNode::Product);
            for child in tree.iter_children(node, network.graph.graph.as_ref()) {
                let child = project_expression(network, tree, child, graph, slot_copies)?;
                add_incidence(graph, header, child);
            }
            Some(header)
        }
        NetworkNode::Op(NetworkOp::Power(power)) if *power > 0 => {
            let child = tree
                .iter_children(node, network.graph.graph.as_ref())
                .next()?;
            let header = graph.add_node(TensorGraphNode::Product);
            for _ in 0..*power {
                let copy = project_expression(network, tree, child, graph, slot_copies)?;
                add_incidence(graph, header, copy);
            }
            Some(header)
        }
        NetworkNode::Op(NetworkOp::Function(_)) => {
            Some(graph.add_node(TensorGraphNode::Opaque(node)))
        }
        NetworkNode::Op(NetworkOp::Neg) => {
            let child = tree
                .iter_children(node, network.graph.graph.as_ref())
                .next()?;
            project_expression(network, tree, child, graph, slot_copies)
        }
        // Sums stay factorized, and function arguments stay opaque. An
        // automorphism confined to one summand does not multiply its parent,
        // but independent factors remain analyzable.
        NetworkNode::Op(NetworkOp::Sum) => Some(graph.add_node(TensorGraphNode::Opaque(node))),
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => add_tensor(
            network,
            node,
            &network.store.tensors[*index],
            tree,
            graph,
            slot_copies,
        ),
        NetworkNode::Leaf(_) => Some(graph.add_node(TensorGraphNode::Opaque(node))),
        NetworkNode::Op(NetworkOp::Power(_)) => None,
    }
}

fn add_tensor<Aind: AbsInd + ParseableAind>(
    network: &SymbolicNet<Aind>,
    network_node: NodeIndex,
    tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    tree: &SimpleTraversalTree<ParentChildStore<()>>,
    graph: &mut TensorGraph<Aind>,
    slot_copies: &mut SlotsByHedge<Aind>,
) -> Option<usize> {
    let slots = tensor
        .structure
        .external_structure_iter()
        .collect::<Vec<_>>();
    let color = tensor_color(tensor, &slots)?;
    let head = color.head;
    let unordered = head.is_symmetric() || head.is_antisymmetric() || head.is_cyclesymmetric();
    let header = graph.add_node(TensorGraphNode::Tensor(color));
    let mut network_slots = network
        .graph
        .graph
        .iter_crown(network_node)
        .filter(|hedge| matches!(network.graph.graph[[hedge]], NetworkEdge::Slot(_)))
        .collect::<Vec<_>>();
    network_slots.sort_unstable_by_key(|hedge| network.graph.slot_order[hedge.0]);
    if network_slots.len() != slots.len() {
        return None;
    }
    let power_path = tree
        .ancestor_iter_node(network_node, network.graph.graph.as_ref())
        .skip(1)
        .filter(|ancestor| {
            matches!(
                network.graph.graph[*ancestor],
                NetworkNode::Op(NetworkOp::Power(1..))
            )
        })
        .collect::<Vec<_>>();

    let mut cycle_slots = Vec::with_capacity(slots.len());
    for (position, (slot, network_slot)) in slots.iter().copied().zip(network_slots).enumerate() {
        let vertex = graph.add_node(TensorGraphNode::Slot(SlotColor::External(slot)));
        let visible_position = if unordered { 0 } else { position };
        graph
            .add_edge(
                header,
                vertex,
                true,
                HiddenData::new((visible_position, None), position),
            )
            .unwrap();
        slot_copies
            .entry(network_slot)
            .or_insert_with(|| SlotCopies {
                slot,
                power_path: power_path.clone(),
                vertices: vec![],
            })
            .vertices
            .push(vertex);
        cycle_slots.push(vertex);
    }

    if head.is_cyclesymmetric() && cycle_slots.len() > 1 {
        for (&left, &right) in cycle_slots
            .iter()
            .zip(cycle_slots.iter().cycle().skip(1))
            .take(cycle_slots.len())
        {
            graph
                .add_edge(
                    left,
                    right,
                    true,
                    HiddenData::new((CYCLE_SYMMETRY_MARKER, None), 0),
                )
                .unwrap();
        }
    }
    Some(header)
}

fn tensor_color<Aind: AbsInd + ParseableAind>(
    tensor: &SymbolicTensor<OrderedStructure<LibraryRep, Aind>>,
    slots: &[LibrarySlot<Aind>],
) -> Option<TensorColor> {
    let AtomView::Fun(function) = tensor.expression.as_view() else {
        return None;
    };
    let mut slot_atoms = slots
        .iter()
        .map(IsAbstractSlot::to_atom)
        .collect::<Vec<_>>();
    let mut arguments = Vec::new();
    for argument in function.iter() {
        if let Some(position) = slot_atoms
            .iter()
            .position(|slot| slot.as_view() == argument)
        {
            let _ = slot_atoms.swap_remove(position);
        } else {
            arguments.push(argument.to_owned());
        }
    }
    if !slot_atoms.is_empty() {
        return None;
    }

    Some(TensorColor {
        head: function.get_symbol(),
        arguments,
    })
}

fn add_incidence<Aind: AbsInd>(graph: &mut TensorGraph<Aind>, parent: usize, child: usize) {
    graph
        .add_edge(parent, child, true, HiddenData::new((0, None), 0))
        .unwrap();
}

fn connect_slots<Aind: AbsInd>(
    network: &SymbolicNet<Aind>,
    graph: &mut TensorGraph<Aind>,
    slot_copies: &SlotsByHedge<Aind>,
) -> Option<()> {
    // Projected tensor ports already identify every relevant network edge.
    // Visiting the whole network for each small region repeats unrelated sums
    // and opaque leaves. Keep the stored source/sink orientation, including an
    // edge whose only projected endpoint is its sink.
    for &hedge in slot_copies.keys() {
        #[cfg(test)]
        PROJECTED_SLOT_VISITS.with(|visits| visits.set(visits.get() + 1));
        let HedgePair::Paired { source, sink } = network.graph.graph.hedge_pair(hedge) else {
            continue;
        };
        if hedge == sink && slot_copies.contains_key(&source) {
            continue;
        }
        let NetworkEdge::Slot(slot) = network.graph.graph[[&source]] else {
            continue;
        };
        let source_copies = slot_copies.get(&source);
        let sink_copies = slot_copies.get(&sink);
        let group = if slot.rep_name().is_dual() {
            slot.rep().dual()
        } else {
            slot.rep()
        };

        match (source_copies, sink_copies) {
            (Some(source), Some(sink))
                if source.vertices.len() == sink.vertices.len()
                    && (source.vertices.len() == 1 || source.power_path == sink.power_path) =>
            {
                for (&source_vertex, &sink_vertex) in source.vertices.iter().zip(&sink.vertices) {
                    connect_slot_pair(
                        graph,
                        source_vertex,
                        sink_vertex,
                        source.slot.rep(),
                        sink.slot.rep(),
                        &group,
                    );
                }
            }
            (Some(copies), None) if is_power_node(network, sink) => {
                let [left, right] = copies.vertices.as_slice() else {
                    return None;
                };
                connect_slot_pair(
                    graph,
                    *left,
                    *right,
                    copies.slot.rep(),
                    copies.slot.rep(),
                    &group,
                );
            }
            (None, Some(copies)) if is_power_node(network, source) => {
                let [left, right] = copies.vertices.as_slice() else {
                    return None;
                };
                connect_slot_pair(
                    graph,
                    *left,
                    *right,
                    copies.slot.rep(),
                    copies.slot.rep(),
                    &group,
                );
            }
            // An opaque subtree can own the other endpoint. Keeping this
            // endpoint fixed is conservative and still permits independent factors.
            (Some(_), None) | (None, Some(_)) => {}
            (None, None) => {}
            _ => return None,
        }
    }
    if slot_copies.iter().any(|(hedge, copies)| {
        copies.vertices.len() != 1 && network.graph.graph.as_ref().is_identity(*hedge)
    }) {
        return None;
    }
    Some(())
}

fn is_power_node<Aind: AbsInd>(network: &SymbolicNet<Aind>, hedge: Hedge) -> bool {
    matches!(
        network.graph.graph[network.graph.graph.node_id(hedge)],
        NetworkNode::Op(NetworkOp::Power(_))
    )
}

fn connect_slot_pair<Aind: AbsInd>(
    graph: &mut TensorGraph<Aind>,
    source: usize,
    sink: usize,
    source_rep: Representation<LibraryRep>,
    sink_rep: Representation<LibraryRep>,
    group: &Representation<LibraryRep>,
) {
    graph.set_node_data(
        source,
        TensorGraphNode::Slot(SlotColor::Internal(source_rep)),
    );
    graph.set_node_data(sink, TensorGraphNode::Slot(SlotColor::Internal(sink_rep)));
    graph
        .add_edge(source, sink, false, HiddenData::new((0, Some(*group)), 0))
        .unwrap();
}

fn generator_is_odd<Aind: AbsInd>(graph: &TensorGraph<Aind>, cycles: &[Vec<usize>]) -> bool {
    let mut automorphism = (0..graph.nodes().len()).collect::<Vec<_>>();
    for cycle in cycles {
        for (source, target) in cycle
            .iter()
            .zip(cycle.iter().cycle().skip(1))
            .take(cycle.len())
        {
            automorphism[*source] = *target;
        }
    }

    let mut odd = false;
    for node_index in 0..graph.nodes().len() {
        let TensorGraphNode::Tensor(color) = &graph.node(node_index).data else {
            continue;
        };
        if !color.head.is_antisymmetric() {
            continue;
        }

        let target_header = automorphism[node_index];
        let mut induced_permutation = vec![];
        for &edge_index in &graph.node(node_index).edges {
            let edge = graph.edge(edge_index);
            if !edge.directed || edge.vertices.0 != node_index {
                continue;
            }

            let target_slot = automorphism[edge.vertices.1];
            let target_position = graph
                .node(target_header)
                .edges
                .iter()
                .map(|edge_index| graph.edge(*edge_index))
                .find(|target_edge| {
                    target_edge.directed && target_edge.vertices == (target_header, target_slot)
                })
                .expect("an automorphism preserves tensor-slot incidence")
                .data
                .hidden;
            induced_permutation.push((edge.data.hidden, target_position));
        }

        induced_permutation.sort_unstable_by_key(|(source, _)| *source);
        odd ^= Permutation::from_map(
            induced_permutation
                .into_iter()
                .map(|(_, target)| target)
                .collect(),
        )
        .sign()
            < 0;
    }
    odd
}

#[cfg(test)]
mod typed_tests {
    use super::*;
    use crate::tensor::SymbolicNetParse;
    use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS};
    use crate::test_support;
    use spenso::network::parsing::ParseSettings;
    use spenso::structure::partial::{PartialIndex, PartialStructureExt};
    use spenso::structure::representation::Minkowski;
    use spenso::{tensor_symbol, vector_symbol};
    use std::collections::HashMap;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::{AtomCore, FunctionBuilder};
    use symbolica::{function, symbol};

    #[test]
    fn signed_zero_projection_visits_only_incident_slots() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let [a, b, c] = [86401, 86402, 86403].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let cycle = spenso::antisym!(a.clone(), b.clone())
            * spenso::antisym!(b, c.clone())
            * spenso::antisym!(c, a);

        for size in [0, 128] {
            let spectator = if size == 0 {
                Atom::one()
            } else {
                Atom::add_many((0..size).map(|i| {
                    let slot = rep
                        .slot::<AbstractIndex, _>(AbstractIndex::Normal(86500 + i))
                        .to_atom();
                    spenso::tensor!(signed_projection_spectator_left, slot.clone())
                        * spenso::tensor!(signed_projection_spectator_right, slot)
                }))
            };
            let expression = &cycle * spectator;
            let network = expression
                .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap();
            if size != 0 {
                assert!(network.graph.graph.iter_edges().count() > 2 * size);
            }
            let tree = network.graph.expr_tree().cast();
            let root = network.graph.graph.node_id(network.graph.head());
            PROJECTED_SLOT_VISITS.with(|visits| visits.set(0));
            assert!(has_odd_automorphism(&network, &tree, root));
            // The closed three-cycle has six tensor ports. None of the
            // hundreds of edges inside the opaque spectator is revisited.
            assert_eq!(PROJECTED_SLOT_VISITS.with(|visits| visits.get()), 6);
        }
    }

    #[test]
    fn projection_without_antisymmetric_tensors_skips_slot_connections() {
        crate::test_support::test_initialize();
        let a = spenso::mink!(4, signed_projection_plain_a);
        let b = spenso::mink!(4, signed_projection_plain_b);
        let expression =
            spenso::sym!(a.clone(), b.clone()) * spenso::tensor!(signed_projection_plain, a, b);
        let network = expression
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        let tree = network.graph.expr_tree().cast();
        let root = network.graph.graph.node_id(network.graph.head());
        PROJECTED_SLOT_VISITS.with(|visits| visits.set(0));
        assert!(!has_odd_automorphism(&network, &tree, root));
        assert_eq!(PROJECTED_SLOT_VISITS.with(|visits| visits.get()), 0);
    }

    #[test]
    fn incident_slot_projection_keeps_opaque_sum_endpoints_external() {
        crate::test_support::test_initialize();
        let a = spenso::lor!(4, signed_projection_opaque_a);
        let b = spenso::lor!(4, signed_projection_opaque_b);
        let sum = spenso::sym!(spenso::dind!(a.clone()), spenso::dind!(b.clone()))
            + spenso::tensor!(
                signed_projection_opaque_ordered,
                spenso::dind!(a.clone()),
                spenso::dind!(b.clone())
            );
        let expression = spenso::antisym!(a, b) * &sum;
        let mut network = expression
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        // The symmetric arm would vanish after distribution, but the ordered
        // arm can have a nonzero antisymmetric part. An opaque sum must not
        // identify its parent's distinct boundary ports with internal slots.
        assert!(!remove_antisymmetric_zero_terms(&mut network));
        let canonical = expression
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        assert!(!canonical.is_zero());
        assert!(matches!(canonical.as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| matches!(factor, AtomView::Add(_)))));
        assert!(
            canonical
                .list_dangling::<AbstractIndex>()
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn typed_canonize_keeps_data_of_opaque_untagged_calls() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::from(spenso::structure::representation::Euclidean {}).new_rep(2);
        let slot = |index: AbstractIndex| rep.slot::<AbstractIndex, _>(index).to_atom();
        let [i, j] = [87_101, 87_102].map(|index| slot(AbstractIndex::Normal(index)));
        let [d0, d1] = [0, 1].map(|index| slot(AbstractIndex::Dummy(index)));
        let call =
            |head: Symbol, argument: &Atom| FunctionBuilder::new(head).add_arg(argument).finish();
        let [u0, u1, w] = ["u0", "u1", "w"].map(|name| symbol!(&format!("typed_opaque_{name}")));
        let [p, q, r] = ["p", "q", "r"].map(|name| {
            spenso::network::tags::SPENSO_TAG
                .rank_one_tensor_symbol(&format!("typed_opaque_{name}"))
        });
        let x = Atom::var(symbol!("typed_opaque_x"));
        // Untagged calls are opaque scalars: the names inside are data. They
        // stay verbatim, they never pair with a tagged index, and tagged pairs
        // are still relabelled.
        let cases = [
            (call(u0, &i) * call(u1, &i), call(u0, &i) * call(u1, &i)),
            (
                call(u0, &i).pow(Atom::num(2)),
                call(u0, &i).pow(Atom::num(2)),
            ),
            (
                (&x + call(u0, &call(p, &i)) * call(u1, &call(q, &i))) * call(p, &j) * call(q, &j),
                (&x + call(u0, &call(p, &i)) * call(u1, &call(q, &i)))
                    * call(p, &d0)
                    * call(q, &d0),
            ),
            (
                call(w, &call(p, &i)) * call(q, &j) * call(r, &j),
                call(w, &call(p, &i)) * call(q, &d0) * call(r, &d0),
            ),
            // The free index `i` of `q(i)` must not pair with the data `i`.
            (
                call(w, &call(p, &i)) * call(q, &i),
                call(w, &call(p, &i)) * call(q, &i),
            ),
            (
                call(w, &call(p, &i)) * call(q, &i) * call(p, &j) * call(r, &j),
                call(w, &call(p, &i)) * call(q, &i) * call(p, &d0) * call(r, &d0),
            ),
        ];
        let _ = &d1;
        // Several opaque calls in inputs that differ only by tagged dummy
        // names canonicalize identically: placeholders follow the literals.
        let k = slot(AbstractIndex::Normal(87_103));
        let pair = |index: &Atom| {
            call(w, &call(p, &i)) * call(u0, &call(q, &i)) * call(p, index) * call(r, index)
                + call(u1, &call(r, &i)) * call(q, index) * call(r, index)
        };
        assert_eq!(
            SymbolicTensor::infer(pair(&j)).unwrap().canonize().unwrap(),
            SymbolicTensor::infer(pair(&k)).unwrap().canonize().unwrap()
        );
        for (source, expected) in cases {
            let typed = SymbolicTensor::infer(source.clone()).unwrap();
            let canonical = typed.canonize().unwrap();
            assert_eq!(canonical.expression(), &expected, "{source}");
            assert_eq!(canonical.structure(), typed.structure(), "{source}");
            assert_eq!(canonical.canonize().unwrap(), canonical, "{source}");
            // No placeholder escapes into the result.
            assert!(
                !canonical
                    .expression()
                    .get_all_symbols(true)
                    .iter()
                    .any(|symbol| symbol.get_name().contains("canonical_opaque")),
                "{source}"
            );
        }
    }

    #[test]
    fn typed_canonize_preserves_open_color_duality() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::from(crate::representations::ColorFundamental {}).new_rep(3);
        let slots = [87001, 87002, 87003, 87004]
            .map(|index| rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index)));
        let direct = SymbolicTensor::infer(
            spenso::g!(slots[0].to_atom(), slots[1].dual().to_atom())
                * spenso::g!(slots[2].to_atom(), slots[3].dual().to_atom()),
        )
        .unwrap();
        let crossed = SymbolicTensor::infer(
            spenso::g!(slots[0].to_atom(), slots[3].dual().to_atom())
                * spenso::g!(slots[2].to_atom(), slots[1].dual().to_atom()),
        )
        .unwrap();
        for source in [&direct, &crossed] {
            let result = source.canonize().unwrap();
            assert_eq!(result.structure(), source.structure());
            assert_eq!(result.expression(), source.expression());
            assert_eq!(result.canonize().unwrap(), result);
        }
        assert_ne!(
            direct.canonize().unwrap().expression(),
            crossed.canonize().unwrap().expression()
        );
    }

    #[test]
    fn typed_canonize_keeps_dual_endpoints_when_relabeling_color_dummies() {
        crate::test_support::test_initialize();
        let fundamental = LibraryRep::from(crate::representations::ColorFundamental {}).new_rep(3);
        let adjoint = LibraryRep::from(crate::representations::ColorAdjoint {}).new_rep(8);
        let left = fundamental.slot::<AbstractIndex, _>(AbstractIndex::Normal(87101));
        let right = fundamental
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(87102))
            .dual();
        let mut canonical = None;
        for label in [87103, 87104] {
            let a = adjoint.slot::<AbstractIndex, _>(AbstractIndex::Normal(label));
            let internal = fundamental.slot::<AbstractIndex, _>(AbstractIndex::Normal(label));
            let source = SymbolicTensor::infer(
                crate::color_t!(a.to_atom(), left.to_atom(), internal.dual().to_atom())
                    * crate::color_t!(a.to_atom(), internal.to_atom(), right.to_atom()),
            )
            .unwrap();
            let result = source.canonize().unwrap();
            assert_eq!(result.structure(), source.structure());
            assert_eq!(result.canonize().unwrap(), result);
            if let Some(expected) = &canonical {
                assert_eq!(result.expression(), expected);
            }
            canonical = Some(result.expression().clone());
        }
    }

    #[test]
    fn canonical_generic_leaves_retain_raw_ports_without_inference() {
        test_support::test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let [a, b, internal] = [86111, 86112, 86113].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let expression = function!(
            tensor_symbol!("canonical_generic_loop"),
            &b,
            &internal,
            &a,
            &internal
        );
        let expected = expression
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap();
        let expected_execution = expected
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap()
            .simple_execute::<()>()
            .unwrap();
        let source = SymbolicTensor::infer(expression).unwrap();
        let source_leaves = source.admitted_leaf_interfaces().lock().unwrap().clone();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let canonical = source.canonize().unwrap();
        let network = canonical
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), source.rank());
        let executed = network.simple_execute::<()>().unwrap();
        let second = canonical.canonize().unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
        assert_eq!(canonical.expression(), &expected);
        assert_eq!(canonical.structure(), source.structure());
        assert_eq!(executed, expected_execution);
        assert_eq!(second, canonical);
        assert_eq!(
            *source.admitted_leaf_interfaces().lock().unwrap(),
            source_leaves
        );
    }

    #[test]
    fn typed_canonize_preserves_declared_order_and_factorized_spectators() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let slots = [86101, 86102, 86103]
            .map(|index| rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index)));
        let t = function!(
            spenso::tensor_symbol!("typed_canonize_T"),
            slots[0].to_atom(),
            slots[1].to_atom()
        );
        let p = function!(
            spenso::vector_symbol!("typed_canonize_p"),
            slots[2].to_atom()
        );
        let q = function!(
            spenso::vector_symbol!("typed_canonize_q"),
            slots[2].to_atom()
        );
        let spectator = (Atom::var(symbol!("typed_canonize_x")) + Atom::one()).pow(7);
        let expression = t * p * q * &spectator;
        let layout = PartialStructure::from_logical_slots(
            [slots[1], slots[0]].map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind))),
        );
        let source = SymbolicTensor::from_produced_parts(expression.clone(), layout).unwrap();
        let result = source.canonize().unwrap();
        assert_eq!(
            result.expression(),
            &expression
                .canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap()
        );
        assert_eq!(result.structure(), source.structure());
        assert!(result.expression().contains(&spectator));
        assert_eq!(result.canonize().unwrap(), result);
    }

    #[test]
    fn typed_canonize_preserves_auto_component_assignment_and_typed_zero() {
        crate::test_support::test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let head = spenso::tensor_symbol!("typed_canonize_auto_T");
        let source =
            SymbolicTensor::infer(function!(head, rep.to_symbolic([]), rep.to_symbolic([])))
                .unwrap();
        let result = source.canonize().unwrap();
        let bindings = HashMap::from([
            (0, AbstractIndex::Normal(86201)),
            (1, AbstractIndex::Normal(86202)),
        ]);
        assert_eq!(result.structure(), source.structure());
        assert_eq!(
            result.materialize_interface_ports(&bindings).unwrap(),
            source.materialize_interface_ports(&bindings).unwrap()
        );
        let zero = source.with_rewritten_expression(Atom::Zero).unwrap();
        assert_eq!(zero.canonize().unwrap(), zero);
    }

    #[test]
    fn typed_canonize_preserves_compliant_callback_interface() {
        test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        // Dummy indices serialize as ordinary named symbols; parsing that atom
        // yields Symbol, not the in-memory Dummy enum variant.
        let canonical_slot = LibraryRep::from(Minkowski {})
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Dummy(0))
            .to_atom();
        let head = tensor_symbol!(
            "typed_canonize_callback_T",
            norm = move |view, out| {
                if let AtomView::Fun(function) = view
                    && function.iter().any(|arg| arg == canonical_slot.as_view())
                {
                    seen.fetch_add(1, Ordering::SeqCst);
                    **out = Atom::num(2)
                        * FunctionBuilder::new(tensor_symbol!("canonize_compliant_result"))
                            .add_args(function.iter())
                            .finish();
                }
            }
        );
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let a = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(86301))
            .to_atom();
        let b = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(86302))
            .to_atom();
        let expression = FunctionBuilder::new(head).add_arg(&a).add_arg(&b).finish()
            * function!(vector_symbol!("typed_canonize_callback_p"), a);
        let source = SymbolicTensor::infer(expression).unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        let canonical = source.canonize().unwrap();
        assert_eq!(canonical.structure(), source.structure());
        assert!(
            source
                .with_checked_expression(canonical.expression().clone())
                .is_ok()
        );
        assert_eq!(
            SymbolicTensor::infer(canonical.expression().clone())
                .unwrap()
                .structure(),
            canonical.structure()
        );
        assert!(calls.load(Ordering::SeqCst) > 0);
    }
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "normalization callback contract violated in canonize")]
    fn canonize_checks_callback_contract() {
        test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = Arc::clone(&calls);
        // Dummy indices serialize as ordinary named symbols; parsing that atom
        // yields Symbol, not the in-memory Dummy enum variant.
        let canonical_slot = LibraryRep::from(Minkowski {})
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Dummy(0))
            .to_atom();
        let head = tensor_symbol!(
            "typed_canonize_callback_T_violation",
            norm = move |view, out| {
                if let AtomView::Fun(function) = view
                    && function.iter().any(|arg| arg == canonical_slot.as_view())
                {
                    seen.fetch_add(1, Ordering::SeqCst);
                    **out = Atom::one();
                }
            }
        );
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let a = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(86301))
            .to_atom();
        let b = rep
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(86302))
            .to_atom();
        let expression = FunctionBuilder::new(head).add_arg(&a).add_arg(&b).finish()
            * function!(vector_symbol!("typed_canonize_callback_p_violation"), a);
        let source = SymbolicTensor::infer(expression).unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        source.canonize().unwrap();
    }
}
