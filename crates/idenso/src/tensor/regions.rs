//! Interfaces retained at admission and composed by known arithmetic producers.

#[cfg(test)]
use std::cell::Cell;
use std::hash::{Hash, Hasher};
use std::ptr;
use std::slice::Iter;
use std::sync::{Arc, Mutex, OnceLock, Weak};

use ahash::{AHashMap, AHashSet, AHasher};
use spenso::network::library::symbolic::ETS;
use spenso::network::parsing::LeafInterfaceCache;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::LibraryRep;
use spenso::structure::slot::{IsAbstractSlot, SlotMatcher};
use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder};

use super::composition;
use super::inference::{InterfaceInference, TensorInferenceError};
use super::{SymbolicTensor, TensorProofs};
use crate::color::simplify::ColorAlgebraSimplifier;

#[cfg(test)]
thread_local! {
    pub(crate) static REGION_DERIVATIONS: Cell<usize> = const { Cell::new(0) };
}

/// Exact serialized bytes, shared at expression-child boundaries. Headers and
/// interstitial bytes remain part of equality, including normalization flags.
#[derive(Debug)]
enum AtomKeyPart {
    Bytes(Arc<[u8]>),
    Child(Arc<AtomKey>),
}

#[derive(Debug)]
struct AtomKey {
    fingerprint: u64,
    bytes: usize,
    parts: Vec<AtomKeyPart>,
    // Compact powers may synthesize views outside their serialized byte range.
    // Keep those semantic children reachable without using them for equality.
    detached_children: Vec<Arc<Self>>,
}

struct AtomKeyBytes<'a> {
    pending: Vec<Iter<'a, AtomKeyPart>>,
}

impl<'a> AtomKeyBytes<'a> {
    fn new(key: &'a AtomKey) -> Self {
        Self {
            pending: vec![key.parts.iter()],
        }
    }
}

impl<'a> Iterator for AtomKeyBytes<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.pending.last_mut()?.next() {
                Some(AtomKeyPart::Bytes(bytes)) => return Some(bytes),
                Some(AtomKeyPart::Child(child)) => self.pending.push(child.parts.iter()),
                None => {
                    self.pending.pop();
                }
            }
        }
    }
}

impl PartialEq for AtomKey {
    fn eq(&self, other: &Self) -> bool {
        if ptr::eq(self, other) {
            return true;
        }
        if self.fingerprint != other.fingerprint || self.bytes != other.bytes {
            return false;
        }
        let mut left = AtomKeyBytes::new(self);
        let mut right = AtomKeyBytes::new(other);
        let mut a = &[][..];
        let mut b = &[][..];
        loop {
            if a.is_empty() {
                a = left.next().unwrap_or_default();
            }
            if b.is_empty() {
                b = right.next().unwrap_or_default();
            }
            if a.is_empty() || b.is_empty() {
                return a.is_empty() && b.is_empty();
            }
            let count = a.len().min(b.len());
            if a[..count] != b[..count] {
                return false;
            }
            a = &a[count..];
            b = &b[count..];
        }
    }
}

impl Eq for AtomKey {}

impl AtomKey {
    fn children(&self) -> impl Iterator<Item = &Arc<Self>> {
        self.parts
            .iter()
            .filter_map(|part| match part {
                AtomKeyPart::Child(child) => Some(child),
                AtomKeyPart::Bytes(_) => None,
            })
            .chain(&self.detached_children)
    }

    fn fingerprint(value: &[u8]) -> u64 {
        let mut hasher = AHasher::default();
        value.hash(&mut hasher);
        hasher.finish()
    }

    fn matches(&self, value: &[u8]) -> bool {
        if let [AtomKeyPart::Bytes(bytes)] = self.parts.as_slice() {
            return bytes.as_ref() == value;
        }
        if self.bytes != value.len() {
            return false;
        }
        let mut offset = 0;
        for bytes in AtomKeyBytes::new(self) {
            if value.get(offset..offset + bytes.len()) != Some(bytes) {
                return false;
            }
            offset += bytes.len();
        }
        offset == value.len()
    }
}

#[derive(Debug, Default)]
struct AtomKeys {
    buckets: AHashMap<u64, Vec<Weak<AtomKey>>>,
    sweep_at: usize,
}

impl AtomKeys {
    /// Index only this payload's shared keys. Each carrier owns its table so
    /// later producers cannot leave dead control blocks on an earlier carrier.
    fn from_root(root: Option<&Arc<AtomKey>>) -> Self {
        let mut result = Self::default();
        let mut pending = root.into_iter().collect::<Vec<_>>();
        let mut visited = AHashSet::new();
        while let Some(key) = pending.pop() {
            if !visited.insert(Arc::as_ptr(key)) {
                continue;
            }
            result
                .buckets
                .entry(key.fingerprint)
                .or_default()
                .push(Arc::downgrade(key));
            pending.extend(key.children());
        }
        result.sweep_at = result.buckets.len().saturating_mul(2).max(1024);
        result
    }

    fn publish(&mut self, key: AtomKey) -> Arc<AtomKey> {
        let key = Arc::new(key);
        let bucket = self.buckets.entry(key.fingerprint).or_default();
        bucket.retain(|key| key.strong_count() != 0);
        bucket.push(Arc::downgrade(&key));
        if self.buckets.len() > self.sweep_at.max(1024) {
            self.buckets.retain(|_, keys| {
                keys.retain(|key| key.strong_count() != 0);
                !keys.is_empty()
            });
            self.sweep_at = self.buckets.len().saturating_mul(2).max(1024);
        }
        key
    }

    fn intern(&mut self, value: AtomView<'_>, fingerprint: u64) -> Arc<AtomKey> {
        let bytes = value.get_data();
        if let Some(key) = self.buckets.get(&fingerprint).and_then(|bucket| {
            bucket
                .iter()
                .filter_map(Weak::upgrade)
                .find(|key| key.matches(bytes))
        }) {
            return key;
        }
        let start = bytes.as_ptr() as usize;
        let mut offset = 0;
        let mut parts = Vec::new();
        for child in value.children() {
            let child_bytes = child.get_data();
            let Some(child_offset) = (child_bytes.as_ptr() as usize).checked_sub(start) else {
                parts.clear();
                break;
            };
            let Some(end) = child_offset.checked_add(child_bytes.len()) else {
                parts.clear();
                break;
            };
            if child_offset < offset || end > bytes.len() {
                // Compact powers may synthesize a child view. Their complete
                // small encoding remains an exact leaf rather than a range.
                parts.clear();
                break;
            }
            if child_offset != offset {
                parts.push(AtomKeyPart::Bytes(Arc::from(&bytes[offset..child_offset])));
            }
            parts.push(AtomKeyPart::Child(
                self.intern(child, AtomKey::fingerprint(child_bytes)),
            ));
            offset = end;
        }
        let detached_children = if parts.is_empty() {
            parts.push(AtomKeyPart::Bytes(Arc::from(bytes)));
            value
                .children()
                .map(|child| self.intern(child, AtomKey::fingerprint(child.get_data())))
                .collect()
        } else {
            if offset < bytes.len() {
                parts.push(AtomKeyPart::Bytes(Arc::from(&bytes[offset..])));
            }
            Vec::new()
        };
        self.publish(AtomKey {
            fingerprint,
            bytes: bytes.len(),
            parts,
            detached_children,
        })
    }
}

#[derive(Debug, Default)]
struct RegionInterfaces {
    buckets: AHashMap<u64, Vec<(Arc<AtomKey>, PartialStructure)>>,
}

impl RegionInterfaces {
    fn get(&self, bytes: &[u8], fingerprint: u64) -> Option<(Arc<AtomKey>, PartialStructure)> {
        self.buckets
            .get(&fingerprint)?
            .iter()
            .find_map(|(key, interface)| {
                key.matches(bytes)
                    .then(|| (Arc::clone(key), interface.clone()))
            })
    }

    fn insert(&mut self, key: Arc<AtomKey>, interface: PartialStructure) {
        let bucket = self.buckets.entry(key.fingerprint).or_default();
        if let Some((_, existing)) = bucket.iter_mut().find(|(existing, _)| *existing == key) {
            *existing = interface;
        } else {
            bucket.push((key, interface));
        }
    }

    fn interface(&self, key: &AtomKey) -> Option<&PartialStructure> {
        self.buckets
            .get(&key.fingerprint)?
            .iter()
            .find_map(|(candidate, interface)| (candidate.as_ref() == key).then_some(interface))
    }
}

/// The current payload owns its admitted regions. A producer shares keys for
/// unchanged subtrees, but never keeps predecessor carriers or obsolete roots.
#[derive(Debug, Default)]
pub(crate) struct AdmittedInterfaces {
    observed: RegionInterfaces,
    derived: Mutex<RegionInterfaces>,
    keys: Mutex<AtomKeys>,
    root: Option<Arc<AtomKey>>,
}

#[derive(Clone, Copy)]
enum MappingOrder {
    BeforeChildren,
    TensorArguments,
    AfterChildren,
}

struct RegionMapping<'a, F, G> {
    order: MappingOrder,
    source: &'a SymbolicTensor<PartialStructure>,
    leaf_map: F,
    interface_map: G,
    interfaces: AdmittedInterfaces,
    source_leaves: LeafInterfaceCache,
    leaves: LeafInterfaceCache,
    written_slots: AHashSet<Atom>,
    mapped_slots: AHashMap<Atom, Atom>,
}

impl<F, G> RegionMapping<'_, F, G> {
    fn map<'a, E>(&mut self, value: AtomView<'a>) -> Result<AtomOrView<'a>, E>
    where
        F: for<'b> FnMut(AtomView<'b>) -> Result<Option<Atom>, E>,
        G: FnMut(&PartialStructure) -> Result<PartialStructure, TensorInferenceError>,
        E: From<TensorInferenceError>,
    {
        let replacement = match self.order {
            MappingOrder::BeforeChildren | MappingOrder::TensorArguments => (self.leaf_map)(value)?,
            MappingOrder::AfterChildren => None,
        };
        let mut mapped_arguments = None;
        let mut producer_changed_head = false;
        let mut result = if let Some(replacement) = replacement {
            AtomOrView::Atom(replacement)
        } else {
            match value {
                AtomView::Num(_) | AtomView::Var(_) => AtomOrView::View(value),
                _ => {
                    let mut matcher = SlotMatcher::default();
                    let children = value
                        .children()
                        .map(|child| {
                            if matches!(self.order, MappingOrder::TensorArguments)
                                && let AtomView::Fun(function) = value
                                && !InterfaceInference::explicit_index_argument(
                                    function.get_symbol(),
                                    child,
                                    &mut matcher,
                                )
                            {
                                Ok(AtomOrView::View(child))
                            } else {
                                self.map(child)
                            }
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    if children
                        .iter()
                        .all(|child| matches!(child, AtomOrView::View(_)))
                    {
                        AtomOrView::View(value)
                    } else {
                        let expression = match value {
                            AtomView::Fun(function) => {
                                let expression = FunctionBuilder::new(function.get_symbol())
                                    .add_args(children.iter().map(AtomCore::as_atom_view))
                                    .finish();
                                mapped_arguments = Some(children);
                                expression
                            }
                            AtomView::Add(_) => Atom::add_many(children),
                            AtomView::Mul(_) => Atom::mul_many(children),
                            AtomView::Pow(_) => {
                                let mut children = children.into_iter();
                                let base = children.next().expect("a power has a base");
                                let exponent = children.next().expect("a power has an exponent");
                                base.as_view().pow(exponent.as_view())
                            }
                            _ => unreachable!("leaf nodes are handled before child mapping"),
                        };
                        AtomOrView::Atom(expression)
                    }
                }
            }
        };
        if matches!(self.order, MappingOrder::AfterChildren)
            && let Some(replacement) = (self.leaf_map)(result.as_view())?
        {
            producer_changed_head = match (result.as_view(), replacement.as_view()) {
                (AtomView::Fun(source), AtomView::Fun(target)) => {
                    source.get_symbol() != target.get_symbol()
                }
                _ => false,
            };
            result = AtomOrView::Atom(replacement);
        }
        if self.written_slots.contains(value.get_data()) {
            self.mapped_slots
                .insert(value.to_owned(), result.as_view().to_owned());
        }
        if matches!(result, AtomOrView::View(_)) {
            let source_leaves = self
                .source_leaves
                .lock()
                .expect("leaf interface cache is not poisoned");
            if let Some((key, facts)) = source_leaves.get_key_value(value.get_data()) {
                self.leaves
                    .lock()
                    .expect("leaf interface cache is not poisoned")
                    .insert(key.clone(), facts.clone());
            }
            return Ok(result);
        }
        if let Some(interface) = self.source.admitted_interface(value) {
            self.interfaces
                .insert(result.as_view(), (self.interface_map)(&interface)?);
        }
        if let AtomView::Fun(source) = value {
            let previous = self
                .source_leaves
                .lock()
                .expect("leaf interface cache is not poisoned")
                .get(value.get_data())
                .cloned();
            if let Some((tensorial, slots)) = previous
                && let Some(mapped) = slots
                    .iter()
                    .map(|slot| self.mapped_slots.get(slot).cloned())
                    .collect::<Option<Vec<_>>>()
            {
                let supplied = mapped_arguments.as_ref().map_or_else(
                    || source.iter().collect::<Vec<_>>(),
                    |arguments| arguments.iter().map(AtomCore::as_atom_view).collect(),
                );
                if let Some((target, raw)) = TensorProofs::raw_ports_after_argument_permutation(
                    &supplied,
                    result.as_view(),
                    &mapped,
                ) && (source.get_symbol() == target.get_symbol() || producer_changed_head)
                {
                    self.leaves
                        .lock()
                        .expect("leaf interface cache is not poisoned")
                        .insert(target.as_view().to_owned(), (tensorial, raw));
                }
                if source.get_symbol().has_tag(&SPENSO_TAG.rank1)
                    && let [port] = mapped.as_slice()
                    && let Some(interface) = self.source.admitted_interface(value)
                {
                    let interface = (self.interface_map)(&interface)?;
                    let mut pending = vec![result.as_view()];
                    let mut matcher = SlotMatcher::default();
                    while let Some(node) = pending.pop() {
                        match node {
                            AtomView::Add(sum) => pending.extend(sum.iter()),
                            AtomView::Mul(product) => pending.extend(product.iter()),
                            AtomView::Pow(power) => pending.push(power.get_base_exp().0),
                            AtomView::Fun(function)
                                if function.get_symbol().has_tag(&SPENSO_TAG.rank1)
                                    && matcher.vector_argument(function)
                                        == Some(port.as_view()) =>
                            {
                                // The producer carries the same single port to
                                // every emitted vector. Constructor semantics
                                // verify that witness; no output shape is inferred.
                                self.interfaces.insert(node, interface.clone());
                                self.leaves
                                    .lock()
                                    .expect("leaf interface cache is not poisoned")
                                    .insert(node.to_owned(), (true, vec![port.clone()]));
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
        Ok(result)
    }
}

impl TensorProofs {
    /// Symmetry normalization may permute the constructor's known arguments and
    /// factor out a numerical sign. Match its surviving leaf exactly and reorder
    /// only the supplied direct ports; the caller verifies the constructor head.
    pub(crate) fn raw_ports_after_argument_permutation<'a>(
        supplied: &[AtomView<'_>],
        target: AtomView<'a>,
        ports: &[Atom],
    ) -> Option<(FunView<'a>, Vec<Atom>)> {
        let target = match target {
            AtomView::Fun(function) => function,
            AtomView::Mul(product) if product.has_coefficient() && product.get_nargs() == 2 => {
                match product
                    .iter()
                    .find(|factor| !matches!(factor, AtomView::Num(_)))?
                {
                    AtomView::Fun(function) => function,
                    _ => return None,
                }
            }
            _ => return None,
        };
        if target.iter().eq(supplied.iter().copied()) {
            return Some((target, ports.to_vec()));
        }
        let mut remaining = supplied.to_vec();
        for argument in target.iter() {
            let position = remaining.iter().position(|&source| source == argument)?;
            remaining.swap_remove(position);
        }
        if !remaining.is_empty() {
            return None;
        }
        let mut remaining = ports.to_vec();
        let mut reordered = Vec::with_capacity(ports.len());
        for argument in target.iter() {
            if let Some(position) = remaining.iter().position(|port| port.as_view() == argument) {
                reordered.push(remaining.remove(position));
            }
        }
        remaining.is_empty().then_some((target, reordered))
    }
}

impl AdmittedInterfaces {
    fn retain_leaves(&self, sources: &[LeafInterfaceCache]) -> LeafInterfaceCache {
        let mut pending = self.root.iter().cloned().collect::<Vec<_>>();
        let mut visited = AHashSet::new();
        let mut reachable: AHashMap<u64, Vec<Arc<AtomKey>>> = AHashMap::new();
        while let Some(key) = pending.pop() {
            if !visited.insert(Arc::as_ptr(&key)) {
                continue;
            }
            pending.extend(key.children().map(Arc::clone));
            reachable.entry(key.fingerprint).or_default().push(key);
        }
        let result = LeafInterfaceCache::default();
        let mut retained = result.lock().expect("leaf interface cache is not poisoned");
        for source in sources {
            let source = source.lock().expect("leaf interface cache is not poisoned");
            for (atom, facts) in source.iter() {
                if retained.contains_key(atom.as_view().get_data()) {
                    continue;
                }
                let bytes = atom.as_view().get_data();
                if reachable
                    .get(&AtomKey::fingerprint(bytes))
                    .is_some_and(|keys| keys.iter().any(|key| key.matches(bytes)))
                {
                    retained.insert(atom.clone(), facts.clone());
                }
            }
        }
        drop(retained);
        result
    }

    pub(crate) fn from_region(value: AtomView<'_>, interface: PartialStructure) -> Self {
        let mut result = Self::default();
        result.insert(value, interface);
        result
    }

    pub(crate) fn produced(
        value: AtomView<'_>,
        interface: Option<PartialStructure>,
        previous: Option<&Arc<Self>>,
    ) -> Arc<Self> {
        let mut result = Self {
            keys: Mutex::new(AtomKeys::from_root(
                previous.and_then(|source| source.root.as_ref()),
            )),
            ..Self::default()
        };
        result.insert_root(value, interface);
        if let Some(previous) = previous {
            result.import_reachable([previous.as_ref()]);
        }
        result.keys = Mutex::new(AtomKeys::from_root(result.root.as_ref()));
        Arc::new(result)
    }

    pub(crate) fn insert(&mut self, value: AtomView<'_>, interface: PartialStructure) {
        self.insert_root(value, Some(interface));
    }

    fn insert_root(&mut self, value: AtomView<'_>, interface: Option<PartialStructure>) {
        if !value.is_zero() {
            let key = self
                .keys
                .lock()
                .expect("regional key interner is not poisoned")
                .intern(value, AtomKey::fingerprint(value.get_data()));
            self.root = Some(Arc::clone(&key));
            if let Some(interface) = interface {
                self.observed.insert(key, interface);
            }
        }
    }

    fn local(&self, bytes: &[u8], fingerprint: u64) -> Option<(Arc<AtomKey>, PartialStructure)> {
        self.observed.get(bytes, fingerprint).or_else(|| {
            self.derived
                .lock()
                .expect("derived interface cache is not poisoned")
                .get(bytes, fingerprint)
        })
    }

    fn get(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        self.local(value.get_data(), AtomKey::fingerprint(value.get_data()))
            .map(|(_, interface)| interface)
    }

    fn retain(&self, value: AtomView<'_>, interface: &PartialStructure) {
        if !value.is_zero() {
            let key = self
                .keys
                .lock()
                .expect("regional key interner is not poisoned")
                .intern(value, AtomKey::fingerprint(value.get_data()));
            self.derived
                .lock()
                .expect("derived interface cache is not poisoned")
                .insert(key, interface.clone());
        }
    }

    /// Copy facts whose exact keys belong to the current payload. Walking
    /// interned child edges does not read tensor syntax or infer boundaries.
    /// Existing facts win, followed by sources in their given order.
    fn import_reachable<'a>(&mut self, sources: impl IntoIterator<Item = &'a Self>) {
        let Some(root) = &self.root else { return };
        let mut pending = vec![Arc::clone(root)];
        let mut visited = AHashSet::new();
        let mut reachable = Vec::new();
        while let Some(key) = pending.pop() {
            if !visited.insert(Arc::as_ptr(&key)) {
                continue;
            }
            pending.extend(key.children().map(Arc::clone));
            if self.observed.interface(&key).is_none() {
                reachable.push(key);
            }
        }
        for source in sources {
            let derived = source
                .derived
                .lock()
                .expect("derived interface cache is not poisoned");
            reachable.retain(|key| {
                if let Some(interface) = source
                    .observed
                    .interface(key)
                    .or_else(|| derived.interface(key))
                {
                    self.observed.insert(Arc::clone(key), interface.clone());
                    false
                } else {
                    true
                }
            });
            if reachable.is_empty() {
                break;
            }
        }
    }

    fn inherit(sources: Vec<Arc<Self>>) -> Arc<Self> {
        let first = &sources[0];
        let mut result = Self {
            keys: Mutex::new(AtomKeys::from_root(first.root.as_ref())),
            root: first.root.clone(),
            ..Self::default()
        };
        result.import_reachable(sources.iter().map(Arc::as_ref));
        Arc::new(result)
    }
}

impl TensorProofs {
    /// Whether these facts belong to this exact current payload, including its
    /// normalization flags. This does not infer or derive any interface.
    pub(crate) fn owns_expression(&self, value: AtomView<'_>) -> bool {
        self.admitted_interfaces
            .get()
            .and_then(|interfaces| interfaces.root.as_ref())
            .is_some_and(|root| root.matches(value.get_data()))
    }

    /// Publish a new authoritative root and retain only its reachable regional
    /// and leaf facts. Input carriers remain independently usable.
    pub(crate) fn publish_interfaces<'a>(
        &mut self,
        value: AtomView<'_>,
        interface: impl Into<Option<PartialStructure>>,
        sources: impl IntoIterator<Item = &'a Self>,
    ) {
        self.admitted_interfaces = OnceLock::from(AdmittedInterfaces::produced(
            value,
            interface.into(),
            self.admitted_interfaces.get(),
        ));
        self.inherit_interfaces(sources);
    }

    /// Fetch only previously supplied facts, without deriving any new region.
    pub(crate) fn retained_interface(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        self.admitted_interfaces.get()?.get(value)
    }

    /// Record the interface supplied by a leaf producer. Publication retains
    /// this fact only while its exact expression belongs to the output tree.
    pub(crate) fn retain_interface(&self, value: AtomView<'_>, interface: PartialStructure) {
        if matches!(value, AtomView::Fun(function)
            if composition::is_composite_head(function.get_symbol(), &SPENSO_TAG)
                || function.get_symbol() == SPENSO_TAG.dot)
        {
            self.leaf_interfaces
                .get_or_init(Default::default)
                .lock()
                .expect("leaf interface cache is not poisoned")
                .insert(
                    value.to_owned(),
                    (
                        true,
                        interface
                            .logical_slots()
                            .into_iter()
                            .map(composition::port_atom)
                            .collect(),
                    ),
                );
        }
        self.admitted_interfaces
            .get_or_init(|| Arc::new(AdmittedInterfaces::default()))
            .retain(value, &interface);
    }

    /// Compose a newly constructed metric from its admitted endpoint roles.
    /// Explicit slots survive; representations create positional axes; vector
    /// operands are consumed inside the metric. This never examines its output.
    pub(crate) fn metric_interface(
        &self,
        arguments: [AtomView<'_>; 2],
    ) -> Result<PartialStructure, TensorInferenceError> {
        let mut matcher = SlotMatcher::default();
        let mut logical = Vec::new();
        for argument in arguments {
            if let Ok(slot) = matcher.parse::<LibraryRep, AbstractIndex>(argument) {
                logical.push(slot.rep().slot(PartialIndex::Explicit(slot.aind())));
            } else if let Ok(representation) = matcher.parse_representation::<LibraryRep>(argument)
            {
                logical.push(representation.slot(PartialIndex::open(logical.len())));
            } else if !matches!(argument, AtomView::Num(_) | AtomView::Var(_))
                && matcher.concrete_component(argument).is_none()
                && !matches!(argument, AtomView::Fun(operand) if operand.get_symbol().has_tag(&SPENSO_TAG.rank1))
            {
                let interface = if let Some(interface) = self.admitted_interface(argument) {
                    interface
                } else if let AtomView::Fun(operand) = argument
                    && composition::is_tensor_leaf_head(operand.get_symbol())
                    && operand.get_nargs() == 1
                    && let Some(port) = operand.iter().next()
                    && matcher.parse_representation::<LibraryRep>(port).is_ok()
                {
                    // A gamma pairing may expose a generic one-port vector
                    // formerly consumed inside the word. The metric supplies
                    // its single compact endpoint and consumes this entire axis.
                    self.retain_leaf_ports(argument, vec![port.to_owned()])?;
                    continue;
                } else {
                    return Err(TensorInferenceError::invalid(format!(
                        "metric operand has no retained spectator interface for {argument}"
                    )));
                };
                let (_, spectators) =
                    InterfaceInference::consume_inner_product_axis(&interface, None)?;
                logical.extend(spectators.logical_slots());
            }
        }
        InterfaceInference::merge_explicit_interface_sequence(&[
            PartialStructure::from_logical_slots(logical),
        ])
    }

    pub(crate) fn inherit_interfaces<'a>(&mut self, sources: impl IntoIterator<Item = &'a Self>) {
        let mut leaves = self
            .leaf_interfaces
            .get()
            .cloned()
            .into_iter()
            .collect::<Vec<_>>();
        let mut predecessors = self
            .admitted_interfaces
            .get()
            .cloned()
            .into_iter()
            .collect::<Vec<_>>();
        for source in sources {
            if let Some(source_leaves) = source.leaf_interfaces.get()
                && !leaves
                    .iter()
                    .any(|existing| Arc::ptr_eq(existing, source_leaves))
            {
                leaves.push(Arc::clone(source_leaves));
            }
            if let Some(source) = source.admitted_interfaces.get()
                && !predecessors
                    .iter()
                    .any(|existing| Arc::ptr_eq(existing, source))
            {
                predecessors.push(Arc::clone(source));
            }
        }
        self.admitted_interfaces = match predecessors.len() {
            0 => OnceLock::new(),
            1 => OnceLock::from(predecessors.pop().expect("one retained admission cache")),
            _ => OnceLock::from(AdmittedInterfaces::inherit(predecessors)),
        };
        if let Some(interfaces) = self.admitted_interfaces.get() {
            self.leaf_interfaces = OnceLock::from(interfaces.retain_leaves(&leaves));
        }
    }

    /// Look up an admitted region, or compose arithmetic from admitted children.
    /// Unknown function syntax is never inspected to rediscover its ports.
    pub(crate) fn admitted_interface(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        if let Some(interface) = self
            .admitted_interfaces
            .get()
            .and_then(|interfaces| interfaces.get(value))
        {
            return Some(interface);
        }
        #[cfg(test)]
        REGION_DERIVATIONS.with(|count| count.set(count.get() + 1));
        let leaf = self.leaf_interfaces.get().and_then(|leaves| {
            leaves
                .lock()
                .expect("leaf interface cache is not poisoned")
                .get(value.get_data())
                .cloned()
        });
        if let Some((tensorial, ports)) = leaf {
            if tensorial {
                // Materializers already supplied these exact raw ports. Merge
                // their incidences without inspecting the emitted leaf again.
                self.retain_leaf_ports(value, ports).ok()?;
                return self.retained_interface(value);
            }
            let scalar = PartialStructure::from_logical_slots([]);
            self.retain_interface(value, scalar.clone());
            return Some(scalar);
        }
        let interface = if self.observations.get().is_some_and(|observed| {
            observed
                .region(value)
                .is_some_and(|region| region.certifies_scalar_interface_syntax())
        }) {
            Some(PartialStructure::from_logical_slots([]))
        } else {
            match value {
                AtomView::Num(_) | AtomView::Var(_) => {
                    Some(PartialStructure::from_logical_slots([]))
                }
                AtomView::Add(sum) => {
                    let mut terms = sum.iter();
                    let interface = self.admitted_interface(terms.next()?)?;
                    terms
                        .all(|term| {
                            self.admitted_interface(term).is_some_and(|other| {
                                InterfaceInference::additive_interfaces_match(&interface, &other)
                            })
                        })
                        .then_some(interface)
                }
                AtomView::Mul(product) => {
                    let interfaces = product
                        .iter()
                        .map(|factor| self.admitted_interface(factor))
                        .collect::<Option<Vec<_>>>()?;
                    InterfaceInference::merge_explicit_interface_sequence(&interfaces).ok()
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let interface = self.admitted_interface(base)?;
                    InterfaceInference::power_interface(&interface, exponent, || value.to_string())
                        .ok()
                }
                AtomView::Fun(function)
                    if [SPENSO_TAG.chain, SPENSO_TAG.trace].contains(&function.get_symbol()) =>
                {
                    let mut args = function.iter();
                    let header = args.next()?;
                    let mut factors = Vec::new();
                    if function.get_symbol() == SPENSO_TAG.chain {
                        factors.push(self.metric_interface([header, args.next()?]).ok()?);
                    }
                    for factor in args {
                        factors.push(self.admitted_interface(factor)?);
                    }
                    InterfaceInference::merge_explicit_interface_sequence(&factors).ok()
                }
                // Brackets and the symmetry projectors carry their operands'
                // combined boundary. Unary broadcast heads preserve it.
                AtomView::Fun(function)
                    if composition::is_composite_head(function.get_symbol(), &SPENSO_TAG) =>
                {
                    let operands = function
                        .iter()
                        .map(|factor| self.admitted_interface(factor))
                        .collect::<Option<Vec<_>>>()?;
                    InterfaceInference::merge_explicit_interface_sequence(&operands).ok()
                }
                AtomView::Fun(function)
                    if function.get_symbol().is_scalar()
                        || function.get_symbol() == SPENSO_TAG.dot =>
                {
                    Some(PartialStructure::from_logical_slots([]))
                }
                AtomView::Fun(function)
                    if function.get_symbol() == ETS.metric && function.get_nargs() == 2 =>
                {
                    let mut arguments = function.iter();
                    self.metric_interface([arguments.next()?, arguments.next()?])
                        .ok()
                }
                AtomView::Fun(function) => {
                    if let Some(ports) = InterfaceInference::constructor_ports(
                        function,
                        &mut SlotMatcher::default(),
                        false,
                    ) {
                        self.retain_leaf_ports(value, ports).ok()?;
                        self.retained_interface(value)
                    } else {
                        ColorAlgebraSimplifier::produced_interface(value)
                    }
                }
            }
        }?;
        if let Some(cache) = self.admitted_interfaces.get() {
            cache.retain(value, &interface);
        }
        Some(interface)
    }

    /// Open a known graph boundary for a local rewrite. The enclosing carrier
    /// keeps its conservative construction history when this payload returns;
    /// regional results must not be exported as independently admitted values.
    pub(crate) fn region_with_interface(
        &self,
        value: AtomView<'_>,
        interface: PartialStructure,
    ) -> Option<SymbolicTensor<PartialStructure>> {
        let occurrences = Arc::clone(self.occurrences.get()?);
        let mut region =
            SymbolicTensor::from_validated_parts(value.to_owned(), interface, occurrences);
        region.proofs.inherit_interfaces([self]);
        Some(region)
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Transform payload and retained regional facts in the same traversal.
    /// Rebuilt nodes run their normalizers; metadata follows the producer's
    /// slot map and never comes from observing the rebuilt payload.
    pub(crate) fn map_admitted_expression(
        &self,
        leaf_map: impl for<'a> FnMut(AtomView<'a>) -> Result<Option<Atom>, TensorInferenceError>,
        interface_map: impl FnMut(&PartialStructure) -> Result<PartialStructure, TensorInferenceError>,
    ) -> Result<(Atom, TensorProofs), TensorInferenceError> {
        self.map_producer_draft(self.expression.as_view(), leaf_map, interface_map)
    }

    /// Finish a producer's intermediate syntax with the source's root contract
    /// and its retained facts for each rewritten region. The intermediate
    /// syntax need not itself be a valid tensor constructor.
    pub(crate) fn map_producer_draft(
        &self,
        expression: AtomView<'_>,
        leaf_map: impl for<'a> FnMut(AtomView<'a>) -> Result<Option<Atom>, TensorInferenceError>,
        interface_map: impl FnMut(&PartialStructure) -> Result<PartialStructure, TensorInferenceError>,
    ) -> Result<(Atom, TensorProofs), TensorInferenceError> {
        self.map_expression(
            expression,
            leaf_map,
            interface_map,
            MappingOrder::BeforeChildren,
        )
    }

    pub(crate) fn map_admitted_expression_bottom_up<E>(
        &self,
        node_map: impl for<'a> FnMut(AtomView<'a>) -> Result<Option<Atom>, E>,
        interface_map: impl FnMut(&PartialStructure) -> Result<PartialStructure, TensorInferenceError>,
    ) -> Result<(Atom, TensorProofs), E>
    where
        E: From<TensorInferenceError>,
    {
        self.map_expression(
            self.expression.as_view(),
            node_map,
            interface_map,
            MappingOrder::AfterChildren,
        )
    }

    /// Map structural indices while keeping tensor metadata and dimensions opaque.
    pub(crate) fn map_admitted_indices(
        &self,
        slot_map: impl for<'a> FnMut(AtomView<'a>) -> Result<Option<Atom>, TensorInferenceError>,
        interface_map: impl FnMut(&PartialStructure) -> Result<PartialStructure, TensorInferenceError>,
    ) -> Result<(Atom, TensorProofs), TensorInferenceError> {
        self.map_expression(
            self.expression.as_view(),
            slot_map,
            interface_map,
            MappingOrder::TensorArguments,
        )
    }

    fn map_expression<E>(
        &self,
        expression: AtomView<'_>,
        leaf_map: impl for<'a> FnMut(AtomView<'a>) -> Result<Option<Atom>, E>,
        interface_map: impl FnMut(&PartialStructure) -> Result<PartialStructure, TensorInferenceError>,
        order: MappingOrder,
    ) -> Result<(Atom, TensorProofs), E>
    where
        E: From<TensorInferenceError>,
    {
        let previous = self.proofs.admitted_interfaces.get();
        let source_leaves = self
            .proofs
            .leaf_interfaces
            .get()
            .cloned()
            .unwrap_or_default();
        let written_slots = source_leaves
            .lock()
            .expect("leaf interface cache is not poisoned")
            .values()
            .flat_map(|(_, slots)| slots.iter().cloned())
            .collect();
        let mut mapping = RegionMapping {
            order,
            source: self,
            leaf_map,
            interface_map,
            interfaces: AdmittedInterfaces {
                keys: Mutex::new(AtomKeys::from_root(
                    previous.and_then(|source| source.root.as_ref()),
                )),
                ..AdmittedInterfaces::default()
            },
            source_leaves,
            leaves: LeafInterfaceCache::default(),
            written_slots,
            mapped_slots: AHashMap::new(),
        };
        let expression = mapping.map(expression)?.into_owned();
        // Keep the declared root order even when normalization returns an
        // already observed child, then retain only still-reachable facts.
        mapping.interfaces.insert(
            expression.as_view(),
            (mapping.interface_map)(&self.structure)?,
        );
        let interfaces = AdmittedInterfaces::inherit(
            [Some(Arc::new(mapping.interfaces)), previous.cloned()]
                .into_iter()
                .flatten()
                .collect(),
        );
        let leaves = interfaces.retain_leaves(&[mapping.leaves, mapping.source_leaves]);
        Ok((
            expression,
            TensorProofs {
                admitted_interfaces: OnceLock::from(interfaces),
                leaf_interfaces: OnceLock::from(leaves),
                ..TensorProofs::default()
            },
        ))
    }

    pub(crate) fn admitted_interface(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        if value == self.expression.as_view() {
            Some(self.structure.clone())
        } else {
            self.proofs.admitted_interface(value)
        }
    }

    #[cfg(test)]
    pub(crate) fn admitted_region(
        &self,
        value: AtomView<'_>,
    ) -> Option<SymbolicTensor<PartialStructure>> {
        self.proofs
            .region_with_interface(value, self.admitted_interface(value)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tensor::SymbolicNetParse;
    use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS, SCOPE_VALIDATIONS};
    use crate::test_support;
    use spenso::network::parsing::ParseSettings;
    use spenso::shadowing::Concretize;
    use spenso::structure::representation::LibrarySlot;
    use spenso::structure::{Canonicalized, NamedStructure, OrderedStructure, TensorStructure};
    use spenso::{g, mink, tensor_symbol, vector_symbol};
    use std::collections::HashMap;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::FunctionBuilder;
    use symbolica::{function, symbol};

    #[derive(Debug, Default)]
    struct KeyStorage {
        nodes: usize,
        edges: usize,
        payload_bytes: usize,
        records: usize,
        interner_buckets: usize,
        interner_records: usize,
        dead_keys: usize,
    }

    impl AdmittedInterfaces {
        fn storage(&self) -> KeyStorage {
            let mut result = KeyStorage::default();
            let interner = self.keys.lock().unwrap();
            result.interner_buckets = interner.buckets.len();
            result.interner_records = interner.buckets.values().map(Vec::len).sum();
            result.dead_keys = interner
                .buckets
                .values()
                .flatten()
                .filter(|key| key.strong_count() == 0)
                .count();
            drop(interner);
            let mut keys = self
                .observed
                .buckets
                .values()
                .flatten()
                .map(|(key, _)| Arc::clone(key))
                .collect::<Vec<_>>();
            keys.extend(
                self.derived
                    .lock()
                    .unwrap()
                    .buckets
                    .values()
                    .flatten()
                    .map(|(key, _)| Arc::clone(key)),
            );
            result.records = keys.len();
            let mut visited_keys = AHashSet::new();
            let mut visited_bytes = AHashSet::new();
            while let Some(key) = keys.pop() {
                if !visited_keys.insert(Arc::as_ptr(&key)) {
                    continue;
                }
                result.nodes += 1;
                result.edges += key.detached_children.len();
                keys.extend(key.detached_children.iter().map(Arc::clone));
                for part in &key.parts {
                    match part {
                        AtomKeyPart::Bytes(bytes) => {
                            if visited_bytes.insert(bytes.as_ptr()) {
                                result.payload_bytes += bytes.len();
                            }
                        }
                        AtomKeyPart::Child(child) => {
                            result.edges += 1;
                            keys.push(Arc::clone(child));
                        }
                    }
                }
            }
            result
        }
    }

    #[test]
    fn exact_key_collisions_and_normalization_flags_are_checked() {
        let head = symbol!("regional_key_collision");
        let first = function!(head, 1);
        let second = function!(head, 2);
        let mut keys = AtomKeys::default();
        let a = keys.intern(first.as_view(), 0);
        let b = keys.intern(second.as_view(), 0);
        assert!(!Arc::ptr_eq(&a, &b));
        assert!(a.matches(first.as_view().get_data()));
        assert!(!a.matches(second.as_view().get_data()));
        assert!(Arc::ptr_eq(&a, &keys.intern(first.as_view(), 0)));

        let mut raw = Atom::new();
        raw.to_fun(head);
        let normalized = FunctionBuilder::new(head).finish();
        assert_ne!(raw.as_view().get_data(), normalized.as_view().get_data());
        let raw_key = keys.intern(raw.as_view(), 0);
        let normalized_key = keys.intern(normalized.as_view(), 0);
        assert!(!Arc::ptr_eq(&raw_key, &normalized_key));
        assert!(!raw_key.matches(normalized.as_view().get_data()));
    }

    #[test]
    fn nested_region_keys_store_each_subtree_payload_once() {
        let head = symbol!("regional_key_nested");
        let depth = 128;
        let mut expression = Atom::one();
        for _ in 0..depth {
            expression = function!(head, expression);
        }
        let mut interfaces = AdmittedInterfaces::default();
        expression.visitor(&mut |value| {
            interfaces.insert(value, PartialStructure::from_logical_slots([]));
            true
        });
        let storage = interfaces.storage();
        assert_eq!(storage.payload_bytes, expression.as_view().get_byte_size());
        assert_eq!(storage.nodes, depth + 1);
        assert_eq!(storage.edges, depth);
    }

    #[test]
    fn successive_producer_roots_share_a_large_scalar_spectator() {
        test_support::test_initialize();
        let head = symbol!("regional_key_large_scalar"; Scalar);
        let spectator = Atom::add_many((0..2048).map(|index| function!(head, index)));
        let mut tensor = SymbolicTensor::infer(spectator.clone()).unwrap();
        let initial = tensor.proofs.admitted_interfaces.get().unwrap().storage();
        let repetitions = 48;
        for coefficient in 2..=repetitions + 1 {
            tensor = tensor
                .with_produced_expression(Atom::num(coefficient) * &spectator)
                .unwrap();
        }
        let storage = tensor.proofs.admitted_interfaces.get().unwrap().storage();
        assert!(storage.payload_bytes <= initial.payload_bytes + repetitions * 64);
        assert!(storage.nodes <= initial.nodes + repetitions * 4);
        assert!(storage.edges <= initial.edges + repetitions * 4);
    }

    #[test]
    fn growing_sums_release_obsolete_roots_and_child_edges() {
        test_support::test_initialize();
        let head = tensor_symbol!("regional_key_growing_sum");
        let slot = mink!(4, 71913);
        let mut tensor = SymbolicTensor::infer(function!(head, 0, &slot)).unwrap();
        for term in 1..256 {
            let previous = Arc::downgrade(tensor.proofs.admitted_interfaces.get().unwrap());
            let next = SymbolicTensor::infer(function!(head, term, &slot)).unwrap();
            tensor = tensor.try_add(&next).unwrap();
            assert!(
                previous.upgrade().is_none(),
                "a result must not own predecessor carriers"
            );
        }
        let storage = tensor.proofs.admitted_interfaces.get().unwrap().storage();
        let mut reference = AdmittedInterfaces::default();
        reference.insert(tensor.expression.as_view(), tensor.structure.clone());
        let current = reference.storage();
        assert_eq!(storage.nodes, current.nodes);
        assert_eq!(storage.edges, current.edges);
        assert_eq!(storage.payload_bytes, current.payload_bytes);
        assert_eq!(storage.interner_buckets, current.interner_buckets);
        assert_eq!(storage.interner_records, current.nodes);
        assert_eq!(storage.dead_keys, 0);
        assert!(storage.records <= current.nodes, "{storage:?}");
        assert!(storage.edges <= 1024, "{storage:?}");
        let leaves = tensor.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
        let AtomView::Add(sum) = tensor.expression.as_view() else {
            panic!("distinct terms stay a sum");
        };
        let current_leaf_bytes = sum.iter().map(|term| term.get_byte_size()).sum::<usize>();
        assert_eq!(leaves.len(), sum.get_nargs());
        assert_eq!(
            leaves
                .keys()
                .map(|leaf| leaf.as_view().get_byte_size())
                .sum::<usize>(),
            current_leaf_bytes
        );
        assert!(sum.iter().all(|term| leaves.contains_key(term.get_data())));
        let remaining = sum.iter().next().unwrap().to_owned();
        drop(leaves);
        let previous_root = Arc::downgrade(
            tensor
                .proofs
                .admitted_interfaces
                .get()
                .unwrap()
                .root
                .as_ref()
                .unwrap(),
        );
        let reduced = tensor.with_produced_expression(remaining).unwrap();
        assert_eq!(
            tensor
                .proofs
                .admitted_interfaces
                .get()
                .unwrap()
                .storage()
                .interner_records,
            storage.interner_records,
            "publishing a result must not mutate the source interner",
        );
        drop(tensor);
        assert!(previous_root.upgrade().is_none());
        let storage = reduced.proofs.admitted_interfaces.get().unwrap().storage();
        assert_eq!(storage.interner_records, storage.nodes, "{storage:?}");
        assert!(storage.interner_buckets <= storage.nodes, "{storage:?}");
        assert_eq!(storage.dead_keys, 0, "{storage:?}");
        assert!(storage.nodes < current.nodes, "{storage:?}");
    }

    #[test]
    fn producer_root_order_wins_over_inherited_order() {
        test_support::test_initialize();
        let source = SymbolicTensor::infer(function!(
            tensor_symbol!("regional_ordered_tensor"),
            mink!(4, 71808),
            mink!(4, 71809)
        ))
        .unwrap();
        let mut slots = source.structure.logical_slots();
        slots.reverse();
        let reversed = PartialStructure::from_logical_slots(slots);
        let reordered = source
            .with_produced_parts(source.expression.clone(), reversed.clone(), |slot| slot)
            .unwrap();
        assert_eq!(reordered.structure, reversed);
        assert_eq!(
            reordered
                .proofs
                .retained_interface(reordered.expression.as_view()),
            Some(reversed)
        );
    }

    #[test]
    fn mapped_leaf_caches_keep_only_current_keys() {
        test_support::test_initialize();
        let head = tensor_symbol!("regional_mapped_leaf_cache");
        let x = symbol!("regional_mapped_leaf_x");
        let y = symbol!("regional_mapped_leaf_y");
        let source = SymbolicTensor::infer(function!(head, Atom::var(x), mink!(4, 71810))).unwrap();
        let original = source.proofs.leaf_interfaces.get().unwrap();
        let (expression, mapped) = source
            .map_admitted_expression(
                |value| {
                    Ok(matches!(value, AtomView::Var(var) if var.get_symbol() == x)
                        .then(|| Atom::var(y)))
                },
                |interface| Ok(interface.clone()),
            )
            .unwrap();
        let leaves = mapped.leaf_interfaces.get().unwrap();
        assert!(!Arc::ptr_eq(original, leaves));
        let leaves = leaves.lock().unwrap();
        assert!(!leaves.contains_key(source.expression.as_view().get_data()));
        assert!(leaves.contains_key(expression.as_view().get_data()));
        assert!(
            leaves
                .keys()
                .map(|atom| atom.as_view().get_byte_size())
                .sum::<usize>()
                <= expression.as_view().get_byte_size()
        );
    }

    #[test]
    fn named_constructors_keep_written_ports_without_inference() {
        test_support::test_initialize();
        let head = tensor_symbol!("regional_constructor_written_order");
        let written = [mink!(4, 71916), mink!(4, 71915)];
        let slots = written
            .iter()
            .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()).unwrap())
            .collect::<Vec<_>>();
        let named = Canonicalized::from_iter(slots).map_canonical(|structure| NamedStructure {
            structure,
            global_name: Some(head),
            additional_args: Some(Vec::<Atom>::new()),
        });
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let logical = SymbolicTensor::<OrderedStructure>::from_canonicalized(&named).unwrap();
        let concrete: SymbolicTensor<OrderedStructure> = named
            .canonical()
            .clone()
            .concretize_logical(named.layout())
            .unwrap();
        let canonical = SymbolicTensor::<OrderedStructure>::from_named(named.canonical()).unwrap();
        for value in [&logical, &concrete, &canonical] {
            let AtomView::Fun(function) = value.expression.as_view() else {
                panic!("constructor preserves the tensor head");
            };
            let ports = value.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
            let retained = &ports.get(value.expression.as_view().get_data()).unwrap().1;
            assert_eq!(
                retained
                    .iter()
                    .map(AtomCore::as_atom_view)
                    .collect::<Vec<_>>(),
                function.iter().collect::<Vec<_>>()
            );
            assert!(value.proofs.owns_expression(value.expression.as_view()));
            assert_eq!(retained.len(), 2);
        }
        assert_eq!(
            logical.expression,
            function!(head, &written[0], &written[1])
        );
        assert_eq!(concrete.expression, logical.expression);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn renamed_generic_contraction_retains_regions_and_network_leaves() {
        test_support::test_initialize();
        let a = tensor_symbol!("regional_rename_a");
        let b = tensor_symbol!("regional_rename_b");
        let old = mink!(4, 71811);
        let new = mink!(4, 71812);
        let source = SymbolicTensor::infer(function!(a, &old) * function!(b, &old)).unwrap();
        INFERENCE_CALLS.with(|calls| calls.set(0));
        let renamed = source
            .renamed_indices(&HashMap::from([(
                AbstractIndex::Normal(71811),
                AbstractIndex::Normal(71812),
            )]))
            .unwrap();
        for head in [a, b] {
            let factor = function!(head, &new);
            assert_eq!(renamed.admitted_region(factor.as_view()).unwrap().rank(), 1);
            assert!(
                renamed
                    .proofs
                    .leaf_interfaces
                    .get()
                    .unwrap()
                    .lock()
                    .unwrap()
                    .contains_key(factor.as_view().get_data())
            );
        }
        let network = renamed
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.graph.dangling_indices().is_empty());
        assert!(renamed.contract(Default::default()).unwrap().is_scalar());
        INFERENCE_CALLS.with(|calls| assert_eq!(calls.get(), 0));
        // Looking through the result must not evict the original tensor's facts.
        let original = source.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
        assert!(original.contains_key(function!(a, &old).as_view().get_data()));
        assert!(!original.contains_key(function!(a, &new).as_view().get_data()));
    }

    #[test]
    fn producer_vector_sums_transport_ports_through_generic_scalar_metadata() {
        test_support::test_initialize();
        let q = vector_symbol!("regional_route_q");
        let k = vector_symbol!("regional_route_k");
        let p = vector_symbol!("regional_route_p");
        let port = mink!(4, 71831);
        let tensor = tensor_symbol!("regional_route_generic");
        let owner = symbol!("regional_route_owner"; Scalar);
        let source = SymbolicTensor::infer(
            function!(tensor, function!(owner, function!(q, 0)), &port) * function!(q, 0, &port),
        )
        .unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let routed = source
            .map_produced_expression(&HashMap::new(), |node| {
                let AtomView::Fun(function) = node else {
                    return Ok(None);
                };
                if function.get_symbol() != q {
                    return Ok(None);
                }
                let arguments = function.iter().skip(1).collect::<Vec<_>>();
                Ok::<_, TensorInferenceError>(Some(
                    FunctionBuilder::new(k)
                        .add_arg(0)
                        .add_args(&arguments)
                        .finish()
                        - FunctionBuilder::new(p)
                            .add_arg(0)
                            .add_args(arguments)
                            .finish(),
                ))
            })
            .unwrap();
        let expected = function!(
            tensor,
            function!(owner, function!(k, 0) - function!(p, 0)),
            &port
        ) * (function!(k, 0, &port) - function!(p, 0, &port));
        assert_eq!(routed.expression, expected);
        let cache = routed.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
        for leaf in [
            function!(
                tensor,
                function!(owner, function!(k, 0) - function!(p, 0)),
                &port
            ),
            function!(k, 0, &port),
            function!(p, 0, &port),
        ] {
            assert_eq!(cache.get(&leaf).unwrap().1, vec![port.clone()]);
        }
        drop(cache);
        let network = routed
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.graph.dangling_indices().is_empty());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn symmetric_generic_renaming_reorders_the_carried_raw_port_witness() {
        test_support::test_initialize();
        let head = symbol!("regional_symmetric_rename"; Symmetric; tags = [&SPENSO_TAG.tensor]);
        let first = mink!(4, 71841);
        let second = mink!(4, 71843);
        let source = SymbolicTensor::infer(function!(head, &first, &second)).unwrap();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let renamed = source
            .renamed_indices(&HashMap::from([
                (AbstractIndex::Normal(71841), AbstractIndex::Normal(71847)),
                (AbstractIndex::Normal(71843), AbstractIndex::Normal(71845)),
            ]))
            .unwrap();
        let AtomView::Fun(function) = renamed.expression.as_view() else {
            panic!("symmetric tensor remains a constructor")
        };
        let written = function
            .iter()
            .map(|port| port.to_owned())
            .collect::<Vec<_>>();
        assert_eq!(
            renamed
                .proofs
                .leaf_interfaces
                .get()
                .unwrap()
                .lock()
                .unwrap()
                .get(&renamed.expression)
                .unwrap()
                .1,
            written
        );
        renamed
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn producer_owner_and_index_binding_reuses_generic_leaf_facts() {
        test_support::test_initialize();
        let a = tensor_symbol!("regional_produced_owner_a");
        let b = tensor_symbol!("regional_produced_owner_b");
        let owner = symbol!("regional_produced_owner"; Scalar);
        let old_owner = function!(owner, 0);
        let new_owner = function!(owner, 7);
        let old = mink!(4, 71813);
        let new = mink!(4, 71814);
        let source =
            SymbolicTensor::infer(function!(a, &old_owner, &old) * function!(b, &old)).unwrap();
        let expected = function!(a, &new_owner, &new) * function!(b, &new);
        INFERENCE_CALLS.with(|calls| calls.set(0));
        let result = source
            .map_produced_expression(
                &HashMap::from([(AbstractIndex::Normal(71813), AbstractIndex::Normal(71814))]),
                |node| {
                    Ok::<_, TensorInferenceError>(if node == old.as_view() {
                        Some(new.clone())
                    } else if node == old_owner.as_view() {
                        Some(new_owner.clone())
                    } else {
                        None
                    })
                },
            )
            .unwrap();
        assert_eq!(result.expression(), &expected);
        let leaf = function!(a, &new_owner, &new);
        assert_eq!(result.admitted_region(leaf.as_view()).unwrap().rank(), 1);
        assert!(
            result
                .proofs
                .leaf_interfaces
                .get()
                .unwrap()
                .lock()
                .unwrap()
                .contains_key(leaf.as_view().get_data())
        );
        let network = result
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert!(network.graph.dangling_indices().is_empty());
        INFERENCE_CALLS.with(|calls| assert_eq!(calls.get(), 0));
    }

    #[test]
    fn arithmetic_regions_retain_their_own_boundaries_without_observation() {
        test_support::test_initialize();
        let a = function!(tensor_symbol!("retained_region_a"), mink!(4, 71801));
        let b = function!(tensor_symbol!("retained_region_b"), mink!(4, 71802));
        let left = SymbolicTensor::infer(a.clone()).unwrap();
        let right = SymbolicTensor::infer(b.clone()).unwrap();
        let product = left.outer(&right).unwrap();
        let square = SymbolicTensor::infer(a.pow(2)).unwrap();
        INFERENCE_CALLS.with(|calls| calls.set(0));
        SCOPE_VALIDATIONS.with(|calls| calls.set(0));

        assert_eq!(product.rank(), 2);
        assert_eq!(product.admitted_region(a.as_view()).unwrap().rank(), 1);
        assert_eq!(product.admitted_region(b.as_view()).unwrap().rank(), 1);
        assert!(square.is_scalar());
        assert_eq!(square.admitted_region(a.as_view()).unwrap().rank(), 1);
        let derived = Atom::num(3) * &a * &b;
        assert_eq!(
            product
                .admitted_interface(derived.as_view())
                .unwrap()
                .canonical()
                .order(),
            2
        );
        let unknown = function!(symbol!("retained_region_unknown"), &a);
        assert!(product.admitted_interface(unknown.as_view()).is_none());
        INFERENCE_CALLS.with(|calls| assert_eq!(calls.get(), 0));
        SCOPE_VALIDATIONS.with(|calls| assert_eq!(calls.get(), 0));
    }

    #[test]
    fn shared_region_cache_never_assigns_a_nominal_interface_to_literal_zero() {
        test_support::test_initialize();
        let a = SymbolicTensor::infer(function!(
            tensor_symbol!("retained_region_zero"),
            mink!(4, 71803)
        ))
        .unwrap();
        let zero = a.try_sub(&a).unwrap();
        assert_eq!(
            zero.admitted_interface(Atom::Zero.as_view()),
            Some(a.structure.clone())
        );
        assert!(
            zero.proofs
                .admitted_interface(Atom::Zero.as_view())
                .unwrap()
                .logical_slots()
                .is_empty()
        );
    }

    #[test]
    fn generated_metric_regions_are_derived_once_and_shared() {
        test_support::test_initialize();
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let metric = g!(mink!(4, 71804), mink!(4, 71805));
        REGION_DERIVATIONS.with(|count| count.set(0));
        let first = source.admitted_interface(metric.as_view()).unwrap();
        assert_eq!(first.logical_slots().len(), 2);
        REGION_DERIVATIONS.with(|count| assert_eq!(count.get(), 1));
        let shared = source.clone();
        assert_eq!(shared.admitted_interface(metric.as_view()), Some(first));
        REGION_DERIVATIONS.with(|count| assert_eq!(count.get(), 1));
    }

    #[test]
    fn routed_vector_constructors_preserve_contractible_scalar_regions() {
        test_support::test_initialize();
        let [k, p, q, r, s] = ["k", "p", "q", "r", "s"]
            .map(|name| SPENSO_TAG.rank_one_tensor_symbol(&format!("regional_route::{name}")));
        let port = mink!(4, 71806);
        let source = SymbolicTensor::infer(
            (function!(k, &port) + function!(s, &port)) * function!(p, &port),
        )
        .unwrap();
        let q_component = function!(q, &port);
        let routed = source
            .with_produced_expression(
                (&q_component + function!(r, &port) + function!(s, &port)) * function!(p, &port),
            )
            .unwrap();
        assert_eq!(
            routed
                .admitted_interface(q_component.as_view())
                .unwrap()
                .logical_slots()
                .len(),
            1
        );
        let contracted = routed.contract(Default::default()).unwrap();
        assert!(contracted.is_scalar());
        let expected = Atom::add_many(
            [q, r, s].map(|head| g!(function!(head, mink!(4)), function!(p, mink!(4)))),
        );
        assert_eq!(contracted.expression(), &expected);
        assert_eq!(contracted.contract(Default::default()).unwrap(), contracted);
    }

    #[test]
    fn mapped_regions_follow_the_single_payload_normalization() {
        test_support::test_initialize();
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let head = tensor_symbol!(
            "regional_map_callback",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let x = symbol!("regional_map_x");
        let y = symbol!("regional_map_y");
        let source =
            SymbolicTensor::infer(Atom::num(2) * function!(head, Atom::var(x), mink!(4, 71807)))
                .unwrap();
        calls.store(0, Ordering::Relaxed);
        let expected = source.expression.replace_map(|node, _, output| {
            if matches!(node, AtomView::Var(variable) if variable.get_symbol() == x) {
                **output = Atom::var(y);
            }
        });
        let expected_calls = calls.load(Ordering::Relaxed);
        assert_eq!(expected_calls, 1);
        calls.store(0, Ordering::Relaxed);
        let (expression, proofs) = source
            .map_admitted_expression(
                |node| {
                    Ok(
                        matches!(node, AtomView::Var(variable) if variable.get_symbol() == x)
                            .then(|| Atom::var(y)),
                    )
                },
                |interface| Ok(interface.clone()),
            )
            .unwrap();
        assert_eq!(expression, expected);
        assert_eq!(calls.load(Ordering::Relaxed), expected_calls);
        assert_eq!(
            proofs.admitted_interface(expression.as_view()),
            Some(source.structure.clone())
        );
    }
}
