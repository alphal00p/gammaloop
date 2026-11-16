//! Interfaces retained at admission and composed by known arithmetic producers.

use std::borrow::Borrow;
#[cfg(test)]
use std::cell::Cell;
#[cfg(test)]
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::mem;
use std::ops::Range;
use std::sync::{Arc, Mutex, OnceLock};

use ahash::{AHashMap, AHashSet};
use spenso::network::library::symbolic::ETS;
use spenso::network::parsing::LeafInterfaceCache;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::TensorStructure;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructure, PartialStructureExt};
use spenso::structure::representation::LibraryRep;
use spenso::structure::slot::{IsAbstractSlot, SlotMatcher};
use symbolica::atom::representation::FunView;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol};

use super::composition;
use super::inference::{InterfaceInference, TensorInferenceError};
use super::raw_leaves::RawLeafInterfaces;
use super::simplification::observation::DomainObservations;
use super::{SymbolicTensor, TensorProofs};
use crate::color::simplify::ColorAlgebraSimplifier;

#[cfg(test)]
thread_local! {
    pub(crate) static REGION_DERIVATIONS: Cell<usize> = const { Cell::new(0) };
    static INDEXED_KEYS: Cell<usize> = const { Cell::new(0) };
    static REGION_PROBES: Cell<usize> = const { Cell::new(0) };
}

/// A fact refers to the exact serialized bytes of a subexpression. All ordinary
/// regions of a published tensor share one allocation for the current payload.
#[derive(Clone, Debug)]
pub(super) struct RegionKey {
    payload: Arc<Atom>,
    range: Range<usize>,
}

impl RegionKey {
    pub(super) fn owned(value: AtomView<'_>) -> Self {
        Self {
            payload: Arc::new(value.to_owned()),
            range: 0..value.get_data().len(),
        }
    }

    /// Compact encodings can synthesize a child outside the payload.
    pub(super) fn within(payload: &Arc<Atom>, value: AtomView<'_>) -> Self {
        let storage = payload.as_view().get_data();
        let bytes = value.get_data();
        if let Some(start) = (bytes.as_ptr() as usize).checked_sub(storage.as_ptr() as usize)
            && let Some(end) = start.checked_add(bytes.len())
            && end <= storage.len()
        {
            Self {
                payload: Arc::clone(payload),
                range: start..end,
            }
        } else {
            Self::owned(value)
        }
    }

    pub(super) fn child(&self, range: &Range<usize>) -> Self {
        debug_assert!(range.end <= self.range.len());
        Self {
            payload: Arc::clone(&self.payload),
            range: self.range.start + range.start..self.range.start + range.end,
        }
    }

    #[cfg(test)]
    pub(super) fn owner(&self) -> &Arc<Atom> {
        &self.payload
    }

    pub(super) fn bytes(&self) -> &[u8] {
        &self.payload.as_view().get_data()[self.range.clone()]
    }
}

impl Borrow<[u8]> for RegionKey {
    fn borrow(&self) -> &[u8] {
        self.bytes()
    }
}

impl Hash for RegionKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bytes().hash(state);
    }
}

impl PartialEq for RegionKey {
    fn eq(&self, other: &Self) -> bool {
        self.bytes() == other.bytes()
    }
}

impl Eq for RegionKey {}

type RegionInterfaces = AHashMap<RegionKey, Arc<PartialStructure>>;

#[derive(Clone, Debug)]
struct RegionPayload {
    atom: Arc<Atom>,
    // During admission, observations borrow the input. This address is only
    // an offset hint into our owned copy. No pointer is dereferenced, and the
    // complete bytes are checked before a range is used.
    source_start: usize,
}

impl RegionPayload {
    fn key(&self, value: AtomView<'_>) -> RegionKey {
        let bytes = value.get_data();
        let payload = self.atom.as_view().get_data();
        for start in [payload.as_ptr() as usize, self.source_start] {
            if let Some(offset) = (bytes.as_ptr() as usize).checked_sub(start)
                && let Some(end) = offset.checked_add(bytes.len())
                && payload.get(offset..end) == Some(bytes)
            {
                return RegionKey {
                    payload: Arc::clone(&self.atom),
                    range: offset..end,
                };
            }
        }
        // Compact powers can synthesize child views outside the root bytes.
        // Explicitly retained off-root derivations also need standalone keys.
        RegionKey::owned(value)
    }
}

/// Projected subexpressions borrow the output's one immutable payload. A scalar
/// or a repeated factor therefore needs only a byte range and its shared boundary.
#[derive(Debug)]
struct ProjectedRegion {
    start: u32,
    length: u32,
    boundary: Arc<PartialStructure>,
}

impl ProjectedRegion {
    fn bytes<'a>(&self, payload: &'a [u8]) -> &'a [u8] {
        let start = self.start as usize;
        &payload[start..start + self.length as usize]
    }
}

/// The root and its reachable facts belong to this payload. The independent
/// values of unchanged interfaces can be shared without retaining an ancestor.
#[derive(Debug, Default)]
pub(crate) struct AdmittedInterfaces {
    root: Option<RegionPayload>,
    interface: Option<Arc<PartialStructure>>,
    observed: RegionInterfaces,
    projected: Box<[ProjectedRegion]>,
    derived: OnceLock<Mutex<RegionInterfaces>>,
}

#[derive(Clone, Copy)]
enum MappingOrder {
    BeforeChildren,
    VectorComponents(Symbol),
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
            MappingOrder::VectorComponents(head) if value.get_symbol() == Some(head) => {
                (self.leaf_map)(value)?
            }
            MappingOrder::VectorComponents(_) | MappingOrder::AfterChildren => None,
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
                    .map(|slot| match self.order {
                        // The momentum producer changes only the component owner
                        // and coefficient. Its original written port stays fixed.
                        MappingOrder::VectorComponents(_) => Some(slot.clone()),
                        _ => self.mapped_slots.get(slot).cloned(),
                    })
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
                if let MappingOrder::VectorComponents(head) = self.order
                    && source.get_symbol() == crate::dirac::AGS.gamma
                    && let Some(arguments) = &mapped_arguments
                {
                    self.carry_routed_gamma(head, source, arguments, &mapped, result.as_view());
                }
                let mut matcher = SlotMatcher::default();
                let component = match mapped.as_slice() {
                    [port] => Some(port.as_view()),
                    [] if matches!(self.order, MappingOrder::VectorComponents(head) if head == source.get_symbol()) =>
                    {
                        // A fixed component has no raw graph port. Carry its
                        // scalar declaration only to that same component.
                        matcher.vector_argument(source).filter(|&argument| {
                            matcher
                                .concrete_component(argument)
                                .is_some_and(|index| index.is_ok())
                        })
                    }
                    _ => None,
                };
                if source.get_symbol().has_tag(&SPENSO_TAG.rank1)
                    && let Some(component) = component
                    && let Some(interface) = self.source.admitted_interface(value)
                {
                    let interface = (self.interface_map)(&interface)?;
                    let mut pending = vec![result.as_view()];
                    while let Some(node) = pending.pop() {
                        match node {
                            AtomView::Add(sum) => pending.extend(sum.iter()),
                            AtomView::Mul(product) => pending.extend(product.iter()),
                            AtomView::Pow(power) => pending.push(power.get_base_exp().0),
                            AtomView::Fun(function)
                                if function.get_symbol().has_tag(&SPENSO_TAG.rank1)
                                    && matcher.vector_argument(function) == Some(component) =>
                            {
                                // Constructor semantics verify the producer's
                                // port or fixed coordinate on every emitted vector.
                                self.interfaces.insert(node, interface.clone());
                                self.leaves
                                    .lock()
                                    .expect("leaf interface cache is not poisoned")
                                    .insert(node.to_owned(), (tensorial, mapped.clone()));
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    /// Linearity pulls sums and scalar weights out of a slashed momentum.
    /// Check each emitted gamma against the exact replacement vector and both
    /// supplied spinor operands before carrying its written incidences.
    fn carry_routed_gamma(
        &mut self,
        vector_head: Symbol,
        source: FunView<'_>,
        arguments: &[AtomOrView<'_>],
        ports: &[Atom],
        result: AtomView<'_>,
    ) {
        fn visit_factors(value: AtomView<'_>, visit: &mut impl FnMut(FunView<'_>)) {
            match value {
                AtomView::Add(_) | AtomView::Mul(_) => {
                    for child in value.children() {
                        visit_factors(child, visit);
                    }
                }
                AtomView::Fun(function) => visit(function),
                _ => {}
            }
        }

        let mut matcher = SlotMatcher::default();
        let Some((position, old_port)) = source.iter().enumerate().find_map(|(index, argument)| {
            let AtomView::Fun(vector) = argument else {
                return None;
            };
            let port = matcher.vector_argument(vector)?;
            (vector.get_symbol() == vector_head
                && matcher.parse_representation::<LibraryRep>(port).is_ok())
            .then_some((index, port))
        }) else {
            return;
        };
        let Some(interface) = self.source.admitted_interface(source.as_view()) else {
            return;
        };
        let mut vectors = AHashSet::new();
        visit_factors(arguments[position].as_view(), &mut |vector| {
            if vector.get_symbol() == vector_head
                && matcher.vector_argument(vector) == Some(old_port)
            {
                vectors.insert(vector.as_view().to_owned());
            }
        });
        let mut leaves = self
            .leaves
            .lock()
            .expect("leaf interface cache is not poisoned");
        visit_factors(result, &mut |gamma| {
            let node = gamma.as_view();
            if gamma.get_symbol() == source.get_symbol()
                && let Some(vector) = gamma.iter().find(|arg| vectors.contains(arg.get_data()))
            {
                let mut supplied = arguments
                    .iter()
                    .map(AtomCore::as_atom_view)
                    .collect::<Vec<_>>();
                supplied[position] = vector;
                if let Some((_, written)) =
                    TensorProofs::raw_ports_after_argument_permutation(&supplied, node, ports)
                {
                    self.interfaces.insert(node, interface.clone());
                    leaves.insert(node.to_owned(), (true, written));
                }
            }
        });
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

/// Start of a subexpression whose complete bytes lie in the shared output payload.
pub(super) fn payload_offset(payload: &[u8], child: &[u8]) -> Option<u32> {
    let offset = (child.as_ptr() as usize).checked_sub(payload.as_ptr() as usize)?;
    checked_payload_offset(payload.len(), offset, child.len())
}

// A compact record needs both a complete subview and a representable start.
fn checked_payload_offset(payload_length: usize, offset: usize, length: usize) -> Option<u32> {
    if offset.checked_add(length)? > payload_length {
        return None;
    }
    u32::try_from(offset).ok()
}

/// Small selections from an enclosing admission probe the source directly.
/// Grow an output index once there are enough distinct nodes to use it.
pub(super) fn output_positions<'a>(
    root: AtomView<'a>,
    try_direct: bool,
) -> (Vec<AtomView<'a>>, Option<AHashMap<&'a [u8], usize>>) {
    const DIRECT_LIMIT: usize = 8;
    let mut nodes = vec![root];
    let mut positions = (!try_direct).then(|| AHashMap::from_iter([(root.get_data(), 0)]));
    let mut cursor = 0;
    while cursor < nodes.len() {
        let node = nodes[cursor];
        let mut add = |child: AtomView<'a>| {
            let bytes = child.get_data();
            if let Some(index) = &mut positions {
                index.entry(bytes).or_insert_with(|| {
                    nodes.push(child);
                    nodes.len() - 1
                });
            } else if !nodes.iter().any(|previous| previous.get_data() == bytes) {
                nodes.push(child);
                if nodes.len() > DIRECT_LIMIT {
                    positions = Some(
                        nodes
                            .iter()
                            .enumerate()
                            .map(|(i, node)| (node.get_data(), i))
                            .collect(),
                    );
                }
            }
        };
        match node {
            AtomView::Fun(function) => function.iter().for_each(add),
            AtomView::Mul(product) => product.iter().for_each(add),
            AtomView::Add(sum) => sum.iter().for_each(add),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                add(base);
                add(exponent);
            }
            AtomView::Num(_) | AtomView::Var(_) => {}
        }
        cursor += 1;
    }
    #[cfg(test)]
    if positions.is_some() {
        INDEXED_KEYS.with(|count| count.set(count.get() + nodes.len()));
    }
    (nodes, positions)
}

impl AdmittedInterfaces {
    pub(crate) fn for_expression(value: AtomView<'_>) -> Self {
        Self {
            root: (!value.is_zero()).then(|| RegionPayload {
                atom: Arc::new(value.to_owned()),
                source_start: value.get_data().as_ptr() as usize,
            }),
            ..Self::default()
        }
    }

    pub(crate) fn from_region(value: AtomView<'_>, interface: PartialStructure) -> Self {
        let mut result = Self::for_expression(value);
        result.insert(value, interface);
        result
    }

    fn owns(&self, bytes: &[u8]) -> bool {
        self.root
            .as_ref()
            .is_some_and(|root| root.atom.as_view().get_data() == bytes)
    }

    fn observed(&self, bytes: &[u8]) -> Option<&Arc<PartialStructure>> {
        if self.owns(bytes) {
            self.interface.as_ref()
        } else {
            self.observed
                .get(bytes)
                .or_else(|| self.projected_interface(bytes))
        }
    }

    fn projected_interface(&self, bytes: &[u8]) -> Option<&Arc<PartialStructure>> {
        if self.projected.is_empty() {
            return None;
        }
        let payload = self.root.as_ref()?.atom.as_view().get_data();
        let position = self
            .projected
            .binary_search_by(|region| {
                (region.length as usize)
                    .cmp(&bytes.len())
                    .then_with(|| region.bytes(payload).cmp(bytes))
            })
            .ok()?;
        Some(&self.projected[position].boundary)
    }

    fn key(&self, value: AtomView<'_>) -> RegionKey {
        self.root
            .as_ref()
            .map_or_else(|| RegionKey::owned(value), |root| root.key(value))
    }

    pub(crate) fn insert(&mut self, value: AtomView<'_>, interface: PartialStructure) {
        if value.is_zero() {
            return;
        }
        if self.owns(value.get_data()) {
            self.interface = Some(Arc::new(interface));
        } else {
            let key = self.key(value);
            self.observed.insert(key, Arc::new(interface));
        }
    }

    pub(crate) fn declare_root(
        &mut self,
        value: AtomView<'_>,
        interface: PartialStructure,
        leaves: LeafInterfaceCache,
    ) -> RawLeafInterfaces {
        if self.owns(value.get_data()) {
            self.insert(value, interface);
            return leaves.into();
        }
        // Product normalization can erase an admitted bracket. Keep the
        // surviving factors' raw ports and release the erased bracket's entry.
        let previous = mem::replace(self, Self::from_region(value, interface));
        self.project(&[Arc::new(previous)], &[leaves.into()])
    }

    fn get(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        let bytes = value.get_data();
        if let Some(interface) = self.observed(bytes) {
            return Some(interface.as_ref().clone());
        }
        self.derived
            .get()?
            .lock()
            .expect("derived interface cache is not poisoned")
            .get(bytes)
            .map(|interface| interface.as_ref().clone())
    }

    fn retain(&self, value: AtomView<'_>, interface: &PartialStructure) {
        if !value.is_zero() {
            self.derived
                .get_or_init(Default::default)
                .lock()
                .expect("derived interface cache is not poisoned")
                .insert(self.key(value), Arc::new(interface.clone()));
        }
    }

    /// Join by exact output bytes. Source precedence applies independently to
    /// region and raw-port facts, and mutexes are acquired one source at a time.
    fn project(
        &mut self,
        regions: &[Arc<Self>],
        leaves: &[RawLeafInterfaces],
    ) -> RawLeafInterfaces {
        let Some(root) = &self.root else {
            return LeafInterfaceCache::default().into();
        };
        if regions.is_empty() && leaves.is_empty() {
            return LeafInterfaceCache::default().into();
        }
        let payload = Arc::clone(&root.atom);
        let (nodes, mut positions) =
            output_positions(payload.as_view(), regions.len() <= 2 && leaves.len() <= 1);
        let retained = RawLeafInterfaces::project(&payload, &nodes, positions.as_ref(), leaves);
        let mut projected = Vec::new();
        let output = payload.as_view().get_data();
        let mut direct = Vec::new();
        if let Some(index) = &mut positions {
            if self.interface.is_some() {
                index.remove(output);
            }
        } else {
            direct.extend(usize::from(self.interface.is_some())..nodes.len());
        }
        let mut keep = |position: usize, boundary: &Arc<PartialStructure>| {
            if position == 0 {
                self.interface = Some(Arc::clone(boundary));
                return;
            }
            let value = nodes[position];
            let bytes = value.get_data();
            if let Some(start) = payload_offset(output, bytes)
                && let Ok(length) = u32::try_from(bytes.len())
            {
                projected.push(ProjectedRegion {
                    start,
                    length,
                    boundary: Arc::clone(boundary),
                });
            } else {
                // A compact power can synthesize its exponent outside the payload.
                // Large atoms retain the same unrestricted exact-key representation.
                let key = self.key(value);
                self.observed
                    .entry(key)
                    .or_insert_with(|| Arc::clone(boundary));
            }
        };
        for source in regions {
            if positions
                .as_ref()
                .map_or(direct.is_empty(), |index| index.is_empty())
            {
                break;
            }
            let derived = source.derived.get().map(|derived| {
                derived
                    .lock()
                    .expect("derived interface cache is not poisoned")
            });
            let mut probe = |bytes: &[u8], position: usize| {
                #[cfg(test)]
                REGION_PROBES.with(|probes| probes.set(probes.get() + 1));
                if let Some(boundary) = source
                    .observed(bytes)
                    .or_else(|| derived.as_ref()?.get(bytes))
                {
                    keep(position, boundary);
                    false
                } else {
                    true
                }
            };
            let Some(positions) = &mut positions else {
                direct.retain(|&position| probe(nodes[position].get_data(), position));
                continue;
            };
            let count = usize::from(source.interface.is_some())
                + source.observed.len()
                + source.projected.len()
                + derived.as_ref().map_or(0, |entries| entries.len());
            if count < positions.len() {
                // A product often inherits many small factors. Read each factor's
                // records once, instead of probing it for every output subexpression.
                let mut take = |bytes: &[u8], boundary: &Arc<PartialStructure>| {
                    #[cfg(test)]
                    REGION_PROBES.with(|probes| probes.set(probes.get() + 1));
                    if let Some(position) = positions.remove(bytes) {
                        keep(position, boundary);
                    }
                };
                if let Some(root) = &source.root
                    && let Some(boundary) = &source.interface
                {
                    take(root.atom.as_view().get_data(), boundary);
                }
                for (key, boundary) in &source.observed {
                    take(key.bytes(), boundary);
                }
                if let Some(root) = &source.root {
                    let bytes = root.atom.as_view().get_data();
                    for region in &source.projected {
                        take(region.bytes(bytes), &region.boundary);
                    }
                }
                if let Some(entries) = &derived {
                    for (key, boundary) in entries.iter() {
                        take(key.bytes(), boundary);
                    }
                }
            } else {
                // Selecting a leaf from a large input has the opposite cost.
                positions.retain(|bytes, position| probe(bytes, *position));
            }
        }
        projected.sort_unstable_by(|a, b| {
            a.length
                .cmp(&b.length)
                .then_with(|| a.bytes(output).cmp(b.bytes(output)))
        });
        self.projected = projected.into_boxed_slice();
        retained
    }
}

impl TensorProofs {
    /// Reuse an admitted payload for the simplifier's scan and its scope keys.
    pub(super) fn observation_payload(&self, value: AtomView<'_>) -> Arc<Atom> {
        self.admitted_interfaces
            .get()
            .and_then(|interfaces| interfaces.root.as_ref())
            .filter(|root| root.atom.as_view() == value)
            .map_or_else(|| Arc::new(value.to_owned()), |root| Arc::clone(&root.atom))
    }

    /// Whether these facts belong to this exact current payload, including its
    /// normalization flags. This does not infer or derive any interface.
    pub(crate) fn owns_expression(&self, value: AtomView<'_>) -> bool {
        self.admitted_interfaces
            .get()
            .is_some_and(|interfaces| interfaces.owns(value.get_data()))
    }

    /// Publish a new authoritative root and retain only its reachable regional
    /// and leaf facts. Input carriers remain independently usable.
    pub(crate) fn publish_interfaces<'a>(
        &mut self,
        value: AtomView<'_>,
        interface: impl Into<Option<PartialStructure>>,
        sources: impl IntoIterator<Item = &'a Self>,
    ) {
        self.project_interfaces(Some((value, interface.into())), sources);
    }

    /// An operation on typed axes may change the declaration while preserving
    /// the exact expression and its written parser ports.
    pub(super) fn update_root_interface(
        &mut self,
        value: AtomView<'_>,
        interface: &PartialStructure,
    ) {
        if value.is_zero() || self.owns_published_interface(value, interface) {
            return;
        }
        if let Some(owner) = self.admitted_interfaces.get_mut().and_then(Arc::get_mut)
            && owner.owns(value.get_data())
        {
            owner.interface = Some(Arc::new(interface.clone()));
        } else {
            self.publish_interfaces(value, interface.clone(), []);
        }
    }

    /// Replace only the root of a private projection. Callback validation may
    /// derive child boundaries from its raw leaves after the initial projection.
    /// Promote those reachable facts without walking the expression a second time.
    pub(super) fn set_projected_boundary(
        &mut self,
        value: AtomView<'_>,
        boundary: PartialStructure,
    ) {
        let owner = self
            .admitted_interfaces
            .get_mut()
            .and_then(Arc::get_mut)
            .expect("a freshly projected interface owner is private");
        debug_assert!(owner.owns(value.get_data()));
        if let Some(derived) = owner.derived.take() {
            let mut derived = derived
                .into_inner()
                .expect("derived interface cache is not poisoned");
            derived.remove(value.get_data());
            for (region, interface) in derived {
                if owner.projected_interface(region.bytes()).is_none() {
                    owner.observed.entry(region).or_insert(interface);
                }
            }
        }
        owner.interface = Some(Arc::new(boundary));
    }

    /// Fetch only previously supplied facts, without deriving any new region.
    pub(crate) fn retained_interface(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        self.admitted_interfaces.get()?.get(value)
    }

    /// A planning leaf may refer back to an immutable published boundary.
    /// Derived facts can be replaced while the plan remains alive.
    pub(super) fn published_interface_matches(
        &self,
        value: AtomView<'_>,
        interface: &PartialStructure,
    ) -> bool {
        self.published_interface(value) == Some(interface)
    }

    pub(super) fn published_interface(&self, value: AtomView<'_>) -> Option<&PartialStructure> {
        self.admitted_interfaces
            .get()?
            .observed(value.get_data())
            .map(Arc::as_ref)
    }

    /// Record the interface supplied by a leaf producer. Publication retains
    /// this fact only while its exact expression belongs to the output tree.
    pub(crate) fn retain_interface(&self, value: AtomView<'_>, interface: PartialStructure) {
        self.record_interface(value, &interface, true);
    }

    /// An enclosing sum or parser storage declaration only supplies the free
    /// boundary. Preserve a composite's independently established written ports.
    pub(super) fn supply_interface(&self, value: AtomView<'_>, interface: &PartialStructure) {
        self.record_interface(value, interface, false);
    }

    fn record_interface(
        &self,
        value: AtomView<'_>,
        interface: &PartialStructure,
        replace_ports: bool,
    ) {
        if matches!(value, AtomView::Fun(function)
            if composition::is_composite_head(function.get_symbol(), &SPENSO_TAG)
                || function.get_symbol() == SPENSO_TAG.dot)
        {
            let mut leaves = self
                .leaf_interfaces
                .get_or_init(Default::default)
                .lock()
                .expect("leaf interface cache is not poisoned");
            if replace_ports
                || !leaves
                    .get(value.get_data())
                    .is_some_and(|(positive, _)| *positive)
            {
                let ports = interface
                    .logical_slots()
                    .into_iter()
                    .map(composition::port_atom)
                    .collect();
                leaves.insert(value.to_owned(), (true, ports));
            }
        }
        self.retain_boundary(value, interface);
    }

    pub(super) fn retain_boundary(&self, value: AtomView<'_>, interface: &PartialStructure) {
        self.admitted_interfaces
            .get_or_init(|| Arc::new(AdmittedInterfaces::default()))
            .retain(value, interface);
    }

    /// Compose a newly constructed metric from its admitted endpoint roles.
    /// Explicit slots survive; representations create positional axes; vector
    /// operands are consumed inside the metric. This never examines its output.
    pub(crate) fn metric_interface(
        &self,
        arguments: [AtomView<'_>; 2],
    ) -> Result<PartialStructure, TensorInferenceError> {
        self.metric_interface_in::<true>(arguments, None)
    }

    fn metric_interface_in<const RECORD: bool>(
        &self,
        arguments: [AtomView<'_>; 2],
        observations: Option<&DomainObservations>,
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
                let interface =
                    if let Some(interface) = self.interface_in::<RECORD>(argument, observations) {
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
                        if RECORD {
                            self.retain_leaf_ports(argument, vec![port.to_owned()])?;
                        }
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
        self.project_interfaces(None, sources);
    }

    fn project_interfaces<'a>(
        &mut self,
        output: Option<(AtomView<'_>, Option<PartialStructure>)>,
        sources: impl IntoIterator<Item = &'a Self>,
    ) {
        let mut leaves = Vec::new();
        if self.leaf_interfaces.is_present() {
            leaves.push(self.leaf_interfaces.clone());
        }
        let mut predecessors = self
            .admitted_interfaces
            .get()
            .cloned()
            .into_iter()
            .collect::<Vec<_>>();
        for source in sources {
            if source.leaf_interfaces.is_present()
                && !leaves
                    .iter()
                    .any(|existing| existing.same_source(&source.leaf_interfaces))
            {
                leaves.push(source.leaf_interfaces.clone());
            }
            if let Some(source) = source.admitted_interfaces.get()
                && !predecessors
                    .iter()
                    .any(|existing| Arc::ptr_eq(existing, source))
            {
                predecessors.push(Arc::clone(source));
            }
        }
        let mut result = if let Some((value, interface)) = output {
            // Copies of the same expression share its immutable bytes. Keep
            // writable derived and raw tables private to the new producer.
            let owner = predecessors
                .iter()
                .find(|source| source.owns(value.get_data()));
            let mut result = owner.map_or_else(
                || AdmittedInterfaces::for_expression(value),
                |source| AdmittedInterfaces {
                    root: source.root.as_ref().map(|root| RegionPayload {
                        atom: Arc::clone(&root.atom),
                        source_start: value.get_data().as_ptr() as usize,
                    }),
                    ..Default::default()
                },
            );
            if let Some(interface) = interface {
                result.insert(value, interface);
            }
            result
        } else if let Some(first) = predecessors.first() {
            let mut result = AdmittedInterfaces {
                root: first.root.clone(),
                ..Default::default()
            };
            if predecessors.len() == 1 {
                self.leaf_interfaces = result.project(&[], &leaves);
                self.admitted_interfaces = OnceLock::from(Arc::clone(first));
                return;
            }
            result
        } else {
            return;
        };
        let retained = result.project(&predecessors, &leaves);
        self.admitted_interfaces = OnceLock::from(Arc::new(result));
        self.leaf_interfaces = retained;
    }

    /// Look up an admitted region, or compose arithmetic from admitted children.
    /// Unknown function syntax is never inspected to rediscover its ports.
    pub(crate) fn admitted_interface(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        self.interface_in::<true>(value, None)
    }

    /// Inspect a candidate rewrite against this admission. Candidate expressions
    /// can be discarded, so their constructor and region facts stay local.
    pub(crate) fn inspect_interface(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        self.interface_in::<false>(value, None)
    }

    fn interface_in<const RECORD: bool>(
        &self,
        value: AtomView<'_>,
        observations: Option<&DomainObservations>,
    ) -> Option<PartialStructure> {
        if let Some(interface) = self
            .admitted_interfaces
            .get()
            .and_then(|interfaces| interfaces.get(value))
        {
            return Some(interface);
        }
        #[cfg(test)]
        REGION_DERIVATIONS.with(|count| count.set(count.get() + 1));
        let leaf = self.leaf_interfaces.with_fact(value, |tensorial, ports| {
            if tensorial {
                Self::boundary_from_leaf_ports(ports)
            } else {
                Ok(PartialStructure::from_logical_slots([]))
            }
        });
        if let Some(boundary) = leaf {
            // Raw ports remain in written order, including internal pairs.
            // A read may memoize their free boundary without changing those ports.
            let boundary = boundary.ok()?;
            if RECORD {
                self.retain_boundary(value, &boundary);
            }
            return Some(boundary);
        }
        let child = |child| self.interface_in::<RECORD>(child, observations);
        let interface = if observations.is_some_and(|observed| {
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
                    let interface = child(terms.next()?)?;
                    terms
                        .all(|term| {
                            child(term).is_some_and(|other| {
                                InterfaceInference::additive_interfaces_match(&interface, &other)
                            })
                        })
                        .then_some(interface)
                }
                AtomView::Mul(product) => {
                    let interfaces = product.iter().map(child).collect::<Option<Vec<_>>>()?;
                    InterfaceInference::merge_explicit_interface_sequence(&interfaces).ok()
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let interface = child(base)?;
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
                        factors.push(
                            self.metric_interface_in::<RECORD>(
                                [header, args.next()?],
                                observations,
                            )
                            .ok()?,
                        );
                    }
                    for factor in args {
                        factors.push(child(factor)?);
                    }
                    InterfaceInference::merge_explicit_interface_sequence(&factors).ok()
                }
                // Brackets and the symmetry projectors carry their operands'
                // combined boundary. Unary broadcast heads preserve it.
                AtomView::Fun(function)
                    if composition::is_composite_head(function.get_symbol(), &SPENSO_TAG) =>
                {
                    let operands = function.iter().map(child).collect::<Option<Vec<_>>>()?;
                    InterfaceInference::merge_explicit_interface_sequence(&operands).ok()
                }
                // Analytic scalar coefficients emitted by integration, such as
                // log(mu^2), retain the known scalar boundary of their argument.
                AtomView::Fun(function)
                    if function.get_nargs() == 1
                        && [
                            Symbol::LOG,
                            Symbol::EXP,
                            Symbol::SIN,
                            Symbol::COS,
                            Symbol::SQRT,
                        ]
                        .contains(&function.get_symbol()) =>
                {
                    child(function.get(0)).filter(|value| value.canonical().is_scalar())
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
                    self.metric_interface_in::<RECORD>(
                        [arguments.next()?, arguments.next()?],
                        observations,
                    )
                    .ok()
                }
                AtomView::Fun(function) => {
                    if let Some(ports) = InterfaceInference::constructor_ports(
                        function,
                        &mut SlotMatcher::default(),
                        false,
                    ) {
                        if RECORD {
                            self.retain_leaf_ports(value, ports).ok()?;
                            self.retained_interface(value)
                        } else {
                            Self::boundary_from_leaf_ports(ports.iter().map(Atom::as_view)).ok()
                        }
                    } else {
                        ColorAlgebraSimplifier::produced_interface(value)
                    }
                }
            }
        }?;
        if RECORD && let Some(cache) = self.admitted_interfaces.get() {
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

    pub(crate) fn map_vector_nodes(
        &self,
        head: Symbol,
        mut route: impl FnMut(AtomView<'_>) -> Option<Atom>,
    ) -> Result<(Atom, TensorProofs), TensorInferenceError> {
        self.map_expression(
            self.expression.as_view(),
            |node| Ok(route(node)),
            |interface| Ok(interface.clone()),
            MappingOrder::VectorComponents(head),
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
        let source_leaves = self.proofs.leaf_interfaces.read_cache().unwrap_or_default();
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
            interfaces: AdmittedInterfaces::default(),
            source_leaves,
            leaves: LeafInterfaceCache::default(),
            written_slots,
            mapped_slots: AHashMap::new(),
        };
        let expression = mapping.map(expression)?.into_owned();
        // Keep the declared root order even when normalization returns an
        // already observed child, then retain only still-reachable facts.
        let mut interfaces = AdmittedInterfaces::from_region(
            expression.as_view(),
            (mapping.interface_map)(&self.structure)?,
        );
        let sources = [Some(Arc::new(mapping.interfaces)), previous.cloned()]
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
        let leaves = interfaces.project(
            &sources,
            &[mapping.leaves.into(), mapping.source_leaves.into()],
        );
        let interfaces = Arc::new(interfaces);
        Ok((
            expression,
            TensorProofs {
                admitted_interfaces: OnceLock::from(interfaces),
                leaf_interfaces: leaves,
                ..TensorProofs::default()
            },
        ))
    }

    pub(crate) fn admitted_interface(&self, value: AtomView<'_>) -> Option<PartialStructure> {
        if value == self.expression.as_view() {
            Some(self.structure.clone())
        } else {
            self.proofs
                .interface_in::<true>(value, self.work.observations().map(AsRef::as_ref))
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
    use crate::tensor::regions;
    use crate::test_support;
    use spenso::network::parsing::ParseSettings;
    use spenso::shadowing::Concretize;
    use spenso::structure::representation::LibrarySlot;
    use spenso::structure::{Canonicalized, NamedStructure, OrderedStructure, TensorStructure};
    use spenso::{g, mink, tensor_symbol, vector_symbol};
    use std::hash::BuildHasherDefault;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use symbolica::atom::FunctionBuilder;
    use symbolica::{function, symbol};

    #[derive(Debug)]
    struct KeyStorage {
        payloads: usize,
        payload_bytes: usize,
        records: usize,
    }

    impl AdmittedInterfaces {
        fn storage(&self) -> KeyStorage {
            let derived = self.derived.get().map(|derived| derived.lock().unwrap());
            let mut payloads = AHashMap::new();
            if let Some(root) = &self.root {
                payloads.insert(Arc::as_ptr(&root.atom), root.atom.as_view().get_byte_size());
            }
            for key in self
                .observed
                .keys()
                .chain(derived.iter().flat_map(|map| map.keys()))
            {
                payloads.insert(
                    Arc::as_ptr(&key.payload),
                    key.payload.as_view().get_byte_size(),
                );
            }
            KeyStorage {
                payloads: payloads.len(),
                payload_bytes: payloads.values().sum(),
                records: self.observed.len()
                    + self.projected.len()
                    + usize::from(self.interface.is_some())
                    + derived.as_ref().map_or(0, |map| map.len()),
            }
        }
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn compact_payload_offsets_leave_large_or_noncontained_ranges_to_owned_keys() {
        let last_offset = u32::MAX as usize;
        assert_eq!(
            regions::checked_payload_offset(last_offset + 2, last_offset, 2),
            Some(u32::MAX)
        );
        assert_eq!(
            regions::checked_payload_offset(last_offset + 2, last_offset + 1, 1),
            None
        );
        assert_eq!(regions::checked_payload_offset(5, 4, 2), None);
        assert_eq!(
            regions::checked_payload_offset(usize::MAX, usize::MAX, 1),
            None
        );
    }

    #[test]
    fn exact_key_collisions_and_normalization_flags_are_checked() {
        #[derive(Default)]
        struct Collision;
        impl Hasher for Collision {
            fn finish(&self) -> u64 {
                0
            }
            fn write(&mut self, _: &[u8]) {}
        }
        let head = symbol!("regional_key_collision");
        let first = function!(head, 1);
        let second = function!(head, 2);
        let mut raw = Atom::new();
        raw.to_fun(head);
        let normalized = FunctionBuilder::new(head).finish();
        assert_ne!(raw.as_view().get_data(), normalized.as_view().get_data());

        let mut keys = HashMap::<_, _, BuildHasherDefault<Collision>>::default();
        for (position, atom) in [&first, &second, &raw, &normalized].into_iter().enumerate() {
            let carrier = AdmittedInterfaces::for_expression(atom.as_view());
            keys.insert(carrier.key(atom.as_view()), position);
        }
        for (position, atom) in [&first, &second, &raw, &normalized].into_iter().enumerate() {
            assert_eq!(keys.get(atom.as_view().get_data()), Some(&position));
        }
    }

    #[test]
    fn projecting_independent_factors_reads_the_smaller_side_of_each_join() {
        let head = symbol!("projected_factor");
        let shared = function!(symbol!("projected_argument"), 997);
        let factors: Vec<_> = (1..=32).map(|i| function!(head, i, &shared)).collect();
        let scalar = PartialStructure::from_logical_slots([]);
        let sources: Vec<_> = factors
            .iter()
            .map(|factor| {
                Arc::new(AdmittedInterfaces::from_region(
                    factor.as_view(),
                    scalar.clone(),
                ))
            })
            .collect();
        let first_boundary = Arc::clone(sources[0].interface.as_ref().unwrap());
        let input = Atom::mul_many(&factors);
        INDEXED_KEYS.with(|count| count.set(0));
        REGION_PROBES.with(|count| count.set(0));
        let mut product = AdmittedInterfaces::from_region(input.as_view(), scalar);
        product.project(&sources, &[]);
        assert_eq!(REGION_PROBES.with(Cell::get), sources.len());
        // Root, 32 factors, 32 labels, and the shared two-node argument.
        assert_eq!(INDEXED_KEYS.with(Cell::get), 67);
        assert_eq!(product.projected.len(), factors.len());
        assert!(product.observed.is_empty());
        assert!(Arc::ptr_eq(
            product
                .projected_interface(factors[0].as_view().get_data())
                .unwrap(),
            &first_boundary
        ));
        let product = Arc::new(product);
        INDEXED_KEYS.with(|count| count.set(0));
        let mut extracted = AdmittedInterfaces::for_expression(factors[0].as_view());
        extracted.project(&[Arc::clone(&product)], &[]);
        assert_eq!(INDEXED_KEYS.with(Cell::get), 0);
        assert!(Arc::ptr_eq(
            extracted.interface.as_ref().unwrap(),
            &first_boundary
        ));
        assert!(extracted.get(factors[1].as_view()).is_none());
    }

    #[test]
    fn publishing_an_observed_leaf_skips_indexing_its_enclosing_expression() {
        let leaf = Atom::var(symbol!("indexed_leaf"));
        let enclosing = function!(symbol!("indexed_parent"), &leaf, 42);
        let scalar = PartialStructure::from_logical_slots([]);
        let mut previous = AdmittedInterfaces::from_region(enclosing.as_view(), scalar.clone());
        previous.insert(leaf.as_view(), scalar.clone());
        let previous = Arc::new(previous);

        INDEXED_KEYS.with(|count| count.set(0));
        let mut published = AdmittedInterfaces::for_expression(leaf.as_view());
        published.project(&[previous], &[]);
        assert_eq!(INDEXED_KEYS.with(Cell::get), 0);
        assert_eq!(published.get(leaf.as_view()), Some(scalar));
        assert!(published.get(enclosing.as_view()).is_none());
    }

    #[test]
    fn admission_releases_the_raw_ports_of_a_normalized_away_bracket() {
        test_support::test_initialize();
        let head = tensor_symbol!("regional_bracket_factor");
        let [mu, nu] = [mink!(4, 71914), mink!(4, 71915)];
        let left = function!(head, &mu);
        let right = function!(head, &nu);
        let bracket = spenso::bracket!(&left, &right);
        let tensor = SymbolicTensor::infer(bracket.clone()).unwrap();
        assert_eq!(tensor.expression, &left * &right);
        let leaves = tensor.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
        assert!(!leaves.contains_key(bracket.as_view().get_data()));
        for (leaf, port) in [(left, mu), (right, nu)] {
            assert_eq!(
                leaves.get(leaf.as_view().get_data()),
                Some(&(true, vec![port]))
            );
        }
    }

    #[test]
    fn nested_region_keys_store_each_subtree_payload_once() {
        let head = symbol!("regional_key_nested");
        let depth = 128;
        let mut expression = Atom::one();
        for _ in 0..depth {
            expression = function!(head, expression);
        }
        let mut interfaces = AdmittedInterfaces::for_expression(expression.as_view());
        expression.visitor(&mut |value| {
            interfaces.insert(value, PartialStructure::from_logical_slots([]));
            true
        });
        let storage = interfaces.storage();
        assert_eq!(storage.payload_bytes, expression.as_view().get_byte_size());
        assert_eq!(storage.records, depth + 1);
        assert_eq!(storage.payloads, 1);
    }

    #[test]
    fn collected_coefficient_keeps_self_loop_ports_and_releases_the_input_payload() {
        test_support::test_initialize();
        let port = mink!(4, 71823);
        let spectator = function!(tensor_symbol!("coefficient_self_loop"), &port, &port);
        let selected = function!(tensor_symbol!("coefficient_selection"), mink!(4, 71824));
        let input = SymbolicTensor::infer(&selected * &spectator).unwrap();
        let owner = input.proofs.admitted_interfaces.get().unwrap();
        let payload = Arc::downgrade(&owner.root.as_ref().unwrap().atom);
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let rows = input
            .coefficient_list_with(|leaf| leaf == selected.as_view())
            .unwrap();
        drop(input);
        assert!(payload.upgrade().is_none());
        let [(monomial, coefficient)] = rows.as_slice() else {
            panic!("expected one selected term and its coefficient");
        };
        assert_eq!(monomial.expression, selected);
        assert_eq!(coefficient.expression, spectator);
        assert!(coefficient.is_scalar());
        assert_eq!(
            coefficient.proofs.leaf_interfaces.fact(spectator.as_view()),
            Some((true, vec![port.clone(), port]))
        );
        assert!(
            coefficient
                .proofs
                .leaf_interfaces
                .fact(selected.as_view())
                .is_none()
        );
        coefficient
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn private_producers_share_the_payload_but_own_their_generated_facts() {
        test_support::test_initialize();
        let leaf = function!(tensor_symbol!("private_context"), mink!(4, 71822));
        let source = SymbolicTensor::infer(leaf).unwrap();
        let (local, _) = source.producer_context();
        let source_owner = source.proofs.admitted_interfaces.get().unwrap();
        let local_owner = local.admitted_interfaces.get().unwrap();
        assert!(Arc::ptr_eq(
            &source_owner.root.as_ref().unwrap().atom,
            &local_owner.root.as_ref().unwrap().atom,
        ));
        let emitted = function!(tensor_symbol!("private_emission"), mink!(4, 71824));
        local
            .retain_leaf_ports(emitted.as_view(), vec![mink!(4, 71824)])
            .unwrap();
        assert!(
            source
                .proofs
                .retained_interface(emitted.as_view())
                .is_none()
        );
        assert!(local.retained_interface(emitted.as_view()).is_some());
        assert!(
            !source
                .admitted_leaf_interfaces()
                .lock()
                .unwrap()
                .contains_key(emitted.as_view().get_data())
        );
    }

    #[test]
    fn derived_keys_after_inheritance_share_the_live_expression_payload() {
        let leaf = Atom::var(symbol!("regional_inherited_scalar"));
        let expression = function!(symbol!("regional_inherited_wrapper"), &leaf, 2);
        let scalar = PartialStructure::from_logical_slots([]);
        let mut source = TensorProofs::default();
        source.publish_interfaces(leaf.as_view(), scalar.clone(), []);
        let mut output = TensorProofs::default();
        output.publish_interfaces(expression.as_view(), scalar.clone(), []);
        output.inherit_interfaces([&source]);
        expression.visitor(&mut |view| {
            output.retain_interface(view, scalar.clone());
            true
        });
        assert_eq!(
            output.admitted_interfaces.get().unwrap().storage().payloads,
            1
        );
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
        assert!(storage.payloads <= initial.payloads + repetitions * 4);
        assert!(storage.records <= initial.records + repetitions * 4);
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
        assert_eq!(
            storage.payload_bytes,
            tensor.expression.as_view().get_byte_size()
        );
        assert_eq!(storage.payloads, 1);
        let mut distinct = AHashSet::new();
        tensor.expression.visitor(&mut |node| {
            distinct.insert(node.to_owned());
            true
        });
        assert!(storage.records <= distinct.len(), "{storage:?}");
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
                .map(|root| &root.atom)
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
                .records,
            storage.records,
            "publishing a result must not mutate the source facts",
        );
        drop(tensor);
        assert!(previous_root.upgrade().is_none());
        let storage = reduced.proofs.admitted_interfaces.get().unwrap().storage();
        assert_eq!(storage.payloads, 1);
        assert_eq!(
            storage.payload_bytes,
            reduced.expression.as_view().get_byte_size()
        );
        assert!(storage.records < 256);
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
    fn routed_slash_carries_written_spinor_ports_to_each_linear_term() {
        test_support::test_initialize();
        let q = vector_symbol!("regional_slash_q");
        let [left, right] = [crate::bis!(4, 71842), crate::bis!(4, 71841)];
        let component = |owner: i64| function!(q, owner, mink!(4));
        let slash = |vector: &Atom| function!(crate::dirac::AGS.gamma, &left, &right, vector);
        let original = component(0);
        let source = SymbolicTensor::infer_indexed(slash(&original)).unwrap();
        let scale = Atom::var(symbol!("regional_slash_scale"; Scalar));
        INFERENCE_CALLS.with(|calls| calls.set(0));
        INFERENCE_VISITS.with(|calls| calls.set(0));
        let routed = source
            .map_produced_vector_components(q, |node| {
                (node == original.as_view()).then(|| &scale * component(1) - component(2))
            })
            .unwrap();
        assert_eq!(
            routed.expression(),
            &(&scale * slash(&component(1)) - slash(&component(2)))
        );
        let routed = SymbolicTensor::<PartialStructure>::from(routed);
        for owner in [1, 2] {
            let leaf = slash(&component(owner));
            let (_, ports) = routed.proofs.leaf_interfaces.get().unwrap().lock().unwrap()
                [leaf.as_view().get_data()]
            .clone();
            assert_eq!(ports, [left.clone(), right.clone()]);
            assert_eq!(routed.admitted_region(leaf.as_view()).unwrap().rank(), 2);
        }
        routed
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(INFERENCE_CALLS.with(Cell::get), 0);
        assert_eq!(INFERENCE_VISITS.with(Cell::get), 0);
    }

    #[test]
    fn routed_fixed_components_keep_empty_ports_and_scalar_owners_stay_opaque() {
        test_support::test_initialize();
        let q = vector_symbol!("regional_fixed_q");
        let fixed = function!(spenso::structure::abstract_index::AIND_SYMBOLS.cind, 0);
        let scalar_owner = symbol!("regional_fixed_owner"; Scalar);
        let old = function!(q, 1, &fixed);
        let new = function!(q, function!(scalar_owner, function!(q, 1)), &fixed);
        let metadata = function!(scalar_owner, function!(q, 2));
        let source = SymbolicTensor::infer_indexed(&old * &metadata).unwrap();
        let routed = source
            .map_produced_vector_components(q, |node| {
                if node == old.as_view() {
                    Some(new.clone())
                } else {
                    assert_eq!(node, function!(q, 2).as_view());
                    Some(function!(q, 3))
                }
            })
            .unwrap();
        assert_eq!(
            routed.expression(),
            &(&new * function!(scalar_owner, function!(q, 3)))
        );
        let routed = SymbolicTensor::<PartialStructure>::from(routed);
        assert!(routed.is_scalar());
        let (_, written) =
            &routed.proofs.leaf_interfaces.get().unwrap().lock().unwrap()[new.as_view().get_data()];
        assert!(written.is_empty());
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
    fn analytic_scalar_weights_preserve_new_tensor_regions() {
        test_support::test_initialize();
        let scale = Atom::var(symbol!("analytic_weight_scale"));
        let vector = function!(tensor_symbol!("analytic_weight_leaf"), mink!(4, 71849));
        let tensor = SymbolicTensor::infer(vector.clone()).unwrap();
        let logarithm = function!(Symbol::LOG, Atom::num(4) / &scale);
        let weighted = &vector * logarithm;
        INFERENCE_CALLS.with(|count| count.set(0));
        assert_eq!(
            tensor.proofs.inspect_interface(weighted.as_view()),
            Some(tensor.structure.clone())
        );
        assert!(
            tensor
                .proofs
                .inspect_interface(function!(Symbol::LOG, &vector).as_view())
                .is_none()
        );
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    }

    #[test]
    fn inspecting_a_rewrite_does_not_publish_its_temporary_regions_or_ports() {
        test_support::test_initialize();
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let head = SPENSO_TAG.rank_one_tensor_symbol("inspected::q");
        let vector = function!(head, mink!(4, 71804));
        let product = &vector * Atom::num(2);
        let interface = source.proofs.inspect_interface(product.as_view()).unwrap();
        assert_eq!(interface.logical_slots().len(), 1);
        assert!(
            source
                .proofs
                .retained_interface(product.as_view())
                .is_none()
        );
        assert!(source.proofs.retained_interface(vector.as_view()).is_none());
        assert!(
            !source
                .admitted_leaf_interfaces()
                .lock()
                .unwrap()
                .contains_key(vector.as_view().get_data())
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
