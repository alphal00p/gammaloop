//! Storage for a constructor's admitted written ports, including closed pairs.

#[cfg(test)]
use std::cell::Cell;
use std::collections::HashMap;
use std::mem;
use std::slice::Iter;
use std::sync::{Arc, Mutex, MutexGuard, OnceLock};

use super::regions;
use ahash::AHashMap;
use spenso::network::parsing::LeafInterfaceCache;
use symbolica::atom::representation::ListIterator;
use symbolica::atom::{Atom, AtomCore, AtomView};

type LeafFact = (bool, Vec<Atom>);

#[derive(Clone, Default)]
pub(crate) struct RawLeafInterfaces(OnceLock<LeafSource>);

#[derive(Clone)]
enum LeafSource {
    Writable(LeafInterfaceCache),
    Published(Arc<PublishedLeaves>),
}

struct PublishedLeaves {
    payload: Arc<Atom>,
    // Sorted by exact expression bytes, independently of their position in the payload.
    direct: Box<[DirectLeaf]>,
    // Parser boundaries can include constructed ports or require a different order.
    exceptional: HashMap<Atom, LeafFact>,
    // Only an explicitly writable user of this table needs to promote it.
    writable: OnceLock<LeafInterfaceCache>,
}

struct DirectLeaf {
    offset: u32,
    length: u16,
    arguments: u16,
}

impl DirectLeaf {
    fn bytes<'a>(&self, payload: &'a [u8]) -> &'a [u8] {
        let start = self.offset as usize;
        &payload[start..start + usize::from(self.length)]
    }
}

#[derive(Clone, Copy)]
enum Ports<'a> {
    Arguments(u16),
    Supplied(bool, &'a [Atom]),
}

/// Written order is independent of the leaf's free boundary and includes
/// repeated slots. A projected mask borrows the constructor's own arguments.
pub(super) enum WrittenPorts<'a> {
    Supplied(Iter<'a, Atom>),
    Arguments {
        remaining: u16,
        arguments: ListIterator<'a>,
    },
}

impl<'a> Iterator for WrittenPorts<'a> {
    type Item = AtomView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Supplied(ports) => ports.next().map(Atom::as_view),
            Self::Arguments {
                remaining,
                arguments,
            } => {
                if *remaining == 0 {
                    return None;
                }
                let skip = remaining.trailing_zeros();
                *remaining >>= skip + 1;
                arguments.nth(skip as usize)
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = match self {
            Self::Supplied(ports) => ports.len(),
            Self::Arguments { remaining, .. } => remaining.count_ones() as usize,
        };
        (n, Some(n))
    }
}

impl ExactSizeIterator for WrittenPorts<'_> {}

#[cfg(test)]
thread_local! {
    static OWNED_DIRECT_READS: Cell<usize> = const { Cell::new(0) };
}

impl<'a> Ports<'a> {
    // The high bit distinguishes a negative fact from an admitted zero-port leaf.
    const NEGATIVE: u16 = u16::MAX;

    fn arguments(self, node: AtomView<'_>) -> Option<u16> {
        let (positive, ports) = match self {
            Self::Arguments(mask) => return Some(mask),
            Self::Supplied(positive, ports) => (positive, ports),
        };
        if !positive {
            return ports.is_empty().then_some(Self::NEGATIVE);
        }
        if ports.is_empty() {
            return Some(0);
        }
        let AtomView::Fun(function) = node else {
            return None;
        };
        let mut written = ports.iter().peekable();
        let mut mask = 0;
        for (position, argument) in function.iter().enumerate() {
            if written.next_if(|port| port.as_view() == argument).is_some() {
                if position >= u16::BITS as usize - 1 {
                    return None;
                }
                mask |= 1 << position;
                if written.peek().is_none() {
                    return Some(mask);
                }
            }
        }
        None
    }

    fn written<'b>(self, node: AtomView<'b>) -> (bool, WrittenPorts<'b>)
    where
        'a: 'b,
    {
        match self {
            Self::Supplied(positive, ports) => (positive, WrittenPorts::Supplied(ports.iter())),
            Self::Arguments(mask @ (Self::NEGATIVE | 0)) => {
                (mask == 0, WrittenPorts::Supplied([].iter()))
            }
            Self::Arguments(remaining) => {
                let AtomView::Fun(function) = node else {
                    unreachable!("recorded argument mask belongs to a function")
                };
                (
                    true,
                    WrittenPorts::Arguments {
                        remaining,
                        arguments: function.iter(),
                    },
                )
            }
        }
    }

    fn to_owned(self, node: AtomView<'_>) -> LeafFact {
        if let Self::Supplied(positive, ports) = self {
            return (positive, ports.to_vec());
        }
        #[cfg(test)]
        OWNED_DIRECT_READS.with(|count| count.set(count.get() + 1));
        let (positive, ports) = self.written(node);
        (positive, ports.map(|port| port.to_owned()).collect())
    }
}

impl PublishedLeaves {
    fn key(&self, entry: &DirectLeaf) -> &[u8] {
        entry.bytes(self.payload.as_view().get_data())
    }

    fn ports(&self, bytes: &[u8]) -> Option<Ports<'_>> {
        if let Some((positive, ports)) = self.exceptional.get(bytes) {
            return Some(Ports::Supplied(*positive, ports));
        }
        let index = self
            .direct
            .binary_search_by(|entry| self.key(entry).cmp(bytes))
            .ok()?;
        Some(Ports::Arguments(self.direct[index].arguments))
    }

    fn to_cache(&self) -> LeafInterfaceCache {
        let mut facts = self.exceptional.clone();
        facts.reserve(self.direct.len());
        for entry in &self.direct {
            let node = AtomView::from(self.key(entry));
            facts.insert(
                node.to_owned(),
                Ports::Arguments(entry.arguments).to_owned(node),
            );
        }
        Arc::new(Mutex::new(facts))
    }
}

enum ReadLeaves<'a> {
    Table(MutexGuard<'a, HashMap<Atom, LeafFact>>),
    Published(&'a PublishedLeaves),
}

impl LeafSource {
    fn read(&self) -> ReadLeaves<'_> {
        let cache = match self {
            Self::Writable(cache) => cache,
            Self::Published(published) => match published.writable.get() {
                Some(cache) => cache,
                None => return ReadLeaves::Published(published),
            },
        };
        ReadLeaves::Table(cache.lock().expect("leaf interface cache is not poisoned"))
    }
}

impl ReadLeaves<'_> {
    fn len(&self) -> usize {
        match self {
            Self::Table(table) => table.len(),
            Self::Published(value) => value.direct.len() + value.exceptional.len(),
        }
    }

    fn ports(&self, bytes: &[u8]) -> Option<Ports<'_>> {
        match self {
            Self::Table(table) => table
                .get(bytes)
                .map(|(positive, ports)| Ports::Supplied(*positive, ports)),
            Self::Published(value) => value.ports(bytes),
        }
    }

    fn for_each(&self, mut visit: impl FnMut(&[u8], Ports<'_>)) {
        let table = match self {
            Self::Table(table) => table,
            Self::Published(value) => {
                for entry in &value.direct {
                    visit(value.key(entry), Ports::Arguments(entry.arguments));
                }
                &value.exceptional
            }
        };
        for (atom, (positive, ports)) in table.iter() {
            visit(atom.as_view().get_data(), Ports::Supplied(*positive, ports));
        }
    }
}

impl From<LeafInterfaceCache> for RawLeafInterfaces {
    fn from(cache: LeafInterfaceCache) -> Self {
        Self(OnceLock::from(LeafSource::Writable(cache)))
    }
}

impl RawLeafInterfaces {
    pub(super) fn is_present(&self) -> bool {
        self.0.get().is_some()
    }

    pub(super) fn same_source(&self, other: &Self) -> bool {
        match (self.0.get(), other.0.get()) {
            (None, None) => true,
            (Some(LeafSource::Writable(left)), Some(LeafSource::Writable(right))) => {
                Arc::ptr_eq(left, right)
            }
            (Some(LeafSource::Published(left)), Some(LeafSource::Published(right))) => {
                Arc::ptr_eq(left, right)
            }
            _ => false,
        }
    }

    fn writable(source: &LeafSource) -> &LeafInterfaceCache {
        match source {
            LeafSource::Writable(value) => value,
            LeafSource::Published(value) => value.writable.get_or_init(|| value.to_cache()),
        }
    }

    pub(crate) fn get(&self) -> Option<&LeafInterfaceCache> {
        self.0.get().map(Self::writable)
    }

    pub(super) fn get_or_init(
        &self,
        create: impl FnOnce() -> LeafInterfaceCache,
    ) -> &LeafInterfaceCache {
        Self::writable(self.0.get_or_init(|| LeafSource::Writable(create())))
    }

    /// Supply Spenso's public map to a reader for this operation. A read never
    /// leaves an expanded duplicate on a long-lived admitted expression.
    pub(crate) fn read_cache(&self) -> Option<LeafInterfaceCache> {
        Some(match self.0.get()? {
            LeafSource::Writable(value) => Arc::clone(value),
            LeafSource::Published(value) => value
                .writable
                .get()
                .cloned()
                .unwrap_or_else(|| value.to_cache()),
        })
    }

    pub(super) fn inherit_replacement(
        &self,
        output: AtomView<'_>,
        previous: AtomView<'_>,
        other: &Self,
    ) {
        let Some(source) = other.0.get() else { return };
        self.0.get_or_init(|| {
            if output.get_data() == previous.get_data() {
                return source.clone();
            }
            // An independently admitted scalar or an extracted region can
            // arrive without a raw table. Keep only the surviving written
            // declarations and release all storage belonging to the old root.
            let source = source.read();
            let mut retained = HashMap::new();
            if source.len() != 0 {
                output.visitor(&mut |node| {
                    if let Some(ports) = source.ports(node.get_data())
                        && !retained.contains_key(node.get_data())
                    {
                        retained.insert(node.to_owned(), ports.to_owned(node));
                    }
                    true
                });
            }
            LeafSource::Writable(Arc::new(Mutex::new(retained)))
        });
    }

    pub(super) fn contains(&self, node: AtomView<'_>) -> bool {
        self.0
            .get()
            .is_some_and(|source| source.read().ports(node.get_data()).is_some())
    }

    pub(super) fn with_fact<T>(
        &self,
        node: AtomView<'_>,
        read: impl FnOnce(bool, WrittenPorts<'_>) -> T,
    ) -> Option<T> {
        let source = self.0.get()?.read();
        let (positive, ports) = source.ports(node.get_data())?.written(node);
        Some(read(positive, ports))
    }

    /// Read one consistent table while an archive visits its reachable expression keys.
    #[cfg(feature = "bincode")]
    pub(super) fn reader(&self) -> impl Fn(AtomView<'_>) -> Option<LeafFact> + '_ {
        let source = self.0.get().map(LeafSource::read);
        move |node| Some(source.as_ref()?.ports(node.get_data())?.to_owned(node))
    }

    pub(crate) fn fact(&self, node: AtomView<'_>) -> Option<LeafFact> {
        let source = self.0.get()?.read();
        Some(source.ports(node.get_data())?.to_owned(node))
    }

    pub(super) fn project<'a>(
        payload: &Arc<Atom>,
        nodes: &[AtomView<'a>],
        positions: Option<&AHashMap<&[u8], usize>>,
        sources: &[Self],
    ) -> Self {
        // Multi-factor products live on generated diagrams together with their
        // original factors. Store their selected facts in the output's own bytes.
        let compact = sources.len() > 2 && nodes.len() >= 16;
        let mut direct = Vec::new();
        let mut exceptional = HashMap::new();
        let mut selected = vec![false; nodes.len()];
        let bytes = payload.as_view().get_data();
        let mut retain = |position: usize, ports: Ports<'_>| {
            if mem::replace(&mut selected[position], true) {
                return;
            }
            let node = nodes[position];
            let key = node.get_data();
            if (compact || matches!(ports, Ports::Arguments(_)))
                && let Some(offset) = regions::payload_offset(bytes, key)
                && let Ok(length) = u16::try_from(key.len())
                && let Some(arguments) = ports.arguments(node)
            {
                direct.push(DirectLeaf {
                    offset,
                    length,
                    arguments,
                });
                return;
            }
            exceptional.insert(node.to_owned(), ports.to_owned(node));
        };
        for source in sources.iter().filter_map(|owner| owner.0.get()) {
            let source = source.read();
            if let Some(positions) = positions
                && source.len() < nodes.len()
            {
                source.for_each(|bytes, ports| {
                    if let Some(&position) = positions.get(bytes) {
                        retain(position, ports);
                    }
                });
            } else {
                for (position, &node) in nodes.iter().enumerate() {
                    if let Some(ports) = source.ports(node.get_data()) {
                        retain(position, ports);
                    }
                }
            }
        }
        if direct.is_empty() {
            return Self::from(Arc::new(Mutex::new(exceptional)));
        }
        direct.sort_unstable_by(|left, right| left.bytes(bytes).cmp(right.bytes(bytes)));
        Self(OnceLock::from(LeafSource::Published(Arc::new(
            PublishedLeaves {
                payload: Arc::clone(payload),
                direct: direct.into_boxed_slice(),
                exceptional,
                writable: OnceLock::new(),
            },
        ))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "bincode")]
    use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS};
    use crate::tensor::{SymbolicTensor, TensorProofs};
    use crate::test_support;
    #[cfg(feature = "bincode")]
    use bincode::config;
    use spenso::structure::TensorStructure;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::partial::PartialStructure;
    use spenso::structure::representation::{Minkowski, RepName};
    use spenso::structure::slot::{IsAbstractSlot, SlotMatcher};
    use spenso::{mink, tensor_symbol};
    use std::iter;
    use std::slice;
    use symbolica::atom::FunctionBuilder;
    #[cfg(feature = "bincode")]
    use symbolica::state::StateMap;
    use symbolica::{function, symbol};

    fn project(root: Atom, sources: &[RawLeafInterfaces]) -> RawLeafInterfaces {
        let root = Arc::new(root);
        let (nodes, positions) = regions::output_positions(root.as_view(), false);
        RawLeafInterfaces::project(&root, &nodes, positions.as_ref(), sources)
    }

    #[test]
    fn raw_port_projection_preserves_order_multiplicity_and_independent_writes() {
        test_support::test_initialize();
        let head = symbol!("packed_written_ports");
        let left = mink!(4, 97641);
        let right = mink!(4, 97642);
        let leaves = (0..10)
            .map(|label| function!(head, label, &left, &right, &right))
            .collect::<Vec<_>>();
        let root = FunctionBuilder::new(symbol!("packed_written_port_parent"))
            .add_args(&leaves)
            .finish();
        let sources = (0..3)
            .map(|_| LeafInterfaceCache::default())
            .collect::<Vec<_>>();
        let mut expected = HashMap::new();
        for (index, leaf) in leaves.iter().enumerate() {
            let fact = match index {
                0 => (false, vec![]),
                1 => (true, vec![]),
                2 => (true, vec![right.clone(), left.clone()]),
                3 => (true, vec![right.clone(), right.clone()]),
                4 => (true, vec![Atom::var(symbol!("synthetic_port"))]),
                _ => (true, vec![left.clone(), right.clone()]),
            };
            sources[index % 3]
                .lock()
                .unwrap()
                .insert(leaf.clone(), fact.clone());
            expected.insert(leaf.clone(), fact);
        }
        // A later donor has a conflicting fact for this key.
        sources[2]
            .lock()
            .unwrap()
            .insert(leaves[0].clone(), (true, vec![right.clone()]));
        let inputs = sources
            .into_iter()
            .map(RawLeafInterfaces::from)
            .collect::<Vec<_>>();
        let table = project(root.clone(), &inputs);
        let Some(LeafSource::Published(packed)) = table.0.get() else {
            panic!("expected compact multi-factor publication")
        };
        assert_eq!(packed.direct.len(), 8);
        assert_eq!(packed.exceptional.len(), 2);

        // The shallow parser reads both raw incidences, while region lookup
        // contracts the closed pair. Both parse the original written arguments.
        let proofs = TensorProofs {
            leaf_interfaces: table.clone(),
            ..Default::default()
        };
        let pair = leaves[3].as_view();
        let AtomView::Fun(function) = pair else {
            unreachable!()
        };
        let written = function
            .iter()
            .skip(2)
            .map(|argument| argument.get_data().as_ptr())
            .collect::<Vec<_>>();
        OWNED_DIRECT_READS.with(|count| count.set(0));
        let pointers = table
            .with_fact(pair, |_, ports| {
                ports
                    .map(|port| port.get_data().as_ptr())
                    .collect::<Vec<_>>()
            })
            .unwrap();
        assert_eq!(pointers, written);
        let layout = proofs
            .parser_structure::<AbstractIndex>(pair, &mut SlotMatcher::default(), || {
                unreachable!("named indices")
            })
            .unwrap();
        assert_eq!(layout.canonical().order(), 2);
        assert_eq!(
            proofs.inspect_interface(pair).unwrap().canonical().order(),
            0
        );
        assert_eq!(OWNED_DIRECT_READS.with(Cell::get), 0);

        // A small extraction can transfer the enclosing leaf's exact mask.
        // Source order also applies when a writable donor disagrees with it.
        let competing = RawLeafInterfaces::default();
        competing
            .get_or_init(Default::default)
            .lock()
            .unwrap()
            .insert(leaves[3].clone(), (true, vec![left.clone()]));
        let extract = |sources: &[RawLeafInterfaces]| {
            OWNED_DIRECT_READS.with(|count| count.set(0));
            let child = project(leaves[3].clone(), sources);
            assert_eq!(OWNED_DIRECT_READS.with(Cell::get), 0);
            child
        };
        for child in [
            extract(slice::from_ref(&table)),
            extract(&[table.clone(), competing.clone()]),
        ] {
            assert!(matches!(child.0.get(), Some(LeafSource::Published(_))));
            assert_eq!(
                child.with_fact(pair, |positive, ports| {
                    positive && ports.eq([right.as_view(), right.as_view()])
                }),
                Some(true)
            );
        }
        let child = extract(&[competing, table.clone()]);
        assert_eq!(
            child.with_fact(pair, |positive, ports| {
                positive && ports.eq([left.as_view()])
            }),
            Some(true)
        );

        assert_eq!(*table.read_cache().unwrap().lock().unwrap(), expected);
        table.read_cache().unwrap().lock().unwrap().clear();
        assert!(packed.writable.get().is_none());
        for (leaf, fact) in &expected {
            assert_eq!(table.fact(leaf.as_view()).as_ref(), Some(fact));
        }
        let other = project(root, slice::from_ref(&table));
        assert!(packed.writable.get().is_none());
        assert_eq!(*other.read_cache().unwrap().lock().unwrap(), expected);
        table
            .get_or_init(Default::default)
            .lock()
            .unwrap()
            .remove(&leaves[0]);
        assert!(table.fact(leaves[0].as_view()).is_none());
        assert_eq!(other.fact(leaves[0].as_view()), Some((false, Vec::new())));
    }

    #[test]
    fn compact_records_use_eight_bytes_and_large_declarations_keep_their_exact_ports() {
        test_support::test_initialize();
        assert_eq!(mem::size_of::<DirectLeaf>(), 8);
        let head = symbol!("large_written_declaration");
        let slot = Atom::num(43);
        let make_leaf = |padding: usize| {
            FunctionBuilder::new(head)
                .add_args(iter::repeat_n(Atom::one(), padding))
                .add_arg(&slot)
                .finish()
        };
        let ordinary = make_leaf(1);
        let wide = make_leaf(16);
        let long = FunctionBuilder::new(head)
            .add_arg(&slot)
            .add_args((100_000..140_000).map(Atom::num))
            .finish();
        assert!(long.as_view().get_data().len() > usize::from(u16::MAX));
        let root = FunctionBuilder::new(head)
            .add_args([&ordinary, &wide, &long])
            .add_args((0..16).map(Atom::num))
            .finish();
        let inputs = [&ordinary, &wide, &long].map(|atom| {
            RawLeafInterfaces::from(Arc::new(Mutex::new(HashMap::from([(
                atom.clone(),
                (true, vec![slot.clone()]),
            )]))))
        });
        let result = project(root, &inputs);
        let Some(LeafSource::Published(packed)) = result.0.get() else {
            panic!("expected a compact table");
        };
        assert_eq!(packed.direct.len(), 1);
        assert_eq!(packed.exceptional.len(), 2);
        for leaf in [&ordinary, &wide, &long] {
            assert_eq!(
                result.fact(leaf.as_view()),
                Some((true, vec![slot.clone()]))
            );
        }
    }

    #[cfg(feature = "bincode")]
    #[test]
    fn archived_product_retains_packed_raw_ports_without_promoting_its_source() {
        test_support::test_initialize();
        let head = tensor_symbol!("packed_archived_product_leaf");
        let rep = Minkowski {}.new_rep(4).to_lib();
        let mut factors = Vec::new();
        let mut expected = Vec::new();
        for index in 0..8 {
            let [left, right] = [97398 - index * 2, 97399 - index * 2].map(|index| {
                rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                    .to_atom()
            });
            let atom = function!(head, &right, &left);
            factors.push(SymbolicTensor::infer(atom.clone()).unwrap());
            expected.push((atom, vec![right, left]));
        }
        let product = SymbolicTensor::multiply_many(&factors).unwrap();
        let Some(LeafSource::Published(packed)) = product.proofs.leaf_interfaces.0.get() else {
            panic!("expected product to own compact raw declarations");
        };
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let bytes = bincode::encode_to_vec(&product, config::standard()).unwrap();
        let (restored, used): (SymbolicTensor<PartialStructure>, _) =
            bincode::decode_from_slice_with_context(
                &bytes,
                config::standard(),
                StateMap::default(),
            )
            .unwrap();
        assert_eq!(used, bytes.len());
        assert_eq!(restored.structure, product.structure);
        for (leaf, ports) in expected {
            assert_eq!(
                restored.proofs.leaf_interfaces.fact(leaf.as_view()),
                Some((true, ports))
            );
        }
        assert!(packed.writable.get().is_none());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }

    #[test]
    fn replacement_without_a_raw_table_releases_its_predecessor() {
        test_support::test_initialize();
        let head = tensor_symbol!("retained_raw_replacement");
        let representation = Minkowski {}.new_rep(4).to_lib();
        let factors = (97650..97658)
            .map(|index| {
                let port = representation
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                    .to_atom();
                SymbolicTensor::infer(function!(head, port)).unwrap()
            })
            .collect::<Vec<_>>();
        let product = SymbolicTensor::multiply_many(&factors).unwrap();
        let Some(LeafSource::Published(published)) = product.proofs.leaf_interfaces.0.get() else {
            panic!("expected a compact predecessor");
        };
        let owner = Arc::downgrade(published);
        let without_raw = |source: &SymbolicTensor<PartialStructure>| {
            SymbolicTensor::from_validated_parts(
                source.expression.clone(),
                source.structure.clone(),
                Arc::clone(source.proofs.occurrences.get().unwrap()),
            )
        };
        let same = without_raw(&product);
        same.observe_replacement(&product);
        assert!(
            same.proofs
                .leaf_interfaces
                .same_source(&product.proofs.leaf_interfaces)
        );
        drop(same);

        let extracted = without_raw(&factors[0]);
        extracted.observe_replacement(&product);
        let scalar = SymbolicTensor::infer(Atom::num(37)).unwrap();
        scalar.observe_replacement(&product);
        drop(product);
        assert!(owner.upgrade().is_none());
        assert!(
            scalar
                .proofs
                .leaf_interfaces
                .read_cache()
                .unwrap()
                .lock()
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            extracted
                .proofs
                .leaf_interfaces
                .fact(factors[0].expression.as_view()),
            factors[0]
                .proofs
                .leaf_interfaces
                .fact(factors[0].expression.as_view())
        );
        assert!(
            !extracted
                .proofs
                .leaf_interfaces
                .contains(factors[1].expression.as_view())
        );
    }
}
