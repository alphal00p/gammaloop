use std::{
    collections::HashSet,
    hash::{Hash, Hasher},
    ops::Range,
};

use super::{NodeStorage, NodeStorageOps};
use crate::{
    half_edge::{
        builder::HedgeNodeBuilder,
        involution::{EdgeIndex, Hedge, Involution},
        subgraph::{BaseSubgraph, Inclusion, SuBitGraph, SubSetLike},
        swap::{partition_with, Swap},
        HedgeGraph, HedgeGraphError, NodeIndex, NodeVec,
    },
    tree::{
        parent_pointer::{PPNode, ParentPointerStore},
        Forest, RootData, RootId, TreeNodeId,
    },
};

/// Marks a missing hedge or node: an empty crown, no anchor, a live node.
const NONE: u32 = u32::MAX;

/// Per-node record of a [`NodeForest`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Crown {
    /// Smallest hedge of the crown, [`NONE`] when the crown is empty.
    head: u32,
    len: u32,
    /// The node's anchor hedge, [`NONE`] when it has none.
    anchor: u32,
    /// The node this one was identified into, [`NONE`] while it is live.
    absorbed_into: u32,
}

impl Crown {
    const EMPTY: Self = Crown {
        head: NONE,
        len: 0,
        anchor: NONE,
        absorbed_into: NONE,
    };

    fn is_live(&self) -> bool {
        self.absorbed_into == NONE
    }

    /// Iterates over this crown's ring.
    #[inline]
    fn ring<'a>(&self, links: &'a [Link]) -> RingNeighborIter<'a> {
        RingNeighborIter {
            links,
            cur: self.head,
            remaining: self.len,
        }
    }
}

/// Entry `i` of the class table of a [`NodeForest`], which holds two inverse
/// permutations of the table's positions in one allocation: the class of node
/// slot `i` and the node slot of class `i`. Positions at and behind the node
/// count hold the free classes, so a free class's slot is behind the nodes too.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ClassEntry {
    /// The class of node slot `i`, or a free class.
    class: u32,
    /// The node slot of class `i`.
    slot: u32,
}

impl ClassEntry {
    /// Class `i` at slot `i`.
    fn identity(i: u32) -> Self {
        ClassEntry { class: i, slot: i }
    }
}

/// A half-edge's neighbours in its node's ascending ring.
#[derive(Clone, Copy, Debug)]
struct Link {
    next: u32,
    prev: u32,
}

impl Link {
    const NONE: Self = Link {
        next: NONE,
        prev: NONE,
    };
}

/// The per-half-edge arrays of a [`NodeForest`]: each half-edge's class and its
/// ring links.
struct Rings<'a> {
    label: &'a mut [u32],
    links: &'a mut [Link],
}

impl Rings<'_> {
    /// Labels every half-edge of the ring of `crown` with `class`.
    fn relabel(&mut self, crown: &Crown, class: u32) {
        let mut hedge = crown.head;
        for _ in 0..crown.len {
            self.label[hedge as usize] = class;
            hedge = self.links[hedge as usize].next;
        }
    }

    /// Links `hedge`, larger than every member, at the end of the ring of `crown`.
    #[inline]
    fn append(&mut self, crown: &mut Crown, hedge: u32) {
        let (head, tail) = if crown.len == 0 {
            crown.head = hedge;
            (hedge, hedge)
        } else {
            (crown.head, self.links[crown.head as usize].prev)
        };
        crown.len += 1;
        self.links[hedge as usize].next = head;
        self.links[hedge as usize].prev = tail;
        self.links[tail as usize].next = hedge;
        self.links[head as usize].prev = hedge;
    }

    /// The largest member below `hedge` of the ring of `class` starting at
    /// `head`, with `head < hedge` below the ring's largest member; `from` is a
    /// member below `hedge`.
    ///
    /// Walks the ring up from `from` while scanning labels down and up from
    /// `hedge`, so the cost is three times the smallest of the members walked and
    /// the half-edges between `hedge` and either ring neighbour. Every half-edge
    /// labelled `class` must already be in the ring.
    fn predecessor(&self, class: u32, head: u32, hedge: u32, mut from: u32) -> u32 {
        let (mut down, mut up) = (hedge - 1, hedge + 1);
        loop {
            let next = self.links[from as usize].next;
            if next == head || next > hedge {
                return from;
            }
            from = next;
            if self.label[down as usize] == class {
                return down;
            }
            if self.label[up as usize] == class {
                return self.links[up as usize].prev;
            }
            down -= 1;
            up += 1;
        }
    }

    /// Links `hedge`, labelled `class`, at its ascending position into the
    /// non-empty ring starting at `head`, searching from the member `from` below
    /// `hedge` or `head`; returns the new head.
    #[inline]
    fn insert_into(&mut self, class: u32, head: u32, hedge: u32, from: u32) -> u32 {
        let tail = self.links[head as usize].prev;
        let after = if hedge > tail || hedge < head {
            tail
        } else {
            self.predecessor(class, head, hedge, from)
        };
        let before = self.links[after as usize].next;
        self.links[after as usize].next = hedge;
        self.links[hedge as usize].prev = after;
        self.links[hedge as usize].next = before;
        self.links[before as usize].prev = hedge;
        head.min(hedge)
    }

    /// Links `hedge`, labelled `class`, into the ring of `crown`.
    fn insert(&mut self, crown: &mut Crown, class: u32, hedge: u32) {
        crown.head = if crown.len == 0 {
            self.links[hedge as usize].next = hedge;
            self.links[hedge as usize].prev = hedge;
            hedge
        } else {
            self.insert_into(class, crown.head, hedge, crown.head)
        };
        crown.len += 1;
    }

    /// Unlinks `hedge` from the ring of `crown`.
    fn unlink(&mut self, crown: &mut Crown, hedge: u32) {
        if crown.len == 1 {
            crown.head = NONE;
        } else {
            let (prev, next) = (
                self.links[hedge as usize].prev,
                self.links[hedge as usize].next,
            );
            self.links[prev as usize].next = next;
            self.links[next as usize].prev = prev;
            if crown.head == hedge {
                crown.head = next;
            }
        }
        crown.len -= 1;
    }

    /// Relabels the ring of `absorbed` with `class`, the class of `target`, and
    /// merges it into the ring of `target` in ascending order.
    fn merge(&mut self, target: &mut Crown, class: u32, absorbed: &Crown) {
        if target.len == 0 {
            self.relabel(absorbed, class);
            target.head = absorbed.head;
        } else {
            // Absorbed half-edges ascend, so each one's predecessor is at or
            // above the previous one.
            let (mut hedge, mut head, mut from) = (absorbed.head, target.head, target.head);
            for _ in 0..absorbed.len {
                let next = self.links[hedge as usize].next;
                self.label[hedge as usize] = class;
                head = self.insert_into(class, head, hedge, from);
                from = hedge;
                hedge = next;
            }
            target.head = head;
        }
        target.len += absorbed.len;
    }
}

/// A node store that keeps each crown as a circular ring of half-edges.
///
/// Every half-edge stores the class of its node and its ring neighbours; a class
/// table maps classes to node slots, so [`NodeStorageOps::node_id_ref`] costs two
/// loads and moving a node to another slot never relabels its crown. Each ring is
/// ascending from the node's smallest half-edge in every state, also while
/// identifications are pending, so crown order never depends on history. Per node
/// the store keeps its class, the ring head and length, an optional anchor
/// half-edge and whether the node was absorbed by an identification; payloads
/// live in a separate [`NodeVec`].
///
/// # Classes
///
/// Class ids are stable while their node lives, whatever its slot. Builds,
/// deletions and extractions number classes like slots. The classes of the nodes
/// dropped by [`NodeStorageOps::forget_identification_history`] become free and
/// [`NodeStorageOps::extend_mut`] reuses them before minting new ones, so the
/// class table never outgrows the largest node count since the last rebuild.
///
/// # Node counts
///
/// [`NodeStorageOps::identify_nodes`] never renumbers: absorbed nodes keep their
/// slot and payload until [`NodeStorageOps::forget_identification_history`].
/// [`Self::node_bound`] counts every slot and is the size for dense side tables
/// indexed by [`NodeIndex`]; [`Self::n_live`] counts the nodes that are not
/// absorbed. The two agree once history is forgotten.
///
/// # Identification
///
/// The survivor is the requested node with the largest crown, the first such in
/// request order on ties; repeated and absorbed nodes resolve to the node that
/// absorbed them. Absorbed half-edges are relabelled immediately (weighted
/// quick-find, so a half-edge is relabelled `O(log H)` times) and spliced into
/// the survivor's ring in ascending position, taking the survivor's class. The
/// survivor keeps its anchor or adopts the first absorbed one.
///
/// # Compaction
///
/// Forgetting, deleting and extracting renumber nodes with the shared tail-fill
/// kernel [`partition_with`]: holes are filled from the tail, so live nodes may
/// change index. Forgetting moves node records and class table entries only, in
/// `O(1)` per moved node. Moving half-edges or renumbering nodes (delete, extract,
/// extract_nodes, swaps) requires forgotten history. Live nodes without
/// half-edges are kept like any other node.
///
/// # Anchors
///
/// A node's anchor is an optional half-edge of its crown set with
/// [`Self::set_anchor`]; nodes start without one. It follows its half-edge
/// through swaps and compaction and falls back to the smallest remaining crown
/// half-edge when its own half-edge leaves the node.
///
/// # Equality
///
/// Comparison and hashing use the canonical image: each half-edge's node, which
/// nodes are live, and the payloads. Classes, ring links, anchors and the node an
/// absorbed node was merged into are not compared.
#[derive(Clone, Debug)]
pub struct NodeForest<V> {
    /// The class of each half-edge's node.
    label: Vec<u32>,
    links: Vec<Link>,
    crowns: Vec<Crown>,
    classes: Vec<ClassEntry>,
    n_absorbed: u32,
    data: NodeVec<V>,
}

/// Iterates over a [`NodeForest`] crown in ascending half-edge order.
#[derive(Clone, Debug)]
pub struct RingNeighborIter<'a> {
    links: &'a [Link],
    cur: u32,
    remaining: u32,
}

impl Iterator for RingNeighborIter<'_> {
    type Item = Hedge;

    #[inline]
    fn next(&mut self) -> Option<Hedge> {
        if self.remaining == 0 {
            return None;
        }
        let hedge = self.cur;
        self.cur = self.links[hedge as usize].next;
        self.remaining -= 1;
        Some(Hedge(hedge as usize))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining as usize, Some(self.remaining as usize))
    }
}

impl ExactSizeIterator for RingNeighborIter<'_> {}

impl<V> NodeStorage for NodeForest<V> {
    type NodeData = V;
    type Storage<N> = NodeForest<N>;
    type Neighbors = SuBitGraph;
    type NeighborsIter<'a>
        = RingNeighborIter<'a>
    where
        Self: 'a;
}

const HISTORY_RETAINED: &str = "moving half-edges or renumbering nodes of a NodeForest \
     requires forget_identification_history() after identify_nodes()";

fn index_u32(index: usize) -> u32 {
    u32::try_from(index)
        .ok()
        .filter(|&i| i != NONE)
        .expect("NodeForest indices must fit below u32::MAX")
}

impl<V> NodeForest<V> {
    /// Builds the store from each half-edge's node, with live nodes without
    /// anchors and each node's class equal to its slot.
    fn from_labels(mut label: Vec<u32>, data: NodeVec<V>) -> Self {
        let n_hedges = label.len();
        let n_nodes: NodeIndex = data.len();
        let mut links = vec![Link::NONE; n_hedges];
        let mut crowns = vec![Crown::EMPTY; n_nodes.0];
        let mut rings = Rings {
            label: &mut label,
            links: &mut links,
        };
        for hedge in 0..n_hedges {
            let node = rings.label[hedge];
            rings.append(&mut crowns[node as usize], hedge as u32);
        }
        NodeForest {
            label,
            links,
            crowns,
            classes: (0..n_nodes.0 as u32).map(ClassEntry::identity).collect(),
            n_absorbed: 0,
            data,
        }
    }

    #[inline]
    fn ring(&self, node: usize) -> RingNeighborIter<'_> {
        self.crowns[node].ring(&self.links)
    }

    /// The node slot of `hedge`.
    #[inline]
    fn node_of(&self, hedge: usize) -> u32 {
        let class = self.label[hedge] as usize;
        debug_assert!(class < self.classes.len());
        // SAFETY: labels enter through `from_labels`, which indexes the node
        // records, as many as classes, with each of them, or are copied from the
        // class table, which never shrinks; so every label is below its length.
        unsafe { self.classes.get_unchecked(class).slot }
    }

    /// Each half-edge's node slot: the labels of the canonical image.
    fn nodes(&self) -> impl ExactSizeIterator<Item = u32> + '_ {
        self.label
            .iter()
            .map(|&class| self.classes[class as usize].slot)
    }

    /// Number of node slots, live or absorbed: the bound of every [`NodeIndex`]
    /// and the size for dense side tables.
    pub fn node_bound(&self) -> usize {
        self.crowns.len()
    }

    /// Number of nodes not absorbed by a pending identification.
    pub fn n_live(&self) -> usize {
        self.crowns.len() - self.n_absorbed as usize
    }

    /// The anchor half-edge of `node`, if it has one.
    pub fn anchor(&self, node: NodeIndex) -> Option<Hedge> {
        let anchor = self.crowns[node.0].anchor;
        (anchor != NONE).then_some(Hedge(anchor as usize))
    }

    /// Sets or clears the anchor of `node` in `O(1)`.
    ///
    /// # Panics
    ///
    /// If `anchor` is not in the crown of `node`.
    pub fn set_anchor(&mut self, node: NodeIndex, anchor: Option<Hedge>) {
        self.crowns[node.0].anchor = match anchor {
            Some(hedge) => {
                assert_eq!(
                    self.node_of(hedge.0) as usize,
                    node.0,
                    "anchor {hedge} is not in the crown of {node}"
                );
                hedge.0 as u32
            }
            None => NONE,
        };
    }

    /// The live node that `node` was, possibly transitively, identified into,
    /// compressing the path on the way.
    fn live_of(&mut self, node: u32) -> u32 {
        let mut live = node;
        while self.crowns[live as usize].absorbed_into != NONE {
            live = self.crowns[live as usize].absorbed_into;
        }
        let mut node = node;
        while node != live {
            let next = self.crowns[node as usize].absorbed_into;
            self.crowns[node as usize].absorbed_into = live;
            node = next;
        }
        live
    }

    /// The half-edge arrays, borrowed apart from the node records.
    fn rings(&mut self) -> (Rings<'_>, &mut [Crown]) {
        let rings = Rings {
            label: &mut self.label,
            links: &mut self.links,
        };
        (rings, &mut self.crowns)
    }

    /// Moves the crown of the live node `absorbed` into the ring of `survivor`.
    fn absorb(&mut self, absorbed: u32, survivor: u32) {
        let crown = std::mem::replace(
            &mut self.crowns[absorbed as usize],
            Crown {
                absorbed_into: survivor,
                ..Crown::EMPTY
            },
        );
        self.n_absorbed += 1;
        let class = self.classes[survivor as usize].class;
        let (mut rings, crowns) = self.rings();
        let target = &mut crowns[survivor as usize];
        rings.merge(target, class, &crown);
        if target.anchor == NONE {
            target.anchor = crown.anchor;
        }
    }

    fn count_absorbed(&self) -> u32 {
        self.crowns.iter().filter(|crown| !crown.is_live()).count() as u32
    }

    /// Panics while identifications are pending: absorbed nodes point at their
    /// survivors by index, so half-edges and nodes must not move.
    fn require_forgotten(&self) {
        assert_eq!(self.n_absorbed, 0, "{HISTORY_RETAINED}");
    }

    /// Swaps two node slots with their classes, crowns, anchors, states and
    /// payloads in `O(1)`: half-edges keep their class.
    fn swap_slots(&mut self, a: usize, b: usize) {
        self.crowns.swap(a, b);
        let (class_a, class_b) = (self.classes[a].class, self.classes[b].class);
        self.classes[a].class = class_b;
        self.classes[b].class = class_a;
        self.classes[class_a as usize].slot = b as u32;
        self.classes[class_b as usize].slot = a as u32;
        self.data.swap(NodeIndex(a), NodeIndex(b));
    }

    /// The order of the half-edges after moving those selected by `moved` behind
    /// the others with the edge store's swaps, and the first moved position.
    fn hedge_order(&self, moved: impl Fn(Hedge) -> bool) -> (Vec<u32>, usize) {
        let n_hedges = self.label.len();
        let mut order: Vec<u32> = (0..n_hedges as u32).collect();
        let left = partition_with(
            &mut order,
            0,
            n_hedges,
            |_, position| !moved(Hedge(position)),
            |order, a, b| order.swap(a, b),
        );
        (order, left)
    }

    /// Partitions the positions `start..` of `order`, the old node at each
    /// position, so that the nodes accepted by `keep` come first, permuting the
    /// payloads along; returns the split point.
    fn partition_nodes(
        &mut self,
        order: &mut Vec<u32>,
        start: usize,
        keep: impl Fn(u32) -> bool,
    ) -> usize {
        let end = order.len();
        partition_with(
            &mut (order, &mut self.data),
            start,
            end,
            |(order, _), position| keep(order[position]),
            |(order, data), a, b| {
                order.swap(a, b);
                data.swap(NodeIndex(a), NodeIndex(b));
            },
        )
    }

    /// Moves the half-edges selected by `moved` behind the others with the edge
    /// store's swaps and orders the nodes as `[untouched | split | inside]`, as
    /// the other stores do. Returns the renumbering, the first moved half-edge
    /// and the starts of the split and inside node ranges.
    fn partition_for_removal(
        &mut self,
        moved: impl Fn(Hedge) -> bool,
    ) -> (Compaction, usize, usize, usize) {
        self.require_forgotten();
        let (order, left) = self.hedge_order(moved);
        let mut n_moved = vec![0; self.crowns.len()];
        for &hedge in &order[left..] {
            n_moved[self.node_of(hedge as usize) as usize] += 1;
        }
        let split_node: Vec<bool> = (self.crowns.iter().zip(&n_moved))
            .map(|(crown, &n_moved)| n_moved < crown.len)
            .collect();
        let mut node_order = (0..self.crowns.len() as u32).collect();
        let split = self.partition_nodes(&mut node_order, 0, |node| n_moved[node as usize] == 0);
        let inside = self.partition_nodes(&mut node_order, split, |node| split_node[node as usize]);
        (Compaction::new(order, node_order), left, split, inside)
    }

    fn with_data<V2>(&self, data: NodeVec<V2>) -> NodeForest<V2> {
        NodeForest {
            label: self.label.clone(),
            links: self.links.clone(),
            crowns: self.crowns.clone(),
            classes: self.classes.clone(),
            n_absorbed: self.n_absorbed,
            data,
        }
    }

    fn into_data<V2>(self, data: NodeVec<V2>) -> NodeForest<V2> {
        NodeForest {
            label: self.label,
            links: self.links,
            crowns: self.crowns,
            classes: self.classes,
            n_absorbed: self.n_absorbed,
            data,
        }
    }

    /// Whether node `i` is live, for each node: part of the canonical image.
    fn live_mask(&self) -> impl Iterator<Item = bool> + '_ {
        self.crowns.iter().map(Crown::is_live)
    }

    /// Rebuilds a store from its canonical parts, the input of every
    /// deserializer: each half-edge's node, and each node's anchor and merge
    /// target, [`NONE`] for none. Rejects inconsistent parts.
    #[cfg(any(feature = "serde", feature = "rkyv"))]
    fn from_parts(
        label: Vec<u32>,
        anchor: Vec<u32>,
        absorbed_into: Vec<u32>,
        data: Vec<V>,
    ) -> Result<Self, String> {
        let n_nodes = data.len();
        if anchor.len() != n_nodes || absorbed_into.len() != n_nodes {
            return Err("one anchor and state per node expected".to_owned());
        }
        if label.iter().any(|&node| node as usize >= n_nodes) {
            return Err("half-edge labelled with a missing node".to_owned());
        }
        let mut store = NodeForest::from_labels(label, NodeVec(data));
        for (crown, (anchor, absorbed_into)) in store
            .crowns
            .iter_mut()
            .zip(anchor.into_iter().zip(absorbed_into))
        {
            crown.anchor = anchor;
            crown.absorbed_into = absorbed_into;
        }
        store.n_absorbed = store.count_absorbed();
        store.check_nodes().map_err(|error| error.to_string())?;
        Ok(store)
    }
}

impl<V: PartialEq> PartialEq for NodeForest<V> {
    fn eq(&self, other: &Self) -> bool {
        self.nodes().eq(other.nodes())
            && self.live_mask().eq(other.live_mask())
            && self.data == other.data
    }
}

impl<V: Eq> Eq for NodeForest<V> {}

impl<V: PartialOrd> PartialOrd for NodeForest<V> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.nodes().cmp(other.nodes()) {
            std::cmp::Ordering::Equal => {}
            ordering => return Some(ordering),
        }
        match self.live_mask().cmp(other.live_mask()) {
            std::cmp::Ordering::Equal => {}
            ordering => return Some(ordering),
        }
        self.data.partial_cmp(&other.data)
    }
}

impl<V: Ord> Ord for NodeForest<V> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.nodes()
            .cmp(other.nodes())
            .then_with(|| self.live_mask().cmp(other.live_mask()))
            .then_with(|| self.data.cmp(&other.data))
    }
}

impl<V: Hash> Hash for NodeForest<V> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.label.len().hash(state);
        for node in self.nodes() {
            node.hash(state);
        }
        for live in self.live_mask() {
            live.hash(state);
        }
        self.data.hash(state);
    }
}

impl<V> Swap<Hedge> for NodeForest<V> {
    fn len(&self) -> Hedge {
        Hedge(self.label.len())
    }

    fn is_zero_length(&self) -> bool {
        self.label.is_empty()
    }

    /// Exchanges two half-edges, repositioning each in its new ring; anchors
    /// follow their half-edge. Each insertion costs at most three times the
    /// smaller of the crown size and the index gap to a new ring neighbour.
    fn swap(&mut self, a: Hedge, b: Hedge) {
        self.require_forgotten();
        if a == b {
            return;
        }
        let (a, b) = (a.0 as u32, b.0 as u32);
        let (class_a, class_b) = (self.label[a as usize], self.label[b as usize]);
        let (node_a, node_b) = (
            self.classes[class_a as usize].slot,
            self.classes[class_b as usize].slot,
        );
        if node_a != node_b {
            let (mut rings, crowns) = self.rings();
            rings.unlink(&mut crowns[node_a as usize], a);
            rings.unlink(&mut crowns[node_b as usize], b);
            rings.label[a as usize] = class_b;
            rings.label[b as usize] = class_a;
            rings.insert(&mut crowns[node_a as usize], class_a, b);
            rings.insert(&mut crowns[node_b as usize], class_b, a);
        }
        for node in [node_a, node_b] {
            let anchor = &mut self.crowns[node as usize].anchor;
            if *anchor == a {
                *anchor = b;
            } else if *anchor == b {
                *anchor = a;
            }
            if node_a == node_b {
                break;
            }
        }
    }
}

impl<V> Swap<NodeIndex> for NodeForest<V> {
    fn len(&self) -> NodeIndex {
        NodeIndex(self.crowns.len())
    }

    fn is_zero_length(&self) -> bool {
        self.crowns.is_empty()
    }

    fn swap(&mut self, a: NodeIndex, b: NodeIndex) {
        self.require_forgotten();
        if a != b {
            self.swap_slots(a.0, b.0);
        }
    }
}

impl<V> NodeStorageOps for NodeForest<V> {
    type OpStorage<N> = NodeForest<N>;
    type Base = SuBitGraph;

    fn extend(mut self, other: Self) -> Self {
        self.extend_mut(other);
        self
    }

    fn new_nodevec<'a, V2>(
        &'a self,
        mut node_map: impl FnMut(NodeIndex, Self::NeighborsIter<'a>, &'a Self::NodeData) -> V2,
    ) -> NodeVec<V2> {
        self.data
            .iter()
            .map(|(node, data)| node_map(node, self.ring(node.0), data))
            .collect()
    }

    /// Appends `other` behind `self`, shifting its half-edges and nodes; pending
    /// identifications on either side are kept. Appended nodes take the free
    /// classes of `self` first.
    fn extend_mut(&mut self, other: Self) {
        let hedge_shift = index_u32(self.label.len());
        let node_shift = index_u32(self.crowns.len());
        // The shifted indices must stay below `NONE` as well.
        index_u32(self.label.len() + other.label.len());
        let n_nodes = index_u32(self.crowns.len() + other.crowns.len());
        for class in self.classes.len() as u32..n_nodes {
            self.classes.push(ClassEntry::identity(class));
        }
        let classes = &self.classes[node_shift as usize..n_nodes as usize];
        self.label.extend(
            other
                .label
                .iter()
                .map(|&class| classes[other.classes[class as usize].slot as usize].class),
        );
        let shift = |index: u32, by: u32| if index == NONE { NONE } else { index + by };
        self.links.extend(other.links.iter().map(|link| Link {
            next: link.next + hedge_shift,
            prev: link.prev + hedge_shift,
        }));
        self.crowns.extend(other.crowns.iter().map(|crown| Crown {
            head: shift(crown.head, hedge_shift),
            len: crown.len,
            anchor: shift(crown.anchor, hedge_shift),
            absorbed_into: shift(crown.absorbed_into, node_shift),
        }));
        self.n_absorbed += other.n_absorbed;
        self.data.extend(other.data);
    }

    fn extract<S: SubSetLike<Base = Self::Base>, V2>(
        &mut self,
        subgraph: &S,
        mut split_node: impl FnMut(&Self::NodeData) -> V2,
        mut owned_node: impl FnMut(Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        let (renumbering, left, split, inside) =
            self.partition_for_removal(|hedge| subgraph.includes(&hedge));
        let mut data: NodeVec<V2> = (split..inside)
            .map(|node| split_node(&self.data[NodeIndex(node)]))
            .collect();
        data.0.extend(
            self.data
                .split_off(NodeIndex(inside))
                .into_iter()
                .map(|(_, node)| owned_node(node)),
        );
        let extracted =
            renumbering.side(self, left..self.label.len(), split..self.crowns.len(), data);
        let kept = std::mem::replace(&mut self.data, NodeVec::new());
        *self = renumbering.side(self, 0..left, 0..inside, kept);
        extracted
    }

    fn extract_nodes(&mut self, nodes: impl IntoIterator<Item = NodeIndex>) -> (Self::Base, Self) {
        self.require_forgotten();
        let nodes: HashSet<NodeIndex> = nodes.into_iter().collect();
        let extracted = SuBitGraph::from_hedge_iter(
            nodes.iter().flat_map(|node| self.ring(node.0)),
            self.label.len(),
        );
        let (order, left) = self.hedge_order(|hedge| extracted.includes(&hedge));
        let mut node_order = (0..self.crowns.len() as u32).collect();
        // Selected nodes move out with their half-edges, isolated ones included.
        let kept = self.partition_nodes(&mut node_order, 0, |node| {
            !nodes.contains(&NodeIndex(node as usize))
        });
        let renumbering = Compaction::new(order, node_order);
        let data = self.data.split_off(NodeIndex(kept));
        let store = renumbering.side(self, left..self.label.len(), kept..self.crowns.len(), data);
        let remaining = std::mem::replace(&mut self.data, NodeVec::new());
        *self = renumbering.side(self, 0..left, 0..kept, remaining);
        (extracted, store)
    }

    fn delete<S: SubSetLike<Base = Self::Base>>(&mut self, subgraph: &S) {
        let (renumbering, left, _, inside) =
            self.partition_for_removal(|hedge| subgraph.includes(&hedge));
        self.data.truncate(inside);
        let kept = std::mem::replace(&mut self.data, NodeVec::new());
        *self = renumbering.side(self, 0..left, 0..inside, kept);
    }

    /// Merges `nodes` into the one with the largest crown, the first such in
    /// request order on ties, which receives `node_data_merge`.
    ///
    /// Repeated and absorbed nodes resolve to the live node that absorbed them.
    /// No index changes; absorbed nodes keep their slot and payload. Costs the
    /// relabelled half-edges plus their ring insertions, each bounded by the
    /// survivor's crown members or half-edge indices it skips.
    fn identify_nodes(
        &mut self,
        nodes: &[NodeIndex],
        node_data_merge: Self::NodeData,
    ) -> NodeIndex {
        let (first, rest) = nodes
            .split_first()
            .expect("identify_nodes needs at least one node");
        let mut survivor = self.live_of(first.0 as u32);
        for node in rest {
            let node = self.live_of(node.0 as u32);
            if self.crowns[node as usize].len > self.crowns[survivor as usize].len {
                survivor = node;
            }
        }
        for node in nodes {
            let node = self.live_of(node.0 as u32);
            if node != survivor {
                self.absorb(node, survivor);
            }
        }
        self.data[NodeIndex(survivor as usize)] = node_data_merge;
        NodeIndex(survivor as usize)
    }

    fn is_absorbed(&self, node: NodeIndex) -> bool {
        !self.crowns[node.0].is_live()
    }

    /// Appends a live node with an empty crown. It takes the first free class,
    /// which the class table keeps at the new slot.
    fn add_node(&mut self, node_data: Self::NodeData) -> NodeIndex {
        let node = index_u32(self.crowns.len());
        if self.classes.len() == node as usize {
            self.classes.push(ClassEntry::identity(node));
        }
        self.crowns.push(Crown::EMPTY);
        self.data.push(node_data);
        NodeIndex(node as usize)
    }

    /// Fills the absorbed slots from the tail with [`partition_with`], moving node
    /// records and classes but no half-edge labels, and returns the absorbed
    /// payloads. The absorbed nodes' classes become free.
    fn forget_identification_history(&mut self) -> NodeVec<Self::NodeData> {
        if self.n_absorbed == 0 {
            return NodeVec::new();
        }
        let n = self.crowns.len();
        let live = partition_with(
            self,
            0,
            n,
            |store, node| store.crowns[node].is_live(),
            |store, a, b| store.swap_slots(a, b),
        );
        self.crowns.truncate(live);
        self.n_absorbed = 0;
        self.data.split_off(NodeIndex(live))
    }

    /// Each crown becomes a tree rooted at its smallest half-edge with the other
    /// half-edges as children; absorbed nodes become empty roots.
    fn to_forest<U, H>(
        &self,
        map_data: impl Fn(&Self::NodeData) -> U,
    ) -> Forest<U, ParentPointerStore<H>> {
        let mut nodes: Vec<Option<PPNode<H>>> = std::iter::repeat_with(|| None)
            .take(self.label.len())
            .collect();
        let mut roots = Vec::with_capacity(self.crowns.len());
        for (node, data) in self.data.iter() {
            let crown = &self.crowns[node.0];
            for hedge in self.ring(node.0) {
                nodes[hedge.0] = Some(if hedge.0 as u32 == crown.head {
                    PPNode::dataless_root(RootId(node.0))
                } else {
                    PPNode::dataless_child(TreeNodeId(crown.head as usize))
                });
            }
            roots.push(RootData {
                root_id: if crown.len == 0 {
                    TreeNodeId::EMPTY
                } else {
                    TreeNodeId(crown.head as usize)
                },
                data: map_data(data),
            });
        }
        Forest {
            nodes: nodes
                .into_iter()
                .map(|node| node.expect("crowns partition the half-edges"))
                .collect(),
            roots,
        }
    }

    fn build<I: IntoIterator<Item = HedgeNodeBuilder<Self::NodeData>>>(
        nodes: I,
        n_hedges: usize,
    ) -> Self {
        Self::build_with_mapping(nodes, n_hedges, std::convert::identity)
    }

    fn build_with_mapping<I: IntoIterator<Item = HedgeNodeBuilder<ND>>, ND>(
        nodes: I,
        n_hedges: usize,
        mut map_data: impl FnMut(ND) -> Self::NodeData,
    ) -> Self {
        index_u32(n_hedges);
        let mut label = vec![NONE; n_hedges];
        let nodes = nodes.into_iter();
        let mut data = NodeVec::with_capacity(nodes.size_hint().0);
        for (node, builder) in nodes.enumerate() {
            let node = index_u32(node);
            for hedge in &builder.hedges {
                assert!(
                    label[hedge.0] == NONE || label[hedge.0] == node,
                    "hedge incidence must form a partition"
                );
                label[hedge.0] = node;
            }
            data.push(map_data(builder.data));
        }
        assert!(
            label.iter().all(|&node| node != NONE),
            "hedge incidence must cover all half-edges"
        );
        Self::from_labels(label, data)
    }

    /// Appends the new largest half-edge at the end of the ring of `source`.
    fn add_dangling_edge_mut(&mut self, source: NodeIndex) -> Result<(), HedgeGraphError> {
        let crown = self.crowns.get(source.0).ok_or(HedgeGraphError::NoNode)?;
        if !crown.is_live() {
            return Err(HedgeGraphError::InvalidNode(source));
        }
        let hedge = index_u32(self.label.len());
        self.label.push(self.classes[source.0].class);
        self.links.push(Link::NONE);
        let (mut rings, crowns) = self.rings();
        rings.append(&mut crowns[source.0], hedge);
        Ok(())
    }

    fn random(sources: &[Self::Neighbors], sinks: &[Self::Neighbors]) -> Self
    where
        Self::NodeData: Default,
    {
        let n_hedges = sources.first().map_or(0, |source| source.size());
        let mut label = vec![NONE; n_hedges];
        for (node, crown) in sources.iter().chain(sinks).enumerate() {
            for hedge in crown.included_iter() {
                label[hedge.0] = index_u32(node);
            }
        }
        let data = (0..sources.len() + sinks.len())
            .map(|_| V::default())
            .collect();
        Self::from_labels(label, data)
    }

    fn drain(self) -> impl Iterator<Item = (NodeIndex, Self::NodeData)> {
        self.data.into_iter()
    }

    fn iter(&self) -> impl Iterator<Item = (NodeIndex, &Self::NodeData)> {
        self.data.iter()
    }

    fn check_and_set_nodes(&mut self) -> Result<(), HedgeGraphError> {
        Ok(())
    }

    /// Checks that the class table is a permutation covering the node slots, that
    /// the rings ascend, close and partition the half-edges consistently with the
    /// labels, that absorbed nodes have no crown and lead to live nodes, and that
    /// anchors lie in their crowns.
    fn check_nodes(&self) -> Result<(), HedgeGraphError> {
        let partition = HedgeGraphError::NodesDoNotPartition;
        let n_nodes = self.crowns.len();
        let n_classes = self.classes.len();
        if n_classes < n_nodes
            || (self.classes.iter().enumerate()).any(|(at, entry)| {
                (self.classes.get(entry.class as usize)).map(|entry| entry.slot) != Some(at as u32)
            })
        {
            return Err(partition("class table is not a permutation".to_owned()));
        }
        if (self.label.iter()).any(|&class| {
            class as usize >= n_classes || self.classes[class as usize].slot as usize >= n_nodes
        }) {
            return Err(partition(
                "half-edge labelled with a missing node".to_owned(),
            ));
        }
        let mut covered = 0;
        for (node, crown) in self.crowns.iter().enumerate() {
            if !crown.is_live() {
                if crown.len != 0 || crown.anchor != NONE {
                    return Err(partition(format!("absorbed node {node} has a crown")));
                }
                let mut live = crown.absorbed_into;
                for _ in 0..n_nodes {
                    match self.crowns.get(live as usize) {
                        Some(next) if next.is_live() => break,
                        Some(next) => live = next.absorbed_into,
                        None => return Err(HedgeGraphError::InvalidNode(NodeIndex(node))),
                    }
                }
                if !self.crowns[live as usize].is_live() {
                    return Err(partition(format!(
                        "absorbed node {node} reaches no live node"
                    )));
                }
                continue;
            }
            let mut previous = None;
            let mut hedge = crown.head;
            for _ in 0..crown.len {
                let at = Hedge(hedge as usize);
                if hedge as usize >= self.label.len() {
                    return Err(HedgeGraphError::InvalidHedge(at));
                }
                if self.node_of(hedge as usize) as usize != node {
                    return Err(HedgeGraphError::NodeDoesNotContainHedge(
                        NodeIndex(node),
                        at,
                    ));
                }
                if previous.is_some_and(|previous| previous >= hedge) {
                    return Err(partition(format!(
                        "crown of {node} does not ascend at {at}"
                    )));
                }
                let next = self.links[hedge as usize].next;
                if self.links[next as usize].prev != hedge {
                    return Err(partition(format!("broken ring link after {at}")));
                }
                previous = Some(hedge);
                hedge = next;
            }
            if hedge != crown.head {
                return Err(partition(format!("crown of {node} does not close")));
            }
            covered += crown.len as usize;
            if crown.anchor != NONE
                && (crown.anchor as usize >= self.label.len()
                    || self.node_of(crown.anchor as usize) as usize != node)
            {
                return Err(partition(format!("anchor of {node} is outside its crown")));
            }
        }
        if covered != self.label.len() {
            return Err(partition("crowns do not cover the half-edges".to_owned()));
        }
        if self.count_absorbed() != self.n_absorbed {
            return Err(partition("absorbed node count is stale".to_owned()));
        }
        Ok(())
    }

    fn map_data_ref_graph<'a, E, V2, H>(
        &'a self,
        graph: &'a HedgeGraph<E, Self::NodeData, H, Self>,
        mut node_map: impl FnMut(
            &'a HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> V2,
    ) -> Self::OpStorage<V2> {
        let data = self
            .data
            .iter()
            .map(|(node, data)| node_map(graph, self.ring(node.0), data))
            .collect();
        self.with_data(data)
    }

    fn map_data_ref_mut_graph<'a, V2>(
        &'a mut self,
        mut node_map: impl FnMut(Self::NeighborsIter<'a>, &'a mut Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        let (links, crowns) = (&self.links, &self.crowns);
        let data = self
            .data
            .iter_mut()
            .map(|(node, data)| node_map(crowns[node.0].ring(links), data))
            .collect();
        NodeForest {
            label: self.label.clone(),
            links: self.links.clone(),
            crowns: self.crowns.clone(),
            classes: self.classes.clone(),
            n_absorbed: self.n_absorbed,
            data,
        }
    }

    fn map_data_ref_graph_result<'a, E, V2, H, Er>(
        &'a self,
        graph: &'a HedgeGraph<E, Self::NodeData, H, Self>,
        mut node_map: impl FnMut(
            &'a HedgeGraph<E, Self::NodeData, H, Self>,
            Self::NeighborsIter<'a>,
            &'a Self::NodeData,
        ) -> Result<V2, Er>,
    ) -> Result<Self::OpStorage<V2>, Er> {
        let data = self
            .data
            .iter()
            .map(|(node, data)| node_map(graph, self.ring(node.0), data))
            .collect::<Result<_, Er>>()?;
        Ok(self.with_data(data))
    }

    fn map_data_graph<'a, V2>(
        mut self,
        involution: &'a Involution<EdgeIndex>,
        mut f: impl FnMut(&'a Involution<EdgeIndex>, NodeIndex, Self::NodeData) -> V2,
    ) -> Self::OpStorage<V2> {
        let data = std::mem::replace(&mut self.data, NodeVec::new())
            .into_iter()
            .map(|(node, data)| f(involution, node, data))
            .collect();
        self.into_data(data)
    }

    fn map_data_graph_result<'a, V2, Err>(
        mut self,
        involution: &'a Involution<EdgeIndex>,
        mut f: impl FnMut(&'a Involution<EdgeIndex>, NodeIndex, Self::NodeData) -> Result<V2, Err>,
    ) -> Result<Self::OpStorage<V2>, Err> {
        let data = std::mem::replace(&mut self.data, NodeVec::new())
            .into_iter()
            .map(|(node, data)| f(involution, node, data))
            .collect::<Result<_, Err>>()?;
        Ok(self.into_data(data))
    }

    fn iter_nodes(
        &self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &Self::NodeData)> {
        self.data
            .iter()
            .map(|(node, data)| (node, self.ring(node.0), data))
    }

    fn iter_nodes_mut(
        &mut self,
    ) -> impl Iterator<Item = (NodeIndex, Self::NeighborsIter<'_>, &mut Self::NodeData)> {
        let (links, crowns) = (&self.links, &self.crowns);
        self.data
            .iter_mut()
            .map(move |(node, data)| (node, crowns[node.0].ring(links), data))
    }

    #[inline]
    fn node_id_ref(&self, hedge: Hedge) -> NodeIndex {
        NodeIndex(self.node_of(hedge.0) as usize)
    }

    #[inline]
    fn get_neighbor_iterator(&self, node_id: NodeIndex) -> Self::NeighborsIter<'_> {
        self.ring(node_id.0)
    }

    fn get_node_data(&self, node_id: NodeIndex) -> &Self::NodeData {
        &self.data[node_id]
    }

    fn get_node_data_mut(&mut self, node_id: NodeIndex) -> &mut Self::NodeData {
        &mut self.data[node_id]
    }
}

/// The renumbering of one compaction: the old half-edge and node at each new
/// position, and their inverses.
struct Compaction {
    order: Vec<u32>,
    position: Vec<u32>,
    node_order: Vec<u32>,
    new_id: Vec<u32>,
}

impl Compaction {
    fn new(order: Vec<u32>, node_order: Vec<u32>) -> Self {
        let inverse = |order: &[u32]| {
            let mut inverse = vec![NONE; order.len()];
            for (at, &old) in order.iter().enumerate() {
                inverse[old as usize] = at as u32;
            }
            inverse
        };
        Compaction {
            position: inverse(&order),
            new_id: inverse(&node_order),
            order,
            node_order,
        }
    }

    /// The store holding the half-edges at `hedges` and the nodes at `nodes` of
    /// `store` after renumbering, with payloads `data`. Anchors follow their
    /// half-edge or fall back to the smallest crown half-edge on this side.
    fn side<V, V2>(
        &self,
        store: &NodeForest<V>,
        hedges: Range<usize>,
        nodes: Range<usize>,
        data: NodeVec<V2>,
    ) -> NodeForest<V2> {
        let shift = nodes.start as u32;
        let label = self.order[hedges.clone()]
            .iter()
            .map(|&old| self.new_id[store.node_of(old as usize) as usize] - shift)
            .collect();
        let mut side = NodeForest::from_labels(label, data);
        for (crown, &old) in side.crowns.iter_mut().zip(&self.node_order[nodes]) {
            let anchor = store.crowns[old as usize].anchor;
            if anchor != NONE {
                let at = self.position[anchor as usize] as usize;
                crown.anchor = if hedges.contains(&at) {
                    (at - hedges.start) as u32
                } else {
                    crown.head
                };
            }
        }
        side
    }
}

#[cfg(feature = "serde")]
mod serialization {
    use super::{NodeForest, NONE};

    /// The canonical parts of a [`NodeForest`]; rings are rebuilt on input.
    #[derive(serde::Deserialize)]
    #[serde(rename = "NodeForest")]
    pub(super) struct Parts<V> {
        pub(super) label: Vec<u32>,
        pub(super) anchor: Vec<Option<u32>>,
        pub(super) absorbed_into: Vec<Option<u32>>,
        pub(super) data: Vec<V>,
    }

    #[derive(serde::Serialize)]
    #[serde(rename = "NodeForest")]
    struct PartsRef<'a, V> {
        label: Vec<u32>,
        anchor: Vec<Option<u32>>,
        absorbed_into: Vec<Option<u32>>,
        data: &'a [V],
    }

    fn optional(index: u32) -> Option<u32> {
        (index != NONE).then_some(index)
    }

    /// The inverse of [`optional`], rejecting the reserved index [`NONE`].
    fn stored(index: Option<u32>) -> Result<u32, String> {
        match index {
            Some(NONE) => Err(format!("index {NONE} is reserved")),
            index => Ok(index.unwrap_or(NONE)),
        }
    }

    impl<V> NodeForest<V> {
        /// Each half-edge's node, and each node's anchor and merge target.
        pub(super) fn stored_parts(&self) -> (Vec<u32>, Vec<Option<u32>>, Vec<Option<u32>>) {
            let (anchor, absorbed_into) = self
                .crowns
                .iter()
                .map(|crown| (optional(crown.anchor), optional(crown.absorbed_into)))
                .unzip();
            (self.nodes().collect(), anchor, absorbed_into)
        }

        /// Rebuilds a store from its stored parts, rejecting inconsistent ones.
        pub(super) fn from_stored(parts: Parts<V>) -> Result<Self, String> {
            let indices = |indices: Vec<Option<u32>>| -> Result<Vec<u32>, String> {
                indices.into_iter().map(stored).collect()
            };
            NodeForest::from_parts(
                parts.label,
                indices(parts.anchor)?,
                indices(parts.absorbed_into)?,
                parts.data,
            )
        }
    }

    impl<V: serde::Serialize> serde::Serialize for NodeForest<V> {
        fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            let (label, anchor, absorbed_into) = self.stored_parts();
            PartsRef {
                label,
                anchor,
                absorbed_into,
                data: &self.data.0,
            }
            .serialize(serializer)
        }
    }

    impl<'de, V: serde::Deserialize<'de>> serde::Deserialize<'de> for NodeForest<V> {
        fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            NodeForest::from_stored(Parts::deserialize(deserializer)?)
                .map_err(serde::de::Error::custom)
        }
    }
}

#[cfg(feature = "bincode")]
mod binary {
    use bincode::{
        de::{BorrowDecoder, Decoder},
        enc::Encoder,
        error::{DecodeError, EncodeError},
        BorrowDecode, Decode, Encode,
    };

    use super::{serialization::Parts, NodeForest};

    impl<V: Encode> Encode for NodeForest<V> {
        fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
            let (label, anchor, absorbed_into) = self.stored_parts();
            label.encode(encoder)?;
            anchor.encode(encoder)?;
            absorbed_into.encode(encoder)?;
            self.data.0.encode(encoder)
        }
    }

    fn rebuild<V>(parts: Parts<V>) -> Result<NodeForest<V>, DecodeError> {
        NodeForest::from_stored(parts).map_err(DecodeError::OtherString)
    }

    impl<Context, V: Decode<Context>> Decode<Context> for NodeForest<V> {
        fn decode<D: Decoder<Context = Context>>(decoder: &mut D) -> Result<Self, DecodeError> {
            rebuild(Parts {
                label: Decode::decode(decoder)?,
                anchor: Decode::decode(decoder)?,
                absorbed_into: Decode::decode(decoder)?,
                data: Decode::decode(decoder)?,
            })
        }
    }

    impl<'de, Context, V: BorrowDecode<'de, Context>> BorrowDecode<'de, Context> for NodeForest<V> {
        fn borrow_decode<D: BorrowDecoder<'de, Context = Context>>(
            decoder: &mut D,
        ) -> Result<Self, DecodeError> {
            rebuild(Parts {
                label: BorrowDecode::borrow_decode(decoder)?,
                anchor: BorrowDecode::borrow_decode(decoder)?,
                absorbed_into: BorrowDecode::borrow_decode(decoder)?,
                data: BorrowDecode::borrow_decode(decoder)?,
            })
        }
    }
}

#[cfg(feature = "rkyv")]
mod archive {
    use std::{error::Error, ptr::addr_of};

    use rkyv::{
        bytecheck::{CheckBytes, StructCheckError},
        out_field,
        ser::{ScratchSpace, Serializer},
        validation::ArchiveContext,
        vec::{ArchivedVec, VecResolver},
        Archive, Archived, Deserialize, Fallible, Serialize,
    };

    use super::NodeForest;
    use crate::half_edge::{involution::Hedge, NodeIndex};

    /// The archived canonical image of a [`NodeForest`]: each half-edge's node,
    /// each node's anchor and merge target (`u32::MAX` for none) and the payloads.
    /// Validation and deserialization rebuild the rings, so validation rejects
    /// archives that do not describe a consistent store.
    #[repr(C)]
    pub struct ArchivedNodeForest<V: Archive> {
        label: ArchivedVec<u32>,
        anchor: ArchivedVec<u32>,
        absorbed_into: ArchivedVec<u32>,
        data: ArchivedVec<Archived<V>>,
    }

    impl<V: Archive> ArchivedNodeForest<V> {
        /// The node of `hedge`.
        pub fn node_id(&self, hedge: Hedge) -> NodeIndex {
            NodeIndex(self.label[hedge.0] as usize)
        }

        /// Each node's payload, absorbed nodes included.
        pub fn node_data(&self) -> &[Archived<V>] {
            &self.data
        }

        fn to_store<W>(&self, data: Vec<W>) -> Result<NodeForest<W>, String> {
            NodeForest::from_parts(
                self.label.to_vec(),
                self.anchor.to_vec(),
                self.absorbed_into.to_vec(),
                data,
            )
        }
    }

    pub struct NodeForestResolver {
        label: VecResolver,
        anchor: VecResolver,
        absorbed_into: VecResolver,
        data: VecResolver,
    }

    impl<V: Archive> Archive for NodeForest<V> {
        type Archived = ArchivedNodeForest<V>;
        type Resolver = NodeForestResolver;

        unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
            let n_nodes = self.crowns.len();
            let (fp, fo) = out_field!(out.label);
            ArchivedVec::resolve_from_len(self.label.len(), pos + fp, resolver.label, fo);
            let (fp, fo) = out_field!(out.anchor);
            ArchivedVec::resolve_from_len(n_nodes, pos + fp, resolver.anchor, fo);
            let (fp, fo) = out_field!(out.absorbed_into);
            ArchivedVec::resolve_from_len(n_nodes, pos + fp, resolver.absorbed_into, fo);
            let (fp, fo) = out_field!(out.data);
            ArchivedVec::resolve_from_len(n_nodes, pos + fp, resolver.data, fo);
        }
    }

    impl<V: Serialize<S>, S: ScratchSpace + Serializer + ?Sized> Serialize<S> for NodeForest<V> {
        fn serialize(&self, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
            let (anchor, absorbed_into) = (
                self.crowns.iter().map(|crown| crown.anchor),
                self.crowns.iter().map(|crown| crown.absorbed_into),
            );
            Ok(NodeForestResolver {
                label: ArchivedVec::serialize_from_iter::<u32, _, _, _>(self.nodes(), serializer)?,
                anchor: ArchivedVec::serialize_from_iter::<u32, _, _, _>(anchor, serializer)?,
                absorbed_into: ArchivedVec::serialize_from_iter::<u32, _, _, _>(
                    absorbed_into,
                    serializer,
                )?,
                data: ArchivedVec::serialize_from_slice(&self.data.0, serializer)?,
            })
        }
    }

    impl<V: Archive, D: Fallible + ?Sized> Deserialize<NodeForest<V>, D> for ArchivedNodeForest<V>
    where
        Archived<V>: Deserialize<V, D>,
    {
        /// # Panics
        ///
        /// If the archive was not validated and is inconsistent.
        fn deserialize(&self, deserializer: &mut D) -> Result<NodeForest<V>, D::Error> {
            let data = (self.data.iter())
                .map(|data| data.deserialize(deserializer))
                .collect::<Result<_, _>>()?;
            Ok(self
                .to_store(data)
                .expect("validated archives hold consistent node forests"))
        }
    }

    #[derive(Debug, thiserror::Error)]
    #[error("{0}")]
    struct Inconsistent(String);

    impl<V: Archive, C: ArchiveContext + ?Sized> CheckBytes<C> for ArchivedNodeForest<V>
    where
        Archived<V>: CheckBytes<C>,
        C::Error: Error,
    {
        type Error = StructCheckError;

        unsafe fn check_bytes<'a>(
            value: *const Self,
            context: &mut C,
        ) -> Result<&'a Self, Self::Error> {
            fn field<E: Error + 'static>(
                field_name: &'static str,
            ) -> impl FnOnce(E) -> StructCheckError {
                move |error| StructCheckError {
                    field_name,
                    inner: Box::new(error),
                }
            }
            ArchivedVec::check_bytes(addr_of!((*value).label), context).map_err(field("label"))?;
            ArchivedVec::check_bytes(addr_of!((*value).anchor), context)
                .map_err(field("anchor"))?;
            ArchivedVec::check_bytes(addr_of!((*value).absorbed_into), context)
                .map_err(field("absorbed_into"))?;
            ArchivedVec::check_bytes(addr_of!((*value).data), context).map_err(field("data"))?;
            let value = &*value;
            value
                .to_store(vec![(); value.data.len()])
                .map_err(|error| field("label")(Inconsistent(error)))?;
            Ok(value)
        }
    }

    #[cfg(test)]
    mod tests {
        use super::super::tests::{crowns, store};
        use super::*;
        use crate::half_edge::nodestore::NodeStorageOps;

        #[test]
        fn rkyv_round_trips_and_rejects_inconsistent_archives() {
            let mut store = store(&[&[1, 3], &[0, 4], &[2], &[]]);
            store.set_anchor(NodeIndex(1), Some(Hedge(4)));
            store.identify_nodes(&[NodeIndex(1), NodeIndex(2)], 10);
            let mut bytes = rkyv::to_bytes::<_, 256>(&store).unwrap();
            let archived = rkyv::check_archived_root::<NodeForest<usize>>(&bytes).unwrap();
            for hedge in (0..5).map(Hedge) {
                assert_eq!(archived.node_id(hedge), store.node_id_ref(hedge));
            }
            let data: Vec<_> = archived.node_data().iter().map(|&d| d as usize).collect();
            assert_eq!(data, [0, 10, 2, 3]);
            let label = archived.label.as_ptr() as usize - bytes.as_ptr() as usize;

            let back: NodeForest<usize> = rkyv::from_bytes(&bytes).unwrap();
            assert_eq!(back, store);
            assert_eq!(crowns(&back), crowns(&store));
            assert_eq!(back.anchor(NodeIndex(1)), Some(Hedge(4)));
            assert!(back.is_absorbed(NodeIndex(2)));

            // Half-edge 0 labelled with a missing node passes the byte checks
            // but not the store's.
            bytes[label..label + 4].copy_from_slice(&7u32.to_ne_bytes());
            assert!(rkyv::check_archived_root::<NodeForest<usize>>(&bytes).is_err());
        }
    }
}

#[cfg(feature = "rkyv")]
pub use archive::{ArchivedNodeForest, NodeForestResolver};

#[cfg(test)]
mod tests {
    use std::{
        collections::{hash_map::DefaultHasher, BTreeSet},
        panic::{catch_unwind, AssertUnwindSafe},
    };

    use super::*;
    use crate::half_edge::{subgraph::ModifySubSet, swap::Swap, NodeVec};

    /// A store whose node `i` has payload `i` and the half-edges `parts[i]`.
    pub(super) fn store(parts: &[&[usize]]) -> NodeForest<usize> {
        let n_hedges = parts.iter().flat_map(|part| part.iter()).count();
        NodeForest::build(
            parts
                .iter()
                .enumerate()
                .map(|(data, hedges)| HedgeNodeBuilder {
                    data,
                    hedges: hedges.iter().copied().map(Hedge).collect(),
                }),
            n_hedges,
        )
    }

    fn crown<V>(store: &NodeForest<V>, node: usize) -> Vec<usize> {
        store
            .get_neighbor_iterator(NodeIndex(node))
            .map(|hedge| hedge.0)
            .collect()
    }

    pub(super) fn crowns<V>(store: &NodeForest<V>) -> Vec<Vec<usize>> {
        (0..store.node_bound())
            .map(|node| crown(store, node))
            .collect()
    }

    fn payloads<V: Clone>(data: NodeVec<V>) -> Vec<V> {
        data.into_iter().map(|(_, data)| data).collect()
    }

    #[cfg(feature = "serde")]
    type Graph = crate::half_edge::HedgeGraph<usize, usize, usize, NodeForest<usize>>;

    /// A triangle with a dangling edge and a pending identification.
    #[cfg(feature = "serde")]
    fn graph() -> Graph {
        use crate::half_edge::{builder::HedgeGraphBuilder, involution::Flow};
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<_> = (0..4).map(|i| builder.add_node(i)).collect();
        for i in 0..3 {
            builder.add_edge(
                nodes[i].add_data(2 * i),
                nodes[(i + 1) % 3].add_data(2 * i + 1),
                i,
                true,
            );
        }
        builder.add_external_edge(nodes[3].add_data(6), 3, false, Flow::Sink);
        let mut graph: Graph = builder.build();
        graph.identify_nodes(&[nodes[3], nodes[0]], 10);
        graph
    }

    fn subgraph(n_hedges: usize, hedges: &[usize]) -> SuBitGraph {
        let mut subgraph = SuBitGraph::empty(n_hedges);
        for &hedge in hedges {
            subgraph.add(Hedge(hedge));
        }
        subgraph
    }

    #[test]
    fn build_orders_crowns_and_labels_half_edges() {
        let store = store(&[&[3, 0, 5], &[], &[4, 1, 2]]);
        store.check_nodes().unwrap();
        assert_eq!(crowns(&store), [vec![0, 3, 5], vec![], vec![1, 2, 4]]);
        assert_eq!(store.get_neighbor_iterator(NodeIndex(2)).len(), 3);
        let nodes: Vec<_> = (0..6).map(|h| store.node_id_ref(Hedge(h)).0).collect();
        assert_eq!(nodes, [0, 2, 2, 0, 2, 0]);
        assert_eq!((store.n_live(), store.node_bound()), (3, 3));
        assert_eq!(store.anchor(NodeIndex(0)), None);
    }

    #[test]
    #[should_panic(expected = "partition")]
    fn build_rejects_shared_half_edges() {
        store(&[&[0, 1], &[1]]);
    }

    #[test]
    fn identify_keeps_the_largest_crown_and_renumbers_nothing() {
        let mut store = store(&[&[0, 5], &[1, 2, 3, 6, 9], &[4, 7, 8], &[]]);
        let n = NodeIndex;
        assert_eq!(store.identify_nodes(&[n(0), n(2)], 10), n(2));
        assert_eq!(crown(&store, 2), [0, 4, 5, 7, 8]);
        assert!(crown(&store, 0).is_empty() && store.is_absorbed(n(0)));
        assert_eq!(store.node_id_ref(Hedge(5)), n(2));
        assert_eq!((store.n_live(), store.node_bound()), (3, 4));
        // Absorbed nodes keep their payload; `0` resolves to the five-hedge node 2,
        // which ties with node 1, so the first requested node survives.
        assert_eq!(store.get_node_data(n(0)), &0);
        assert_eq!(store.identify_nodes(&[n(1), n(0), n(1)], 11), n(1));
        assert_eq!(crown(&store, 1), (0..10).collect::<Vec<_>>());
        // A repeated isolated node only takes the new payload.
        assert_eq!(store.identify_nodes(&[n(3), n(3)], 12), n(3));
        assert_eq!(store.get_node_data(n(3)), &12);
        // An absorbed node resolves to the node that absorbed it last.
        assert_eq!(store.identify_nodes(&[n(3), n(2)], 13), n(1));
        assert_eq!((store.n_live(), store.node_bound()), (1, 4));
        store.check_nodes().unwrap();

        let removed = payloads(store.forget_identification_history());
        assert_eq!(
            removed.into_iter().collect::<BTreeSet<_>>(),
            BTreeSet::from([0, 10, 12])
        );
        assert_eq!(crowns(&store), [(0..10).collect::<Vec<_>>()]);
        assert_eq!(store.get_node_data(n(0)), &13);
        store.check_nodes().unwrap();
    }

    #[test]
    fn crowns_stay_ascending_while_merges_interleave() {
        // Node `i` holds every half-edge congruent to `i` modulo 5, and the hub
        // holds the last ones, so absorbed half-edges land inside and below it.
        let parts: Vec<Vec<usize>> = (0..5).map(|i| (i..40).step_by(5).collect()).collect();
        let mut hub: Vec<usize> = (40..60).collect();
        hub.push(2);
        let mut parts: Vec<&[usize]> = parts.iter().map(Vec::as_slice).collect();
        let without_two: Vec<usize> = parts[2].iter().copied().filter(|&h| h != 2).collect();
        parts[2] = &without_two;
        parts.push(&hub);
        let mut store = store(&parts);
        let hub = NodeIndex(5);
        for node in [4, 0, 2, 3, 1] {
            assert_eq!(store.identify_nodes(&[NodeIndex(node), hub], 0), hub);
            store.check_nodes().unwrap();
        }
        assert_eq!(crown(&store, 5), (0..60).collect::<Vec<_>>());
    }

    #[test]
    fn forget_fills_absorbed_slots_from_the_tail() {
        let mut store = store(&[&[0], &[1, 2], &[], &[3], &[]]);
        store.identify_nodes(&[NodeIndex(1), NodeIndex(0)], 10);
        store.identify_nodes(&[NodeIndex(3), NodeIndex(2)], 11);
        assert_eq!(payloads(store.forget_identification_history()), [2, 0]);
        // Slot 0 takes the isolated node 4; the others keep their index.
        let data: Vec<_> = store.iter().map(|(_, data)| *data).collect();
        assert_eq!(data, [4, 10, 11]);
        assert_eq!(crowns(&store), [vec![], vec![0, 1, 2], vec![3]]);
        assert_eq!(store.node_id_ref(Hedge(3)), NodeIndex(2));
        store.check_nodes().unwrap();
        assert!(store.forget_identification_history().0.is_empty());
    }

    /// A hub at the highest index moves into the first hole with its class, so
    /// none of its half-edges is relabelled; the image matches a fresh build.
    #[test]
    fn forget_moves_a_hub_without_relabelling_its_crown() {
        let mut moved = store(&[&[0], &[1], &[2], &[3, 4, 5, 6]]);
        let hub = NodeIndex(3);
        moved.set_anchor(hub, Some(Hedge(5)));
        moved.identify_nodes(&[NodeIndex(0), hub], 10);
        moved.identify_nodes(&[hub, NodeIndex(2)], 11);
        let labels = moved.label.clone();
        assert_eq!(payloads(moved.forget_identification_history()), [2, 0]);
        assert_eq!(moved.label, labels);
        // Classes 2 and 0 of the dropped nodes are free behind the node slots.
        assert_eq!(
            moved
                .classes
                .iter()
                .map(|entry| entry.class)
                .collect::<Vec<_>>(),
            [3, 1, 2, 0]
        );
        assert_eq!((moved.n_live(), moved.node_bound()), (2, 2));
        assert_eq!(crowns(&moved), [vec![0, 2, 3, 4, 5, 6], vec![1]]);
        assert_eq!(moved.node_id_ref(Hedge(4)), NodeIndex(0));
        assert_eq!(moved.anchor(NodeIndex(0)), Some(Hedge(5)));
        moved.check_nodes().unwrap();
        let mut fresh = store(&[&[0, 2, 3, 4, 5, 6], &[1]]);
        *fresh.get_node_data_mut(NodeIndex(0)) = 11;
        assert_eq!(moved, fresh);
        assert_eq!(hash(&moved), hash(&fresh));
        // Serialization stores nodes, not classes.
        #[cfg(feature = "serde")]
        {
            fresh.set_anchor(NodeIndex(0), Some(Hedge(5)));
            assert_eq!(
                serde_json::to_string(&moved).unwrap(),
                serde_json::to_string(&fresh).unwrap()
            );
        }
    }

    /// Appended nodes reuse the classes that forgetting freed, so repeated
    /// identify, forget and append cycles keep the class table bounded, also
    /// when the appended store has moved classes of its own.
    #[test]
    fn classes_are_reused_across_identify_forget_and_append() {
        let mut store = store(&[&[0, 1], &[2], &[3]]);
        for round in 0..20 {
            // Its hub sits last and moves into slot 0 when forgetting.
            let mut other = super::tests::store(&[&[0], &[3], &[1, 2]]);
            other.identify_nodes(&[NodeIndex(0), NodeIndex(2)], 100);
            other.forget_identification_history();
            assert_eq!(other.classes[0], ClassEntry { class: 2, slot: 2 });
            let shift = store.node_bound();
            let n_hedges = store.label.len();
            store.extend_mut(other);
            store.check_nodes().unwrap();
            assert!(
                store.classes.len() <= 5,
                "round {round}: {:?}",
                store.classes
            );
            assert_eq!(crown(&store, shift), [n_hedges, n_hedges + 1, n_hedges + 2]);
            assert_eq!(crown(&store, shift + 1), [n_hedges + 3]);
            // Merge both appended nodes into the node of half-edge 0 and forget.
            let first = store.node_id_ref(Hedge(0));
            store.identify_nodes(&[NodeIndex(shift), first], round);
            store.identify_nodes(&[first, NodeIndex(shift + 1)], round);
            store.forget_identification_history();
            store.check_nodes().unwrap();
            assert_eq!(store.node_bound(), 3);
            let merged = store.node_id_ref(Hedge(0)).0;
            let mut expected = vec![0, 1];
            expected.extend(4..store.label.len());
            assert_eq!(crown(&store, merged), expected);
            assert_eq!(crown(&store, store.node_id_ref(Hedge(2)).0), [2]);
            assert_eq!(crown(&store, store.node_id_ref(Hedge(3)).0), [3]);
        }
    }

    #[test]
    fn anchors_follow_their_half_edge_and_fall_back_to_the_smallest() {
        let mut store = store(&[&[0, 3], &[1, 4], &[2, 5], &[6]]);
        store.set_anchor(NodeIndex(1), Some(Hedge(4)));
        store.set_anchor(NodeIndex(2), Some(Hedge(5)));
        // The survivor without an anchor adopts the first absorbed one.
        store.identify_nodes(&[NodeIndex(0), NodeIndex(2), NodeIndex(1)], 10);
        assert_eq!(store.anchor(NodeIndex(0)), Some(Hedge(5)));
        store.forget_identification_history();
        store.set_anchor(NodeIndex(1), Some(Hedge(6)));
        assert_eq!(crowns(&store), [vec![0, 1, 2, 3, 4, 5], vec![6]]);

        Swap::<Hedge>::swap(&mut store, Hedge(5), Hedge(6));
        assert_eq!(store.anchor(NodeIndex(0)), Some(Hedge(6)));
        assert_eq!(store.anchor(NodeIndex(1)), Some(Hedge(5)));
        assert_eq!(crowns(&store), [vec![0, 1, 2, 3, 4, 6], vec![5]]);

        // Extracting an anchor leaves the smallest remaining half-edge behind.
        let extracted = store.extract(&subgraph(7, &[2, 6]), |d| *d, |d| d);
        let node = |store: &NodeForest<usize>, payload| {
            store.iter().find(|(_, data)| **data == payload).unwrap().0
        };
        assert_eq!(crowns(&store), [vec![2], vec![0, 1, 3, 4]]);
        assert_eq!(store.anchor(node(&store, 10)), Some(Hedge(0)));
        assert_eq!(store.anchor(node(&store, 3)), Some(Hedge(2)));
        assert_eq!(crown(&extracted, 0), [0, 1]);
        assert_eq!(extracted.anchor(NodeIndex(0)), Some(Hedge(1)));
        store.delete(&subgraph(5, &[0]));
        assert_eq!(crowns(&store), [vec![2], vec![0, 1, 3]]);
        assert_eq!(store.anchor(node(&store, 10)), Some(Hedge(0)));
        store.delete(&subgraph(4, &[0, 1, 3]));
        assert_eq!(crowns(&store), [vec![0]]);
        assert_eq!(store.anchor(NodeIndex(0)), Some(Hedge(0)));
        store.check_nodes().unwrap();
        extracted.check_nodes().unwrap();
    }

    #[test]
    #[should_panic(expected = "not in the crown")]
    fn anchors_must_lie_in_their_crown() {
        store(&[&[0], &[1]]).set_anchor(NodeIndex(0), Some(Hedge(1)));
    }

    #[test]
    fn moving_half_edges_or_nodes_requires_forgotten_history() {
        let mut pending = store(&[&[0, 1], &[2], &[3]]);
        pending.identify_nodes(&[NodeIndex(1), NodeIndex(2)], 10);
        let sub = subgraph(4, &[0]);
        for operation in 0..5 {
            let run = |store: &mut NodeForest<usize>| match operation {
                0 => store.delete(&sub),
                1 => drop(store.extract(&sub, |d| *d, |d| d)),
                2 => drop(store.extract_nodes([NodeIndex(0)])),
                3 => Swap::<Hedge>::swap(store, Hedge(0), Hedge(3)),
                _ => Swap::<NodeIndex>::swap(store, NodeIndex(0), NodeIndex(1)),
            };
            let mut store = pending.clone();
            assert!(catch_unwind(AssertUnwindSafe(|| run(&mut store))).is_err());
            let mut store = pending.clone();
            store.forget_identification_history();
            run(&mut store);
            store.check_nodes().unwrap();
        }
    }

    #[test]
    fn dangling_edges_need_a_live_node() {
        let mut store = store(&[&[0], &[1], &[]]);
        store.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 10);
        let before = store.clone();
        assert!(matches!(
            store.add_dangling_edge_mut(NodeIndex(3)),
            Err(HedgeGraphError::NoNode)
        ));
        assert!(matches!(
            store.add_dangling_edge_mut(NodeIndex(1)),
            Err(HedgeGraphError::InvalidNode(NodeIndex(1)))
        ));
        assert_eq!(
            (&store.label, crowns(&store)),
            (&before.label, crowns(&before))
        );
        store.add_dangling_edge_mut(NodeIndex(2)).unwrap();
        store.add_dangling_edge_mut(NodeIndex(0)).unwrap();
        assert_eq!(crowns(&store), [vec![0, 1, 3], vec![], vec![2]]);
        store.check_nodes().unwrap();
    }

    /// Appended nodes take the free classes left by forgetting first, then mint
    /// new ones, and accept half-edges like any other node.
    #[test]
    fn added_nodes_take_free_classes_first() {
        let mut store = store(&[&[0], &[1], &[2]]);
        store.identify_nodes(&[NodeIndex(0), NodeIndex(1), NodeIndex(2)], 10);
        store.forget_identification_history();
        assert_eq!(store.classes.len(), 3);
        let added: Vec<_> = (20..24).map(|data| store.add_node(data)).collect();
        assert_eq!(added, (1..5).map(NodeIndex).collect::<Vec<_>>());
        assert_eq!(store.classes.len(), 5);
        store.check_nodes().unwrap();
        store.add_dangling_edge_mut(added[1]).unwrap();
        store.add_dangling_edge_mut(added[3]).unwrap();
        assert_eq!(
            crowns(&store),
            [vec![0, 1, 2], vec![], vec![3], vec![], vec![4]]
        );
        let data: Vec<_> = store.iter().map(|(_, d)| *d).collect();
        assert_eq!(data, [10, 20, 21, 22, 23]);
        store.check_nodes().unwrap();
    }

    #[test]
    fn appending_shifts_pending_identifications() {
        let mut first = store(&[&[0], &[1]]);
        first.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 10);
        let mut second = store(&[&[0, 2], &[1], &[]]);
        second.set_anchor(NodeIndex(1), Some(Hedge(1)));
        second.identify_nodes(&[NodeIndex(1), NodeIndex(0)], 11);
        first.extend_mut(second);
        first.check_nodes().unwrap();
        assert_eq!((first.n_live(), first.node_bound()), (3, 5));
        assert_eq!(
            crowns(&first),
            [vec![0, 1], vec![], vec![2, 3, 4], vec![], vec![]]
        );
        assert_eq!(first.anchor(NodeIndex(2)), Some(Hedge(3)));
        // The appended absorbed node still resolves to its survivor.
        assert_eq!(
            first.identify_nodes(&[NodeIndex(4), NodeIndex(3)], 12),
            NodeIndex(2)
        );
        assert_eq!(payloads(first.forget_identification_history()).len(), 3);
        first.check_nodes().unwrap();
    }

    /// Removal orders nodes as `[untouched | split | inside]` with the shared
    /// tail-fill kernel, moves the selected half-edges behind the others with the
    /// edge store's swaps and keeps isolated nodes.
    #[test]
    fn removal_renumbers_with_the_shared_kernel() {
        let parts: &[&[usize]] = &[&[0, 7], &[1, 2], &[], &[3, 8], &[4, 5], &[6], &[9]];
        let selection = subgraph(10, &[1, 2, 4, 8, 9]);
        let payloads =
            |store: &NodeForest<usize>| store.iter().map(|(_, d)| *d).collect::<Vec<_>>();
        let mut deleted = store(parts);
        deleted.delete(&selection);
        deleted.check_nodes().unwrap();
        assert_eq!(payloads(&deleted), [0, 5, 2, 3, 4]);
        let nodes: Vec<_> = (0..5)
            .map(|hedge| deleted.node_id_ref(Hedge(hedge)).0)
            .collect();
        assert_eq!(nodes, [0, 0, 1, 3, 4]);

        let mut remaining = store(parts);
        let extracted = remaining.extract(&selection, |d| d + 100, |d| d);
        assert_eq!(payloads(&remaining), payloads(&deleted));
        // Split copies come first, then the nodes wholly inside.
        assert_eq!(payloads(&extracted), [100 + 3, 104, 1, 6]);
        assert_eq!(crowns(&extracted), [vec![3], vec![0], vec![1, 2], vec![4]]);
        extracted.check_nodes().unwrap();
        remaining.check_nodes().unwrap();
    }

    #[test]
    fn extract_nodes_moves_selected_nodes_with_isolated_ones() {
        let mut store = store(&[&[0, 3], &[], &[1, 2], &[4]]);
        let (selected, extracted) = store.extract_nodes([NodeIndex(1), NodeIndex(2)]);
        assert_eq!(
            selected.included_iter().map(|h| h.0).collect::<Vec<_>>(),
            [1, 2]
        );
        let data: Vec<_> = store.iter().map(|(_, d)| *d).collect();
        assert_eq!(data, [0, 3]);
        assert_eq!(crowns(&store), [vec![0, 2], vec![1]]);
        let data: Vec<_> = extracted.iter().map(|(_, d)| *d).collect();
        assert_eq!(data, [2, 1]);
        assert_eq!(crowns(&extracted), [vec![0, 1], vec![]]);
        store.check_nodes().unwrap();
        extracted.check_nodes().unwrap();
    }

    /// Hedge swaps keep crowns ascending and exchange the two half-edges' nodes;
    /// node swaps exchange two nodes' crowns.
    #[test]
    fn hedge_swaps_exchange_nodes_and_keep_crowns_ascending() {
        let parts: &[&[usize]] = &[
            &[0, 5, 9, 13],
            &[1, 2, 14],
            &[3, 4, 6, 7, 8],
            &[10],
            &[11, 12],
        ];
        let mut ring = store(parts);
        // The node of each half-edge.
        let mut model = vec![0; 15];
        for (node, hedges) in parts.iter().enumerate() {
            for &hedge in *hedges {
                model[hedge] = node;
            }
        }
        let expected = |model: &[usize], node: usize| -> Vec<usize> {
            (0..15).filter(|&hedge| model[hedge] == node).collect()
        };
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        for _ in 0..200 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let (a, b) = (state as usize % 15, (state >> 32) as usize % 15);
            Swap::<Hedge>::swap(&mut ring, Hedge(a), Hedge(b));
            model.swap(a, b);
            ring.check_nodes().unwrap();
            for node in 0..5 {
                assert_eq!(crown(&ring, node), expected(&model, node));
            }
        }
        Swap::<NodeIndex>::swap(&mut ring, NodeIndex(0), NodeIndex(4));
        for node in &mut model {
            *node = match *node {
                0 => 4,
                4 => 0,
                node => node,
            };
        }
        for (hedge, &node) in model.iter().enumerate() {
            assert_eq!(ring.node_id_ref(Hedge(hedge)).0, node);
        }
        ring.check_nodes().unwrap();
    }

    fn hash(store: &NodeForest<usize>) -> u64 {
        let mut hasher = DefaultHasher::new();
        store.hash(&mut hasher);
        hasher.finish()
    }

    #[test]
    fn equality_uses_the_canonical_image() {
        let mut a = store(&[&[0, 2], &[1], &[3]]);
        let mut b = a.clone();
        b.set_anchor(NodeIndex(0), Some(Hedge(2)));
        assert_eq!(a, b);
        assert_eq!(hash(&a), hash(&b));
        // The same partition reached through different merges compares equal.
        a.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 7);
        a.identify_nodes(&[NodeIndex(0), NodeIndex(2)], 9);
        b.identify_nodes(&[NodeIndex(0), NodeIndex(2)], 7);
        b.identify_nodes(&[NodeIndex(1), NodeIndex(0)], 9);
        assert_eq!(a, b);
        assert_eq!(hash(&a), hash(&b));
        assert_eq!(a.cmp(&b), std::cmp::Ordering::Equal);
        let c = store(&[&[0, 2], &[1], &[3]]);
        assert_ne!(a, c);
        assert!(a < c);
    }

    #[test]
    fn equality_compares_which_nodes_are_live() {
        // Absorbing node 1, which has no half-edges, keeps labels and payloads.
        let unmerged = store(&[&[0], &[]]);
        let mut merged = unmerged.clone();
        merged.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 0);
        assert!(merged.nodes().eq(unmerged.nodes()));
        assert_eq!(merged.data, unmerged.data);
        assert_ne!(merged, unmerged);
        assert_eq!(merged.cmp(&unmerged), std::cmp::Ordering::Less);
        assert!(merged < unmerged);
        assert_ne!(hash(&merged), hash(&unmerged));
    }

    #[test]
    fn to_forest_roots_each_crown_at_its_smallest_half_edge() {
        let mut store = store(&[&[1, 3], &[0], &[2], &[]]);
        store.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 10);
        let forest = store.to_forest::<usize, ()>(|d| *d);
        forest.validate_structure().unwrap();
        assert_eq!(forest.roots[0].root_id, TreeNodeId(0));
        assert!(forest.roots[1].root_id.is_empty() && forest.roots[3].root_id.is_empty());
        for hedge in 0..4 {
            let node = forest.root(TreeNodeId(hedge));
            assert_eq!(NodeIndex(node.0), store.node_id_ref(Hedge(hedge)));
        }
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serde_round_trips_pending_identifications_and_anchors() {
        let mut store = store(&[&[1, 3], &[0, 4], &[2], &[]]);
        store.set_anchor(NodeIndex(1), Some(Hedge(4)));
        store.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 10);
        let json = serde_json::to_string(&store).unwrap();
        let back: NodeForest<usize> = serde_json::from_str(&json).unwrap();
        back.check_nodes().unwrap();
        assert_eq!(back, store);
        assert_eq!(crowns(&back), crowns(&store));
        assert_eq!(back.anchor(NodeIndex(0)), Some(Hedge(4)));
        assert!(back.is_absorbed(NodeIndex(1)));
        assert_eq!(
            back.clone().identify_nodes(&[NodeIndex(1)], 0),
            NodeIndex(0)
        );
        let graph = graph();
        let json = serde_json::to_string(&graph).unwrap();
        assert_eq!(serde_json::from_str::<Graph>(&json).unwrap(), graph);

        for invalid in [
            // A half-edge of a missing node.
            r#"{"label":[0,5],"anchor":[null],"absorbed_into":[null],"data":[0]}"#,
            // An absorbed node with a crown.
            r#"{"label":[0,1],"anchor":[null,null],"absorbed_into":[null,0],"data":[0,1]}"#,
            // An anchor outside the crown.
            r#"{"label":[0,1],"anchor":[1,null],"absorbed_into":[null,null],"data":[0,1]}"#,
            // Absorbed nodes merged into each other.
            r#"{"label":[],"anchor":[null,null],"absorbed_into":[1,0],"data":[0,1]}"#,
            // The reserved index u32::MAX as an anchor or a merge target.
            r#"{"label":[0],"anchor":[4294967295],"absorbed_into":[null],"data":[0]}"#,
            r#"{"label":[0],"anchor":[null,null],"absorbed_into":[null,4294967295],"data":[0,1]}"#,
        ] {
            assert!(
                serde_json::from_str::<NodeForest<usize>>(invalid).is_err(),
                "{invalid}"
            );
        }
    }

    #[cfg(feature = "bincode")]
    #[test]
    fn bincode_round_trips() {
        let mut store = store(&[&[1, 3], &[0, 4], &[2], &[]]);
        store.set_anchor(NodeIndex(1), Some(Hedge(4)));
        store.identify_nodes(&[NodeIndex(0), NodeIndex(1)], 10);
        let config = bincode::config::standard();
        let bytes = bincode::encode_to_vec(&store, config).unwrap();
        let (back, _): (NodeForest<usize>, _) = bincode::decode_from_slice(&bytes, config).unwrap();
        assert_eq!(back, store);
        assert_eq!(crowns(&back), crowns(&store));
        assert_eq!(back.anchor(NodeIndex(0)), Some(Hedge(4)));
        let graph = graph();
        let bytes = bincode::encode_to_vec(&graph, config).unwrap();
        let (back, _): (Graph, _) = bincode::decode_from_slice(&bytes, config).unwrap();
        assert_eq!(back, graph);
    }

    /// Random identifications of possibly absorbed or repeated nodes, dangling
    /// edges and forgetting, against a union-find model with crown sizes.
    #[test]
    fn random_identifications_match_a_union_find() {
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut random = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        let find = |parent: &[usize], mut node: usize| {
            while parent[node] != node {
                node = parent[node];
            }
            node
        };
        for _ in 0..300 {
            let n_nodes = 1 + random(8);
            let mut owner: Vec<usize> = (0..random(16)).map(|_| random(n_nodes)).collect();
            let mut store = NodeForest::build(
                (0..n_nodes).map(|node| HedgeNodeBuilder {
                    data: node,
                    hedges: (0..owner.len())
                        .filter(|&h| owner[h] == node)
                        .map(Hedge)
                        .collect(),
                }),
                owner.len(),
            );
            let mut parent: Vec<usize> = (0..n_nodes).collect();
            let mut size: Vec<usize> = (0..n_nodes)
                .map(|node| owner.iter().filter(|&&o| o == node).count())
                .collect();
            for step in 0..16 {
                match random(5) {
                    0..=2 => {
                        let requested: Vec<_> =
                            (0..=random(4)).map(|_| random(parent.len())).collect();
                        let mut survivor = find(&parent, requested[0]);
                        for &node in &requested {
                            let node = find(&parent, node);
                            if size[node] > size[survivor] {
                                survivor = node;
                            }
                        }
                        for &node in &requested {
                            let node = find(&parent, node);
                            if node != survivor {
                                size[survivor] += size[node];
                                parent[node] = survivor;
                            }
                        }
                        let requested: Vec<_> = requested.into_iter().map(NodeIndex).collect();
                        assert_eq!(
                            store.identify_nodes(&requested, 100 + step),
                            NodeIndex(survivor)
                        );
                    }
                    3 => {
                        let node = random(parent.len());
                        match store.clone().add_dangling_edge(NodeIndex(node)) {
                            Ok(extended) => {
                                assert_eq!(find(&parent, node), node);
                                store = extended;
                                owner.push(node);
                                size[node] += 1;
                            }
                            Err(HedgeGraphError::InvalidNode(_)) => {
                                assert_ne!(find(&parent, node), node)
                            }
                            Err(error) => panic!("{error}"),
                        }
                    }
                    _ => {
                        let survivors: BTreeSet<_> = (0..parent.len())
                            .filter(|&node| find(&parent, node) == node)
                            .map(|node| *store.get_node_data(NodeIndex(node)))
                            .collect();
                        let expected: Vec<_> = owner
                            .iter()
                            .map(|&node| *store.get_node_data(NodeIndex(find(&parent, node))))
                            .collect();
                        store.forget_identification_history();
                        let nodes = store.iter().map(|(_, data)| *data);
                        assert_eq!(nodes.collect::<BTreeSet<_>>(), survivors);
                        owner = (0..owner.len())
                            .map(|h| store.node_id_ref(Hedge(h)).0)
                            .collect();
                        let actual: Vec<_> = owner
                            .iter()
                            .map(|&node| *store.get_node_data(NodeIndex(node)))
                            .collect();
                        assert_eq!(actual, expected);
                        parent = (0..store.node_bound()).collect();
                        size = (0..store.node_bound())
                            .map(|node| crown(&store, node).len())
                            .collect();
                    }
                }
                store.check_nodes().unwrap();
                for (hedge, &node) in owner.iter().enumerate() {
                    assert_eq!(
                        store.node_id_ref(Hedge(hedge)),
                        NodeIndex(find(&parent, node))
                    );
                }
                let live = (0..parent.len()).filter(|&n| find(&parent, n) == n).count();
                assert_eq!(store.n_live(), live);
            }
        }
    }
}
