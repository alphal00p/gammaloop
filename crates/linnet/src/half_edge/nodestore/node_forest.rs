use std::{
    hash::{Hash, Hasher},
    ops::Range,
};

use super::{Checkpoint, NodeFate, Partition, Topology};
use crate::half_edge::{
    involution::Hedge,
    remap::{index_u32, partition_order, positions, HedgeRemap, NodeRemap},
    HedgeGraphError, NodeIndex,
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

/// A [`Topology`] that keeps each crown as a circular ring of half-edges.
///
/// Every half-edge stores the class of its node and its ring neighbours; a class
/// table maps classes to node slots, so [`Topology::node_of`] costs two loads
/// and moving a node to another slot never relabels its crown. Each ring is
/// ascending from the node's smallest half-edge in every state, also while
/// identifications are pending, so crown order never depends on history. Per node
/// the store keeps its class, the ring head and length, an optional anchor
/// half-edge and whether the node was absorbed by an identification.
///
/// # Classes
///
/// Class ids are stable while their node lives, whatever its slot. Builds,
/// partitions and permutations number classes like slots. The classes of the
/// nodes dropped by [`Topology::compact`] become free and [`Topology::extend`]
/// and [`Topology::push_node`] reuse them before minting new ones, so the class
/// table never outgrows the largest node count since the last rebuild.
///
/// # Identification
///
/// Absorbed half-edges are relabelled immediately (weighted quick-find, so a
/// half-edge is relabelled `O(log H)` times) and spliced into the survivor's
/// ring in ascending position, taking the survivor's class.
///
/// # Compaction
///
/// Compacting moves node records and class table entries only, in `O(1)` per
/// moved node. Partitions and permutations rebuild the rings.
///
/// # Checkpoints
///
/// While a [`Checkpoint`] is open, edits are journaled: identifications record
/// the absorbed crowns, so rolling one back costs its relabelled half-edges.
/// Path compression is suspended while journaling.
///
/// # Equality
///
/// Comparison and hashing use the canonical image: each half-edge's node and
/// which nodes are live. Classes, ring links, anchors, the node an absorbed node
/// was merged into and open checkpoints are not compared.
#[derive(Clone, Debug)]
pub struct NodeForest {
    /// The class of each half-edge's node.
    label: Vec<u32>,
    links: Vec<Link>,
    crowns: Vec<Crown>,
    classes: Vec<ClassEntry>,
    n_absorbed: u32,
    journal: Option<Box<Journal>>,
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

/// The edits recorded since the outermost open checkpoint.
#[derive(Clone, Debug, Default)]
struct Journal {
    /// The log length at each open checkpoint, innermost last.
    marks: Vec<usize>,
    log: Vec<Undo>,
    /// The crowns of the nodes absorbed by logged identifications.
    hedges: Vec<u32>,
}

/// How to undo one journaled edit.
#[derive(Clone, Copy, Debug)]
enum Undo {
    PushNode {
        minted_class: bool,
    },
    Attach {
        node: u32,
    },
    SetAnchor {
        node: u32,
        anchor: u32,
    },
    /// `absorbed`, with crown `crown` whose half-edges are logged from
    /// `hedges` on, merged into `survivor`, whose crown was `survivor_crown`.
    Absorb {
        absorbed: u32,
        crown: Crown,
        survivor: u32,
        survivor_crown: Crown,
        hedges: usize,
    },
}

impl NodeForest {
    /// Builds the store from each half-edge's node, with `n_nodes` live nodes
    /// without anchors and each node's class equal to its slot.
    fn from_label_vec(mut label: Vec<u32>, n_nodes: usize) -> Self {
        let n_hedges = label.len();
        let mut links = vec![Link::NONE; n_hedges];
        let mut crowns = vec![Crown::EMPTY; n_nodes];
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
            classes: (0..index_u32(n_nodes)).map(ClassEntry::identity).collect(),
            n_absorbed: 0,
            journal: None,
        }
    }

    #[inline]
    fn ring(&self, node: usize) -> RingNeighborIter<'_> {
        self.crowns[node].ring(&self.links)
    }

    /// The node slot of `hedge`.
    #[inline]
    fn slot_of(&self, hedge: usize) -> u32 {
        let class = self.label[hedge] as usize;
        debug_assert!(class < self.classes.len());
        // SAFETY: labels enter through `from_label_vec`, which indexes the node
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

    /// The live node that `node` was, possibly transitively, identified into.
    fn find(&self, mut node: u32) -> u32 {
        while self.crowns[node as usize].absorbed_into != NONE {
            node = self.crowns[node as usize].absorbed_into;
        }
        node
    }

    /// Like [`Self::find`], compressing the path unless edits are journaled.
    fn live_of(&mut self, node: u32) -> u32 {
        let live = self.find(node);
        if self.journal.is_none() {
            let mut node = node;
            while node != live {
                let next = self.crowns[node as usize].absorbed_into;
                self.crowns[node as usize].absorbed_into = live;
                node = next;
            }
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
        if let Some(journal) = &mut self.journal {
            let undo = Undo::Absorb {
                absorbed,
                crown,
                survivor,
                survivor_crown: self.crowns[survivor as usize],
                hedges: journal.hedges.len(),
            };
            journal.log.push(undo);
            journal
                .hedges
                .extend(crown.ring(&self.links).map(|h| h.0 as u32));
        }
        self.n_absorbed += 1;
        let class = self.classes[survivor as usize].class;
        let (mut rings, crowns) = self.rings();
        let target = &mut crowns[survivor as usize];
        rings.merge(target, class, &crown);
        if target.anchor == NONE {
            target.anchor = crown.anchor;
        }
    }

    /// Undoes one journaled edit; `hedges` holds the crowns logged so far.
    fn undo(&mut self, undo: Undo, hedges: &mut Vec<u32>) {
        match undo {
            Undo::PushNode { minted_class } => {
                self.crowns.pop();
                if minted_class {
                    self.classes.pop();
                }
            }
            Undo::Attach { node } => {
                let hedge = self.label.len() as u32 - 1;
                let (mut rings, crowns) = self.rings();
                rings.unlink(&mut crowns[node as usize], hedge);
                self.label.pop();
                self.links.pop();
            }
            Undo::SetAnchor { node, anchor } => self.crowns[node as usize].anchor = anchor,
            Undo::Absorb {
                absorbed,
                crown,
                survivor,
                survivor_crown,
                hedges: start,
            } => {
                let class = self.classes[absorbed as usize].class;
                let (mut rings, crowns) = self.rings();
                for &hedge in &hedges[start..] {
                    rings.unlink(&mut crowns[survivor as usize], hedge);
                }
                debug_assert_eq!(crowns[survivor as usize].len, survivor_crown.len);
                crowns[survivor as usize] = survivor_crown;
                let restored = &mut crowns[absorbed as usize];
                *restored = Crown {
                    anchor: crown.anchor,
                    ..Crown::EMPTY
                };
                for &hedge in &hedges[start..] {
                    rings.label[hedge as usize] = class;
                    rings.append(restored, hedge);
                }
                debug_assert_eq!((restored.head, restored.len), (crown.head, crown.len));
                hedges.truncate(start);
                self.n_absorbed -= 1;
            }
        }
    }

    /// Journals an edit while a checkpoint is open.
    #[inline]
    fn record(&mut self, undo: impl FnOnce() -> Undo) {
        if let Some(journal) = &mut self.journal {
            journal.log.push(undo());
        }
    }

    #[cfg(any(feature = "serde", feature = "bincode", feature = "rkyv", test))]
    fn count_absorbed(&self) -> u32 {
        self.crowns.iter().filter(|crown| !crown.is_live()).count() as u32
    }

    /// Fails while a checkpoint is open.
    fn unjournaled(&self) -> Result<(), HedgeGraphError> {
        match self.journal {
            Some(_) => Err(HedgeGraphError::CheckpointOpen),
            None => Ok(()),
        }
    }

    /// Fails while a checkpoint is open or identifications are pending: absorbed
    /// nodes point at their survivors by index, so half-edges and nodes must not
    /// move.
    fn settled(&self) -> Result<(), HedgeGraphError> {
        self.unjournaled()?;
        match self.n_absorbed {
            0 => Ok(()),
            _ => Err(HedgeGraphError::PendingIdentifications),
        }
    }

    /// Swaps two node slots with their classes, crowns, anchors and states in
    /// `O(1)`: half-edges keep their class.
    fn swap_slots(&mut self, a: usize, b: usize) {
        self.crowns.swap(a, b);
        let (class_a, class_b) = (self.classes[a].class, self.classes[b].class);
        self.classes[a].class = class_b;
        self.classes[b].class = class_a;
        self.classes[class_a as usize].slot = b as u32;
        self.classes[class_b as usize].slot = a as u32;
    }

    /// Whether node `i` is live, for each node: part of the canonical image.
    fn live_mask(&self) -> impl Iterator<Item = bool> + '_ {
        self.crowns.iter().map(Crown::is_live)
    }

    /// Checks the representation behind the [`Topology`] contract: the class
    /// table is a permutation covering the node slots, labels resolve to nodes,
    /// rings close and match the labels, absorbed nodes lead to live ones and
    /// the absorbed count is current.
    #[cfg(any(feature = "serde", feature = "bincode", feature = "rkyv", test))]
    fn check_rings(&self) -> Result<(), HedgeGraphError> {
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
        for (node, crown) in self.crowns.iter().enumerate() {
            if !crown.is_live() {
                // A chain of merges longer than the node count is a cycle.
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
            let mut hedge = crown.head;
            for _ in 0..crown.len {
                if hedge as usize >= self.label.len() {
                    return Err(HedgeGraphError::InvalidHedge(Hedge(hedge as usize)));
                }
                let next = self.links[hedge as usize].next;
                if self.links[next as usize].prev != hedge {
                    return Err(partition(format!("broken ring link after {hedge}")));
                }
                hedge = next;
            }
            if hedge != crown.head {
                return Err(partition(format!("crown of {node} does not close")));
            }
        }
        if self.count_absorbed() != self.n_absorbed {
            return Err(partition("absorbed node count is stale".to_owned()));
        }
        Ok(())
    }

    /// Rebuilds a store from its canonical parts, the input of every
    /// deserializer: each half-edge's node, and each node's anchor and merge
    /// target, [`NONE`] for none. Rejects inconsistent parts.
    #[cfg(any(feature = "serde", feature = "bincode", feature = "rkyv"))]
    fn from_parts(
        label: Vec<u32>,
        anchor: Vec<u32>,
        absorbed_into: Vec<u32>,
    ) -> Result<Self, String> {
        let n_nodes = anchor.len();
        if absorbed_into.len() != n_nodes {
            return Err("one anchor and state per node expected".to_owned());
        }
        if label.iter().any(|&node| node as usize >= n_nodes) {
            return Err("half-edge labelled with a missing node".to_owned());
        }
        let mut store = NodeForest::from_label_vec(label, n_nodes);
        for (crown, (anchor, absorbed_into)) in store
            .crowns
            .iter_mut()
            .zip(anchor.into_iter().zip(absorbed_into))
        {
            crown.anchor = anchor;
            crown.absorbed_into = absorbed_into;
        }
        store.n_absorbed = store.count_absorbed();
        (store.check_rings())
            .and_then(|()| store.check())
            .map_err(|error| error.to_string())?;
        Ok(store)
    }
}

impl PartialEq for NodeForest {
    fn eq(&self, other: &Self) -> bool {
        self.nodes().eq(other.nodes()) && self.live_mask().eq(other.live_mask())
    }
}

impl Eq for NodeForest {}

impl PartialOrd for NodeForest {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NodeForest {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.nodes()
            .cmp(other.nodes())
            .then_with(|| self.live_mask().cmp(other.live_mask()))
    }
}

impl Hash for NodeForest {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.label.len().hash(state);
        for node in self.nodes() {
            node.hash(state);
        }
        for live in self.live_mask() {
            live.hash(state);
        }
    }
}

impl Topology for NodeForest {
    type Crown<'a> = RingNeighborIter<'a>;

    fn n_hedges(&self) -> usize {
        self.label.len()
    }

    fn node_bound(&self) -> usize {
        self.crowns.len()
    }

    fn n_live(&self) -> usize {
        self.crowns.len() - self.n_absorbed as usize
    }

    #[inline]
    fn node_of(&self, hedge: Hedge) -> NodeIndex {
        NodeIndex(self.slot_of(hedge.0) as usize)
    }

    #[inline]
    fn crown(&self, node: NodeIndex) -> RingNeighborIter<'_> {
        self.ring(node.0)
    }

    fn is_live(&self, node: NodeIndex) -> bool {
        self.crowns[node.0].is_live()
    }

    fn representative(&self, node: NodeIndex) -> NodeIndex {
        NodeIndex(self.find(node.0 as u32) as usize)
    }

    fn anchor(&self, node: NodeIndex) -> Option<Hedge> {
        let anchor = self.crowns[node.0].anchor;
        (anchor != NONE).then_some(Hedge(anchor as usize))
    }

    /// Sets or clears the anchor of `node` in `O(1)`.
    fn set_anchor(&mut self, node: NodeIndex, anchor: Option<Hedge>) {
        let anchor = match anchor {
            Some(hedge) => {
                assert_eq!(
                    self.slot_of(hedge.0) as usize,
                    node.0,
                    "anchor {hedge} is not in the crown of {node}"
                );
                hedge.0 as u32
            }
            None => NONE,
        };
        let previous = std::mem::replace(&mut self.crowns[node.0].anchor, anchor);
        self.record(|| Undo::SetAnchor {
            node: node.0 as u32,
            anchor: previous,
        });
    }

    fn from_labels(labels: impl ExactSizeIterator<Item = NodeIndex>, node_bound: usize) -> Self {
        index_u32(labels.len());
        index_u32(node_bound);
        let label = labels
            .map(|node| {
                assert!(
                    node.0 < node_bound,
                    "half-edge labelled with a missing node"
                );
                node.0 as u32
            })
            .collect();
        Self::from_label_vec(label, node_bound)
    }

    /// Appends a live node with an empty crown. It takes the first free class,
    /// which the class table keeps at the new slot.
    fn push_node(&mut self) -> NodeIndex {
        let node = index_u32(self.crowns.len());
        let minted_class = self.classes.len() == node as usize;
        if minted_class {
            self.classes.push(ClassEntry::identity(node));
        }
        self.crowns.push(Crown::EMPTY);
        self.record(|| Undo::PushNode { minted_class });
        NodeIndex(node as usize)
    }

    /// Appends the new largest half-edge at the end of the ring of `node`.
    fn attach(&mut self, node: NodeIndex) -> Result<Hedge, HedgeGraphError> {
        let crown = self.crowns.get(node.0).ok_or(HedgeGraphError::NoNode)?;
        if !crown.is_live() {
            return Err(HedgeGraphError::InvalidNode(node));
        }
        let hedge = index_u32(self.label.len());
        self.label.push(self.classes[node.0].class);
        self.links.push(Link::NONE);
        let (mut rings, crowns) = self.rings();
        rings.append(&mut crowns[node.0], hedge);
        self.record(|| Undo::Attach {
            node: node.0 as u32,
        });
        Ok(Hedge(hedge as usize))
    }

    /// Appended nodes take the free classes of `self` first.
    fn extend(&mut self, other: Self) -> Result<(), HedgeGraphError> {
        self.unjournaled()?;
        other.unjournaled()?;
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
        Ok(())
    }

    /// Costs the relabelled half-edges plus their ring insertions, each bounded
    /// by the survivor's crown members or half-edge indices it skips.
    fn identify(&mut self, nodes: &[NodeIndex]) -> NodeIndex {
        let (first, rest) = nodes
            .split_first()
            .expect("identify needs at least one node");
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
        NodeIndex(survivor as usize)
    }

    /// Fills the absorbed slots from the tail with the shared kernel, moving node
    /// records and classes but no half-edge labels. The absorbed nodes' classes
    /// become free.
    fn compact(&mut self) -> Result<NodeRemap, HedgeGraphError> {
        self.unjournaled()?;
        let n = self.crowns.len();
        if self.n_absorbed == 0 {
            return Ok(NodeRemap::identity(n));
        }
        let remap = NodeRemap::compaction(
            self,
            n,
            |store, node| store.crowns[node].is_live(),
            |store, node| store.crowns[node].absorbed_into as usize,
            |store, hole, live| store.swap_slots(hole, live),
        );
        self.crowns.truncate(remap.new_bound());
        self.n_absorbed = 0;
        Ok(remap)
    }

    fn partition_hedges(
        &mut self,
        moved: impl Fn(Hedge) -> bool,
        fate: impl Fn(NodeIndex, usize, usize) -> NodeFate,
        build_moved: bool,
    ) -> Result<Partition<Self>, HedgeGraphError> {
        self.settled()?;
        let (n_hedges, n_nodes) = (self.label.len(), self.crowns.len());
        let mut order: Vec<u32> = (0..n_hedges as u32).collect();
        let n_staying = partition_order(&mut order, 0, |hedge| !moved(Hedge(hedge as usize)));
        let mut n_moved = vec![0; n_nodes];
        for &hedge in &order[n_staying..] {
            n_moved[self.slot_of(hedge as usize) as usize] += 1;
        }
        let fates: Vec<NodeFate> = (self.crowns.iter().zip(&n_moved).enumerate())
            .map(|(node, (crown, &n_moved))| {
                let fate = fate(
                    NodeIndex(node),
                    (crown.len - n_moved) as usize,
                    n_moved as usize,
                );
                debug_assert!(
                    (fate != NodeFate::Kept || n_moved == 0)
                        && (fate != NodeFate::Moved || n_moved == crown.len),
                    "node {node} cannot keep all its half-edges on one side"
                );
                fate
            })
            .collect();
        let mut node_order: Vec<u32> = (0..n_nodes as u32).collect();
        let split = partition_order(&mut node_order, 0, |node| {
            fates[node as usize] == NodeFate::Kept
        });
        let inside = partition_order(&mut node_order, split, |node| {
            fates[node as usize] == NodeFate::Split
        });
        let compaction = Compaction::new(order, node_order);
        let moved_side =
            build_moved.then(|| compaction.side(self, n_staying..n_hedges, split..n_nodes));
        *self = compaction.side(self, 0..n_staying, 0..inside);
        Ok(Partition {
            moved: moved_side,
            hedges: HedgeRemap::bijection(compaction.position),
            n_staying,
            nodes: NodeRemap::bijection(compaction.new_id),
            split,
            inside,
        })
    }

    fn permute(&mut self, hedges: &HedgeRemap, nodes: &NodeRemap) -> Result<(), HedgeGraphError> {
        self.settled()?;
        let (n_hedges, n_nodes) = (self.label.len(), self.crowns.len());
        assert!(
            (
                hedges.old_bound(),
                hedges.new_bound(),
                nodes.old_bound(),
                nodes.new_bound()
            ) == (n_hedges, n_hedges, n_nodes, n_nodes),
            "permutations must cover every half-edge and node"
        );
        let new_hedge = |hedge: usize| hedges.get(Hedge(hedge)).expect("a bijection").0;
        let new_node = |node: usize| nodes.get(NodeIndex(node)).expect("a bijection").0;
        let mut label = vec![NONE; n_hedges];
        for hedge in 0..n_hedges {
            label[new_hedge(hedge)] = new_node(self.slot_of(hedge) as usize) as u32;
        }
        let mut permuted = Self::from_label_vec(label, n_nodes);
        for (node, crown) in self.crowns.iter().enumerate() {
            if crown.anchor != NONE {
                permuted.crowns[new_node(node)].anchor = new_hedge(crown.anchor as usize) as u32;
            }
        }
        *self = permuted;
        Ok(())
    }
}

impl Checkpoint for NodeForest {
    fn checkpoint(&mut self) {
        let journal = self.journal.get_or_insert_default();
        journal.marks.push(journal.log.len());
    }

    fn rollback(&mut self) {
        let mut journal = self.journal.take().expect("no checkpoint is open");
        let mark = journal.marks.pop().expect("no checkpoint is open");
        while journal.log.len() > mark {
            let undo = journal.log.pop().expect("logged edits");
            self.undo(undo, &mut journal.hedges);
        }
        if !journal.marks.is_empty() {
            self.journal = Some(journal);
        }
    }

    fn commit(&mut self) {
        let journal = self.journal.as_mut().expect("no checkpoint is open");
        journal.marks.pop();
        if journal.marks.is_empty() {
            self.journal = None;
        }
    }

    fn has_checkpoint(&self) -> bool {
        self.journal.is_some()
    }
}

/// The renumbering of one partition: the old half-edge and node at each new
/// position, and their inverses.
struct Compaction {
    order: Vec<u32>,
    position: Vec<u32>,
    node_order: Vec<u32>,
    new_id: Vec<u32>,
}

impl Compaction {
    fn new(order: Vec<u32>, node_order: Vec<u32>) -> Self {
        Compaction {
            position: positions(&order),
            new_id: positions(&node_order),
            order,
            node_order,
        }
    }

    /// The store holding the half-edges at `hedges` and the nodes at `nodes` of
    /// `store` after renumbering. Anchors follow their half-edge or fall back to
    /// the smallest crown half-edge on this side.
    fn side(&self, store: &NodeForest, hedges: Range<usize>, nodes: Range<usize>) -> NodeForest {
        let shift = nodes.start as u32;
        let label = self.order[hedges.clone()]
            .iter()
            .map(|&old| self.new_id[store.slot_of(old as usize) as usize] - shift)
            .collect();
        let mut side = NodeForest::from_label_vec(label, nodes.len());
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

#[cfg(any(feature = "serde", feature = "bincode"))]
mod serialization {
    use super::{NodeForest, NONE};

    /// The canonical parts of a [`NodeForest`]; rings are rebuilt on input.
    #[cfg_attr(feature = "serde", derive(serde::Deserialize))]
    #[cfg_attr(feature = "serde", serde(rename = "NodeForest"))]
    pub(super) struct Parts {
        pub(super) label: Vec<u32>,
        pub(super) anchor: Vec<Option<u32>>,
        pub(super) absorbed_into: Vec<Option<u32>>,
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

    impl NodeForest {
        /// Each half-edge's node, and each node's anchor and merge target.
        pub(super) fn stored_parts(&self) -> Parts {
            let (anchor, absorbed_into) = self
                .crowns
                .iter()
                .map(|crown| (optional(crown.anchor), optional(crown.absorbed_into)))
                .unzip();
            Parts {
                label: self.nodes().collect(),
                anchor,
                absorbed_into,
            }
        }

        /// Rebuilds a store from its stored parts, rejecting inconsistent ones.
        pub(super) fn from_stored(parts: Parts) -> Result<Self, String> {
            let indices = |indices: Vec<Option<u32>>| -> Result<Vec<u32>, String> {
                indices.into_iter().map(stored).collect()
            };
            NodeForest::from_parts(
                parts.label,
                indices(parts.anchor)?,
                indices(parts.absorbed_into)?,
            )
        }
    }

    #[cfg(feature = "serde")]
    impl serde::Serialize for NodeForest {
        fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            use serde::ser::SerializeStruct;
            let parts = self.stored_parts();
            let mut state = serializer.serialize_struct("NodeForest", 3)?;
            state.serialize_field("label", &parts.label)?;
            state.serialize_field("anchor", &parts.anchor)?;
            state.serialize_field("absorbed_into", &parts.absorbed_into)?;
            state.end()
        }
    }

    #[cfg(feature = "serde")]
    impl<'de> serde::Deserialize<'de> for NodeForest {
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

    impl Encode for NodeForest {
        fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
            let parts = self.stored_parts();
            parts.label.encode(encoder)?;
            parts.anchor.encode(encoder)?;
            parts.absorbed_into.encode(encoder)
        }
    }

    fn rebuild(parts: Parts) -> Result<NodeForest, DecodeError> {
        NodeForest::from_stored(parts).map_err(DecodeError::OtherString)
    }

    impl<Context> Decode<Context> for NodeForest {
        fn decode<D: Decoder<Context = Context>>(decoder: &mut D) -> Result<Self, DecodeError> {
            rebuild(Parts {
                label: Decode::decode(decoder)?,
                anchor: Decode::decode(decoder)?,
                absorbed_into: Decode::decode(decoder)?,
            })
        }
    }

    impl<'de, Context> BorrowDecode<'de, Context> for NodeForest {
        fn borrow_decode<D: BorrowDecoder<'de, Context = Context>>(
            decoder: &mut D,
        ) -> Result<Self, DecodeError> {
            rebuild(Parts {
                label: BorrowDecode::borrow_decode(decoder)?,
                anchor: BorrowDecode::borrow_decode(decoder)?,
                absorbed_into: BorrowDecode::borrow_decode(decoder)?,
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
    use crate::half_edge::{involution::Hedge, nodestore::NodeStore, NodeIndex, NodeVec};

    /// The archived canonical image of a [`NodeForest`]: each half-edge's node
    /// and each node's anchor and merge target (`u32::MAX` for none).
    /// Validation and deserialization rebuild the rings, so validation rejects
    /// archives that do not describe a consistent store.
    #[repr(C)]
    pub struct ArchivedNodeForest {
        label: ArchivedVec<u32>,
        anchor: ArchivedVec<u32>,
        absorbed_into: ArchivedVec<u32>,
    }

    impl ArchivedNodeForest {
        /// The node of `hedge`.
        pub fn node_id(&self, hedge: Hedge) -> NodeIndex {
            NodeIndex(self.label[hedge.0] as usize)
        }

        /// The number of node slots.
        pub fn node_bound(&self) -> usize {
            self.anchor.len()
        }

        fn to_store(&self) -> Result<NodeForest, String> {
            NodeForest::from_parts(
                self.label.to_vec(),
                self.anchor.to_vec(),
                self.absorbed_into.to_vec(),
            )
        }
    }

    pub struct NodeForestResolver {
        label: VecResolver,
        anchor: VecResolver,
        absorbed_into: VecResolver,
    }

    impl Archive for NodeForest {
        type Archived = ArchivedNodeForest;
        type Resolver = NodeForestResolver;

        unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
            let n_nodes = self.crowns.len();
            let (fp, fo) = out_field!(out.label);
            ArchivedVec::resolve_from_len(self.label.len(), pos + fp, resolver.label, fo);
            let (fp, fo) = out_field!(out.anchor);
            ArchivedVec::resolve_from_len(n_nodes, pos + fp, resolver.anchor, fo);
            let (fp, fo) = out_field!(out.absorbed_into);
            ArchivedVec::resolve_from_len(n_nodes, pos + fp, resolver.absorbed_into, fo);
        }
    }

    impl<S: ScratchSpace + Serializer + ?Sized> Serialize<S> for NodeForest {
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
            })
        }
    }

    impl<D: Fallible + ?Sized> Deserialize<NodeForest, D> for ArchivedNodeForest {
        /// # Panics
        ///
        /// If the archive was not validated and is inconsistent.
        fn deserialize(&self, _: &mut D) -> Result<NodeForest, D::Error> {
            Ok(self
                .to_store()
                .expect("validated archives hold consistent node forests"))
        }
    }

    #[derive(Debug, thiserror::Error)]
    #[error("{0}")]
    struct Inconsistent(String);

    fn field<E: Error + 'static>(field_name: &'static str) -> impl FnOnce(E) -> StructCheckError {
        move |error| StructCheckError {
            field_name,
            inner: Box::new(error),
        }
    }

    impl<C: ArchiveContext + ?Sized> CheckBytes<C> for ArchivedNodeForest
    where
        C::Error: Error,
    {
        type Error = StructCheckError;

        unsafe fn check_bytes<'a>(
            value: *const Self,
            context: &mut C,
        ) -> Result<&'a Self, Self::Error> {
            ArchivedVec::check_bytes(addr_of!((*value).label), context).map_err(field("label"))?;
            ArchivedVec::check_bytes(addr_of!((*value).anchor), context)
                .map_err(field("anchor"))?;
            ArchivedVec::check_bytes(addr_of!((*value).absorbed_into), context)
                .map_err(field("absorbed_into"))?;
            let value = &*value;
            value
                .to_store()
                .map_err(|error| field("label")(Inconsistent(error)))?;
            Ok(value)
        }
    }

    /// The archive of a [`NodeStore`] over a [`NodeForest`]: the topology and
    /// one payload per node slot, which validation checks.
    #[repr(C)]
    pub struct ArchivedNodeStore<V: Archive> {
        topology: ArchivedNodeForest,
        data: ArchivedVec<Archived<V>>,
    }

    impl<V: Archive> ArchivedNodeStore<V> {
        /// The node of `hedge`.
        pub fn node_id(&self, hedge: Hedge) -> NodeIndex {
            self.topology.node_id(hedge)
        }

        /// Each node slot's payload, absorbed nodes included.
        pub fn node_data(&self) -> &[Archived<V>] {
            &self.data
        }
    }

    pub struct NodeStoreResolver {
        topology: NodeForestResolver,
        data: VecResolver,
    }

    impl<V: Archive> Archive for NodeStore<NodeForest, V> {
        type Archived = ArchivedNodeStore<V>;
        type Resolver = NodeStoreResolver;

        unsafe fn resolve(&self, pos: usize, resolver: Self::Resolver, out: *mut Self::Archived) {
            let (fp, fo) = out_field!(out.topology);
            self.topology.resolve(pos + fp, resolver.topology, fo);
            let (fp, fo) = out_field!(out.data);
            ArchivedVec::resolve_from_len(self.data.0.len(), pos + fp, resolver.data, fo);
        }
    }

    impl<V: Serialize<S>, S: ScratchSpace + Serializer + ?Sized> Serialize<S>
        for NodeStore<NodeForest, V>
    {
        fn serialize(&self, serializer: &mut S) -> Result<Self::Resolver, S::Error> {
            Ok(NodeStoreResolver {
                topology: self.topology.serialize(serializer)?,
                data: ArchivedVec::serialize_from_slice(&self.data.0, serializer)?,
            })
        }
    }

    impl<V: Archive, D: Fallible + ?Sized> Deserialize<NodeStore<NodeForest, V>, D>
        for ArchivedNodeStore<V>
    where
        Archived<V>: Deserialize<V, D>,
    {
        fn deserialize(&self, deserializer: &mut D) -> Result<NodeStore<NodeForest, V>, D::Error> {
            let data = (self.data.iter())
                .map(|data| data.deserialize(deserializer))
                .collect::<Result<_, _>>()?;
            Ok(NodeStore::new(
                self.topology.deserialize(deserializer)?,
                NodeVec(data),
            ))
        }
    }

    impl<V: Archive, C: ArchiveContext + ?Sized> CheckBytes<C> for ArchivedNodeStore<V>
    where
        Archived<V>: CheckBytes<C>,
        C::Error: Error,
    {
        type Error = StructCheckError;

        unsafe fn check_bytes<'a>(
            value: *const Self,
            context: &mut C,
        ) -> Result<&'a Self, Self::Error> {
            ArchivedNodeForest::check_bytes(addr_of!((*value).topology), context)
                .map_err(field("topology"))?;
            ArchivedVec::check_bytes(addr_of!((*value).data), context).map_err(field("data"))?;
            let value = &*value;
            if value.data.len() != value.topology.node_bound() {
                return Err(field("data")(Inconsistent(
                    "one payload per node slot expected".to_owned(),
                )));
            }
            Ok(value)
        }
    }

    #[cfg(test)]
    mod tests {
        use super::super::tests::{crowns, store};
        use super::*;
        use crate::half_edge::nodestore::Topology;

        #[test]
        fn rkyv_round_trips_and_rejects_inconsistent_archives() {
            let mut store = store(&[&[1, 3], &[0, 4], &[2], &[]]);
            store.set_anchor(NodeIndex(1), Some(Hedge(4)));
            store.identify(&[NodeIndex(1), NodeIndex(2)], 10);
            let mut bytes = rkyv::to_bytes::<_, 256>(&store).unwrap();
            let archived =
                rkyv::check_archived_root::<NodeStore<NodeForest, usize>>(&bytes).unwrap();
            for hedge in (0..5).map(Hedge) {
                assert_eq!(archived.node_id(hedge), store.node_id(hedge));
            }
            let data: Vec<_> = archived.node_data().iter().map(|&d| d as usize).collect();
            assert_eq!(data, [0, 10, 2, 3]);
            let label = archived.topology.label.as_ptr() as usize - bytes.as_ptr() as usize;

            let back: NodeStore<NodeForest, usize> = rkyv::from_bytes(&bytes).unwrap();
            assert_eq!(back, store);
            assert_eq!(crowns(&back), crowns(&store));
            assert_eq!(back.anchor(NodeIndex(1)), Some(Hedge(4)));
            assert!(!back.topology().is_live(NodeIndex(2)));

            // Half-edge 0 labelled with a missing node passes the byte checks
            // but not the store's.
            bytes[label..label + 4].copy_from_slice(&7u32.to_ne_bytes());
            assert!(rkyv::check_archived_root::<NodeStore<NodeForest, usize>>(&bytes).is_err());
        }
    }
}

#[cfg(feature = "rkyv")]
pub use archive::{ArchivedNodeForest, ArchivedNodeStore};

#[cfg(test)]
mod tests {
    use std::collections::{hash_map::DefaultHasher, BTreeSet};

    use super::*;
    use crate::{
        half_edge::{
            builder::HedgeNodeBuilder,
            nodestore::NodeStore,
            subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
        },
        permutation::Permutation,
        tree::TreeNodeId,
    };

    type Store = NodeStore<NodeForest, usize>;

    /// A store whose node `i` has payload `i` and the half-edges `parts[i]`.
    pub(super) fn store(parts: &[&[usize]]) -> Store {
        let n_hedges = parts.iter().flat_map(|part| part.iter()).count();
        NodeStore::build(
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

    fn crown(store: &Store, node: usize) -> Vec<usize> {
        store.crown(NodeIndex(node)).map(|hedge| hedge.0).collect()
    }

    pub(super) fn crowns(store: &Store) -> Vec<Vec<usize>> {
        (0..store.node_bound())
            .map(|node| crown(store, node))
            .collect()
    }

    fn payloads(store: &Store) -> Vec<usize> {
        store.iter().map(|(_, data)| *data).collect()
    }

    fn images<I: From<usize> + Into<usize>>(
        remap: &crate::half_edge::remap::Remap<I>,
    ) -> Vec<Option<usize>> {
        remap.iter().map(|(_, new)| new.map(Into::into)).collect()
    }

    fn subgraph(n_hedges: usize, hedges: &[usize]) -> SuBitGraph {
        let mut subgraph = SuBitGraph::empty(n_hedges);
        for &hedge in hedges {
            subgraph.add(Hedge(hedge));
        }
        subgraph
    }

    /// The renumbering that exchanges `a` and `b` among `n` indices.
    fn transposition<I: From<usize> + Into<usize>>(
        n: usize,
        a: usize,
        b: usize,
    ) -> crate::half_edge::remap::Remap<I> {
        let mut map: Vec<usize> = (0..n).collect();
        map.swap(a, b);
        crate::half_edge::remap::Remap::from_permutation(&Permutation::from_map(map))
    }

    #[cfg(feature = "serde")]
    type Graph = crate::half_edge::HedgeGraph<usize, usize, usize, NodeForest>;

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

    #[test]
    fn build_orders_crowns_and_labels_half_edges() {
        let store = store(&[&[3, 0, 5], &[], &[4, 1, 2]]);
        store.check_nodes().unwrap();
        store.topology.check_rings().unwrap();
        assert_eq!(crowns(&store), [vec![0, 3, 5], vec![], vec![1, 2, 4]]);
        assert_eq!(store.crown(NodeIndex(2)).len(), 3);
        let nodes: Vec<_> = (0..6).map(|h| store.node_id(Hedge(h)).0).collect();
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
        assert_eq!(store.identify(&[n(0), n(2)], 10), n(2));
        assert_eq!(crown(&store, 2), [0, 4, 5, 7, 8]);
        assert!(crown(&store, 0).is_empty() && !store.is_live(n(0)));
        assert_eq!(store.node_id(Hedge(5)), n(2));
        assert_eq!((store.n_live(), store.node_bound()), (3, 4));
        // Absorbed nodes keep their payload; `0` resolves to the five-hedge node 2,
        // which ties with node 1, so the first requested node survives.
        assert_eq!(store.get(n(0)), &0);
        assert_eq!(store.identify(&[n(1), n(0), n(1)], 11), n(1));
        assert_eq!(crown(&store, 1), (0..10).collect::<Vec<_>>());
        // A repeated isolated node only takes the new payload.
        assert_eq!(store.identify(&[n(3), n(3)], 12), n(3));
        assert_eq!(store.get(n(3)), &12);
        // An absorbed node resolves to the node that absorbed it last.
        assert_eq!(store.topology().representative(n(2)), n(1));
        assert_eq!(store.identify(&[n(3), n(2)], 13), n(1));
        assert_eq!((store.n_live(), store.node_bound()), (1, 4));
        store.check_nodes().unwrap();

        let remap = store.compact();
        // Every node maps to the survivor, which is the only origin.
        assert_eq!(images(&remap), [Some(0); 4]);
        assert_eq!(
            (0..4)
                .filter(|&old| remap.is_merged(n(old)))
                .collect::<Vec<_>>(),
            [0, 2, 3]
        );
        assert_eq!(crowns(&store), [(0..10).collect::<Vec<_>>()]);
        assert_eq!(store.get(n(0)), &13);
        store.check_nodes().unwrap();
        store.topology.check_rings().unwrap();
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
            assert_eq!(store.identify(&[NodeIndex(node), hub], 0), hub);
            store.check_nodes().unwrap();
        }
        assert_eq!(crown(&store, 5), (0..60).collect::<Vec<_>>());
    }

    #[test]
    fn compaction_fills_absorbed_slots_from_the_tail() {
        let mut store = store(&[&[0], &[1, 2], &[], &[3], &[]]);
        store.identify(&[NodeIndex(1), NodeIndex(0)], 10);
        store.identify(&[NodeIndex(3), NodeIndex(2)], 11);
        let remap = store.compact();
        // Slot 0 takes the isolated node 4, slot 2 the survivor 3; absorbed
        // nodes follow their survivor.
        assert_eq!(
            images(&remap),
            [Some(1), Some(1), Some(2), Some(2), Some(0)]
        );
        assert_eq!(remap.new_bound(), 3);
        assert_eq!(payloads(&store), [4, 10, 11]);
        assert_eq!(crowns(&store), [vec![], vec![0, 1, 2], vec![3]]);
        assert_eq!(store.node_id(Hedge(3)), NodeIndex(2));
        store.check_nodes().unwrap();
        assert!(store.compact().is_identity());
    }

    /// A hub at the highest index moves into the first hole with its class, so
    /// none of its half-edges is relabelled; the image matches a fresh build.
    #[test]
    fn compaction_moves_a_hub_without_relabelling_its_crown() {
        let mut moved = store(&[&[0], &[1], &[2], &[3, 4, 5, 6]]);
        let hub = NodeIndex(3);
        moved.set_anchor(hub, Some(Hedge(5)));
        moved.identify(&[NodeIndex(0), hub], 10);
        moved.identify(&[hub, NodeIndex(2)], 11);
        let labels = moved.topology.label.clone();
        let remap = moved.compact();
        assert_eq!(images(&remap), [Some(0), Some(1), Some(0), Some(0)]);
        assert_eq!(moved.topology.label, labels);
        // Classes 2 and 0 of the dropped nodes are free behind the node slots.
        assert_eq!(
            (moved.topology.classes.iter())
                .map(|entry| entry.class)
                .collect::<Vec<_>>(),
            [3, 1, 2, 0]
        );
        assert_eq!((moved.n_live(), moved.node_bound()), (2, 2));
        assert_eq!(crowns(&moved), [vec![0, 2, 3, 4, 5, 6], vec![1]]);
        assert_eq!(moved.node_id(Hedge(4)), NodeIndex(0));
        assert_eq!(moved.anchor(NodeIndex(0)), Some(Hedge(5)));
        moved.check_nodes().unwrap();
        moved.topology.check_rings().unwrap();
        let mut fresh = store(&[&[0, 2, 3, 4, 5, 6], &[1]]);
        *fresh.get_mut(NodeIndex(0)) = 11;
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

    /// Appended nodes reuse the classes that compaction freed, so repeated
    /// identify, compact and append cycles keep the class table bounded, also
    /// when the appended store has moved classes of its own.
    #[test]
    fn classes_are_reused_across_identify_compact_and_append() {
        let mut store = store(&[&[0, 1], &[2], &[3]]);
        for round in 0..20 {
            // Its hub sits last and moves into slot 0 when compacting.
            let mut other = super::tests::store(&[&[0], &[3], &[1, 2]]);
            other.identify(&[NodeIndex(0), NodeIndex(2)], 100);
            other.compact();
            assert_eq!(other.topology.classes[0], ClassEntry { class: 2, slot: 2 });
            let shift = store.node_bound();
            let n_hedges = store.n_hedges();
            store.extend(other);
            store.check_nodes().unwrap();
            store.topology.check_rings().unwrap();
            let classes = &store.topology.classes;
            assert!(classes.len() <= 5, "round {round}: {classes:?}");
            assert_eq!(crown(&store, shift), [n_hedges, n_hedges + 1, n_hedges + 2]);
            assert_eq!(crown(&store, shift + 1), [n_hedges + 3]);
            // Merge both appended nodes into the node of half-edge 0 and compact.
            let first = store.node_id(Hedge(0));
            store.identify(&[NodeIndex(shift), first], round);
            store.identify(&[first, NodeIndex(shift + 1)], round);
            store.compact();
            store.check_nodes().unwrap();
            assert_eq!(store.node_bound(), 3);
            let merged = store.node_id(Hedge(0)).0;
            let mut expected = vec![0, 1];
            expected.extend(4..store.n_hedges());
            assert_eq!(crown(&store, merged), expected);
            assert_eq!(crown(&store, store.node_id(Hedge(2)).0), [2]);
            assert_eq!(crown(&store, store.node_id(Hedge(3)).0), [3]);
        }
    }

    #[test]
    fn anchors_follow_their_half_edge_and_fall_back_to_the_smallest() {
        let mut store = store(&[&[0, 3], &[1, 4], &[2, 5], &[6]]);
        store.set_anchor(NodeIndex(1), Some(Hedge(4)));
        store.set_anchor(NodeIndex(2), Some(Hedge(5)));
        // The survivor without an anchor adopts the first absorbed one.
        store.identify(&[NodeIndex(0), NodeIndex(2), NodeIndex(1)], 10);
        assert_eq!(store.anchor(NodeIndex(0)), Some(Hedge(5)));
        store.compact();
        store.set_anchor(NodeIndex(1), Some(Hedge(6)));
        assert_eq!(crowns(&store), [vec![0, 1, 2, 3, 4, 5], vec![6]]);

        store
            .permute(&transposition(7, 5, 6), &NodeRemap::identity(2))
            .unwrap();
        assert_eq!(store.anchor(NodeIndex(0)), Some(Hedge(6)));
        assert_eq!(store.anchor(NodeIndex(1)), Some(Hedge(5)));
        assert_eq!(crowns(&store), [vec![0, 1, 2, 3, 4, 6], vec![5]]);

        // Extracting an anchor leaves the smallest remaining half-edge behind.
        let (extracted, _) = store.extract(&subgraph(7, &[2, 6]), |d| *d, |d| d).unwrap();
        let node =
            |store: &Store, payload| store.iter().find(|(_, data)| **data == payload).unwrap().0;
        assert_eq!(crowns(&store), [vec![2], vec![0, 1, 3, 4]]);
        assert_eq!(store.anchor(node(&store, 10)), Some(Hedge(0)));
        assert_eq!(store.anchor(node(&store, 3)), Some(Hedge(2)));
        assert_eq!(crown(&extracted, 0), [0, 1]);
        assert_eq!(extracted.anchor(NodeIndex(0)), Some(Hedge(1)));
        store.delete(&subgraph(5, &[0])).unwrap();
        assert_eq!(crowns(&store), [vec![2], vec![0, 1, 3]]);
        assert_eq!(store.anchor(node(&store, 10)), Some(Hedge(0)));
        store.delete(&subgraph(4, &[0, 1, 3])).unwrap();
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
    fn renumbering_requires_compacted_identifications() {
        let mut pending = store(&[&[0, 1], &[2], &[3]]);
        pending.identify(&[NodeIndex(1), NodeIndex(2)], 10);
        let sub = subgraph(4, &[0]);
        for operation in 0..5 {
            let run = |store: &mut Store| match operation {
                0 => store.delete(&sub).map(drop),
                1 => store.extract(&sub, |d| *d, |d| d).map(drop),
                2 => store.extract_nodes([NodeIndex(0)]).map(drop),
                // Compaction leaves two nodes.
                3 => store.permute(&transposition(4, 0, 3), &NodeRemap::identity(2)),
                _ => store.permute(&HedgeRemap::identity(4), &transposition(2, 0, 1)),
            };
            let mut store = pending.clone();
            assert!(matches!(
                run(&mut store),
                Err(HedgeGraphError::PendingIdentifications)
            ));
            assert_eq!(store, pending);
            store.compact();
            run(&mut store).unwrap();
            store.check_nodes().unwrap();
        }
    }

    #[test]
    fn dangling_edges_need_a_live_node() {
        let mut store = store(&[&[0], &[1], &[]]);
        store.identify(&[NodeIndex(0), NodeIndex(1)], 10);
        let before = store.clone();
        assert!(matches!(
            store.attach(NodeIndex(3)),
            Err(HedgeGraphError::NoNode)
        ));
        assert!(matches!(
            store.attach(NodeIndex(1)),
            Err(HedgeGraphError::InvalidNode(NodeIndex(1)))
        ));
        assert_eq!(
            (&store.topology.label, crowns(&store)),
            (&before.topology.label, crowns(&before))
        );
        assert_eq!(store.attach(NodeIndex(2)).unwrap(), Hedge(2));
        assert_eq!(store.attach(NodeIndex(0)).unwrap(), Hedge(3));
        assert_eq!(crowns(&store), [vec![0, 1, 3], vec![], vec![2]]);
        store.check_nodes().unwrap();
    }

    /// Appended nodes take the free classes left by compaction first, then mint
    /// new ones, and accept half-edges like any other node.
    #[test]
    fn added_nodes_take_free_classes_first() {
        let mut store = store(&[&[0], &[1], &[2]]);
        store.identify(&[NodeIndex(0), NodeIndex(1), NodeIndex(2)], 10);
        store.compact();
        assert_eq!(store.topology.classes.len(), 3);
        let added: Vec<_> = (20..24).map(|data| store.add_node(data)).collect();
        assert_eq!(added, (1..5).map(NodeIndex).collect::<Vec<_>>());
        assert_eq!(store.topology.classes.len(), 5);
        store.check_nodes().unwrap();
        store.topology.check_rings().unwrap();
        store.attach(added[1]).unwrap();
        store.attach(added[3]).unwrap();
        assert_eq!(
            crowns(&store),
            [vec![0, 1, 2], vec![], vec![3], vec![], vec![4]]
        );
        assert_eq!(payloads(&store), [10, 20, 21, 22, 23]);
        store.check_nodes().unwrap();
    }

    #[test]
    fn appending_shifts_pending_identifications() {
        let mut first = store(&[&[0], &[1]]);
        first.identify(&[NodeIndex(0), NodeIndex(1)], 10);
        let mut second = store(&[&[0, 2], &[1], &[]]);
        second.set_anchor(NodeIndex(1), Some(Hedge(1)));
        second.identify(&[NodeIndex(1), NodeIndex(0)], 11);
        first.extend(second);
        first.check_nodes().unwrap();
        first.topology.check_rings().unwrap();
        assert_eq!((first.n_live(), first.node_bound()), (3, 5));
        assert_eq!(
            crowns(&first),
            [vec![0, 1], vec![], vec![2, 3, 4], vec![], vec![]]
        );
        assert_eq!(first.anchor(NodeIndex(2)), Some(Hedge(3)));
        // The appended absorbed node still resolves to its survivor.
        assert_eq!(
            first.identify(&[NodeIndex(4), NodeIndex(3)], 12),
            NodeIndex(2)
        );
        let remap = first.compact();
        assert_eq!(remap.new_bound(), 2);
        assert_eq!(
            (0..5)
                .filter(|&old| remap.is_merged(NodeIndex(old)))
                .count(),
            3
        );
        first.check_nodes().unwrap();
    }

    /// Removal orders nodes as `[untouched | split | inside]` with the shared
    /// tail-fill kernel, moves the selected half-edges behind the others with the
    /// edge store's swaps and keeps isolated nodes.
    #[test]
    fn removal_renumbers_with_the_shared_kernel() {
        let parts: &[&[usize]] = &[&[0, 7], &[1, 2], &[], &[3, 8], &[4, 5], &[6], &[9]];
        let selection = subgraph(10, &[1, 2, 4, 8, 9]);
        let mut deleted = store(parts);
        let renumbering = deleted.delete(&selection).unwrap();
        deleted.check_nodes().unwrap();
        assert_eq!(payloads(&deleted), [0, 5, 2, 3, 4]);
        let nodes: Vec<_> = (0..5)
            .map(|hedge| deleted.node_id(Hedge(hedge)).0)
            .collect();
        assert_eq!(nodes, [0, 0, 1, 3, 4]);
        // Each surviving half-edge and node reports where it went.
        assert_eq!(
            images(&renumbering.hedges),
            [
                Some(0),
                None,
                None,
                Some(3),
                None,
                Some(4),
                Some(2),
                Some(1),
                None,
                None
            ]
        );
        assert_eq!(
            images(&renumbering.nodes),
            [Some(0), None, Some(2), Some(3), Some(4), Some(1), None]
        );

        let mut remaining = store(parts);
        let (extracted, split) = remaining.extract(&selection, |d| d + 100, |d| d).unwrap();
        assert_eq!(payloads(&remaining), payloads(&deleted));
        assert_eq!(split.kept, renumbering);
        // Split copies come first, then the nodes wholly inside.
        assert_eq!(payloads(&extracted), [100 + 3, 104, 1, 6]);
        assert_eq!(crowns(&extracted), [vec![3], vec![0], vec![1, 2], vec![4]]);
        assert_eq!(
            images(&split.moved.nodes),
            [None, Some(2), None, Some(0), Some(1), None, Some(3)]
        );
        assert_eq!(split.moved.hedges.get(Hedge(8)), Some(Hedge(3)));
        extracted.check_nodes().unwrap();
        remaining.check_nodes().unwrap();
    }

    #[test]
    fn extract_nodes_moves_selected_nodes_with_isolated_ones() {
        let mut store = store(&[&[0, 3], &[], &[1, 2], &[4]]);
        let (selected, extracted, split) =
            store.extract_nodes([NodeIndex(1), NodeIndex(2)]).unwrap();
        assert_eq!(
            selected.included_iter().map(|h| h.0).collect::<Vec<_>>(),
            [1, 2]
        );
        assert_eq!(payloads(&store), [0, 3]);
        assert_eq!(crowns(&store), [vec![0, 2], vec![1]]);
        assert_eq!(payloads(&extracted), [2, 1]);
        assert_eq!(crowns(&extracted), [vec![0, 1], vec![]]);
        assert_eq!(images(&split.kept.nodes), [Some(0), None, None, Some(1)]);
        assert_eq!(images(&split.moved.nodes), [None, Some(1), Some(0), None]);
        store.check_nodes().unwrap();
        extracted.check_nodes().unwrap();
    }

    /// Permutations rebuild the crowns ascending and move payloads and anchors
    /// with their nodes and half-edges.
    #[test]
    fn permutations_move_half_edges_nodes_and_anchors() {
        let parts: &[&[usize]] = &[
            &[0, 5, 9, 13],
            &[1, 2, 14],
            &[3, 4, 6, 7, 8],
            &[10],
            &[11, 12],
        ];
        let mut permuted = store(parts);
        permuted.set_anchor(NodeIndex(2), Some(Hedge(6)));
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut random = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        // The node of each half-edge, and the half-edge and node of each original.
        let mut model: Vec<usize> = vec![0; 15];
        for (node, hedges) in parts.iter().enumerate() {
            for &hedge in *hedges {
                model[hedge] = node;
            }
        }
        for _ in 0..50 {
            let mut hedge_map: Vec<usize> = (0..15).collect();
            let mut node_map: Vec<usize> = (0..5).collect();
            for i in (1..15).rev() {
                hedge_map.swap(i, random(i + 1));
            }
            for i in (1..5).rev() {
                node_map.swap(i, random(i + 1));
            }
            let anchor = permuted.anchor(permuted.node_id(Hedge(0))).map(|h| h.0);
            let payload_of_zero = *permuted.get(permuted.node_id(Hedge(0)));
            let hedges = HedgeRemap::from_permutation(&Permutation::from_map(hedge_map.clone()));
            let nodes = NodeRemap::from_permutation(&Permutation::from_map(node_map.clone()));
            permuted.permute(&hedges, &nodes).unwrap();
            let mut next = vec![0; 15];
            for (old, &node) in model.iter().enumerate() {
                next[hedge_map[old]] = node_map[node];
            }
            model = next;
            permuted.check_nodes().unwrap();
            permuted.topology.check_rings().unwrap();
            for (hedge, &node) in model.iter().enumerate() {
                assert_eq!(permuted.node_id(Hedge(hedge)).0, node);
            }
            let zero = permuted.node_id(Hedge(hedge_map[0]));
            assert_eq!(*permuted.get(zero), payload_of_zero);
            assert_eq!(
                permuted.anchor(zero).map(|h| h.0),
                anchor.map(|a| hedge_map[a])
            );
        }
    }

    fn hash(store: &Store) -> u64 {
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
        a.identify(&[NodeIndex(0), NodeIndex(1)], 7);
        a.identify(&[NodeIndex(0), NodeIndex(2)], 9);
        b.identify(&[NodeIndex(0), NodeIndex(2)], 7);
        b.identify(&[NodeIndex(1), NodeIndex(0)], 9);
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
        merged.identify(&[NodeIndex(0), NodeIndex(1)], 0);
        assert!(merged.topology.nodes().eq(unmerged.topology.nodes()));
        assert_eq!(payloads(&merged), payloads(&unmerged));
        assert_ne!(merged, unmerged);
        assert_eq!(merged.cmp(&unmerged), std::cmp::Ordering::Less);
        assert!(merged < unmerged);
        assert_ne!(hash(&merged), hash(&unmerged));
    }

    #[test]
    fn to_forest_roots_each_crown_at_its_smallest_half_edge() {
        let mut store = store(&[&[1, 3], &[0], &[2], &[]]);
        store.identify(&[NodeIndex(0), NodeIndex(1)], 10);
        let forest = store.to_forest::<usize, ()>(|d| *d);
        forest.validate_structure().unwrap();
        assert_eq!(forest.roots[0].root_id, TreeNodeId(0));
        assert!(forest.roots[1].root_id.is_empty() && forest.roots[3].root_id.is_empty());
        for hedge in 0..4 {
            let node = forest.root(TreeNodeId(hedge));
            assert_eq!(NodeIndex(node.0), store.node_id(Hedge(hedge)));
        }
    }

    /// Rolling back undoes appended nodes and half-edges, anchors and
    /// identifications, innermost checkpoint first; committing keeps them.
    #[test]
    fn checkpoints_roll_back_and_commit_in_place_edits() {
        let mut topology = store(&[&[0, 4], &[1, 5, 6], &[2], &[3]]).topology;
        topology.set_anchor(NodeIndex(0), Some(Hedge(4)));
        topology.identify(&[NodeIndex(2), NodeIndex(3)]);
        let original = topology.clone();
        let snapshot = |t: &NodeForest| {
            let crowns: Vec<Vec<Hedge>> = (0..t.node_bound())
                .map(|node| t.crown(NodeIndex(node)).collect())
                .collect();
            let anchors: Vec<_> = (0..t.node_bound())
                .map(|node| t.anchor(NodeIndex(node)))
                .collect();
            let representatives: Vec<_> = (0..t.node_bound())
                .map(|node| t.representative(NodeIndex(node)))
                .collect();
            (
                t.clone(),
                crowns,
                anchors,
                representatives,
                t.classes.clone(),
            )
        };
        let before = snapshot(&topology);

        topology.checkpoint();
        assert!(topology.has_checkpoint());
        let added = topology.push_node();
        topology.attach(added).unwrap();
        topology.attach(NodeIndex(1)).unwrap();
        // The absorbed node 3 stands for its survivor 2.
        assert_eq!(
            topology.identify(&[NodeIndex(0), NodeIndex(3)]),
            NodeIndex(0)
        );
        topology.checkpoint();
        assert_eq!(topology.identify(&[added, NodeIndex(1)]), NodeIndex(1));
        topology.set_anchor(NodeIndex(1), Some(Hedge(7)));
        let inner = snapshot(&topology);
        topology.checkpoint();
        topology.identify(&[NodeIndex(1), NodeIndex(0)]);
        topology.rollback();
        assert_eq!(snapshot(&topology), inner);
        topology.check().unwrap();
        topology.check_rings().unwrap();
        // While checkpoints are open, renumbering operations refuse.
        for result in [
            topology.compact().map(drop),
            topology
                .partition_hedges(|_| false, NodeFate::of_removal, false)
                .map(drop),
            topology.extend(original.clone()),
            topology.permute(&HedgeRemap::identity(9), &NodeRemap::identity(5)),
        ] {
            assert!(matches!(result, Err(HedgeGraphError::CheckpointOpen)));
        }
        topology.commit();
        topology.rollback();
        assert!(!topology.has_checkpoint());
        assert_eq!(snapshot(&topology), before);
        topology.check().unwrap();
        topology.check_rings().unwrap();

        topology.checkpoint();
        let added = topology.push_node();
        topology.attach(added).unwrap();
        topology.identify(&[added, NodeIndex(2)]);
        topology.commit();
        assert!(!topology.has_checkpoint());
        assert_eq!(
            topology
                .crown(NodeIndex(2))
                .map(|h| h.0)
                .collect::<Vec<_>>(),
            [2, 3, 7]
        );
        assert_eq!(topology.compact().unwrap().new_bound(), 3);
        topology.check().unwrap();
    }

    /// Identifications journaled under a checkpoint leave the chains of nodes
    /// absorbed before it uncompressed, so rolling back restores where they
    /// lead.
    #[test]
    fn rollback_restores_chains_absorbed_before_the_checkpoint() {
        let labels = [0, 1, 1, 2, 2, 2].map(NodeIndex);
        let mut topology = NodeForest::from_labels(labels.into_iter(), 3);
        assert_eq!(
            topology.identify(&[NodeIndex(1), NodeIndex(0)]),
            NodeIndex(1)
        );
        topology.checkpoint();
        assert_eq!(
            topology.identify(&[NodeIndex(2), NodeIndex(1)]),
            NodeIndex(2)
        );
        // Resolving 0 walks 0 -> 1 -> 2 and must not point 0 at 2.
        assert_eq!(
            topology.identify(&[NodeIndex(0), NodeIndex(2)]),
            NodeIndex(2)
        );
        topology.rollback();
        assert_eq!(topology.representative(NodeIndex(0)), NodeIndex(1));
        topology.check().unwrap();
        let remap = topology.compact().unwrap();
        assert_eq!(remap.get(NodeIndex(0)), remap.get(NodeIndex(1)));
        assert_ne!(remap.get(NodeIndex(0)), remap.get(NodeIndex(2)));
    }

    /// Random edits rolled back restore the exact state, classes included.
    #[test]
    fn random_checkpointed_edits_roll_back_exactly() {
        let mut state = 0x853c_49e6_748f_ea9b_u64;
        let mut random = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound as u64) as usize
        };
        for _ in 0..200 {
            let n_nodes = 1 + random(6);
            let labels: Vec<_> = (0..random(12))
                .map(|_| NodeIndex(random(n_nodes)))
                .collect();
            let mut topology = NodeForest::from_labels(labels.into_iter(), n_nodes);
            // Identifications left pending build chains of absorbed nodes that
            // the checkpoint must restore.
            for _ in 0..random(4) {
                topology.identify(&[NodeIndex(random(n_nodes)), NodeIndex(random(n_nodes))]);
            }
            if random(2) == 0 {
                topology.compact().unwrap();
            }
            let representatives = |t: &NodeForest| -> Vec<_> {
                (0..t.node_bound())
                    .map(|node| t.representative(NodeIndex(node)))
                    .collect()
            };
            let before_representatives = representatives(&topology);
            let before = (
                topology.clone(),
                topology.label.clone(),
                topology.classes.clone(),
            );
            topology.checkpoint();
            for _ in 0..random(10) {
                let bound = topology.node_bound();
                match random(4) {
                    0 => {
                        topology.push_node();
                    }
                    1 => {
                        let _ = topology.attach(NodeIndex(random(bound)));
                    }
                    2 => {
                        let node = topology.representative(NodeIndex(random(bound)));
                        let crown: Vec<_> = topology.crown(node).collect();
                        let anchor = (!crown.is_empty()).then(|| crown[random(crown.len())]);
                        topology.set_anchor(node, anchor);
                    }
                    _ => {
                        let nodes: Vec<_> = (0..1 + random(3))
                            .map(|_| NodeIndex(random(bound)))
                            .collect();
                        topology.identify(&nodes);
                    }
                }
                topology.check().unwrap();
                topology.check_rings().unwrap();
            }
            topology.rollback();
            topology.check_rings().unwrap();
            assert_eq!(
                (
                    topology.clone(),
                    topology.label.clone(),
                    topology.classes.clone()
                ),
                before
            );
            for node in (0..topology.node_bound()).map(NodeIndex) {
                assert_eq!(topology.anchor(node), before.0.anchor(node));
                assert!(topology.crown(node).eq(before.0.crown(node)));
            }
            assert_eq!(representatives(&topology), before_representatives);
        }
    }

    #[cfg(feature = "serde")]
    #[test]
    fn serde_round_trips_pending_identifications_and_anchors() {
        let mut store = store(&[&[1, 3], &[0, 4], &[2], &[]]);
        store.set_anchor(NodeIndex(1), Some(Hedge(4)));
        store.identify(&[NodeIndex(0), NodeIndex(1)], 10);
        let json = serde_json::to_string(&store).unwrap();
        let back: Store = serde_json::from_str(&json).unwrap();
        back.check_nodes().unwrap();
        assert_eq!(back, store);
        assert_eq!(crowns(&back), crowns(&store));
        assert_eq!(back.anchor(NodeIndex(0)), Some(Hedge(4)));
        assert!(!back.is_live(NodeIndex(1)));
        assert_eq!(back.clone().identify(&[NodeIndex(1)], 0), NodeIndex(0));
        let graph = graph();
        let json = serde_json::to_string(&graph).unwrap();
        assert_eq!(serde_json::from_str::<Graph>(&json).unwrap(), graph);

        let store_json = |label: &str, anchor: &str, absorbed_into: &str, data: &str| {
            format!(
                r#"{{"topology":{{"label":{label},"anchor":{anchor},"absorbed_into":{absorbed_into}}},"data":{data}}}"#
            )
        };
        for invalid in [
            // A half-edge of a missing node.
            store_json("[0,5]", "[null]", "[null]", "[0]"),
            // An absorbed node with a crown.
            store_json("[0,1]", "[null,null]", "[null,0]", "[0,1]"),
            // An anchor outside the crown.
            store_json("[0,1]", "[1,null]", "[null,null]", "[0,1]"),
            // Absorbed nodes merged into each other.
            store_json("[]", "[null,null]", "[1,0]", "[0,1]"),
            // The reserved index u32::MAX as an anchor or a merge target.
            store_json("[0]", "[4294967295]", "[null]", "[0]"),
            store_json("[0]", "[null,null]", "[null,4294967295]", "[0,1]"),
            // One payload too few.
            store_json("[0]", "[null,null]", "[null,null]", "[0]"),
        ] {
            assert!(
                serde_json::from_str::<Store>(&invalid).is_err(),
                "{invalid}"
            );
        }
    }

    #[cfg(feature = "bincode")]
    #[test]
    fn bincode_round_trips() {
        let mut store = store(&[&[1, 3], &[0, 4], &[2], &[]]);
        store.set_anchor(NodeIndex(1), Some(Hedge(4)));
        store.identify(&[NodeIndex(0), NodeIndex(1)], 10);
        let config = bincode::config::standard();
        let bytes = bincode::encode_to_vec(&store, config).unwrap();
        let (back, _): (Store, _) = bincode::decode_from_slice(&bytes, config).unwrap();
        assert_eq!(back, store);
        assert_eq!(crowns(&back), crowns(&store));
        assert_eq!(back.anchor(NodeIndex(0)), Some(Hedge(4)));
        #[cfg(feature = "serde")]
        {
            let graph = graph();
            let bytes = bincode::encode_to_vec(&graph, config).unwrap();
            let (back, _): (Graph, _) = bincode::decode_from_slice(&bytes, config).unwrap();
            assert_eq!(back, graph);
        }
    }

    /// Random identifications of possibly absorbed or repeated nodes, dangling
    /// edges and compaction, against a union-find model with crown sizes.
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
            let mut store: Store = NodeStore::build(
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
                        assert_eq!(store.identify(&requested, 100 + step), NodeIndex(survivor));
                    }
                    3 => {
                        let node = random(parent.len());
                        let before = store.clone();
                        match store.attach(NodeIndex(node)) {
                            Ok(hedge) => {
                                assert_eq!(find(&parent, node), node);
                                assert_eq!(hedge, Hedge(owner.len()));
                                owner.push(node);
                                size[node] += 1;
                            }
                            Err(HedgeGraphError::InvalidNode(_)) => {
                                assert_ne!(find(&parent, node), node);
                                assert_eq!(store, before);
                            }
                            Err(error) => panic!("{error}"),
                        }
                    }
                    _ => {
                        let survivors: BTreeSet<_> = (0..parent.len())
                            .filter(|&node| find(&parent, node) == node)
                            .map(|node| *store.get(NodeIndex(node)))
                            .collect();
                        let expected: Vec<_> = owner
                            .iter()
                            .map(|&node| *store.get(NodeIndex(find(&parent, node))))
                            .collect();
                        let remap = store.compact();
                        for (old, new) in remap.iter() {
                            let survivor = find(&parent, old.0);
                            assert_eq!(new, remap.get(NodeIndex(survivor)));
                            assert_eq!(remap.is_merged(old), survivor != old.0);
                        }
                        let nodes = store.iter().map(|(_, data)| *data);
                        assert_eq!(nodes.collect::<BTreeSet<_>>(), survivors);
                        owner = (0..owner.len())
                            .map(|h| store.node_id(Hedge(h)).0)
                            .collect();
                        let actual: Vec<_> = owner
                            .iter()
                            .map(|&node| *store.get(NodeIndex(node)))
                            .collect();
                        assert_eq!(actual, expected);
                        parent = (0..store.node_bound()).collect();
                        size = (0..store.node_bound())
                            .map(|node| crown(&store, node).len())
                            .collect();
                    }
                }
                store.check_nodes().unwrap();
                store.topology.check_rings().unwrap();
                for (hedge, &node) in owner.iter().enumerate() {
                    assert_eq!(store.node_id(Hedge(hedge)), NodeIndex(find(&parent, node)));
                }
                let live = (0..parent.len()).filter(|&n| find(&parent, n) == n).count();
                assert_eq!(store.n_live(), live);
            }
        }
    }
}
