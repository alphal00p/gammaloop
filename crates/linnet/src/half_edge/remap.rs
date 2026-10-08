//! Typed renumberings of half-edges and nodes.
//!
//! Every operation that moves half-edges or renumbers nodes reports where each
//! old index went as a [`Remap`]: [`super::HedgeGraph::compact`],
//! [`super::HedgeGraph::delete_hedges`], [`super::HedgeGraph::extract`],
//! [`super::HedgeGraph::extract_nodes`] and [`super::HedgeGraph::contract_subgraph`]
//! return them, and [`super::HedgeGraph::permute_hedges`] and
//! [`super::HedgeGraph::permute_nodes`] accept them. All of them derive from the
//! shared tail-fill kernel [`partition_with`].

use std::{marker::PhantomData, ops::Range};

use super::{involution::Hedge, swap::partition_with, NodeIndex};
use crate::permutation::Permutation;

/// Flags an old index that maps to a new index it is not the origin of.
const MERGED: u32 = 1 << 31;
/// Marks an old index without a new one.
const NONE: u32 = u32::MAX;

/// Where each old index went: a dense table from old to new indices.
///
/// An old index may have no new one (it was removed or moved elsewhere) and,
/// for nodes, several old indices may share one new index: the nodes merged
/// into a survivor by an identification map to the survivor's new index. Each
/// new index has exactly one *origin*, the old index whose record, and payload,
/// moved there; merged indices are not origins. Tables hold fewer than `2^31`
/// indices.
#[derive(Clone)]
pub struct Remap<I> {
    /// Each old index's position, flagged with [`MERGED`], or [`NONE`]; the
    /// new indices are the positions in `start..start + bound`, shifted down.
    /// Splitting a partition into its sides only moves this window.
    new: Vec<u32>,
    start: u32,
    bound: usize,
    id: PhantomData<fn() -> I>,
}

/// Where each old node went.
pub type NodeRemap = Remap<NodeIndex>;
/// Where each old half-edge went.
pub type HedgeRemap = Remap<Hedge>;

impl<I> Remap<I> {
    /// The new index of `entry`, and whether its old index is the origin.
    #[inline]
    fn resolve(&self, entry: u32) -> Option<(u32, bool)> {
        let new = (entry & !MERGED).wrapping_sub(self.start);
        (entry != NONE && (new as usize) < self.bound).then_some((new, entry & MERGED == 0))
    }

    /// The new index of each old index, and whether it is the origin.
    fn resolved(&self) -> impl ExactSizeIterator<Item = Option<(u32, bool)>> + '_ {
        self.new.iter().map(|&entry| self.resolve(entry))
    }
}

impl<I> PartialEq for Remap<I> {
    fn eq(&self, other: &Self) -> bool {
        self.bound == other.bound && self.resolved().eq(other.resolved())
    }
}

impl<I> Eq for Remap<I> {}

impl<I> std::fmt::Debug for Remap<I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let new: Vec<_> = self.resolved().map(|new| new.map(|(new, _)| new)).collect();
        let merged = (self.resolved().enumerate())
            .filter(|(_, new)| new.is_some_and(|(_, origin)| !origin))
            .map(|(old, _)| old);
        f.debug_struct("Remap")
            .field("new", &new)
            .field("merged", &merged.collect::<Vec<_>>())
            .field("bound", &self.bound)
            .finish()
    }
}

impl<I: From<usize> + Into<usize>> Remap<I> {
    fn checked(new: Vec<u32>, bound: usize) -> Self {
        assert!(
            new.len() < MERGED as usize,
            "renumberings hold fewer than 2^31 indices"
        );
        Remap {
            new,
            start: 0,
            bound,
            id: PhantomData,
        }
    }

    /// The table of a bijection from the new index of every old one.
    pub(crate) fn bijection(new: Vec<u32>) -> Self {
        let bound = new.len();
        Self::checked(new, bound)
    }

    /// The table that keeps every one of `n` indices in place.
    pub fn identity(n: usize) -> Self {
        Self::bijection((0..index_u32(n)).collect())
    }

    /// The renumbering that [`crate::half_edge::swap::Swap::permute`] applies
    /// with `permutation`: old index `i` moves to `permutation.map()[i]`.
    pub fn from_permutation(permutation: &Permutation) -> Self {
        Self::bijection(
            permutation
                .map()
                .iter()
                .map(|&new| index_u32(new))
                .collect(),
        )
    }

    /// Compacts `n` slots of `ctx` with the shared kernel: live slots fill the
    /// holes of absorbed ones from the tail, calling `swap(ctx, hole, live)`
    /// for each, and absorbed slots merge into the new index of the live slot
    /// their chain of `parent`s, the slots they were merged into, ends at.
    ///
    /// Beyond the kernel's scan, this fills the table and makes three passes
    /// over the absorbed slots.
    pub(crate) fn compaction<C: ?Sized>(
        ctx: &mut C,
        n: usize,
        is_live: impl Fn(&C, usize) -> bool,
        parent: impl Fn(&C, usize) -> usize,
        mut swap: impl FnMut(&mut C, usize, usize),
    ) -> Self {
        // Live slots hold their own index until the kernel moves them. Absorbed
        // slots hold their parent, flagged, read while their record is still
        // in place: holes just before the kernel fills them, and the absorbed
        // slots behind the live ones after it, where only these stay put.
        let mut remap = Self::checked((0..index_u32(n)).collect(), 0);
        let live = partition_with(
            &mut (&mut *ctx, &mut remap.new),
            0,
            n,
            |(ctx, _), slot| is_live(ctx, slot),
            |(ctx, new), hole, moved| {
                new[hole] = MERGED | parent(ctx, hole) as u32;
                swap(ctx, hole, moved);
                new[moved] = hole as u32;
            },
        );
        remap.bound = live;
        // Each slot behind the live ones is absorbed, or moved and names the
        // absorbed hole it filled.
        let absorbed = |new: &[u32], slot: usize| match new[slot] & MERGED {
            0 => new[slot] as usize,
            _ => slot,
        };
        for slot in live..n {
            if remap.new[slot] as usize == slot {
                remap.new[slot] = MERGED | parent(ctx, slot) as u32;
            }
        }
        // Point every absorbed slot at the old index of its live root, then at
        // the root's new index; flagged entries hold old indices until then.
        for slot in live..n {
            let slot = absorbed(&remap.new, slot);
            let mut root = (remap.new[slot] & !MERGED) as usize;
            while remap.new[root] & MERGED != 0 {
                root = (remap.new[root] & !MERGED) as usize;
            }
            remap.new[slot] = MERGED | root as u32;
        }
        for slot in live..n {
            let slot = absorbed(&remap.new, slot);
            let root = (remap.new[slot] & !MERGED) as usize;
            remap.new[slot] = MERGED | remap.new[root];
        }
        remap
    }

    /// The table sending each old index to its new index, if any, and whether
    /// it is the origin there.
    #[cfg(test)]
    pub(crate) fn from_entries(
        entries: impl IntoIterator<Item = Option<(usize, bool)>>,
        bound: usize,
    ) -> Self {
        let new = (entries.into_iter())
            .map(|entry| match entry {
                Some((new, origin)) => index_u32(new) | if origin { 0 } else { MERGED },
                None => NONE,
            })
            .collect();
        Self::checked(new, bound)
    }

    /// The new index of `old`, if it has one.
    #[inline]
    pub fn get(&self, old: I) -> Option<I> {
        let (new, _) = self.resolve(self.new[old.into()])?;
        Some(I::from(new as usize))
    }

    /// The number of old indices.
    pub fn old_bound(&self) -> usize {
        self.new.len()
    }

    /// The number of new indices.
    pub fn new_bound(&self) -> usize {
        self.bound
    }

    /// Whether `old` maps to a new index whose origin is another old index.
    pub fn is_merged(&self, old: I) -> bool {
        (self.resolve(self.new[old.into()])).is_some_and(|(_, origin)| !origin)
    }

    /// Whether every index keeps its place.
    pub fn is_identity(&self) -> bool {
        self.bound == self.new.len()
            && (self.resolved().enumerate()).all(|(old, new)| new == Some((old as u32, true)))
    }

    /// Each old index with its new one.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (I, Option<I>)> + '_ {
        (self.resolved().enumerate())
            .map(|(old, new)| (I::from(old), new.map(|(new, _)| I::from(new as usize))))
    }

    /// The renumbering of `self` followed by `next`, whose old indices are the
    /// new indices of `self`; an index is merged if either step merged it.
    ///
    /// # Panics
    ///
    /// If `next` does not start from the new indices of `self`.
    pub fn then(&self, next: &Self) -> Self {
        assert_eq!(
            self.bound,
            next.new.len(),
            "composed renumberings must chain"
        );
        let new = (self.resolved())
            .map(|first| {
                let (middle, origin) = first?;
                let (new, next_origin) = next.resolve(next.new[middle as usize])?;
                Some(new | if origin && next_origin { 0 } else { MERGED })
            })
            .map(|entry| entry.unwrap_or(NONE))
            .collect();
        Self::checked(new, next.bound)
    }

    /// The permutation this bijection applies, for stores that permute by swaps.
    ///
    /// # Panics
    ///
    /// If some old index has no new index of its own.
    pub(crate) fn to_permutation(&self) -> Permutation {
        let map = (self.resolved())
            .map(|new| match new {
                Some((new, true)) => new as usize,
                _ => panic!("only bijections are permutations"),
            })
            .collect();
        Permutation::from_map(map)
    }

    /// The new indices in `range`, renumbered from its start; the others are
    /// dropped. Splits a bijection over `[kept | moved]` into its sides in
    /// constant time.
    pub(crate) fn window(mut self, range: Range<usize>) -> Self {
        assert!(
            range.end <= self.bound,
            "windows lie inside the new indices"
        );
        self.start += range.start as u32;
        self.bound = range.len();
        self
    }

    /// Moves the value at each origin to its new index and drops the others.
    ///
    /// # Panics
    ///
    /// If `values` does not hold one value per old index.
    pub(crate) fn apply<T>(&self, values: &mut Vec<T>) {
        assert_eq!(values.len(), self.new.len(), "one value per old index");
        // Complete the table to a permutation that sends the dropped values
        // behind the new indices, then follow its cycles.
        let mut dropped = self.bound as u32;
        let mut target: Vec<u32> = (self.resolved())
            .map(|new| match new {
                Some((new, true)) => new,
                _ => {
                    dropped += 1;
                    dropped - 1
                }
            })
            .collect();
        for start in 0..target.len() {
            while target[start] as usize != start {
                let to = target[start] as usize;
                values.swap(start, to);
                target.swap(start, to);
            }
        }
        values.truncate(self.bound);
    }

    /// Like [`Self::apply`] for a bijection without window, which it borrows to
    /// mark the visited cycles instead of allocating.
    pub(crate) fn permute_in_place<T>(&mut self, values: &mut [T]) {
        assert!(
            self.start == 0 && self.bound == self.new.len() && values.len() == self.bound,
            "one value per index of a bijection"
        );
        for start in 0..self.new.len() {
            if self.new[start] & MERGED != 0 {
                continue;
            }
            // Rotate the cycle through `start`: each swap settles one value.
            let mut next = self.new[start];
            self.new[start] |= MERGED;
            while next as usize != start {
                values.swap(start, next as usize);
                let after = self.new[next as usize];
                self.new[next as usize] |= MERGED;
                next = after;
            }
        }
        for entry in &mut self.new {
            *entry &= !MERGED;
        }
    }

    /// Like [`Self::apply`] for a [`Self::compaction`], whose origins below the
    /// new bound stay: each origin behind it fills the hole of a dropped value,
    /// in time linear in the values behind the bound.
    pub(crate) fn fill_holes<T>(&self, values: &mut Vec<T>) {
        assert_eq!(values.len(), self.new.len(), "one value per old index");
        debug_assert!(
            (self.resolved().enumerate().take(self.bound))
                .all(|(old, new)| new.is_none_or(|(new, origin)| !origin || new as usize == old)),
            "only compactions fill holes"
        );
        for (old, &entry) in self.new.iter().enumerate().skip(self.bound) {
            if let Some((new, true)) = self.resolve(entry) {
                values.swap(old, new as usize);
            }
        }
        values.truncate(self.bound);
    }
}

/// How an operation renumbered half-edges and nodes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Renumbering {
    pub hedges: HedgeRemap,
    pub nodes: NodeRemap,
}

/// The renumberings into the two graphs an extraction leaves: the half-edges
/// and nodes that stay, and those moved out. A node crossed by the selection
/// boundary appears on both sides.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SplitRenumbering {
    pub kept: Renumbering,
    pub moved: Renumbering,
}

/// Converts an index to `u32`. Indices stay below `2^31 - 1`, so remaps can
/// flag merges in the top bit and keep `u32::MAX` free as a sentinel.
pub(crate) fn index_u32(index: usize) -> u32 {
    u32::try_from(index)
        .ok()
        .filter(|&i| i < NONE & !MERGED)
        .expect("indices must fit below 2^31 - 1")
}

/// Partitions the positions `start..` of `order`, the old index at each
/// position, with the shared kernel so that those accepted by `keep` come
/// first; returns the split point.
pub(crate) fn partition_order(
    order: &mut [u32],
    start: usize,
    keep: impl Fn(u32) -> bool,
) -> usize {
    let end = order.len();
    partition_with(
        order,
        start,
        end,
        |order, position| keep(order[position]),
        |order, a, b| order.swap(a, b),
    )
}

/// The inverse of `order`: each old index's position.
pub(crate) fn positions(order: &[u32]) -> Vec<u32> {
    let mut position = vec![NONE; order.len()];
    for (at, &old) in order.iter().enumerate() {
        position[old as usize] = at as u32;
    }
    position
}

#[cfg(test)]
mod tests {
    use super::*;

    fn images<I: From<usize> + Into<usize>>(remap: &Remap<I>) -> Vec<Option<usize>> {
        remap.iter().map(|(_, new)| new.map(Into::into)).collect()
    }

    /// Slots 0 and 3 are absorbed into 1 and 4, which stay live with 2.
    fn compacted() -> (NodeRemap, Vec<usize>) {
        let mut slots: Vec<usize> = (0..5).collect();
        let parent = [1, 1, 2, 4, 4];
        let remap = NodeRemap::compaction(
            &mut slots,
            5,
            |_, slot| parent[slot] == slot,
            |_, slot| parent[slot],
            |slots, hole, live| slots.swap(hole, live),
        );
        (remap, slots)
    }

    #[test]
    fn compaction_fills_holes_from_the_tail_and_merges_into_survivors() {
        let (remap, slots) = compacted();
        assert_eq!(&slots[..3], [4, 1, 2]);
        assert_eq!(
            images(&remap),
            [Some(1), Some(1), Some(2), Some(0), Some(0)]
        );
        let merged: Vec<_> = (0..5)
            .filter(|&old| remap.is_merged(NodeIndex(old)))
            .collect();
        assert_eq!((merged, remap.new_bound()), (vec![0, 3], 3));
        let mut values = vec!['a', 'b', 'c', 'd', 'e'];
        remap.fill_holes(&mut values);
        assert_eq!(values, ['e', 'b', 'c']);
        let mut general = vec!['a', 'b', 'c', 'd', 'e'];
        remap.apply(&mut general);
        assert_eq!(general, values);
    }

    /// Chains of absorbed slots resolve through filled holes, slots still in
    /// place and moved live slots alike.
    #[test]
    fn compaction_follows_chains_to_their_live_root() {
        let mut slots: Vec<usize> = (0..6).collect();
        let parent = [2, 1, 3, 4, 4, 5];
        let remap = NodeRemap::compaction(
            &mut slots,
            6,
            |_, slot| parent[slot] == slot,
            |_, slot| parent[slot],
            |slots, hole, live| slots.swap(hole, live),
        );
        assert_eq!(&slots[..3], [5, 1, 4]);
        assert_eq!(
            images(&remap),
            [Some(2), Some(1), Some(2), Some(2), Some(2), Some(0)]
        );
        let merged: Vec<_> = (0..6)
            .filter(|&old| remap.is_merged(NodeIndex(old)))
            .collect();
        assert_eq!(merged, [0, 2, 3]);
    }

    #[test]
    fn apply_moves_origins_and_drops_the_rest() {
        // A permutation of the first three values; the last two are dropped.
        let map = NodeRemap::bijection(vec![2, 0, 1, 3, 4]).window(0..3);
        let mut values = vec!['a', 'b', 'c', 'd', 'e'];
        map.apply(&mut values);
        assert_eq!(values, ['b', 'c', 'a']);
        assert_eq!(map.get(NodeIndex(3)), None);
    }

    #[test]
    fn composition_follows_both_steps_and_keeps_merges() {
        let (first, _) = compacted();
        // Drop the new node 1, keep 0 and 2 as 0 and 1.
        let second = NodeRemap::bijection(vec![0, 2, 1]).window(0..2);
        let both = first.then(&second);
        assert_eq!(images(&both), [None, None, Some(1), Some(0), Some(0)]);
        assert!(both.is_merged(NodeIndex(3)) && !both.is_merged(NodeIndex(4)));
        assert!(!both.is_merged(NodeIndex(0)));
        assert_eq!(both.new_bound(), 2);
        assert!(NodeRemap::identity(3).is_identity());
        assert!(!first.is_identity());
    }

    #[test]
    fn windows_split_a_bijection_into_sides() {
        let all = HedgeRemap::bijection(vec![3, 0, 2, 1]);
        let kept = all.clone().window(0..2);
        let moved = all.window(2..4);
        assert_eq!(kept.get(Hedge(1)), Some(Hedge(0)));
        assert_eq!(kept.get(Hedge(0)), None);
        assert_eq!(moved.get(Hedge(0)), Some(Hedge(1)));
        assert_eq!((kept.new_bound(), moved.new_bound()), (2, 2));
    }

    #[test]
    fn permutations_match_swap_permute() {
        use crate::half_edge::swap::Swap;
        let permutation = Permutation::from_map(vec![2, 0, 3, 1]);
        let mut swapped: crate::half_edge::NodeVec<usize> = vec![10, 11, 12, 13].into();
        swapped.permute(&permutation);
        let remap = NodeRemap::from_permutation(&permutation);
        let mut applied = vec![10, 11, 12, 13];
        remap.apply(&mut applied);
        assert_eq!(swapped.0, applied);
        assert_eq!(remap.to_permutation(), permutation);
        let mut permuted = vec![10, 11, 12, 13];
        let mut borrowed = remap.clone();
        borrowed.permute_in_place(&mut permuted);
        assert_eq!((permuted, borrowed), (applied, remap));
    }
}
