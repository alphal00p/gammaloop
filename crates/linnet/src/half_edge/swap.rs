use std::ops::{AddAssign, SubAssign};

use crate::permutation::Permutation;

use super::{
    hedgevec::SmartEdgeVec,
    involution::{EdgeIndex, Hedge},
    nodestore::NodeStorageOps,
    HedgeGraph, NodeIndex,
};

#[allow(clippy::len_without_is_empty)]
pub trait Swap<Id> {
    // type Item;

    fn swap(&mut self, i: Id, j: Id);

    fn len(&self) -> Id;

    fn is_zero_length(&self) -> bool;

    fn permute(&mut self, perm: &Permutation)
    where
        Id: From<usize>,
    {
        let trans = perm.transpositions();

        // println!("trans:{trans:?}");

        for (i, j) in trans.into_iter().rev() {
            self.swap(Id::from(i), Id::from(j));
        }
    }

    // fn filter(&self, id: &Id, filter: &impl Fn(&Id, &Self::Item) -> bool) -> bool;

    /// Moves every position accepted by `keep` in front of the rejected ones and
    /// returns the split point, swapping as [`partition_with`] does.
    fn partition(&mut self, mut keep: impl FnMut(&Id) -> bool) -> Id
    where
        Id: From<usize> + AddAssign + PartialOrd + SubAssign + Copy,
    {
        self.partition_from(Id::from(0), |_, id| keep(&id))
    }

    /// Partitions `start..len` like [`Swap::partition`], leaving `..start` in
    /// place; `keep` may inspect the current state.
    fn partition_from(&mut self, start: Id, keep: impl FnMut(&Self, Id) -> bool) -> Id
    where
        Id: From<usize> + AddAssign + PartialOrd + SubAssign + Copy,
    {
        let end = self.len();
        partition_with(self, start, end, keep, |s, i, j| s.swap(i, j))
    }
}

/// The two-pointer compaction behind every delete, extract and forget.
///
/// Partitions `start..end` of `ctx` so that the positions accepted by `keep` come
/// first and returns the split point. `left` scans up over kept positions; at a
/// rejected `left`, `right` steps down and the two are swapped only if `right` is
/// kept, so the holes are filled from the tail. Stores compacted with the same
/// predicate stay aligned because they all swap in this order.
///
/// `keep` sees a position only before any swap moves it, so a predicate on the
/// original layout and one on the current state agree.
// Without the hint rustc keeps the kernel out of line, costing 1-4% more
// instructions per delete/extract than the per-store loops it replaced.
#[inline]
pub fn partition_with<C: ?Sized, Id>(
    ctx: &mut C,
    start: Id,
    end: Id,
    mut keep: impl FnMut(&C, Id) -> bool,
    mut swap: impl FnMut(&mut C, Id, Id),
) -> Id
where
    Id: From<usize> + AddAssign + PartialOrd + SubAssign + Copy,
{
    let one = Id::from(1);
    let (mut left, mut right) = (start, end);
    while left < right {
        if keep(ctx, left) {
            left += one;
        } else {
            right -= one;
            if keep(ctx, right) {
                swap(ctx, left, right);
                left += one;
            }
        }
    }
    left
}
impl<E, V, H, N: NodeStorageOps<NodeData = V>> Swap<Hedge> for HedgeGraph<E, V, H, N> {
    fn len(&self) -> Hedge {
        self.hedge_data.len()
    }

    // type Item = H;

    // fn filter(&self, id: &Hedge, filter: &impl Fn(&Hedge, &Self::Item) -> bool) -> bool {
    //     filter(id, &self[*id])
    // }

    fn is_zero_length(&self) -> bool {
        self.hedge_data.is_zero_length()
    }

    fn swap(&mut self, i: Hedge, j: Hedge) {
        // println!("Swapping {i:?} with {j:?}");
        self.hedge_data.swap(i, j);
        // println!("nodeswap");
        self.node_store.swap(i, j);
        // println!("edgeswap");
        self.edge_store.swap(i, j);
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> Swap<EdgeIndex> for HedgeGraph<E, V, H, N> {
    // type Item = (E, HedgePair);

    // fn filter(&self, id: &EdgeIndex, filter: &impl Fn(&EdgeIndex, &Self::Item) -> bool) -> bool {
    //     <SmartEdgeVec<E> as Swap<EdgeIndex>>::filter(&self.edge_store, id, filter)
    // }

    fn swap(&mut self, i: EdgeIndex, j: EdgeIndex) {
        self.edge_store.swap(i, j);
    }

    fn is_zero_length(&self) -> bool {
        <SmartEdgeVec<E> as Swap<EdgeIndex>>::is_zero_length(&self.edge_store)
    }

    fn len(&self) -> EdgeIndex {
        self.edge_store.len()
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> Swap<NodeIndex> for HedgeGraph<E, V, H, N> {
    fn swap(&mut self, i: NodeIndex, j: NodeIndex) {
        self.node_store.swap(i, j);
    }

    // type Item = <N as Swap<NodeIndex>>::Item;

    // fn filter(&self, id: &NodeIndex, filter: &impl Fn(&NodeIndex, &Self::Item) -> bool) -> bool {
    //     <N as Swap<NodeIndex>>::filter(&self.node_store, id, filter)
    // }

    fn is_zero_length(&self) -> bool {
        <N as Swap<NodeIndex>>::is_zero_length(&self.node_store)
    }

    fn len(&self) -> NodeIndex {
        self.node_store.len()
    }
}

#[cfg(test)]
mod test {
    use super::{partition_with, Swap};
    use crate::{
        dot,
        half_edge::{NodeIndex, NodeVec},
        parser::DotGraph,
    };

    #[test]
    fn swap() {
        let graph: DotGraph = dot!(digraph {
            edge [label = "test"]
            node [label = "test"]
            in [style=invis]
            A
            B
            in -> A:3
            in -> A
            in -> A
            2->3 [dir=back]
            3->2
            in -> A [label = "override"]
            A -> in
        })
        .unwrap();

        println!("{}", graph.dot_of(&graph.full_filter()));
        println!("{}", graph.debug_dot())
    }

    /// Partitions `keep` from `start`, returning the split, the swaps and the result.
    fn run(keep: &[bool], start: usize) -> (usize, Vec<(usize, usize)>, Vec<bool>) {
        let mut values = keep.to_vec();
        let mut swaps = vec![];
        let split = partition_with(
            &mut values,
            start,
            keep.len(),
            |v, i| v[i],
            |v, i, j| {
                v.swap(i, j);
                swaps.push((i, j));
            },
        );
        (split, swaps, values)
    }

    #[test]
    fn kernel_edge_cases() {
        assert_eq!(run(&[], 0), (0, vec![], vec![]));
        assert_eq!(run(&[true; 3], 0), (3, vec![], vec![true; 3]));
        assert_eq!(run(&[false; 3], 0), (0, vec![], vec![false; 3]));
        assert_eq!(run(&[false], 1), (1, vec![], vec![false]));
    }

    #[test]
    fn kernel_fills_holes_from_the_tail() {
        let (t, f) = (true, false);
        assert_eq!(
            run(&[f, t, f, t, t], 0),
            (3, vec![(0, 4), (2, 3)], vec![t, t, t, f, f])
        );
        // The tail is scanned past rejected positions before swapping.
        assert_eq!(
            run(&[f, t, t, f, f], 0),
            (2, vec![(0, 2)], vec![t, t, f, f, f])
        );
    }

    #[test]
    fn kernel_leaves_positions_before_start() {
        let (t, f) = (true, false);
        assert_eq!(
            run(&[f, f, t, f, t], 2),
            (4, vec![(3, 4)], vec![f, f, t, t, f])
        );
    }

    /// Every keep pattern up to length 8 from every start: the result is a
    /// partition, and a predicate on the original positions agrees with one on
    /// the current values.
    #[test]
    fn kernel_partitions_and_reads_positions_before_moving_them() {
        for n in 0..=8 {
            for mask in 0..1u32 << n {
                let keep: Vec<bool> = (0..n).map(|i| mask >> i & 1 == 1).collect();
                for start in 0..=n {
                    let (split, swaps, values) = run(&keep, start);
                    assert_eq!(values[..start], keep[..start]);
                    assert!(values[start..split].iter().all(|&k| k));
                    assert!(values[split..].iter().all(|&k| !k));

                    let mut positional = keep.clone();
                    let mut positional_swaps = vec![];
                    let positional_split = partition_with(
                        &mut positional,
                        start,
                        n,
                        |_, i| keep[i],
                        |v, i, j| {
                            v.swap(i, j);
                            positional_swaps.push((i, j));
                        },
                    );
                    assert_eq!((positional_split, positional_swaps), (split, swaps));
                }
            }
        }
    }

    /// The three-way node split used by delete and extract: untouched, split, inside.
    #[test]
    fn partition_from_splits_three_ways() {
        let mut nodes: NodeVec<u8> = vec![2, 0, 1, 2, 0, 1].into();
        let split = nodes.partition_from(NodeIndex(0), |v, n| v[n] == 0);
        let inside = nodes.partition_from(split, |v, n| v[n] == 1);
        assert_eq!((split, inside), (NodeIndex(2), NodeIndex(4)));
        let order: Vec<u8> = nodes.into_iter().map(|(_, v)| v).collect();
        assert_eq!(order, [0, 0, 1, 1, 2, 2]);
    }
}
