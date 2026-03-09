use std::ops::{AddAssign, SubAssign};

use crate::permutation::Permutation;

use super::{
    hedgevec::SmartEdgeVec,
    involution::{EdgeIndex, Hedge},
    nodestore::NodeStorageOps,
    HedgeGraph, NodeIndex,
};

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

    fn partition(&mut self, filter: impl Fn(&Id) -> bool) -> Id
    where
        Id: From<usize> + AddAssign + PartialOrd + SubAssign + Copy,
    {
        let mut left = Id::from(0);
        let one = Id::from(1);
        let mut extracted = self.len();

        while left < extracted {
            if filter(&left) {
                //left is in the right place
                left += one;
            } else {
                //left needs to be swapped
                extracted -= one;
                if filter(&extracted) {
                    // println!("{extracted}<=>{left}");
                    //only with an extracted that is in the wrong spot
                    self.swap(left, extracted);
                    left += one;
                }
            }
        }
        left
    }
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
    use crate::{
        dot,
        half_edge::nodestore::NodeStorageVec,
        parser::{DotGraph, DotVertexData},
    };

    #[test]
    fn swap() {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = dot!(digraph {
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
}
