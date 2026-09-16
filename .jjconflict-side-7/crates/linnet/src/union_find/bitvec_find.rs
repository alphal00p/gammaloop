use std::ops::{Index, IndexMut};

use bitvec::vec::BitVec;

use crate::half_edge::{
    involution::Hedge,
    nodestore::{BitVecNeighborIter, NodeStorage},
    NodeIndex,
};

use super::{ParentPointer, SetIndex, UnionFind};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LightIndex(pub usize);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HeavyIndex(pub usize);

pub type HLindex = HeavyLight<HeavyIndex, LightIndex>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HeavyLight<H, L> {
    Light(L),
    Heavy(H),
}

pub struct BitFilterData<U> {
    pub data: Option<U>,
    pub filter: BitVec,
}

pub struct UnionFindBitFilterHL<H, L> {
    pub inner: UnionFind<HLindex>,
    pub heavy_data: Vec<BitFilterData<H>>,
    pub light_data: Vec<BitFilterData<L>>,
}
impl<H, L> UnionFindBitFilterHL<H, L> {
    pub fn new_heavy(data: Vec<H>) -> Self {
        let n = data.len();
        // assert_eq!(n, data.len());

        let light_data = vec![];
        let mut heavy_data = vec![];
        let mut pointers = Vec::with_capacity(data.len());

        for (i, d) in data.into_iter().enumerate() {
            let mut filter: BitVec = BitVec::repeat(false, n);
            filter.set(i, true);
            let index = HeavyLight::Heavy(HeavyIndex(i));
            pointers.push(index);
            heavy_data.push(BitFilterData {
                data: Some(d),
                filter,
            });
        }

        let inner = UnionFind::new(pointers);

        Self {
            inner,
            heavy_data,
            light_data,
        }
    }

    pub fn new_light(data: Vec<L>) -> Self {
        let n = data.len();
        // assert_eq!(n, data.len());

        let mut light_data = vec![];
        let heavy_data = vec![];
        let mut pointers = Vec::with_capacity(n);

        for (i, d) in data.into_iter().enumerate() {
            let mut filter: BitVec = BitVec::repeat(false, n);
            filter.set(i, true);
            let index = HeavyLight::Light(LightIndex(i));
            pointers.push(index);
            light_data.push(BitFilterData {
                data: Some(d),
                filter,
            });
        }

        let inner = UnionFind::new(pointers);

        Self {
            inner,
            heavy_data,
            light_data,
        }
    }

    pub fn new(data_enum: Vec<HeavyLight<H, L>>) -> Self {
        let n = data_enum.len();
        // assert_eq!(n, data_enum.len());

        let mut light_data = vec![];
        let mut heavy_data = vec![];
        let mut pointers = Vec::with_capacity(n);

        for (i, d) in data_enum.into_iter().enumerate() {
            let mut filter: BitVec = BitVec::repeat(false, n);
            filter.set(i, true);
            let index = match d {
                HeavyLight::Heavy(h) => {
                    let id = HeavyLight::Heavy(HeavyIndex(heavy_data.len()));
                    heavy_data.push(BitFilterData {
                        data: Some(h),
                        filter: filter.into(),
                    });
                    id
                }
                HeavyLight::Light(l) => {
                    let id = HeavyLight::Light(LightIndex(light_data.len()));
                    light_data.push(BitFilterData {
                        data: Some(l),
                        filter: filter.into(),
                    });
                    id
                }
            };

            pointers.push(index);
        }

        let inner = UnionFind::new(pointers);

        Self {
            inner,
            heavy_data,
            light_data,
        }
    }

    /// Finds the representative of the set containing the element at ParentPointer x.
    pub fn find(&self, x: ParentPointer) -> ParentPointer {
        self.inner.find(x)
    }

    pub fn find_from_heavy(&self, h: HeavyIndex) -> ParentPointer {
        self.find(Hedge(self[&h].filter.iter_ones().next().unwrap()))
    }

    pub fn find_from_light(&self, l: LightIndex) -> ParentPointer {
        self.find(Hedge(self[&l].filter.iter_ones().next().unwrap()))
    }

    /// Returns a reference to the BitVec filter for the set containing the element at ParentPointer x.
    pub fn find_index(&self, x: ParentPointer) -> &HLindex {
        self.inner.find_data(x)
    }

    pub fn get(&self, hl_pointer: HLindex) -> HeavyLight<&H, &L> {
        match hl_pointer {
            HeavyLight::Heavy(h) => HeavyLight::Heavy(&self[h]),
            HeavyLight::Light(l) => HeavyLight::Light(&self[l]),
        }
    }

    pub fn find_data(&self, x: ParentPointer) -> HeavyLight<&H, &L> {
        let ptr = self.find_index(x);
        self.get(*ptr)
    }

    pub fn union<FH, FL>(
        &mut self,
        x: ParentPointer,
        y: ParentPointer,
        merge_h: FH,
        merge_l: FL,
    ) -> ParentPointer
    where
        FL: FnOnce(L, L) -> L,
        FH: FnOnce(H, H) -> H,
    {
        let mut loser_idx = None;

        let closure = |winner: HLindex, loser: HLindex| {
            match (winner, loser) {
                //-------------------------------------------------------
                // Case 1) Heavy vs. Heavy
                //-------------------------------------------------------
                (HeavyLight::Heavy(wi), HeavyLight::Heavy(li)) => {
                    let widx = wi.0;
                    let lidx = li.0;

                    // 1) Merge data
                    let wval = self.heavy_data[widx]
                        .data
                        .take()
                        .expect("Winner missing heavy data?");
                    let lval = self.heavy_data[lidx]
                        .data
                        .take()
                        .expect("Loser missing heavy data?");

                    self.heavy_data[widx].data = Some(merge_h(wval, lval));
                    if widx < lidx {
                        let (wslice, lslice) = self.heavy_data.split_at_mut(lidx);
                        wslice[widx].filter |= &lslice[0].filter;
                    } else {
                        let (lslice, wslice) = self.heavy_data.split_at_mut(widx);
                        wslice[0].filter |= &lslice[lidx].filter;
                    }
                }

                //-------------------------------------------------------
                // Case 2) Light vs. Light
                //-------------------------------------------------------
                (HeavyLight::Light(wi), HeavyLight::Light(li)) => {
                    let widx = wi.0;
                    let lidx = li.0;

                    // 1) Merge data
                    let wval = self.light_data[widx]
                        .data
                        .take()
                        .expect("Winner missing light data?");
                    let lval = self.light_data[lidx]
                        .data
                        .take()
                        .expect("Loser missing light data?");
                    self.light_data[widx].data = Some(merge_l(wval, lval));

                    // 2) Merge filters
                    if widx < lidx {
                        let (wslice, lslice) = self.light_data.split_at_mut(lidx);
                        wslice[widx].filter |= &lslice[0].filter;
                    } else {
                        let (lslice, wslice) = self.light_data.split_at_mut(widx);
                        wslice[0].filter |= &lslice[lidx].filter;
                    }
                }

                //-------------------------------------------------------
                // Case 3) Mismatch (can't unify a Heavy with a Light)
                //-------------------------------------------------------
                _ => panic!("Cannot unify a Heavy root with a Light root!"),
            }
            loser_idx = Some(loser);
            winner
        };
        let out = self.inner.union(x, y, closure);
        match loser_idx {
            Some(HeavyLight::Heavy(l)) => {
                if l.0 + 1 < self.heavy_data.len() {
                    self.heavy_data.swap_remove(l.0);
                    let set_index = self.inner.find_data_index(self.find_from_heavy(l));
                    self.inner[set_index] = HeavyLight::Heavy(l);
                } else {
                    self.heavy_data.pop();
                }
            }
            Some(HeavyLight::Light(l)) => {
                if l.0 + 1 < self.light_data.len() {
                    self.light_data.swap_remove(l.0);
                    let set_index = self.inner.find_data_index(self.find_from_light(l));
                    self.inner[set_index] = HeavyLight::Light(l);
                } else {
                    self.light_data.pop();
                }
            }
            _ => panic!("Cannot unify a Heavy root with a Light root!"),
        };
        out
    }

    // fn root_hedge(&self, node: NodeIndex) -> Hedge {
    //     self.inner[&SetIndex(node.0)].root_pointer
    // }
}

impl<H, L> Index<HeavyIndex> for UnionFindBitFilterHL<H, L> {
    type Output = H;
    fn index(&self, index: HeavyIndex) -> &Self::Output {
        self.heavy_data[index.0]
            .data
            .as_ref()
            .expect("missing heavy?")
    }
}

impl<H, L> IndexMut<HeavyIndex> for UnionFindBitFilterHL<H, L> {
    fn index_mut(&mut self, index: HeavyIndex) -> &mut Self::Output {
        self.heavy_data[index.0]
            .data
            .as_mut()
            .expect("missing heavy?")
    }
}

impl<H, L> Index<&HeavyIndex> for UnionFindBitFilterHL<H, L> {
    type Output = BitFilterData<H>;
    fn index(&self, index: &HeavyIndex) -> &Self::Output {
        &self.heavy_data[index.0]
    }
}

impl<H, L> IndexMut<&HeavyIndex> for UnionFindBitFilterHL<H, L> {
    fn index_mut(&mut self, index: &HeavyIndex) -> &mut Self::Output {
        &mut self.heavy_data[index.0]
    }
}

impl<H, L> Index<LightIndex> for UnionFindBitFilterHL<H, L> {
    type Output = L;
    fn index(&self, index: LightIndex) -> &Self::Output {
        self.light_data[index.0]
            .data
            .as_ref()
            .expect("missing light?")
    }
}

impl<H, L> IndexMut<LightIndex> for UnionFindBitFilterHL<H, L> {
    fn index_mut(&mut self, index: LightIndex) -> &mut Self::Output {
        self.light_data[index.0]
            .data
            .as_mut()
            .expect("missing light?")
    }
}

impl<H, L> Index<&LightIndex> for UnionFindBitFilterHL<H, L> {
    type Output = BitFilterData<L>;
    fn index(&self, index: &LightIndex) -> &Self::Output {
        &self.light_data[index.0]
    }
}

impl<H, L> IndexMut<&LightIndex> for UnionFindBitFilterHL<H, L> {
    fn index_mut(&mut self, index: &LightIndex) -> &mut Self::Output {
        &mut self.light_data[index.0]
    }
}

impl<H, L> NodeStorage for UnionFindBitFilterHL<H, L> {
    type Storage<N> = UnionFindBitFilterHL<N, L>;
    type NodeData = H;
    type Neighbors = BitVec;
    type NeighborsIter<'a>
        = BitVecNeighborIter<'a>
    where
        Self: 'a;
}

// impl<H, L> NodeStorageOps for UnionFindBitFilterHL<H, L> {
//     type OpStorage<A> = Self::Storage<A>;
//     fn hedge_len(&self) -> usize {
//         self.inner.n_elements()
//     }

//     fn identify_nodes(
//         &mut self,
//         nodes: &[NodeIndex],
//         node_data_merge: Self::NodeData,
//     ) -> NodeIndex {
//         let hedges = nodes
//             .iter()
//             .map(|n| self.root_hedge(*n))
//             .collect::<Vec<_>>();

//         let n_init_root = hedges[0];
//         for n in hedges.iter().skip(1) {
//             self.inner.union(n_init_root, *n, |a, _| a);
//         }
//         self.inner.replace_set_data_of(n_init_root, |mut a| {
//             a.data = node_data_merge;
//             a
//         });
//         NodeIndex(self.nodes.find_data_index(n_init_root).0)
//     }

//     fn to_forest<U>(
//         &self,
//         map_data: impl Fn(&Self::NodeData) -> U,
//     ) -> crate::tree::Forest<U, crate::tree::parent_pointer::ParentPointerStore<()>> {
//         Forest {
//             roots: self
//                 .nodes
//                 .set_data
//                 .iter()
//                 .map(|a| RootData {
//                     root_id: a.root_pointer.into(),
//                     data: map_data(&a.data.as_ref().unwrap().data),
//                 })
//                 .collect(),
//             nodes: self
//                 .nodes
//                 .nodes
//                 .iter()
//                 .map(|a| {
//                     let ad = a.get();
//                     let n = match &ad {
//                         super::UFNode::Child(c) => PPNode::child((), (*c).into()),
//                         super::UFNode::Root { set_data_idx, .. } => {
//                             PPNode::root((), RootId(set_data_idx.0))
//                         }
//                     };

//                     a.set(ad);

//                     n
//                 })
//                 .collect(),
//         }
//     }

//     fn node_len(&self) -> usize {
//         self.nodes.n_sets()
//     }

//     fn set_hedge_data(&mut self, _hedge: crate::half_edge::involution::Hedge, _nodeid: NodeIndex) {
//         panic!("should not need to set")
//     }

//     fn check_and_set_nodes(&mut self) -> Result<(), crate::half_edge::HedgeGraphError> {
//         Ok(())
//     }

//     fn map_data_ref_graph<'a, E, V2>(
//         &'a self,
//         graph: &'a crate::half_edge::HedgeGraph<E, Self::NodeData, Self>,
//         mut node_map: impl FnMut(
//             &'a crate::half_edge::HedgeGraph<E, Self::NodeData, Self>,
//             &'a HedgeNode,
//             &'a Self::NodeData,
//         ) -> V2,
//     ) -> Self::Storage<V2> {
//         UnionFindNodeStore {
//             nodes: self.nodes.map_set_data_ref(|n| HedgeNodeStore {
//                 data: node_map(graph, &n.node, &n.data),
//                 node: n.node.clone(),
//             }),
//         }
//     }

//     fn get_node(&self, node_id: NodeIndex) -> &HedgeNode {
//         &self.nodes[SetIndex(node_id.0)].node
//     }

//     fn iter_nodes(&self) -> impl Iterator<Item = (&HedgeNode, &Self::NodeData)> {
//         self.nodes.iter_set_data().map(|(_, d)| (&d.node, &d.data))
//     }

//     fn iter_nodes_mut(&mut self) -> impl Iterator<Item = (&HedgeNode, &mut Self::NodeData)> {
//         self.nodes
//             .iter_set_data_mut()
//             .map(|(_, d)| (&d.node, &mut d.data))
//     }

//     fn node_id_ref(&self, hedge: crate::half_edge::involution::Hedge) -> NodeIndex {
//         NodeIndex(self.nodes.find_data_index(hedge).0)
//     }

//     fn get_node_data_mut(&mut self, node_id: NodeIndex) -> &mut Self::NodeData {
//         &mut self.nodes[SetIndex(node_id.0)].data
//     }

//     fn iter_node_id(&self) -> impl Iterator<Item = NodeIndex> {
//         self.nodes.iter_set_data().map(|(i, _)| NodeIndex(i.0))
//     }

//     fn get_node_data(&self, node_id: NodeIndex) -> &Self::NodeData {
//         &self.nodes[SetIndex(node_id.0)].data
//     }

//     fn map_data_graph<'a, V2>(
//         self,
//         involution: &crate::half_edge::involution::Involution<
//             crate::half_edge::involution::EdgeIndex,
//         >,
//         mut f: impl FnMut(
//             &crate::half_edge::involution::Involution<crate::half_edge::involution::EdgeIndex>,
//             &HedgeNode,
//             NodeIndex,
//             Self::NodeData,
//         ) -> V2,
//     ) -> Self::Storage<V2> {
//         UnionFindNodeStore {
//             nodes: self.nodes.map_set_data(|i, n| HedgeNodeStore {
//                 data: f(involution, &n.node, NodeIndex(i.0), n.data),
//                 node: n.node.clone(),
//             }),
//         }
//     }

//     fn extend(mut self, other: Self) -> Self {
//         self.nodes.extend(other.nodes);
//         self
//     }

//     fn iter(&self) -> impl Iterator<Item = (NodeIndex, &Self::NodeData)> {
//         self.nodes
//             .iter_set_data()
//             .map(|(i, d)| (NodeIndex(i.0), &d.data))
//     }

//     fn build<
//         I: IntoIterator<Item = crate::half_edge::builder::HedgeNodeBuilder<Self::NodeData>>,
//     >(
//         nodes: I,
//         n_hedges: usize,
//     ) -> Self {
//         NodeStorageVec::build(nodes, n_hedges).into()
//     }

//     fn drain(self) -> impl Iterator<Item = (NodeIndex, Self::NodeData)> {
//         self.nodes
//             .drain_set_data()
//             .map(|(s, d)| (NodeIndex(s.0), d.data))
//     }

//     fn random(sources: &[HedgeNode], sinks: &[HedgeNode]) -> Self
//     where
//         Self::NodeData: Default,
//     {
//         NodeStorageVec::random(sources, sinks).into()
//     }

//     fn add_dangling_edge(
//         mut self,
//         source: NodeIndex,
//     ) -> Result<Self, crate::half_edge::HedgeGraphError> {
//         let setid = SetIndex(source.0);
//         let _ = self.nodes.add_child(setid);
//         self.nodes.iter_set_data_mut().for_each(|(s, n)| {
//             if s == setid {
//                 n.node.hairs.push(true);
//             } else {
//                 n.node.hairs.push(false);
//             }
//             n.node.internal_graph.filter.push(false);
//         });
//         Ok(self)
//     }
// }
