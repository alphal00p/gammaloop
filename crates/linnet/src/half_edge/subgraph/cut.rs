use std::ops::{Range, RangeFrom, RangeInclusive, RangeTo, RangeToInclusive};

use super::{Cycle, Inclusion, SubSetIter, SubSetLike, SubSetOps};
use crate::half_edge::hedgevec::Accessors;
use crate::half_edge::involution::{EdgeIndex, HedgePair};
use crate::half_edge::nodestore::NodeStorageOps;
use crate::half_edge::subgraph::{ModifySubSet, SuBitGraph, SubGraphLike};
use crate::half_edge::EdgeAccessors;
use crate::half_edge::{
    involution::SignOrZero, EdgeData, Flow, Hedge, HedgeGraph, InvolutiveMapping, NodeStorage,
    Orientation, PowersetIterator,
};
use crate::parser::DotEdgeData;
use std::cmp::Ordering;
use std::{
    fmt::{Display, Formatter},
    hash::Hash,
};
use thiserror::Error;

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents an oriented cut in a graph.
///
/// A cut is a partition of the graph's half-edges into two sets, typically
/// defining a boundary. This `OrientedCut` explicitly stores half-edges
/// considered to be on the "left" and "right" sides of this boundary.
/// The orientation implies a direction across the cut.
///
/// For a given full edge (composed of two half-edges, say `h1` and `h2` where
/// `h2 = inv(h1)`), if `h1` is in `left`, then `h2` should be in `right`,
/// and vice-versa, for the cut to be consistent along that edge.
pub struct OrientedCut {
    /// A bitmask representing the set of half-edges on the "left" side of the cut.
    /// In a scattering diagram context, these might be the edges one would drag
    /// to the right to pass through the cut.
    pub left: SuBitGraph,
    /// A bitmask representing the set of half-edges on the "right" side of the cut.
    /// These are typically theinvolutional pairs of the half-edges in the `left` set.
    /// In a scattering diagram context, these might be the edges one would drag
    /// to the left.
    ///
    /// It's generally expected that `right` contains `inv(h)` for every `h` in `left`.
    pub right: SuBitGraph,
}

impl Display for OrientedCut {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for h in self.left.union(&self.right).included_iter() {
            if self.left.includes(&h) {
                write!(f, "+{}", h.0)?;
            }
            if self.right.includes(&h) {
                write!(f, "-{}", h.0)?;
            }
        }

        Ok(())
    }
}

/// Errors that can occur during operations related to graph cuts.
#[derive(Debug, Error)]
pub enum CutError {
    #[error("Invalid edge: The specified edge is not suitable for the cut operation.")]
    /// Indicates an edge is invalid for the context of a cut operation.
    InvalidEdge,
    #[error("Invalid orientation: The orientation specified or found is not valid for the cut operation.")]
    /// Indicates an orientation is invalid in the context of a cut.
    InvalidOrientation,
    #[error("Cut edge has already been set: The edge is already part of the cut in an incompatible way.")]
    /// Attempted to define a cut edge that was already defined, or its inverse was.
    CutEdgeAlreadySet,
    #[error("Cut edge is identity: An unpaired (identity) half-edge cannot form part of a separating cut in this context.")]
    /// An attempt was made to use an identity (unpaired) half-edge in a cut where a paired edge was expected.
    CutEdgeIsIdentity,
}

impl OrientedCut {
    /// disregards identity edges
    pub fn from_underlying_coerce<E, V, H, N: NodeStorageOps<NodeData = V>>(
        cut: SuBitGraph,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Result<Self, CutError> {
        let mut right = graph.empty_subgraph::<SuBitGraph>();

        for i in cut.included_iter() {
            let invh = graph.inv(i);
            if cut.includes(&invh) {
                return Err(CutError::CutEdgeAlreadySet);
            } else if invh == i {
                right.sub(i);
            }
            right.add(invh);
        }

        cut.subtract(&right);
        Ok(OrientedCut { left: cut, right })
    }

    /// Errors for identity edges
    pub fn from_underlying_strict<E, V, H, N: NodeStorageOps<NodeData = V>>(
        cut: SuBitGraph,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Result<Self, CutError> {
        let mut right = graph.empty_subgraph::<SuBitGraph>();

        for i in cut.included_iter() {
            let invh = graph.inv(i);
            if cut.includes(&invh) {
                return Err(CutError::CutEdgeAlreadySet);
            } else if invh == i {
                return Err(CutError::CutEdgeIsIdentity);
            }
            right.add(invh);
        }
        Ok(OrientedCut { left: cut, right })
    }

    pub fn iter_edges_flow<'a, E, V, H, N: NodeStorageOps<NodeData = V>>(
        &'a self,
        graph: &'a HedgeGraph<E, V, H, N>,
    ) -> impl Iterator<Item = (Flow, EdgeIndex, EdgeData<&'a E>)> {
        graph.iter_edges_of(&self.left).filter_map(|(pair, b, c)| {
            if let HedgePair::Split { split, .. } = pair {
                Some((split, b, c))
            } else {
                None
            }
        })
    }

    /// takes all non-cut edges and gives all possible signs to them.
    ///
    /// If C is a set of edges (pairs of half-edges), then take all S in Pset(C), and put them to the left of the cut. All other edges are put to the right.
    pub fn all_initial_state_cuts<E, V, H, N: NodeStorageOps<NodeData = V>>(
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Vec<Self> {
        let mut all_cuts = Vec::new();

        for c in graph.non_cut_edges() {
            if c.is_empty() {
                continue;
            }
            let mut all_sources = graph.empty_subgraph::<SuBitGraph>();

            for h in c.included_iter() {
                match graph.edge_store.inv_full(h) {
                    InvolutiveMapping::Identity { .. } => {
                        panic!("cut edge is identity")
                    }
                    InvolutiveMapping::Source { .. } => {
                        all_sources.add(h);
                    }
                    InvolutiveMapping::Sink { .. } => {}
                }
            }

            let n_cut_edges: u8 = all_sources.n_included().try_into().unwrap();

            let pset = PowersetIterator::new(n_cut_edges); //.unchecked_sub(1)

            for i in pset {
                let mut left = graph.empty_subgraph::<SuBitGraph>();
                for (j, h) in all_sources.included_iter().enumerate() {
                    // if let Some(j) = j.checked_sub(1) {
                    if i[Hedge(j)] {
                        left.add(h);
                    } else {
                        left.add(graph.inv(h));
                    }
                }
                all_cuts.push(Self::from_underlying_strict(left, graph).unwrap());
            }
        }

        all_cuts
    }

    pub fn winding_number(&self, cycle: &Cycle) -> i32 {
        let mut winding_number = 0;

        for h in cycle.filter.included_iter() {
            winding_number += SignOrZero::from(self.relative_orientation(h)) * 1;
        }

        winding_number
    }
    pub fn iter_edges<'a, E, V, H, N: NodeStorageOps<NodeData = V>>(
        &'a self,
        graph: &'a HedgeGraph<E, V, H, N>,
    ) -> impl Iterator<Item = (Orientation, EdgeData<&'a E>)> {
        self.left
            .included_iter()
            .map(|i| (self.orientation(i, graph), graph.get_edge_data_full(i)))
    }

    pub fn iter_edges_idx_relative(&self) -> impl Iterator<Item = (Hedge, Orientation)> + '_ {
        self.left
            .included_iter()
            .map(|i| (i, self.relative_orientation(i)))
    }

    pub fn iter_left_hedges(&self) -> impl Iterator<Item = Hedge> + '_ {
        self.left.included_iter()
    }

    pub fn iter_right_hedges(&self) -> impl Iterator<Item = Hedge> + '_ {
        self.right.included_iter()
    }

    /// - If the left and right are aligned=> panic
    /// - If in left set => Default
    /// - If in right set => Reversed
    /// - If in neither => Undirected
    pub fn relative_orientation(&self, i: Hedge) -> Orientation {
        match (self.left.includes(&i), self.right.includes(&i)) {
            (true, true) => panic!("Both left and right are included in the reference"),
            (true, false) => Orientation::Default,
            (false, true) => Orientation::Reversed,
            (false, false) => Orientation::Undirected,
        }
    }

    /// essentially tells you if the edge is in the cut or not ([Orientation::Undirected])
    ///
    /// If it is in the cut, and the left set contains the source hedge [Orientation::Default] or the sink hedge [Orientation::Reversed].
    ///
    /// Equivalently tells you if the source hedge is in the left side of the cut [Orientation::Default] or the right side [Orientation::Reversed].
    pub fn get_from_pair(&self, pair: HedgePair) -> Orientation {
        match pair {
            HedgePair::Paired { source, sink } => {
                debug_assert!(
                    (self.left.includes(&source) && !self.left.includes(&sink))
                        || (!self.left.includes(&source) && self.left.includes(&sink))
                );
                if self.left.includes(&source) {
                    // debug_assert!(self.right.includes(&sink));
                    Orientation::Default
                } else if self.left.includes(&sink) {
                    // debug_assert!(self.right.includes(&source));
                    Orientation::Reversed
                } else {
                    Orientation::Undirected
                }
            }
            HedgePair::Split {
                source,
                sink,
                split,
            } => {
                debug_assert!(
                    (self.left.includes(&source) && !self.left.includes(&sink))
                        || (!self.left.includes(&source) && self.left.includes(&sink))
                );
                match split {
                    Flow::Sink => {
                        if self.left.includes(&sink) {
                            // debug_assert!(self.right.includes(&source));
                            Orientation::Reversed
                        } else {
                            Orientation::Undirected
                        }
                    }
                    Flow::Source => {
                        if self.left.includes(&source) {
                            // debug_assert!(self.right.includes(&sink));
                            Orientation::Default
                        } else {
                            Orientation::Undirected
                        }
                    }
                }
            }
            HedgePair::Unpaired { hedge, .. } => {
                if self.left.includes(&hedge) {
                    Orientation::Default
                } else if self.right.includes(&hedge) {
                    Orientation::Reversed
                } else {
                    Orientation::Undirected
                }
            }
        }
    }

    /// Set the left cut containing the [Flow] hedge.
    pub fn set(&mut self, pair: HedgePair, flow: Flow) {
        match pair {
            HedgePair::Paired { source, sink } => match flow {
                Flow::Source => {
                    self.left.add(source);
                    self.right.sub(source);
                    self.left.sub(sink);
                    self.right.add(sink);
                }
                Flow::Sink => {
                    self.left.add(sink);
                    self.right.sub(sink);
                    self.left.sub(source);
                    self.right.add(source);
                }
            },
            HedgePair::Unpaired { hedge, .. } => match flow {
                Flow::Source => {
                    self.left.add(hedge);
                    self.right.sub(hedge);
                }
                Flow::Sink => {
                    self.left.sub(hedge);
                    self.right.add(hedge);
                }
            },
            _ => {}
        }
    }

    /// essentially tells you if the edge is in the cut or not ([Orientation::Undirected])
    ///
    /// If it is in the cut, and the left set contains the source hedge [Orientation::Default] or the sink hedge [Orientation::Reversed].
    ///
    /// Equivalently tells you if the source hedge is in the left side of the cut [Orientation::Default] or the right side [Orientation::Reversed].
    pub fn orientation<E, V, H, N: NodeStorage<NodeData = V>>(
        &self,
        i: Hedge,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Orientation {
        let pair = graph.edge_store.pair(i);
        self.get_from_pair(pair)
    }

    pub fn cut_edge<E>(&self, data: E, pair: HedgePair, id: EdgeIndex) -> PossiblyCutEdge<E> {
        let mut edge = PossiblyCutEdge::uncut(data, id);

        let orientation = self.get_from_pair(pair);
        match orientation {
            Orientation::Default => {
                // println!("Source:{id:?}");
                edge.cut(Flow::Source);
            }
            Orientation::Reversed => {
                // println!("Sink:{id:?}");
                edge.cut(Flow::Sink);
            }
            Orientation::Undirected => {}
        }

        edge
    }
}

impl Inclusion<Hedge> for OrientedCut {
    fn includes(&self, other: &Hedge) -> bool {
        self.left.includes(other)
    }
    fn intersects(&self, other: &Hedge) -> bool {
        self.left.intersects(other)
    }
}

impl Inclusion<HedgePair> for OrientedCut {
    fn includes(&self, other: &HedgePair) -> bool {
        self.left.includes(other)
    }
    fn intersects(&self, other: &HedgePair) -> bool {
        self.left.intersects(other)
    }
}

impl Inclusion<SuBitGraph> for OrientedCut {
    fn includes(&self, other: &SuBitGraph) -> bool {
        self.left.includes(other)
    }

    fn intersects(&self, other: &SuBitGraph) -> bool {
        self.left.intersects(other)
    }
}

impl Inclusion<OrientedCut> for OrientedCut {
    fn includes(&self, other: &OrientedCut) -> bool {
        self.left.includes(&other.left)
    }

    fn intersects(&self, other: &OrientedCut) -> bool {
        self.left.intersects(&other.left)
    }
}
impl Inclusion<Range<Hedge>> for OrientedCut {
    fn includes(&self, other: &Range<Hedge>) -> bool {
        (other.start.0..other.end.0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &Range<Hedge>) -> bool {
        (other.start.0..other.end.0).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeTo<Hedge>> for OrientedCut {
    fn includes(&self, other: &RangeTo<Hedge>) -> bool {
        (0..other.end.0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeTo<Hedge>) -> bool {
        (0..other.end.0).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeToInclusive<Hedge>> for OrientedCut {
    fn includes(&self, other: &RangeToInclusive<Hedge>) -> bool {
        (0..=other.end.0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeToInclusive<Hedge>) -> bool {
        (0..=other.end.0).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeFrom<Hedge>> for OrientedCut {
    fn includes(&self, other: &RangeFrom<Hedge>) -> bool {
        (other.start.0..).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeFrom<Hedge>) -> bool {
        (other.start.0..).any(|a| self.includes(&Hedge(a)))
    }
}

impl Inclusion<RangeInclusive<Hedge>> for OrientedCut {
    fn includes(&self, other: &RangeInclusive<Hedge>) -> bool {
        (other.start().0..=other.end().0).all(|a| self.includes(&Hedge(a)))
    }

    fn intersects(&self, other: &RangeInclusive<Hedge>) -> bool {
        (other.start().0..=other.end().0).any(|a| self.includes(&Hedge(a)))
    }
}

impl SubGraphLike for OrientedCut {
    fn nedges<E, V, H, N: NodeStorage<NodeData = V>>(
        &self,
        _graph: &HedgeGraph<E, V, H, N>,
    ) -> usize {
        self.n_included()
    }
    fn hairs(&self, node: impl Iterator<Item = Hedge>) -> SuBitGraph {
        self.left.hairs(node)
    }
}

impl SubSetLike for OrientedCut {
    type Base = SuBitGraph;

    type BaseIter<'a> = SubSetIter<'a>;

    fn size(&self) -> usize {
        self.left.size()
    }

    fn has_greater(&self, hedge: Hedge) -> bool {
        self.left.has_greater(hedge) || self.right.has_greater(hedge)
    }

    fn join_mut(&mut self, other: Self) {
        self.left.join_mut(other.left);
        self.right.join_mut(other.right);
    }

    fn included(&self) -> &SuBitGraph {
        self.left.included()
    }

    fn included_iter(&self) -> Self::BaseIter<'_> {
        self.left.included_iter()
    }

    fn n_included(&self) -> usize {
        self.left.n_included()
    }
    fn empty(size: usize) -> Self {
        OrientedCut {
            left: SuBitGraph::empty(size),
            right: SuBitGraph::empty(size),
        }
    }

    fn is_empty(&self) -> bool {
        self.left.is_empty()
    }

    fn string_label(&self) -> String {
        self.left.string_label()
    }
    fn from_base62(label: &str, size: usize) -> Option<Self> {
        None
    }
}

impl OrientedCut {
    /// Take the graph and split it along the cut, putting the cut orientation, and original edge index as additional data.
    pub fn to_owned_graph_ref<E, V, H, N: NodeStorageOps<NodeData = V>>(
        self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> HedgeGraph<PossiblyCutEdge<&E>, &V, &H, N::OpStorage<&V>> {
        let mut new_graph = graph.map_data_ref(
            |_, _, v| v,
            |_, i, _, e| e.map(|d| PossiblyCutEdge::uncut(d, i)),
            |_, h| h,
        );
        for h in self.iter_left_hedges() {
            new_graph[[&h]].cut(Flow::Source);
            let data = EdgeData::new(new_graph[[&h]].reverse(), new_graph.orientation(h));
            let invh = new_graph.inv(h);
            new_graph.split_edge(invh, data).unwrap();
        }

        new_graph
    }

    /// Take the graph and split it along the cut, putting the cut orientation, and original edge index as additional data.
    pub fn to_owned_graph<E, V, H, N: NodeStorageOps<NodeData = V>>(
        self,
        graph: HedgeGraph<E, V, H, N>,
    ) -> HedgeGraph<PossiblyCutEdge<E>, V, H, N::OpStorage<V>> {
        let mut new_graph = graph.map(
            |_, _, v| v,
            |i, _, h, _, e| {
                e.map(|d| {
                    let h = i[h.any_hedge()];
                    PossiblyCutEdge::uncut(d, h)
                })
            },
            |_, h| h,
        );
        for h in self.iter_left_hedges() {
            new_graph[[&h]].cut(Flow::Source);
            let data = EdgeData::new(
                new_graph[[&h]].duplicate_without_data().reverse(),
                new_graph.orientation(h),
            );
            let invh = new_graph.inv(h);
            new_graph.split_edge(invh, data).unwrap();
        }

        new_graph
    }
}

pub type CutGraph<E, V, H, N> = HedgeGraph<PossiblyCutEdge<E>, V, H, N>;

impl<E, V, H, N: NodeStorageOps<NodeData = V>> CutGraph<E, V, H, N> {
    pub fn split(&mut self) {
        let cut = self.cut();
        for h in cut.iter_left_hedges() {
            self[[&h]].cut(Flow::Source);
            let data = EdgeData::new(
                self[[&h]].duplicate_without_data().reverse(),
                self.orientation(h),
            );
            let invh = self.inv(h);
            self.split_edge(invh, data).unwrap();
        }
    }

    pub fn split_clone(&mut self)
    where
        E: Clone,
    {
        let cut = self.cut();
        for h in cut.iter_left_hedges() {
            self[[&h]].cut(Flow::Source);
            let data = EdgeData::new(self[[&h]].clone().reverse(), self.orientation(h));
            let invh = self.inv(h);
            self.split_edge(invh, data).unwrap();
        }
    }

    pub fn split_copy(&mut self)
    where
        E: Copy,
    {
        let cut = self.cut();
        for h in cut.iter_left_hedges() {
            self[[&h]].cut(Flow::Source);
            let data = EdgeData::new(self[[&h]].reverse(), self.orientation(h));
            let invh = self.inv(h);
            self.split_edge(invh, data).unwrap();
        }
    }

    pub fn round_trip_split(&mut self) {
        self.glue_back_strict();
        self.split();
    }

    pub fn round_trip_glue(&mut self) {
        self.split();
        self.glue_back_strict();
    }

    pub fn cut(&self) -> OrientedCut {
        let mut cut = OrientedCut::empty(self.n_hedges());
        for (h, _, d) in self.iter_edges() {
            match d.data.flow {
                Orientation::Default => cut.set(h, Flow::Source),
                Orientation::Reversed => cut.set(h, Flow::Sink),
                Orientation::Undirected => {}
            }
        }

        cut
    }

    pub fn debug_cut_dot(&self) -> String {
        self.dot_impl(
            &self.full_filter(),
            "",
            &|_| None,
            &|a| Some(format!("label=\"{}\"", a.label())),
            &|_| Some("".to_string()),
        )
    }

    /// glues_back a cut graph. only matches edges with the same cut data, reversed cut flow, and compatible orientation and flow.
    pub fn glue_back_strict(&mut self) {
        self.sew(
            |lf, ld, rf, rd| {
                if lf == -rf {
                    if ld.orientation == rd.orientation {
                        ld.data.matches(rd.data)
                    } else {
                        false
                    }
                } else {
                    false
                }
            },
            |lf, ld, rf, rd| {
                let lo: Flow = ld.data.flow.try_into().unwrap();
                let ro: Flow = rd.data.flow.try_into().unwrap();
                // println!("lf{lf:?},lo{lo:?},rf{rf:?},ro{ro:?}");
                debug_assert_eq!(lo, -ro);
                debug_assert_eq!(lf, -rf);

                let mut data = ld.data.merge(rd.data).unwrap();

                let orientation = ld.orientation;
                match (lf, lo) {
                    // A source hedge on the right of the cut
                    // Means that edge data needs to say Orientation::Reversed,
                    // so we need a relative difference in the flow
                    (Flow::Source, Flow::Sink) => {
                        data.cut(Flow::Sink);
                        (Flow::Source, EdgeData::new(data, orientation))
                    }
                    // A sink hedge on the right of the cut
                    // Means that edge data needs to say Orientation::Default
                    // so we need an alignment in the flow
                    (Flow::Sink, Flow::Source) => {
                        data.cut(Flow::Sink);
                        (Flow::Sink, EdgeData::new(data, orientation))
                    }
                    // A source hedge on the left of the cut
                    // Means that edge data needs to say Orientation::Default
                    // so we need an alignment in the flow
                    (Flow::Source, Flow::Source) => {
                        data.cut(Flow::Source);
                        (Flow::Source, EdgeData::new(data, orientation))
                    }
                    // A sink hedge on the left of the cut
                    // Means that edge data needs to say Orientation::Reversed
                    // so we need a relative difference in the flow
                    (Flow::Sink, Flow::Sink) => {
                        data.cut(Flow::Source);
                        (Flow::Sink, EdgeData::new(data, orientation))
                    }
                }
            },
        )
        .unwrap()
    }

    /// glues_back a cut graph. only matches edges with the same cut data, reversed cut flow, and compatible orientation and without flow matching.
    pub fn glue_back_lenient(&mut self) {
        self.sew(
            |lf, ld, rf, rd| {
                if lf == -rf {
                    if ld.orientation == rd.orientation {
                        ld.data.matches(rd.data)
                    } else {
                        false
                    }
                } else if ld.orientation == rd.orientation.reverse() {
                    ld.data.matches(rd.data)
                } else {
                    false
                }
            },
            |lf, ld, _, rd| {
                let lo: Flow = ld.data.flow.try_into().unwrap();
                let ro: Flow = rd.data.flow.try_into().unwrap();
                debug_assert_eq!(lo, -ro);
                // debug_assert_eq!(lf, -rf);

                let mut data = ld.data.merge(rd.data).unwrap();

                let orientation = ld.orientation;
                match (lf, lo) {
                    // A source hedge on the right of the cut
                    // Means that edge data needs to say Orientation::Reversed,
                    // so we need a relative difference in the flow
                    (Flow::Source, Flow::Sink) => {
                        data.cut(Flow::Sink);
                        (Flow::Source, EdgeData::new(data, orientation))
                    }
                    // A sink hedge on the right of the cut
                    // Means that edge data needs to say Orientation::Default
                    // so we need an alignment in the flow
                    (Flow::Sink, Flow::Source) => {
                        data.cut(Flow::Sink);
                        (Flow::Sink, EdgeData::new(data, orientation))
                    }
                    // A source hedge on the left of the cut
                    // Means that edge data needs to say Orientation::Default
                    // so we need an alignment in the flow
                    (Flow::Source, Flow::Source) => {
                        data.cut(Flow::Source);
                        (Flow::Source, EdgeData::new(data, orientation))
                    }
                    // A sink hedge on the left of the cut
                    // Means that edge data needs to say Orientation::Reversed
                    // so we need a relative difference in the flow
                    (Flow::Sink, Flow::Sink) => {
                        data.cut(Flow::Source);
                        (Flow::Sink, EdgeData::new(data, orientation))
                    }
                }
            },
        )
        .unwrap()
    }
}

#[derive(Debug, Eq, Clone, Copy)]
/// Represents an edge that may or may not be part of a cut, and if it is,
/// stores its orientation relative to that cut.
///
/// This struct is typically used as the edge data type `E` in a `HedgeGraph`
/// (often aliased as `CutGraph`) when the graph has been "split" along an
/// [`OrientedCut`]. It allows each edge to remember its original `EdgeIndex`
/// and how it relates to the cut.
///
/// # Type Parameters
///
/// - `E`: The original custom data type of the edge before being potentially marked as cut.
pub struct PossiblyCutEdge<E> {
    /// The original custom data of the edge. This is `None` if the edge was
    /// created as a conceptual part of a split without original data (e.g., the
    /// "other side" of a split external edge).
    data: Option<E>,
    /// The orientation of this edge relative to a cut.
    /// - [`Orientation::Undirected`]: The edge is not part of the cut.
    /// - [`Orientation::Default`]: The edge crosses the cut in one direction (e.g., "left to right", or "source side").
    /// - [`Orientation::Reversed`]: The edge crosses the cut in the opposite direction (e.g., "right to left", or "sink side").
    flow: Orientation,
    /// The original [`EdgeIndex`] of this edge in the graph before any cut operations.
    /// This helps in identifying the edge across different graph representations or
    /// after operations like splitting and gluing.
    pub index: EdgeIndex,
}
impl<E> From<PossiblyCutEdge<E>> for DotEdgeData
where
    DotEdgeData: From<E>,
{
    fn from(value: PossiblyCutEdge<E>) -> Self {
        let mut statements = DotEdgeData::empty();
        if let Some(data) = value.data {
            let data_statements = DotEdgeData::from(data);
            statements.extend(data_statements);
        }

        statements.add_statement(
            "cut_flow",
            match value.flow {
                Orientation::Default => "aligned",
                Orientation::Reversed => "reversed",
                Orientation::Undirected => "uncut",
            },
        );

        let edge_id: usize = value.index.into();
        statements.add_statement("edge_id", edge_id.to_string());

        statements
    }
}

impl<E: TryFrom<DotEdgeData>> TryFrom<DotEdgeData> for PossiblyCutEdge<E> {
    type Error = String;
    fn try_from(dot_edge_data: DotEdgeData) -> Result<Self, Self::Error> {
        let flow = dot_edge_data
            .statements
            .get("cut_flow")
            .ok_or("Missing 'cut_flow' attribute")?;

        let flow = match flow.as_str() {
            "aligned" => Orientation::Default,
            "reversed" => Orientation::Reversed,
            "uncut" => Orientation::Undirected,
            _ => return Err("Invalid 'cut_flow' value".to_string()),
        };

        let edge_id = dot_edge_data
            .statements
            .get("edge_id")
            .ok_or("Missing 'edge_id' attribute")?;
        let edge_id: usize = edge_id
            .parse()
            .map_err(|_| "Invalid 'edge_id' value".to_string())?;

        let data = dot_edge_data.try_into().ok();

        Ok(PossiblyCutEdge {
            data,
            flow,
            index: edge_id.into(),
        })
    }
}

impl<E: Hash> Hash for PossiblyCutEdge<E> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        (&self.data, self.flow).hash(state);
    }
}

impl<E: PartialOrd> PartialOrd for PossiblyCutEdge<E> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        (&self.data, self.flow).partial_cmp(&(&other.data, other.flow))
    }
}

impl<E: Ord> Ord for PossiblyCutEdge<E> {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.data, self.flow).cmp(&(&other.data, other.flow))
    }
}

impl<E: PartialEq> PartialEq for PossiblyCutEdge<E> {
    fn eq(&self, other: &Self) -> bool {
        (&self.data, self.flow).eq(&(&other.data, other.flow))
    }
}

impl<E> PossiblyCutEdge<E> {
    pub fn edge_data(&self) -> &E {
        self.data.as_ref().unwrap()
    }

    pub fn set_data(&mut self, data: E) {
        self.data = Some(data);
    }

    pub fn map<F, T>(self, f: F) -> PossiblyCutEdge<T>
    where
        F: FnOnce(E) -> T,
    {
        PossiblyCutEdge {
            data: self.data.map(f),
            flow: self.flow,
            index: self.index,
        }
    }

    pub fn as_ref(&self) -> PossiblyCutEdge<&E> {
        PossiblyCutEdge {
            data: self.data.as_ref(),
            flow: self.flow,
            index: self.index,
        }
    }

    pub fn duplicate_without_data(&self) -> Self {
        Self {
            data: None,
            flow: self.flow,
            index: self.index,
        }
    }

    pub fn reverse_mut(&mut self) {
        self.flow = self.flow.reverse();
    }

    pub fn reverse(self) -> Self {
        Self {
            data: self.data,
            flow: self.flow.reverse(),
            index: self.index,
        }
    }

    pub fn matches(&self, other: &Self) -> bool {
        self.flow == other.flow.reverse() && self.index == other.index
    }

    pub fn merge(self, other: Self) -> Option<Self> {
        if self.matches(&other) {
            Some(Self {
                data: Some(self.data.or(other.data)?),
                flow: self.flow,
                index: self.index,
            })
        } else {
            None
        }
    }

    pub fn flow(&self) -> Option<Flow> {
        self.flow.try_into().ok()
    }

    pub fn is_cut(&self) -> bool {
        !matches!(self.flow, Orientation::Undirected)
    }

    pub fn label(&self) -> String {
        let mut label = format!("{}", self.index);

        match self.flow {
            Orientation::Default => label.push_str("(left)"),
            Orientation::Reversed => label.push_str("(right)"),
            _ => {}
        }
        label
    }

    pub fn uncut(data: E, index: EdgeIndex) -> Self {
        Self {
            data: Some(data),
            flow: Orientation::Undirected,
            index,
        }
    }

    pub fn cut(&mut self, flow: Flow) {
        self.flow = flow.into();
    }
}

#[cfg(test)]
mod tests {
    use crate::{dot, parser::DotGraph};

    use super::*;

    #[test]
    fn test_roundtrip() {
        let gr: DotGraph = dot!(
            digraph {
              num = "1";
              overall_factor = "(AutG(1))^(-1)*InternalFermionLoopSign(-1)";
              0[int_id=V_71];
              1[int_id=V_71];
              2[int_id=V_98];
              3[int_id=V_98];

              ext0 [style=invis];
              2:0-> ext0 [id=0 is_cut=0 label="e-"];
              ext1 [style=invis];
              ext1-> 3:3 [id=6 is_cut=0 label="e+"];
              ext2 [style=invis];
              2:2-> ext2 [id=1 dir=back is_cut=2 label="e+"];
              ext3 [style=invis];
              ext3-> 3:1 [id=7 dir=back is_cut=2 label="e-"];
              0:4-> 1:5 [id=4 dir=back  label="d~"];
              0:6-> 1:7 [id=5  label="d"];
              0:8-> 3:9 [id=2 dir=none  label="a"];
              1:10-> 2:11 [id=3 dir=none  label="a"];
            }
        )
        .unwrap();

        let mut a = gr.graph.map(
            |_, _, v| v,
            |inv, _, _, ei, e| {
                e.map(|a| {
                    if let Some(is_cut) = a.get::<_, usize>("is_cut") {
                        let is_cut_h = Hedge(is_cut.unwrap());
                        let mut flow = inv.flow(is_cut_h);
                        let eid = inv[is_cut_h];
                        if eid != ei {
                            flow = -flow;
                        }
                        let mut e = PossiblyCutEdge::uncut(1, eid);
                        e.cut(flow);
                        e
                    } else {
                        PossiblyCutEdge::uncut(1, ei)
                    }
                })
            },
            |_, h| h,
        );

        // println!("{}", a.debug_cut_dot());
        a.glue_back_strict();

        fn from(e: &PossiblyCutEdge<i32>) -> DotEdgeData {
            let mut data = DotEdgeData::empty();

            data.add_statement(
                "cut_flow",
                match e.flow {
                    Orientation::Default => "aligned",
                    Orientation::Reversed => "reversed",
                    Orientation::Undirected => "uncut",
                },
            );

            let edge_id: usize = e.index.into();
            data.add_statement("edge_id", edge_id.to_string());
            data
        }

        let mut out = String::new();
        a.dot_serialize_fmt(&mut out, (), &|h| h.clone(), &from, &|v| v.clone())
            .unwrap();

        insta::assert_snapshot!(out,@r#"
        digraph {
          0[int_id="V_71"];
          1[int_id="V_71"];
          2[int_id="V_98"];
          3[int_id="V_98"];

          2:0	-> 3:3	 [id=0  cut_flow="aligned" edge_id="0"];
          2:2	-> 3:1	 [id=1 dir=back  cut_flow="aligned" edge_id="1"];
          0:8	-> 3:9	 [id=2 dir=none  cut_flow="uncut" edge_id="2"];
          1:10	-> 2:11	 [id=3 dir=none  cut_flow="uncut" edge_id="3"];
          0:4	-> 1:5	 [id=4 dir=back  cut_flow="uncut" edge_id="4"];
          0:6	-> 1:7	 [id=5  cut_flow="uncut" edge_id="5"];
        }
        "#);

        let olda = a.clone();

        a.round_trip_glue();

        assert_eq!(a, olda);

        let gr: DotGraph = dot!(
            digraph {
              num = "1";
              overall_factor = "(AutG(1))^(-1)";
              0[int_id=V_74];
              1[int_id=V_74];
              2[int_id=V_74];
              3[int_id=V_74];
              4[int_id=V_71];
              5[int_id=V_71];
              ext0 [style=invis];
              ext0-> 5:0 [id=9 dir=none is_cut=0 is_dummy=false particle="a"];
              0:1-> 1:2 [id=1 dir=back  is_dummy=false particle="d~"];
              0:3-> 1:4 [id=2  is_dummy=false particle="d"];
              0:5-> 3:6 [id=3 dir=none  is_dummy=false particle="g"];
              1:7-> 2:8 [id=4 dir=none  is_dummy=false particle="g"];
              2:9-> 4:10 [id=5 dir=back  is_dummy=false particle="d~"];
              2:11-> 5:12 [id=6  is_dummy=false particle="d"];
              3:13-> 4:14 [id=7  is_dummy=false particle="d"];
              3:15-> 5:16 [id=0 dir=back  is_dummy=false particle="d~"];
              ext9 [style=invis];
              4:17-> ext9 [id=8 dir=none is_cut=0 is_dummy=false particle="a"];
              }
        )
        .unwrap();

        let mut a = gr.graph.map(
            |_, _, v| v,
            |inv, _, _, ei, e| {
                e.map(|a| {
                    if let Some(is_cut) = a.get::<_, usize>("is_cut") {
                        let is_cut_val = is_cut.unwrap();
                        let is_cut_h = Hedge(is_cut_val);
                        let mut flow = inv.flow(is_cut_h);
                        let eid = inv[is_cut_h];
                        if eid != ei {
                            //a split edge with is_cut not set at the hedge id of it..
                            flow = -flow;
                        }
                        let mut e = PossiblyCutEdge::uncut(1, eid);
                        e.cut(flow);
                        e
                    } else {
                        PossiblyCutEdge::uncut(1, ei)
                    }
                })
            },
            |_, h| h,
        );

        // println!("{}", a.debug_cut_dot());
        a.glue_back_strict();

        let mut out = String::new();
        a.dot_serialize_fmt(&mut out, (), &|h| h.clone(), &from, &|v| v.clone())
            .unwrap();

        insta::assert_snapshot!(out,@r#"
        digraph {
          0[int_id="V_74"];
          1[int_id="V_74"];
          2[int_id="V_74"];
          3[int_id="V_74"];
          4[int_id="V_71"];
          5[int_id="V_71"];

          3:15	-> 5:16	 [id=0 dir=back  cut_flow="uncut" edge_id="0"];
          0:1	-> 1:2	 [id=1 dir=back  cut_flow="uncut" edge_id="1"];
          0:3	-> 1:4	 [id=2  cut_flow="uncut" edge_id="2"];
          0:5	-> 3:6	 [id=3 dir=none  cut_flow="uncut" edge_id="3"];
          1:7	-> 2:8	 [id=4 dir=none  cut_flow="uncut" edge_id="4"];
          2:9	-> 4:10	 [id=5 dir=back  cut_flow="uncut" edge_id="5"];
          2:11	-> 5:12	 [id=6  cut_flow="uncut" edge_id="6"];
          3:13	-> 4:14	 [id=7  cut_flow="uncut" edge_id="7"];
          4:17	-> 5:0	 [id=8 dir=none  cut_flow="aligned" edge_id="9"];
        }
        "#);

        let olda = a.clone();

        a.round_trip_glue();
        let mut out = String::new();
        a.dot_serialize_fmt(&mut out, (), &|h| h.clone(), &from, &|v| v.clone())
            .unwrap();

        let mut out2 = String::new();
        olda.dot_serialize_fmt(&mut out, (), &|h| h.clone(), &from, &|v| v.clone())
            .unwrap();

        assert_eq!(a, olda, "{}\nvs\n{}", out, out2);
    }
}
