use std::{hash::Hash, num::TryFromIntError};

use ahash::AHashSet;
use bitvec::order::Lsb0;
use indenter::CodeFormatter;
use std::fmt::Write;

use super::{
    involution::HedgePair, nodestore::NodeStorageOps, GVEdgeAttrs, Hedge, HedgeGraph,
    PowersetIterator,
};

const BASE62_ALPHABET: &[u8; 62] =
    b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

/// Defines operations that can be performed between subgraphs or collections of subgraphs.
///
/// This trait provides utilities for common set-like operations such as unions,
/// intersections (implicitly via other methods like `subtract`), and symmetric differences,
/// often applied pairwise to collections of subgraphs or iteratively.
///
/// # Type Parameters
/// - `Other`: The type of the other subgraph to perform operations with. Defaults to `Self`.
pub trait SubSetOps<ID = Hedge, Other: SubSetLike<ID> = Self>: SubSetLike<ID> {
    // type Other = Self; // Potentially for future use if more varied inter-type operations are needed.

    /// Applies a pairwise operation `op` between all elements of `left` (a mutable set)
    /// and all elements of `right` (a slice), adding results to `left`.
    /// Returns `true` if `left` was modified.
    fn all_pairwise_ops(
        left: &mut AHashSet<Self>,
        right: &[Self],
        op: &impl Fn(&Self, &Self) -> Self,
    ) -> bool {
        Self::all_pairwise_ops_filter_map(left, right, op, &|x| Some(x))
    }

    fn all_pairwise_ops_filter_map(
        left: &mut AHashSet<Self>,
        right: &[Self],
        op: &impl Fn(&Self, &Self) -> Self,
        filter_map: &impl Fn(Self) -> Option<Self>,
    ) -> bool {
        let mut added = false;
        let mut new = AHashSet::new();
        for l in left.iter() {
            for r in right {
                if let Some(n) = filter_map(op(l, r)) {
                    new.insert(n);
                }
            }
        }
        for n in new.drain() {
            if left.insert(n) {
                added = true;
            }
        }
        added
    }

    fn all_pairwise_unions(left: &mut AHashSet<Self>, right: &[Self]) -> bool {
        Self::all_pairwise_ops(left, right, &|l, r| l.union(r))
    }

    fn all_pairwise_sym_diff(left: &mut AHashSet<Self>, right: &[Self]) -> bool {
        Self::all_pairwise_ops(left, right, &|l, r| l.sym_diff(r))
    }

    fn all_ops_iterative_filter_map(
        set: &[Self],
        op: &impl Fn(&Self, &Self) -> Self,
        filter_map: &impl Fn(Self) -> Option<Self>,
    ) -> AHashSet<Self> {
        let mut s: AHashSet<_> = set.iter().cloned().collect();
        while Self::all_pairwise_ops_filter_map(&mut s, set, op, filter_map) {}
        s.drain().filter_map(filter_map).collect()
    }

    fn all_unions_iterative(set: &[Self]) -> AHashSet<Self> {
        Self::all_ops_iterative_filter_map(set, &|a, b| a.union(b), &|a| Some(a))
    }

    fn all_sym_diff_iterative(set: &[Self]) -> AHashSet<Self> {
        Self::all_ops_iterative_filter_map(set, &|a, b| a.sym_diff(b), &|a| Some(a))
    }

    fn all_op_powerset_filter_map(
        set: &[Self],
        op: impl Fn(&mut Self, &Self),
        filter_map: &impl Fn(Self) -> Option<Self>,
    ) -> Result<AHashSet<Self>, TryFromIntError> {
        let mut s = AHashSet::new();
        let mut pset = PowersetIterator::<usize>::new(set.len().try_into()?);

        pset.next().unwrap(); //Skip the empty set

        for i in pset {
            let mut ones = i.included_iter();

            let mut union = set[ones.next().unwrap()].clone();

            for o in ones {
                op(&mut union, &set[o]);
            }

            if let Some(union) = filter_map(union) {
                s.insert(union);
            }
        }

        Ok(s)
    }

    fn all_unions_powerset_filter_map(
        set: &[Self],
        filter_map: &impl Fn(Self) -> Option<Self>,
    ) -> Result<AHashSet<Self>, TryFromIntError> {
        Self::all_op_powerset_filter_map(set, |l, r| l.union_with(r), filter_map)
    }

    fn all_sym_diff_powerset(
        set: &[Self],
        filter_map: &impl Fn(Self) -> Option<Self>,
    ) -> Result<AHashSet<Self>, TryFromIntError> {
        Self::all_op_powerset_filter_map(set, |l, r| l.sym_diff_with(r), filter_map)
    }

    fn n_op(n: usize, set: &[Self], op: &impl Fn(&Self, &Self) -> Self) -> AHashSet<Self> {
        if n == 0 {
            AHashSet::new()
        } else {
            let mut s = Self::n_op(n - 1, set, op);
            Self::all_pairwise_ops(&mut s, set, op);
            s
        }
    }

    fn intersect_with(&mut self, other: &Self);
    fn union_with(&mut self, other: &Self);
    fn union_with_iter(&mut self, other: impl Iterator<Item = ID>);
    fn sym_diff_with(&mut self, other: &Self);
    fn empty_intersection(&self, other: &Self) -> bool;
    fn empty_union(&self, other: &Self) -> bool;
    fn intersection(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.intersect_with(other);
        new
    }
    fn union(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.union_with(other);

        new
    }
    fn sym_diff(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.sym_diff_with(other);
        new
    }

    fn subtract_with(&mut self, other: &Self);
    fn subtract(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.subtract_with(other);
        new
    }
}

pub trait SubGraphOps {
    fn complement<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Self;
}

/// A trait for checking inclusion or intersection relationships between a subgraph
/// and another entity `T`.
///
/// This is used to determine if a subgraph contains certain elements (like specific hedges,
/// ranges of hedges, or other subgraphs) or overlaps with them.
///
/// # Type Parameters
/// - `T`: The type of the entity to check for inclusion or intersection with.
pub trait Inclusion<T> {
    /// Checks if this subgraph fully includes/contains the `other` entity.
    fn includes(&self, other: &T) -> bool;
    /// Checks if this subgraph has any common elements with the `other` entity.
    fn intersects(&self, other: &T) -> bool;
}

/// An iterator over [`Hedge`]s included in a subgraph.
///
/// This struct typically wraps an iterator over set bits in a `BitVec`
/// (like `bitvec::slice::IterOnes`) and maps the `usize` indices to `Hedge` types.
#[derive(Clone, Debug)]
pub struct SubSetIter<'a, ID = Hedge> {
    /// The underlying iterator that yields `usize` indices, which are then mapped to `Hedge`.
    iter: std::iter::Map<bitvec::slice::IterOnes<'a, usize, Lsb0>, fn(usize) -> ID>,
}

impl<ID> DoubleEndedIterator for SubSetIter<'_, ID> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.iter.next_back()
    }
}

impl<ID> ExactSizeIterator for SubSetIter<'_, ID> {
    fn len(&self) -> usize {
        self.iter.len()
    }
}

impl<ID> Iterator for SubSetIter<'_, ID> {
    type Item = ID;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }
}

pub trait BaseSubgraph:
    SubGraphLike + ModifySubSet<HedgePair> + ModifySubSet<Hedge> + SubSetOps<Hedge>
{
    fn from_filter<E, V, H, N: NodeStorageOps<NodeData = V>, F: FnMut(&E) -> bool>(
        graph: &HedgeGraph<E, V, H, N>,
        filter: F,
    ) -> Self;

    /// Creates a base subgraph containing all hedges yielded by the iterator.
    /// `len` specifies the total number of possible hedges in the graph context.
    fn from_hedge_iter<I: Iterator<Item = Hedge>>(iter: I, len: usize) -> Self;
}

/// Defines methods for modifying a subgraph by adding or removing elements.
///
/// # Type Parameters
/// - `Index`: The type of element to add or remove (e.g., [`Hedge`], [`HedgePair`]).
pub trait ModifySubSet<Index> {
    /// Adds an element `index` to the subgraph.
    fn add(&mut self, index: Index);

    /// Removes an element `index` from the subgraph.
    fn sub(&mut self, index: Index);
}

pub trait SubGraphLike: SubSetLike<Hedge> + Inclusion<HedgePair> {
    /// maximal graph that contains all nodes of the subgraph
    fn covers<E, V, H, N: NodeStorageOps<NodeData = V>, S: SubSetOps<Hedge> + SubSetLike<Hedge>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> S {
        let mut covering = graph.empty_subgraph::<S>();
        for i in self.included_iter() {
            covering.union_with_iter(graph.neighbors(i))
        }
        covering
    }

    fn background_color(&self, hedge_pair: Option<HedgePair>) -> Option<String> {
        let color = "gray".to_string();

        if let Some(p) = hedge_pair {
            if let HedgePair::Split { .. } = p {
                Some(color)
            } else {
                None
            }
        } else {
            Some(color)
        }
    }

    /// Number of full edges included in the subgraph
    fn nedges<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> usize; //not counting unpaired hedges

    fn dot_fmt<W: std::fmt::Write, E, V, H, N: NodeStorageOps<NodeData = V>, Str: AsRef<str>>(
        &self,
        writer: &mut W,
        graph: &HedgeGraph<E, V, H, N>,
        graph_info: Str,
        hedge_attr: &impl Fn(&H) -> Option<String>,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> Result<(), std::fmt::Error> {
        writeln!(writer, "digraph {{")?;
        let mut writer = CodeFormatter::new(writer, "  ");
        writer.indent(1);
        write!(
            writer,
            "\nnode\t [shape=circle,height=0.1,label=\"\"];\noverlap = \"scale\";\nlayout = \"neato\";",
        )?;
        write!(writer, "\n{}", graph_info.as_ref())?;

        for (n, _, v) in graph.iter_nodes_of(self) {
            if let Some(a) = node_attr(v) {
                write!(writer, "\n{}\t [{}];", n.0, a)?;
            }
        }

        for (hedge_pair, eid, data) in graph.iter_edges() {
            let subgraph_pair = hedge_pair.with_subgraph(self);

            let attr = GVEdgeAttrs {
                color: self.background_color(subgraph_pair),
                label: None,
                other: edge_attr(data.data),
            };
            if let Some(p) = subgraph_pair {
                let attr = p.fill_color(attr);

                p.add_data_of(graph, self).dot_fmt(
                    &mut writer,
                    graph,
                    eid,
                    hedge_attr,
                    |a| a.to_string(),
                    data.orientation,
                    attr,
                )?;
            } else {
                let attr = hedge_pair.fill_color(attr);

                hedge_pair.add_data_of(graph, self).dot_fmt(
                    &mut writer,
                    graph,
                    eid,
                    hedge_attr,
                    |a| a.to_string(),
                    data.orientation,
                    attr,
                )?;
            }
        }
        writeln!(writer, "}}")?;
        Ok(())
    }

    fn dot_io<W: std::io::Write, E, V, H, N: NodeStorageOps<NodeData = V>, Str: AsRef<str>>(
        &self,
        writer: &mut W,
        graph: &HedgeGraph<E, V, H, N>,
        graph_info: Str,
        hedge_attr: &impl Fn(&H) -> Option<String>,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> Result<(), std::io::Error> {
        writeln!(writer, "digraph {{")?;
        writeln!(
            writer,
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";",
        )?;
        writeln!(writer, "{}", graph_info.as_ref())?;

        for (n, _, v) in graph.iter_nodes_of(self) {
            if let Some(a) = node_attr(v) {
                writeln!(writer, "  {} [{}];", n.0, a)?;
            }
        }

        for (hedge_pair, eid, data) in graph.iter_edges() {
            let subgraph_pair = hedge_pair.with_subgraph(self);

            let attr = GVEdgeAttrs {
                color: self.background_color(subgraph_pair),
                label: None,
                other: edge_attr(data.data),
            };
            write!(writer, "  ")?;
            if let Some(p) = subgraph_pair {
                let attr = p.fill_color(attr);

                p.add_data_of(graph, self).dot_io(
                    writer,
                    graph,
                    eid,
                    hedge_attr,
                    |a| a.to_string(),
                    data.orientation,
                    attr,
                )?;
            } else {
                let attr = hedge_pair.fill_color(attr);
                hedge_pair.add_data_of(graph, self).dot_io(
                    writer,
                    graph,
                    eid,
                    hedge_attr,
                    |a| a.to_string(),
                    data.orientation,
                    attr,
                )?;
            }
        }
        writeln!(writer, "}}")?;
        Ok(())
    }

    fn hairs(&self, node: impl Iterator<Item = Hedge>) -> SuBitGraph {
        let mut hairs = SuBitGraph::empty(self.size());
        for h in node {
            if self.includes(&h) {
                hairs.add(h)
            }
        }
        hairs
    }
}

/// The core trait defining the properties and behaviors of a subgraph representation.
///
/// A subgraph is fundamentally a selection of half-edges from a larger [`HedgeGraph`].
/// This trait provides methods for querying the subgraph's contents, its relationship
/// to the parent graph, and performing common operations.
///
/// Implementors must also satisfy `Clone`, `Eq`, `Hash`, and various [`Inclusion`] bounds.
pub trait SubSetLike<ID = Hedge>:
    Clone
    + Eq
    + Hash
    + Inclusion<Self>
    + Inclusion<Self::Base>
    + Inclusion<std::ops::Range<ID>>
    + Inclusion<std::ops::RangeToInclusive<ID>>
    + Inclusion<std::ops::RangeInclusive<ID>>
    + Inclusion<std::ops::RangeTo<ID>>
    + Inclusion<std::ops::RangeFrom<ID>>
    + Inclusion<ID>
{
    type Base: SubSetLike<ID>;
    type BaseIter<'a>: Iterator<Item = ID> + DoubleEndedIterator
    where
        Self: 'a,
        ID: 'a;

    /// Contains an id with index >= hedge
    fn has_greater(&self, hedge: ID) -> bool {
        self.intersects(&(hedge..))
    }

    /// Contains an id with index < hedge
    fn has_lesser(&self, hedge: ID) -> bool {
        self.intersects(&(..hedge))
    }

    ///joins two subsets into one, (this is not union)
    fn join(mut self, other: Self) -> Self {
        self.join_mut(other);
        self
    }

    /// Appends all incuded half-edges at the end of other
    fn join_mut(&mut self, other: Self);

    fn string_label(&self) -> String;
    fn from_base62(label: &str, size: usize) -> Option<Self>;
    fn included_iter(&self) -> Self::BaseIter<'_>;
    // fn included_iter_(&self) -> Self::BaseIter<'_>;
    // SubGraphHedgeIter {
    //     SubGraphHedgeIter {
    //         iter: self.included().iter_ones().map(Hedge),
    //     }
    // }
    /// Returns a simple Self::Base of all included hedges
    fn included(&self) -> &Self::Base;

    /// Number of half-edges in the graph this is a subgraph of
    fn size(&self) -> usize;
    /// Number of half-edges included in the subgraph
    fn n_included(&self) -> usize;

    fn empty(size: usize) -> Self;
    fn is_empty(&self) -> bool;
}

pub mod contracted;
pub use contracted::ContractedSubGraph;
pub mod cut;
pub use cut::OrientedCut;
pub mod cycle;
pub use cycle::Cycle;
pub mod internal;
pub use internal::InternalSubGraph;
pub mod node;
pub use node::HedgeNode;
pub mod empty;
pub use empty::Empty;
pub mod full;
pub use full::FullOrEmpty;
pub mod subset;
pub use subset::SuBitGraph;
