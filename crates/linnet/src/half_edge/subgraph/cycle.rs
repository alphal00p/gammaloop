use std::{num::TryFromIntError, ops::Add};

use ahash::AHashSet;

use crate::half_edge::{
    nodestore::NodeStorageOps,
    subgraph::{ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
    Hedge, HedgeGraph, PowersetIterator,
};

use super::{Inclusion, InternalSubGraph};

#[derive(Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents a cycle in a graph with a specific orientation or traversal direction.
///
/// While a [`Cycle`] represents the set of edges forming a cycle, a `SignedCycle`
/// implies a particular directed path along these edges. This can be important
/// in contexts where the direction of traversal matters.
pub struct SignedCycle {
    /// A bitmask indicating the set of half-edges that constitute this signed cycle,
    /// potentially implying a specific order or direction of traversal.
    pub filter: SuBitGraph,
    /// An optional count, possibly representing the number of times the cycle
    /// traverses its edges or some other property related to its structure.
    pub loop_count: Option<usize>,
}

impl SignedCycle {
    pub fn from_cycle<V, E, H, N: NodeStorageOps<NodeData = V>>(
        cycle: Cycle,
        according_to: Hedge,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Option<Self> {
        if !cycle.is_circuit(graph) {
            return None;
        }

        if !cycle.filter.includes(&according_to) {
            return None;
        }

        let mut filter = graph.empty_subgraph::<SuBitGraph>();

        let mut current_hedge = according_to;

        loop {
            if filter.includes(&current_hedge) {
                break;
            }
            filter.add(current_hedge);

            current_hedge = graph.inv(
                graph
                    .iter_crown(graph.node_id(current_hedge))
                    .find(|h| cycle.filter.includes(h) && (*h != current_hedge))?,
            );
        }

        Some(SignedCycle {
            filter,
            loop_count: cycle.loop_count,
        })
    }
}

#[derive(Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
/// Represents a cycle in a graph, defined as a set of half-edges.
///
/// A cycle is a path in a graph that starts and ends at the same vertex and
/// does not revisit intermediate vertices. This struct typically stores the
/// half-edges forming the cycle as a bitmask.
///
/// Operations on cycles, such as addition (symmetric difference), are common
/// in cycle basis theory.
pub struct Cycle {
    /// A bitmask representing the set of half-edges that form this cycle.
    /// The specific half-edges included define the path of the cycle.
    pub filter: SuBitGraph,
    /// An optional count, which might represent properties like the number of
    /// times the cycle traverses its edges (e.g., for fundamental cycles, this
    /// would often be 1).
    pub loop_count: Option<usize>,
}
impl Cycle {
    pub fn internal_graph<E, V, H, N: NodeStorageOps<NodeData = V>>(
        self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_pessimist(self.filter, graph)
    }
    pub fn new_unchecked(filter: SuBitGraph) -> Self {
        Self {
            filter,
            loop_count: None,
        }
    }
    pub fn is_circuit<E, V, H, N: NodeStorageOps<NodeData = V>>(
        &self,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> bool {
        for (_, c, _) in graph.iter_nodes_of(&self.filter) {
            let adgacent = c.filter(|a| self.filter.includes(a));
            let count = adgacent.count();
            if count != 2 {
                return false;
            }
        }
        if graph.count_connected_components(&self.filter) > 1 {
            return false;
        }
        true
    }
    pub fn new_circuit<E, V, H, N: NodeStorageOps<NodeData = V>>(
        filter: SuBitGraph,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Option<Self> {
        let circuit = Self {
            filter,
            loop_count: Some(1),
        };
        if circuit.is_circuit(graph) {
            Some(circuit)
        } else {
            None
        }
    }
    pub fn new<E, V, H, N: NodeStorageOps<NodeData = V>>(
        filter: SuBitGraph,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> Option<Self> {
        for (_, c, _) in graph.iter_nodes_of(&filter) {
            let adgacent = c.filter(|a| filter.includes(a));
            if adgacent.count() % 2 == 1 {
                return None;
            }
        }

        Some(Self {
            filter,
            loop_count: None,
        })
    }

    pub fn all_sum_powerset_filter_map(
        set: &[Self],
        filter_map: &impl Fn(Self) -> Option<Self>,
    ) -> Result<AHashSet<Self>, TryFromIntError> {
        let mut s = AHashSet::new();
        let mut pset = PowersetIterator::<usize>::new(set.len().try_into()?);

        pset.next().unwrap(); //Skip the empty set

        for i in pset {
            let mut ones = i.included_iter();

            let mut union = set[ones.next().unwrap()].clone();

            for o in ones {
                union = &union + &set[o];
            }

            if let Some(union) = filter_map(union) {
                s.insert(union);
            }
        }

        Ok(s)
    }
}

impl Add<&Cycle> for &Cycle {
    type Output = Cycle;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, other: &Cycle) -> Cycle {
        Cycle {
            filter: self.filter.sym_diff(&other.filter),
            loop_count: None,
        }
    }
}
