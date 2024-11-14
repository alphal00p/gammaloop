use std::{num::TryFromIntError, ops::Add};

use ahash::AHashSet;
use bitvec::vec::BitVec;
use serde::{Deserialize, Serialize};

use crate::graph::half_edge::{HedgeGraph, PowersetIterator};

use super::{InternalSubGraph, SubGraphOps};

#[derive(Clone, Debug, Serialize, Deserialize, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Cycle {
    pub filter: BitVec,
    pub loop_count: Option<usize>,
}
impl Cycle {
    pub fn internal_graph<N, E>(self, graph: &HedgeGraph<N, E>) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_pessimist(self.filter, graph)
    }
    pub fn new_unchecked(filter: BitVec) -> Self {
        Self {
            filter,
            loop_count: None,
        }
    }
    pub fn is_circuit<E, V>(&self, graph: &HedgeGraph<E, V>) -> bool {
        for e in graph.iter_egde_node(&self.filter) {
            let adgacent = self.filter.intersection(&e.hairs);
            if adgacent.count_ones() != 2 {
                return false;
            }
        }
        if graph.count_connected_components(&self.filter) > 1 {
            return false;
        }
        true
    }
    pub fn new_circuit<E, V>(filter: BitVec, graph: &HedgeGraph<E, V>) -> Option<Self> {
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
    pub fn new<E, V>(filter: BitVec, graph: &HedgeGraph<E, V>) -> Option<Self> {
        for e in graph.iter_egde_node(&filter) {
            let adgacent = filter.intersection(&e.hairs);
            if adgacent.count_ones() % 2 == 1 {
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
        let mut pset = PowersetIterator::new(set.len().try_into()?);

        pset.next().unwrap(); //Skip the empty set

        for i in pset {
            let mut ones = i.iter_ones();

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

    fn add(self, other: &Cycle) -> Cycle {
        Cycle {
            filter: self.filter.clone() ^ other.filter.clone(),
            loop_count: None,
        }
    }
}
