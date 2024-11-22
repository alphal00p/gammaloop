use std::ops::Add;

use bitvec::vec::BitVec;
use serde::{Deserialize, Serialize};

use crate::graph::half_edge::HedgeGraph;

use super::SubGraphOps;

#[derive(Clone, Debug, Serialize, Deserialize, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct Cycle {
    pub filter: BitVec,
}
impl Cycle {
    pub fn new_unchecked(filter: BitVec) -> Self {
        Self { filter }
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
        let circuit = Self { filter };
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

        Some(Self { filter })
    }
}

impl Add<&Cycle> for &Cycle {
    type Output = Cycle;

    fn add(self, other: &Cycle) -> Cycle {
        Cycle {
            filter: self.filter.clone() ^ other.filter.clone(),
        }
    }
}
