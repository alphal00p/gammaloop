use std::hash::{Hash, Hasher};

use linnet::half_edge::{
    HedgeGraph,
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike},
};
use tracing::{debug, error};

use crate::{
    graph::{LMBext, LoopMomentumBasis, cuts::CutSet},
    uv::{ApproximationType, UltravioletGraph},
};

#[derive(Clone, Debug, Eq)]
pub struct Spinney {
    pub subgraph: InternalSubGraph,
    components: Vec<SuBitGraph>,
    pub dod: i32,
    pub lmb: LoopMomentumBasis,
    pub renormalization_scheme: ApproximationType,
    max_comp_loop_count: usize,
}

impl Spinney {
    pub fn compatible_with(&self, cut: &CutSet) -> bool {
        self.renormalization_scheme == ApproximationType::VacuumLimit
            || !self.subgraph.filter.intersects(&cut.union)
    }

    pub fn empty<E, V, H, G: AsRef<HedgeGraph<E, V, H>> + LMBext + ?Sized>(g: &G) -> Self {
        Self {
            subgraph: InternalSubGraph::empty(g.as_ref().n_hedges()),
            components: vec![],
            dod: 0,
            lmb: g.empty_lmb(),
            renormalization_scheme: ApproximationType::MUV,
            max_comp_loop_count: 0,
        }
    }

    pub fn new<E, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>> + ?Sized>(
        subgraph: InternalSubGraph,
        g: &G,
        lmb: &LoopMomentumBasis,
    ) -> Self {
        Self::with_scheme(subgraph, g, lmb, ApproximationType::MUV)
    }

    pub fn with_scheme<E, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>> + ?Sized>(
        subgraph: InternalSubGraph,
        g: &G,
        lmb: &LoopMomentumBasis,
        renormalization_scheme: ApproximationType,
    ) -> Self {
        let components = g.as_ref().connected_components(&subgraph);
        let max_comp_loop_count = components
            .iter()
            .map(|component| g.as_ref().cyclotomatic_number(component))
            .max()
            .unwrap_or(0);
        let lmb = g.compatible_sub_lmb(&subgraph, g.dummy_less_full_crown(&subgraph), lmb);

        let dod = g.dod(&subgraph);

        if dod != g.local_dod(&subgraph) {
            error!(
                "dod mismatch: global={} local={} of graph:\n{}",
                dod,
                g.local_dod(&subgraph),
                g.dod(&subgraph),
            );
        }

        debug!(
            dod = %dod,
            graph = %g.dot_lmb_of(&subgraph,&lmb),
            renormalization_scheme = %renormalization_scheme,
            string_label =%subgraph.string_label(),
            "Constucted spinney with scheme: {:?}",
            renormalization_scheme
        );

        Self {
            components,
            dod,
            lmb,
            subgraph,
            renormalization_scheme,
            max_comp_loop_count,
        }
    }

    pub fn filter(&self) -> &SuBitGraph {
        &self.subgraph.filter
    }

    pub fn max_comp_loop_count(&self) -> usize {
        self.max_comp_loop_count
    }

    pub fn n_components(&self) -> usize {
        self.components.len()
    }
}

impl Hash for Spinney {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.subgraph.hash(state);
    }
}

impl PartialEq for Spinney {
    fn eq(&self, other: &Self) -> bool {
        self.subgraph == other.subgraph
    }
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl PartialOrd for Spinney {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.subgraph == other.subgraph {
            Some(std::cmp::Ordering::Equal)
        } else if self.subgraph.includes(&other.subgraph) {
            Some(std::cmp::Ordering::Greater)
        } else if other.subgraph.includes(&self.subgraph) {
            Some(std::cmp::Ordering::Less)
        } else {
            None
        }
    }
}

impl Ord for Spinney {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.filter().cmp(other.filter())
    }
}
