use spenso::network::library::TensorLibraryData;
use symbolica::atom::Atom;

use crate::{
    graph::BareEdge,
    model::{ArcParticle, ArcPropagator},
};

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Edge {
    // #[bincode(with_serde)]
    pub name: String,
    // pub edge_type: EdgeType,
    pub propagator: ArcPropagator,
    pub particle: ArcParticle,
    pub color_num: Atom,
    pub spin_num: Atom,
    pub dod: i32,
    // #[bincode(with_serde)]
    // pub internal_index: Vec<AbstractIndex>,
}

impl Edge {
    pub fn n_dummies(&self) -> usize {
        5
    }
}

impl From<BareEdge> for Edge {
    fn from(value: BareEdge) -> Self {
        Self {
            name: value.name.into(),
            propagator: ArcPropagator(value.propagator),
            particle: value.particle,
            color_num: Atom::one(),
            spin_num: Atom::one(),
            dod: -2,
        }
    }
}
