use spenso::{
    structure::{representation::Euclidean, OrderedStructure},
    tensors::parametric::ParamTensor,
};

use crate::{model::ArcVertexRule, numerator::aind::Aind, GammaLoopContext};

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Vertex {
    // #[bincode(with_serde)]
    pub label: String,
    pub vertex_rule: ArcVertexRule,
    pub num_spin: ParamTensor<OrderedStructure<Euclidean, Aind>>,
    pub num_color: ParamTensor<OrderedStructure<Euclidean, Aind>>,
    pub dod: i32,
}

impl Vertex {}
