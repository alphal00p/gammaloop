use itertools::Itertools;
use linnet::{
    half_edge::{
        involution::Orientation, nodestore::BitVecNeighborIter, EdgeAccessors, HedgeGraph,
    },
    parser::DotVertexData,
};
use spenso::{
    contraction::Contract,
    structure::{
        concrete_index::FlatIndex, representation::Euclidean, HasStructure, OrderedStructure,
        TensorStructure,
    },
    tensors::{data::GetTensorData, parametric::ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore},
    parse,
};

use crate::{
    model::{ArcParticle, ArcVertexRule, Model},
    numerator::aind::Aind,
    GammaLoopContext,
};
use color_eyre::Result;
use eyre::eyre;

use super::{
    edge::ParseEdge,
    hedge_data::ParseHedge,
    parse::{StripParse, ToQuoted},
};

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

impl From<&Vertex> for DotVertexData {
    fn from(value: &Vertex) -> Self {
        let mut v = DotVertexData::empty();
        v.add_statement("label", value.label.clone());
        v.add_statement("int_id", value.vertex_rule.name.as_str());
        v.add_statement("dod", value.dod);

        if value.num_color.size().unwrap() > 1 {
            v.add_statement(
                "num",
                value
                    .num_color
                    .contract(&value.num_spin)
                    .unwrap()
                    .scalar()
                    .unwrap()
                    .to_quoted(),
            );
        } else {
            v.add_statement(
                "num",
                value
                    .num_spin
                    .get_owned_linear(FlatIndex::from(0))
                    .unwrap()
                    .to_quoted(),
            );
            v.add_statement(
                "color_num",
                value
                    .num_color
                    .get_owned_linear(FlatIndex::from(0))
                    .unwrap()
                    .to_quoted(),
            );
        }

        v
    }
}
impl Vertex {}

#[derive(Debug, Clone)]
pub struct ParseVertex {
    pub label: Option<String>,
    pub vertex_rule: ArcVertexRule,
    pub spin_num: Option<Atom>,
    pub dod: Option<i32>,
    pub color_num: Option<Atom>,
}

impl ParseVertex {
    pub fn with_spin_num(mut self, spin_num: Atom) -> Self {
        self.spin_num = Some(spin_num);
        self
    }

    pub fn with_color_num(mut self, color_num: Atom) -> Self {
        self.color_num = Some(color_num);
        self
    }

    pub fn with_label(mut self, label: String) -> Self {
        self.label = Some(label);
        self
    }
}

impl From<ArcVertexRule> for ParseVertex {
    fn from(vertex_rule: ArcVertexRule) -> Self {
        ParseVertex {
            label: None,
            vertex_rule,
            dod: None,
            spin_num: None,
            color_num: None,
        }
    }
}

impl ParseVertex {
    pub fn parse<'a>(
        model: &'a Model,
        auto_detect_vertex_rule: bool,
    ) -> impl FnMut(
        &'a HedgeGraph<ParseEdge, &'a DotVertexData, ParseHedge>,
        BitVecNeighborIter<'a>,
        &'a &'a DotVertexData,
    ) -> Result<Self> {
        move |g, n, v| {
            let label = v.name().map(|id| id.to_string());

            let dod = v.get::<_, i32>("dod").transpose()?;

            let spin_num = v
                .get::<_, String>("num")
                .transpose()?
                .map(|a| a.strip_parse());
            let color_num = v
                .get::<_, String>("color_num")
                .transpose()?
                .map(|a| a.strip_parse());

            if let Some(n) = v.get::<_, String>("int_id") {
                let vertex_rule = model.get_vertex_rule(n.unwrap());

                Ok(ParseVertex {
                    dod,
                    label,
                    vertex_rule,
                    spin_num,
                    color_num,
                })
            } else if auto_detect_vertex_rule {
                let mut particles: Vec<ArcParticle> = n
                    .map(|h| {
                        let eid = g[&h];
                        let particle = match g.orientation(h).relative_to(g.flow(h)) {
                            Orientation::Reversed => g[eid].particle.get_anti_particle(model),
                            _ => g[eid].particle.clone(),
                        };
                        particle
                    })
                    .collect();
                particles.sort();

                let res = model.particle_set_to_vertex_rules_map.get(&particles);
                if let Some(res) = res {
                    if res.len() == 1 {
                        Ok(ParseVertex {
                            dod,
                            label,
                            vertex_rule: res[0].clone(),
                            spin_num,
                            color_num,
                        })
                    } else {
                        Err(eyre!("Multiple vertex rules for {:?}", particles))
                    }
                } else {
                    let particles = particles.iter().map(|p| p.name.as_str()).collect_vec();
                    Err(eyre!(
                        "Failed to find vertex rule for particles: {:?}",
                        particles
                    ))
                }
            } else {
                Err(eyre!("Vertex rule not supplied"))
            }
        }
    }
}
