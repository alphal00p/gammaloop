use itertools::Itertools;
use linnet::{
    half_edge::{nodestore::BitVecNeighborIter, HedgeGraph, NodeIndex},
    parser::DotVertexData,
};

use symbolica::atom::Atom;

use crate::{
    feyngen::diagram_generator::NodeColorWithVertexRule,
    model::{ArcParticle, ArcVertexRule, Model},
    GammaLoopContext,
};
use color_eyre::Result;
use eyre::{eyre, Context};

use super::{
    edge::ParseEdge,
    global::ParseData,
    hedge_data::ParseHedge,
    parse::{extract_oriented_particles_from_vertex_hedges, StripParse, ToQuoted},
};

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Vertex {
    // #[bincode(with_serde)]
    pub name: String,
    pub vertex_rule: Option<ArcVertexRule>,
    pub num: Atom,
    // pub num_spin: ParamTensor<OrderedStructure<Euclidean, Aind>>,
    // pub num_color: ParamTensor<OrderedStructure<Euclidean, Aind>>,
    pub dod: i32,
}

impl Vertex {
    pub(crate) fn get_num(&self) -> Atom {
        self.num.clone()
    }
}

impl From<&Vertex> for DotVertexData {
    fn from(value: &Vertex) -> Self {
        let mut v = DotVertexData::empty();
        v.name = Some(value.name.clone());
        if let Some(vertex_rule) = &value.vertex_rule {
            v.add_statement("int_id", vertex_rule.name.as_str());
        }
        v.add_statement("dod", value.dod);
        v.add_statement("num", value.num.to_quoted());
        v
    }
}
impl Vertex {}

#[derive(Debug, Clone)]
pub struct ParseVertex {
    pub name: Option<String>,
    // pub strict: bool,
    pub vertex_rule: Option<ArcVertexRule>,
    pub num: Option<Atom>,
    pub dod: Option<i32>,
}

impl From<&ParseVertex> for DotVertexData {
    fn from(value: &ParseVertex) -> Self {
        let mut v = DotVertexData::empty();
        if let Some(name) = &value.name {
            v.name = Some(name.clone());
        }
        if let Some(vertex_rule) = &value.vertex_rule {
            v.add_statement("int_id", vertex_rule.name.as_str());
        }
        if let Some(dod) = value.dod {
            v.add_statement("dod", dod);
        }
        if let Some(num) = &value.num {
            v.add_statement("num", num.to_quoted());
        }
        v
    }
}

impl ParseVertex {
    pub fn with_num(mut self, num: Atom) -> Self {
        self.num = Some(num);
        self
    }

    pub fn with_label(mut self, label: String) -> Self {
        self.name = Some(label);
        self
    }
}

impl From<ArcVertexRule> for ParseVertex {
    fn from(vertex_rule: ArcVertexRule) -> Self {
        ParseVertex {
            name: None,
            // strict: false,
            vertex_rule: Some(vertex_rule),
            dod: None,
            num: None,
        }
    }
}

pub trait ParticleEdge {
    fn particle(&self) -> Option<ArcParticle>;
    fn is_dummy(&self) -> bool;
}

impl From<&NodeColorWithVertexRule> for ParseVertex {
    fn from(value: &NodeColorWithVertexRule) -> Self {
        value.vertex_rule.clone().into()
    }
}
impl ParseVertex {
    pub(crate) fn parse<'a>(
        model: &'a Model,
        parse_data: &'a ParseData,
    ) -> impl FnMut(
        &'a HedgeGraph<ParseEdge, &'a DotVertexData, ParseHedge>,
        BitVecNeighborIter<'a>,
        &'a &'a DotVertexData,
    ) -> Result<Self> {
        move |g, n, v| {
            let name = v.name().map(|id| id.to_string());

            let dod = v
                .get::<_, String>("dod")
                .transpose()
                .with_context(|| format!("Error parsing vertex dod"))?
                .map(|a| a.strip_parse())
                .transpose()?;

            if let Some(num) = v.get::<_, String>("num") {
                let num = num?;
                Ok(ParseVertex {
                    dod,
                    name,
                    vertex_rule: None,
                    num: Some(num.strip_parse().with_context(|| {
                        format!(
                            "Error parsing vertex num {num} of graph {}",
                            parse_data.name
                        )
                    })?),
                })
            } else if let Some(n) = v.get::<_, String>("int_id") {
                let vertex_rule = Some(model.get_vertex_rule(n.unwrap()));

                Ok(ParseVertex {
                    dod,
                    name,
                    vertex_rule,
                    num: None,
                })
            } else {
                let node_id = n
                    .clone()
                    .next()
                    .map(|h| g.node_id(h))
                    .unwrap_or(NodeIndex(0));
                let mut particles = extract_oriented_particles_from_vertex_hedges(g, n, model);
                particles.sort();

                let res = model.particle_set_to_vertex_rules_map.get(&particles);
                if let Some(res) = res {
                    if res.len() == 1 {
                        Ok(ParseVertex {
                            dod,
                            name,
                            vertex_rule: Some(res[0].clone()),
                            num: None,
                        })
                    } else {
                        Err(eyre!("Multiple vertex rules for {:?}", particles))
                    }
                } else {
                    let particles = particles.iter().map(|p| p.name.as_str()).collect_vec();
                    Err(eyre!(
                        "Failed to find vertex rule for particles: {:?} for node {node_id} in graph {}",
                        particles,parse_data.name,
                    ))
                }
            }
        }
    }
}
