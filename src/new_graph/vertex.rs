use itertools::Itertools;
use linnet::{
    half_edge::{
        involution::Orientation, nodestore::BitVecNeighborIter, EdgeAccessors, HedgeGraph,
        NodeIndex,
    },
    parser::DotVertexData,
};

use symbolica::atom::Atom;

use crate::{
    model::{ArcParticle, ArcVertexRule, Model},
    GammaLoopContext,
};
use color_eyre::Result;
use eyre::eyre;

use super::{
    edge::ParseEdge,
    global::ParseData,
    hedge_data::ParseHedge,
    parse::{StripParse, ToQuoted},
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
        v.add_statement("name", value.name.clone());
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

impl ParseVertex {
    pub(crate) fn with_num(mut self, num: Atom) -> Self {
        self.num = Some(num);
        self
    }

    pub(crate) fn with_label(mut self, label: String) -> Self {
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

            let dod = v.get::<_, i32>("dod").transpose()?;

            // let strict = v.get::<_, bool>("strict").transpose()?.unwrap_or(false);

            if let Some(num) = v.get::<_, String>("num") {
                let num = num?;
                // println!("Parsed with num:{}", num);
                Ok(ParseVertex {
                    dod,
                    // strict,
                    name,
                    vertex_rule: None,
                    num: Some(<String as StripParse<Atom>>::strip_parse(&num)),
                })
            } else if let Some(n) = v.get::<_, String>("int_id") {
                let vertex_rule = Some(model.get_vertex_rule(n.unwrap()));

                Ok(ParseVertex {
                    dod,
                    // strict,
                    name,
                    vertex_rule,
                    num: None,
                })
            } else {
                let mut node_id = NodeIndex(0);
                let mut particles: Vec<ArcParticle> = n
                    .filter_map(|h| {
                        node_id = g.node_id(h);
                        let eid = g[&h];
                        if g[eid].is_dummy {
                            return None;
                        }

                        // println!(
                        //     "{:?}{:?}{:?}{}",
                        //     g.flow(h),
                        //     g.orientation(h),
                        //     g.orientation(h).reverse().relative_to(g.flow(h)),
                        //     g[eid].particle.name
                        // );
                        let particle = match g.orientation(h).relative_to(g.flow(h)) {
                            Orientation::Reversed => {
                                g[eid].particle.particle()?.get_anti_particle(model)
                            }
                            _ => g[eid].particle.particle()?.clone(),
                        };
                        Some(particle)
                    })
                    .collect();
                particles.sort();

                let res = model.particle_set_to_vertex_rules_map.get(&particles);
                if let Some(res) = res {
                    if res.len() == 1 {
                        Ok(ParseVertex {
                            dod,
                            // strict,
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
