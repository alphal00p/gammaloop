use linnet::{
    half_edge::{
        involution::{EdgeData, EdgeIndex, Flow, HedgePair},
        HedgeGraph,
    },
    parser::{DotEdgeData, DotHedgeData, DotVertexData},
};
use spenso::network::library::TensorLibraryData;
use symbolica::atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder};

use crate::{
    graph::BareEdge,
    model::{ArcParticle, ArcPropagator, Model},
    momentum_sample::LoopIndex,
    numerator::aind::NewAind,
    utils::GS,
};

use color_eyre::Result;
use eyre::eyre;

use super::{
    parse::{StripParse, ToQuoted},
    NumHedgeData,
};

#[derive(Debug, Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Edge {
    // #[bincode(with_serde)]
    pub name: String,
    // pub edge_type: EdgeType,
    pub propagator: ArcPropagator,
    pub particle: ArcParticle,
    pub num: Atom,
    // pub spin_num: Atom,
    pub dod: i32,
    pub is_dummy: bool, // #[bincode(with_serde)]
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
            // color_num: Atom::one(),
            num: Atom::one(),
            dod: -2,
            is_dummy: false,
        }
    }
}

impl From<&Edge> for DotEdgeData {
    fn from(value: &Edge) -> Self {
        let mut e = DotEdgeData::empty();
        e.add_statement("name", value.name.clone());
        e.add_statement("particle", format!("\"{}\"", value.particle.name));
        e.add_statement("dod", value.dod);
        e.add_statement("num", value.num.to_quoted());
        // e.add_statement("color_num", value.color_num.to_quoted());
        e.add_statement("is_dummy", value.is_dummy);
        e
    }
}

#[derive(Debug, Clone)]
pub struct ParseEdge {
    pub label: Option<String>,
    pub particle: ArcParticle,
    pub dod: Option<i32>,
    pub is_dummy: bool,
    pub lmb_id: Option<LoopIndex>,
    pub num: Option<Atom>,
    // pub color_num: Option<sAtom>,
}

impl ParseEdge {
    pub fn new(particle: ArcParticle) -> Self {
        ParseEdge {
            is_dummy: false,
            label: None,
            particle,
            dod: None,
            lmb_id: None,
            num: None,
        }
    }

    pub fn is_dummy(mut self) -> Self {
        self.is_dummy = true;
        self
    }

    pub fn with_label(mut self, label: String) -> Self {
        self.label = Some(label);
        self
    }

    pub fn with_dod(mut self, dod: i32) -> Self {
        self.dod = Some(dod);
        self
    }

    pub fn with_lmb_id(mut self, lmb_id: LoopIndex) -> Self {
        self.lmb_id = Some(lmb_id);
        self
    }

    pub fn with_num(mut self, num: Atom) -> Self {
        self.num = Some(num);
        self
    }
}

impl ParseEdge {
    pub fn localize_ainds(atom: impl AtomCore, eid: EdgeIndex, hedge_pair: HedgePair) -> Atom {
        let a = atom
            .replace(GS.edgeid)
            .with(Atom::num(usize::from(eid) as i64))
            .replace_map(|term, ctx, out| {
                if let AtomView::Fun(f) = term {
                    if f.get_symbol() == GS.edgeid {
                        if f.get_nargs() == 1 {
                            if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                if let Ok(u) = u16::try_from(i) {
                                    *out = eid.aind(u).into();
                                    return true;
                                }
                            }
                        }
                    }
                }
                false
            });

        match hedge_pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => a
                .replace(GS.sink_id)
                .with(Atom::num(sink.0 as i64))
                .replace_map(|term, ctx, out| {
                    if let AtomView::Fun(f) = term {
                        if f.get_symbol() == GS.sink_id {
                            if f.get_nargs() == 1 {
                                if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                    if let Ok(u) = u16::try_from(i) {
                                        *out = sink.aind(u).into();
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                    false
                })
                .replace(GS.source_id)
                .with(Atom::num(source.0 as i64))
                .replace_map(|term, ctx, out| {
                    if let AtomView::Fun(f) = term {
                        if f.get_symbol() == GS.source_id {
                            if f.get_nargs() == 1 {
                                if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                    if let Ok(u) = u16::try_from(i) {
                                        *out = source.aind(u).into();
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                    false
                }),
            HedgePair::Unpaired { hedge, flow } => match flow {
                Flow::Source => a
                    .replace(GS.source_id)
                    .with(Atom::num(hedge.0 as i64))
                    .replace_map(|term, ctx, out| {
                        if let AtomView::Fun(f) = term {
                            if f.get_symbol() == GS.source_id {
                                if f.get_nargs() == 1 {
                                    if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                        if let Ok(u) = u16::try_from(i) {
                                            *out = hedge.aind(u).into();
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                        false
                    }),
                Flow::Sink => a
                    .replace(GS.sink_id)
                    .with(Atom::num(hedge.0 as i64))
                    .replace_map(|term, ctx, out| {
                        if let AtomView::Fun(f) = term {
                            if f.get_symbol() == GS.sink_id {
                                if f.get_nargs() == 1 {
                                    if let Ok(i) = i64::try_from(f.iter().next().unwrap()) {
                                        if let Ok(u) = u16::try_from(i) {
                                            *out = hedge.aind(u).into();
                                            return true;
                                        }
                                    }
                                }
                            }
                        }
                        false
                    }),
            },
        }
    }
    pub fn parse<'a>(
        model: &'a Model,
    ) -> impl FnMut(
        &'a HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData>,
        EdgeIndex,
        HedgePair,
        EdgeData<&'a DotEdgeData>,
    ) -> Result<EdgeData<Self>> {
        |graph: &'a HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData>,
         eid: EdgeIndex,
         p: HedgePair,
         e_data: EdgeData<&'a DotEdgeData>| {
            e_data.map_result(|e| {
                let label = e.get::<_, String>("label").transpose()?;

                let lmb_id: Option<LoopIndex> = e
                    .get::<_, usize>("lmb_id")
                    .transpose()?
                    .map(LoopIndex::from);

                let dod = e.get::<_, i32>("dod").transpose()?;
                let is_dummy = e.get::<_, bool>("is_dummy").transpose()?.unwrap_or(false);

                let num = e.get::<_, String>("num").transpose()?.map(|a| {
                    Self::localize_ainds(<String as StripParse<Atom>>::strip_parse(&a), eid, p)
                });

                let particle = if let Some(v) = e.get::<_, isize>("pdg") {
                    model.get_particle_from_pdg(v?)
                } else if let Some(v) = e.get::<_, String>("particle") {
                    let pname = v?;
                    let pname = pname
                        .as_str()
                        .strip_prefix('"')
                        .unwrap_or(&pname)
                        .strip_suffix('"')
                        .unwrap_or(&pname);
                    model.get_particle(pname)
                } else {
                    return Err(eyre!("no pdg or name found for edge"));
                };

                Ok(ParseEdge {
                    is_dummy,
                    dod,
                    lmb_id,
                    particle,
                    num,
                    label,
                })
            })
        }
    }
}
