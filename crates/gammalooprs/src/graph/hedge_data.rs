use itertools::Itertools;
use linnet::{
    half_edge::involution::{Flow, Hedge},
    parser::DotHedgeData,
};
use serde::{Deserialize, Serialize};
use spenso::structure::{
    OrderedStructure, PermutedStructure, TensorStructure,
    representation::{LibraryRep, LibrarySlot},
    slot::IsAbstractSlot,
};
use symbolica::atom::{Atom, FunctionBuilder};

use crate::{graph::edge::PossibleParticle, numerator::aind::Aind};

use super::{Autogen, parse::ParseGraph};

use color_eyre::Result;

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct HedgeIndices {
    pub vertex_indices: OrderedStructure<LibraryRep, Aind>,
    pub edge_indices: OrderedStructure<LibraryRep, Aind>,
}

impl HedgeIndices {
    pub(crate) fn new(edge_indices: OrderedStructure<LibraryRep, Aind>) -> Self {
        Self {
            vertex_indices: edge_indices.clone().dual(),
            edge_indices,
        }
    }
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct NumIndices {
    pub color_indices: HedgeIndices,
    pub spin_indices: HedgeIndices,
}

impl NumIndices {
    fn from_particle(particle: &PossibleParticle, flow: Flow, h: Hedge) -> Self {
        let creps = particle.color_reps(flow);
        let sreps = particle.spin_reps();

        let mut init = 0;
        let mut last = None;
        let color_structure: PermutedStructure<_> = creps
            .external_reps_iter()
            .map(|r| {
                if let Some(l) = last {
                    if l != r {
                        last = Some(r);
                        init = 0;
                    } else {
                        init += 1;
                    }
                } else {
                    last = Some(r);
                    init = 0;
                }
                r.slot(Aind::Hedge(h.0 as u16, init))
            })
            .collect();

        let spin_structure: PermutedStructure<_> = sreps
            .external_reps_iter()
            .map(|r| {
                if let Some(l) = last {
                    if l != r {
                        last = Some(r);
                        init = 0;
                    } else {
                        init += 1;
                    }
                } else {
                    last = Some(r);
                    init = 0;
                }
                r.slot(Aind::Hedge(h.0 as u16, init))
            })
            .collect();

        NumIndices {
            color_indices: HedgeIndices::new(color_structure.structure),
            spin_indices: HedgeIndices::new(spin_structure.structure),
        }
    }

    pub(crate) fn parse<'a>(
        graph: &'a ParseGraph,
    ) -> impl FnMut(Hedge, &'a ParseHedgeData) -> Self {
        |h, _| {
            let eid = graph[&h];
            let flow = graph.flow(h);

            Self::from_particle(&graph[eid].particle, flow, h)
        }
    }
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct HedgeData {
    pub num_indices: NumIndices,
    pub ufo_order: Autogen<u8>,
}

impl HedgeData {
    pub(crate) fn edge_color_slots<'a>(&'a self) -> impl Iterator<Item = LibrarySlot<Aind>> + 'a {
        self.num_indices
            .color_indices
            .edge_indices
            .external_structure_iter()
    }

    pub(crate) fn edge_spin_slots<'a>(&'a self) -> impl Iterator<Item = LibrarySlot<Aind>> + 'a {
        self.num_indices
            .spin_indices
            .edge_indices
            .external_structure_iter()
    }

    pub(crate) fn polarization(&self, mut builder: FunctionBuilder) -> Atom {
        for s in self.edge_spin_slots() {
            builder = builder.add_arg(s.to_atom())
        }
        builder.finish()
    }

    pub(crate) fn color_kronekers(&self, other: &Self) -> Atom {
        let mut color = Atom::num(1);
        for (i, j) in self.edge_color_slots().zip_eq(other.edge_color_slots()) {
            if i.rep().matches(&j.rep()) {
                color *= j.rep().id(i.aind, j.aind);
            } else {
                panic!("Should be the same rep found:{} and {}", i, j)
            }
        }

        color
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ParseHedgeData {
    pub ufo_order: Option<u8>,
}

impl ParseHedgeData {
    pub(crate) fn parse<'a>() -> impl FnMut((Hedge, &'a DotHedgeData)) -> Result<Self> {
        |(_, h)| {
            let Some(statement) = h
                .statement
                .as_deref()
                .map(str::trim)
                .filter(|s| !s.is_empty())
            else {
                return Ok(ParseHedgeData::default());
            };

            Ok(json5::from_str(statement)?)
        }
    }
}

impl From<&HedgeData> for DotHedgeData {
    fn from(value: &HedgeData) -> Self {
        let payload = ParseHedgeData {
            ufo_order: (!value.ufo_order.autogenerated).then_some(value.ufo_order.value),
        };

        let statement = if payload.ufo_order.is_none() {
            None
        } else {
            Some(json5::to_string(&payload).expect("serializing hedge payload should not fail"))
        };

        DotHedgeData::from(statement)
    }
}

impl From<&ParseHedgeData> for DotHedgeData {
    fn from(value: &ParseHedgeData) -> Self {
        let statement = if value.ufo_order.is_none() {
            None
        } else {
            Some(json5::to_string(value).expect("serializing hedge payload should not fail"))
        };

        DotHedgeData::from(statement)
    }
}

#[cfg(test)]
mod tests {
    use linnet::{half_edge::involution::Hedge, parser::DotHedgeData};

    use super::ParseHedgeData;

    #[test]
    fn parse_json5_hedge_payload() {
        let dot = DotHedgeData::from(Some("{ufo_order: 2}".to_string()));
        let parsed = ParseHedgeData::parse()((Hedge(0), &dot)).unwrap();

        assert_eq!(parsed.ufo_order, Some(2));
    }

    #[test]
    fn serialize_parse_hedge_data_to_json5_statement() {
        let dot = DotHedgeData::from(Some("{ufo_order: 1}".to_string()));
        let parsed = ParseHedgeData::parse()((Hedge(0), &dot)).unwrap();

        let serialized: DotHedgeData = (&parsed).into();
        let statement = serialized.statement.expect("statement should be present");

        assert!(statement.contains("ufo_order"));
    }
}
