use itertools::Itertools;
use linnet::{
    half_edge::involution::{Flow, Hedge},
    parser::DotHedgeData,
};
use serde::{Deserialize, Serialize};
use spenso::structure::{
    Canonicalized, OrderedStructure, TensorStructure,
    representation::{LibraryRep, LibrarySlot, RepName},
    slot::IsAbstractSlot,
};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::atom::AtomCore;
use symbolica::atom::{Atom, FunctionBuilder};

use crate::{graph::edge::PossibleParticle, model::Model, numerator::aind::Aind};

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

#[derive(Clone, Copy)]
pub(crate) enum NumeratorIndexSource {
    /// FeynKit localizes positive UFO indices starting at one.
    Feynkit,
    /// A finalized artifact owns its literal represented indices.
    Literal,
}

impl NumIndices {
    fn from_particle(
        particle: &PossibleParticle,
        model: &Model,
        flow: Flow,
        h: Hedge,
        source: NumeratorIndexSource,
        labels: &BTreeMap<(usize, LibraryRep), BTreeSet<u16>>,
    ) -> Result<Self> {
        let creps = particle.color_reps(model, flow);
        let sreps = particle.spin_reps(model);
        let structure = |representations: Vec<
            spenso::structure::representation::Representation<LibraryRep>,
        >|
         -> Result<OrderedStructure<LibraryRep, Aind>> {
            let mut counts = BTreeMap::new();
            for representation in &representations {
                *counts.entry(representation.rep.base()).or_insert(0usize) += 1;
            }
            let mut ordinals = BTreeMap::new();
            let mut slots = Vec::new();
            for representation in representations {
                let family = representation.rep.base();
                let ordinal = ordinals.entry(family).or_insert(0usize);
                let local = match source {
                    NumeratorIndexSource::Feynkit => Some(u16::try_from(*ordinal + 1)?),
                    NumeratorIndexSource::Literal => match labels.get(&(h.0, family)) {
                        None => None,
                        Some(indices) => {
                            eyre::ensure!(
                                indices.len() == counts[&family],
                                "ambiguous literal tensor interface at half-edge {} for {}: expected {} index labels, found {:?}",
                                h.0,
                                family,
                                counts[&family],
                                indices
                            );
                            indices.iter().nth(*ordinal).copied()
                        }
                    },
                };
                if let Some(local) = local {
                    slots.push(representation.slot(Aind::Hedge(u16::try_from(h.0)?, local)));
                }
                *ordinal += 1;
            }
            let structure: Canonicalized<_> = slots.into_iter().collect();
            Ok(structure.into_canonical())
        };
        Ok(NumIndices {
            color_indices: HedgeIndices::new(structure(creps.external_reps_iter().collect())?),
            spin_indices: HedgeIndices::new(structure(sreps.external_reps_iter().collect())?),
        })
    }

    pub(crate) fn parse<'a>(
        graph: &'a ParseGraph,
        model: &'a Model,
        source: NumeratorIndexSource,
    ) -> impl FnMut(Hedge, &'a ParseHedgeData) -> Result<Self> {
        let mut labels = BTreeMap::<(usize, LibraryRep), BTreeSet<u16>>::new();
        if matches!(source, NumeratorIndexSource::Literal) {
            // Each label explicitly identifies its physical half-edge and
            // representation. No momentum inference or expression expansion
            // participates in this interface certificate.
            let mut inspect = |atom: &Atom| {
                atom.visitor(&mut |part| {
                    if let Ok(slot) = LibrarySlot::<Aind>::try_from(part)
                        && let Aind::Hedge(hedge, local) = slot.aind
                    {
                        labels
                            .entry((usize::from(hedge), slot.rep().rep.base()))
                            .or_default()
                            .insert(local);
                        return false;
                    }
                    true
                });
            };
            for (_, _, vertex) in graph.iter_nodes() {
                if let Some(atom) = &vertex.num {
                    inspect(atom);
                }
            }
            for (_, _, edge) in graph.iter_edges() {
                if let Some(atom) = &edge.data.num {
                    inspect(atom);
                }
            }
            inspect(&graph.global_data.num);
            inspect(&graph.global_data.overall_factor);
            if let Some(atom) = &graph.global_data.projectors {
                inspect(atom);
            }
        }
        move |h, _| {
            let eid = graph[&h];
            let flow = graph.flow(h);
            Self::from_particle(&graph[eid].particle, model, flow, h, source, &labels)
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
    pub(crate) fn edge_color_slots(&self) -> impl Iterator<Item = LibrarySlot<Aind>> + '_ {
        self.num_indices
            .color_indices
            .edge_indices
            .external_structure_iter()
    }

    pub(crate) fn color_kronekers(&self, other: &Self) -> Atom {
        self.edge_color_slots()
            .zip_eq(other.edge_color_slots())
            .map(|(left, right)| {
                assert!(
                    left.rep().matches(&right.rep()),
                    "color representations must match"
                );
                right.rep().id(left.aind, right.aind)
            })
            .product()
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
        (&payload).into()
    }
}

impl From<&ParseHedgeData> for DotHedgeData {
    fn from(value: &ParseHedgeData) -> Self {
        let Some(order) = value.ufo_order else {
            return DotHedgeData::default();
        };

        // Linnet quotes hedge statements when writing `source=`/`sink=`. Keep
        // the compact JSON5 payload unquoted here so serialized finalized
        // artifacts can be parsed again without relying on inference.
        DotHedgeData::default().with_statement(format!("{{ufo_order:{order}}}"))
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

        assert_eq!(statement, "{ufo_order:1}");
    }
}
