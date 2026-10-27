#[cfg(test)]
use itertools::Itertools;
use linnet::half_edge::involution::{Flow, Hedge};
use spenso::structure::representation::{LibraryRep, LibrarySlot, RepName};
use spenso::structure::slot::IsAbstractSlot;
use spenso::structure::{Canonicalized, OrderedStructure, TensorStructure};
use std::collections::BTreeMap;
use symbolica::atom::{Atom, FunctionBuilder};

use crate::graph::edge::PossibleParticle;
use crate::model::Model;
use crate::numerator::aind::Aind;

use super::Autogen;
use super::parse::ParseGraph;

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
    /// The indices of half-edge `h` in FeynKit's frame: the k-th index of each
    /// representation family is `hedge(h,k)`, counted from one.
    fn from_particle(
        particle: &PossibleParticle,
        model: &Model,
        flow: Flow,
        h: Hedge,
    ) -> Result<Self> {
        let creps = particle.color_reps(model, flow);
        let sreps = particle.spin_reps(model);
        let structure = |representations: Vec<
            spenso::structure::representation::Representation<LibraryRep>,
        >|
         -> Result<OrderedStructure<LibraryRep, Aind>> {
            let mut ordinals = BTreeMap::new();
            let mut slots = Vec::new();
            for representation in representations {
                let ordinal = ordinals.entry(representation.rep.base()).or_insert(0usize);
                *ordinal += 1;
                slots.push(
                    representation.slot(Aind::Hedge(u16::try_from(h.0)?, u16::try_from(*ordinal)?)),
                );
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
    ) -> impl FnMut(Hedge, &'a ParseHedgeData) -> Result<Self> {
        move |h, _| {
            let eid = graph[&h];
            let flow = graph.flow(h);
            Self::from_particle(&graph[eid].particle, model, flow, h)
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
    #[cfg(test)]
    pub(crate) fn edge_color_slots(&self) -> impl Iterator<Item = LibrarySlot<Aind>> + '_ {
        self.num_indices
            .color_indices
            .edge_indices
            .external_structure_iter()
    }

    #[cfg(test)]
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

/// The interaction slot FeynKit assigned to a half-edge.
#[derive(Debug, Clone)]
pub struct ParseHedgeData {
    pub ufo_order: u8,
}
