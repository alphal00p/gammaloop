use linnet::half_edge::involution::{Flow, Hedge};
use spenso::structure::{
    Canonicalized, OrderedStructure, TensorStructure,
    representation::{LibraryRep, LibrarySlot},
    slot::IsAbstractSlot,
};
use symbolica::atom::{Atom, FunctionBuilder};

use crate::{graph::edge::PossibleParticle, model::Model, numerator::aind::Aind};

use super::{Autogen, parse::ParseGraph};

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
    fn from_particle(particle: &PossibleParticle, model: &Model, flow: Flow, h: Hedge) -> Self {
        let creps = particle.color_reps(model, flow);
        let sreps = particle.spin_reps(model);

        let mut init = 0;
        let mut last = None;
        let color_structure: Canonicalized<_> = creps
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

        let spin_structure: Canonicalized<_> = sreps
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
            color_indices: HedgeIndices::new(color_structure.into_canonical()),
            spin_indices: HedgeIndices::new(spin_structure.into_canonical()),
        }
    }

    pub(crate) fn parse<'a>(
        graph: &'a ParseGraph,
        model: &'a Model,
    ) -> impl FnMut(Hedge, &'a ParseHedgeData) -> Self {
        |h, _| {
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
