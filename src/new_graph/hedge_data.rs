use itertools::Itertools;
use linnet::{
    half_edge::{
        involution::{Hedge, Orientation},
        EdgeAccessors,
    },
    parser::DotHedgeData,
};
use spenso::structure::{
    representation::{LibraryRep, LibrarySlot},
    slot::IsAbstractSlot,
    OrderedStructure, PermutedStructure, TensorStructure,
};
use symbolica::atom::{Atom, FunctionBuilder};

use crate::numerator::aind::Aind;

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
    pub(crate) fn parse<'a>(graph: &'a ParseGraph) -> impl FnMut(Hedge, &'a ParseHedge) -> Self {
        |h, _| {
            let eid = graph[&h];
            let flow = graph.flow(h);
            let orientation = graph.orientation(h);

            let flow = match orientation {
                Orientation::Default => flow,
                Orientation::Reversed => -flow,
                Orientation::Undirected => flow,
            };

            let creps = graph[eid].particle.color_reps(flow);
            let sreps = graph[eid].particle.spin_reps();

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
    }
}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct NumHedgeData {
    pub num_indices: NumIndices,
    pub node_order: u8,
}

impl NumHedgeData {
    pub(crate) fn edge_slots<'a>(&'a self) -> impl Iterator<Item = LibrarySlot<Aind>> + 'a {
        self.edge_color_slots().chain(self.edge_spin_slots())
    }

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

    pub(crate) fn vertex_slots<'a>(&'a self) -> impl Iterator<Item = LibrarySlot<Aind>> + 'a {
        self.vertex_color_slots().chain(self.vertex_spin_slots())
    }

    pub(crate) fn vertex_color_slots<'a>(&'a self) -> impl Iterator<Item = LibrarySlot<Aind>> + 'a {
        self.num_indices
            .color_indices
            .vertex_indices
            .external_structure_iter()
    }

    pub(crate) fn vertex_spin_slots<'a>(&'a self) -> impl Iterator<Item = LibrarySlot<Aind>> + 'a {
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

        return color;
    }
}

#[derive(Debug, Clone, Default)]
pub struct ParseHedge {
    hedge_id: Option<usize>,
}

impl ParseHedge {
    pub(crate) fn parse<'a>() -> impl FnMut((Hedge, &'a DotHedgeData)) -> Result<Self> {
        |(i, h)| {
            let hedge_id = h
                .statement
                .as_ref()
                .map(|s| s.parse::<usize>().ok())
                .flatten();
            Ok(ParseHedge { hedge_id })
        }
    }
}

impl From<&NumHedgeData> for DotHedgeData {
    fn from(value: &NumHedgeData) -> Self {
        let h = DotHedgeData::from(Some(value.node_order.to_string()));
        h
    }
}
