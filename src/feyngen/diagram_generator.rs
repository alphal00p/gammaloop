use log::debug;

use super::{FeynGenError, FeynGenOptions};
use crate::{graph::BareGraph, model::Model};
use symbolica::graph::Graph as SymbolicaGraph;

pub struct FeynGen {
    pub options: FeynGenOptions,
}

impl FeynGen {
    pub fn new(options: FeynGenOptions) -> Self {
        Self { options }
    }

    pub fn generate(&self, model: &Model) -> Result<Vec<BareGraph>, FeynGenError> {
        debug!(
            "Generating Feynman diagrams for model {} and process:\n{}",
            model.name, self.options
        );

        let vertex_signatures = model
            .vertex_rules
            .iter()
            .map(|v| {
                v.particles
                    .iter()
                    .map(|p| p.name.as_str())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let mut external_edges = self
            .options
            .initial_pdgs
            .iter()
            .enumerate()
            .map(|(i_initial, pdg)| {
                (
                    i_initial + 1,
                    model
                        .get_particle_from_pdg(*pdg as isize)
                        .clone()
                        .name
                        .clone(),
                )
            })
            .collect::<Vec<_>>();
        external_edges.extend(
            self.options
                .final_pdgs
                .iter()
                .enumerate()
                .map(|(i_final, pdg)| {
                    (
                        self.options.initial_pdgs.len() + i_final + 1,
                        model
                            .get_particle_from_pdg(*pdg as isize)
                            .clone()
                            .name
                            .clone(),
                    )
                })
                .collect::<Vec<_>>(),
        );
        debug!("external_edges = {:?}", external_edges);
        debug!("vertex_signatures = {:?}", vertex_signatures);

        let graphs = SymbolicaGraph::generate(
            &external_edges
                .iter()
                .map(|(i, name)| (*i, name.as_str()))
                .collect::<Vec<_>>(),
            &vertex_signatures,
            None,
            Some(self.options.loop_count_range.1),
            self.options.filters.get_max_bridge(),
            self.options.filters.allow_tadpoles(),
        );

        debug!("Symbolica generated {} graphs", graphs.len());

        Ok(vec![])
    }
}
