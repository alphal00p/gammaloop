use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetOps};
use symbolica::atom::Atom;

use crate::{
    cff::CutCFFIndex,
    debug_tags,
    graph::{Graph, GraphThreeDSource},
    utils::GS,
    uv::{DeferredIntegrands, Integrands, UltravioletGraph, approx::local_4d::Full4dCts},
};

use super::{Local3DCts, Localizer};

impl Localizer<'_> {
    #[cfg(test)]
    pub(crate) fn project_4d(
        self,
        source: &Full4dCts,
        graph: &mut Graph,
        current_spinney: &SuBitGraph,
    ) -> Result<Local3DCts> {
        let indices = self.cutset.residue_selector.generate_allowed_keys();
        let mut projected = DeferredIntegrands::from_indices(indices.iter().copied());
        let mut exact_cff_generation_cache =
            crate::cff::generation::ExactCffGenerationCache::default();
        let reduced = graph
            .full_filter()
            .subtract(current_spinney)
            .subtract(&graph.initial_state_cut);
        let outside_numerator = graph
            .numerator(&reduced, current_spinney)
            .get_single_atom()
            .expect("Graph numerator should be available")
            * graph.global_atom();

        struct ExactProjectionTerm {
            analysis_numerator: Atom,
            active_denominators: Vec<crate::graph::FourDDenominator>,
            residual_factor: Atom,
            frozen_localizer: Atom,
        }

        let options = self.orientation.cff_options(graph);
        let uv_edges = graph
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge_id, edge)| {
                (pair.is_paired() && !edge.data.is_dummy && current_spinney.includes(&pair))
                    .then_some(edge_id)
            })
            .collect::<Vec<_>>();
        let mut exact_terms = Vec::new();

        for sector in source.sectors() {
            let frozen_localizer = sector
                .frozen_lmbs()
                .iter()
                .fold(Atom::one(), |product, lmb| {
                    product * GS.localizing_integrand(lmb)
                });
            if !sector.active_components.is_empty() {
                // Typed Taylor sectors carry their own component frames. They
                // must never fall through to the unframed whole-source oracle,
                // whose production LMB can restore a loop direction already
                // demoted by the local Taylor operation.
                let active =
                    self.project_factorized_taylor_sector(graph, sector, &outside_numerator)?;
                for index in &indices {
                    projected.push(*index, &active * &frozen_localizer)?;
                }
                continue;
            }
            for term in sector.physical_terms()? {
                // Analyze and later map the numerator owned by this exact 4D term
                // together with the factors grown outside its owning spinney. Keep
                // the product factorized throughout both operations.
                let analysis_numerator = &term.numerator * &outside_numerator;
                let mut active_denominators = Vec::new();
                let mut residual_factor = Atom::one();
                for denominator in term.denominators {
                    let is_uv = current_spinney.includes(&graph[&denominator.source_edge].1);
                    if denominator.depends_on_loop(graph, is_uv)? {
                        active_denominators.push(denominator);
                    } else {
                        // Keep each exact occurrence. This deliberately does not
                        // group by source edge: equal edge IDs can carry distinct
                        // full expressions after a 4D Taylor expansion.
                        residual_factor /= denominator.full_expr;
                    }
                }
                if active_denominators.is_empty() {
                    if !indices.contains(&CutCFFIndex::new_all_none()) {
                        continue;
                    }
                    // A denominatorless exact source has no local energy residue
                    // and therefore no orientation selector. Its complete
                    // numerator must be independent of every loop energy; the
                    // empty exact source proves that before setting its inactive
                    // temporal placeholders to zero.
                    let exact_source = GraphThreeDSource::from_exact_denominators(graph, &[])
                        .map_err(|error| {
                            eyre!("could not build denominatorless exact 4D source: {error}")
                        })?;
                    // Initial-cut and tree-edge energies are external to the empty
                    // CFF source. The mapper below still projects their complete
                    // external-affine physical energy into the factorized
                    // numerator, while `residual_factor` owns every corresponding
                    // denominator occurrence.
                    let cff_external_edges = graph
                        .iter_edges_of(&graph.initial_state_cut)
                        .chain(graph.iter_edges_of(&graph.tree_edges))
                        .map(|(_, edge_id, _)| edge_id)
                        .collect::<Vec<_>>();
                    let mapper = exact_source
                        .exact_source_energy_mapper()
                        .expect("empty exact source has an owned parent-energy mapper");
                    let candidates = mapper.equivalent_energy_candidates(std::iter::empty())?;
                    let plan = graph
                        .plan_numerator_energy_assignment_in_atom_excluding(
                            &analysis_numerator,
                            cff_external_edges,
                            &candidates,
                        )
                        .map_err(|error| {
                            eyre!(
                                "denominatorless exact 4D term has an active numerator EMR energy but no denominator residue can carry it: {error}"
                            )
                        })?;
                    let mapped_numerator = mapper.map_planned_numerator(&[], &[], &plan)?;
                    projected.push(
                        CutCFFIndex::new_all_none(),
                        mapped_numerator * &residual_factor * &frozen_localizer,
                    )?;
                    continue;
                }

                exact_terms.push(ExactProjectionTerm {
                    analysis_numerator,
                    active_denominators,
                    residual_factor,
                    frozen_localizer: frozen_localizer.clone(),
                });
            }
        }

        // Plan every term before generating any CFF. Canonically equivalent
        // UV topologies can require different sparse numerator supports; the
        // cache joins those requirements by algebraic energy channel and
        // supplies one minimax capacity envelope to their shared generation.
        for term in &exact_terms {
            let exact_source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
                graph,
                &term.active_denominators,
                uv_edges.iter().copied(),
            )?;
            graph.register_3d_expression_for_4d_term(
                &exact_source,
                &options,
                &term.analysis_numerator,
                &mut exact_cff_generation_cache,
            )?;
        }

        for term in exact_terms {
            debug_tags!(#generation, #uv, #cff, #term, #inspect;
                graph = %graph.name,
                denominator_count = term.active_denominators.len(),
                log.numerator = term.analysis_numerator,
                file.denominators = ?term.active_denominators,
                "Projecting exact 4D counterterm sector"
            );

            // The exact source CFF resolves contracted directions internally.
            // Its orientations are an explicit source-local sum, so a
            // production orientation pattern is deliberately never applied.
            let (cff, _contract_subgraph) = graph.cff_from_4d_denominators_in_uv_edges(
                &term.active_denominators,
                uv_edges.iter().copied(),
                self.cutset,
                &options,
                &term.analysis_numerator,
                Some(&mut exact_cff_generation_cache),
            )?;
            self.orientation
                .record_energy_degree_bound_report(&cff.energy_degree_bound_report);
            let production_prefactor = Atom::num(cff.production_prefactor_factor());

            for (index, cff_term) in cff.terms {
                for orientation in &cff_term.orientations {
                    let mapped_numerator = cff_term
                        .map_exact_source_numerator(&orientation.orientation)
                        .map_err(|error| {
                            eyre!(
                                "{error}; exact 4D term active denominators are {:?}",
                                term.active_denominators
                            )
                        })?;
                    projected.push(
                        index,
                        orientation.expression.clone()
                            * &production_prefactor
                            * mapped_numerator
                            * &term.residual_factor
                            * &term.frozen_localizer,
                    )?;
                }
            }
        }

        debug_tags!(#generation, #uv, #cff, #summary;
            distinct_exact_topologies = exact_cff_generation_cache.len(),
            "Reused canonically equivalent exact 4D CFF topologies"
        );

        Ok(Local3DCts::from_projected(Integrands::from_deferred(
            projected,
        )))
    }
}
