use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::atom::Atom;
use three_dimensional_reps::CffGenerationContext;

use crate::{
    cff::{
        CutCFF, CutCFFIndex,
        expression::{OrientationExpression, OrientationID, OrientationSelector},
        orientations::GraphOrientation,
        surface::LinearEnergyExpr,
    },
    graph::{Graph, cuts::CutSet},
    settings::global::OrientationPattern,
    utils::GS,
    uv::approx::OrientationProjection,
};

use super::{Localizer, OrientationIntegrandBranch, OrientationIntegrands, SourceSelectorHosting};

impl<'a> Localizer<'a> {
    pub(crate) fn new(cutset: &'a CutSet, orientation: OrientationProjection<'a>) -> Self {
        Self {
            cutset,
            orientation,
            source_selector_hosting: SourceSelectorHosting::PhysicalPrefix,
        }
    }

    /// Projected local-4D pieces carry independent CFF sums. Their temporary
    /// production host never becomes a selector and cannot discard a
    /// source-local residue map merely because it has no complete physical
    /// extension.
    pub(in crate::uv::approx) fn with_independent_source_sum(mut self) -> Self {
        self.source_selector_hosting = SourceSelectorHosting::IndependentSum;
        self
    }

    fn cff_contract_subgraph(self, graph: &Graph, to_contract: &SuBitGraph) -> SuBitGraph {
        to_contract
            .union(&graph.tree_edges)
            .subtract(&graph.initial_state_cut)
    }

    fn cff(
        self,
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        analysis_numerator: &Atom,
        generation_context: CffGenerationContext,
    ) -> Result<(CutCFF, SuBitGraph)> {
        let contract_subgraph = self.cff_contract_subgraph(graph, to_contract);
        let mut options = self.orientation.cff_options(graph);
        options.cff_generation_context = generation_context;
        // Exact projection applies the user pattern to full production maps.
        // Contracted edges are undirected in the reduced CFF and cannot be
        // filtered against a pattern that still constrains those edges.
        let unfiltered = OrientationPattern::default();
        let orientation_pattern = if self.orientation.exact_orientations().is_some() {
            &unfiltered
        } else {
            self.orientation.orientation_pattern
        };
        // Only production exact maps own numerator-energy capacity. Coarse
        // projectors retain their denominator-only CFF and identity numerator
        // mapping, as used by legacy exports and isolated UV tests.
        let analysis_numerator = self
            .orientation
            .exact_orientations()
            .is_some()
            .then_some(analysis_numerator);
        let cff = if to_contract.is_empty() {
            if let Some(root_expression) = self.orientation.root_expression() {
                graph.cff_from_production_expression(
                    root_expression,
                    self.cutset,
                    orientation_pattern,
                )?
            } else {
                graph.cff(
                    &contract_subgraph,
                    self.cutset,
                    orientation_pattern,
                    &options,
                    analysis_numerator,
                )?
            }
        } else {
            graph.cff(
                &contract_subgraph,
                self.cutset,
                orientation_pattern,
                &options,
                analysis_numerator,
            )?
        };
        self.orientation
            .record_energy_degree_bound_report(&cff.energy_degree_bound_report);
        Ok((cff, contract_subgraph))
    }

    pub(super) fn exact_representatives(
        self,
        graph: &Graph,
        reduced: &OrientationExpression,
        contract_subgraph: &SuBitGraph,
    ) -> Result<Vec<OrientationID>> {
        let production = self
            .orientation
            .exact_orientations()
            .expect("exact representatives are only requested for an exact projector");
        let contracted_edges = graph.paired_edges(contract_subgraph);
        // A UV source's contracted-edge orientations resolve its internal
        // energy residues; they are not outer production-map directions.
        let mut explicit_reduced_orientation = reduced.data.orientation.clone();
        for edge in &contracted_edges {
            if explicit_reduced_orientation
                .iter()
                .any(|(explicit_edge, _)| explicit_edge == *edge)
            {
                explicit_reduced_orientation[*edge] = Orientation::Undirected;
            }
        }
        let surviving_edges = graph
            .as_ref()
            .iter_edges()
            .filter_map(|(pair, edge_id, edge)| {
                (pair.is_paired() && !edge.data.is_dummy && !contracted_edges.contains(&edge_id))
                    .then_some(edge_id)
            })
            .collect::<Vec<_>>();

        let exact_matches = production
            .iter_enumerated()
            .filter(|(_, full)| {
                full.data
                    .orientation
                    .is_compatible_with(&explicit_reduced_orientation)
                    && surviving_edges.iter().all(|edge_id| {
                        let edge = usize::from(*edge_id);
                        full.edge_energy_map
                            .get(edge)
                            .cloned()
                            .unwrap_or_default()
                            .canonical()
                            == reduced
                                .edge_energy_map
                                .get(edge)
                                .cloned()
                                .unwrap_or_default()
                                .canonical()
                    })
            })
            .map(|(id, _)| id)
            .collect::<Vec<_>>();

        if exact_matches.is_empty() {
            let compatible = production
                .iter_enumerated()
                .filter(|(_, full)| {
                    full.data
                        .orientation
                        .is_compatible_with(&explicit_reduced_orientation)
                })
                .collect::<Vec<_>>();
            let first_mismatch = compatible.first().map(|(id, full)| {
                let mismatches = surviving_edges
                    .iter()
                    .filter_map(|edge_id| {
                        let edge = usize::from(*edge_id);
                        let full_energy = full
                            .edge_energy_map
                            .get(edge)
                            .cloned()
                            .unwrap_or_default()
                            .canonical();
                        let reduced_energy = reduced
                            .edge_energy_map
                            .get(edge)
                            .cloned()
                            .unwrap_or_default()
                            .canonical();
                        (full_energy != reduced_energy).then_some((
                            *edge_id,
                            full_energy,
                            reduced_energy,
                        ))
                    })
                    .collect::<Vec<_>>();
                (
                    *id,
                    full.data.label.as_deref(),
                    full.data.numerator_map_index,
                    full.variants
                        .iter()
                        .filter_map(|variant| variant.origin.as_deref())
                        .collect::<Vec<_>>(),
                    mismatches,
                )
            });
            Err(eyre!(
                "no production energy map exactly extends normalized reduced map {} after contracting edges {:?}; reduced provenance (label, numerator map, origins): ({:?}, {:?}, {:?}); {} production orientations match the surviving directions; first mismatch (orientation, label, numerator map, origins, edge/production/reduced): {first_mismatch:#?}",
                GS.orientation_delta(&explicit_reduced_orientation),
                contracted_edges,
                reduced.data.label,
                reduced.data.numerator_map_index,
                reduced
                    .variants
                    .iter()
                    .filter_map(|variant| variant.origin.as_deref())
                    .collect::<Vec<_>>(),
                compatible.len(),
            ))
        } else {
            Ok(exact_matches
                .into_iter()
                .filter(|id| {
                    self.orientation.orientation_pattern.filter_orientation(
                        &production
                            .get(*id)
                            .expect("exact representative belongs to the production map set")
                            .data
                            .orientation,
                    )
                })
                .collect())
        }
    }

    pub(crate) fn source_selector_representatives(
        self,
        graph: &Graph,
        reduced: &OrientationExpression,
        contract_subgraph: &SuBitGraph,
    ) -> Result<Vec<OrientationID>> {
        let production = self
            .orientation
            .exact_orientations()
            .expect("source selectors are only requested for an exact projector");
        let contracted_edges = graph.paired_edges(contract_subgraph);
        // The reduced source map owns numerator energies. Production maps only
        // host its surviving physical graph-edge directions. Exact denominator
        // occurrences may append synthetic directions, but those are algebraic
        // contour pieces rather than selector coordinates and are deliberately
        // excluded from this physical prefix. Contracted and initial-cut edges
        // cannot constrain the outer production host. No LMB or energy-frame
        // reconstruction is involved at this selector boundary.
        let mut explicit_reduced_orientation = EdgeVec::from_iter(
            reduced
                .data
                .orientation
                .iter()
                .take(graph.underlying.n_edges())
                .map(|(_, orientation)| *orientation),
        );
        for edge in &contracted_edges {
            if explicit_reduced_orientation
                .iter()
                .any(|(explicit_edge, _)| explicit_edge == *edge)
            {
                explicit_reduced_orientation[*edge] = Orientation::Undirected;
            }
        }
        // Initial-state cut energies are external to every UV-child contour.
        // Their exact source values still map the numerator, but they cannot
        // constrain the production orientation which hosts that child.
        for (_, edge, _) in graph.iter_edges_of(&graph.initial_state_cut) {
            explicit_reduced_orientation[edge] = Orientation::Undirected;
        }
        let compatible = production
            .iter_enumerated()
            .filter(|(_, full)| {
                full.data
                    .orientation
                    .is_compatible_with(&explicit_reduced_orientation)
            })
            .map(|(id, _)| id)
            .collect::<Vec<_>>();

        if compatible.is_empty() {
            return Err(eyre!(
                "no valid global orientation matches reduced orientation {}",
                GS.orientation_delta(&explicit_reduced_orientation)
            ));
        }

        Ok(compatible
            .into_iter()
            .filter(|id| {
                self.orientation.orientation_pattern.filter_orientation(
                    &production
                        .get(*id)
                        .expect("source selector belongs to the production map set")
                        .data
                        .orientation,
                )
            })
            .collect())
    }

    fn coarse_representatives(
        self,
        reduced_orientation: &EdgeVec<Orientation>,
    ) -> Result<Vec<OrientationID>> {
        let ids = self.orientation.orientation_ids();
        if ids.len() == 1 && self.orientation.orientation(ids[0]).is_none() {
            return Ok(ids);
        }
        let compatible = ids
            .into_iter()
            .filter(|id| {
                self.orientation
                    .orientation(*id)
                    .is_some_and(|orientation| orientation.is_compatible_with(reduced_orientation))
            })
            .collect::<Vec<_>>();

        if compatible.is_empty() {
            Err(eyre!(
                "no valid global orientation matches reduced orientation {}",
                GS.orientation_delta(reduced_orientation)
            ))
        } else {
            Ok(compatible)
        }
    }

    /// Assign a reduced CFF map to one complete production-map key. Ordinary
    /// 3D evaluation carries that key as the opaque sparse selector; explicit-
    /// sum evaluation keeps the same reduced residue once with the key used as
    /// mapping metadata only.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn localized_orientation_terms(
        self,
        graph: &Graph,
        reduced: &OrientationExpression,
        reduced_expression: &Atom,
        contract_subgraph: &SuBitGraph,
        internal_edges: &[EdgeIndex],
        valid_production_ids: Option<&BTreeSet<OrientationID>>,
        production_orientation_id: Option<OrientationID>,
        source_edge_energy_map: Option<&[LinearEnergyExpr]>,
    ) -> Result<Vec<(OrientationID, Atom)>> {
        if let Some(id) = production_orientation_id {
            if valid_production_ids.is_some_and(|valid| !valid.contains(&id)) {
                return Ok(Vec::new());
            }
            // A stored root residue is already diagonal in its complete
            // production-map key. Keep that key in the sparse branch sidecar:
            // the UV Taylor operators act on the branch body and therefore see
            // the exact factorized equivalent of `sigma(id) * body` without
            // mistaking the coarser physical edge directions for the selector.
            // Reduced/new CFF terms have no stored production ID and continue
            // through representative reconstruction below.
            return Ok(self
                .orientation
                .orientation(id)
                .filter(|orientation| {
                    self.orientation
                        .orientation_pattern
                        .filter_orientation(orientation)
                })
                .map(|_| vec![(id, reduced_expression.clone())])
                .unwrap_or_default());
        }
        let candidate_representatives = if self.orientation.exact_orientations().is_some() {
            if source_edge_energy_map.is_none() {
                self.exact_representatives(graph, reduced, contract_subgraph)?
            } else {
                match self.source_selector_representatives(graph, reduced, contract_subgraph) {
                    Ok(representatives) if !representatives.is_empty() => representatives,
                    Ok(_) | Err(_)
                        if self.source_selector_hosting
                            == SourceSelectorHosting::IndependentSum =>
                    {
                        // A compatible physical-prefix host is the most stable
                        // mapping frame. Some generalized numerator samples do
                        // not define a complete production orientation at all;
                        // the projected local-4D lane must still keep them once
                        // under a permitted deterministic bookkeeping host.
                        self.orientation.orientation_ids()
                    }
                    Ok(representatives) => representatives,
                    Err(error) => return Err(error),
                }
            }
        } else {
            self.coarse_representatives(&reduced.data.orientation)?
        };
        let mut representatives = candidate_representatives
            .iter()
            .copied()
            .filter(|id| valid_production_ids.is_none_or(|valid| valid.contains(id)))
            .collect::<Vec<_>>();
        if representatives.is_empty()
            && source_edge_energy_map.is_some()
            && (self.source_selector_hosting == SourceSelectorHosting::IndependentSum
                || self.orientation.explicit_orientation_sum_only)
        {
            // A selected cut can exclude every otherwise compatible physical
            // host. An independently summed source—or any source after the
            // orientation selectors have explicitly been summed—still owns
            // this residue through its exact energy map. Use a remaining
            // cut-valid ID only as deterministic bookkeeping rather than
            // dropping that residue-map sample.
            representatives = self
                .orientation
                .orientation_ids()
                .into_iter()
                .filter(|id| valid_production_ids.is_none_or(|valid| valid.contains(id)))
                .collect();
        }
        let mut selector_edges = graph
            .as_ref()
            .iter_edges()
            .filter_map(|(pair, edge_id, edge)| {
                (pair.is_paired() && !edge.data.is_dummy).then_some(edge_id)
            })
            .collect::<Vec<_>>();
        for edge_id in internal_edges {
            if !selector_edges.contains(edge_id) {
                selector_edges.push(*edge_id);
            }
        }
        let representative_score = |id: &OrientationID| {
            self.orientation.orientation(*id).map(|orientation| {
                // Generalized numerator sampling can add under-resolved maps
                // to the CFF carrier. They own numerator values, but are not
                // additional physical orientation channels, so resolve as many
                // selector directions as possible before applying the ordinary
                // deterministic representative ordering.
                let directed_count = selector_edges
                    .iter()
                    .filter(|edge_id| orientation[**edge_id] != Orientation::Undirected)
                    .count();
                (directed_count, orientation.score(&selector_edges))
            })
        };
        if self.orientation.explicit_orientation_sum_only {
            let Some(representative) = representatives
                .iter()
                .copied()
                .max_by_key(representative_score)
            else {
                return Ok(Vec::new());
            };
            return Ok(vec![(representative, reduced_expression.clone())]);
        }
        let representatives = representatives
            .iter()
            .copied()
            .max_by_key(representative_score)
            .into_iter()
            .collect::<Vec<_>>();
        // Exact production maps keep their complete selector key in the sparse
        // sidecar. Coarse diagnostic/export projectors have no such key and
        // retain their historical physical-theta selector in the atom.
        let exact_map_selector = self.orientation.exact_orientations().is_some();
        Ok(representatives
            .into_iter()
            .map(|representative| {
                let selector = if exact_map_selector {
                    Atom::one()
                } else {
                    self.orientation
                        .orientation(representative)
                        .map(|orientation| orientation.orientation_thetas())
                        .unwrap_or_else(Atom::one)
                };
                (representative, reduced_expression.clone() * selector)
            })
            .collect())
    }

    #[cfg(test)]
    pub(super) fn localized_orientation_term(
        self,
        reduced_expression: &Atom,
        reduced_orientation: &EdgeVec<Orientation>,
        internal_edges: &[EdgeIndex],
    ) -> Result<Atom> {
        let representative = self
            .coarse_representatives(reduced_orientation)?
            .into_iter()
            .max_by_key(|id| {
                self.orientation
                    .orientation(*id)
                    .map(|orientation| orientation.score(internal_edges))
            });
        let reduced_selector = reduced_orientation.orientation_thetas();
        Ok(representative
            .map(|id| {
                reduced_expression.clone()
                    * reduced_selector
                    * self
                        .orientation
                        .orientation(id)
                        .map(|orientation| {
                            orientation.internal_orientation_selector(internal_edges)
                        })
                        .unwrap_or_else(Atom::one)
            })
            .unwrap_or(Atom::Zero))
    }

    pub(crate) fn projected_cff(
        self,
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        analysis_numerator: &Atom,
        generation_context: CffGenerationContext,
    ) -> Result<OrientationIntegrands> {
        let (cff, contract_subgraph) =
            self.cff(graph, to_contract, analysis_numerator, generation_context)?;
        // A generalized source map can have several resolved production
        // extensions, but only some of them support the selected Cutkosky
        // residue. Restrict selector hosts to those admissible production IDs,
        // separately for every raised order. The root CFF is already the
        // authoritative selected source; reduced CFFs must consult that root.
        let has_selected_residue = self.cutset.residue_selector.lu.is_some()
            || self.cutset.residue_selector.left_th_cut.is_some()
            || self.cutset.residue_selector.right_th_cut.is_some();
        let valid_production_ids = if has_selected_residue && to_contract.is_empty() {
            Some(
                cff.terms
                    .iter()
                    .map(|(index, term)| {
                        let ids = term
                            .orientations
                            .iter()
                            .filter(|orientation| !orientation.orientation.variants.is_empty())
                            .filter_map(|orientation| orientation.production_orientation_id)
                            .collect::<BTreeSet<_>>();
                        (*index, ids)
                    })
                    .collect::<BTreeMap<_, _>>(),
            )
        } else if has_selected_residue
            && let Some(root_expression) = self.orientation.root_expression()
        {
            let unfiltered = OrientationPattern::default();
            Some(
                graph
                    .cff_from_production_expression(root_expression, self.cutset, &unfiltered)?
                    .terms
                    .into_iter()
                    .map(|(index, term)| {
                        let ids = term
                            .orientations
                            .into_iter()
                            .filter(|orientation| !orientation.orientation.variants.is_empty())
                            .filter_map(|orientation| orientation.production_orientation_id)
                            .collect::<BTreeSet<_>>();
                        (index, ids)
                    })
                    .collect::<BTreeMap<_, _>>(),
            )
        } else {
            None
        };
        let no_valid_production_ids = BTreeSet::new();
        let indices = cff.terms.keys().copied().collect::<Vec<_>>();
        let production_prefactor = Atom::num(cff.production_prefactor_factor());
        let ids = self.orientation.orientation_ids();
        if ids.is_empty() {
            return Err(eyre!(
                "orientation pattern selects no production energy maps"
            ));
        }
        #[allow(clippy::type_complexity)]
        let mut terms: Vec<(
            OrientationID,
            Option<Vec<LinearEnergyExpr>>,
            BTreeMap<CutCFFIndex, Atom>,
        )> = Vec::new();
        let internal_edges = if self.orientation.exact_orientations().is_some() {
            // Exact full maps must also recover selectors for graph tree edges
            // contracted by the CFF source, not only the explicit UV subgraph.
            graph.paired_edges(&contract_subgraph)
        } else {
            // Preserve the ordinary coarse-localization convention.
            graph.paired_edges(to_contract)
        };
        for (index, cff_term) in cff.terms {
            let valid_production_ids = valid_production_ids.as_ref().map(|valid_by_index| {
                valid_by_index
                    .get(&index)
                    .unwrap_or(&no_valid_production_ids)
            });
            for reduced in cff_term.orientations {
                // A generalized carrier map can be under-resolved even for a
                // stored root residue. Keep that branch-owned numerator map
                // separate from the resolved production orientation selected
                // solely to partition the ordinary runtime orientation sum.
                let source_edge_energy_map = self
                    .orientation
                    .exact_orientations()
                    .is_some()
                    .then(|| reduced.orientation.edge_energy_map.clone());
                let reduced_expression = &reduced.expression * &production_prefactor;
                let localized = self.localized_orientation_terms(
                    graph,
                    &reduced.orientation,
                    &reduced_expression,
                    &contract_subgraph,
                    &internal_edges,
                    valid_production_ids,
                    reduced.production_orientation_id,
                    source_edge_energy_map.as_deref(),
                )?;
                for (id, expression) in localized {
                    if !terms.iter().any(|(selector_id, energy_map, _)| {
                        *selector_id == id && energy_map == &source_edge_energy_map
                    }) {
                        terms.push((
                            id,
                            source_edge_energy_map.clone(),
                            indices.iter().map(|index| (*index, Atom::Zero)).collect(),
                        ));
                    }
                    *terms
                        .iter_mut()
                        .find(|(selector_id, energy_map, _)| {
                            *selector_id == id && energy_map == &source_edge_energy_map
                        })
                        .and_then(|(_, _, integrands)| integrands.get_mut(&index))
                        .expect("all projected CFF branch keys were initialized") += expression;
                }
            }
        }

        Ok(OrientationIntegrands(
            terms
                .into_iter()
                .map(|(selector_id, source_edge_energy_map, integrands)| {
                    OrientationIntegrandBranch {
                        selector_id,
                        source_edge_energy_map,
                        integrands: integrands.into_iter().collect(),
                    }
                })
                .collect(),
        ))
    }

    pub(crate) fn map_numerator(
        self,
        graph: &Graph,
        orientation_id: OrientationID,
        source_edge_energy_map: Option<&[LinearEnergyExpr]>,
        numerator: &Atom,
    ) -> Result<Atom> {
        self.orientation
            .map_numerator(graph, orientation_id, source_edge_energy_map, numerator)
    }

    pub(crate) fn uses_exact_maps(self) -> bool {
        self.orientation.exact_orientations().is_some()
    }

    /// Materialize the opaque selector carried by an exact production-map
    /// branch. Explicit sums and coarse projectors have no runtime map key.
    #[cfg(test)]
    pub(crate) fn residue_map_key_selector(self, id: OrientationID) -> Atom {
        if self.uses_exact_maps() && !self.orientation.explicit_orientation_sum_only {
            id.atom()
        } else {
            Atom::one()
        }
    }
}
