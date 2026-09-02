use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::CffGenerationContext;

use crate::{
    cff::{
        CutCFF, CutCFFIndex,
        expression::{OrientationExpression, OrientationID, OrientationSelector},
        orientations::GraphOrientation,
        surface::{GammaLoopLinearEnergyExpr, LinearEnergyExpr},
    },
    graph::{Graph, GraphThreeDSource, LoopMomentumBasis, cuts::CutSet},
    momentum::SignOrZero,
    settings::global::OrientationPattern,
    utils::GS,
    uv::approx::OrientationProjection,
};

use super::{Localizer, OrientationIntegrandBranch, OrientationIntegrands};

impl<'a> Localizer<'a> {
    pub(crate) fn new(cutset: &'a CutSet, orientation: OrientationProjection<'a>) -> Self {
        Self {
            cutset,
            orientation,
        }
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

    #[allow(clippy::type_complexity)]
    pub(crate) fn source_selector_representatives(
        self,
        graph: &Graph,
        source_lmb: &LoopMomentumBasis,
        reduced: &OrientationExpression,
        contract_subgraph: &SuBitGraph,
        physical_source_frame: Option<(&[(EdgeIndex, Atom)], &[EdgeIndex])>,
    ) -> Result<Vec<OrientationID>> {
        let production = self
            .orientation
            .exact_orientations()
            .expect("source selectors are only requested for an exact projector");
        if let Some((physical_source_energies, physical_boundary_edges)) = physical_source_frame {
            {
                // A generalized production map owns two distinct coordinate
                // systems: `loop_energy_map` retains the physical contour, while
                // `edge_energy_map` may contain synthetic +/-M/zero samples used
                // only to evaluate a higher-rank numerator. Reconstruct the former
                // in the same contracted source frame used for production CFF
                // generation before comparing physical selector hosts.
                let production_contract_subgraph =
                    graph.tree_edges.subtract(&graph.initial_state_cut);
                let production_contract_edges = graph.paired_edges(&production_contract_subgraph);
                let production_source = GraphThreeDSource::new(graph, &production_contract_edges)?;
                let physical_edges = physical_boundary_edges.to_vec();
                let production_physical_energy_maps = production
                .iter()
                .map(|full| {
                    physical_edges
                        .iter()
                        .map(|edge| {
                            let coordinates = production_source
                                .reconstructible_outer_loop_coordinates(*edge)
                                .ok_or_else(|| {
                                    eyre!(
                                        "production edge {} cannot be reconstructed in the contracted physical source frame",
                                        usize::from(*edge),
                                    )
                                })?;
                            if coordinates.len() != full.loop_energy_map.len() {
                                return Err(eyre!(
                                    "production edge {} has {} outer coordinates for {} physical loop-energy maps",
                                    usize::from(*edge),
                                    coordinates.len(),
                                    full.loop_energy_map.len(),
                                ));
                            }
                            let mut energy = coordinates.iter().zip(&full.loop_energy_map).fold(
                                LinearEnergyExpr::zero(),
                                |sum, (coefficient, loop_energy)| {
                                    sum + loop_energy.clone().scale_rational(coefficient.clone())
                                },
                            );
                            let signature = &graph.loop_momentum_basis.edge_signatures[*edge];
                            for (external_edge, sign) in graph
                                .loop_momentum_basis
                                .ext_edges
                                .iter()
                                .zip(&signature.external)
                            {
                                let coefficient = match sign {
                                    SignOrZero::Zero => 0,
                                    SignOrZero::Plus => 1,
                                    SignOrZero::Minus => -1,
                                };
                                energy = energy
                                    + LinearEnergyExpr::external(*external_edge, coefficient);
                            }
                            Ok((*edge, energy.canonical().to_atom_gs(&[])))
                        })
                        .collect::<Result<BTreeMap<_, _>>>()
                })
                .collect::<Result<Vec<_>>>()?;
                let mut has_physical_pole_carrier = false;
                let compatible = production
                    .iter_enumerated()
                    .zip(&production_physical_energy_maps)
                    .filter(|((_, full), full_physical_energies)| {
                        let mut saw_pole_carrier = false;
                        let matches_pole_carriers =
                            physical_source_energies.iter().all(|(edge, energy)| {
                                let physical_energy = physical_boundary_edges.iter().fold(
                                    energy.clone(),
                                    |physical_energy, boundary_edge| {
                                        let production_energy =
                                            full_physical_energies[boundary_edge].clone();
                                        physical_energy
                                            .replace(GS.emr_mom(*boundary_edge, GS.cind(0)))
                                            .with(production_energy.clone())
                                            .replace(crate::utils::external_energy_atom_from_index(
                                                *boundary_edge,
                                            ))
                                            .with(production_energy)
                                    },
                                );
                                let positive_pole = GS.ose(*edge);
                                let is_positive = (physical_energy.clone() - &positive_pole)
                                    .expand()
                                    .is_zero();
                                let is_negative = (physical_energy.clone() + &positive_pole)
                                    .expand()
                                    .is_zero();
                                if !is_positive && !is_negative {
                                    // A dependent source edge is not the residue carrier
                                    // of this term. Its affine energy belongs in the CFF
                                    // denominator and cannot constrain a global line
                                    // orientation, which is represented by +/-OSE(edge).
                                    return true;
                                }
                                saw_pole_carrier = true;
                                // The source-local map owns the numerator
                                // coordinates. A production map only hosts
                                // the corresponding exact residue-map key, so its
                                // affine energy need not reproduce the source
                                // pole coordinate.  Compare only the physical
                                // pole direction with that sector.
                                matches!(
                                    (is_positive, is_negative, full.data.orientation[*edge]),
                                    (true, false, Orientation::Default)
                                        | (false, true, Orientation::Reversed)
                                )
                            });
                        has_physical_pole_carrier |= saw_pole_carrier;
                        saw_pole_carrier && matches_pole_carriers
                    })
                    .map(|((id, _), _)| id)
                    .collect::<Vec<_>>();
                if has_physical_pole_carrier {
                    if compatible.is_empty() {
                        return Err(eyre!(
                            "no production orientation matches the exact physical pole-carrier energies {:?} after substituting physical boundary edges {:?}",
                            physical_source_energies,
                            physical_boundary_edges,
                        ));
                    }
                    return Ok(compatible
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
                        .collect());
                }
            }
            // This residue has no simple physical pole among the term's
            // singleton hard channels. Its poles therefore belong to a
            // repeated hard channel, whose synthetic occurrence directions
            // must be summed before selecting one physical production host.
            // Fall through to physical-edge-prefix compatibility below.
        }
        let contracted_edges = graph.paired_edges(contract_subgraph);
        // The reduced source map owns numerator energies, while its loop lift
        // still owns the residue direction in the exact LMB which generated
        // that source. A generalized contact can set all physical sampling
        // energies to zero without erasing the surviving pole: reconstruct only
        // exact +/-OSE directions from that LMB so the two opposite lower
        // residues are not assigned to the same residue-map-key host.
        // Production maps then extend those directions into full sectors; the
        // zero source map remains authoritative for numerator evaluation.
        // Exact denominator occurrences append synthetic orientation entries.
        // They describe algebraic contour pieces, not extra physical graph
        // directions, and only their complete sum has the powered-denominator
        // normalization.  Restrict production-host compatibility to the
        // original graph edge prefix; the complete exact source map remains
        // authoritative for numerator sampling below.
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
        let external_energy_map = source_lmb
            .ext_edges
            .iter()
            .map(|edge| LinearEnergyExpr::external(*edge, 1))
            .collect::<Vec<_>>();
        for edge in graph
            .as_ref()
            .iter_edges()
            .filter_map(|(pair, edge, data)| {
                (pair.is_paired() && !data.data.is_dummy && !contracted_edges.contains(&edge))
                    .then_some(edge)
            })
        {
            let signature = &source_lmb.edge_signatures[edge];
            if signature.internal.len() != reduced.loop_energy_map.len() {
                continue;
            }
            let Some(edge_energy) = signature
                .try_compute_momentum(&reduced.loop_energy_map, &external_energy_map)
                .map(LinearEnergyExpr::canonical)
            else {
                continue;
            };
            if explicit_reduced_orientation[edge] != Orientation::Undirected {
                continue;
            }
            if edge_energy == LinearEnergyExpr::ose(edge, 1) {
                explicit_reduced_orientation[edge] = Orientation::Default;
            } else if edge_energy == LinearEnergyExpr::ose(edge, -1) {
                explicit_reduced_orientation[edge] = Orientation::Reversed;
            }
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
        source_lmb: &LoopMomentumBasis,
        reduced: &OrientationExpression,
        reduced_expression: &Atom,
        contract_subgraph: &SuBitGraph,
        internal_edges: &[EdgeIndex],
        valid_production_ids: Option<&BTreeSet<OrientationID>>,
        production_orientation_id: Option<OrientationID>,
        source_edge_energy_map: Option<&[LinearEnergyExpr]>,
        physical_source_frame: Option<(&[(EdgeIndex, Atom)], &[EdgeIndex])>,
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
            if source_edge_energy_map.is_some() {
                self.source_selector_representatives(
                    graph,
                    source_lmb,
                    reduced,
                    contract_subgraph,
                    physical_source_frame,
                )?
            } else {
                self.exact_representatives(graph, reduced, contract_subgraph)?
            }
        } else {
            self.coarse_representatives(&reduced.data.orientation)?
        };
        let representatives = candidate_representatives
            .iter()
            .copied()
            .filter(|id| valid_production_ids.is_none_or(|valid| valid.contains(id)))
            .collect::<Vec<_>>();
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
        // The generated source's loop map is the physical contour authority.
        // Generalized numerator interpolation changes only `edge_energy_map`,
        // so all of its +/-/zero samples must inherit the selector obtained by
        // lifting this unchanged loop map through the same contracted source
        // coordinates that were used during CFF generation.
        let contracted_source = self
            .orientation
            .exact_orientations()
            .is_some()
            .then(|| GraphThreeDSource::new(graph, &internal_edges))
            .transpose()?;
        let physical_source_edges = contracted_source.as_ref().map(|source| {
            graph
                .as_ref()
                .iter_edges()
                .filter_map(|(pair, edge_id, edge)| {
                    (pair.is_paired()
                        && !edge.data.is_dummy
                        && !internal_edges.contains(&edge_id)
                        && source
                            .reconstructible_outer_loop_coordinates(edge_id)
                            .is_some())
                    .then_some(edge_id)
                })
                .collect::<Vec<_>>()
        });
        let physical_boundary_edges = contracted_source.as_ref().map(|_| {
            let mut edges = graph
                .dummy_stripped_external_flows_of(&contract_subgraph)
                .included_iter()
                .map(|hedge| graph.underlying[&hedge])
                .collect::<Vec<_>>();
            edges.sort_unstable();
            edges.dedup();
            edges
        });

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
                let physical_source_energies = contracted_source
                    .as_ref()
                    .zip(physical_source_edges.as_ref())
                    .map(|(source, edges)| {
                        edges
                            .iter()
                            .map(|edge_id| {
                                let coordinates = source
                                    .reconstructible_outer_loop_coordinates(*edge_id)
                                    .expect("physical source edges were filtered above");
                                if coordinates.len() != reduced.orientation.loop_energy_map.len() {
                                    return Err(eyre!(
                                        "contracted edge {} has {} outer coordinates for {} generated loop-energy maps",
                                        usize::from(*edge_id),
                                        coordinates.len(),
                                        reduced.orientation.loop_energy_map.len(),
                                    ));
                                }
                                let mut energy = coordinates
                                    .iter()
                                    .zip(&reduced.orientation.loop_energy_map)
                                    .fold(LinearEnergyExpr::zero(), |sum, (coefficient, loop_energy)| {
                                        sum + loop_energy
                                            .clone()
                                            .scale_rational(coefficient.clone())
                                    });
                                let signature = &graph.loop_momentum_basis.edge_signatures[*edge_id];
                                for (external_edge, sign) in graph
                                    .loop_momentum_basis
                                    .ext_edges
                                    .iter()
                                    .zip(&signature.external)
                                {
                                    let coefficient = match sign {
                                        SignOrZero::Zero => 0,
                                        SignOrZero::Plus => 1,
                                        SignOrZero::Minus => -1,
                                    };
                                    energy = energy
                                        + LinearEnergyExpr::external(*external_edge, coefficient);
                                }
                                Ok((*edge_id, energy.canonical().to_atom_gs(&[])))
                            })
                            .collect::<Result<Vec<_>>>()
                    })
                    .transpose()?;
                let physical_source_frame = physical_source_energies
                    .as_deref()
                    .filter(|energies| !energies.is_empty())
                    .zip(physical_boundary_edges.as_deref());
                let reduced_expression = &reduced.expression * &production_prefactor;
                let localized = self.localized_orientation_terms(
                    graph,
                    &graph.loop_momentum_basis,
                    &reduced.orientation,
                    &reduced_expression,
                    &contract_subgraph,
                    &internal_edges,
                    valid_production_ids,
                    reduced.production_orientation_id,
                    source_edge_energy_map.as_deref(),
                    physical_source_frame,
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
