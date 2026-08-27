use std::{collections::BTreeMap, ops::Neg, sync::LazyLock};

use eyre::eyre;

use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, HedgePair, Orientation},
    subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function,
    id::Replacement,
    symbol,
};

use crate::{
    cff::{
        CutCFF, CutCFFIndex,
        expression::{OrientationExpression, OrientationID, OrientationSelector},
        orientations::GraphOrientation,
        surface::LinearEnergyExpr,
    },
    debug_tags,
    graph::{Graph, GraphThreeDSource, LMBext, LoopMomentumBasis, cuts::CutSet},
    settings::global::OrientationPattern,
    utils::{GS, W_},
    uv::{
        ApproximationType, DeferredIntegrands, Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike, OrientationProjection, UVCtx, integrated::IntegratedCts,
            local_4d::Full4dCts,
        },
        marker::{UvMarker, UvOperation},
        uv_graph::UVE,
    },
};
use color_eyre::Result;

pub(crate) struct FrozenActiveCt {
    pub active: OrientationIntegrands,
    pub frozen_integrands: Integrands,
}

impl FrozenActiveCt {
    pub(crate) fn combine(self) -> Result<OrientationIntegrands> {
        self.active.zip_mul_unmapped(&self.frozen_integrands)
    }
}

/// Residue integrands grouped by the selector and exact energy map that own
/// every numerator factor in the term. A missing source map uses the selected
/// production map; a present source map remains authoritative while production
/// IDs only partition theta sectors.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct OrientationIntegrandBranch {
    selector_id: OrientationID,
    source_edge_energy_map: Option<Vec<LinearEnergyExpr>>,
    integrands: Integrands,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct OrientationIntegrands(Vec<OrientationIntegrandBranch>);

impl From<Integrands> for OrientationIntegrands {
    fn from(integrands: Integrands) -> Self {
        Self(vec![OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: None,
            integrands,
        }])
    }
}

impl OrientationIntegrands {
    fn from_ids_and_indices(
        ids: impl IntoIterator<Item = OrientationID>,
        indices: &[CutCFFIndex],
    ) -> Self {
        Self(
            ids.into_iter()
                .map(|selector_id| OrientationIntegrandBranch {
                    selector_id,
                    source_edge_energy_map: None,
                    integrands: indices.iter().map(|index| (*index, Atom::Zero)).collect(),
                })
                .collect(),
        )
    }

    pub(crate) fn zip_add(&self, other: &Self) -> Result<Self> {
        // Deferred integrands expose nonzero evaluator calls through `iter`, so
        // this predicate only prunes compact branches proven to be zero.
        let is_zero = |branch: &OrientationIntegrandBranch| {
            branch.integrands.iter().all(|(_, atom)| atom.is_zero())
        };
        let fallback_zero = self
            .0
            .iter()
            .chain(&other.0)
            .find(|branch| is_zero(branch))
            .cloned();
        let mut branches = self
            .0
            .iter()
            .filter(|branch| !is_zero(branch))
            .cloned()
            .collect::<Vec<_>>();
        for right in other.0.iter().filter(|branch| !is_zero(branch)) {
            if let Some(left) = branches.iter_mut().find(|left| {
                left.selector_id == right.selector_id
                    && left.source_edge_energy_map == right.source_edge_energy_map
            }) {
                left.integrands = left.integrands.clone().zip_add(right.integrands.clone())?;
            } else {
                branches.push(right.clone());
            }
        }
        if branches.is_empty() {
            branches.extend(fallback_zero);
        }
        Ok(Self(branches))
    }

    fn zip_mul_unmapped(&self, other: &Integrands) -> Result<Self> {
        self.0
            .iter()
            .map(|branch| {
                Ok(OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.zip_mul(other)?,
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
    }

    pub(crate) fn map(&self, mut f: impl FnMut(&Atom) -> Atom) -> Self {
        Self(
            self.0
                .iter()
                .map(|branch| OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.map(&mut f),
                })
                .collect(),
        )
    }

    fn fallible_map(
        &self,
        mut f: impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        self.0
            .iter()
            .map(|branch| {
                Ok(OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.fallible_map(|atom| {
                        f(
                            branch.selector_id,
                            branch.source_edge_energy_map.as_deref(),
                            atom,
                        )
                    })?,
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
    }

    #[cfg(test)]
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.0.iter().flat_map(|branch| branch.integrands.iter())
    }

    pub(crate) fn iter_orientations(
        &self,
    ) -> impl Iterator<Item = (OrientationID, Option<&[LinearEnergyExpr]>, &Integrands)> {
        self.0.iter().map(|branch| {
            (
                branch.selector_id,
                branch.source_edge_energy_map.as_deref(),
                &branch.integrands,
            )
        })
    }
}

impl FromIterator<(OrientationID, Integrands)> for OrientationIntegrands {
    fn from_iter<T: IntoIterator<Item = (OrientationID, Integrands)>>(iter: T) -> Self {
        Self(
            iter.into_iter()
                .map(|(selector_id, integrands)| OrientationIntegrandBranch {
                    selector_id,
                    source_edge_energy_map: None,
                    integrands,
                })
                .collect(),
        )
    }
}

impl Neg for OrientationIntegrands {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(
            self.0
                .into_iter()
                .map(|branch| OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map,
                    integrands: -branch.integrands,
                })
                .collect(),
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Localizer<'a> {
    cutset: &'a CutSet,
    orientation: OrientationProjection<'a>,
}

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
    ) -> Result<(CutCFF, SuBitGraph)> {
        let contract_subgraph = self.cff_contract_subgraph(graph, to_contract);
        let options = self.orientation.cff_options(graph);
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

    fn exact_representatives(
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

    fn source_selector_representatives(
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
        // extend its surviving directions into full theta sectors.
        let mut explicit_reduced_orientation = reduced.data.orientation.clone();
        for edge in &contracted_edges {
            if explicit_reduced_orientation
                .iter()
                .any(|(explicit_edge, _)| explicit_edge == *edge)
            {
                explicit_reduced_orientation[*edge] = Orientation::Undirected;
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

    /// Assign a reduced CFF map to its production representation. Ordinary 3D
    /// evaluation keeps the reduced selector and one deterministic integrated-
    /// subgraph selector. Explicit-sum evaluation keeps each reduced residue
    /// once without orientation selectors; its ID is then map metadata only.
    #[allow(clippy::too_many_arguments)]
    fn localized_orientation_terms(
        self,
        graph: &Graph,
        reduced: &OrientationExpression,
        reduced_expression: &Atom,
        contract_subgraph: &SuBitGraph,
        internal_edges: &[EdgeIndex],
        production_orientation_id: Option<OrientationID>,
        source_edge_energy_map: Option<&[LinearEnergyExpr]>,
    ) -> Result<Vec<(OrientationID, Atom)>> {
        if let Some(id) = production_orientation_id
            && source_edge_energy_map.is_none()
        {
            // A stored root residue is already diagonal in its production map.
            // Its outer orientation key owns the numerator map, while its full
            // theta selector still partitions ordinary parametric evaluation.
            // Explicit-sum evaluators set all orientation inputs to zero, which
            // activates every theta factor and retains the complete sum.
            return Ok(self
                .orientation
                .orientation(id)
                .filter(|orientation| {
                    self.orientation
                        .orientation_pattern
                        .filter_orientation(orientation)
                })
                .map(|orientation| {
                    vec![(
                        id,
                        if self.orientation.explicit_orientation_sum_only {
                            reduced_expression.clone()
                        } else {
                            reduced_expression.clone() * orientation.orientation_thetas()
                        },
                    )]
                })
                .unwrap_or_default());
        }
        let representatives = if self.orientation.exact_orientations().is_some() {
            if source_edge_energy_map.is_some() {
                self.source_selector_representatives(graph, reduced, contract_subgraph)?
            } else {
                self.exact_representatives(graph, reduced, contract_subgraph)?
            }
        } else {
            self.coarse_representatives(&reduced.data.orientation)?
        };
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
        let Some(representative) = representatives.into_iter().max_by_key(|id| {
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
        }) else {
            return Ok(Vec::new());
        };
        if self.orientation.explicit_orientation_sum_only {
            return Ok(vec![(representative, reduced_expression.clone())]);
        }
        // The reduced term fixes every surviving physical direction. Any
        // remaining under-resolved direction is either contracted or belongs
        // to a generalized numerator-sampling contact map. Host the complete
        // branch in one resolved production orientation so ordinary runtime
        // orientation summation sees it exactly once; the branch-owned source
        // map still evaluates its numerator at the original sampling point.
        let selector = self
            .orientation
            .orientation(representative)
            .map(|orientation| orientation.orientation_thetas())
            .unwrap_or_else(Atom::one);
        Ok(vec![(
            representative,
            reduced_expression.clone() * selector,
        )])
    }

    #[cfg(test)]
    fn localized_orientation_term(
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

    fn projected_cff(
        self,
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        analysis_numerator: &Atom,
    ) -> Result<OrientationIntegrands> {
        let (cff, contract_subgraph) = self.cff(graph, to_contract, analysis_numerator)?;
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
                for (id, expression) in self.localized_orientation_terms(
                    graph,
                    &reduced.orientation,
                    &reduced_expression,
                    &contract_subgraph,
                    &internal_edges,
                    reduced.production_orientation_id,
                    source_edge_energy_map.as_deref(),
                )? {
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

    pub(crate) fn localize<S: ForestNodeLike>(
        &self,
        expr: &Atom,
        graph: &mut Graph,
        integrated_node: &S,
    ) -> Result<FrozenActiveCt> {
        if expr.is_zero() {
            let indices = self.cutset.residue_selector.generate_allowed_keys();
            return Ok(FrozenActiveCt {
                active: OrientationIntegrands::from_ids_and_indices(
                    self.orientation.orientation_ids(),
                    &indices,
                ),
                frozen_integrands: indices
                    .into_iter()
                    .map(|index| (index, Atom::one()))
                    .collect(),
            });
        }

        let to_contract = integrated_node.subgraph();
        let integrated_loop_count = graph.n_loops(to_contract);
        // Keep finite addbacks for nested multi-loop entries. Integrated expressions carry
        // their forest-composition signs where needed, so dropping the localized finite
        // representative here removes the Tint(T(...)) terms.
        let finite_ct = expr.clone();
        debug_tags!(#generation, #profile, #uv, #integrated, #local, #summary;
            stage = "localize_integrated_ct_forest_overlap",
            integrated_loop_count,
            "Applied integrated UV forest-overlap addback rule"
        );
        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        // CFF capacity belongs to the complete factorized numerator of this
        // reduced term. The integrated finite CT replaces the numerator of
        // the contracted spinney; the remaining graph and global factors are
        // still grown later by their existing owners.
        let reduced = graph
            .full_filter()
            .subtract(to_contract)
            .subtract(&graph.initial_state_cut);
        let outside_numerator = graph
            .numerator(&reduced, to_contract)
            .get_single_atom()
            .expect("Graph numerator should be available")
            * graph.global_atom();
        let analysis_numerator = expr * &outside_numerator;

        let localizing_integrand = GS.localizing_integrand(integrated_node.lmb());
        let active = self
            .projected_cff(graph, to_contract, &analysis_numerator)?
            .fallible_map(|orientation_id, source_edge_energy_map, localized| {
                let localized = localized * &fourddenoms;
                let localized_cff_byte_size = localized.as_view().get_byte_size();
                let mapped_finite_ct = self.orientation.map_numerator(
                    graph,
                    orientation_id,
                    source_edge_energy_map,
                    &finite_ct,
                )?;
                let active_ct = localized * mapped_finite_ct;
                let localized_ct = &active_ct * &localizing_integrand;
                debug_tags!(#generation, #profile, #uv, #integrated, #local, #term, #summary;
                    stage = "localize_integrated_ct_term",
                    localized_cff_byte_size,
                    active_ct_byte_size = active_ct.as_view().get_byte_size(),
                    localized_ct_byte_size = localized_ct.as_view().get_byte_size(),
                    "Integrated UV CT localization size checkpoint"
                );
                Ok(active_ct)
            })?;
        let frozen_integrands = self
            .cutset
            .residue_selector
            .generate_allowed_keys()
            .into_iter()
            .map(|index| (index, localizing_integrand.clone()))
            .collect();

        Ok(FrozenActiveCt {
            active,
            frozen_integrands,
        })
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
}

pub(crate) struct Local3DApproximation<'a> {
    localizer: Localizer<'a>,
    graph: &'a mut Graph,
    settings: &'a UVgenerationSettings,
}

impl<'a> Local3DApproximation<'a> {
    pub(crate) fn new(
        localizer: Localizer<'a>,
        graph: &'a mut Graph,
        settings: &'a UVgenerationSettings,
    ) -> Self {
        Self {
            localizer,
            graph,
            settings,
        }
    }

    pub(crate) fn run<S: ForestNodeLike, M: ForestNodeLike>(
        self,
        local: &Local3DCts,
        integrated: &IntegratedCts,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Local3DCts> {
        let integrated_t = self.localizer.localize(
            &integrated.physical_finite_counterterm_atom(),
            self.graph,
            given,
        )?;
        let ctx = UVCtx::new(self.graph, self.settings);
        let marker = UvMarker::new(ctx.settings);

        if let Some(active_sectors) = local.active_sectors() {
            let reduced_subgraph = current.reduced_subgraph(given);
            let mut next_sectors = Vec::with_capacity(active_sectors.len() + 1);

            for (active_subgraph, integrands) in active_sectors {
                let active_subgraph = active_subgraph.union(&reduced_subgraph);
                let integrands =
                    -integrands.fallible_map(Local3DLoopRescaling::FullSubgraph.map(
                        &ctx,
                        self.localizer.orientation,
                        current,
                        given,
                        Some(active_subgraph.clone()),
                    ))?;
                next_sectors.push((active_subgraph, integrands));
            }

            let integrated = -(integrated_t
                .active
                .fallible_map(Local3DLoopRescaling::ReducedSubgraph.map(
                    &ctx,
                    self.localizer.orientation,
                    current,
                    given,
                    Some(reduced_subgraph.clone()),
                ))?
                .zip_mul_unmapped(&integrated_t.frozen_integrands)?);
            next_sectors.push((reduced_subgraph, integrated));

            return Local3DCts::from_active_sectors(next_sectors)?.map(|atom| {
                Ok(marker.apply(
                    UvOperation::Approx,
                    marker_current.subgraph(),
                    marker_given.subgraph(),
                    atom,
                ))
            });
        }

        let local =
            -(local.map_orientations(full(&ctx, self.localizer.orientation, current, given))?);
        let integrated = -(integrated_t
            .active
            .fallible_map(reduced(&ctx, self.localizer.orientation, current, given))?
            .zip_mul_unmapped(&integrated_t.frozen_integrands)?);

        local.zip_add(&integrated)?.map(|atom| {
            Ok(marker.apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }

    pub(crate) fn run_local<S: ForestNodeLike, M: ForestNodeLike>(
        self,
        local: &Local3DCts,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Local3DCts> {
        let ctx = UVCtx::new(self.graph, self.settings);
        let reduced_subgraph = current.reduced_subgraph(given);
        let active_sectors = match local.active_sectors() {
            Some(active_sectors) => active_sectors
                .iter()
                .map(|(active_subgraph, integrands)| {
                    let active_subgraph = active_subgraph.union(&reduced_subgraph);
                    // Retain the full mask for descendants, but this replay step acts only
                    // on the part of that mask covered by its component-local path.
                    let rescaled_subgraph = active_subgraph.intersection(current.subgraph());
                    Ok((
                        active_subgraph,
                        -integrands.fallible_map(Local3DLoopRescaling::FullSubgraph.map(
                            &ctx,
                            self.localizer.orientation,
                            current,
                            given,
                            Some(rescaled_subgraph),
                        ))?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
            None => vec![(
                reduced_subgraph.clone(),
                -local
                    .integrands()
                    .fallible_map(Local3DLoopRescaling::FullSubgraph.map(
                        &ctx,
                        self.localizer.orientation,
                        current,
                        given,
                        Some(reduced_subgraph),
                    ))?,
            )],
        };
        Local3DCts::from_active_sectors(active_sectors)?.map(|atom| {
            Ok(UvMarker::new(ctx.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }

    pub(crate) fn run_integrated<S: ForestNodeLike, I: ForestNodeLike, M: ForestNodeLike>(
        self,
        integrated: &IntegratedCts,
        integrated_node: &I,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<Local3DCts> {
        let integrated = self.localizer.localize(
            &integrated.physical_finite_counterterm_atom(),
            self.graph,
            integrated_node,
        )?;
        let ctx = UVCtx::new(self.graph, self.settings);
        let active_subgraph = current.reduced_subgraph(given);
        let integrated = -(integrated
            .active
            .fallible_map(Local3DLoopRescaling::ReducedSubgraph.map(
                &ctx,
                self.localizer.orientation,
                current,
                given,
                Some(active_subgraph.clone()),
            ))?
            .zip_mul_unmapped(&integrated.frozen_integrands)?);

        Local3DCts::from_active_sectors(vec![(active_subgraph, integrated)])?.map(|atom| {
            Ok(UvMarker::new(ctx.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                atom,
            ))
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Local3DCts {
    representation: Local3DRepresentation,
    // A disconnected join can leave a different set of loop variables active
    // in each local/integrated cross term.
    active_sectors: Option<Vec<(SuBitGraph, OrientationIntegrands)>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Local3DRepresentation {
    Direct(OrientationIntegrands),
    Projected(Integrands),
}

impl From<Integrands> for Local3DCts {
    fn from(integrands: Integrands) -> Self {
        Self {
            representation: Local3DRepresentation::Direct(integrands.into()),
            active_sectors: None,
        }
    }
}

impl Neg for Local3DCts {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            representation: match self.representation {
                Local3DRepresentation::Direct(integrands) => {
                    Local3DRepresentation::Direct(integrands.map(|a| a.neg()))
                }
                Local3DRepresentation::Projected(integrands) => {
                    Local3DRepresentation::Projected(-integrands)
                }
            },
            active_sectors: self.active_sectors.map(|sectors| {
                sectors
                    .into_iter()
                    .map(|(active, integrands)| (active, -integrands))
                    .collect()
            }),
        }
    }
}

impl Local3DCts {
    pub(crate) fn zip_add(&self, other: &OrientationIntegrands) -> Result<Self> {
        if self.active_sectors.is_some() {
            return Err(eyre!(
                "an unlabelled local term cannot be added to active UV sectors"
            ));
        }
        let Local3DRepresentation::Direct(integrands) = &self.representation else {
            return Err(eyre!(
                "production-orientation terms cannot be added to projected 4D CFF terms"
            ));
        };
        Ok(Self {
            representation: Local3DRepresentation::Direct(integrands.zip_add(other)?),
            active_sectors: None,
        })
    }

    pub(crate) fn integrands(&self) -> &OrientationIntegrands {
        let Local3DRepresentation::Direct(integrands) = &self.representation else {
            panic!("projected 4D CFF terms have no production orientation IDs")
        };
        integrands
    }

    pub(crate) fn projected_integrands(&self) -> Result<&Integrands> {
        match &self.representation {
            Local3DRepresentation::Projected(integrands) => Ok(integrands),
            Local3DRepresentation::Direct(_) => Err(eyre!(
                "direct 3D terms are grouped by production orientation IDs"
            )),
        }
    }

    pub(crate) fn active_sectors(&self) -> Option<&[(SuBitGraph, OrientationIntegrands)]> {
        self.active_sectors.as_deref()
    }

    pub(crate) fn from_active_sectors<I: Into<OrientationIntegrands>>(
        active_sectors: Vec<(SuBitGraph, I)>,
    ) -> Result<Self> {
        let active_sectors = active_sectors
            .into_iter()
            .map(|(active, integrands)| (active, integrands.into()))
            .collect::<Vec<_>>();
        let mut sectors = active_sectors.iter();
        let mut integrands = sectors
            .next()
            .ok_or_else(|| eyre!("active UV counterterm sectors cannot be empty"))?
            .1
            .clone();
        for (_, sector) in sectors {
            integrands = integrands.zip_add(sector)?;
        }

        Ok(Self {
            representation: Local3DRepresentation::Direct(integrands),
            active_sectors: Some(active_sectors),
        })
    }

    pub(crate) fn map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        if let Some(active_sectors) = &self.active_sectors {
            let active_sectors = active_sectors
                .iter()
                .map(|(active, integrands)| {
                    Ok((
                        active.clone(),
                        integrands.fallible_map(|_, _, atom| f(atom))?,
                    ))
                })
                .collect::<Result<_>>()?;
            Self::from_active_sectors(active_sectors)
        } else {
            let representation = match &self.representation {
                Local3DRepresentation::Direct(integrands) => {
                    Local3DRepresentation::Direct(integrands.fallible_map(|_, _, atom| f(atom))?)
                }
                Local3DRepresentation::Projected(integrands) => {
                    Local3DRepresentation::Projected(integrands.fallible_map(f)?)
                }
            };
            Ok(Self {
                representation,
                active_sectors: None,
            })
        }
    }

    fn map_orientations(
        &self,
        mut f: impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        if let Some(active_sectors) = &self.active_sectors {
            let active_sectors = active_sectors
                .iter()
                .map(|(active, integrands)| Ok((active.clone(), integrands.fallible_map(&mut f)?)))
                .collect::<Result<_>>()?;
            Self::from_active_sectors(active_sectors)
        } else {
            let Local3DRepresentation::Direct(integrands) = &self.representation else {
                return Err(eyre!(
                    "projected 4D CFF terms cannot be remapped to production orientations"
                ));
            };
            Ok(Self {
                representation: Local3DRepresentation::Direct(integrands.fallible_map(f)?),
                active_sectors: None,
            })
        }
    }

    pub(crate) fn root(graph: &mut Graph, localizer: Localizer<'_>) -> Result<Self> {
        let analysis_numerator = graph.production_numerator_atom_for_full_3d_expression();
        let to_contract = graph.empty_subgraph::<SuBitGraph>();
        let cff = localizer.projected_cff(graph, &to_contract, &analysis_numerator)?;

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        Ok(Local3DCts {
            representation: Local3DRepresentation::Direct(cff.map(|atom| atom * &fourddenoms)),
            active_sectors: None,
        })
    }

    fn from_projected(integrands: Integrands) -> Self {
        Self {
            representation: Local3DRepresentation::Projected(integrands),
            active_sectors: None,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum Local3DLoopRescaling {
    FullSubgraph,
    ReducedSubgraph,
}

static OSE_FOR_LOCAL_3D_SERIES: LazyLock<Symbol> = LazyLock::new(|| {
    symbol!(
        "gammalooprs::OSE_for_local_3d_series",
        der = |_, arg, out| {
            if arg == 2 {
                **out = Atom::num(1);
            } else {
                **out = Atom::Zero;
            }
        }
    )
});

// #[debug_instrument(
//     current = %current.log_display(),
//     given = %given.log_display(),
//     reduced,
// )]
#[allow(clippy::too_many_arguments)]
pub(crate) fn t_tilde<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    orientation_id: OrientationID,
    source_edge_energy_map: Option<&[LinearEnergyExpr]>,
    current: &S,
    given: &S,
    cff: &Atom,
    active_subgraph: Option<&SuBitGraph>,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    let graph = ctx.graph;
    let settings = ctx.settings;
    let reduced = current.reduced_subgraph(given);
    let rescaled_subgraph = active_subgraph.unwrap_or_else(|| current.subgraph());
    let lmb_id = lmb
        .loop_edges
        .first()
        .copied()
        .unwrap_or_else(|| current.lmb_id());

    // split numerator momenta into OSEs and spatial parts
    let mut reps = Vec::new();
    for (p, eid, e) in graph.iter_edges_of(rescaled_subgraph) {
        if p.is_paired() {
            let e_mass = e.data.mass_atom();
            reps.push(GS.split_mom_pattern(eid, lmb_id, e_mass, settings.inner_products));
        }
    }

    let numerator = graph
        .numerator(&reduced, given.subgraph())
        .get_single_atom()
        .unwrap();
    let mut numerator = orientation
        .map_numerator(graph, orientation_id, source_edge_energy_map, &numerator)?
        .replace_multiple(&reps);

    // rescale the external momenta in the added numerator subgraph
    for e in &lmb.ext_edges {
        // println!("Rescale {}", e);
        numerator = numerator
            .replace(GS.emr_vec_index(*e, W_.x___))
            .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
    }

    let mut atomarg = cff * numerator;

    // add data for OSE computation and add an explicit sqrt
    for (p, ei, e) in graph.iter_edges_of(rescaled_subgraph) {
        let eid = usize::from(ei) as i64;
        if p.is_paired() {
            // set energies from inner_t on-shell
            atomarg = atomarg.replace(function!(GS.energy, eid)).with(GS.ose(ei));

            let e_mass = e.data.mass_atom();
            atomarg = atomarg.replace(GS.ose(ei)).with(GS.ose_full(
                ei,
                lmb_id,
                e_mass,
                None,
                settings.inner_products,
            ));
        }
    }

    atomarg = atomarg.replace_multiple(&reps);

    let mom_reps = graph.replacement_impl(
        |e, loops, externals| {
            Replacement::new(
                GS.emr_vec
                    .call_args([Atom::num(usize::from(e)), Atom::var(W_.x___)])
                    .to_pattern(),
                (loops
                    .replace(function!(GS.emr_vec, W_.x_))
                    .allow_new_wildcards_on_rhs(true)
                    .with(
                        FunctionBuilder::new(GS.emr_vec)
                            .add_arg(W_.x_)
                            .add_args([W_.x___])
                            .finish(),
                    )
                    + externals * GS.rescale)
                    .to_pattern(),
            )
        },
        &reduced,
        lmb,
        GS.emr_vec,
        GS.emr_vec,
        &[],
        &[W_.x___],
        HedgePair::is_paired,
        true,
    );

    atomarg = atomarg.replace_multiple(&mom_reps);
    atomarg = atomarg
        .replace(function!(GS.ose, W_.a___))
        .with(function!(*OSE_FOR_LOCAL_3D_SERIES, W_.a___));

    let a = atomarg.series(GS.rescale, Atom::Zero, -1).unwrap();

    let mut a = a
        .to_atom()
        .replace(function!(
            Symbol::DERIVATIVE,
            0,
            1,
            *OSE_FOR_LOCAL_3D_SERIES,
            W_.y___
        ))
        .with(Atom::num(1))
        .replace(function!(
            Symbol::DERIVATIVE,
            W_.x___,
            *OSE_FOR_LOCAL_3D_SERIES,
            W_.y___
        ))
        .with(Atom::num(0));
    a = a
        .replace(function!(*OSE_FOR_LOCAL_3D_SERIES, W_.a___))
        .with(function!(GS.ose, W_.a___));
    a = a.replace(GS.rescale).with(Atom::num(1));
    Ok(a)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn start<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    orientation_id: OrientationID,
    source_edge_energy_map: Option<&[LinearEnergyExpr]>,
    current: &S,
    given: &S,
    cff: &Atom,
    active_subgraph: Option<&SuBitGraph>,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    let graph = ctx.graph;
    let settings = ctx.settings;
    let reduced = current.reduced_subgraph(given);
    let rescaled_subgraph = active_subgraph.unwrap_or_else(|| current.subgraph());
    let lmb_id = lmb
        .loop_edges
        .first()
        .copied()
        .unwrap_or_else(|| current.lmb_id());
    let numerator = graph
        .numerator(&reduced, given.subgraph())
        .get_single_atom()
        .unwrap();
    let mut atomarg = cff
        * orientation.map_numerator(graph, orientation_id, source_edge_energy_map, &numerator)?;
    debug_tags!(#generation, #profile, #uv, #local, #trace;
        stage = "local_3d_start_initial",
        byte_size = atomarg.as_view().get_byte_size(),
        file.expr = %atomarg,
        "Local 3D start expression checkpoint"
    );
    // println!("CFF: {}", cff);

    // add data for OSE computation and add an explicit sqrt
    for (p, ei, e) in graph.iter_edges_of(rescaled_subgraph) {
        let eid = usize::from(ei) as i64;
        if p.is_paired() {
            // set energies from inner_t on-shell
            atomarg = atomarg.replace(function!(GS.energy, eid)).with(GS.ose(ei));

            let e_mass = e.data.mass_atom();
            atomarg = atomarg.replace(GS.ose(ei)).with(GS.ose_full(
                ei,
                lmb_id,
                e_mass,
                None,
                settings.inner_products,
            ));
        }
    }
    debug_tags!(#generation, #profile, #uv, #local, #trace;
        stage = "local_3d_start_after_ose_full",
        byte_size = atomarg.as_view().get_byte_size(),
        file.expr = %atomarg,
        "Local 3D start expression checkpoint"
    );

    // split numerator momenta into OSEs and spatial parts
    let mut reps = Vec::new();
    for (p, eid, e) in graph.iter_edges_of(rescaled_subgraph) {
        if p.is_paired() {
            let e_mass = e.data.mass_atom();
            let rep = GS.split_mom_pattern(eid, lmb_id, e_mass, settings.inner_products);
            debug_tags!(#uv, #local, #momentum, #trace;
                stage = "local_3d_start_split_mom_pattern",
                split_rep = %rep,
                "Local 3D start momentum split"
            );
            reps.push(rep);
        }
    }
    let atomarg = atomarg.replace_multiple(&reps);
    debug_tags!(#generation, #profile, #uv, #local, #trace;
        stage = "local_3d_start_output",
        byte_size = atomarg.as_view().get_byte_size(),
        file.expr = %atomarg,
        "Local 3D start expression checkpoint"
    );
    Ok(atomarg)
}

fn full<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    current: &S,
    given: &S,
) -> impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom> {
    Local3DLoopRescaling::FullSubgraph.map(ctx, orientation, current, given, None)
}
fn reduced<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    current: &S,
    given: &S,
) -> impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom> {
    Local3DLoopRescaling::ReducedSubgraph.map(ctx, orientation, current, given, None)
}
impl Local3DLoopRescaling {
    // #[debug_instrument(
    //     current = %current.log_display(),
    //     given = %given.log_display(),
    //     reduced,
    // )]
    pub(crate) fn t<S: super::ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        integrand: &Atom,
        active_subgraph: Option<&SuBitGraph>,
        lmb: &LoopMomentumBasis,
    ) -> Result<Atom> {
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);

        // only apply replacements for edges in the reduced graph
        let mom_reps = graph.uv_spatial_wrapped_replacement(&reduced, lmb, &[W_.x___]);
        for m in &mom_reps {
            debug_tags!(#uv,#momentum,#trace;mom_rep=%m,"Mom rep");
        }

        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_input",
            byte_size = integrand.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );
        let mut atomarg = integrand.replace_multiple(&mom_reps);
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_momentum_replacements",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        // Rescale every loop momentum still active in this sector, including
        // cycles expanded by earlier local operations.
        for e in &lmb.loop_edges {
            // println!("Rescale {}", e);
            atomarg = atomarg
                .replace(GS.emr_vec_index(*e, W_.x___))
                .with(GS.emr_vec_index(*e, W_.x___) * GS.rescale);
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_loop_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        // (re-)expand OSEs from the subgraph only
        for eid in lmb.loop_edges.iter() {
            let eid = eid.0 as i64;
            // rescale the whole OSE so that the function itself has no poles during the expansion
            atomarg = atomarg.replace(function!(GS.ose, eid, W_.prop_)).with(
                function!(
                    GS.ose,
                    eid,
                    (GS.m_uv_expansion * GS.m_uv_expansion * GS.rescale * GS.rescale + W_.prop_
                        - GS.m_uv_expansion * GS.m_uv_expansion)
                        / GS.rescale
                        / GS.rescale
                ) * GS.rescale
                    * GS.rescale,
            )
        }
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_ose_rescale",
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        atomarg = (atomarg * self.measure_scaling(ctx, current, given, active_subgraph))
            .replace(GS.rescale)
            .with(Atom::num(1) / GS.rescale);
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_before_series",
            loop_edges = ?lmb.loop_edges,
            byte_size = atomarg.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        debug_tags!(#uv, #local, #before_series; log.expr = atomarg, "Before series in t");

        let series = atomarg.series(GS.rescale, Atom::Zero, 0).unwrap();
        let series_atom = series.to_atom();
        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_after_series",
            byte_size = series_atom.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );

        debug_tags!(#uv, #local; expr = %series, "After series in t");
        let a = series_atom.replace(GS.rescale).with(Atom::num(1));

        debug_tags!(#generation, #profile, #uv, #local, #summary;
            stage = "local_3d_t_output",
            byte_size = a.as_view().get_byte_size(),
            "Local 3D T size checkpoint"
        );
        debug_tags!(#uv, #local; log.expr = a, "Local 3D approximation");
        Ok(a)
    }

    // #[instrument(skip(self, ctx, current, given))]
    fn map<S: super::ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        orientation: OrientationProjection<'_>,
        current: &S,
        given: &S,
        active_subgraph: Option<SuBitGraph>,
    ) -> impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom> {
        move |orientation_id, source_edge_energy_map, integrand| {
            let active_lmb = active_subgraph
                .as_ref()
                .filter(|active| *active != current.subgraph())
                .map(|active| {
                    ctx.graph.try_compatible_sub_lmb(
                        active,
                        ctx.graph.dummy_less_full_crown(active),
                        current.lmb(),
                    )
                })
                .transpose()?;
            let lmb = active_lmb.as_ref().unwrap_or_else(|| current.lmb());

            match current.renormalization_scheme() {
                ApproximationType::MUV | ApproximationType::PolePart => {
                    let started = start(
                        ctx,
                        orientation,
                        orientation_id,
                        source_edge_energy_map,
                        current,
                        given,
                        integrand,
                        active_subgraph.as_ref(),
                        lmb,
                    )?;
                    crate::debug_tags!(#generation, #profile, #uv, #local, #summary;
                        stage = "local_3d_kernel_after_start",
                        input_byte_size = integrand.as_view().get_byte_size(),
                        output_byte_size = started.as_view().get_byte_size(),
                        "Local 3D kernel size checkpoint"
                    );
                    self.t(ctx, current, given, &started, active_subgraph.as_ref(), lmb)
                }
                ApproximationType::IR => {
                    let t_tilde = t_tilde(
                        ctx,
                        orientation,
                        orientation_id,
                        source_edge_energy_map,
                        current,
                        given,
                        integrand,
                        active_subgraph.as_ref(),
                        lmb,
                    )?;
                    Ok(self.t(
                        ctx,
                        current,
                        given,
                        &start(
                            ctx,
                            orientation,
                            orientation_id,
                            source_edge_energy_map,
                            current,
                            given,
                            integrand,
                            active_subgraph.as_ref(),
                            lmb,
                        )?,
                        active_subgraph.as_ref(),
                        lmb,
                    )? + &t_tilde
                        - self.t(ctx, current, given, &t_tilde, active_subgraph.as_ref(), lmb)?)
                }
                ApproximationType::VaccuumLimit => Err(eyre!("Not yet implemented VaccuumLimit")),
                ApproximationType::OS => Err(eyre!("Not yet implemented OS")),
                ApproximationType::Unsubtracted => {
                    panic!("should have been kept out of the wood");
                }
            }
        }
    }

    fn measure_scaling<S: super::ForestNodeLike>(
        self,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        active_subgraph: Option<&SuBitGraph>,
    ) -> Atom {
        let n_rescaled_loops = match active_subgraph {
            Some(active_subgraph) => ctx.graph.n_loops(active_subgraph),
            None => match self {
                Local3DLoopRescaling::FullSubgraph => ctx.graph.n_loops(current.subgraph()),
                Local3DLoopRescaling::ReducedSubgraph => {
                    ctx.graph.n_loops(current.subgraph()) - ctx.graph.n_loops(given.subgraph())
                }
            },
        };

        Atom::var(GS.rescale).pow(3 * n_rescaled_loops as i64)
    }
}

#[cfg(test)]
mod tests {
    use super::{OrientationIntegrandBranch, OrientationIntegrands};
    use crate::{
        cff::{
            CutCFFIndex,
            expression::{
                GammaLoopOrientationExpression, OrientationData, OrientationExpression,
                OrientationID,
            },
            orientations::GraphOrientation,
            surface::LinearEnergyExpr,
        },
        dot,
        graph::{FeynmanGraph, FourDDenominator, Graph, cuts::CutSet, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        utils::GS,
        uv::{
            ApproximationType, Spinney, UltravioletGraph,
            approx::{
                OrientationProjection, Rooted,
                local_3d::Localizer,
                local_4d::{Full4dCts, Local4dCts},
            },
            hedge_poset::OwnedForestNode,
        },
    };
    use color_eyre::Result;
    use eyre::eyre;
    use linnet::half_edge::{
        involution::{EdgeIndex, EdgeVec, Orientation},
        subgraph::{InternalSubGraph, SubSetOps},
    };
    use std::{collections::BTreeMap, sync::OnceLock};
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, FunctionBuilder},
        function,
    };
    use typed_index_collections::TiVec;

    static TWO_EDGE_GRAPH: OnceLock<Graph> = OnceLock::new();

    fn orientation(value: i8) -> Orientation {
        match value {
            1 => Orientation::Default,
            -1 => Orientation::Reversed,
            0 => Orientation::Undirected,
            _ => panic!("invalid orientation encoding"),
        }
    }

    fn edgevec(values: impl IntoIterator<Item = i8>) -> EdgeVec<Orientation> {
        EdgeVec::from_iter(values.into_iter().map(orientation))
    }

    fn edges(values: impl IntoIterator<Item = usize>) -> Vec<EdgeIndex> {
        values.into_iter().map(EdgeIndex).collect()
    }

    fn energy_map(
        orientation: EdgeVec<Orientation>,
        edge_energy_map: Vec<LinearEnergyExpr>,
    ) -> OrientationExpression {
        OrientationExpression {
            data: OrientationData::new(orientation),
            loop_energy_map: Vec::new(),
            edge_energy_map,
            variants: Vec::new(),
        }
    }

    fn two_edge_graph() -> Result<Graph> {
        Ok(TWO_EDGE_GRAPH
            .get_or_init(|| {
                test_initialise().expect("test model initialization succeeds");
                dot!(
                    digraph G {
                        edge [particle="scalar_1"];
                        node [num=1];
                        a -> b [id=0];
                        a -> b [id=1];
                    },
                    "scalars"
                )
                .expect("two-edge projection graph parses")
            })
            .clone())
    }

    #[test]
    fn production_emr_map_cancels_one_powered_denominator() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let edge = EdgeIndex(0);
        let on_shell_energy_squared = (1..=3).fold(
            graph.underlying[edge].particle.mass_atom().pow(2),
            |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
            },
        );
        let numerator = GS.emr_mom(edge, GS.cind(0)).pow(2) - &on_shell_energy_squared;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&numerator),
        )?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(&production, &options, &pattern, false),
        );
        let contract = graph.empty_subgraph();
        let powered = localizer.projected_cff(&mut graph, &contract, &numerator)?;
        let mut explicit_sum = Atom::Zero;
        for (selector_id, source_map, integrands) in powered.iter_orientations() {
            let mapped = localizer.map_numerator(&graph, selector_id, source_map, &numerator)?;
            explicit_sum += integrands
                .iter()
                .fold(Atom::Zero, |sum, (_, term)| sum + term * &mapped);
        }
        // The two parallel physical edges carry opposite routing signs. Their
        // acyclic production sectors are therefore (+,-) and (-,+); (+,+) is
        // not a production orientation. The under-resolved cancellation
        // contact is deterministically hosted once in the former sector.
        let mut selected_production = edgevec([1, -1]).select(&explicit_sum);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
        let denominators = [
            FourDDenominator {
                source_edge: edge,
                momentum: momentum.clone(),
                mass_squared: mass_squared.clone(),
                full_expr: Atom::var(symbolica::symbol!("local_3d_test::first")),
            },
            FourDDenominator {
                source_edge: EdgeIndex(1),
                momentum: -momentum,
                mass_squared: graph.underlying[EdgeIndex(1)].particle.mass_atom().pow(2),
                full_expr: Atom::var(symbolica::symbol!("local_3d_test::second")),
            },
        ];
        let (exact_powered, _) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        // A projected post-4D source owns its complete orientation sum. Once
        // the numerator cancels an occurrence, that provenance slot can be
        // undirected and must not be used to filter the surviving contact.
        let mut exact_powered_sum = exact_powered
            .terms
            .values()
            .flat_map(|term| {
                term.orientations.iter().map(|orientation| {
                    Ok::<_, eyre::Error>(
                        &orientation.expression
                            * term.map_exact_source_numerator(&orientation.orientation)?,
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .sum::<Atom>()
            * Atom::num(exact_powered.production_prefactor_factor());
        let (exact_lower, _) =
            graph.cff_from_4d_denominators(&denominators[1..], &cutset, &options, &Atom::one())?;
        // The exact one-denominator source enumerates both contour directions.
        // Compare the remaining-edge direction selected by the (+,-)
        // production sector. Denominator-routing signs and physical selector
        // signs are separate conventions.
        let remaining_lower_occurrence = EdgeIndex(graph.underlying.n_edges());
        let mut exact_lower_sum = exact_lower
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .find(|orientation| {
                orientation.orientation.data.orientation[remaining_lower_occurrence]
                    == Orientation::Default
            })
            .expect("the lower source has the selected remaining-edge direction")
            .expression
            .clone()
            * Atom::num(exact_lower.production_prefactor_factor());
        let mass = graph.underlying[edge].particle.mass_atom();
        for expression in [
            &mut selected_production,
            &mut exact_powered_sum,
            &mut exact_lower_sum,
        ] {
            *expression = expression.replace(mass.clone()).with(Atom::one());
            for spatial_index in 1..=3 {
                *expression = expression
                    .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                    .with(Atom::Zero);
            }
            *expression = expression
                .replace(GS.ose(EdgeIndex(0)))
                .with(Atom::one())
                .replace(GS.ose(EdgeIndex(1)))
                .with(Atom::one());
        }
        let exact_difference = (&exact_powered_sum - &exact_lower_sum).together();
        assert!(
            exact_difference.is_zero(),
            "summed exact-source powered and selected lower channels differ: powered={exact_powered_sum}, lower={exact_lower_sum}, difference={exact_difference}"
        );
        let production_difference = (&selected_production - &exact_lower_sum).together();
        assert!(
            production_difference.is_zero(),
            "production EMR mapping does not cancel the powered denominator: production={selected_production}, lower={exact_lower_sum}, difference={production_difference}"
        );
        Ok(())
    }

    #[test]
    fn orientation_term_keeps_external_selectors_and_adds_internal_ones() {
        let reduced_expression = function!(GS.ose, 0);
        let reduced_orientation = edgevec([1, 0, -1]);
        let valid = vec![edgevec([1, 1, -1]), edgevec([1, -1, -1])];
        let pat = OrientationPattern::default();
        let cutset = CutSet::empty(1);
        let loc = Localizer::new(&cutset, OrientationProjection::new(&valid, &pat));
        let localized = loc
            .localized_orientation_term(&reduced_expression, &reduced_orientation, &edges([1]))
            .unwrap();

        let expected = reduced_expression
            * GS.sign_theta(GS.sign(EdgeIndex(0)))
            * GS.sign_theta(-GS.sign(EdgeIndex(2)))
            * GS.sign_theta(GS.sign(EdgeIndex(1)));
        assert_eq!(localized, expected);
    }

    #[test]
    fn exact_projection_distinguishes_affine_maps_with_the_same_coarse_orientation() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 2),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = energy_map(
            edgevec([1, 0]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::zero(),
            ],
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );

        assert_eq!(
            localizer.exact_representatives(
                &graph,
                &reduced,
                &graph.get_edge_subgraph(EdgeIndex(1)),
            )?,
            vec![OrientationID(0)]
        );
        Ok(())
    }

    #[test]
    fn exact_projection_errors_when_no_affine_map_matches() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        )]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = energy_map(
            edgevec([1, 0]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 2),
                LinearEnergyExpr::zero(),
            ],
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );

        let error = localizer
            .exact_representatives(&graph, &reduced, &graph.get_edge_subgraph(EdgeIndex(1)))
            .expect_err("a different surviving affine map must not use a coarse fallback");
        assert!(
            error
                .to_string()
                .contains("no production energy map exactly extends")
        );
        Ok(())
    }

    #[test]
    fn source_energy_map_owns_factorized_contact_branch() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
            energy_map(
                edgevec([1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = energy_map(
            edgevec([1, 0]),
            vec![LinearEnergyExpr::uniform_scale(2), LinearEnergyExpr::zero()],
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let contract = graph.get_edge_subgraph(EdgeIndex(1));

        let strict_error = localizer
            .exact_representatives(&graph, &reduced, &contract)
            .expect_err("the source map must not weaken strict production-map matching");
        assert!(
            strict_error
                .to_string()
                .contains("no production energy map exactly extends")
        );
        let localized = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            Some(&reduced.edge_energy_map),
        )?;
        let outer_selector = GS.sign_theta(GS.sign(EdgeIndex(0)));
        assert_eq!(localized.len(), 1);
        assert_eq!(localized[0].0, OrientationID(0));
        assert_eq!(
            localized[0].1,
            &outer_selector * GS.sign_theta(GS.sign(EdgeIndex(1)))
        );

        let factor_a = Atom::var(symbolica::symbol!("source_factor_a"));
        let factor_b = Atom::var(symbolica::symbol!("source_factor_b"));
        let edge_energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let numerator = (edge_energy.clone() + &factor_a) * (edge_energy + &factor_b);
        let mapped = localizer.map_numerator(
            &graph,
            OrientationID(0),
            Some(&reduced.edge_energy_map),
            &numerator,
        )?;
        assert_eq!(
            localizer.map_numerator(
                &graph,
                OrientationID(1),
                Some(&reduced.edge_energy_map),
                &numerator,
            )?,
            mapped,
            "production IDs partition theta sectors but do not remap source-owned numerators"
        );
        let sampled_energy = Atom::num(2) * Atom::var(GS.numerator_sampling_scale);
        assert_eq!(
            mapped,
            (sampled_energy.clone() + &factor_a) * (sampled_energy + &factor_b)
        );
        let AtomView::Mul(factors) = mapped.as_view() else {
            panic!("the source-mapped numerator must remain factorized")
        };
        assert_eq!(
            factors
                .iter()
                .filter(|factor| matches!(factor, AtomView::Add(_)))
                .count(),
            2
        );
        assert_ne!(
            localizer.map_numerator(&graph, OrientationID(0), None, &numerator)?,
            mapped
        );
        let coarse = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
        assert!(
            coarse
                .map_numerator(
                    &graph,
                    OrientationID(0),
                    Some(&reduced.edge_energy_map),
                    &numerator,
                )
                .expect_err("a coarse projector cannot consume an exact source map")
                .to_string()
                .contains("cannot be carried by a coarse orientation projector")
        );

        let index = CutCFFIndex::new_all_none();
        let source = OrientationIntegrands(vec![OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: Some(reduced.edge_energy_map.clone()),
            integrands: [(index, Atom::num(3))].into_iter().collect(),
        }]);
        let distinct_map = vec![LinearEnergyExpr::zero(), LinearEnergyExpr::zero()];
        let distinct = OrientationIntegrands(vec![OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: Some(distinct_map.clone()),
            integrands: [(index, Atom::num(5))].into_iter().collect(),
        }]);
        let production_zero = OrientationIntegrands(vec![OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: None,
            integrands: [(index, Atom::Zero)].into_iter().collect(),
        }]);
        let union = source.zip_add(&distinct)?.zip_add(&production_zero)?;
        assert_eq!(union.iter_orientations().count(), 2);
        assert!(
            union
                .iter_orientations()
                .all(|(_, source_map, _)| source_map.is_some()),
            "a zero production-map identity branch is pruned beside source-owned terms"
        );
        assert!(
            union
                .iter_orientations()
                .any(|(_, map, _)| map == Some(reduced.edge_energy_map.as_slice()))
        );
        assert!(
            union
                .iter_orientations()
                .any(|(_, map, _)| map == Some(distinct_map.as_slice()))
        );

        let same_source = OrientationIntegrands(vec![OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: Some(reduced.edge_energy_map.clone()),
            integrands: [(index, Atom::num(7))].into_iter().collect(),
        }]);
        let merged = source.zip_add(&same_source)?;
        assert_eq!(merged.iter_orientations().count(), 1);
        assert_eq!(
            merged.iter().next(),
            Some((&index, &Atom::num(10))),
            "only identical selector/source-map branches are coalesced"
        );
        Ok(())
    }

    #[test]
    fn exact_projection_skips_extensions_excluded_by_the_full_pattern() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
            energy_map(
                edgevec([-1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), -1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = energy_map(
            edgevec([1, 0]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::zero(),
            ],
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::from_user_pattern("(-1,-1)")?;
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let contract = graph.get_edge_subgraph(EdgeIndex(1));

        assert!(
            localizer
                .exact_representatives(&graph, &reduced, &contract)?
                .is_empty(),
            "an exact extension excluded by the full pattern is a zero contribution, not a missing-map error"
        );
        assert!(
            localizer
                .localized_orientation_terms(
                    &graph,
                    &reduced,
                    &Atom::one(),
                    &contract,
                    &[EdgeIndex(1)],
                    None,
                    None,
                )?
                .is_empty(),
            "an excluded extension contributes zero without normalization"
        );
        Ok(())
    }

    #[test]
    fn exact_cff_defers_contracted_edge_patterns_to_full_projection() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
            energy_map(
                edgevec([-1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), -1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::from_user_pattern("(1,-1)")?;
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let to_contract = graph.get_edge_subgraph(EdgeIndex(1));
        let analysis_numerator = graph.production_numerator_atom_for_full_3d_expression();
        let projected = localizer.projected_cff(&mut graph, &to_contract, &analysis_numerator)?;

        assert_eq!(
            projected
                .iter_orientations()
                .map(|(id, _, _)| id)
                .collect::<Vec<_>>(),
            vec![OrientationID(0)]
        );
        assert!(projected.iter().any(|(_, atom)| !atom.is_zero()));
        Ok(())
    }

    #[test]
    fn coarse_cff_keeps_denominator_only_capacity() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
        let to_contract = graph.empty_subgraph();
        let unsupported_numerator = (GS.emr_mom(EdgeIndex(0), GS.cind(0)) + Atom::one()).pow(-1);

        let projected =
            localizer.projected_cff(&mut graph, &to_contract, &unsupported_numerator)?;

        assert!(projected.iter().any(|(_, atom)| !atom.is_zero()));
        Ok(())
    }

    #[test]
    fn exact_localization_capacity_excludes_replaced_spinney_numerator() -> Result<()> {
        // Share model initialization with the other graph-backed tests in this module.
        let _ = two_edge_graph()?;
        let mut graph: Graph = dot!(
            digraph replaced_spinney_numerator {
                edge [pdg=1000 num=1 mass=1]
                node [num=1]
                a -> b [id=0 lmb_id=0 num="1/(Q(0,spenso::cind(0))+1)"]
                a -> b [id=1]
            },
            "scalars"
        )?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let contracted = graph.tree_edges.subtract(&graph.initial_state_cut);
        let contract_edges = graph.paired_edges(&contracted);
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph
            .generate_3d_expression_for_integrand(&contract_edges, &canonization, &options, None)?
            .expression
            .orientations;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let subgraph =
            InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph.as_ref());
        let integrated_node = OwnedForestNode {
            spinney: Spinney::with_scheme(
                subgraph,
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                0,
            )
            .expect("the complete bubble has a compatible loop basis"),
            topo_order: 0,
        };

        let localized = localizer.localize(&Atom::one(), &mut graph, &integrated_node)?;

        assert!(localized.active.iter().any(|(_, atom)| !atom.is_zero()));
        Ok(())
    }

    #[test]
    fn exact_localization_maps_finite_ct_instead_of_leaving_it_unmapped() -> Result<()> {
        // Share model initialization with the other graph-backed tests in this module.
        let _ = two_edge_graph()?;
        let mut graph: Graph = dot!(
            digraph finite_ct_map {
                edge [particle="scalar_1"];
                node [num=1];
                external [style=invis];
                external -> A:0 [id=3];
                C:1 -> external [id=4];
                A -> B [id=0];
                A -> B [id=1];
                B -> C [id=2];
            },
            "scalars"
        )?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let contracted = graph.tree_edges.subtract(&graph.initial_state_cut);
        let contract_edges = graph.paired_edges(&contracted);
        assert!(contract_edges.contains(&EdgeIndex(2)));
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph
            .generate_3d_expression_for_integrand(&contract_edges, &canonization, &options, None)?
            .expression
            .orientations;
        let (orientation_id, finite_ct, mapped_finite_ct) = production
            .iter_enumerated()
            .flat_map(|(orientation_id, orientation)| {
                orientation
                    .edge_energy_map
                    .iter()
                    .enumerate()
                    .filter(|(edge, energy)| {
                        contract_edges.contains(&EdgeIndex(*edge))
                            && !energy.external_terms.is_empty()
                    })
                    .map(move |(edge, _)| (orientation_id, orientation, EdgeIndex(edge)))
            })
            .find_map(|(orientation_id, orientation, edge)| {
                let finite_ct = GS.emr_mom(edge, GS.cind(0));
                let replacements = orientation.energy_replacements_gs(&graph);
                let mapped_finite_ct = finite_ct.replace_multiple(&replacements);
                (mapped_finite_ct != finite_ct).then_some((
                    orientation_id,
                    finite_ct,
                    mapped_finite_ct,
                ))
            })
            .ok_or_else(|| {
                eyre!(
                    "the finite-CT fixture must contain a contracted affine map that changes an edge energy"
                )
            })?;
        let pattern =
            OrientationPattern::from_orientation(&production[orientation_id].data.orientation);
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let integrated_node = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };

        let baseline = localizer.localize(&Atom::one(), &mut graph, &integrated_node)?;
        let localized = localizer.localize(&finite_ct, &mut graph, &integrated_node)?;
        let baseline = baseline
            .active
            .iter_orientations()
            .find_map(|(id, _, integrands)| (id == orientation_id).then_some(integrands))
            .expect("the selected production map has a baseline branch");
        let localized = localized
            .active
            .iter_orientations()
            .find_map(|(id, _, integrands)| (id == orientation_id).then_some(integrands))
            .expect("the selected production map has a finite-CT branch");

        assert_eq!(localized, &baseline.map(|atom| atom * &mapped_finite_ct));
        assert_ne!(localized, &baseline.map(|atom| atom * &finite_ct));
        Ok(())
    }

    #[test]
    fn contracted_exact_extensions_follow_evaluator_orientation_mode() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
            energy_map(
                edgevec([1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = energy_map(
            edgevec([1, 0]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::zero(),
            ],
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let contract = graph.get_edge_subgraph(EdgeIndex(1));
        let localized = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            None,
        )?;
        let expected_default =
            GS.sign_theta(GS.sign(EdgeIndex(0))) * GS.sign_theta(GS.sign(EdgeIndex(1)));

        assert_eq!(localized, vec![(OrientationID(0), expected_default)]);

        let explicit_localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let explicit = explicit_localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            None,
        )?;
        assert_eq!(explicit, vec![(OrientationID(0), Atom::one())]);
        Ok(())
    }

    #[test]
    fn stored_root_residue_keeps_its_production_orientation_diagonal() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
            energy_map(
                edgevec([1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = energy_map(
            edgevec([1, 0]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::zero(),
            ],
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let contract = graph.get_edge_subgraph(EdgeIndex(1));

        let localized = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            Some(OrientationID(1)),
            None,
        )?;

        assert_eq!(localized.len(), 1);
        assert_eq!(localized[0].0, OrientationID(1));
        assert_eq!(
            localized[0].1,
            GS.sign_theta(GS.sign(EdgeIndex(0))) * GS.sign_theta(-GS.sign(EdgeIndex(1)))
        );
        assert_eq!(
            production[OrientationID(0)]
                .data
                .orientation
                .select(localized[0].1.as_view()),
            Atom::Zero,
            "a stored root residue must not leak into another runtime orientation"
        );
        assert_eq!(
            production[OrientationID(1)]
                .data
                .orientation
                .select(localized[0].1.as_view()),
            Atom::one()
        );

        let explicit_localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let explicit = explicit_localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            Some(OrientationID(1)),
            None,
        )?;
        assert_eq!(explicit, vec![(OrientationID(1), Atom::one())]);
        Ok(())
    }

    #[test]
    fn stored_generalized_root_map_uses_one_resolved_orientation_selector() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
            energy_map(
                edgevec([1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
            energy_map(
                edgevec([1, 0]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::zero(),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = production[OrientationID(2)].clone();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let contract = graph.empty_subgraph();

        let localized = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[],
            Some(OrientationID(2)),
            Some(&reduced.edge_energy_map),
        )?;

        assert_eq!(
            localized,
            vec![(
                OrientationID(0),
                GS.sign_theta(GS.sign(EdgeIndex(0))) * GS.sign_theta(GS.sign(EdgeIndex(1))),
            )],
            "an under-resolved numerator sampling map must live in one resolved physical orientation",
        );

        let explicit_localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        assert_eq!(
            explicit_localizer.localized_orientation_terms(
                &graph,
                &reduced,
                &Atom::one(),
                &contract,
                &[],
                Some(OrientationID(2)),
                Some(&reduced.edge_energy_map),
            )?,
            vec![(OrientationID(0), Atom::one())],
            "an explicit orientation sum keeps the generalized branch exactly once without a selector",
        );
        Ok(())
    }

    #[test]
    fn exact_typed_4d_projection_keeps_source_local_explicit_orientation_sum() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
            energy_map(
                edgevec([-1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), -1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::from_user_pattern("(1,-1)")?;
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let cograph = graph.full_filter();
        let source = Full4dCts::from_coefficient(&Atom::one(), &graph, &cograph);
        let root_spinney = graph.empty_subgraph();
        let projected = localizer.project_4d(&source, &mut graph, &root_spinney)?;
        let projected = projected.projected_integrands()?;
        let root = CutCFFIndex::new_all_none();
        let source_local_terms = projected
            .deferred_terms(&root)
            .expect("the projected root keeps its source-local CFF bodies");
        let opposite_pattern = OrientationPattern::from_user_pattern("(-1,1)")?;
        let opposite = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &opposite_pattern, false),
        )
        .project_4d(&source, &mut graph, &root_spinney)?;

        assert!(
            source_local_terms.iter().any(|atom| !atom.is_zero()),
            "the exact source-local orientation sum must retain a nonzero body"
        );
        assert_eq!(
            opposite.projected_integrands()?,
            projected,
            "production orientation patterns must not filter an exact source-local orientation sum"
        );
        assert!(
            projected
                .materialize()
                .iter()
                .any(|(_, atom)| !atom.is_zero())
        );
        Ok(())
    }

    #[test]
    fn typed_4d_projection_localizes_each_frozen_loop_sector_once() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
        let cograph = graph.full_filter();
        let root_spinney = graph.empty_subgraph();
        let frozen_lmb = graph.loop_momentum_basis.clone();
        let localizing_integrand = GS.localizing_integrand(&frozen_lmb);
        assert!(!localizing_integrand.is_one());

        let baseline_source = Full4dCts::from_coefficient(&Atom::one(), &graph, &cograph);
        let frozen_source = Full4dCts::from_frozen_coefficient(
            &Atom::var(GS.integrated_loop_scale).pow(4),
            &graph,
            &cograph,
            &frozen_lmb,
        );
        let baseline = localizer.project_4d(&baseline_source, &mut graph, &root_spinney)?;
        let frozen = localizer.project_4d(&frozen_source, &mut graph, &root_spinney)?;
        let baseline = baseline.projected_integrands()?.materialize();
        let frozen = frozen.projected_integrands()?.materialize();
        let expected = baseline.map(|atom| atom * &localizing_integrand);

        assert_eq!(frozen, expected);
        assert!(frozen.iter().all(|(_, atom)| {
            !atom.contains_symbol(GS.integrated_loop_scale)
                && atom.contains_symbol(GS.localizing_integrand)
        }));
        Ok(())
    }

    #[test]
    fn typed_4d_projection_maps_finite_ct_instead_of_leaving_it_unmapped() -> Result<()> {
        // Share model initialization with the other graph-backed tests in this module.
        let _ = two_edge_graph()?;
        let mut graph: Graph = dot!(
            digraph finite_ct_map {
                edge [particle="scalar_1"];
                node [num=1];
                external [style=invis];
                external -> A:0 [id=3];
                C:1 -> external [id=4];
                A -> B [id=0];
                A -> B [id=1];
                B -> C [id=2];
            },
            "scalars"
        )?;
        let options = graph.denominator_only_cff_3d_expression_options();
        let contracted = graph.tree_edges.subtract(&graph.initial_state_cut);
        let contract_edges = graph.paired_edges(&contracted);
        assert!(contract_edges.contains(&EdgeIndex(2)));
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph
            .generate_3d_expression_for_integrand(&contract_edges, &canonization, &options, None)?
            .expression
            .orientations;
        let (orientation_id, finite_ct) = production
            .iter_enumerated()
            .flat_map(|(orientation_id, orientation)| {
                orientation
                    .edge_energy_map
                    .iter()
                    .enumerate()
                    .filter(|(edge, energy)| {
                        contract_edges.contains(&EdgeIndex(*edge))
                            && !energy.external_terms.is_empty()
                    })
                    .map(move |(edge, _)| (orientation_id, orientation, EdgeIndex(edge)))
            })
            .find_map(|(orientation_id, orientation, edge)| {
                let finite_ct = GS.emr_mom(edge, GS.cind(0));
                let replacements = orientation.energy_replacements_gs(&graph);
                let mapped_finite_ct = finite_ct.replace_multiple(&replacements);
                (mapped_finite_ct != finite_ct).then_some((orientation_id, finite_ct))
            })
            .ok_or_else(|| {
                eyre!(
                    "the finite-CT fixture must contain a contracted affine map that changes an edge energy"
                )
            })?;
        let pattern =
            OrientationPattern::from_orientation(&production[orientation_id].data.orientation);
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let cograph = graph.full_filter().subtract(&graph.initial_state_cut);
        let baseline_source = Full4dCts::from_coefficient(&Atom::one(), &graph, &cograph);
        let localized_source = Full4dCts::from_coefficient(&finite_ct, &graph, &cograph);
        let root_spinney = graph.empty_subgraph();

        let mut baseline_terms = baseline_source.terms()?;
        assert_eq!(baseline_terms.len(), 1);
        let mut active_denominators = Vec::new();
        for denominator in baseline_terms
            .pop()
            .expect("the baseline has one exact 4D term")
            .denominators
        {
            if denominator.depends_on_loop(&graph, false)? {
                active_denominators.push(denominator);
            }
        }
        let (source_cff, _) =
            graph.cff_from_4d_denominators(&active_denominators, &cutset, &options, &finite_ct)?;
        let mapped_finite_cts = source_cff
            .terms
            .iter()
            .map(|(index, term)| {
                let mapped = term
                    .orientations
                    .iter()
                    .map(|orientation| term.map_exact_source_numerator(&orientation.orientation))
                    .collect::<Result<Vec<_>>>()?;
                Ok((*index, mapped))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;

        let baseline = localizer.project_4d(&baseline_source, &mut graph, &root_spinney)?;
        let localized = localizer.project_4d(&localized_source, &mut graph, &root_spinney)?;
        let baseline = baseline.projected_integrands()?;
        let localized = localized.projected_integrands()?;
        let mut mapped_any = false;
        for (index, mapped_finite_cts) in mapped_finite_cts {
            let baseline_terms = baseline
                .deferred_terms(&index)
                .expect("the exact source has a baseline branch");
            let localized_terms = localized
                .deferred_terms(&index)
                .expect("the exact source has a finite-CT branch");
            assert_eq!(baseline_terms.len(), mapped_finite_cts.len());
            assert_eq!(localized_terms.len(), mapped_finite_cts.len());
            for ((baseline, localized), mapped_finite_ct) in baseline_terms
                .iter()
                .zip(localized_terms)
                .zip(mapped_finite_cts)
            {
                assert_eq!(localized, &(baseline * &mapped_finite_ct));
                if mapped_finite_ct != finite_ct {
                    mapped_any = true;
                    assert_ne!(localized, &(baseline * &finite_ct));
                }
            }
        }
        assert!(
            mapped_any,
            "the source-local finite-CT fixture must exercise a nontrivial affine energy map"
        );
        Ok(())
    }

    #[test]
    fn typed_4d_projection_analyzes_each_term_owned_numerator() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
        let cograph = graph.full_filter();
        let non_polynomial_term = (GS.emr_mom(EdgeIndex(0), GS.cind(0)) + Atom::one()).pow(-1);
        let source = Full4dCts::from_coefficient(&non_polynomial_term, &graph, &cograph);
        let root_spinney = graph.empty_subgraph();

        let error = localizer
            .project_4d(&source, &mut graph, &root_spinney)
            .expect_err("the exact term numerator must determine its own CFF support");
        assert!(
            error
                .to_string()
                .contains("could not analyze numerator in physical EMR energy variables")
        );
        Ok(())
    }

    #[test]
    fn typed_4d_projection_does_not_reanalyze_the_current_spinney_numerator() -> Result<()> {
        let mut graph: Graph = dot!(
            digraph current_owned_numerator {
                edge [pdg=1000 num=1 mass=1]
                node [num=1]
                a -> b [id=0 lmb_id=0 num="1/(Q(0,spenso::cind(0))+1)"]
                a -> b [id=1]
            },
            "scalars"
        )?;
        let current_spinney = graph.get_edge_subgraph(EdgeIndex(0));
        let reduced = graph
            .full_filter()
            .subtract(&current_spinney)
            .subtract(&graph.initial_state_cut);
        assert_eq!(
            graph
                .numerator(&reduced, &current_spinney)
                .get_single_atom()
                .expect("outside-spinney numerator is available"),
            Atom::one(),
            "the fixture must own its non-polynomial numerator only inside the current spinney"
        );

        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
        let source = Full4dCts::with_cograph(&Local4dCts::root(), &graph, &reduced);
        let projected = localizer.project_4d(&source, &mut graph, &current_spinney)?;
        let projected = projected.projected_integrands()?.materialize();

        assert!(projected.iter().any(|(_, atom)| !atom.is_zero()));
        Ok(())
    }

    #[test]
    fn contracted_uv_source_directions_use_production_selectors() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = vec![
            energy_map(
                edgevec([1, 1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), 1),
                    LinearEnergyExpr::ose(EdgeIndex(1), 1),
                ],
            ),
            energy_map(
                edgevec([-1, -1]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), -1),
                    LinearEnergyExpr::ose(EdgeIndex(1), -1),
                ],
            ),
        ]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let reduced = energy_map(
            edgevec([-1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), -1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        );
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let contract = graph.full_filter();
        let localized = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(0), EdgeIndex(1)],
            None,
            None,
        )?;

        assert_eq!(localized.len(), 1);
        assert_eq!(localized[0].0, OrientationID(0));
        assert_eq!(
            localized[0].1,
            GS.sign_theta(GS.sign(EdgeIndex(0))) * GS.sign_theta(GS.sign(EdgeIndex(1)))
        );
        Ok(())
    }

    #[test]
    fn typed_4d_projection_retains_external_tree_denominators() -> Result<()> {
        let _ = two_edge_graph()?;
        let mut graph: Graph = dot!(digraph external_tree {
            edge [particle="scalar_1"]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=4]
            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
            c -> d [id=3]
            d -> outgoing [id=5]
        }, "scalars")?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
        let cograph = graph.full_filter().subtract(&graph.initial_state_cut);
        let source = Full4dCts::with_cograph(&Local4dCts::root(), &graph, &cograph);
        let root_spinney = graph.empty_subgraph();

        let projected = localizer.project_4d(&source, &mut graph, &root_spinney)?;
        let projected = projected.projected_integrands()?.materialize();
        assert!(projected.iter().any(|(_, atom)| !atom.is_zero()));
        assert!(
            projected
                .iter()
                .all(|(_, atom)| !atom.contains_symbol(GS.den)),
            "tree propagators must remain exact scalar factors, not 4D wrapper atoms"
        );
        Ok(())
    }

    #[test]
    fn typed_4d_projection_supports_pure_trees() -> Result<()> {
        let _ = two_edge_graph()?;
        let mut graph: Graph = dot!(digraph pure_tree {
            edge [particle="scalar_1"]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
        }, "scalars")?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
        let cograph = graph.full_filter();
        let source = Full4dCts::with_cograph(&Local4dCts::root(), &graph, &cograph);
        let root_spinney = graph.empty_subgraph();

        let projected = localizer.project_4d(&source, &mut graph, &root_spinney)?;
        let root = CutCFFIndex::new_all_none();
        let terms = projected
            .projected_integrands()?
            .deferred_terms(&root)
            .expect("a pure tree has one denominatorless exact source term");
        assert_eq!(terms.len(), 1);
        assert!(!terms[0].is_zero());
        assert!(!terms[0].contains_symbol(GS.den));
        Ok(())
    }
}
