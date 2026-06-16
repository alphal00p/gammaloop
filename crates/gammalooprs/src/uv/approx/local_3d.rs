use std::{collections::BTreeMap, ops::Neg, sync::LazyLock};

use eyre::eyre;

use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, HedgePair, Orientation},
    subgraph::{SuBitGraph, SubSetOps},
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
    },
    debug_tags,
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet},
    settings::global::OrientationPattern,
    utils::{GS, W_},
    uv::{
        ApproximationType, Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{ForestNodeLike, OrientationProjection, UVCtx, integrated::IntegratedCts},
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

/// Residue integrands grouped by the one production energy map that owns every
/// numerator factor in the term. Operations zip equal IDs and never form a
/// Cartesian product of energy maps.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct OrientationIntegrands(BTreeMap<OrientationID, Integrands>);

impl From<Integrands> for OrientationIntegrands {
    fn from(integrands: Integrands) -> Self {
        Self(BTreeMap::from([(OrientationID(0), integrands)]))
    }
}

impl OrientationIntegrands {
    fn from_ids_and_indices(
        ids: impl IntoIterator<Item = OrientationID>,
        indices: &[CutCFFIndex],
    ) -> Self {
        Self(
            ids.into_iter()
                .map(|id| {
                    (
                        id,
                        indices.iter().map(|index| (*index, Atom::Zero)).collect(),
                    )
                })
                .collect(),
        )
    }

    fn checked_zip(
        &self,
        other: &Self,
        mut map: impl FnMut(OrientationID, &Integrands, &Integrands) -> Result<Integrands>,
    ) -> Result<Self> {
        if self.0.len() != other.0.len() {
            return Err(eyre!(
                "energy-map integrands have different map counts: {} and {}",
                self.0.len(),
                other.0.len()
            ));
        }
        self.0
            .iter()
            .map(|(id, left)| {
                let right = other
                    .0
                    .get(id)
                    .ok_or_else(|| eyre!("right integrands are missing energy map {}", id.0))?;
                Ok((*id, map(*id, left, right)?))
            })
            .collect::<Result<BTreeMap<_, _>>>()
            .map(Self)
    }

    pub(crate) fn zip_add(&self, other: &Self) -> Result<Self> {
        self.checked_zip(other, |_, left, right| left.zip_add(right))
    }

    fn zip_mul_unmapped(&self, other: &Integrands) -> Result<Self> {
        self.0
            .iter()
            .map(|(id, integrands)| Ok((*id, integrands.zip_mul(other)?)))
            .collect::<Result<BTreeMap<_, _>>>()
            .map(Self)
    }

    pub(crate) fn map(&self, mut f: impl FnMut(&Atom) -> Atom) -> Self {
        Self(
            self.0
                .iter()
                .map(|(id, integrands)| (*id, integrands.map(&mut f)))
                .collect(),
        )
    }

    fn fallible_map(
        &self,
        mut f: impl FnMut(OrientationID, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        self.0
            .iter()
            .map(|(id, integrands)| Ok((*id, integrands.fallible_map(|atom| f(*id, atom))?)))
            .collect::<Result<BTreeMap<_, _>>>()
            .map(Self)
    }

    #[cfg(test)]
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.0.values().flat_map(Integrands::iter)
    }

    pub(crate) fn iter_orientations(&self) -> impl Iterator<Item = (OrientationID, &Integrands)> {
        self.0.iter().map(|(id, integrands)| (*id, integrands))
    }
}

impl FromIterator<(OrientationID, Integrands)> for OrientationIntegrands {
    fn from_iter<T: IntoIterator<Item = (OrientationID, Integrands)>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl Neg for OrientationIntegrands {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(
            self.0
                .into_iter()
                .map(|(id, integrands)| (id, -integrands))
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

    fn cff(self, graph: &mut Graph, to_contract: &SuBitGraph) -> Result<(CutCFF, SuBitGraph)> {
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
        let analysis_numerator = self
            .orientation
            .exact_orientations()
            .is_some()
            .then(|| graph.production_numerator_atom_for_full_3d_expression());
        let cff = graph.cff(
            &contract_subgraph,
            self.cutset,
            orientation_pattern,
            &options,
            analysis_numerator.as_ref(),
        )?;
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
                    .is_compatible_with(&reduced.data.orientation)
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
                        .is_compatible_with(&reduced.data.orientation)
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
                            (full_energy != reduced_energy)
                                .then_some((*edge_id, full_energy, reduced_energy))
                        })
                        .collect::<Vec<_>>();
                    (
                        id,
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
                "no production energy map exactly extends reduced map {} after contracting edges {:?}; reduced provenance (label, numerator map, origins): ({:?}, {:?}, {:?}); {} production orientations are direction-compatible; first mismatch (orientation, label, numerator map, origins, edge/production/reduced): {first_mismatch:#?}",
                GS.orientation_delta(&reduced.data.orientation),
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

    fn coarse_representatives(
        self,
        reduced_orientation: &EdgeVec<Orientation>,
        internal_edges: &[EdgeIndex],
        average: bool,
    ) -> Result<Vec<OrientationID>> {
        let ids = self.orientation.orientation_ids();
        if ids.len() == 1 && self.orientation.orientation(ids[0]).is_none() {
            return Ok(ids);
        }
        let mut compatible = ids
            .into_iter()
            .filter(|id| {
                self.orientation
                    .orientation(*id)
                    .is_some_and(|orientation| orientation.is_compatible_with(reduced_orientation))
            })
            .collect::<Vec<_>>();

        if !average {
            compatible = compatible
                .into_iter()
                .max_by_key(|id| {
                    self.orientation
                        .orientation(*id)
                        .expect("coarse representative has an orientation")
                        .score(internal_edges)
                })
                .into_iter()
                .collect();
        }

        if compatible.is_empty() {
            Err(eyre!(
                "no valid global orientation matches reduced orientation {}",
                GS.orientation_delta(reduced_orientation)
            ))
        } else {
            Ok(compatible)
        }
    }

    /// Assign a reduced CFF map to exact production map IDs. The reduced
    /// selector remains attached to the denominator, while selectors for the
    /// integrated subgraph come from its full-map extensions.
    fn localized_orientation_terms(
        self,
        graph: &Graph,
        reduced: &OrientationExpression,
        reduced_expression: &Atom,
        contract_subgraph: &SuBitGraph,
        internal_edges: &[EdgeIndex],
        average: bool,
    ) -> Result<Vec<(OrientationID, Atom)>> {
        let representatives = if self.orientation.exact_orientations().is_some() {
            self.exact_representatives(graph, reduced, contract_subgraph)?
        } else {
            self.coarse_representatives(&reduced.data.orientation, internal_edges, average)?
        };
        if representatives.is_empty() {
            return Ok(Vec::new());
        }
        // Integrated localization splits a reduced term evenly among its full
        // extensions. Root CFF selectors already partition unity and retain
        // every full extension without this normalization.
        let normalization = if average {
            Atom::num(representatives.len() as i64)
        } else {
            Atom::one()
        };
        let reduced_selector = reduced.data.orientation.orientation_thetas();

        Ok(representatives
            .into_iter()
            .map(|id| {
                let internal_selector = self
                    .orientation
                    .orientation(id)
                    .map(|orientation| orientation.internal_orientation_selector(internal_edges))
                    .unwrap_or_else(Atom::one);
                (
                    id,
                    reduced_expression.clone() * &reduced_selector * internal_selector
                        / &normalization,
                )
            })
            .collect())
    }

    #[cfg(test)]
    fn localized_orientation_term(
        self,
        reduced_expression: &Atom,
        reduced_orientation: &EdgeVec<Orientation>,
        internal_edges: &[EdgeIndex],
        average: bool,
    ) -> Result<Atom> {
        let representatives =
            self.coarse_representatives(reduced_orientation, internal_edges, average)?;
        let normalization = if average {
            Atom::num(representatives.len() as i64)
        } else {
            Atom::one()
        };
        let reduced_selector = reduced_orientation.orientation_thetas();
        Ok(representatives
            .into_iter()
            .map(|id| {
                reduced_expression.clone()
                    * &reduced_selector
                    * self
                        .orientation
                        .orientation(id)
                        .map(|orientation| {
                            orientation.internal_orientation_selector(internal_edges)
                        })
                        .unwrap_or_else(Atom::one)
                    / &normalization
            })
            .reduce(|left, right| left + right)
            .unwrap_or(Atom::Zero))
    }

    fn projected_cff(
        self,
        graph: &mut Graph,
        to_contract: &SuBitGraph,
        average: bool,
    ) -> Result<OrientationIntegrands> {
        let (cff, contract_subgraph) = self.cff(graph, to_contract)?;
        let indices = cff.terms.keys().copied().collect::<Vec<_>>();
        let ids = self.orientation.orientation_ids();
        if ids.is_empty() {
            return Err(eyre!(
                "orientation pattern selects no production energy maps"
            ));
        }
        let mut terms = ids
            .iter()
            .flat_map(|id| indices.iter().map(move |index| ((*id, *index), Atom::Zero)))
            .collect::<BTreeMap<_, _>>();
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
                for (id, expression) in self.localized_orientation_terms(
                    graph,
                    &reduced.orientation,
                    &reduced.expression,
                    &contract_subgraph,
                    &internal_edges,
                    average,
                )? {
                    *terms.get_mut(&(id, index)).ok_or_else(|| {
                        eyre!("projected CFF returned unselected map {}", id.0)
                    })? += expression;
                }
            }
        }

        Ok(OrientationIntegrands(
            ids.into_iter()
                .map(|id| {
                    let integrands = indices
                        .iter()
                        .map(|index| {
                            (
                                *index,
                                terms
                                    .remove(&(id, *index))
                                    .expect("all projected CFF keys were initialized"),
                            )
                        })
                        .collect();
                    (id, integrands)
                })
                .collect(),
        ))
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

        let localizing_integrand = GS.localizing_integrand(integrated_node.lmb());
        let active = self.projected_cff(graph, to_contract, true)?.fallible_map(
            |orientation_id, localized| {
                let localized = localized * &fourddenoms;
                let localized_cff_byte_size = localized.as_view().get_byte_size();
                let mapped_finite_ct =
                    self.orientation
                        .map_numerator(graph, orientation_id, &finite_ct)?;
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
            },
        )?;
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
        numerator: &Atom,
    ) -> Result<Atom> {
        self.orientation
            .map_numerator(graph, orientation_id, numerator)
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
    integrands: OrientationIntegrands,
    // A disconnected join can leave a different set of loop variables active
    // in each local/integrated cross term.
    active_sectors: Option<Vec<(SuBitGraph, OrientationIntegrands)>>,
}

impl From<Integrands> for Local3DCts {
    fn from(integrands: Integrands) -> Self {
        Self {
            integrands: integrands.into(),
            active_sectors: None,
        }
    }
}

impl Neg for Local3DCts {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            integrands: self.integrands.map(|a| a.neg()),
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
        Ok(Self {
            integrands: self.integrands.zip_add(other)?,
            active_sectors: None,
        })
    }

    pub(crate) fn integrands(&self) -> &OrientationIntegrands {
        &self.integrands
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
            integrands,
            active_sectors: Some(active_sectors),
        })
    }

    pub(crate) fn map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        if let Some(active_sectors) = &self.active_sectors {
            let active_sectors = active_sectors
                .iter()
                .map(|(active, integrands)| {
                    Ok((active.clone(), integrands.fallible_map(|_, atom| f(atom))?))
                })
                .collect::<Result<_>>()?;
            Self::from_active_sectors(active_sectors)
        } else {
            Ok(Self {
                integrands: self.integrands.fallible_map(|_, atom| f(atom))?,
                active_sectors: None,
            })
        }
    }

    fn map_orientations(
        &self,
        mut f: impl FnMut(OrientationID, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        if let Some(active_sectors) = &self.active_sectors {
            let active_sectors = active_sectors
                .iter()
                .map(|(active, integrands)| Ok((active.clone(), integrands.fallible_map(&mut f)?)))
                .collect::<Result<_>>()?;
            Self::from_active_sectors(active_sectors)
        } else {
            Ok(Self {
                integrands: self.integrands.fallible_map(f)?,
                active_sectors: None,
            })
        }
    }

    pub(crate) fn root(graph: &mut Graph, localizer: Localizer<'_>) -> Result<Self> {
        let cff = localizer.projected_cff(graph, &graph.empty_subgraph::<SuBitGraph>(), false)?;

        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        Ok(Local3DCts {
            integrands: cff.map(|atom| atom * &fourddenoms),
            active_sectors: None,
        })
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
pub(crate) fn t_tilde<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    orientation_id: OrientationID,
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
        .map_numerator(graph, orientation_id, &numerator)?
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

pub(crate) fn start<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    orientation_id: OrientationID,
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
    let mut atomarg = cff * orientation.map_numerator(graph, orientation_id, &numerator)?;
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
) -> impl FnMut(OrientationID, &Atom) -> Result<Atom> {
    Local3DLoopRescaling::FullSubgraph.map(ctx, orientation, current, given, None)
}
fn reduced<S: super::ForestNodeLike>(
    ctx: &UVCtx<'_>,
    orientation: OrientationProjection<'_>,
    current: &S,
    given: &S,
) -> impl FnMut(OrientationID, &Atom) -> Result<Atom> {
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
    ) -> impl FnMut(OrientationID, &Atom) -> Result<Atom> {
        move |orientation_id, integrand| {
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
    use crate::{
        cff::{
            expression::{
                GammaLoopOrientationExpression, OrientationData, OrientationExpression,
                OrientationID,
            },
            surface::LinearEnergyExpr,
        },
        dot,
        graph::{FeynmanGraph, Graph, cuts::CutSet, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        utils::GS,
        uv::{
            Spinney,
            approx::{OrientationProjection, local_3d::Localizer},
            hedge_poset::OwnedForestNode,
        },
    };
    use color_eyre::Result;
    use eyre::eyre;
    use linnet::half_edge::{
        involution::{EdgeIndex, EdgeVec, Orientation},
        subgraph::SubSetOps,
    };
    use std::sync::OnceLock;
    use symbolica::{
        atom::{Atom, AtomCore},
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
    fn orientation_term_keeps_external_selectors_and_adds_internal_ones() {
        let reduced_expression = function!(GS.ose, 0);
        let reduced_orientation = edgevec([1, 0, -1]);
        let valid = vec![edgevec([1, 1, -1]), edgevec([1, -1, -1])];
        let pat = OrientationPattern::default();
        let cutset = CutSet::empty(1);
        let loc = Localizer::new(&cutset, OrientationProjection::new(&valid, &pat));
        let localized = loc
            .localized_orientation_term(
                &reduced_expression,
                &reduced_orientation,
                &edges([1]),
                false,
            )
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
            OrientationProjection::exact(&production, &options, &pattern),
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
            OrientationProjection::exact(&production, &options, &pattern),
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
            OrientationProjection::exact(&production, &options, &pattern),
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
                    true,
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
            OrientationProjection::exact(&production, &options, &pattern),
        );
        let to_contract = graph.get_edge_subgraph(EdgeIndex(1));
        let projected = localizer.projected_cff(&mut graph, &to_contract, true)?;

        assert_eq!(
            projected
                .iter_orientations()
                .map(|(id, _)| id)
                .collect::<Vec<_>>(),
            vec![OrientationID(0)]
        );
        assert!(projected.iter().any(|(_, atom)| !atom.is_zero()));
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
            .generate_3d_expression_for_integrand(
                &contract_edges,
                &canonization,
                &options,
                None,
                false,
            )?
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
            OrientationProjection::exact(&production, &options, &pattern),
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
            .find_map(|(id, integrands)| (id == orientation_id).then_some(integrands))
            .expect("the selected production map has a baseline branch");
        let localized = localized
            .active
            .iter_orientations()
            .find_map(|(id, integrands)| (id == orientation_id).then_some(integrands))
            .expect("the selected production map has a finite-CT branch");

        assert_eq!(localized, &baseline.map(|atom| atom * &mapped_finite_ct));
        assert_ne!(localized, &baseline.map(|atom| atom * &finite_ct));
        Ok(())
    }

    #[test]
    fn contracted_exact_extensions_are_unweighted_at_root_and_averaged_when_localized() -> Result<()>
    {
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
            OrientationProjection::exact(&production, &options, &pattern),
        );
        let contract = graph.get_edge_subgraph(EdgeIndex(1));
        let root = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            false,
        )?;
        let localized = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            true,
        )?;
        let expected_default =
            GS.sign_theta(GS.sign(EdgeIndex(0))) * GS.sign_theta(GS.sign(EdgeIndex(1)));
        let expected_reversed =
            GS.sign_theta(GS.sign(EdgeIndex(0))) * GS.sign_theta(-GS.sign(EdgeIndex(1)));

        assert_eq!(
            root.iter().map(|(id, _)| *id).collect::<Vec<_>>(),
            vec![OrientationID(0), OrientationID(1)]
        );
        assert_eq!(root[0].1, expected_default);
        assert_eq!(root[1].1, expected_reversed);
        assert_eq!(localized[0].1, root[0].1.clone() / 2);
        assert_eq!(localized[1].1, root[1].1.clone() / 2);
        Ok(())
    }
}
