use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Neg,
};

use eyre::eyre;

use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore},
    symbol,
};
use three_dimensional_reps::CffGenerationContext;

use crate::{
    cff::{
        CutCFF, CutCFFIndex,
        expression::{OrientationExpression, OrientationID, OrientationSelector},
        orientations::GraphOrientation,
        surface::{GammaLoopLinearEnergyExpr, LinearEnergyExpr},
    },
    debug_tags,
    graph::{Graph, GraphThreeDSource, LoopMomentumBasis, cuts::CutSet},
    momentum::SignOrZero,
    settings::global::OrientationPattern,
    utils::GS,
    uv::{
        Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{ForestNodeLike, OrientationProjection, Rooted, direct_3d::Direct3dCts},
    },
};
use color_eyre::Result;

#[cfg(test)]
use crate::uv::{DeferredIntegrands, approx::local_4d::Full4dCts};
#[cfg(test)]
use linnet::half_edge::subgraph::Inclusion;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct FrozenActiveCt {
    pub active: OrientationIntegrands,
    pub frozen_integrands: Integrands,
    /// The coordinate frame in which this sector's still-active Taylor
    /// coefficient was formed. Direct complete-CFF sectors have no such
    /// independent 4D frame.
    pub active_lmb: Option<LoopMomentumBasis>,
}

impl FrozenActiveCt {
    pub(crate) fn combine(&self) -> Result<OrientationIntegrands> {
        self.active.zip_mul_unmapped(&self.frozen_integrands)
    }
}

impl From<OrientationIntegrands> for FrozenActiveCt {
    fn from(active: OrientationIntegrands) -> Self {
        let frozen_integrands = active
            .0
            .first()
            .map(|branch| {
                branch
                    .integrands
                    .iter()
                    .map(|(index, _)| (*index, Atom::one()))
                    .collect()
            })
            .unwrap_or_else(Integrands::root);
        Self {
            active,
            frozen_integrands,
            active_lmb: None,
        }
    }
}

impl From<Integrands> for FrozenActiveCt {
    fn from(integrands: Integrands) -> Self {
        OrientationIntegrands::from(integrands).into()
    }
}

impl Neg for FrozenActiveCt {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            active: -self.active,
            frozen_integrands: self.frozen_integrands,
            active_lmb: self.active_lmb,
        }
    }
}

/// Residue integrands grouped by the selector and exact energy map that own
/// every numerator factor in the term. A missing source map uses the selected
/// production map; a present source map remains authoritative while production
/// IDs only partition exact residue-map-key hosts.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(super) struct OrientationIntegrandBranch {
    pub(super) selector_id: OrientationID,
    pub(super) source_edge_energy_map: Option<Vec<LinearEnergyExpr>>,
    pub(super) integrands: Integrands,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct OrientationIntegrands(pub(super) Vec<OrientationIntegrandBranch>);

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

    pub(crate) fn zip_mul_unmapped(&self, other: &Integrands) -> Result<Self> {
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

    pub(crate) fn multiply_mapped(
        &self,
        mut map: impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>) -> Result<Atom>,
    ) -> Result<Self> {
        self.0
            .iter()
            .map(|branch| {
                let mapped = map(branch.selector_id, branch.source_edge_energy_map.as_deref())?;
                Ok(OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.map(|atom| atom * &mapped),
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
    }

    /// Multiply branches hosted by the same production selector. A branch
    /// absent on either side contributes zero, so independently projected
    /// factors retain only their selector intersection. The mapper sees the
    /// complete active factor only after its host has been selected.
    pub(crate) fn zip_mul_mapped_factor(
        &self,
        factor: &Self,
        mut map: impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        self.0
            .iter()
            .filter_map(|outer| {
                let matching = factor
                    .0
                    .iter()
                    .filter(|inner| inner.selector_id == outer.selector_id)
                    .collect::<Vec<_>>();
                (!matching.is_empty()).then(|| {
                    let mut product = outer.integrands.zero_like();
                    for inner in matching {
                        let mapped = inner.integrands.fallible_map(|atom| {
                            map(
                                outer.selector_id,
                                outer.source_edge_energy_map.as_deref(),
                                atom,
                            )
                        })?;
                        product = product.zip_add(outer.integrands.zip_mul(&mapped)?)?;
                    }
                    Ok(OrientationIntegrandBranch {
                        selector_id: outer.selector_id,
                        source_edge_energy_map: outer.source_edge_energy_map.clone(),
                        integrands: product,
                    })
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
    }

    /// A factorized additive projection for coefficient diagnostics. Repeated
    /// selector hosts of the same active expression are included only once.
    #[cfg(test)]
    pub(crate) fn factorized_sum(&self) -> Atom {
        let mut distinct = Vec::new();
        for atom in self
            .0
            .iter()
            .flat_map(|branch| branch.integrands.iter().map(|(_, atom)| atom))
        {
            if !distinct.contains(atom) {
                distinct.push(atom.clone());
            }
        }
        distinct
            .into_iter()
            .fold(Atom::Zero, |sum, atom| sum + atom)
    }

    /// Keep independently evaluated production branches algebraically
    /// independent while deriving one conservative outer-CFF capacity.
    /// Branch tags are analysis-only scalar coefficients: they neither expand
    /// the factorized atoms nor enter the mapped production numerator.
    pub(crate) fn factorized_capacity_envelope(&self) -> Atom {
        self.0
            .iter()
            .flat_map(|branch| branch.integrands.iter().map(|(_, atom)| atom))
            .filter(|atom| !atom.is_zero())
            .enumerate()
            .fold(Atom::Zero, |sum, (branch, atom)| {
                let tag = Atom::var(symbol!(format!(
                    "__gammaloop_outer_cff_capacity_branch_{branch}"
                )));
                sum + tag * atom
            })
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
    pub(super) cutset: &'a CutSet,
    pub(super) orientation: OrientationProjection<'a>,
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
    fn localized_orientation_terms(
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
                active_lmb: None,
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
            .projected_cff(
                graph,
                to_contract,
                &analysis_numerator,
                CffGenerationContext::Standalone,
            )?
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
                    integrated_node = %integrated_node.log_display(),
                    contracted = %to_contract.string_label(),
                    reduced = %reduced.string_label(),
                    residue_map_key = orientation_id.0,
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
            active_lmb: None,
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

    /// Materialize the opaque selector carried by an exact production-map
    /// branch. Explicit sums and coarse projectors have no runtime map key.
    pub(crate) fn residue_map_key_selector(self, id: OrientationID) -> Atom {
        if self.uses_exact_maps() && !self.orientation.explicit_orientation_sum_only {
            id.atom()
        } else {
            Atom::one()
        }
    }
}

pub(crate) struct Local3DApproximation<'a> {
    pub(super) localizer: Localizer<'a>,
    pub(super) graph: &'a mut Graph,
    pub(super) settings: &'a UVgenerationSettings,
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
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Local3DCts {
    /// The typed direct lane. Its complete generalized residue keys remain
    /// opaque while the Taylor forest is built.
    Direct(Direct3dCts),
    /// Projected local-4D Taylor coefficients. They omit the untouched outer
    /// CFF, which is attached only during final assembly.
    Projected4d(Vec<(SuBitGraph, FrozenActiveCt)>),
    /// Exact whole-expression projection remains a diagnostic-only route.
    #[cfg(test)]
    Projected(Integrands),
}

impl Neg for Local3DCts {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Self::Direct(direct) => Self::Direct(-direct),
            Self::Projected4d(sectors) => Self::Projected4d(
                sectors
                    .into_iter()
                    .map(|(active, integrands)| (active, -integrands))
                    .collect(),
            ),
            #[cfg(test)]
            Self::Projected(integrands) => Self::Projected(-integrands),
        }
    }
}

impl Local3DCts {
    pub(crate) fn direct(&self) -> Result<&Direct3dCts> {
        match self {
            Self::Direct(direct) => Ok(direct),
            _ => Err(eyre!(
                "projected local-4D coefficients are not direct local-3D counterterms"
            )),
        }
    }

    pub(crate) fn projected_4d_sectors(&self) -> Result<&[(SuBitGraph, FrozenActiveCt)]> {
        match self {
            Self::Projected4d(sectors) => Ok(sectors),
            Self::Direct(_) => Err(eyre!(
                "direct complete-CFF terms are not projected local-4D coefficients"
            )),
            #[cfg(test)]
            Self::Projected(_) => Err(eyre!(
                "diagnostic whole-expression projections are not projected local-4D coefficients"
            )),
        }
    }

    #[cfg(test)]
    pub(crate) fn projected_integrands(&self) -> Result<&Integrands> {
        match self {
            Self::Projected(integrands) => Ok(integrands),
            _ => Err(eyre!(
                "direct 3D terms are grouped by production orientation IDs"
            )),
        }
    }

    #[cfg(test)]
    pub(crate) fn map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        let map_sectors = |sectors: &[(SuBitGraph, FrozenActiveCt)],
                           f: &mut F|
         -> Result<Vec<(SuBitGraph, FrozenActiveCt)>> {
            sectors
                .iter()
                .map(|(active, integrands)| {
                    Ok((
                        active.clone(),
                        FrozenActiveCt {
                            active: integrands.active.fallible_map(|_, _, atom| f(atom))?,
                            frozen_integrands: integrands.frozen_integrands.clone(),
                            active_lmb: integrands.active_lmb.clone(),
                        },
                    ))
                })
                .collect()
        };
        match self {
            Self::Direct(direct) => Ok(Self::Direct(direct.map(&mut f)?)),
            Self::Projected4d(sectors) => Ok(Self::Projected4d(map_sectors(sectors, &mut f)?)),
            #[cfg(test)]
            Self::Projected(integrands) => Ok(Self::Projected(integrands.fallible_map(f)?)),
        }
    }

    pub(crate) fn root(graph: &mut Graph, localizer: Localizer<'_>) -> Result<Self> {
        Ok(Self::Direct(Direct3dCts::root(graph, localizer)?))
    }

    #[cfg(test)]
    fn from_projected(integrands: Integrands) -> Self {
        Self::Projected(integrands)
    }
}

#[cfg(test)]
mod tests {
    use super::{Local3DCts, OrientationIntegrandBranch, OrientationIntegrands};
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
        graph::{FeynmanGraph, FourDDenominator, Graph, LMBext, cuts::CutSet, parse::IntoGraph},
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
    use std::{
        collections::{BTreeMap, BTreeSet},
        sync::OnceLock,
    };
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, FunctionBuilder},
        function,
    };
    use three_dimensional_reps::CffGenerationContext;
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
        let powered = localizer.projected_cff(
            &mut graph,
            &contract,
            &numerator,
            CffGenerationContext::Standalone,
        )?;
        let mut explicit_sum = Atom::Zero;
        for (selector_id, source_map, integrands) in powered.iter_orientations() {
            let mapped = localizer.map_numerator(&graph, selector_id, source_map, &numerator)?;
            explicit_sum += integrands
                .iter()
                .fold(Atom::Zero, |sum, (_, term)| sum + term * &mapped);
        }
        // The two parallel physical edges carry opposite routing signs. Their
        // acyclic production sectors are therefore (+,-) and (-,+); (+,+) is
        // not a production orientation. The cancellation contacts have zero
        // numerator sampling maps, but their opposite loop lifts still host
        // one matching lower residue in each sector.
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
    fn direct_root_preserves_powered_production_entries() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let edge = EdgeIndex(0);
        let on_shell_energy_squared = (1..=3).fold(
            graph.underlying[edge].particle.mass_atom().pow(2),
            |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
            },
        );
        let numerator = GS.emr_mom(edge, GS.cind(0)).pow(2) - on_shell_energy_squared;
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
        let raw_root = graph.cff_from_production_expression(&production, &cutset, &pattern)?;
        let production_prefactor = Atom::num(raw_root.production_prefactor_factor());
        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );

        // The direct root is the stored production CFF itself. A generalized
        // numerator sample may leave an edge undirected, but that is part of
        // its authoritative contour and cannot trigger representative
        // reconstruction. Both generation modes keep that complete map key as
        // an opaque branch while the direct Taylor forest is replayed.
        for explicit_orientation_sum_only in [false, true] {
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &pattern,
                    explicit_orientation_sum_only,
                ),
            );
            let localized = Local3DCts::root(&mut graph, localizer)?
                .direct()?
                .branches()?;
            let expected_keys = raw_root
                .terms
                .values()
                .flat_map(|term| &term.orientations)
                .map(|raw| {
                    (
                        raw.production_orientation_id
                            .expect("stored production CFF entries retain their original ID"),
                        raw.orientation.edge_energy_map.clone(),
                    )
                })
                .fold(Vec::new(), |mut keys, key| {
                    if !keys.contains(&key) {
                        keys.push(key);
                    }
                    keys
                });
            let actual_keys = localized
                .iter_keys()
                .map(|(key, _)| {
                    (
                        key.selector_host,
                        production.expression.orientations[key.selector_host]
                            .edge_energy_map
                            .clone(),
                    )
                })
                .collect::<Vec<_>>();
            assert_eq!(actual_keys.len(), expected_keys.len());
            assert!(expected_keys.iter().all(|key| actual_keys.contains(key)));
            for (index, term) in &raw_root.terms {
                for raw in &term.orientations {
                    let id = raw
                        .production_orientation_id
                        .expect("stored production CFF entries retain their original ID");
                    let original = &production.expression.orientations[id];
                    assert_eq!(raw.orientation.loop_energy_map, original.loop_energy_map);
                    assert_eq!(raw.orientation.edge_energy_map, original.edge_energy_map);
                    let (key, integrands) = localized
                        .iter_keys()
                        .find(|(key, _)| key.selector_host == id)
                        .expect("the direct root preserves each stored ID and exact energy map");
                    assert_eq!(key.source_edge_energy_map(), None);
                    let actual = integrands
                        .iter()
                        .find_map(|(actual_index, atom)| (actual_index == index).then_some(atom))
                        .expect("the direct root preserves every selected CutCFF index");
                    let expected = &raw.expression * &production_prefactor * &fourddenoms;
                    assert_eq!(actual, &expected);
                }
            }
        }
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
    fn zero_sampling_maps_retain_opposite_loop_lift_sectors() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = [[1, -1], [-1, 1]]
            .into_iter()
            .map(|directions| {
                energy_map(
                    edgevec(directions),
                    directions
                        .into_iter()
                        .enumerate()
                        .map(|(edge, direction)| {
                            LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                        })
                        .collect(),
                )
            })
            .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let surviving_edge = EdgeIndex(1);
        let mut hosts = Vec::new();

        for loop_energy in [
            LinearEnergyExpr::ose(surviving_edge, 1),
            LinearEnergyExpr::ose(surviving_edge, -1),
        ] {
            let mut reduced = energy_map(
                edgevec([0, 0]),
                vec![LinearEnergyExpr::zero(), LinearEnergyExpr::zero()],
            );
            reduced.loop_energy_map = vec![loop_energy];
            let representatives = localizer.source_selector_representatives(
                &graph,
                &graph.loop_momentum_basis,
                &reduced,
                &graph.empty_subgraph(),
                None,
            )?;
            assert_eq!(representatives.len(), 1);
            hosts.push(representatives[0]);
        }

        assert!(hosts.contains(&OrientationID(0)));
        assert!(hosts.contains(&OrientationID(1)));
        assert_ne!(hosts[0], hosts[1]);
        assert_ne!(
            production[hosts[0]].data.orientation[surviving_edge],
            production[hosts[1]].data.orientation[surviving_edge],
            "opposite loop-pole lifts must host opposite surviving-edge sectors"
        );
        Ok(())
    }

    #[test]
    fn source_selector_interprets_loop_lift_in_source_generation_lmb() -> Result<()> {
        let graph = two_edge_graph()?;
        let source_lmb = graph
            .generate_loop_momentum_bases_of(&graph.full_filter())
            .into_iter()
            .find(|lmb| lmb.loop_edges != graph.loop_momentum_basis.loop_edges)
            .ok_or_else(|| eyre!("the selector-hosting fixture needs a non-global LMB"))?;
        let source_edge = *source_lmb
            .loop_edges
            .first()
            .expect("the two-edge graph has one loop");
        let production = [[1, -1], [-1, 1]]
            .into_iter()
            .map(|directions| {
                energy_map(
                    edgevec(directions),
                    directions
                        .into_iter()
                        .enumerate()
                        .map(|(edge, direction)| {
                            LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                        })
                        .collect(),
                )
            })
            .collect::<TiVec<OrientationID, _>>();
        let mut reduced = energy_map(
            edgevec([0, 0]),
            vec![LinearEnergyExpr::zero(), LinearEnergyExpr::zero()],
        );
        reduced.loop_energy_map = vec![LinearEnergyExpr::ose(source_edge, 1)];
        let interpret_source_edge = |lmb: &crate::graph::LoopMomentumBasis| {
            lmb.edge_signatures[source_edge]
                .try_compute_momentum(&reduced.loop_energy_map, &[])
                .expect("the source loop edge depends on its loop coordinate")
                .canonical()
        };
        assert_eq!(
            interpret_source_edge(&source_lmb),
            LinearEnergyExpr::ose(source_edge, 1),
        );
        assert_ne!(
            interpret_source_edge(&graph.loop_momentum_basis),
            LinearEnergyExpr::ose(source_edge, 1),
            "the fixture must distinguish the child source LMB from the graph LMB",
        );

        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let expected_host = production
            .iter_enumerated()
            .find_map(|(id, map)| {
                (map.data.orientation[source_edge] == Orientation::Default).then_some(id)
            })
            .expect("one production sector directs the source carrier forward");
        let valid_hosts = BTreeSet::from([OrientationID(0), OrientationID(1)]);
        let hosted = localizer.localized_orientation_terms(
            &graph,
            &source_lmb,
            &reduced,
            &Atom::one(),
            &graph.empty_subgraph(),
            &[],
            Some(&valid_hosts),
            None,
            Some(&reduced.edge_energy_map),
            None,
        )?;

        assert_eq!(
            hosted,
            vec![(expected_host, Atom::one())],
            "explicit direct-3D replay must retain the child-LMB host as partition metadata",
        );
        Ok(())
    }

    #[test]
    fn nested_powered_contact_keeps_loop_lift_as_provenance_only() -> Result<()> {
        let mut graph = two_edge_graph()?;
        let powered_edge = EdgeIndex(0);
        let on_shell_energy_squared = (1..=3).fold(
            graph.underlying[powered_edge].particle.mass_atom().pow(2),
            |norm_squared, spatial_index| {
                norm_squared + GS.emr_mom(powered_edge, GS.cind(spatial_index)).pow(2)
            },
        );
        let numerator = GS.emr_mom(powered_edge, GS.cind(0)).pow(2) - on_shell_energy_squared;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonization,
            &options,
            Some(&Atom::one()),
        )?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(&production, &options, &pattern, true),
        );
        let source_denominator = || FourDDenominator {
            source_edge: powered_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(powered_edge))
                .finish(),
            mass_squared: graph.underlying[powered_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };
        let denominators = [source_denominator(), source_denominator()];
        let (exact, contract_subgraph) =
            graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
        let internal_edges = graph.paired_edges(&contract_subgraph);
        let term = exact
            .terms
            .values()
            .next()
            .expect("the uncut powered source has one CFF term");
        let mut shared_hosts = None;
        let mut contact_count = 0;
        let mut remainder_count = 0;
        let mut saw_zero_canonical_contact = false;
        for orientation in &term.orientations {
            let is_contact = orientation.orientation.variants.iter().any(|variant| {
                variant
                    .origin
                    .as_deref()
                    .is_some_and(|origin| origin.contains("contact"))
            });
            let is_remainder = orientation.orientation.variants.iter().any(|variant| {
                variant
                    .origin
                    .as_deref()
                    .is_some_and(|origin| origin.contains("remainder"))
            });
            assert_ne!(is_contact, is_remainder);
            let canonical_sample = term.map_exact_source_atom(
                &orientation.orientation,
                &GS.emr_mom(powered_edge, GS.cind(0)),
            )?;
            let physical_loop_lift_energies = term.map_exact_source_physical_loop_lift_energies(
                &orientation.orientation,
                [powered_edge],
            )?;
            let powered_loop_lift = physical_loop_lift_energies
                .iter()
                .find_map(|(edge, energy)| (*edge == powered_edge).then_some(energy));
            let powered_loop_lift = powered_loop_lift
                .expect("the contact owner is reconstructed from the exact source loop lift");
            assert!(
                orientation
                    .orientation
                    .loop_energy_map
                    .iter()
                    .any(|energy| energy.clone().canonical() != LinearEnergyExpr::zero()),
                "each powered-source branch must retain a nonzero residue loop lift"
            );
            assert!(!powered_loop_lift.is_zero());
            if is_contact && canonical_sample.is_zero() {
                saw_zero_canonical_contact = true;
                assert_eq!(powered_loop_lift, &GS.ose(powered_edge));
            }
            let hosts = localizer.source_selector_representatives(
                &graph,
                &graph.loop_momentum_basis,
                &orientation.orientation,
                &contract_subgraph,
                None,
            )?;
            let host_set = hosts.iter().copied().collect::<BTreeSet<_>>();
            if let Some(expected) = &shared_hosts {
                assert_eq!(
                    &host_set, expected,
                    "synthetic contact and remainder directions must not partition physical production hosts"
                );
            } else {
                shared_hosts = Some(host_set);
            }
            if is_contact {
                contact_count += 1;
            } else {
                remainder_count += 1;
            }
            let selected_host = hosts[0];
            let valid_host = BTreeSet::from([selected_host]);
            let hosted = localizer.localized_orientation_terms(
                &graph,
                &graph.loop_momentum_basis,
                &orientation.orientation,
                &Atom::one(),
                &contract_subgraph,
                &internal_edges,
                Some(&valid_host),
                orientation.production_orientation_id,
                Some(&orientation.orientation.edge_energy_map),
                None,
            )?;
            assert_eq!(hosted.len(), 1);
            assert_eq!(hosted[0].0, selected_host);
        }
        assert!(
            saw_zero_canonical_contact,
            "at least one recursive contact must expose the stale zero numerator sample"
        );
        assert_eq!(contact_count, 3);
        assert_eq!(remainder_count, 2);
        assert_eq!(
            shared_hosts
                .expect("the exact powered source has contour pieces")
                .len(),
            2,
            "every algebraic contour piece may be placed under either physical production host; the caller selects one deterministically"
        );
        Ok(())
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
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            None,
            Some(&reduced.edge_energy_map),
            None,
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
            "production IDs partition complete residue-map keys but do not remap source-owned numerators"
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
    fn factorized_products_keep_only_shared_production_hosts() -> Result<()> {
        let index = CutCFFIndex::new_all_none();
        let branch = |selector, value| OrientationIntegrandBranch {
            selector_id: OrientationID(selector),
            source_edge_energy_map: None,
            integrands: [(index, Atom::num(value))].into_iter().collect(),
        };
        let left = OrientationIntegrands(vec![branch(0, 2), branch(1, 3)]);
        let right = OrientationIntegrands(vec![branch(1, 5), branch(2, 7)]);
        let mut mapped_hosts = Vec::new();
        let product = left.zip_mul_mapped_factor(&right, |host, _, factor| {
            mapped_hosts.push(host);
            Ok(factor.clone())
        })?;

        assert_eq!(mapped_hosts, vec![OrientationID(1)]);
        assert_eq!(
            product
                .iter_orientations()
                .map(|(host, _, _)| host)
                .collect::<Vec<_>>(),
            vec![OrientationID(1)],
            "factorized multiplication must neither retain nor synthesize unmatched hosts"
        );
        assert_eq!(product.iter().next(), Some((&index, &Atom::num(15))));
        Ok(())
    }

    fn energy_bounds(graph: &Graph, atom: &Atom) -> Result<Vec<(usize, usize)>> {
        Ok(
            graph.automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
                atom,
                [],
                1,
            )?,
        )
    }

    #[test]
    fn outer_cff_capacity_does_not_cancel_between_selector_branches() -> Result<()> {
        let graph = two_edge_graph()?;
        let energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let factorized = (energy.clone() + Atom::num(1)) * (energy.clone() + Atom::num(2));
        let root = CutCFFIndex::new_all_none();
        let branches = OrientationIntegrands(vec![
            OrientationIntegrandBranch {
                selector_id: OrientationID(0),
                source_edge_energy_map: None,
                integrands: [(root, factorized.clone())].into_iter().collect(),
            },
            OrientationIntegrandBranch {
                selector_id: OrientationID(1),
                source_edge_energy_map: None,
                integrands: [(root, -&factorized + &energy)].into_iter().collect(),
            },
        ]);

        assert_eq!(
            energy_bounds(&graph, &branches.factorized_sum())?,
            vec![(0, 1)]
        );
        assert_eq!(
            energy_bounds(&graph, &branches.factorized_capacity_envelope())?,
            vec![(0, 2)],
            "mutually exclusive selector branches need the maximum of their separate ranks"
        );
        assert_eq!(
            branches
                .iter_orientations()
                .next()
                .and_then(|(_, _, integrands)| integrands.iter().next())
                .map(|(_, atom)| atom),
            Some(&factorized),
            "capacity analysis must not expand or rewrite the stored factorized numerator"
        );
        Ok(())
    }

    #[test]
    fn outer_cff_capacity_does_not_cancel_between_cut_orders() -> Result<()> {
        let graph = two_edge_graph()?;
        let energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let cubic = (energy.clone() + Atom::num(1))
            * (energy.clone() + Atom::num(2))
            * (energy.clone() + Atom::num(3));
        let raised = CutCFFIndex {
            lu_cut_order: Some(1),
            ..CutCFFIndex::new_all_none()
        };
        let branches = OrientationIntegrands(vec![OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: None,
            integrands: [
                (CutCFFIndex::new_all_none(), cubic.clone()),
                (raised, -&cubic + &energy),
            ]
            .into_iter()
            .collect(),
        }]);

        assert_eq!(
            energy_bounds(&graph, &branches.factorized_sum())?,
            vec![(0, 1)]
        );
        assert_eq!(
            energy_bounds(&graph, &branches.factorized_capacity_envelope())?,
            vec![(0, 3)],
            "separately evaluated CutCFFIndex values need the maximum of their separate ranks"
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
                    &graph.loop_momentum_basis,
                    &reduced,
                    &Atom::one(),
                    &contract,
                    &[EdgeIndex(1)],
                    None,
                    None,
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
        let projected = localizer.projected_cff(
            &mut graph,
            &to_contract,
            &analysis_numerator,
            CffGenerationContext::Standalone,
        )?;

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

        let projected = localizer.projected_cff(
            &mut graph,
            &to_contract,
            &unsupported_numerator,
            CffGenerationContext::Standalone,
        )?;

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
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            None,
            None,
            None,
        )?;
        assert_eq!(localized, vec![(OrientationID(0), Atom::one())]);
        assert_eq!(
            localizer.residue_map_key_selector(localized[0].0),
            OrientationID(0).atom(),
            "an exact reduced residue is localized by its complete map key, not by a physical-theta product"
        );

        let explicit_localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let explicit = explicit_localizer.localized_orientation_terms(
            &graph,
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            None,
            None,
            None,
        )?;
        assert_eq!(explicit, vec![(OrientationID(0), Atom::one())]);
        Ok(())
    }

    #[test]
    fn cut_valid_ids_host_one_inner_representative_per_outer_sector() -> Result<()> {
        let graph = two_edge_graph()?;
        let production = [[1, 1], [1, -1], [-1, 1], [-1, -1]]
            .into_iter()
            .map(|directions| {
                energy_map(
                    edgevec(directions),
                    directions
                        .into_iter()
                        .enumerate()
                        .map(|(edge, direction)| {
                            LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                        })
                        .collect(),
                )
            })
            .collect::<TiVec<OrientationID, _>>();
        let reduced = [1, -1].map(|outer_direction| {
            energy_map(
                edgevec([outer_direction, 0]),
                vec![
                    LinearEnergyExpr::ose(EdgeIndex(0), i64::from(outer_direction)),
                    LinearEnergyExpr::zero(),
                ],
            )
        });
        let valid_ids = BTreeSet::from([OrientationID(1), OrientationID(2)]);
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let contract = graph.get_edge_subgraph(EdgeIndex(1));

        for explicit_orientation_sum_only in [false, true] {
            let localizer = Localizer::new(
                &cutset,
                OrientationProjection::exact(
                    &production,
                    &options,
                    &pattern,
                    explicit_orientation_sum_only,
                ),
            );
            let localized = reduced
                .iter()
                .enumerate()
                .map(|(outer_sector, reduced)| {
                    let body =
                        Atom::var(symbolica::symbol!(format!("cut_valid_body_{outer_sector}")));
                    let terms = localizer.localized_orientation_terms(
                        &graph,
                        &graph.loop_momentum_basis,
                        reduced,
                        &body,
                        &contract,
                        &[EdgeIndex(1)],
                        Some(&valid_ids),
                        None,
                        Some(&reduced.edge_energy_map),
                        None,
                    )?;
                    Ok((body, terms))
                })
                .collect::<Result<Vec<_>>>()?;

            assert!(localized.iter().all(|(_, terms)| terms.len() == 1));
            assert_eq!(localized[0].1[0].0, OrientationID(1));
            assert_eq!(localized[1].1[0].0, OrientationID(2));
            assert!(
                localized
                    .iter()
                    .flat_map(|(_, terms)| terms)
                    .all(|(id, _)| {
                        valid_ids.contains(id) && *id != OrientationID(0) && *id != OrientationID(3)
                    })
            );

            for (outer_sector, (body, terms)) in localized.iter().enumerate() {
                let (id, expression) = &terms[0];
                assert_eq!(expression, body);
                if explicit_orientation_sum_only {
                    assert_eq!(localizer.residue_map_key_selector(*id), Atom::one());
                } else {
                    let selector = localizer.residue_map_key_selector(*id);
                    assert_eq!(
                        id.select(selector.as_view()),
                        Atom::one(),
                        "the valid representative owns its complete map-key sector exactly once"
                    );
                    let invalid_id = [OrientationID(0), OrientationID(3)][outer_sector];
                    assert_eq!(
                        invalid_id.select(selector.as_view()),
                        Atom::Zero,
                        "a cut-invalid map key must not receive the body"
                    );
                }
            }
        }
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
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            Some(OrientationID(1)),
            None,
            None,
        )?;

        assert_eq!(localized.len(), 1);
        assert_eq!(localized[0].0, OrientationID(1));
        assert_eq!(localized[0].1, Atom::one());
        let selected = localizer.residue_map_key_selector(localized[0].0);
        assert_eq!(
            OrientationID(0).select(selected.as_view()),
            Atom::Zero,
            "a stored root residue must not leak into another residue-map key"
        );
        assert_eq!(OrientationID(1).select(selected.as_view()), Atom::one());

        let explicit_localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let explicit = explicit_localizer.localized_orientation_terms(
            &graph,
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(1)],
            None,
            Some(OrientationID(1)),
            None,
            None,
        )?;
        assert_eq!(explicit, vec![(OrientationID(1), Atom::one())]);
        Ok(())
    }

    #[test]
    fn stored_generalized_root_map_preserves_its_own_selector() -> Result<()> {
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
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[],
            None,
            Some(OrientationID(2)),
            Some(&reduced.edge_energy_map),
            None,
        )?;

        assert_eq!(
            localized,
            vec![(OrientationID(2), Atom::one())],
            "a stored generalized root map keeps its original ID in the sparse selector sidecar",
        );

        let cut_valid_ids = BTreeSet::from([OrientationID(1), OrientationID(2)]);
        let cut_localized = localizer.localized_orientation_terms(
            &graph,
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[],
            Some(&cut_valid_ids),
            Some(OrientationID(2)),
            Some(&reduced.edge_energy_map),
            None,
        )?;
        assert_eq!(
            cut_localized,
            vec![(OrientationID(2), Atom::one())],
            "cut filtering must preserve a cut-valid stored root ID instead of re-hosting it",
        );
        assert_eq!(
            OrientationID(2).select(
                localizer
                    .residue_map_key_selector(cut_localized[0].0)
                    .as_view(),
            ),
            Atom::one(),
            "the branch's complete map key must select it without resolving its undirected edge",
        );

        let explicit_localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        assert_eq!(
            explicit_localizer.localized_orientation_terms(
                &graph,
                &graph.loop_momentum_basis,
                &reduced,
                &Atom::one(),
                &contract,
                &[],
                None,
                Some(OrientationID(2)),
                Some(&reduced.edge_energy_map),
                None,
            )?,
            vec![(OrientationID(2), Atom::one())],
            "an explicit orientation sum keeps the stored generalized branch exactly once without a selector",
        );
        assert_eq!(
            explicit_localizer.localized_orientation_terms(
                &graph,
                &graph.loop_momentum_basis,
                &reduced,
                &Atom::one(),
                &contract,
                &[],
                Some(&cut_valid_ids),
                Some(OrientationID(2)),
                Some(&reduced.edge_energy_map),
                None,
            )?,
            vec![(OrientationID(2), Atom::one())],
            "the explicit sum retains the stored cut-valid branch metadata without a selector",
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
            &graph.loop_momentum_basis,
            &reduced,
            &Atom::one(),
            &contract,
            &[EdgeIndex(0), EdgeIndex(1)],
            None,
            None,
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
