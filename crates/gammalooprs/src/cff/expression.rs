use std::collections::{BTreeMap, BTreeSet};

pub use three_dimensional_reps::expression::{
    AllOrientations, CFFVariant, OrientationData, OrientationExpression, OrientationID,
    OrientationSelector, RaisedEsurfaceData, RaisedEsurfaceDataView, RaisedEsurfaceGroup,
    RaisedEsurfaceGroupView, RaisedEsurfaceId,
};

use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    function,
    id::Replacement,
};

use super::{
    CutCFFIndex,
    esurface::{self, Esurface},
    hsurface::Hsurface,
    surface::{self, GammaLoopLinearEnergyExpr, LinearEnergyExpr},
};
use crate::{
    graph::{Graph, cuts::CutSet},
    utils::{GS, W_, ose_atom_from_index},
};

pub type ThreeDExpression<O> =
    three_dimensional_reps::expression::ThreeDExpression<O, Esurface, Hsurface>;
pub type CFFExpression<O> =
    three_dimensional_reps::expression::CFFExpression<O, Esurface, Hsurface>;

pub(crate) fn normalize_three_d_expression_cut_support_with_raised_edge_groups<O>(
    expression: &mut ThreeDExpression<O>,
    raised_edge_groups: &[Vec<EdgeIndex>],
) where
    O: From<usize> + Into<usize>,
{
    for orientation in expression.orientations.iter_mut() {
        for variant in &mut orientation.variants {
            variant.denominator_edges = normalize_cut_edge_support_with_raised_edge_groups(
                &variant.denominator_edges,
                raised_edge_groups,
            );
            variant.denominator_edge_support_signs = normalize_cut_edge_support_signs(
                std::mem::take(&mut variant.denominator_edge_support_signs),
                raised_edge_groups,
            );
        }
    }
}

fn normalize_cut_edge_support_signs(
    support_signs: BTreeMap<Vec<EdgeIndex>, i64>,
    raised_edge_groups: &[Vec<EdgeIndex>],
) -> BTreeMap<Vec<EdgeIndex>, i64> {
    support_signs
        .into_iter()
        .fold(BTreeMap::new(), |mut normalized, (support, sign)| {
            let support =
                normalize_cut_edge_support_with_raised_edge_groups(&support, raised_edge_groups);
            *normalized.entry(support).or_insert(1) *= sign;
            normalized
        })
}

pub(crate) fn normalize_cut_edge_support_with_raised_edge_groups(
    edges: &[EdgeIndex],
    raised_edge_groups: &[Vec<EdgeIndex>],
) -> Vec<EdgeIndex> {
    edges
        .iter()
        .map(|edge| {
            raised_edge_groups
                .iter()
                .find(|group| group.contains(edge))
                .and_then(|group| group.first())
                .copied()
                .unwrap_or(*edge)
        })
        .sorted()
        .dedup()
        .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum CutCffResidueAxis {
    RightThreshold,
    LeftThreshold,
    LuCut,
}

impl CutCffResidueAxis {
    fn set_order(self, index: &mut CutCFFIndex, order: usize) {
        match self {
            Self::RightThreshold => index.right_threshold_order = Some(order),
            Self::LeftThreshold => index.left_threshold_order = Some(order),
            Self::LuCut => index.lu_cut_order = Some(order),
        }
    }
}

fn apply_indexed_residue_selection<F>(
    residues: Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>,
    axis: CutCffResidueAxis,
    mut select: F,
) -> Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>
where
    F: FnMut(ThreeDExpression<OrientationID>) -> Vec<ThreeDExpression<OrientationID>>,
{
    residues
        .into_iter()
        .flat_map(|(index, expression)| {
            select(expression)
                .into_iter()
                .enumerate()
                .map(move |(i, residue)| {
                    let mut new_index = index;
                    axis.set_order(&mut new_index, i + 1);
                    // One causal coordinate has been consumed regardless of
                    // which raised-order branch `i` labels. Its energy-factor,
                    // contact, and residue signs remain owned by the selected
                    // CFF variant; the index records only the residue order.
                    (new_index, residue)
                })
        })
        .collect()
}

pub(super) fn select_indexed_cff_residues(
    mut cff: ThreeDExpression<OrientationID>,
    cutset: &CutSet,
    representation: three_dimensional_reps::RepresentationMode,
    raised_edge_groups: &[Vec<EdgeIndex>],
    mut bounds_provider: impl FnMut() -> Result<Vec<(usize, usize)>>,
) -> Result<Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>> {
    if let Some(lu) = cutset.residue_selector.lu.as_ref()
        && (lu.cut_edge_alternatives.is_empty()
            || lu.cut_edge_alternatives.iter().any(Vec::is_empty))
    {
        return Err(eyre::eyre!(
            "LU residue selection requires at least one non-empty physical cut alternative"
        ));
    }
    for group in cutset
        .residue_selector
        .right_th_cut
        .iter()
        .chain(&cutset.residue_selector.left_th_cut)
        .chain(cutset.residue_selector.lu.iter().map(|lu| &lu.raised_group))
    {
        if group.esurface_ids.is_empty() || group.max_occurence == 0 {
            return Err(eyre::eyre!(
                "CFF residue selection requires a nonempty surface group and a positive residue order"
            ));
        }
        if let Some(id) = group
            .esurface_ids
            .iter()
            .find(|id| cff.surfaces.esurface_cache.get(**id).is_none())
        {
            return Err(eyre::eyre!(
                "CFF residue selection references missing E-surface {id:?}"
            ));
        }
    }
    if representation == three_dimensional_reps::RepresentationMode::Ltd {
        use surface::{HybridSurfaceID, HybridSurfaceRef};
        // Raised physical edges share one OSE even when the common surface
        // catalogue names a different occurrence first. Remap energy values
        // only: numerator-map slots still name their original propagators.
        let aliases = raised_edge_groups
            .iter()
            .filter_map(|group| group.first().map(|first| (group, first)))
            .flat_map(|(group, first)| group.iter().map(move |edge| (edge.0, first.0)))
            .collect::<BTreeMap<_, _>>();
        let canonical_edge = |edge: &mut EdgeIndex| {
            if let Some(representative) = aliases.get(&edge.0) {
                *edge = EdgeIndex(*representative);
            }
        };
        for surface in &mut cff.surfaces.esurface_cache {
            *surface =
                Graph::normalize_esurface_with_raised_edge_groups(surface, raised_edge_groups);
        }
        for surface in &mut cff.surfaces.hsurface_cache {
            surface
                .positive_energies
                .iter_mut()
                .for_each(&canonical_edge);
            surface
                .negative_energies
                .iter_mut()
                .for_each(&canonical_edge);
        }
        for surface in &mut cff.surfaces.linear_surface_cache {
            surface.expression = surface.expression.clone().remap_internal_edges(&aliases);
        }
        for orientation in &mut cff.orientations {
            orientation.remap_on_shell_energy_values(&aliases);
        }
        let physical_surfaces = cff
            .surfaces
            .iter_all_surfaces()
            .filter_map(|(id, surface)| {
                let (positive, negative, shift) = match surface {
                    HybridSurfaceRef::Esurface(surface) => {
                        (&surface.energies, &[][..], &surface.external_shift)
                    }
                    HybridSurfaceRef::Hsurface(surface) => (
                        &surface.positive_energies,
                        surface.negative_energies.as_slice(),
                        &surface.external_shift,
                    ),
                    _ => return None,
                };
                let expression = positive
                    .iter()
                    .map(|edge| LinearEnergyExpr::ose(*edge, 1))
                    .chain(negative.iter().map(|edge| LinearEnergyExpr::ose(*edge, -1)))
                    .chain(
                        shift
                            .iter()
                            .map(|(edge, sign)| LinearEnergyExpr::external(*edge, *sign)),
                    )
                    .fold(LinearEnergyExpr::zero(), |sum, term| sum + term);
                Some((id, expression))
            })
            .collect::<BTreeMap<_, _>>();
        let lu_cut_alternatives = cutset.residue_selector.lu.as_ref().map(|lu| {
            lu.cut_edge_alternatives
                .iter()
                .map(|edges| {
                    normalize_cut_edge_support_with_raised_edge_groups(edges, raised_edge_groups)
                })
                .collect::<Vec<_>>()
        });
        let mut energy_degree_bounds: Option<Vec<(usize, usize)>> = None;
        let selections = [
            (
                CutCffResidueAxis::RightThreshold,
                cutset.residue_selector.right_th_cut.as_ref(),
            ),
            (
                CutCffResidueAxis::LeftThreshold,
                cutset.residue_selector.left_th_cut.as_ref(),
            ),
            (
                CutCffResidueAxis::LuCut,
                cutset
                    .residue_selector
                    .lu
                    .as_ref()
                    .map(|lu| &lu.raised_group),
            ),
        ];
        let mut canonical_selections = selections.map(|(axis, group)| {
            (
                axis,
                group.map(|group| {
                    let selected =
                        &physical_surfaces[&HybridSurfaceID::Esurface(group.esurface_ids[0])];
                    esurface::RaisedEsurfaceGroup {
                        esurface_ids: physical_surfaces
                            .iter()
                            .filter_map(|(id, expression)| match id {
                                HybridSurfaceID::Esurface(id) if expression == selected => {
                                    Some(*id)
                                }
                                _ => None,
                            })
                            .collect(),
                        max_occurence: group.max_occurence,
                    }
                }),
            )
        });
        for (_, group) in &mut canonical_selections {
            let Some(group) = group else { continue };
            cff.normalize_single_raising(group);
            group.max_occurence = group.max_occurence.max(
                cff.orientations
                    .iter()
                    .map(|orientation| {
                        orientation.max_effective_denominator_value_count_on_branch(
                            &HybridSurfaceID::Esurface(group.esurface_ids[0]),
                        )
                    })
                    .max()
                    .unwrap_or_default(),
            );
        }
        // Exact denominator sources can carry higher physical poles than the
        // parent graph's CutSet. Preserve the source's effective order while
        // leaving the authoritative CFF selection contract unchanged.
        let selections = canonical_selections
            .each_ref()
            .map(|(axis, group)| (*axis, group.as_ref()));
        let mut residues = vec![(CutCFFIndex::new_all_none(), cff)];
        for (axis, group) in selections {
            let Some(group) = group else { continue };
            let selected =
                physical_surfaces[&HybridSurfaceID::Esurface(group.esurface_ids[0])].clone();
            let other_selected = selections
                .iter()
                .filter(|(other_axis, _)| *other_axis != axis)
                .filter_map(|(_, group)| {
                    group.map(|group| {
                        physical_surfaces[&HybridSurfaceID::Esurface(group.esurface_ids[0])].clone()
                    })
                })
                .collect::<Vec<_>>();
            let mut next = Vec::new();
            for (index, mut expression) in residues {
                if let CutCffResidueAxis::LuCut = axis {
                    expression = expression
                        .restrict_to_cut_alternatives(lu_cut_alternatives.as_ref().unwrap());
                }
                let mut surface_expressions = physical_surfaces.clone();
                surface_expressions.extend(
                    expression
                        .surfaces
                        .linear_surface_cache
                        .iter_enumerated()
                        .map(|(id, surface)| {
                            (HybridSurfaceID::Linear(id), surface.expression.clone())
                        }),
                );
                crate::debug_tags!(#cff, #trace;
                    stage = "ltd_pole_selection_input",
                    axis = ?axis,
                    cut_index = ?index,
                    selected_surface = ?selected,
                    scalar_surfaces = ?surface_expressions,
                    residue_maps = ?expression.orientations,
                    "Selecting an LTD physical pole"
                );
                let coefficients = match expression.select_ltd_residue(
                    &surface_expressions,
                    &selected,
                    &other_selected,
                    energy_degree_bounds.as_deref()) {
                    Ok(coefficients) => coefficients,
                    Err(three_dimensional_reps::generation::GenerationError::LtdResidueRequiresEnergyBounds) => {
                        let bounds = bounds_provider()?;
                        let coefficients = expression.select_ltd_residue(&surface_expressions, &selected, &other_selected, Some(&bounds))?;
                        energy_degree_bounds = Some(bounds);
                        coefficients
                    }
                    Err(error) => return Err(error.into()),
                };
                for (order, coefficient) in coefficients.into_iter().enumerate() {
                    let mut index = index;
                    axis.set_order(&mut index, order + 1);
                    crate::debug_tags!(#cff, #trace;
                        stage = "ltd_pole_selection_coefficient",
                        axis = ?axis,
                        cut_index = ?index,
                        scalar_surfaces = ?coefficient.surfaces.linear_surface_cache,
                        residue_maps = ?coefficient.orientations,
                        "Selected LTD physical-pole coefficient"
                    );
                    next.push((index, coefficient));
                }
            }
            residues = next;
        }
        return Ok(residues);
    }
    // Carry selected causal coordinates through the operation pipeline rather
    // than reconstructing them later from raised-order labels in
    // `CutCFFIndex`. Energy-factor ownership stays with the generated CFF
    // variants and is deliberately independent of this indexing operation.
    let mut residues = vec![(CutCFFIndex::new_all_none(), cff)];

    if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
        residues = apply_indexed_residue_selection(
            residues,
            CutCffResidueAxis::RightThreshold,
            |expression| expression.select_esurface_residue(right_threshold, representation),
        );
    }

    if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
        residues = apply_indexed_residue_selection(
            residues,
            CutCffResidueAxis::LeftThreshold,
            |expression| expression.select_esurface_residue(left_threshold, representation),
        );
    }

    if let Some(lu) = cutset.residue_selector.lu.as_ref() {
        residues =
            apply_indexed_residue_selection(residues, CutCffResidueAxis::LuCut, |expression| {
                // Physical support selection is independent of the 3D
                // representation. Each representation owns the conversion
                // from its signed denominator to the selected pole.
                expression
                    .restrict_to_cut_alternatives(&lu.cut_edge_alternatives)
                    .select_esurface_residue(&lu.raised_group, representation)
            });
    }

    Ok(residues)
}

pub trait GammaLoopCFFVariant {
    fn to_atom_gs(&self) -> Atom;
}

impl GammaLoopCFFVariant for CFFVariant {
    fn to_atom_gs(&self) -> Atom {
        let half_edge_factor = self
            .half_edges
            .iter()
            .map(|edge_id| Atom::num(1) / (Atom::num(2) * ose_atom_from_index(*edge_id)))
            .reduce(|acc, factor| acc * factor)
            .unwrap_or_else(|| Atom::num(1));
        let scale_factor = if self.uniform_scale_power == 0 {
            Atom::num(1)
        } else {
            Atom::num(1)
                / Atom::var(GS.numerator_sampling_scale).pow(self.uniform_scale_power as i64)
        };

        let numerator_surface_factor = self
            .numerator_surfaces
            .iter()
            .map(|surface_id| Atom::from(*surface_id))
            .reduce(|acc, factor| acc * factor)
            .unwrap_or_else(|| Atom::num(1));

        Atom::num(self.prefactor.clone())
            * half_edge_factor
            * scale_factor
            * numerator_surface_factor
            * self.denominator.to_atom_inv()
    }
}

pub trait GammaLoopOrientationExpression {
    fn to_atom_gs(&self) -> Atom;
    fn energy_replacements_gs(&self, graph: &Graph) -> Vec<Replacement>;
}

impl GammaLoopOrientationExpression for OrientationExpression {
    fn to_atom_gs(&self) -> Atom {
        self.variants
            .iter()
            .map(GammaLoopCFFVariant::to_atom_gs)
            .reduce(|acc, atom| acc + atom)
            .unwrap_or_else(Atom::new)
    }

    fn energy_replacements_gs(&self, graph: &Graph) -> Vec<Replacement> {
        energy_map_replacements_gs(
            self.edge_energy_map
                .iter()
                .map(|energy| energy.to_atom_gs(&[])),
            graph,
        )
    }
}

pub(crate) fn energy_map_replacements_gs(
    edge_energy_map: impl IntoIterator<Item = Atom>,
    graph: &Graph,
) -> Vec<Replacement> {
    let edge_energy_map = edge_energy_map.into_iter().collect::<Vec<_>>();
    let mut replacements = Vec::new();
    let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);
    let external_edges = graph
        .underlying
        .iter_edges()
        .filter_map(|(pair, edge_id, _)| (!pair.is_paired()).then_some(edge_id))
        .collect::<BTreeSet<_>>();

    for (edge_id, energy) in edge_energy_map.iter().enumerate() {
        let edge_id = EdgeIndex(edge_id);
        // Remapping a generated internal-edge map into GammaLoop's physical
        // namespace pads unpaired external-edge slots with zero. Those slots
        // are not residue samples: their EMR momenta must remain the original
        // external momenta carried by the factorized numerator.
        if external_edges.contains(&edge_id) {
            continue;
        }
        replacements.push(Replacement::new(
            GS.emr_mom(edge_id, AIND_SYMBOLS.cind.call(Atom::Zero))
                .to_pattern(),
            energy.clone().to_pattern(),
        ));
        replacements.push(Replacement::new(
            GS.emr_mom(edge_id, &mink_index).to_pattern(),
            (GS.emr_vec_index(edge_id, &mink_index) + energy * GS.energy_delta(&mink_index))
                .to_pattern(),
        ));
    }

    for (loop_id, loop_edge_id) in graph.loop_momentum_basis.loop_edges.iter_enumerated() {
        let loop_id = usize::from(loop_id);
        let loop_id_atom = Atom::num(loop_id as i64);
        let energy = edge_energy_map
            .get(usize::from(*loop_edge_id))
            .cloned()
            .unwrap_or_else(Atom::new);
        replacements.push(Replacement::new(
            function!(
                GS.loop_mom,
                loop_id_atom.clone(),
                AIND_SYMBOLS.cind.call(Atom::Zero)
            )
            .to_pattern(),
            energy.clone().to_pattern(),
        ));
        for spatial_index in 1..=3 {
            replacements.push(Replacement::new(
                function!(
                    GS.loop_mom,
                    loop_id_atom.clone(),
                    AIND_SYMBOLS.cind.call(spatial_index)
                )
                .to_pattern(),
                GS.emr_mom(*loop_edge_id, AIND_SYMBOLS.cind.call(spatial_index))
                    .to_pattern(),
            ));
        }
        replacements.push(Replacement::new(
            FunctionBuilder::new(GS.loop_mom)
                .add_arg(loop_id as i64)
                .add_arg(mink_index.as_view())
                .finish()
                .to_pattern(),
            (GS.emr_vec_index(*loop_edge_id, &mink_index) + energy * GS.energy_delta(&mink_index))
                .to_pattern(),
        ));
    }

    replacements
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::{EdgeVec, Orientation};

    use super::*;
    use crate::{
        cff::surface::LinearEnergyExpr, dot, graph::parse::from_dot::IntoGraph,
        initialisation::test_initialise, utils::external_energy_atom_from_index,
    };

    #[test]
    fn affine_energy_map_keeps_two_numerator_factors_under_one_map() -> color_eyre::Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph affine_map {
                edge [num=1 mass=0]
                node [num=1]
                A -> B [id=0]
                A -> B [id=1]
            }
        )?;
        let energy_map = LinearEnergyExpr {
            internal_terms: vec![(EdgeIndex(0), 2.into()), (EdgeIndex(1), (-3).into())],
            external_terms: vec![(EdgeIndex(2), 5.into())],
            uniform_scale_coeff: 7.into(),
            constant: 11.into(),
        }
        .canonical();
        let orientation = OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([
                Orientation::Undirected,
                Orientation::Undirected,
            ])),
            loop_energy_map: Vec::new(),
            edge_energy_map: vec![energy_map.clone(), LinearEnergyExpr::zero()],
            variants: Vec::new(),
        };
        let factor_a = Atom::var(symbolica::symbol!("factor_a"));
        let factor_b = Atom::var(symbolica::symbol!("factor_b"));
        let energy_component = GS.emr_mom(EdgeIndex(0), AIND_SYMBOLS.cind.call(Atom::Zero));
        let numerator =
            (energy_component.clone() + factor_a.clone()) * (energy_component + factor_b.clone());

        let mapped = numerator.replace_multiple(orientation.energy_replacements_gs(&graph));
        let mapped_energy = Atom::num(2) * ose_atom_from_index(EdgeIndex(0))
            - Atom::num(3) * ose_atom_from_index(EdgeIndex(1))
            + Atom::num(5) * external_energy_atom_from_index(EdgeIndex(2))
            + Atom::num(7) * Atom::var(GS.numerator_sampling_scale)
            + Atom::num(11);
        let expected =
            (mapped_energy.clone() + factor_a.clone()) * (mapped_energy.clone() + factor_b.clone());
        assert_eq!(mapped, expected);
        let mixed = (mapped_energy.clone() + factor_a) * (mapped_energy + Atom::num(1) + factor_b);
        assert_ne!(mapped, mixed);
        Ok(())
    }

    #[test]
    fn padded_external_energy_slots_do_not_erase_factorized_external_momenta()
    -> color_eyre::Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph padded_external_energy_slots {
            edge [num=1 mass=0]
            node [num=1]
            ext [style=invis]
            v0;
            v1;
            v2;
            ext -> v0 [id=0]
            ext -> v1 [id=1]
            v2 -> ext [id=2]
            v0 -> v2 [id=3 lmb_id=0]
            v1 -> v0 [id=4]
            v2 -> v1 [id=5]
        })?;
        let mapped_internal_energy = LinearEnergyExpr {
            internal_terms: vec![(EdgeIndex(3), 2.into())],
            external_terms: vec![(EdgeIndex(0), (-1).into())],
            uniform_scale_coeff: 0.into(),
            constant: 0.into(),
        }
        .canonical();
        let mut edge_energy_map = vec![LinearEnergyExpr::zero(); graph.underlying.n_edges()];
        edge_energy_map[3] = mapped_internal_energy.clone();
        let external_energy = GS.emr_mom(EdgeIndex(0), AIND_SYMBOLS.cind.call(Atom::Zero));
        let internal_energy = GS.emr_mom(EdgeIndex(3), AIND_SYMBOLS.cind.call(Atom::Zero));
        let spectator = Atom::var(symbolica::symbol!("spectator"));
        let numerator = external_energy.clone() * (internal_energy + spectator.clone());

        let mapped = numerator.replace_multiple(energy_map_replacements_gs(
            edge_energy_map.iter().map(|energy| energy.to_atom_gs(&[])),
            &graph,
        ));
        assert_eq!(
            mapped,
            external_energy * (mapped_internal_energy.to_atom_gs(&[]) + spectator)
        );
        Ok(())
    }

    #[test]
    fn raised_alias_normalization_merges_cut_support_signs() {
        let support_signs = BTreeMap::from([
            (vec![EdgeIndex(2), EdgeIndex(3)], -1),
            (vec![EdgeIndex(3), EdgeIndex(4)], -1),
        ]);

        assert_eq!(
            normalize_cut_edge_support_signs(support_signs, &[vec![EdgeIndex(2), EdgeIndex(4)]],),
            BTreeMap::from([(vec![EdgeIndex(2), EdgeIndex(3)], 1)])
        );
    }
}
