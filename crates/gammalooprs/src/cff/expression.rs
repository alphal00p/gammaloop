pub use three_dimensional_reps::expression::{
    AllOrientations, CFFVariant, GraphOrientation, OrientationData, OrientationExpression,
    OrientationID, OrientationSelector, RaisedEsurfaceData, RaisedEsurfaceDataView,
    RaisedEsurfaceGroup, RaisedEsurfaceGroupView, RaisedEsurfaceId, ResidualDenominator,
};

use color_eyre::Result;
use eyre::eyre;
use itertools::{EitherOrBoth, Itertools};
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, HedgePair, Orientation};
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder},
    function,
    id::Replacement,
    symbol,
};
use three_dimensional_reps::tree::{NodeId, Tree};

use super::{
    esurface::{Esurface, RaisedEsurfaceGroup as GammaLoopRaisedEsurfaceGroup},
    hsurface::Hsurface,
    surface::{
        EsurfaceID, GammaLoopLinearEnergyExpr, GammaLoopSurfaceCache, HybridSurfaceID,
        LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind, SurfaceOrigin,
    },
};
use crate::{
    graph::{
        Graph,
        cuts::{CutSet, LuLocalSeriesCoordinate, LuResidueSelectionBasis},
    },
    settings::global::{OrientationPattern, ThreeDRepresentation},
    utils::{
        GS, W_, external_energy_atom_from_index, ose_atom_from_index, symbolica_ext::CallSymbol,
    },
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
        }
    }
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

pub(crate) fn remove_ltd_global_contact_completions_from_local_residue(
    residue: &mut ThreeDExpression<OrientationID>,
) {
    // Repeated-signature LTD generation may add CFF lower-sector contact
    // completions so that the standalone, globally evaluated 3D expression is
    // equivalent for higher energy-degree numerators. These completions are not
    // local LTD residues: confluent LTD residues already carry the local
    // numerator derivatives needed for LU and threshold residues. Projecting the
    // global CFF completion again into a local residue double counts the local
    // pole.
    for orientation in &mut residue.orientations {
        orientation.variants.retain(|variant| {
            !variant.origin.as_deref().is_some_and(|origin| {
                origin.starts_with("bounded_degree_quadratic_recursive_contact")
                    || origin.starts_with("bounded_degree_known_factor_cff_contact")
                    || origin.starts_with("bounded_degree_ltd_finite_pole_contact")
            })
        });
    }
}

pub(crate) fn retain_three_d_variants_supporting_any_cut_edges(
    expression: &mut ThreeDExpression<OrientationID>,
    cut_edge_sets: &[Vec<EdgeIndex>],
) {
    if cut_edge_sets.is_empty() {
        return;
    }

    for orientation in &mut expression.orientations {
        orientation.variants.retain(|variant| {
            cut_edge_sets.iter().any(|cut_edges| {
                cut_edges
                    .iter()
                    .all(|cut_edge| variant.denominator_edges.contains(cut_edge))
            })
        });
    }
}

pub(crate) fn prepare_ltd_lu_local_series_expression(
    expression: &mut ThreeDExpression<OrientationID>,
    cutset: &CutSet,
) {
    remove_ltd_global_contact_completions_from_local_residue(expression);
    retain_three_d_variants_supporting_any_cut_edges(
        expression,
        cutset.residue_selector.lu_cut_edge_sets(),
    );
}

pub(crate) fn cutset_has_repeated_lu_pole(cutset: &CutSet) -> bool {
    cutset
        .residue_selector
        .lu_cut()
        .is_some_and(|lu_cut| lu_cut.max_occurence > 1)
}

pub(crate) fn cff_expression_uses_ltd_lu_residue_basis(
    expression: &ThreeDExpression<OrientationID>,
) -> bool {
    expression
        .orientations
        .iter()
        .flat_map(|orientation| &orientation.variants)
        .any(|variant| {
            variant
                .origin
                .as_deref()
                .is_some_and(|origin| origin.starts_with("ltd_confluent"))
        })
}

pub(crate) fn ltd_lu_local_series_prefactor_sign(cutset: &CutSet) -> i64 {
    cutset
        .residue_selector
        .ltd_local_series_residue_prefactor_sign()
}

pub(crate) fn ltd_lu_local_series_coefficients_from_parametric_atom(
    graph_name: &str,
    expression: &ThreeDExpression<OrientationID>,
    atom: Atom,
    lu_cut: &GammaLoopRaisedEsurfaceGroup,
    selected_esurface_signs: &[(EsurfaceID, i64)],
    local_series_coordinates: &[LuLocalSeriesCoordinate],
    _allow_orientation_dependent_surface_family_sign: bool,
) -> Result<Vec<Atom>> {
    ltd_lu_local_series_coefficients(
        graph_name,
        expression,
        lu_cut,
        selected_esurface_signs,
        local_series_coordinates,
        |selected_context| {
            let surface_family_sign = ltd_lu_local_series_surface_family_sign_atom(
                graph_name,
                expression,
                lu_cut,
                selected_context.esurface_id,
                selected_context.localized_external_edge,
                selected_context.localized_external_coeff,
            )?;
            Ok(atom.clone() * surface_family_sign)
        },
    )
}

struct LtdLuLocalSeriesSelectedSurface {
    esurface_id: EsurfaceID,
    localized_external_edge: EdgeIndex,
    localized_external_coeff: i64,
    surface_without_localized_external: Atom,
}

fn ltd_lu_local_series_coefficients(
    graph_name: &str,
    expression: &ThreeDExpression<OrientationID>,
    lu_cut: &GammaLoopRaisedEsurfaceGroup,
    selected_esurface_signs: &[(EsurfaceID, i64)],
    local_series_coordinates: &[LuLocalSeriesCoordinate],
    mut parametric_atom_for_selected_surface: impl FnMut(
        &LtdLuLocalSeriesSelectedSurface,
    ) -> Result<Atom>,
) -> Result<Vec<Atom>> {
    let residue_parameter = symbol!("ltd_lu_residue_parameter");
    let residue_parameter_atom = Atom::var(residue_parameter);
    let max_occurrence = lu_cut.max_occurence;
    let mut coefficients = vec![Atom::Zero; max_occurrence];

    let selected_surfaces = if max_occurrence == 1 {
        selected_esurface_signs.to_vec()
    } else {
        lu_cut
            .esurface_ids
            .iter()
            .copied()
            .map(|esurface_id| (esurface_id, 1))
            .collect_vec()
    };

    for (esurface_id, selected_sign) in &selected_surfaces {
        let esurface = expression.surfaces.esurface_cache[*esurface_id].clone();
        let coordinate = local_series_coordinates
            .iter()
            .find(|coordinate| coordinate.esurface_id == *esurface_id)
            .ok_or_else(|| {
                eyre!(
                    "cannot take an LTD LU residue for graph {graph_name} on E-surface {esurface_id:?}: no local-series coordinate was recorded"
                )
            })?;
        let localized_external_edge = coordinate.external_edge;
        let localized_external_coeff = esurface
            .external_shift
            .iter()
            .filter_map(|(edge_id, sign)| (*edge_id == localized_external_edge).then_some(*sign))
            .sum::<i64>();
        if localized_external_coeff == 0 {
            return Err(eyre!(
                "cannot take an LTD LU residue for graph {graph_name} on E-surface {esurface_id:?}: local-series external coordinate {} is not part of the generated E-surface external shift {:?}",
                usize::from(localized_external_edge),
                esurface.external_shift,
            ));
        }

        let surface_without_localized_external_expr =
            surface_without_localized_external_expr(&esurface, localized_external_edge);
        let surface_without_localized_external =
            surface_without_localized_external_expr.to_atom_gs(&[]);
        let selected_context = LtdLuLocalSeriesSelectedSurface {
            esurface_id: *esurface_id,
            localized_external_edge,
            localized_external_coeff,
            surface_without_localized_external,
        };
        let localized_external_replacement = (Atom::num(*selected_sign)
            * residue_parameter_atom.clone()
            - &selected_context.surface_without_localized_external)
            / Atom::num(localized_external_coeff);
        let localized_atom = parametric_atom_for_selected_surface(&selected_context)?
            .replace(external_energy_atom_from_index(localized_external_edge))
            .with(localized_external_replacement)
            .expand();
        let series = localized_atom
            .series(residue_parameter, Atom::Zero, 0)
            .map_err(|error| {
                eyre!(
                    "failed to extract LTD LU residue for graph {graph_name} on E-surface {esurface_id:?}: {error}"
                )
            })?;
        for occurrence in 1..=max_occurrence {
            let occurrence_i64 = i64::try_from(occurrence).map_err(|_| {
                eyre!(
                    "failed to extract LTD LU residue for graph {graph_name} on E-surface {esurface_id:?}: occurrence does not fit in i64"
                )
            })?;
            coefficients[occurrence - 1] += series.coefficient((-occurrence_i64, 1).into());
        }
    }

    Ok(coefficients)
}

fn surface_without_localized_external_expr(
    esurface: &Esurface,
    localized_external_edge: EdgeIndex,
) -> LinearEnergyExpr {
    let surface_without_localized_external = esurface
        .energies
        .iter()
        .fold(LinearEnergyExpr::zero(), |acc, edge_id| {
            acc + LinearEnergyExpr::ose(*edge_id, 1)
        });
    esurface.external_shift.iter().fold(
        surface_without_localized_external,
        |acc, (external_edge, external_sign)| {
            if *external_edge == localized_external_edge {
                acc
            } else {
                acc + LinearEnergyExpr::external(*external_edge, *external_sign)
            }
        },
    )
}

fn ltd_lu_local_series_surface_family_sign_atom(
    graph_name: &str,
    expression: &ThreeDExpression<OrientationID>,
    lu_cut: &GammaLoopRaisedEsurfaceGroup,
    selected_esurface_id: EsurfaceID,
    localized_external_edge: EdgeIndex,
    localized_external_coeff: i64,
) -> Result<Atom> {
    let selected_surface_ids = lu_cut.esurface_ids.iter().copied().collect_vec();
    let replacement = selected_lu_surface_zero_replacement(
        expression,
        selected_esurface_id,
        localized_external_edge,
        localized_external_coeff,
    )?;

    let mut global_signs = std::collections::BTreeSet::new();
    for orientation in &expression.orientations {
        let mut signs = std::collections::BTreeSet::new();
        for variant in &orientation.variants {
            let selected_numerator_count = variant
                .numerator_surfaces
                .iter()
                .filter(|surface_id| {
                    matches!(
                        surface_id,
                        HybridSurfaceID::Esurface(esurface_id)
                            if selected_surface_ids.contains(esurface_id)
                    )
                })
                .count();

            for chain in denominator_tree_chains_gs(&variant.denominator) {
                let selected_denominator_count = chain
                    .iter()
                    .filter(|surface_id| {
                        matches!(
                            surface_id,
                            HybridSurfaceID::Esurface(esurface_id)
                                if selected_surface_ids.contains(esurface_id)
                        )
                    })
                    .count();
                if selected_denominator_count <= selected_numerator_count {
                    continue;
                }
                signs.insert(localized_variant_surface_family_sign(
                    graph_name,
                    expression,
                    &chain,
                    &variant.numerator_surfaces,
                    &selected_surface_ids,
                    selected_esurface_id,
                    localized_external_edge,
                    &replacement,
                )?);
            }
        }
        global_signs.insert(unique_or_unfactorable_surface_family_sign(&signs));
    }

    Ok(Atom::num(unique_or_unfactorable_surface_family_sign(
        &global_signs,
    )))
}

fn unique_or_unfactorable_surface_family_sign(signs: &std::collections::BTreeSet<i64>) -> i64 {
    if signs.len() == 1 {
        *signs.iter().next().expect("one sign was checked above")
    } else {
        // Branch-dependent signs cannot be pulled in front of the branch-summed
        // local atom. They require per-branch Laurent extraction, so keep the
        // exact atom-level branch signs unchanged at this level.
        1
    }
}

fn selected_lu_surface_zero_replacement(
    expression: &ThreeDExpression<OrientationID>,
    esurface_id: EsurfaceID,
    localized_external_edge: EdgeIndex,
    localized_external_coeff: i64,
) -> Result<LinearEnergyExpr> {
    if localized_external_coeff.unsigned_abs() != 1 {
        return Err(eyre!(
            "cannot build local zero-surface replacement for E-surface {esurface_id:?}: expected a unit localized external coefficient, found {localized_external_coeff}"
        ));
    }
    let esurface = &expression.surfaces.esurface_cache[esurface_id];
    let surface_without_localized_external =
        surface_without_localized_external_expr(esurface, localized_external_edge);
    Ok(surface_without_localized_external.scale(-localized_external_coeff.signum()))
}

fn localized_variant_surface_family_sign(
    graph_name: &str,
    expression: &ThreeDExpression<OrientationID>,
    denominator_chain: &[HybridSurfaceID],
    numerator_surfaces: &[HybridSurfaceID],
    selected_surface_ids: &[EsurfaceID],
    selected_esurface_id: EsurfaceID,
    localized_external_edge: EdgeIndex,
    replacement: &LinearEnergyExpr,
) -> Result<i64> {
    let denominator_sign = denominator_chain
        .iter()
        .filter(|surface_id| !is_selected_lu_surface(**surface_id, selected_surface_ids))
        .map(|surface_id| {
            localized_surface_family_sign(
                graph_name,
                expression,
                *surface_id,
                selected_esurface_id,
                localized_external_edge,
                replacement,
            )
        })
        .product::<Result<i64>>()?;
    let numerator_sign = numerator_surfaces
        .iter()
        .filter(|surface_id| !is_selected_lu_surface(**surface_id, selected_surface_ids))
        .map(|surface_id| {
            localized_surface_family_sign(
                graph_name,
                expression,
                *surface_id,
                selected_esurface_id,
                localized_external_edge,
                replacement,
            )
        })
        .product::<Result<i64>>()?;
    Ok(denominator_sign * numerator_sign)
}

fn is_selected_lu_surface(
    surface_id: HybridSurfaceID,
    selected_surface_ids: &[EsurfaceID],
) -> bool {
    matches!(
        surface_id,
        HybridSurfaceID::Esurface(esurface_id) if selected_surface_ids.contains(&esurface_id)
    )
}

fn localized_surface_family_sign(
    graph_name: &str,
    expression: &ThreeDExpression<OrientationID>,
    surface_id: HybridSurfaceID,
    selected_esurface_id: EsurfaceID,
    localized_external_edge: EdgeIndex,
    replacement: &LinearEnergyExpr,
) -> Result<i64> {
    if surface_id == HybridSurfaceID::Unit {
        return Ok(1);
    }
    let Some(expression_for_surface) = residue_surface_linear_expression(expression, surface_id)
    else {
        return Ok(1);
    };
    let localized = expression_for_surface
        .substitute_external_energy(localized_external_edge, replacement)
        .canonical();
    if localized.is_zero() {
        return Err(eyre!(
            "localizing graph {graph_name} on LU E-surface {selected_esurface_id:?} made spectator surface {surface_id:?} vanish"
        ));
    }
    Ok(find_matching_localized_surface(
        expression,
        &localized,
        localized_external_edge,
        replacement,
    )
    .map(|(_, sign)| sign)
    .unwrap_or(1))
}

fn denominator_tree_chains_gs(tree: &Tree<HybridSurfaceID>) -> Vec<Vec<HybridSurfaceID>> {
    fn walk(
        tree: &Tree<HybridSurfaceID>,
        node_id: NodeId,
        current: &mut Vec<HybridSurfaceID>,
        out: &mut Vec<Vec<HybridSurfaceID>>,
    ) {
        let node = tree.get_node(node_id);
        if node.data != HybridSurfaceID::Unit {
            current.push(node.data);
        }
        if node.children.is_empty() {
            out.push(current.clone());
        } else {
            for child in &node.children {
                walk(tree, *child, current, out);
            }
        }
        if node.data != HybridSurfaceID::Unit {
            current.pop();
        }
    }

    let mut out = Vec::new();
    walk(tree, NodeId::root(), &mut Vec::new(), &mut out);
    out
}

pub(crate) fn select_lu_cut_residue_for_representation<O>(
    expression: ThreeDExpression<O>,
    lu_cut: &GammaLoopRaisedEsurfaceGroup,
    cut_edge_sets: &[Vec<EdgeIndex>],
    ltd_lu_cut_esurface_signs: &[(EsurfaceID, i64)],
    representation: ThreeDRepresentation,
) -> Vec<ThreeDExpression<O>>
where
    O: Clone + From<usize> + Into<usize>,
    usize: From<O>,
{
    let selection_basis = match representation {
        ThreeDRepresentation::Cff => LuResidueSelectionBasis::PositiveEnergyCutkosky,
        ThreeDRepresentation::Ltd => LuResidueSelectionBasis::GeneratedEsurface,
    };
    select_lu_cut_residue_for_basis(
        expression,
        lu_cut,
        cut_edge_sets,
        ltd_lu_cut_esurface_signs,
        selection_basis,
    )
}

pub(crate) fn select_lu_cut_residue_for_basis<O>(
    expression: ThreeDExpression<O>,
    lu_cut: &GammaLoopRaisedEsurfaceGroup,
    cut_edge_sets: &[Vec<EdgeIndex>],
    ltd_lu_cut_esurface_signs: &[(EsurfaceID, i64)],
    selection_basis: LuResidueSelectionBasis,
) -> Vec<ThreeDExpression<O>>
where
    O: Clone + From<usize> + Into<usize>,
    usize: From<O>,
{
    match selection_basis {
        LuResidueSelectionBasis::PositiveEnergyCutkosky => {
            expression.select_esurface_residue_with_cut_edges(lu_cut, cut_edge_sets)
        }
        LuResidueSelectionBasis::GeneratedEsurface => {
            // LTD and confluent-generated CFF expressions expose LU residues in
            // the generated E-surface variables. Keep the selected-denominator
            // orientation signs with that basis; pure CFF expressions are
            // already assembled in the positive-energy Cutkosky convention.
            expression.select_esurface_residue_with_cut_edges_and_esurface_signs(
                lu_cut,
                cut_edge_sets,
                ltd_lu_cut_esurface_signs,
            )
        }
    }
}

pub(crate) fn localize_three_d_expression_on_esurface(
    residue: &mut ThreeDExpression<OrientationID>,
    esurface_id: EsurfaceID,
) -> Result<()> {
    let esurface = residue.surfaces.esurface_cache[esurface_id].clone();
    let Some((localized_external_edge, localized_external_sign)) =
        esurface.external_shift.first().copied()
    else {
        return Ok(());
    };

    if localized_external_sign.unsigned_abs() != 1 {
        return Err(eyre!(
            "cannot localize 3D expression on E-surface {esurface_id:?}: expected a unit external-energy coefficient, found {localized_external_sign}"
        ));
    }

    let replacement = esurface
        .energies
        .iter()
        .fold(LinearEnergyExpr::zero(), |acc, edge_id| {
            acc + LinearEnergyExpr::ose(*edge_id, 1)
        });
    let replacement = esurface
        .external_shift
        .iter()
        .fold(replacement, |acc, (external_edge, external_sign)| {
            if *external_edge == localized_external_edge {
                acc
            } else {
                acc + LinearEnergyExpr::external(*external_edge, *external_sign)
            }
        })
        .scale(-localized_external_sign);

    for orientation in &mut residue.orientations {
        for energy_expr in &mut orientation.loop_energy_map {
            *energy_expr = energy_expr
                .clone()
                .substitute_external_energy(localized_external_edge, &replacement);
        }
        for energy_expr in &mut orientation.edge_energy_map {
            *energy_expr = energy_expr
                .clone()
                .substitute_external_energy(localized_external_edge, &replacement);
        }
    }

    let referenced_surfaces = referenced_residue_surfaces(residue);
    let mut surface_map = std::collections::BTreeMap::new();
    for (surface_id, expression) in residue_surface_linear_expressions(residue) {
        if !referenced_surfaces.contains(&surface_id) {
            continue;
        }
        let localized = expression
            .substitute_external_energy(localized_external_edge, &replacement)
            .canonical();
        if localized.is_zero() {
            return Err(eyre!(
                "localizing on E-surface {esurface_id:?} made spectator surface {surface_id:?} vanish"
            ));
        }
        // LTD threshold localization can turn a subsequent LU E-surface into
        // an already known surface with the opposite convention, often an
        // H-surface. Canonicalize these sign-equivalent surfaces here and carry
        // the sign explicitly through the residue data.
        let localized_id = find_matching_localized_surface(
            residue,
            &localized,
            localized_external_edge,
            &replacement,
        )
        .unwrap_or_else(|| {
            (
                HybridSurfaceID::Linear(intern_localized_linear_surface(residue, localized)),
                1,
            )
        });
        surface_map.insert(surface_id, localized_id);
    }

    for orientation in &mut residue.orientations {
        for variant in &mut orientation.variants {
            let numerator_surface_sign = variant
                .numerator_surfaces
                .iter()
                .filter_map(|surface_id| surface_map.get(surface_id).map(|(_, sign)| *sign))
                .product::<i64>();
            for surface_id in &mut variant.numerator_surfaces {
                if let Some((localized_id, _)) = surface_map.get(surface_id) {
                    *surface_id = *localized_id;
                }
            }
            let localized_denominator_signs = variant
                .denominator
                .iter_nodes()
                .filter_map(|node| surface_map.get(&node.data).copied())
                .collect::<Vec<_>>();
            let denominator_surface_sign = localized_denominator_signs
                .iter()
                .map(|(_, sign)| *sign)
                .product::<i64>();
            variant.denominator.map_mut(|surface_id| {
                if let Some((localized_id, _)) = surface_map.get(surface_id) {
                    *surface_id = *localized_id;
                }
            });
            let localization_sign = numerator_surface_sign * denominator_surface_sign;
            if localization_sign < 0 {
                variant.prefactor *= Atom::num(localization_sign);
            }
            let mut denominator_surface_signs = std::collections::BTreeMap::new();
            for (surface_id, sign) in std::mem::take(&mut variant.denominator_surface_signs) {
                let (localized_id, localization_sign) = surface_map
                    .get(&surface_id)
                    .copied()
                    .unwrap_or((surface_id, 1));
                let total_sign = sign * localization_sign;
                denominator_surface_signs
                    .entry(localized_id)
                    .and_modify(|existing_sign| *existing_sign *= total_sign)
                    .or_insert(total_sign);
            }
            for (localized_id, localization_sign) in localized_denominator_signs {
                if localization_sign < 0 {
                    denominator_surface_signs
                        .entry(localized_id)
                        .and_modify(|existing_sign| *existing_sign *= localization_sign)
                        .or_insert(localization_sign);
                }
            }
            denominator_surface_signs.retain(|_, sign| *sign < 0);
            variant.denominator_surface_signs = denominator_surface_signs;
        }
    }

    Ok(())
}

fn find_matching_localized_surface(
    residue: &ThreeDExpression<OrientationID>,
    localized: &LinearEnergyExpr,
    localized_external_edge: EdgeIndex,
    replacement: &LinearEnergyExpr,
) -> Option<(HybridSurfaceID, i64)> {
    let negative_localized = localized.clone().scale(-1);
    for (id, surface) in residue.surfaces.esurface_cache.iter_enumerated() {
        let expression = esurface_linear_expr(surface)
            .substitute_external_energy(localized_external_edge, replacement)
            .canonical();
        if &expression == localized {
            return Some((HybridSurfaceID::Esurface(id), 1));
        }
        if expression == negative_localized {
            return Some((HybridSurfaceID::Esurface(id), -1));
        }
    }
    for (id, surface) in residue.surfaces.hsurface_cache.iter_enumerated() {
        let expression = hsurface_linear_expr(surface)
            .substitute_external_energy(localized_external_edge, replacement)
            .canonical();
        if &expression == localized {
            return Some((HybridSurfaceID::Hsurface(id), 1));
        }
        if expression == negative_localized {
            return Some((HybridSurfaceID::Hsurface(id), -1));
        }
    }
    for (id, surface) in residue.surfaces.linear_surface_cache.iter_enumerated() {
        let expression = surface
            .expression
            .clone()
            .substitute_external_energy(localized_external_edge, replacement)
            .canonical();
        if &expression == localized {
            return Some((HybridSurfaceID::Linear(id), 1));
        }
        if expression == negative_localized {
            return Some((HybridSurfaceID::Linear(id), -1));
        }
    }
    None
}

fn referenced_residue_surfaces(
    residue: &ThreeDExpression<OrientationID>,
) -> std::collections::BTreeSet<HybridSurfaceID> {
    let mut surfaces = std::collections::BTreeSet::new();
    for orientation in &residue.orientations {
        for variant in &orientation.variants {
            surfaces.extend(variant.numerator_surfaces.iter().copied());
            surfaces.extend(
                variant
                    .denominator
                    .iter_nodes()
                    .map(|node| node.data)
                    .filter(|surface_id| *surface_id != HybridSurfaceID::Unit),
            );
        }
    }
    surfaces
}

fn residue_surface_linear_expressions(
    residue: &ThreeDExpression<OrientationID>,
) -> Vec<(HybridSurfaceID, LinearEnergyExpr)> {
    let esurfaces = residue
        .surfaces
        .esurface_cache
        .iter_enumerated()
        .map(|(id, surface)| (HybridSurfaceID::Esurface(id), esurface_linear_expr(surface)));
    let hsurfaces = residue
        .surfaces
        .hsurface_cache
        .iter_enumerated()
        .map(|(id, surface)| (HybridSurfaceID::Hsurface(id), hsurface_linear_expr(surface)));
    let linear_surfaces =
        residue
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, surface)| {
                (
                    HybridSurfaceID::Linear(id),
                    surface.expression.clone().canonical(),
                )
            });
    esurfaces.chain(hsurfaces).chain(linear_surfaces).collect()
}

fn residue_surface_linear_expression(
    residue: &ThreeDExpression<OrientationID>,
    surface_id: HybridSurfaceID,
) -> Option<LinearEnergyExpr> {
    match surface_id {
        HybridSurfaceID::Unit => Some(LinearEnergyExpr::zero()),
        HybridSurfaceID::Infinite => None,
        HybridSurfaceID::Esurface(id) => residue
            .surfaces
            .esurface_cache
            .get(id)
            .map(esurface_linear_expr),
        HybridSurfaceID::Hsurface(id) => residue
            .surfaces
            .hsurface_cache
            .get(id)
            .map(hsurface_linear_expr),
        HybridSurfaceID::Linear(id) => residue
            .surfaces
            .linear_surface_cache
            .get(id)
            .map(|surface| surface.expression.clone().canonical()),
    }
}

fn esurface_linear_expr(surface: &Esurface) -> LinearEnergyExpr {
    surface
        .energies
        .iter()
        .fold(LinearEnergyExpr::zero(), |acc, edge_id| {
            acc + LinearEnergyExpr::ose(*edge_id, 1)
        })
        + surface
            .external_shift
            .iter()
            .fold(LinearEnergyExpr::zero(), |acc, (edge_id, coeff)| {
                acc + LinearEnergyExpr::external(*edge_id, *coeff)
            })
}

fn hsurface_linear_expr(surface: &Hsurface) -> LinearEnergyExpr {
    let positive = surface
        .positive_energies
        .iter()
        .fold(LinearEnergyExpr::zero(), |acc, edge_id| {
            acc + LinearEnergyExpr::ose(*edge_id, 1)
        });
    let negative = surface
        .negative_energies
        .iter()
        .fold(LinearEnergyExpr::zero(), |acc, edge_id| {
            acc + LinearEnergyExpr::ose(*edge_id, -1)
        });
    let external = surface
        .external_shift
        .iter()
        .fold(LinearEnergyExpr::zero(), |acc, (edge_id, coeff)| {
            acc + LinearEnergyExpr::external(*edge_id, *coeff)
        });
    positive + negative + external
}

fn intern_localized_linear_surface(
    residue: &mut ThreeDExpression<OrientationID>,
    expression: LinearEnergyExpr,
) -> LinearSurfaceID {
    let surface = LinearSurface {
        kind: LinearSurfaceKind::Hsurface,
        expression,
        origin: SurfaceOrigin::Helper,
        numerator_only: false,
    };
    residue
        .surfaces
        .linear_surface_cache
        .position(|existing| existing == &surface)
        .unwrap_or_else(|| {
            residue.surfaces.linear_surface_cache.push(surface);
            LinearSurfaceID(residue.surfaces.linear_surface_cache.len() - 1)
        })
}

pub trait GammaLoopGraphOrientation: GraphOrientation {
    fn orientation_thetas_gs(&self) -> Atom {
        let mut thetas = Atom::num(1);

        for (edge, orientation) in self.orientation() {
            match orientation {
                Orientation::Default => thetas *= GS.sign_theta(GS.sign(edge)),
                Orientation::Reversed => thetas *= GS.sign_theta(-GS.sign(edge)),
                Orientation::Undirected => {}
            }
        }

        thetas
    }

    fn orientation_delta_gs(&self) -> Atom
    where
        Self: Sized,
    {
        GS.orientation_delta(self)
    }

    fn select_gs<'a>(&self, atom: impl Into<AtomOrView<'a>>) -> Atom {
        let theta_reps = vec![
            Replacement::new(GS.sign_theta(Atom::num(1)).to_pattern(), Atom::num(1)),
            Replacement::new(GS.sign_theta(Atom::num(-1)).to_pattern(), Atom::Zero),
        ];

        let mut reps = Vec::new();

        for (edge, orientation) in self.orientation() {
            match orientation {
                Orientation::Default => {
                    reps.push(Replacement::new(GS.sign(edge).to_pattern(), Atom::num(1)));
                }
                Orientation::Reversed => {
                    reps.push(Replacement::new(GS.sign(edge).to_pattern(), Atom::num(-1)));
                }
                Orientation::Undirected => {}
            }
        }

        let orientation = self.orientation();
        atom.into()
            .replace_multiple(&reps)
            .replace_multiple(&theta_reps)
            .replace_map(|term, _ctx, out| {
                if let AtomView::Fun(function) = term
                    && function.get_symbol() == GS.orientation_delta
                {
                    if function
                        .iter()
                        .zip_longest(orientation)
                        .all(|either| match either {
                            EitherOrBoth::Both(arg, (_, orientation)) => {
                                if let Ok(arg) = i64::try_from(arg) {
                                    match orientation {
                                        Orientation::Default => arg >= 0,
                                        Orientation::Reversed => arg <= 0,
                                        Orientation::Undirected => true,
                                    }
                                } else {
                                    false
                                }
                            }
                            EitherOrBoth::Left(_) | EitherOrBoth::Right(_) => false,
                        })
                    {
                        **out = Atom::num(1);
                    } else {
                        **out = Atom::Zero;
                    }
                }
            })
    }
}

impl<T: GraphOrientation> GammaLoopGraphOrientation for T {}

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

        self.prefactor.clone()
            * half_edge_factor
            * scale_factor
            * numerator_surface_factor
            * self.denominator.to_atom_inv()
    }
}

pub trait GammaLoopOrientationExpression {
    fn to_atom_gs(&self) -> Atom;
    fn energy_replacements_gs(&self, graph: &Graph) -> Vec<Replacement>;
    fn numerator_atom_gs(&self, graph: &Graph, numerator: &Atom) -> Atom;
    fn parametric_atom_without_orientation_thetas_gs<E, H>(
        &self,
        graph: &Graph,
        numerator: &Atom,
        surfaces: &three_dimensional_reps::surface::SurfaceCache<E, H>,
    ) -> Atom
    where
        three_dimensional_reps::surface::SurfaceCache<E, H>: GammaLoopSurfaceCache;
    fn parametric_atom_gs<E, H>(
        &self,
        graph: &Graph,
        numerator: &Atom,
        surfaces: &three_dimensional_reps::surface::SurfaceCache<E, H>,
    ) -> Atom
    where
        three_dimensional_reps::surface::SurfaceCache<E, H>: GammaLoopSurfaceCache;
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
        let mut replacements = Vec::new();
        let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);

        for (edge_id, energy_expr) in self.edge_energy_map.iter().enumerate() {
            let edge_id = linnet::half_edge::involution::EdgeIndex(edge_id);
            let energy = energy_expr.to_atom_gs(&[]);
            replacements.push(Replacement::new(
                GS.emr_mom(edge_id, AIND_SYMBOLS.cind.f([Atom::Zero]))
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
            let energy = self
                .edge_energy_map
                .get(usize::from(*loop_edge_id))
                .map(|energy_expr| energy_expr.to_atom_gs(&[]))
                .unwrap_or_else(Atom::new);
            replacements.push(Replacement::new(
                function!(
                    GS.loop_mom,
                    loop_id_atom.clone(),
                    AIND_SYMBOLS.cind.f([Atom::Zero])
                )
                .to_pattern(),
                energy.clone().to_pattern(),
            ));
            for spatial_index in 1..=3 {
                replacements.push(Replacement::new(
                    function!(
                        GS.loop_mom,
                        loop_id_atom.clone(),
                        AIND_SYMBOLS.cind.f([spatial_index])
                    )
                    .to_pattern(),
                    GS.emr_mom(*loop_edge_id, AIND_SYMBOLS.cind.f([spatial_index]))
                        .to_pattern(),
                ));
            }
            replacements.push(Replacement::new(
                FunctionBuilder::new(GS.loop_mom)
                    .add_arg(loop_id as i64)
                    .add_arg(mink_index.as_view())
                    .finish()
                    .to_pattern(),
                (GS.emr_vec_index(*loop_edge_id, &mink_index)
                    + energy * GS.energy_delta(&mink_index))
                .to_pattern(),
            ));
        }

        replacements
    }

    fn numerator_atom_gs(&self, graph: &Graph, numerator: &Atom) -> Atom {
        numerator.replace_multiple(self.energy_replacements_gs(graph))
    }

    fn parametric_atom_without_orientation_thetas_gs<E, H>(
        &self,
        graph: &Graph,
        numerator: &Atom,
        surfaces: &three_dimensional_reps::surface::SurfaceCache<E, H>,
    ) -> Atom
    where
        three_dimensional_reps::surface::SurfaceCache<E, H>: GammaLoopSurfaceCache,
    {
        self.numerator_atom_gs(graph, numerator)
            * surfaces.substitute_energies_gs(&self.to_atom_gs(), &[])
    }

    fn parametric_atom_gs<E, H>(
        &self,
        graph: &Graph,
        numerator: &Atom,
        surfaces: &three_dimensional_reps::surface::SurfaceCache<E, H>,
    ) -> Atom
    where
        three_dimensional_reps::surface::SurfaceCache<E, H>: GammaLoopSurfaceCache,
    {
        self.parametric_atom_without_orientation_thetas_gs(graph, numerator, surfaces)
            * self.orientation_thetas_gs()
    }
}

pub trait GammaLoopThreeDExpression<O> {
    fn to_atom_gs(&self, pattern: &OrientationPattern) -> Atom;
    fn get_orientation_atoms_gs(
        &self,
        pattern: &OrientationPattern,
    ) -> typed_index_collections::TiVec<OrientationID, Atom>;
    fn get_orientation_atom_gs(&self, orientation_id: O) -> Atom;
    fn diagnostic_parametric_atom_with_numerator_gs(
        &self,
        graph: &Graph,
        numerator: &Atom,
        pattern: &OrientationPattern,
    ) -> Atom;
    fn parametric_atom_with_numerator_gs(
        &self,
        graph: &Graph,
        numerator: &Atom,
        pattern: &OrientationPattern,
    ) -> Atom;
    fn diagnostic_parametric_atom_gs(&self, graph: &Graph, pattern: &OrientationPattern) -> Atom;
    fn parametric_atom_gs(&self, graph: &Graph, pattern: &OrientationPattern) -> Atom;
}

impl<O, E, H> GammaLoopThreeDExpression<O>
    for three_dimensional_reps::expression::ThreeDExpression<O, E, H>
where
    O: Clone + From<usize> + Into<usize>,
    three_dimensional_reps::surface::SurfaceCache<E, H>: GammaLoopSurfaceCache,
    usize: From<O>,
{
    fn to_atom_gs(&self, pattern: &OrientationPattern) -> Atom {
        self.orientations
            .iter()
            .filter_map(|orientation| {
                if pattern.filter_orientation(orientation.orientation()) {
                    Some(orientation.to_atom_gs())
                } else {
                    None
                }
            })
            .reduce(|a, b| a + b)
            .unwrap_or_default()
    }

    fn get_orientation_atoms_gs(
        &self,
        pattern: &OrientationPattern,
    ) -> typed_index_collections::TiVec<OrientationID, Atom> {
        self.orientations
            .iter()
            .map(|orientation| {
                if pattern.filter_orientation(orientation.orientation()) {
                    orientation.to_atom_gs()
                } else {
                    Atom::new()
                }
            })
            .collect()
    }

    fn get_orientation_atom_gs(&self, orientation_id: O) -> Atom {
        self.orientations[orientation_id].to_atom_gs()
    }

    fn diagnostic_parametric_atom_with_numerator_gs(
        &self,
        graph: &Graph,
        numerator: &Atom,
        pattern: &OrientationPattern,
    ) -> Atom {
        self.orientations
            .iter()
            .filter_map(|orientation| {
                if pattern.filter_orientation(orientation.orientation()) {
                    Some(orientation.parametric_atom_without_orientation_thetas_gs(
                        graph,
                        numerator,
                        &self.surfaces,
                    ))
                } else {
                    None
                }
            })
            .reduce(|acc, atom| acc + atom)
            .unwrap_or_default()
    }

    fn parametric_atom_with_numerator_gs(
        &self,
        graph: &Graph,
        numerator: &Atom,
        pattern: &OrientationPattern,
    ) -> Atom {
        self.orientations
            .iter()
            .filter_map(|orientation| {
                if pattern.filter_orientation(orientation.orientation()) {
                    Some(orientation.parametric_atom_gs(graph, numerator, &self.surfaces))
                } else {
                    None
                }
            })
            .reduce(|acc, atom| acc + atom)
            .unwrap_or_default()
    }

    fn diagnostic_parametric_atom_gs(&self, graph: &Graph, pattern: &OrientationPattern) -> Atom {
        self.diagnostic_parametric_atom_with_numerator_gs(
            graph,
            &graph.full_numerator_atom(),
            pattern,
        )
    }

    fn parametric_atom_gs(&self, graph: &Graph, pattern: &OrientationPattern) -> Atom {
        self.parametric_atom_with_numerator_gs(graph, &graph.full_numerator_atom(), pattern)
    }
}

pub fn numerator_with_positive_internal_ose_gs(graph: &Graph) -> Atom {
    graph
        .full_numerator_atom()
        .replace_multiple(positive_internal_ose_energy_replacements_gs(graph))
}

pub fn numerator_with_internal_energy_parameters_gs(graph: &Graph) -> Atom {
    graph
        .full_numerator_atom()
        .replace_multiple(internal_energy_parameter_replacements_gs(graph))
}

pub fn internal_energy_parameter_atom_gs(edge_id: EdgeIndex) -> Atom {
    GS.emr_mom(edge_id, AIND_SYMBOLS.cind.f([Atom::Zero]))
}

pub fn positive_internal_ose_energy_replacements_gs(graph: &Graph) -> Vec<Replacement> {
    let mut replacements = Vec::new();
    let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);

    for (pair, edge_id, _) in graph.underlying.iter_edges() {
        if !matches!(pair, HedgePair::Paired { .. }) {
            continue;
        }
        let energy = ose_atom_from_index(edge_id);
        replacements.push(Replacement::new(
            GS.emr_mom(edge_id, AIND_SYMBOLS.cind.f([Atom::Zero]))
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
        let energy = ose_atom_from_index(*loop_edge_id);
        replacements.push(Replacement::new(
            function!(
                GS.loop_mom,
                loop_id_atom.clone(),
                AIND_SYMBOLS.cind.f([Atom::Zero])
            )
            .to_pattern(),
            energy.clone().to_pattern(),
        ));
        for spatial_index in 1..=3 {
            replacements.push(Replacement::new(
                function!(
                    GS.loop_mom,
                    loop_id_atom.clone(),
                    AIND_SYMBOLS.cind.f([spatial_index])
                )
                .to_pattern(),
                GS.emr_mom(*loop_edge_id, AIND_SYMBOLS.cind.f([spatial_index]))
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

pub fn internal_energy_parameter_replacements_gs(graph: &Graph) -> Vec<Replacement> {
    let mut replacements = Vec::new();
    let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);

    for (pair, edge_id, _) in graph.underlying.iter_edges() {
        if !matches!(pair, HedgePair::Paired { .. }) {
            continue;
        }
        let energy = internal_energy_parameter_atom_gs(edge_id);
        replacements.push(Replacement::new(
            GS.emr_mom(edge_id, AIND_SYMBOLS.cind.f([Atom::Zero]))
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
        let energy = internal_energy_parameter_atom_gs(*loop_edge_id);
        replacements.push(Replacement::new(
            function!(
                GS.loop_mom,
                loop_id_atom.clone(),
                AIND_SYMBOLS.cind.f([Atom::Zero])
            )
            .to_pattern(),
            energy.clone().to_pattern(),
        ));
        for spatial_index in 1..=3 {
            replacements.push(Replacement::new(
                function!(
                    GS.loop_mom,
                    loop_id_atom.clone(),
                    AIND_SYMBOLS.cind.f([spatial_index])
                )
                .to_pattern(),
                GS.emr_mom(*loop_edge_id, AIND_SYMBOLS.cind.f([spatial_index]))
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

pub fn orientation_delta_gs(orientation: &EdgeVec<Orientation>) -> Atom {
    orientation.orientation_delta_gs()
}
