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
};

use super::{
    esurface::{Esurface, RaisedEsurfaceGroup as GammaLoopRaisedEsurfaceGroup},
    hsurface::Hsurface,
    surface::{
        EsurfaceID, GammaLoopLinearEnergyExpr, GammaLoopSurfaceCache, HybridSurfaceID,
        LinearEnergyExpr, LinearSurface, LinearSurfaceID, LinearSurfaceKind, SurfaceOrigin,
    },
};
use crate::{
    graph::Graph,
    settings::global::{OrientationPattern, ThreeDRepresentation},
    utils::{GS, W_, ose_atom_from_index, symbolica_ext::CallSymbol},
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
            })
        });
    }
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
    if representation == ThreeDRepresentation::Ltd {
        // GammaLoop Cutkosky cuts are oriented positive-energy constraints.
        // LTD residues are selected from the canonical generated E-surface
        // variables, so both simple and confluent LU cuts need the effective
        // selected-variable signs supplied by the cross-section cut resolver.
        // Any LTD bridge factor that cannot be represented as a selected
        // variable sign is applied outside this selector.
        expression.select_esurface_residue_with_cut_edges_and_esurface_signs(
            lu_cut,
            cut_edge_sets,
            ltd_lu_cut_esurface_signs,
        )
    } else {
        expression.select_esurface_residue_with_cut_edges(lu_cut, cut_edge_sets)
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
