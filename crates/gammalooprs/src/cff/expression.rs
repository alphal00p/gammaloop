pub use three_dimensional_reps::expression::{
    AllOrientations, CFFVariant, GraphOrientation, OrientationData, OrientationExpression,
    OrientationID, OrientationSelector, RaisedEsurfaceData, RaisedEsurfaceDataView,
    RaisedEsurfaceGroup, RaisedEsurfaceGroupView, RaisedEsurfaceId, ResidualDenominator,
};

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
    esurface::Esurface,
    hsurface::Hsurface,
    surface::{GammaLoopLinearEnergyExpr, GammaLoopSurfaceCache},
};
use crate::{
    graph::Graph,
    settings::global::OrientationPattern,
    utils::{GS, W_, ose_atom_from_index, symbolica_ext::CallSymbol},
};

pub type ThreeDExpression<O> =
    three_dimensional_reps::expression::ThreeDExpression<O, Esurface, Hsurface>;
pub type CFFExpression<O> =
    three_dimensional_reps::expression::CFFExpression<O, Esurface, Hsurface>;

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
        numerator.replace_multiple(&self.energy_replacements_gs(graph))
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
        .replace_multiple(&positive_internal_ose_energy_replacements_gs(graph))
}

pub fn numerator_with_internal_energy_parameters_gs(graph: &Graph) -> Atom {
    graph
        .full_numerator_atom()
        .replace_multiple(&internal_energy_parameter_replacements_gs(graph))
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
