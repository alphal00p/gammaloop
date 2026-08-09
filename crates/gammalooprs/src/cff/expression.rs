use std::collections::BTreeMap;

pub use three_dimensional_reps::expression::{
    AllOrientations, CFFVariant, OrientationData, OrientationExpression, OrientationID,
    OrientationSelector, RaisedEsurfaceData, RaisedEsurfaceDataView, RaisedEsurfaceGroup,
    RaisedEsurfaceGroupView, RaisedEsurfaceId,
};

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

use super::{esurface::Esurface, hsurface::Hsurface};
use crate::{
    graph::Graph,
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
            let edge_id = EdgeIndex(edge_id);
            let energy = super::surface::GammaLoopLinearEnergyExpr::to_atom_gs(energy_expr, &[]);
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
            let energy = self
                .edge_energy_map
                .get(usize::from(*loop_edge_id))
                .map(|energy_expr| {
                    super::surface::GammaLoopLinearEnergyExpr::to_atom_gs(energy_expr, &[])
                })
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
                (GS.emr_vec_index(*loop_edge_id, &mink_index)
                    + energy * GS.energy_delta(&mink_index))
                .to_pattern(),
            ));
        }

        replacements
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::{EdgeVec, Orientation};
    use symbolica::atom::AtomView;

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
            internal_terms: vec![(EdgeIndex(0), Atom::num(2)), (EdgeIndex(1), Atom::num(-3))],
            external_terms: vec![(EdgeIndex(2), Atom::num(5))],
            uniform_scale_coeff: Atom::num(7),
            constant: Atom::num(11),
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
        let AtomView::Mul(product) = mapped.as_view() else {
            panic!("mapped numerator should remain a product");
        };
        assert_eq!(
            product
                .iter()
                .filter(|factor| matches!(factor, AtomView::Add(_)))
                .count(),
            2,
        );

        let mixed = (mapped_energy.clone() + factor_a) * (mapped_energy + Atom::num(1) + factor_b);
        assert_ne!(mapped, mixed);
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
