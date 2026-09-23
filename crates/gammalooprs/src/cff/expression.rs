<<<<<<< pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (rebased revision)
use std::collections::{BTreeMap, BTreeSet};
||||||| pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (parents of rebased revision)
#![allow(dead_code)]
=======
//! GammaLoop extensions for FeynKit-owned CFF expressions and orientations.
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)

<<<<<<< pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (rebased revision)
pub use three_dimensional_reps::expression::{
    AllOrientations, CFFVariant, OrientationData, OrientationExpression, OrientationID,
    OrientationSelector, RaisedEsurfaceData, RaisedEsurfaceDataView, RaisedEsurfaceGroup,
    RaisedEsurfaceGroupView, RaisedEsurfaceId,
};
||||||| pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (parents of rebased revision)
use std::{borrow::Borrow, collections::HashMap, fmt::Display};
=======
use std::borrow::Borrow;
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)

<<<<<<< pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (rebased revision)
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski, RepName},
};
||||||| pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (parents of rebased revision)
use crate::{
    cff::{
        esurface::{EsurfaceID, RaisedEsurfaceData, RaisedEsurfaceGroup},
        orientations::GraphOrientation,
    },
    settings::global::OrientationPattern,
    utils::{W_, ose_atom_from_index},
};
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use linnet::half_edge::{
    GVEdgeAttrs, HedgeGraph,
    involution::{EdgeVec, Orientation, SignOrZero},
    nodestore::NodeStorageOps,
};
use serde::{Deserialize, Serialize};
=======
use feynkit_cff::{CffResult, OrientationData, OrientationExpression, OrientationId};
use linnet::half_edge::involution::{EdgeVec, Orientation};
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)
use symbolica::{
<<<<<<< pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (rebased revision)
    atom::{Atom, AtomCore, FunctionBuilder},
    function,
    id::Replacement,
};

use super::{esurface::Esurface, hsurface::Hsurface, surface::GammaLoopLinearEnergyExpr};
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
||||||| pvqotmrp c31c29c1 "lower UFO numerators to native Spenso tensors" (parents of rebased revision)
    atom::{Atom, AtomCore, AtomOrView, Symbol},
    function,
    id::{Pattern, Replacement},
    symbol,
};
use tabled::{builder::Builder, settings::Style};
use typed_index_collections::TiVec;

use super::{generation::SurfaceCache, surface::HybridSurfaceID, tree::Tree};

#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    From,
    Into,
    Hash,
    PartialEq,
    Eq,
    Copy,
    Encode,
    Decode,
    PartialOrd,
    Ord,
)]
pub struct OrientationID(pub usize);

impl GraphOrientation for EdgeVec<Orientation> {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        self
    }
}

impl OrientationID {
    pub fn symbol() -> Symbol {
        symbol!("sigma")
    }

    pub fn atom(self) -> Atom {
        let id: usize = self.into();
        function!(Self::symbol(), id as i64)
    }

    pub fn select<'a>(self, atom: impl Into<AtomOrView<'a>>) -> Atom {
        atom.into()
            .as_view()
            .replace(self.atom())
            .with(Atom::num(1))
            .replace(function!(Self::symbol(), W_.x_))
            .with(Atom::Zero)
    }
}

#[derive(
    Clone, Debug, PartialOrd, Ord, Hash, PartialEq, Eq, Serialize, Deserialize, Encode, Decode,
)]
pub struct OrientationData {
    pub orientation: EdgeVec<Orientation>,
}

impl GraphOrientation for OrientationData {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.orientation
    }
}

impl Display for OrientationData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = Builder::new();

        for (i, item) in self.orientation.iter() {
            table.push_column(&[i.to_string(), SignOrZero::from(*item).to_string()]);
        }
        table.build().with(Style::rounded()).fmt(f)
    }
}

impl OrientationData {
    pub(crate) fn dot<E, V, N: NodeStorageOps<NodeData = V>>(&self, graph: &HedgeGraph<E, V, N>) {
        let mut writer = String::new();
        writer.push_str("digraph {{");

        writer.push_str(
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";",
        );

        for (hedge_pair, id, _) in graph.iter_edges() {
            let attr = GVEdgeAttrs {
                color: None,
                label: None,
                other: None,
            };
            writer.push_str("  ");

            let attr = hedge_pair.fill_color(attr);
            hedge_pair
                .add_data(graph)
                .dot_fmt(
                    &mut writer,
                    graph,
                    id,
                    |_| None,
                    |a| a.to_string(),
                    self.orientation[id],
                    attr,
                )
                .unwrap();
        }
        writer.push_str("}}");
    }

    pub(crate) fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.orientation
            .borrow()
            .into_iter()
            .filter_map(|(edge_index, orientation)| {
                if matches!(orientation, Orientation::Reversed) {
                    let energy_atom = ose_atom_from_index(edge_index);

                    let neg_energy_atom = -&energy_atom;

                    Some(Replacement::new(
                        Pattern::from(energy_atom),
                        Pattern::from(neg_energy_atom),
                    ))
                } else {
                    None
                }
            })
            .collect()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationExpression {
    pub data: OrientationData,
    pub expression: Tree<HybridSurfaceID>,
}
#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct CFFExpression<O>
where
    O: From<usize> + Into<usize>,
{
    pub orientations: TiVec<O, OrientationExpression>,
    pub surfaces: SurfaceCache,
}

impl GraphOrientation for OrientationExpression {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.data.orientation
    }
}

impl<O: Clone + From<usize> + Into<usize>> CFFExpression<O>
where
    usize: From<O>,
{
    pub(crate) fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
            surfaces: SurfaceCache {
                esurface_cache: TiVec::new(),
                hsurface_cache: TiVec::new(),
            },
        }
    }

    pub fn to_atom(&self, pattern: OrientationPattern) -> Atom {
        self.orientations
            .iter()
            .filter_map(|orientation| {
                if pattern.filter(orientation.orientation()) {
                    Some(orientation.expression.to_atom_inv())
                } else {
                    None
                }
            })
            .reduce(|a, b| a + b)
            .unwrap_or_default()
    }

    pub(crate) fn get_orientation_atoms(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<OrientationID, Atom> {
        self.orientations
            .iter()
            .map(|orientation| {
                if pattern.filter(orientation.orientation()) {
                    orientation.expression.to_atom_inv()
                } else {
                    Atom::new()
                }
            })
            .collect()
    }

    pub fn get_orientation_atoms_with_data(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<OrientationID, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = if pattern.filter(orientation.orientation()) {
                    orientation.expression.to_atom_inv()
                } else {
                    Atom::new()
                };
                let data = orientation.data.clone();
                (atom, data)
            })
            .collect()
    }

    pub(crate) fn get_orientation_atom(&self, orientation_id: O) -> Atom {
        self.orientations[orientation_id].expression.to_atom_inv()
    }

    pub(crate) fn num_unfolded_terms(&self) -> usize {
        self.orientations
            .iter()
            .map(|o| o.expression.get_bottom_layer().len())
            .sum()
    }

    pub(crate) fn get_orientations_with_esurface(&self, esurface_id: EsurfaceID) -> Vec<O> {
        self.orientations
            .iter_enumerated()
            .filter_map(|(id, orientation)| {
                if orientation
                    .expression
                    .iter_nodes()
                    .any(|node| node.data == HybridSurfaceID::Esurface(esurface_id))
                {
                    Some(id)
                } else {
                    None
                }
            })
            .collect()
    }

    pub(crate) fn select_esurface_residue(
        mut self,
        raised_esurface_group: &RaisedEsurfaceGroup,
    ) -> Vec<CFFExpression<O>> {
        self.normalize_single_raising(raised_esurface_group);

        let reprentative_esurface_id = raised_esurface_group.esurface_ids[0];

        let mut result = vec![];
        for occurence in 1..=raised_esurface_group.max_occurence {
            let mut new_expression = self.clone();

            for orientation in new_expression.orientations.iter_mut() {
                orientation.expression.keep_branches_with_value_count_mut(
                    &HybridSurfaceID::Esurface(reprentative_esurface_id),
                    occurence,
                );

                orientation
                    .expression
                    .map_mut(|hybrid_surface_id| match hybrid_surface_id {
                        HybridSurfaceID::Esurface(esurface_id)
                            if *esurface_id == reprentative_esurface_id =>
                        {
                            *hybrid_surface_id = HybridSurfaceID::Unit;
                        }
                        _ => (),
                    });
            }

            result.push(new_expression);
        }

        result
    }

    pub(crate) fn normalize_single_raising(&mut self, raised_esurface_group: &RaisedEsurfaceGroup) {
        let reprentative_esurface_id = raised_esurface_group.esurface_ids[0];

        for orientation in self.orientations.iter_mut() {
            orientation
                .expression
                .map_mut(|hybrid_surface_id| match hybrid_surface_id {
                    HybridSurfaceID::Esurface(esurface_id)
                        if raised_esurface_group.esurface_ids.contains(esurface_id) =>
                    {
                        *hybrid_surface_id = HybridSurfaceID::Esurface(reprentative_esurface_id);
                    }
                    _ => (),
                });
        }
    }

    pub(crate) fn normalize_wrt_all_raisings(&mut self, raised_data: &RaisedEsurfaceData) {
        let mut esurface_mappings = HashMap::new();

        for cut_group in raised_data.raised_groups.iter() {
            let esurface_id_of_first = cut_group.esurface_ids[0];

            for esurface_id in cut_group.esurface_ids.iter() {
                esurface_mappings.insert(*esurface_id, esurface_id_of_first);
            }
        }

        for orientation in self.orientations.iter_mut() {
            orientation.expression.map_mut(|hybrid_surface_id| {
                if let HybridSurfaceID::Esurface(esurface_id) = hybrid_surface_id
                    && let Some(normalized_esurface_id) = esurface_mappings.get(esurface_id)
                {
                    *esurface_id = *normalized_esurface_id;
                }
            });
        }
=======
    atom::Atom,
    id::{Pattern, Replacement},
};
use typed_index_collections::TiVec;

use crate::{
    cff::orientations::GraphOrientation, settings::global::OrientationPattern,
    utils::ose_atom_from_index,
};

impl GraphOrientation for OrientationData {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.orientation
    }
}

impl GraphOrientation for feynkit_cff::GraphOrientation {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        self
    }
}

impl GraphOrientation for OrientationExpression {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.data.orientation
    }
}

/// GammaLoop diagnostics and symbolic rewrites for a FeynKit orientation.
pub trait OrientationDataExt {
    fn get_ose_replacements(&self) -> Vec<Replacement>;
}

impl OrientationDataExt for OrientationData {
    fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.orientation
            .borrow()
            .into_iter()
            .filter(|(_, orientation)| **orientation == Orientation::Reversed)
            .map(|(edge, _)| {
                let energy = ose_atom_from_index(edge);
                Replacement::new(Pattern::from(energy.clone()), Pattern::from(-energy))
            })
            .collect()
    }
}

/// GammaLoop's symbolic lowering and residue API over a canonical FeynKit result.
pub trait CffExpressionExt {
    fn to_atom(&self, pattern: OrientationPattern) -> Atom;
    fn get_orientation_atoms(&self, pattern: OrientationPattern) -> TiVec<OrientationId, Atom>;
    fn get_orientation_atoms_with_data(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<OrientationId, (Atom, OrientationData)>;
}

impl CffExpressionExt for CffResult {
    fn to_atom(&self, pattern: OrientationPattern) -> Atom {
        self.orientations
            .iter()
            .filter(|orientation| pattern.filter(orientation.orientation()))
            .map(|orientation| orientation.expression.to_atom_inverse())
            .reduce(|sum, term| sum + term)
            .unwrap_or_default()
    }

    fn get_orientation_atoms(&self, pattern: OrientationPattern) -> TiVec<OrientationId, Atom> {
        self.orientations
            .iter()
            .map(|orientation| {
                if pattern.filter(orientation.orientation()) {
                    orientation.expression.to_atom_inverse()
                } else {
                    Atom::Zero
                }
            })
            .collect()
    }

    fn get_orientation_atoms_with_data(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<OrientationId, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = if pattern.filter(orientation.orientation()) {
                    orientation.expression.to_atom_inverse()
                } else {
                    Atom::Zero
                };
                (atom, orientation.data.clone())
            })
            .collect()
>>>>>>> xzuonylm 3e5f0ddf (rebased revision)
    }
}
