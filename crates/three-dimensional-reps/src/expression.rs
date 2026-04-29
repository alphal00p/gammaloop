use std::{borrow::Borrow, collections::HashMap, fmt::Display};

use bincode_trait_derive::{Decode, Encode};
use itertools::{EitherOrBoth, Itertools};
use linnet::half_edge::{
    GVEdgeAttrs, HedgeGraph,
    involution::{EdgeIndex, EdgeVec, Orientation, SignOrZero},
    nodestore::NodeStorageOps,
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol},
    function,
    id::{Pattern, Replacement},
    symbol,
};
use typed_index_collections::TiVec;

use crate::{
    graph_io::EnergyEdgeIndexMap,
    surface::{EsurfaceID, HybridSurfaceID, LinearEnergyExpr, RationalCoefficient, SurfaceCache},
    symbols::{
        numerator_sampling_scale, orientation_delta, orientation_delta_symbol, ose_atom_from_index,
        sign, sign_theta,
    },
    tree::{Tree, TreeNode},
};

#[derive(
    Debug, Clone, Serialize, Deserialize, Hash, PartialEq, Eq, Copy, Encode, Decode, PartialOrd, Ord,
)]
pub struct OrientationID(pub usize);

impl From<usize> for OrientationID {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<OrientationID> for usize {
    fn from(value: OrientationID) -> Self {
        value.0
    }
}

impl GraphOrientation for EdgeVec<Orientation> {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        self
    }
}

pub trait GraphOrientation {
    fn orientation(&self) -> &EdgeVec<Orientation>;

    fn orientation_thetas(&self) -> Atom {
        let mut thetas = Atom::num(1);

        for (e, h) in self.orientation() {
            match h {
                Orientation::Default => {
                    thetas *= sign_theta(sign(e));
                }
                Orientation::Reversed => {
                    thetas *= sign_theta(-sign(e));
                }
                Orientation::Undirected => {}
            }
        }
        thetas
    }

    fn orientation_delta(&self) -> Atom
    where
        Self: Sized,
    {
        orientation_delta(self)
    }

    fn select<'a>(&self, atom: impl Into<AtomOrView<'a>>) -> Atom {
        let theta_reps = vec![
            Replacement::new(sign_theta(Atom::num(1)).to_pattern(), Atom::num(1)),
            Replacement::new(sign_theta(Atom::num(-1)).to_pattern(), Atom::Zero),
        ];

        let mut reps = Vec::new();

        for (e, h) in self.orientation() {
            match h {
                Orientation::Default => {
                    reps.push(Replacement::new(sign(e).to_pattern(), Atom::num(1)));
                }
                Orientation::Reversed => {
                    reps.push(Replacement::new(sign(e).to_pattern(), Atom::num(-1)));
                }
                Orientation::Undirected => {}
            }
        }

        let orientation = self.orientation();
        atom.into()
            .replace_multiple(&reps)
            .replace_multiple(&theta_reps)
            .replace_map(|term, _ctx, out| {
                if let AtomView::Fun(f) = term
                    && f.get_symbol() == orientation_delta_symbol()
                {
                    if f.iter()
                        .zip_longest(orientation)
                        .all(|either| match either {
                            EitherOrBoth::Both(a, (_, o)) => {
                                if let Ok(a) = i64::try_from(a) {
                                    match o {
                                        Orientation::Default => a >= 0,
                                        Orientation::Reversed => a <= 0,
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
            .replace(function!(Self::symbol(), symbol!("x_")))
            .with(Atom::Zero)
    }
}

#[derive(
    Clone, Debug, PartialOrd, Ord, Hash, PartialEq, Eq, Serialize, Deserialize, Encode, Decode,
)]
pub struct OrientationData {
    pub orientation: EdgeVec<Orientation>,
    pub label: Option<String>,
    pub numerator_map_index: Option<usize>,
}

impl GraphOrientation for OrientationData {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.orientation
    }
}

impl Display for OrientationData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let entries = self
            .orientation
            .iter()
            .map(|(i, item)| format!("{i}:{}", SignOrZero::from(*item)))
            .join(", ");
        f.write_str(&entries)
    }
}

impl OrientationData {
    pub fn new(orientation: EdgeVec<Orientation>) -> Self {
        let label = Some(format_graph_orientation_label(&orientation));
        Self {
            orientation,
            label,
            numerator_map_index: None,
        }
    }

    pub fn dot<E, V, N: NodeStorageOps<NodeData = V>>(&self, graph: &HedgeGraph<E, V, N>) {
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

    pub fn get_ose_replacements(&self) -> Vec<Replacement> {
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
pub struct CFFVariant {
    pub origin: Option<String>,
    pub prefactor: RationalCoefficient,
    pub half_edges: Vec<EdgeIndex>,
    pub uniform_scale_power: usize,
    pub numerator_surfaces: Vec<HybridSurfaceID>,
    pub denominator: Tree<HybridSurfaceID>,
}

impl CFFVariant {
    pub fn pure_cff(denominator: Tree<HybridSurfaceID>) -> Self {
        Self {
            origin: Some("cff".to_string()),
            prefactor: RationalCoefficient::one(),
            half_edges: Vec::new(),
            uniform_scale_power: 0,
            numerator_surfaces: Vec::new(),
            denominator,
        }
    }

    pub fn to_atom(&self) -> Atom {
        let half_edge_factor = self
            .half_edges
            .iter()
            .map(|edge_id| Atom::num(1) / (Atom::num(2) * ose_atom_from_index(*edge_id)))
            .reduce(|acc, factor| acc * factor)
            .unwrap_or_else(|| Atom::num(1));

        let scale_factor = if self.uniform_scale_power == 0 {
            Atom::num(1)
        } else {
            Atom::num(1) / numerator_sampling_scale().pow(self.uniform_scale_power as i64)
        };

        let numerator_surface_factor = self
            .numerator_surfaces
            .iter()
            .map(|surface_id| Atom::from(*surface_id))
            .reduce(|acc, factor| acc * factor)
            .unwrap_or_else(|| Atom::num(1));

        self.prefactor.to_atom()
            * half_edge_factor
            * scale_factor
            * numerator_surface_factor
            * self.denominator.to_atom_inv()
    }

    pub fn remap_energy_edge_indices(&mut self, edge_map: &EnergyEdgeIndexMap) {
        self.half_edges = self
            .half_edges
            .iter()
            .map(|edge_id| {
                EdgeIndex(
                    edge_map
                        .internal
                        .get(&edge_id.0)
                        .copied()
                        .unwrap_or(edge_id.0),
                )
            })
            .collect();
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationExpression {
    pub data: OrientationData,
    pub loop_energy_map: Vec<LinearEnergyExpr>,
    pub edge_energy_map: Vec<LinearEnergyExpr>,
    pub variants: Vec<CFFVariant>,
}

impl OrientationExpression {
    pub fn pure_cff(data: OrientationData, denominator: Tree<HybridSurfaceID>) -> Self {
        let edge_energy_map = data
            .orientation
            .iter()
            .map(|(edge_id, orientation)| match orientation {
                Orientation::Default => LinearEnergyExpr::ose(edge_id, 1),
                Orientation::Reversed => LinearEnergyExpr::ose(edge_id, -1),
                Orientation::Undirected => LinearEnergyExpr::zero(),
            })
            .collect();

        Self {
            data,
            loop_energy_map: Vec::new(),
            edge_energy_map,
            variants: vec![CFFVariant::pure_cff(denominator)],
        }
    }

    pub fn to_atom(&self) -> Atom {
        self.variants
            .iter()
            .map(CFFVariant::to_atom)
            .reduce(|acc, atom| acc + atom)
            .unwrap_or_else(Atom::new)
    }

    pub fn iter_denominator_nodes(&self) -> impl Iterator<Item = &TreeNode<HybridSurfaceID>> {
        self.variants
            .iter()
            .flat_map(|variant| variant.denominator.iter_nodes())
    }

    pub fn for_each_denominator_tree_mut(&mut self, mut f: impl FnMut(&mut Tree<HybridSurfaceID>)) {
        for variant in &mut self.variants {
            f(&mut variant.denominator);
        }
    }

    pub fn num_unfolded_terms(&self) -> usize {
        self.variants
            .iter()
            .map(|variant| variant.denominator.get_bottom_layer().len())
            .sum()
    }

    pub fn max_denominator_value_count_on_branch(&self, value: &HybridSurfaceID) -> usize {
        self.variants
            .iter()
            .map(|variant| variant.denominator.max_value_count_on_branch(value))
            .max()
            .unwrap_or(0)
    }

    pub fn remap_energy_edge_indices(&mut self, edge_map: &EnergyEdgeIndexMap) {
        self.loop_energy_map = self
            .loop_energy_map
            .drain(..)
            .map(|expr| expr.remap_energy_edges(&edge_map.internal, &edge_map.external))
            .collect();

        let mut remapped_edge_energy_map =
            vec![LinearEnergyExpr::zero(); edge_map.orientation_edge_count];
        for (local_edge_id, expr) in self.edge_energy_map.drain(..).enumerate() {
            let target_edge_id = edge_map
                .internal
                .get(&local_edge_id)
                .copied()
                .unwrap_or(local_edge_id);
            if target_edge_id >= remapped_edge_energy_map.len() {
                remapped_edge_energy_map.resize(target_edge_id + 1, LinearEnergyExpr::zero());
            }
            remapped_edge_energy_map[target_edge_id] =
                expr.remap_energy_edges(&edge_map.internal, &edge_map.external);
        }
        self.edge_energy_map = remapped_edge_energy_map;

        let mut remapped_orientation =
            EdgeVec::from_iter((0..self.edge_energy_map.len()).map(|_| Orientation::Undirected));
        for (local_edge_id, orientation) in self.data.orientation.iter() {
            let target_edge_id = edge_map
                .internal
                .get(&local_edge_id.0)
                .copied()
                .unwrap_or(local_edge_id.0);
            if target_edge_id < self.edge_energy_map.len() {
                remapped_orientation[EdgeIndex(target_edge_id)] = *orientation;
            }
        }
        self.data.orientation = remapped_orientation;
        self.data.label = Some(format_graph_orientation_label(&self.data.orientation));

        for variant in &mut self.variants {
            variant.remap_energy_edge_indices(edge_map);
        }
    }
}

impl GraphOrientation for OrientationExpression {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.data.orientation
    }
}

fn format_graph_orientation_label(orientation: &EdgeVec<Orientation>) -> String {
    orientation
        .iter()
        .map(|(_, orientation)| match *orientation {
            Orientation::Default => '+',
            Orientation::Reversed => '-',
            Orientation::Undirected => 'x',
        })
        .collect()
}

pub trait OrientationSelector {
    fn filter_orientation(&self, orientation: &EdgeVec<Orientation>) -> bool;
}

#[derive(Clone, Debug, Default)]
pub struct AllOrientations;

impl OrientationSelector for AllOrientations {
    fn filter_orientation(&self, _orientation: &EdgeVec<Orientation>) -> bool {
        true
    }
}

impl<F> OrientationSelector for F
where
    F: Fn(&EdgeVec<Orientation>) -> bool,
{
    fn filter_orientation(&self, orientation: &EdgeVec<Orientation>) -> bool {
        self(orientation)
    }
}

pub trait RaisedEsurfaceGroupView {
    fn esurface_ids(&self) -> &[EsurfaceID];
    fn max_occurrence(&self) -> usize;
}

pub trait RaisedEsurfaceDataView {
    fn for_each_raised_group(&self, f: impl FnMut(&dyn RaisedEsurfaceGroupView));
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct ThreeDExpression<O, E = (), H = ()>
where
    O: From<usize> + Into<usize>,
{
    pub orientations: TiVec<O, OrientationExpression>,
    pub surfaces: SurfaceCache<E, H>,
}

pub type CFFExpression<O, E = (), H = ()> = ThreeDExpression<O, E, H>;

impl<O: Clone + From<usize> + Into<usize>, E, H> ThreeDExpression<O, E, H>
where
    usize: From<O>,
{
    pub fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
            surfaces: SurfaceCache::new(),
        }
    }

    pub fn to_atom<P: OrientationSelector>(&self, pattern: P) -> Atom {
        self.orientations
            .iter()
            .filter_map(|orientation| {
                if pattern.filter_orientation(orientation.orientation()) {
                    Some(orientation.to_atom())
                } else {
                    None
                }
            })
            .reduce(|a, b| a + b)
            .unwrap_or_default()
    }

    pub fn get_orientation_atoms<P: OrientationSelector>(
        &self,
        pattern: P,
    ) -> TiVec<OrientationID, Atom> {
        self.orientations
            .iter()
            .map(|orientation| {
                if pattern.filter_orientation(orientation.orientation()) {
                    orientation.to_atom()
                } else {
                    Atom::new()
                }
            })
            .collect()
    }

    pub fn get_orientation_atoms_with_data<P: OrientationSelector>(
        &self,
        pattern: P,
    ) -> TiVec<OrientationID, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = if pattern.filter_orientation(orientation.orientation()) {
                    orientation.to_atom()
                } else {
                    Atom::new()
                };
                let data = orientation.data.clone();
                (atom, data)
            })
            .collect()
    }

    pub fn get_orientation_atom(&self, orientation_id: O) -> Atom {
        self.orientations[orientation_id].to_atom()
    }

    pub fn remap_energy_edge_indices(mut self, edge_map: &EnergyEdgeIndexMap) -> Self {
        for orientation in self.orientations.iter_mut() {
            orientation.remap_energy_edge_indices(edge_map);
        }

        for surface in self.surfaces.linear_surface_cache.iter_mut() {
            surface.expression = surface
                .expression
                .clone()
                .remap_energy_edges(&edge_map.internal, &edge_map.external);
        }

        self
    }

    pub fn num_unfolded_terms(&self) -> usize {
        self.orientations
            .iter()
            .map(OrientationExpression::num_unfolded_terms)
            .sum()
    }

    pub fn get_orientations_with_esurface(&self, esurface_id: EsurfaceID) -> Vec<O> {
        self.orientations
            .iter_enumerated()
            .filter_map(|(id, orientation)| {
                if orientation
                    .iter_denominator_nodes()
                    .any(|node| node.data == HybridSurfaceID::Esurface(esurface_id))
                {
                    Some(id)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn select_esurface_residue(
        mut self,
        raised_esurface_group: &impl RaisedEsurfaceGroupView,
    ) -> Vec<ThreeDExpression<O, E, H>>
    where
        E: Clone,
        H: Clone,
    {
        self.normalize_single_raising(raised_esurface_group);

        let representative_esurface_id = raised_esurface_group.esurface_ids()[0];

        let mut result = vec![];
        for occurrence in 1..=raised_esurface_group.max_occurrence() {
            let mut new_expression = self.clone();

            for orientation in new_expression.orientations.iter_mut() {
                orientation.for_each_denominator_tree_mut(|denominator| {
                    denominator.keep_branches_with_value_count_mut(
                        &HybridSurfaceID::Esurface(representative_esurface_id),
                        occurrence,
                    );

                    denominator.map_mut(|hybrid_surface_id| match hybrid_surface_id {
                        HybridSurfaceID::Esurface(esurface_id)
                            if *esurface_id == representative_esurface_id =>
                        {
                            *hybrid_surface_id = HybridSurfaceID::Unit;
                        }
                        _ => (),
                    });
                });
            }

            result.push(new_expression);
        }

        result
    }

    pub fn normalize_single_raising(
        &mut self,
        raised_esurface_group: &impl RaisedEsurfaceGroupView,
    ) {
        let representative_esurface_id = raised_esurface_group.esurface_ids()[0];

        for orientation in self.orientations.iter_mut() {
            orientation.for_each_denominator_tree_mut(|denominator| {
                denominator.map_mut(|hybrid_surface_id| match hybrid_surface_id {
                    HybridSurfaceID::Esurface(esurface_id)
                        if raised_esurface_group.esurface_ids().contains(esurface_id) =>
                    {
                        *hybrid_surface_id = HybridSurfaceID::Esurface(representative_esurface_id);
                    }
                    _ => (),
                });
            });
        }
    }

    pub fn normalize_wrt_all_raisings(&mut self, raised_data: &impl RaisedEsurfaceDataView) {
        let mut esurface_mappings = HashMap::new();

        raised_data.for_each_raised_group(|cut_group| {
            let esurface_id_of_first = cut_group.esurface_ids()[0];

            for esurface_id in cut_group.esurface_ids().iter() {
                esurface_mappings.insert(*esurface_id, esurface_id_of_first);
            }
        });

        for orientation in self.orientations.iter_mut() {
            orientation.for_each_denominator_tree_mut(|denominator| {
                denominator.map_mut(|hybrid_surface_id| {
                    if let HybridSurfaceID::Esurface(esurface_id) = hybrid_surface_id
                        && let Some(normalized_esurface_id) = esurface_mappings.get(esurface_id)
                    {
                        *esurface_id = *normalized_esurface_id;
                    }
                });
            });
        }
    }
}

impl<O, E, H> Display for ThreeDExpression<O, E, H>
where
    O: Clone + Copy + From<usize> + Into<usize>,
    usize: From<O>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "id | orientation | variant | pref | M power | half edges | den surfaces | num surfaces"
        )?;

        for (orientation_id, orientation) in self.orientations.iter_enumerated() {
            let id: usize = orientation_id.into();
            let label =
                orientation.data.label.clone().unwrap_or_else(|| {
                    format_graph_orientation_label(&orientation.data.orientation)
                });

            for (variant_index, variant) in orientation.variants.iter().enumerate() {
                let half_edges = variant
                    .half_edges
                    .iter()
                    .map(|edge_id| edge_id.0.to_string())
                    .join(", ");
                let denominator_surfaces = variant
                    .denominator
                    .iter_nodes()
                    .map(|node| format!("{:?}", node.data))
                    .unique()
                    .join(", ");
                let numerator_surfaces = variant
                    .numerator_surfaces
                    .iter()
                    .map(|surface_id| format!("{surface_id:?}"))
                    .join(", ");
                writeln!(
                    f,
                    "{} | {} | {} | {} | {} | [{}] | [{}] | [{}]",
                    id,
                    label,
                    variant
                        .origin
                        .clone()
                        .unwrap_or_else(|| format!("variant {variant_index}")),
                    variant.prefactor,
                    variant.uniform_scale_power,
                    half_edges,
                    denominator_surfaces,
                    numerator_surfaces
                )?;
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq)]
pub struct RaisedEsurfaceData {
    pub raised_groups: TiVec<RaisedEsurfaceId, RaisedEsurfaceGroup>,
}

impl RaisedEsurfaceDataView for RaisedEsurfaceData {
    fn for_each_raised_group(&self, mut f: impl FnMut(&dyn RaisedEsurfaceGroupView)) {
        for group in &self.raised_groups {
            f(group);
        }
    }
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Encode, Decode, Hash, PartialOrd, Ord,
)]
pub struct RaisedEsurfaceId(pub usize);

impl From<usize> for RaisedEsurfaceId {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<RaisedEsurfaceId> for usize {
    fn from(value: RaisedEsurfaceId) -> Self {
        value.0
    }
}

#[derive(
    Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord,
)]
pub struct RaisedEsurfaceGroup {
    pub esurface_ids: Vec<EsurfaceID>,
    pub max_occurence: usize,
}

impl RaisedEsurfaceGroupView for RaisedEsurfaceGroup {
    fn esurface_ids(&self) -> &[EsurfaceID] {
        &self.esurface_ids
    }

    fn max_occurrence(&self) -> usize {
        self.max_occurence
    }
}
