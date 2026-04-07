use std::fmt::Display;

use crate::graph::{Graph, LmbIndex, LoopMomentumBasis};
use crate::integrands::process::evaluators::SingleOrAllOrientations;
use crate::momentum::{FourMomentum, Polarization, Rotatable, Rotation, SignOrZero, ThreeMomentum};

use crate::momentum::signature::LoopSignature;
use crate::utils::hyperdual_utils::new_constant;
use crate::utils::{F, FloatLike, Length, PrecisionUpgradable};
use crate::{DependentMomentaConstructor, define_index, settings::runtime::kinematic::Externals};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use derive_more::{From, Into};
use eyre::eyre;
use linnet::half_edge::HedgeGraph;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use linnet::half_edge::nodestore::NodeStorageOps;
use linnet::half_edge::subgraph::subset::SubSet;
use linnet::half_edge::subgraph::{
    Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps,
};
use linnet::half_edge::typed_vec::IndexLike;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Index, IndexMut, Sub};
use symbolica::domains::dual::HyperDual;
use tabled::settings::Style;

use typed_index_collections::TiVec;

#[derive(
    From,
    Into,
    Copy,
    Clone,
    Hash,
    Eq,
    Ord,
    PartialEq,
    PartialOrd,
    bincode::Encode,
    bincode::Decode,
    Debug,
    Serialize,
    Deserialize,
)]
pub struct LoopIndex(pub usize);

// #[derive(
//     From,
//     Into,
//     Copy,
//     Clone,
//     Hash,
//     Eq,
//     Ord,
//     PartialEq,
//     PartialOrd,
//     Encode,
//     Decode,
//     Debug,
//     Serialize,
//     Deserialize,
// )]
// pub struct ExternalIndex(pub usize);

define_index!(
    pub struct ExternalIndex;
);

// define_indexed_vec!(
//     EdgeIndex;

//     pub struct ExternalThr;
// );

impl Display for LoopIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Display for ExternalIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(From, Into, Serialize, Deserialize, Clone, PartialEq, Debug, Encode, Decode)]
pub struct LoopMomenta<T>(pub Vec<ThreeMomentum<T>>);

impl<T: Display> Display for LoopMomenta<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table_builder = tabled::builder::Builder::new();
        for (i, mom) in self.0.iter().enumerate() {
            let mut row = Vec::new();
            row.push(i.to_string());
            for m in mom {
                row.push(m.to_string());
            }
            table_builder.push_record(row);
        }

        write!(f, "{}", table_builder.build().with(Style::modern_rounded()))
    }
}

impl<T> Length for LoopMomenta<T> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }
}

pub type Subspace<'a> = Option<&'a [LoopIndex]>; // None means full space

#[derive(Clone, Debug, bincode::Encode, bincode::Decode)]
pub struct SubspaceData {
    subgraph: InternalSubGraph,
    #[bincode(with_serde)]
    complement_subgraph: SuBitGraph,
    lmb: LmbIndex,
    lmb_indices: Vec<LoopIndex>,
}

impl SubspaceData {
    pub(crate) fn is_mergable_with(&self, other: &Self) -> bool {
        self.lmb == other.lmb
            && self
                .lmb_indices
                .iter()
                .all(|idx| !other.lmb_indices.contains(idx))
            && other
                .lmb_indices
                .iter()
                .all(|idx| !self.lmb_indices.contains(idx))
    }

    fn cleaned_filter_pessimist<E, V, H, N: NodeStorageOps<NodeData = V>>(
        mut filter: SuBitGraph,
        graph: &HedgeGraph<E, V, H, N>,
    ) -> InternalSubGraph {
        let mut to_remove = SuBitGraph::empty(filter.size());

        for i in filter.included_iter() {
            if !filter.includes(&graph.inv(i)) {
                to_remove.add(i);
            }
        }
        filter.subtract_with(&to_remove);

        InternalSubGraph {
            filter,
            loopcount: None,
        }
    }

    pub(crate) fn new_with_user_selected_lmb(
        subgraph: SuBitGraph,
        lmb: LmbIndex,
        graph: &Graph,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> Result<Self> {
        let mut subgraph = Self::cleaned_filter_pessimist(subgraph, graph);
        subgraph.set_loopcount(graph);
        let edges_in_subgraph = graph
            .iter_edges_of(&subgraph.filter)
            .map(|ed| ed.1)
            .collect::<Vec<_>>();
        let complement_subgraph = graph.full_filter().subtract(&subgraph.filter);

        let mut lmb_indices = Vec::new();
        for (loop_index, edge_index) in all_lmbs[lmb].loop_edges.iter_enumerated() {
            if edges_in_subgraph.contains(edge_index) {
                lmb_indices.push(loop_index);
            }
        }

        if lmb_indices.len() != subgraph.loopcount.unwrap() {
            return Err(color_eyre::eyre::eyre!(
                "Provided loop momentum basis does not align with subgraph, lmb_indices: {:?}, lmb_edges: {:?}, subgraph loopcount: {}, edges_in_subgraph: {:?}",
                lmb_indices,
                all_lmbs[lmb].loop_edges,
                subgraph.loopcount.unwrap(),
                edges_in_subgraph,
            ));
        }

        Ok(Self {
            subgraph,
            complement_subgraph,
            lmb,
            lmb_indices,
        })
    }

    /// this function chooses the lmb automatically based on the subgraph
    #[allow(dead_code)]
    pub(crate) fn new(
        subgraph: SuBitGraph,
        graph: &Graph,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> Result<Self> {
        let mut subgraph = InternalSubGraph::cleaned_filter_pessimist(subgraph, graph);
        subgraph.set_loopcount(graph);
        let complement_subgraph = graph.full_filter().subtract(&subgraph.filter);

        let edges_in_subgraph = graph
            .iter_edges_of(&subgraph.filter)
            .map(|ed| ed.1)
            .collect::<Vec<_>>();

        for (lmb_index, lmb) in all_lmbs.iter_enumerated() {
            let mut lmb_indices = Vec::new();
            for (loop_index, edge_index) in lmb.loop_edges.iter_enumerated() {
                if edges_in_subgraph.contains(edge_index) {
                    lmb_indices.push(loop_index);
                }
            }

            if lmb_indices.len() == subgraph.loopcount.unwrap() {
                return Ok(Self {
                    subgraph,
                    complement_subgraph,
                    lmb: lmb_index,
                    lmb_indices,
                });
            }
        }

        Err(eyre!("No suitable loop momentum basis found for subgraph"))
    }
    pub(crate) fn get_lmb<'a>(
        &self,
        all_lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> &'a LoopMomentumBasis {
        &all_lmbs[self.lmb]
    }

    pub(crate) fn does_not_contain<'a>(
        &'a self,
        edges: &'a [EdgeIndex],
        graph: &'a Graph,
    ) -> impl Iterator<Item = EdgeIndex> + 'a {
        graph
            .iter_edges_of(&self.complement_subgraph)
            .map(|ed| ed.1)
            .filter(|e| edges.contains(e))
    }

    pub(crate) fn contains<'a>(
        &'a self,
        edges: &'a [EdgeIndex],
        graph: &'a Graph,
    ) -> impl Iterator<Item = EdgeIndex> + 'a {
        graph
            .iter_edges_of(&self.subgraph)
            .map(|ed| ed.1)
            .filter(|e| edges.contains(e))
    }

    pub(crate) fn project_loop_signature<'a>(
        &'a self,
        signature: &'a LoopSignature,
    ) -> impl Iterator<Item = SignOrZero> + 'a {
        signature.iter_enumerated().map(|(loop_index, sign)| {
            if self.lmb_indices.contains(&loop_index) {
                *sign
            } else {
                SignOrZero::Zero
            }
        })
    }

    pub(crate) fn project_loop_signature_filtered<'a>(
        &'a self,
        signature: &'a LoopSignature,
    ) -> impl Iterator<Item = SignOrZero> + 'a {
        signature
            .iter_enumerated()
            .filter_map(move |(loop_index, sign)| {
                if self.lmb_indices.contains(&loop_index) {
                    Some(*sign)
                } else {
                    None
                }
            })
    }

    pub(crate) fn project_complement_loop_signature<'a>(
        &'a self,
        signature: &'a LoopSignature,
    ) -> impl Iterator<Item = SignOrZero> + 'a {
        signature.iter_enumerated().map(|(loop_index, sign)| {
            if !self.lmb_indices.contains(&loop_index) {
                *sign
            } else {
                SignOrZero::Zero
            }
        })
    }

    pub(crate) fn iter_lmb_indices<'a>(&'a self) -> impl Iterator<Item = LoopIndex> + 'a {
        self.lmb_indices.iter().copied()
    }

    pub(crate) fn contains_loop_index(&self, loop_index: LoopIndex) -> bool {
        self.lmb_indices.contains(&loop_index)
    }

    pub(crate) fn loopcount(&self) -> usize {
        self.subgraph.loopcount.unwrap()
    }

    pub(crate) fn as_subspace_simple(&self) -> Subspace<'_> {
        Some(self.lmb_indices.as_slice())
    }
}

// #[comemo::track]
// impl<T> LoopMomenta<T> {}

impl<T> LoopMomenta<T> {
    pub(crate) fn iter(&self) -> std::slice::Iter<'_, ThreeMomentum<T>> {
        self.0.iter()
    }

    pub(crate) fn first(&self) -> Option<&ThreeMomentum<T>> {
        self.0.first()
    }

    pub(crate) fn iter_enumerated(&self) -> impl Iterator<Item = (LoopIndex, &ThreeMomentum<T>)> {
        self.0.iter().enumerate().map(|(i, m)| (LoopIndex(i), m))
    }
}

impl<T> IntoIterator for LoopMomenta<T> {
    type Item = ThreeMomentum<T>;
    type IntoIter = std::vec::IntoIter<ThreeMomentum<T>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<T> FromIterator<ThreeMomentum<T>> for LoopMomenta<T> {
    fn from_iter<I: IntoIterator<Item = ThreeMomentum<T>>>(iter: I) -> Self {
        LoopMomenta(iter.into_iter().collect())
    }
}

impl<T> Index<LoopIndex> for LoopMomenta<T> {
    type Output = ThreeMomentum<T>;

    fn index(&self, index: LoopIndex) -> &Self::Output {
        &self.0[index.0]
    }
}

impl<T> IndexMut<LoopIndex> for LoopMomenta<T> {
    fn index_mut(&mut self, index: LoopIndex) -> &mut Self::Output {
        &mut self.0[index.0]
    }
}

impl<T: FloatLike> LoopMomenta<F<T>> {
    pub(crate) fn hyper_radius_squared(&self, subspace: Subspace) -> F<T> {
        let zero = self.0[0].px.zero();
        match subspace {
            None => self.iter().fold(zero, |acc, x| acc + x.norm_squared()),
            Some(subspace) => subspace
                .iter()
                .fold(zero, |acc, &i| acc + self[i].norm_squared()),
        }
    }

    pub(crate) fn rescale(&self, factor: &F<T>, subspace: Subspace) -> Self {
        match subspace {
            None => LoopMomenta::from_iter(self.iter().map(|k| k * factor)),
            // this branch is wrong
            Some(subspace) => LoopMomenta::from_iter(self.iter_enumerated().map(|(i, k)| {
                if subspace.contains(&i) {
                    k * factor
                } else {
                    k.clone()
                }
            })),
        }
    }

    pub(crate) fn rescale_with_hyper_dual(
        &self,
        factor: &HyperDual<F<T>>,
        subspace: Subspace,
    ) -> LoopMomenta<HyperDual<F<T>>> {
        match subspace {
            None => LoopMomenta::from_iter(self.iter().map(|k| {
                ThreeMomentum::new(
                    new_constant(factor, &k.px) * factor,
                    new_constant(factor, &k.py) * factor,
                    new_constant(factor, &k.pz) * factor,
                )
            })),
            Some(subspace) => LoopMomenta::from_iter(self.iter_enumerated().map(|(i, k)| {
                if subspace.contains(&i) {
                    ThreeMomentum::new(
                        new_constant(factor, &k.px) * factor,
                        new_constant(factor, &k.py) * factor,
                        new_constant(factor, &k.pz) * factor,
                    )
                } else {
                    ThreeMomentum::new(
                        new_constant(factor, &k.px),
                        new_constant(factor, &k.py),
                        new_constant(factor, &k.pz),
                    )
                }
            })),
        }
    }

    pub(crate) fn rotate(&self, rotation: &Rotation) -> Self {
        LoopMomenta::from_iter(self.iter().map(|k| k.rotate(rotation)))
    }

    pub(crate) fn lmb_transform(
        &self,
        from: &LoopMomentumBasis,
        to: &LoopMomentumBasis,
        externals: &ExternalThreeMomenta<F<T>>,
    ) -> Self {
        LoopMomenta::from_iter(
            to.loop_edges
                .iter()
                .map(|e_id| from.edge_signatures[*e_id].compute_momentum(self, externals)),
        )
    }
}

impl<T: FloatLike> LoopMomenta<HyperDual<F<T>>> {
    pub fn lmb_transform(
        &self,
        from: &LoopMomentumBasis,
        to: &LoopMomentumBasis,
        externals: &ExternalThreeMomenta<HyperDual<F<T>>>,
    ) -> Self {
        LoopMomenta::from_iter(to.loop_edges.iter().map(|e_id| {
            from.edge_signatures[*e_id]
                .try_compute_momentum(&self.0, &externals.raw)
                .unwrap()
        }))
    }

    pub fn rescale(&self, factor: &HyperDual<F<T>>, subspace: Subspace) -> Self {
        match subspace {
            None => LoopMomenta::from_iter(self.iter().map(|k| k * factor)),
            // this branch is wrong
            Some(subspace) => LoopMomenta::from_iter(self.iter_enumerated().map(|(i, k)| {
                if subspace.contains(&i) {
                    k * factor
                } else {
                    k.clone()
                }
            })),
        }
    }
}

impl LoopMomenta<F<f64>> {
    pub(crate) fn cast<T: FloatLike>(&self) -> LoopMomenta<F<T>> {
        LoopMomenta::from_iter(self.iter().map(|m| m.map(&|x| F::from_ff64(x))))
    }
}

impl<T: FloatLike> Sub<&LoopMomenta<F<T>>> for &LoopMomenta<F<T>> {
    type Output = LoopMomenta<F<T>>;

    fn sub(self, rhs: &LoopMomenta<F<T>>) -> Self::Output {
        LoopMomenta::from_iter(self.iter().zip(rhs.iter()).map(|(l, r)| l - r))
    }
}

impl<T: FloatLike> Add<&LoopMomenta<F<T>>> for &LoopMomenta<F<T>> {
    type Output = LoopMomenta<F<T>>;

    fn add(self, rhs: &LoopMomenta<F<T>>) -> Self::Output {
        LoopMomenta::from_iter(self.iter().zip(rhs.iter()).map(|(l, r)| l + r))
    }
}

// define_indexed_vec!(
//     MyIdx;
//     pub struct ExternalMomentumtwo<ThreeMomentum<T>>;
// );
pub type ExternalThreeMomenta<T> = TiVec<ExternalIndex, ThreeMomentum<T>>;
// define_indexed_vec!()
pub type ExternalFourMomenta<T> = TiVec<ExternalIndex, FourMomentum<T>>;
// impl<T: Display> Display for ExternalFourMomenta<T> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         let mut table_builder = tabled::builder::Builder::new();
//         for (i, mom) in self.iter_enumerated() {
//             let mut row = Vec::new();
//             row.push(i.to_string());
//             for m in mom {
//                 row.push(m.to_string());
//             }
//             table_builder.push_record(row);
//         }

//         write!(f, "{}", table_builder.build().with(Style::modern_rounded()))
//     }
// }
pub type PolarizationVectors<T> = TiVec<ExternalIndex, Polarization<T>>; // should be the same length as #externals

fn extract_external_spatial<T: Clone>(
    external_four_momenta: &ExternalFourMomenta<T>,
) -> ExternalThreeMomenta<T> {
    ExternalThreeMomenta::from_iter(external_four_momenta.iter().map(|fm| fm.spatial.clone()))
}

#[derive(Debug, Clone)]
pub struct BareMomentumSample<T: FloatLike> {
    pub loop_moms: LoopMomenta<F<T>>,
    pub dual_loop_moms: Option<LoopMomenta<HyperDual<F<T>>>>,
    pub loop_mom_cache_id: usize,
    /// Base cache ID for the fundamental loop momentum configuration (before transformations)
    pub loop_mom_base_cache_id: usize,
    pub external_moms: ExternalFourMomenta<F<T>>,
    pub external_mom_cache_id: usize,
    /// Base cache ID for the fundamental external momentum configuration (before transformations)
    pub external_mom_base_cache_id: usize,
    pub jacobian: F<T>,
    pub orientation: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct MomentumSample<T: FloatLike> {
    pub sample: BareMomentumSample<T>,
    // pub uuid: Uuid,
}

impl<T: FloatLike> Display for MomentumSample<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = tabled::builder::Builder::new();

        table.push_record(["Sample"]);

        table.push_record(["Loop Momenta", "p_x", "p_y", "p_z"]);
        // table

        for (index, loop_mom) in self.loop_moms().0.iter().enumerate() {
            table.push_record([
                index.to_string(),
                loop_mom.px.to_string(),
                loop_mom.py.to_string(),
                loop_mom.pz.to_string(),
            ]);
        }

        table.push_record(["External Momenta", "E", "p_x", "p_y", "p_z"]);
        for (index, external_mom) in self.external_moms().iter_enumerated() {
            table.push_record([
                index.to_string(),
                external_mom.temporal.to_string(),
                external_mom.spatial.px.to_string(),
                external_mom.spatial.py.to_string(),
                external_mom.spatial.pz.to_string(),
            ]);
        }

        table.push_record(["Jacobian".into(), format!("{:+e}", self.sample.jacobian)]);
        table.build().with(Style::rounded()).fmt(f)
    }
}

impl<T: FloatLike> BareMomentumSample<T> {
    #[inline(never)]
    pub(crate) fn new(
        loop_moms: LoopMomenta<F<T>>,
        loop_mom_cache_id: usize,
        external_moms: &Externals,
        external_mom_cache_id: usize,
        jacobian: F<T>,
        dependent_momenta_constructor: DependentMomentaConstructor,
        orientation: Option<usize>,
    ) -> Result<Self> {
        let external_moms = external_moms.get_dependent_externals(dependent_momenta_constructor)?;

        Ok(Self {
            loop_moms,
            dual_loop_moms: None,
            loop_mom_cache_id,
            loop_mom_base_cache_id: loop_mom_cache_id, // Initially same as cache_id
            external_mom_cache_id,
            external_mom_base_cache_id: external_mom_cache_id, // Initially same as cache_id
            external_moms,
            jacobian,
            orientation,
        })
    }

    pub(crate) fn one(&self) -> F<T> {
        if let Some(f) = self.loop_moms.first() {
            f.px.one()
        } else if let Some(f) = self.external_moms.first() {
            f.spatial.px.one()
        } else {
            panic!("No momenta in sample")
        }
    }

    pub(crate) fn zero(&self) -> F<T> {
        if let Some(f) = self.loop_moms.first() {
            f.px.zero()
        } else if let Some(f) = self.external_moms.first() {
            f.spatial.px.zero()
        } else {
            panic!("No momenta in sample")
        }
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> BareMomentumSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        BareMomentumSample {
            loop_mom_cache_id: self.loop_mom_cache_id,
            loop_mom_base_cache_id: self.loop_mom_base_cache_id,
            external_mom_cache_id: self.external_mom_cache_id,
            external_mom_base_cache_id: self.external_mom_base_cache_id,
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::cast).collect(),
            dual_loop_moms: self
                .dual_loop_moms
                .as_ref()
                .map(|_dlm| todo!("make sure the cast works if there are hyperdual momenta")),
            external_moms: self.external_moms.iter().map(FourMomentum::cast).collect(),
            jacobian: self.jacobian.clone().into(),
            orientation: self.orientation,
        }
    }

    pub(crate) fn higher_precision(&self) -> BareMomentumSample<T::Higher>
    where
        T::Higher: FloatLike,
        T::Lower: FloatLike,
    {
        BareMomentumSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::higher).collect(),
            dual_loop_moms: self.dual_loop_moms.as_ref().map(|dlm| {
                LoopMomenta::from_iter(dlm.iter().map(|m| m.clone().map(&|x| x.higher())))
            }),
            external_moms: self
                .external_moms
                .iter()
                .map(FourMomentum::higher)
                .collect(),
            loop_mom_cache_id: self.loop_mom_cache_id,
            loop_mom_base_cache_id: self.loop_mom_base_cache_id,
            external_mom_cache_id: self.external_mom_cache_id,
            external_mom_base_cache_id: self.external_mom_base_cache_id,
            jacobian: self.jacobian.higher(),
            orientation: self.orientation,
        }
    }

    pub(crate) fn lower_precision(&self) -> BareMomentumSample<T::Lower>
    where
        T::Higher: FloatLike,
        T::Lower: FloatLike,
    {
        BareMomentumSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::lower).collect(),
            dual_loop_moms: self.dual_loop_moms.as_ref().map(|dlm| {
                LoopMomenta::from_iter(dlm.iter().map(|m| m.clone().map(&|x| x.lower())))
            }),
            external_moms: self.external_moms.iter().map(FourMomentum::lower).collect(),
            jacobian: self.jacobian.lower(),
            orientation: self.orientation,
            loop_mom_cache_id: self.loop_mom_cache_id,
            loop_mom_base_cache_id: self.loop_mom_base_cache_id,
            external_mom_cache_id: self.external_mom_cache_id,
            external_mom_base_cache_id: self.external_mom_base_cache_id,
        }
    }

    #[inline]
    pub(crate) fn rotate(
        &self,
        rotation: &Rotation,
        loop_mom_cache_id: usize,
        external_mom_cache_id: usize,
    ) -> Self {
        Self {
            loop_moms: self.loop_moms.iter().map(|l| l.rotate(rotation)).collect(),
            dual_loop_moms: self
                .dual_loop_moms
                .as_ref()
                .map(|dlm| LoopMomenta::from_iter(dlm.iter().map(|l| l.rotate(rotation)))),
            external_moms: self
                .external_moms
                .iter()
                .map(|l| l.rotate(rotation))
                .collect(),
            loop_mom_cache_id,
            loop_mom_base_cache_id: self.loop_mom_base_cache_id, // Preserve base cache ID
            external_mom_cache_id,
            external_mom_base_cache_id: self.external_mom_base_cache_id, // Preserve base cache ID
            jacobian: self.jacobian.clone(),
            orientation: self.orientation,
        }
    }

    #[inline]
    #[allow(dead_code)]
    pub(crate) fn rescaled_loop_momenta(&self, factor: &F<T>, subspace: Subspace) -> Self {
        Self {
            loop_moms: self.loop_moms.rescale(factor, subspace),
            dual_loop_moms: self
                .dual_loop_moms
                .as_ref()
                .map(|dlm| dlm.rescale(&new_constant(&dlm[LoopIndex(0)].px, factor), subspace)),
            loop_mom_cache_id: self.loop_mom_cache_id + 1,
            loop_mom_base_cache_id: self.loop_mom_base_cache_id, // Preserve base cache ID
            external_moms: self.external_moms.clone(),
            external_mom_cache_id: self.external_mom_cache_id,
            external_mom_base_cache_id: self.external_mom_base_cache_id, // Preserve base cache ID
            jacobian: self.jacobian.clone(),
            orientation: self.orientation,
        }
    }

    #[inline]
    pub(crate) fn lmb_transform(&self, from: &LoopMomentumBasis, to: &LoopMomentumBasis) -> Self {
        Self {
            loop_moms: self.loop_moms.lmb_transform(
                from,
                to,
                &extract_external_spatial(&self.external_moms),
            ),
            dual_loop_moms: self.dual_loop_moms.as_ref().map(|dlm| {
                let dual_externals = self
                    .external_moms
                    .iter()
                    .map(|four_mom| {
                        ThreeMomentum::new(
                            new_constant(&dlm[LoopIndex(0)].px, &four_mom.spatial.px),
                            new_constant(&dlm[LoopIndex(0)].px, &four_mom.spatial.py),
                            new_constant(&dlm[LoopIndex(0)].px, &four_mom.spatial.pz),
                        )
                    })
                    .collect();
                dlm.lmb_transform(from, to, &dual_externals)
            }),
            loop_mom_cache_id: self.loop_mom_cache_id + 1,
            loop_mom_base_cache_id: self.loop_mom_base_cache_id, // Preserve base cache ID
            external_moms: self.external_moms.clone(),
            external_mom_cache_id: self.external_mom_cache_id,
            external_mom_base_cache_id: self.external_mom_base_cache_id, // Preserve base cache ID
            jacobian: self.jacobian.clone(),
            orientation: self.orientation,
        }
    }
}

impl<T: FloatLike> MomentumSample<T> {
    pub(crate) fn orientations<'a, OID: From<usize> + Copy + IndexLike>(
        &self,
        filter: &'a SubSet<OID>,
        orientations: &'a TiVec<OID, EdgeVec<Orientation>>,
    ) -> SingleOrAllOrientations<'a, OID>
    where
        usize: From<OID>,
    {
        if let Some(o) = self.sample.orientation {
            let id = if filter.is_full() {
                OID::from(o)
            } else {
                filter
                    .included_iter()
                    .nth(o)
                    .unwrap_or_else(|| {
                        panic!(
                            "sampled orientation index {o} must resolve within the filtered orientation subset"
                        )
                    })
            };
            SingleOrAllOrientations::Single {
                id,
                orientation: &orientations[id],
            }
        } else {
            SingleOrAllOrientations::All {
                all: orientations,
                filter,
            }
        }
    }

    // pub(crate) fn uuid(&self) -> Option<Uuid> {
    //     if self.rotated_sample.is_some() {
    //         None
    //     } else {
    //         Some(self.uuid)
    //     }
    // }

    pub(crate) fn loop_moms(&self) -> &LoopMomenta<F<T>> {
        &self.sample.loop_moms
    }

    pub(crate) fn external_moms(&self) -> &ExternalFourMomenta<F<T>> {
        &self.sample.external_moms
    }

    pub(crate) fn jacobian(&self) -> F<T> {
        self.sample.jacobian.clone()
    }

    pub(crate) fn new(
        loop_moms: LoopMomenta<F<T>>,
        loop_mom_cache_id: usize,
        external_moms: &Externals,
        external_mom_cache_id: usize,
        jacobian: F<T>,
        dependent_momenta_constructor: DependentMomentaConstructor,
        orientation: Option<usize>,
    ) -> Result<Self> {
        // info!("New with ext cache id:{external_mom_cache_id}");
        Ok(Self {
            sample: BareMomentumSample::new(
                loop_moms,
                loop_mom_cache_id,
                external_moms,
                external_mom_cache_id,
                jacobian,
                dependent_momenta_constructor,
                orientation,
            )?,
            // uuid: Uuid::new_v4(),
        })
    }

    pub(crate) fn one(&self) -> F<T> {
        self.sample.one()
    }

    pub(crate) fn zero(&self) -> F<T> {
        self.sample.zero()
    }

    /// Cast the sample to a different precision
    #[inline]
    pub(crate) fn cast_sample<T2: FloatLike>(&self) -> MomentumSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        MomentumSample {
            sample: self.sample.cast_sample(),
        }
    }

    pub(crate) fn higher_precision(&self) -> MomentumSample<T::Higher>
    where
        T::Higher: FloatLike + Default,
        T::Lower: FloatLike + Default,
    {
        MomentumSample {
            sample: self.sample.higher_precision(),
        }
    }

    pub(crate) fn lower_precision(&self) -> MomentumSample<T::Lower>
    where
        T::Lower: FloatLike + Default,
        T::Higher: FloatLike + Default,
    {
        MomentumSample {
            sample: self.sample.lower_precision(),
        }
    }

    #[inline]
    pub(crate) fn rotate(
        &self,
        rotation: &Rotation,
        loop_mom_cache_id: usize,
        external_mom_cache_id: usize,
    ) -> Self {
        Self {
            sample: self
                .sample
                .rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
        }
    }

    #[allow(dead_code)]
    #[inline]
    pub(crate) fn rescaled_loop_momenta(&self, factor: &F<T>, subspace: Subspace) -> Self {
        Self {
            sample: self.sample.rescaled_loop_momenta(factor, subspace),
        }
    }

    #[inline]
    pub(crate) fn lmb_transform(&self, from: &LoopMomentumBasis, to: &LoopMomentumBasis) -> Self {
        Self {
            sample: self.sample.lmb_transform(from, to),
        }
    }
}
