use std::collections::BTreeMap;

use linnet::half_edge::{
    HedgeGraph,
    involution::{Flow, HedgePair},
    nodestore::NodeStorageOps,
    subgraph::{OrientedCut, SubGraphLike},
};
use serde::{Deserialize, Serialize};
use typed_index_collections::TiVec;

use crate::{CffError, EdgeId, VertexId};

pub(crate) const MAX_VERTEX_COUNT: usize = u64::BITS as usize;

macro_rules! index_type {
    ($name:ident) => {
        #[derive(
            Clone,
            Copy,
            Debug,
            PartialEq,
            Eq,
            PartialOrd,
            Ord,
            Hash,
            Serialize,
            Deserialize,
            bincode::Encode,
            bincode::Decode,
        )]
        pub struct $name(pub usize);

        impl $name {
            pub const fn index(self) -> usize {
                self.0
            }
        }

        impl From<usize> for $name {
            fn from(value: usize) -> Self {
                Self(value)
            }
        }

        impl From<$name> for usize {
            fn from(value: $name) -> Self {
                value.0
            }
        }
    };
}

index_type!(EnergySurfaceId);
index_type!(HSurfaceId);
index_type!(RaisedEnergySurfaceId);

/// One equivalence class of energy surfaces produced by raised propagators.
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct RaisedEnergySurfaceGroup {
    pub surface_ids: Vec<EnergySurfaceId>,
    pub max_occurrence: usize,
}

/// Structural raised-surface information for one CFF expression.
#[derive(
    Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode,
)]
pub struct RaisedEnergySurfaceData {
    #[bincode(with_serde)]
    pub groups: TiVec<RaisedEnergySurfaceId, RaisedEnergySurfaceGroup>,
}

/// A canonical linear combination of external edge energies.
#[derive(
    Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode,
)]
pub struct ExternalShift(BTreeMap<EdgeId, i64>);

impl ExternalShift {
    pub fn new(terms: impl IntoIterator<Item = (EdgeId, i64)>) -> Result<Self, CffError> {
        let mut shift = Self::default();
        for (edge, coefficient) in terms {
            shift.add(edge, coefficient)?;
        }
        Ok(shift)
    }

    pub fn get(&self, edge: EdgeId) -> i64 {
        self.0.get(&edge).copied().unwrap_or(0)
    }

    pub fn iter(&self) -> impl Iterator<Item = (&EdgeId, &i64)> {
        self.0.iter()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn remove(&mut self, edge: EdgeId) -> Option<i64> {
        self.0.remove(&edge)
    }

    pub(crate) fn add(&mut self, edge: EdgeId, coefficient: i64) -> Result<(), CffError> {
        if coefficient == 0 {
            return Ok(());
        }
        let new_coefficient = self
            .get(edge)
            .checked_add(coefficient)
            .ok_or(CffError::ExternalShiftOverflow { edge })?;
        if new_coefficient == 0 {
            self.0.remove(&edge);
        } else {
            self.0.insert(edge, new_coefficient);
        }
        Ok(())
    }

    pub(crate) fn rewrite(
        &mut self,
        dependent: EdgeId,
        replacement: &Self,
    ) -> Result<(), CffError> {
        let coefficient = self.get(dependent);
        self.0.remove(&dependent);
        for (edge, replacement_coefficient) in replacement.iter() {
            let term = coefficient
                .checked_mul(*replacement_coefficient)
                .ok_or(CffError::ExternalShiftOverflow { edge: *edge })?;
            self.add(*edge, term)?;
        }
        Ok(())
    }
}

impl FromIterator<(EdgeId, i64)> for ExternalShift {
    fn from_iter<T: IntoIterator<Item = (EdgeId, i64)>>(iter: T) -> Self {
        Self::new(iter).expect("external-shift coefficients overflowed i64")
    }
}

impl From<Vec<(EdgeId, i64)>> for ExternalShift {
    fn from(value: Vec<(EdgeId, i64)>) -> Self {
        value.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::{
        HedgeGraph,
        builder::HedgeGraphBuilder,
        involution::{EdgeIndex, Flow, HedgePair, Orientation},
        subgraph::{ModifySubSet, OrientedCut, SuBitGraph},
    };

    use super::*;

    #[test]
    fn duplicate_terms_report_coefficient_overflow() {
        let edge = EdgeId::new(0);
        assert_eq!(
            ExternalShift::new([(edge, i64::MAX), (edge, 1)]),
            Err(CffError::ExternalShiftOverflow { edge })
        );
    }

    #[test]
    fn rewrite_reports_multiplication_overflow() {
        let dependent = EdgeId::new(0);
        let replacement_edge = EdgeId::new(1);
        let mut shift = ExternalShift::new([(dependent, i64::MAX)]).unwrap();
        let replacement = ExternalShift::new([(replacement_edge, 2)]).unwrap();

        assert_eq!(
            shift.rewrite(dependent, &replacement),
            Err(CffError::ExternalShiftOverflow {
                edge: replacement_edge,
            })
        );
    }

    #[test]
    fn cut_side_surface_preserves_original_edges_vertices_and_external_flow() {
        let mut builder = HedgeGraphBuilder::<(), (), ()>::new();
        let left = builder.add_node(());
        let right = builder.add_node(());
        builder.add_external_edge(left, (), Orientation::Undirected, Flow::Sink);
        builder.add_edge(left, right, (), Orientation::Default);
        builder.add_external_edge(right, (), Orientation::Undirected, Flow::Source);
        let graph: HedgeGraph<(), (), ()> = builder.into();

        let mut side: SuBitGraph = graph.empty_subgraph();
        let (_, neighbours, _) = graph
            .iter_nodes()
            .find(|(vertex, _, _)| *vertex == left)
            .unwrap();
        neighbours.for_each(|hedge| side.add(hedge));

        let HedgePair::Paired { source, .. } = graph[&EdgeIndex(1)].1 else {
            panic!("the middle edge must be internal");
        };
        let mut crossing: SuBitGraph = graph.empty_subgraph();
        crossing.add(source);
        let crossing = OrientedCut::from_underlying_strict(crossing, &graph).unwrap();

        let surface = EnergySurface::from_cut_side(&graph, &crossing, &side, None).unwrap();
        assert_eq!(surface.energies, vec![EdgeIndex(1)]);
        assert_eq!(surface.external_shift.get(EdgeIndex(0)), -1);
        assert!(surface.vertex_set.contains(left));
        assert!(!surface.vertex_set.contains(right));
    }
}

/// Set of original vertices represented by a contracted CFF vertex.
#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct VertexSet(u64);

impl VertexSet {
    pub fn singleton(vertex: VertexId) -> Self {
        Self(1_u64 << vertex.index())
    }

    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    pub fn from_usize(vertex: usize) -> Self {
        Self::singleton(VertexId::new(vertex))
    }

    pub const fn dummy() -> Self {
        Self(0)
    }

    pub const fn join(&self, other: &Self) -> Self {
        self.union(*other)
    }

    pub fn contains(self, vertex: VertexId) -> bool {
        vertex.index() < MAX_VERTEX_COUNT && self.0 & (1_u64 << vertex.index()) != 0
    }

    pub fn iter(self) -> impl Iterator<Item = VertexId> {
        (0..MAX_VERTEX_COUNT)
            .filter(move |index| self.0 & (1_u64 << index) != 0)
            .map(VertexId::new)
    }
}

/// An E-surface generated from an ordinary source or sink.
#[derive(Clone, Debug, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct EnergySurface {
    pub energies: Vec<EdgeId>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
}

impl EnergySurface {
    /// Construct the energy surface associated with one side of an oriented cut.
    ///
    /// The crossing edges supply the positive on-shell energies. External shifts
    /// are read either from an explicit initial-state cut or from the unpaired
    /// edges incident on `side`. Original Linnet vertex identities are retained
    /// in the resulting surface.
    pub fn from_cut_side<E, V, H, N, C, S>(
        graph: &HedgeGraph<E, V, H, N>,
        crossing: &C,
        side: &S,
        initial_state_cut: Option<&OrientedCut>,
    ) -> Result<Self, CffError>
    where
        N: NodeStorageOps<NodeData = V>,
        C: SubGraphLike,
        S: SubGraphLike,
    {
        let mut energies = graph
            .iter_edges_of(crossing)
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        energies.sort_unstable();

        let mut external_shift = if let Some(initial_state_cut) = initial_state_cut {
            graph
                .iter_edges_of(initial_state_cut)
                .map(|(_, edge, _)| (edge, -1))
                .collect::<Vec<_>>()
        } else {
            graph
                .iter_edges_of(side)
                .filter_map(|(pair, edge, _)| match pair {
                    HedgePair::Unpaired { flow, .. } => Some((
                        edge,
                        match flow {
                            Flow::Sink => -1,
                            Flow::Source => 1,
                        },
                    )),
                    _ => None,
                })
                .collect::<Vec<_>>()
        };
        external_shift.sort_by_key(|(edge, _)| *edge);

        let mut vertices = graph.iter_nodes_of(side).map(|(vertex, _, _)| vertex);
        let first = vertices.next().ok_or(CffError::EmptyGraph)?;
        if first.0 >= MAX_VERTEX_COUNT {
            return Err(CffError::VertexIdentityOutOfRange {
                vertex: first,
                maximum: MAX_VERTEX_COUNT - 1,
            });
        }
        let mut vertex_set = VertexSet::singleton(first);
        for vertex in vertices {
            if vertex.0 >= MAX_VERTEX_COUNT {
                return Err(CffError::VertexIdentityOutOfRange {
                    vertex,
                    maximum: MAX_VERTEX_COUNT - 1,
                });
            }
            vertex_set = vertex_set.union(VertexSet::singleton(vertex));
        }

        Ok(Self {
            energies,
            external_shift: ExternalShift::new(external_shift)?,
            vertex_set,
        })
    }

    pub(crate) fn cache_equivalent(&self, other: &Self) -> bool {
        self.energies == other.energies && self.external_shift == other.external_shift
    }
}

impl PartialEq for EnergySurface {
    fn eq(&self, other: &Self) -> bool {
        self.cache_equivalent(other)
    }
}

impl Eq for EnergySurface {}

/// An H-surface generated when a subgraph has severed internal edges.
///
/// This is the historical CFF “H-surface” name: in an amplitude subgraph it
/// usually represents an E-surface of the containing graph.
#[derive(Clone, Debug, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct HSurface {
    pub positive_energies: Vec<EdgeId>,
    pub negative_energies: Vec<EdgeId>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
}

impl PartialEq for HSurface {
    fn eq(&self, other: &Self) -> bool {
        self.cache_equivalent(other)
    }
}

impl Eq for HSurface {}

impl HSurface {
    pub(crate) fn cache_equivalent(&self, other: &Self) -> bool {
        // Preserve the CFF cache identity: the signed on-shell energies identify
        // an H-surface; the external shift and originating vertices are metadata.
        self.positive_energies == other.positive_energies
            && self.negative_energies == other.negative_energies
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub enum Surface {
    Energy(EnergySurface),
    H(HSurface),
    Unit,
    Infinite,
}

#[derive(
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub enum SurfaceId {
    Energy(EnergySurfaceId),
    H(HSurfaceId),
    Unit,
    Infinite,
}

/// Surface interning owned by one [`crate::CffGenerator::generate`] call.
#[derive(
    Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode,
)]
pub struct SurfaceCache {
    #[bincode(with_serde)]
    energy_surfaces: TiVec<EnergySurfaceId, EnergySurface>,
    #[bincode(with_serde)]
    h_surfaces: TiVec<HSurfaceId, HSurface>,
}

/// ID translation produced while merging one surface arena into another.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SurfaceIdMap {
    energy: Vec<EnergySurfaceId>,
    h: Vec<HSurfaceId>,
}

impl SurfaceIdMap {
    pub fn remap(&self, id: SurfaceId) -> Option<SurfaceId> {
        match id {
            SurfaceId::Energy(id) => self.energy.get(id.index()).copied().map(SurfaceId::Energy),
            SurfaceId::H(id) => self.h.get(id.index()).copied().map(SurfaceId::H),
            SurfaceId::Unit => Some(SurfaceId::Unit),
            SurfaceId::Infinite => Some(SurfaceId::Infinite),
        }
    }
}

impl SurfaceCache {
    pub fn energy_surfaces(&self) -> &TiVec<EnergySurfaceId, EnergySurface> {
        &self.energy_surfaces
    }

    pub fn h_surfaces(&self) -> &TiVec<HSurfaceId, HSurface> {
        &self.h_surfaces
    }

    pub fn get(&self, id: SurfaceId) -> Option<Surface> {
        match id {
            SurfaceId::Energy(id) => self.energy_surfaces.get(id).cloned().map(Surface::Energy),
            SurfaceId::H(id) => self.h_surfaces.get(id).cloned().map(Surface::H),
            SurfaceId::Unit => Some(Surface::Unit),
            SurfaceId::Infinite => Some(Surface::Infinite),
        }
    }

    pub fn len(&self) -> usize {
        self.energy_surfaces.len() + self.h_surfaces.len()
    }

    pub fn is_empty(&self) -> bool {
        self.energy_surfaces.is_empty() && self.h_surfaces.is_empty()
    }

    pub fn intern(&mut self, surface: Surface) -> SurfaceId {
        match surface {
            Surface::Energy(surface) => {
                if let Some(index) = self
                    .energy_surfaces
                    .iter()
                    .position(|cached| cached.cache_equivalent(&surface))
                {
                    SurfaceId::Energy(EnergySurfaceId(index))
                } else {
                    let id = EnergySurfaceId(self.energy_surfaces.len());
                    self.energy_surfaces.push(surface);
                    SurfaceId::Energy(id)
                }
            }
            Surface::H(surface) => {
                if let Some(index) = self
                    .h_surfaces
                    .iter()
                    .position(|cached| cached.cache_equivalent(&surface))
                {
                    SurfaceId::H(HSurfaceId(index))
                } else {
                    let id = HSurfaceId(self.h_surfaces.len());
                    self.h_surfaces.push(surface);
                    SurfaceId::H(id)
                }
            }
            Surface::Unit => SurfaceId::Unit,
            Surface::Infinite => SurfaceId::Infinite,
        }
    }

    /// Merge `other` into this arena while preserving this arena's existing IDs.
    pub fn merge(&mut self, other: &Self) -> SurfaceIdMap {
        let energy = other
            .energy_surfaces
            .iter()
            .cloned()
            .map(|surface| match self.intern(Surface::Energy(surface)) {
                SurfaceId::Energy(id) => id,
                _ => unreachable!("interning an energy surface returned a different kind"),
            })
            .collect();
        let h = other
            .h_surfaces
            .iter()
            .cloned()
            .map(|surface| match self.intern(Surface::H(surface)) {
                SurfaceId::H(id) => id,
                _ => unreachable!("interning an H-surface returned a different kind"),
            })
            .collect();
        SurfaceIdMap { energy, h }
    }
}

impl From<EnergySurfaceId> for symbolica::atom::Atom {
    fn from(id: EnergySurfaceId) -> Self {
        symbolica::parse!(&format!("η({})", id.index()))
    }
}

impl From<HSurfaceId> for symbolica::atom::Atom {
    fn from(id: HSurfaceId) -> Self {
        symbolica::parse!(&format!("H({})", id.index()))
    }
}

impl From<SurfaceId> for symbolica::atom::Atom {
    fn from(id: SurfaceId) -> Self {
        match id {
            SurfaceId::Energy(id) => id.into(),
            SurfaceId::H(id) => id.into(),
            SurfaceId::Unit => Self::num(1),
            SurfaceId::Infinite => symbolica::parse!("η_inf"),
        }
    }
}

#[cfg(test)]
mod cache_tests {
    use super::*;

    fn energy(edge: usize) -> EnergySurface {
        EnergySurface {
            energies: vec![EdgeId::new(edge)],
            external_shift: ExternalShift::default(),
            vertex_set: VertexSet::default(),
        }
    }

    #[test]
    fn merge_deduplicates_surfaces_and_returns_a_total_remap() {
        let mut target = SurfaceCache::default();
        let existing = target.intern(Surface::Energy(energy(1)));
        let mut source = SurfaceCache::default();
        let duplicate = source.intern(Surface::Energy(energy(1)));
        let distinct = source.intern(Surface::Energy(energy(2)));

        let remap = target.merge(&source);

        assert_eq!(remap.remap(duplicate), Some(existing));
        assert_eq!(
            remap.remap(distinct),
            Some(SurfaceId::Energy(EnergySurfaceId(1)))
        );
        assert_eq!(target.energy_surfaces().len(), 2);
    }
}
