use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{CffError, EdgeId, VertexId};

macro_rules! index_type {
    ($name:ident) => {
        #[derive(
            Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
        )]
        pub struct $name(pub usize);

        impl $name {
            pub const fn index(self) -> usize {
                self.0
            }
        }
    };
}

index_type!(EnergySurfaceId);
index_type!(HSurfaceId);

/// A canonical linear combination of external edge energies.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
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

    pub fn iter(&self) -> impl Iterator<Item = (EdgeId, i64)> + '_ {
        self.0
            .iter()
            .map(|(edge, coefficient)| (*edge, *coefficient))
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
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
                .checked_mul(replacement_coefficient)
                .ok_or(CffError::ExternalShiftOverflow { edge })?;
            self.add(edge, term)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
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
}

/// Set of original vertices represented by a contracted CFF vertex.
#[derive(
    Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
)]
pub struct VertexSet(u64);

impl VertexSet {
    pub(crate) fn singleton(vertex: VertexId) -> Self {
        Self(1_u64 << vertex.index())
    }

    pub(crate) const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    pub fn contains(self, vertex: VertexId) -> bool {
        vertex.index() < VertexId::MAX_COUNT && self.0 & (1_u64 << vertex.index()) != 0
    }

    pub fn iter(self) -> impl Iterator<Item = VertexId> {
        (0..VertexId::MAX_COUNT)
            .filter(move |index| self.0 & (1_u64 << index) != 0)
            .map(VertexId::new)
    }
}

/// An E-surface generated from an ordinary source or sink.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnergySurface {
    pub energies: Vec<EdgeId>,
    pub external_shift: ExternalShift,
    pub vertices: VertexSet,
}

impl EnergySurface {
    pub(crate) fn cache_equivalent(&self, other: &Self) -> bool {
        self.energies == other.energies && self.external_shift == other.external_shift
    }
}

/// An H-surface generated when a subgraph has severed internal edges.
///
/// This is the historical CFF “H-surface” name: in an amplitude subgraph it
/// usually represents an E-surface of the containing graph.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HSurface {
    pub positive_energies: Vec<EdgeId>,
    pub negative_energies: Vec<EdgeId>,
    pub external_shift: ExternalShift,
    pub vertices: VertexSet,
}

impl HSurface {
    pub(crate) fn cache_equivalent(&self, other: &Self) -> bool {
        // Preserve the CFF cache identity: the signed on-shell energies identify
        // an H-surface; the external shift and originating vertices are metadata.
        self.positive_energies == other.positive_energies
            && self.negative_energies == other.negative_energies
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Surface {
    Energy(EnergySurface),
    H(HSurface),
    Unit,
    Infinite,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum SurfaceId {
    Energy(EnergySurfaceId),
    H(HSurfaceId),
    Unit,
    Infinite,
}

/// Surface interning owned by one [`crate::CffGenerator::generate`] call.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SurfaceCache {
    energy_surfaces: Vec<EnergySurface>,
    h_surfaces: Vec<HSurface>,
}

impl SurfaceCache {
    pub fn energy_surfaces(&self) -> &[EnergySurface] {
        &self.energy_surfaces
    }

    pub fn h_surfaces(&self) -> &[HSurface] {
        &self.h_surfaces
    }

    pub fn get(&self, id: SurfaceId) -> Option<Surface> {
        match id {
            SurfaceId::Energy(id) => self
                .energy_surfaces
                .get(id.index())
                .cloned()
                .map(Surface::Energy),
            SurfaceId::H(id) => self.h_surfaces.get(id.index()).cloned().map(Surface::H),
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

    pub(crate) fn intern(&mut self, surface: Surface) -> SurfaceId {
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
}
