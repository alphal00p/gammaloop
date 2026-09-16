use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};
use typed_index_collections::TiVec;

use crate::{
    CffError, EdgeId, EnergySurfaceId, ExpressionTree, GraphOrientation, OrientationId,
    RaisedEnergySurfaceData, RaisedEnergySurfaceGroup, SurfaceCache, SurfaceId, SurfaceIdMap,
};

/// One CFF denominator tree and its acyclic graph orientation.
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct OrientationData {
    pub orientation: GraphOrientation,
}

/// One CFF denominator tree and its acyclic graph orientation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct OrientationExpression {
    pub id: OrientationId,
    pub data: OrientationData,
    pub expression: ExpressionTree<SurfaceId>,
}

impl OrientationExpression {
    /// Unfold the tree into denominator products, one vector per additive term.
    pub fn denominator_products(&self) -> Vec<Vec<SurfaceId>> {
        self.expression
            .root_to_leaf_paths()
            .into_iter()
            .map(|path| path.into_iter().copied().collect())
            .collect()
    }
}

/// Sum of the expressions associated with all accepted acyclic orientations.
#[derive(
    Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode,
)]
pub struct CffExpression {
    #[bincode(with_serde)]
    pub orientations: TiVec<OrientationId, OrientationExpression>,
}

impl CffExpression {
    pub(crate) fn new(orientations: Vec<OrientationExpression>) -> Self {
        Self {
            orientations: orientations.into(),
        }
    }

    pub fn orientations(&self) -> &[OrientationExpression] {
        &self.orientations.raw
    }

    pub fn orientations_mut(&mut self) -> &mut [OrientationExpression] {
        &mut self.orientations.raw
    }

    pub fn orientation(&self, id: OrientationId) -> Option<&OrientationExpression> {
        self.orientations.get(id)
    }

    pub fn unfolded_term_count(&self) -> usize {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.leaves().count())
            .sum()
    }

    pub fn is_empty(&self) -> bool {
        self.orientations.is_empty()
    }

    pub fn retain_orientations(&mut self, mut predicate: impl FnMut(&GraphOrientation) -> bool) {
        self.orientations
            .retain(|orientation| predicate(&orientation.data.orientation));
        for (id, orientation) in self.orientations.iter_mut().enumerate() {
            orientation.id = OrientationId(id);
        }
    }

    pub fn remap_surfaces(&mut self, remap: &SurfaceIdMap) -> Result<(), CffError> {
        let mut missing = None;
        for orientation in &mut self.orientations {
            orientation.expression.map_mut(|surface| {
                if let Some(remapped) = remap.remap(*surface) {
                    *surface = remapped;
                } else {
                    missing = Some(*surface);
                }
            });
        }
        match missing {
            Some(surface) => Err(CffError::Invariant(format!(
                "surface remap does not contain {surface:?}"
            ))),
            None => Ok(()),
        }
    }

    /// Energy-surface IDs actually referenced by this expression.
    pub fn energy_surface_ids(&self) -> BTreeSet<EnergySurfaceId> {
        self.orientations
            .iter()
            .flat_map(|orientation| orientation.expression.nodes())
            .filter_map(|node| match node.data {
                SurfaceId::Energy(id) => Some(id),
                _ => None,
            })
            .collect()
    }

    pub fn normalize_raised_surfaces(&mut self, raised: &RaisedEnergySurfaceData) {
        let representatives = raised
            .groups
            .iter()
            .flat_map(|group| {
                group
                    .surface_ids
                    .first()
                    .into_iter()
                    .flat_map(move |representative| {
                        group
                            .surface_ids
                            .iter()
                            .map(move |surface| (*surface, *representative))
                    })
            })
            .collect::<BTreeMap<_, _>>();
        for orientation in &mut self.orientations {
            orientation.expression.map_mut(|surface| {
                if let SurfaceId::Energy(id) = surface
                    && let Some(representative) = representatives.get(id)
                {
                    *id = *representative;
                }
            });
        }
    }

    /// Select every non-zero residue order for a raised energy-surface group.
    pub fn select_energy_surface_residue(&self, group: &RaisedEnergySurfaceGroup) -> Vec<Self> {
        let Some(representative) = group.surface_ids.first().copied() else {
            return Vec::new();
        };
        let mut normalized = self.clone();
        normalized.normalize_raised_surfaces(&RaisedEnergySurfaceData {
            groups: vec![group.clone()].into(),
        });

        (1..=group.max_occurrence)
            .map(|occurrence| {
                let mut residue = normalized.clone();
                for orientation in &mut residue.orientations {
                    orientation.expression.retain_branches_with_value_count(
                        &SurfaceId::Energy(representative),
                        occurrence,
                    );
                    orientation.expression.map_mut(|surface| {
                        if *surface == SurfaceId::Energy(representative) {
                            *surface = SurfaceId::Unit;
                        }
                    });
                }
                residue
            })
            .collect()
    }

    pub fn max_surface_occurrence(&self, surface: SurfaceId) -> usize {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.max_value_count_on_branch(&surface))
            .max()
            .unwrap_or(0)
    }

    /// Discover equivalence classes induced by replacing raised parallel edges.
    ///
    /// `representative` maps every energy-carrying edge to the canonical edge of
    /// its raised-propagator class. Only surfaces referenced by this expression
    /// participate, even when the arena is shared with another expression.
    pub fn raised_energy_surfaces(
        &self,
        surfaces: &SurfaceCache,
        mut representative: impl FnMut(EdgeId) -> EdgeId,
    ) -> Result<RaisedEnergySurfaceData, CffError> {
        let mut normalized = Vec::new();
        for id in self.energy_surface_ids() {
            let surface = surfaces.energy_surfaces().get(id).ok_or_else(|| {
                CffError::Invariant(format!(
                    "expression references missing energy surface {id:?}"
                ))
            })?;
            let mut energies = surface
                .energies
                .iter()
                .copied()
                .map(&mut representative)
                .collect::<Vec<_>>();
            energies.sort_unstable();
            normalized.push((id, energies, surface.external_shift.clone()));
        }

        let mut groups: Vec<RaisedEnergySurfaceGroup> = Vec::new();
        for (id, energies, shift) in &normalized {
            if let Some(group) = groups.iter_mut().find(|group| {
                let representative = group.surface_ids[0];
                normalized
                    .iter()
                    .any(|(candidate, candidate_energies, candidate_shift)| {
                        *candidate == representative
                            && candidate_energies == energies
                            && candidate_shift == shift
                    })
            }) {
                group.surface_ids.push(*id);
            } else {
                groups.push(RaisedEnergySurfaceGroup {
                    surface_ids: vec![*id],
                    max_occurrence: 0,
                });
            }
        }

        let mut normalized_expression = self.clone();
        normalized_expression.normalize_raised_surfaces(&RaisedEnergySurfaceData {
            groups: groups.clone().into(),
        });
        for group in &mut groups {
            group.max_occurrence = normalized_expression
                .max_surface_occurrence(SurfaceId::Energy(group.surface_ids[0]));
        }
        Ok(RaisedEnergySurfaceData {
            groups: groups.into(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EdgeOrientation, GraphOrientation, NodeId};

    fn expression_with_raised_surfaces() -> CffExpression {
        let first = SurfaceId::Energy(EnergySurfaceId(0));
        let equivalent = SurfaceId::Energy(EnergySurfaceId(1));
        let mut tree = ExpressionTree::from_root(first);
        tree.insert(NodeId::ROOT, equivalent).unwrap();
        CffExpression::new(vec![OrientationExpression {
            id: OrientationId(0),
            data: OrientationData {
                orientation: GraphOrientation::new(),
            },
            expression: tree,
        }])
    }

    #[test]
    fn raised_surface_residue_selects_each_pole_order() {
        let expression = expression_with_raised_surfaces();
        let group = RaisedEnergySurfaceGroup {
            surface_ids: vec![EnergySurfaceId(0), EnergySurfaceId(1)],
            max_occurrence: 2,
        };

        let residues = expression.select_energy_surface_residue(&group);

        assert_eq!(residues.len(), 2);
        assert!(residues[0].orientations()[0].expression.is_empty());
        assert_eq!(
            residues[1].orientations()[0]
                .expression
                .root_to_leaf_paths(),
            vec![vec![&SurfaceId::Unit, &SurfaceId::Unit]]
        );
    }

    #[test]
    fn retaining_orientations_renumbers_the_survivors() {
        let mut expression = expression_with_raised_surfaces();
        let mut second = expression.orientations()[0].clone();
        second.id = OrientationId(1);
        second.data.orientation = GraphOrientation::from_iter([EdgeOrientation::Reversed]);
        expression.orientations.push(second);

        expression.retain_orientations(|orientation| {
            orientation.get(crate::EdgeId::new(0)) == Some(&EdgeOrientation::Reversed)
        });

        assert_eq!(expression.orientations().len(), 1);
        assert_eq!(expression.orientations()[0].id, OrientationId(0));
    }
}
