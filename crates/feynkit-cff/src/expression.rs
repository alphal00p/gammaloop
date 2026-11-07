use std::collections::{BTreeMap, BTreeSet};

use serde::ser::SerializeStruct;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use symbolica::atom::Atom;
use typed_index_collections::TiVec;

use crate::{
    CffError, EdgeId, EnergySurfaceId, ExpressionTree, GraphOrientation, OrientationId,
    OrientationWeight, RaisedEnergySurfaceData, RaisedEnergySurfaceGroup, SurfaceCache, SurfaceId,
    SurfaceIdMap,
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
///
/// Serde writes `weight` in every compact format, including the
/// `bincode(with_serde)` path used by [`CffExpression`], because such formats
/// decode every field by position. Human-readable formats omit a unit weight,
/// and a missing weight reads as the unit weight.
#[derive(Clone, Debug, PartialEq, Eq, bincode::Encode, bincode::Decode)]
pub struct OrientationExpression {
    pub id: OrientationId,
    pub data: OrientationData,
    /// Exact factor multiplying `expression`: unit unless the orientation
    /// directs self-loops.
    pub weight: OrientationWeight,
    pub expression: ExpressionTree<SurfaceId>,
}

/// Field order shared by the serializer and every positional decoder.
#[derive(Deserialize)]
#[serde(rename = "OrientationExpression")]
struct OrientationExpressionFields {
    id: OrientationId,
    data: OrientationData,
    #[serde(default)]
    weight: OrientationWeight,
    expression: ExpressionTree<SurfaceId>,
}

impl Serialize for OrientationExpression {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let omit_weight = serializer.is_human_readable() && self.weight.is_unit();
        let mut fields = serializer
            .serialize_struct("OrientationExpression", if omit_weight { 3 } else { 4 })?;
        fields.serialize_field("id", &self.id)?;
        fields.serialize_field("data", &self.data)?;
        if omit_weight {
            fields.skip_field("weight")?;
        } else {
            fields.serialize_field("weight", &self.weight)?;
        }
        fields.serialize_field("expression", &self.expression)?;
        fields.end()
    }
}

impl<'de> Deserialize<'de> for OrientationExpression {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let fields = OrientationExpressionFields::deserialize(deserializer)?;
        Ok(Self {
            id: fields.id,
            data: fields.data,
            weight: fields.weight,
            expression: fields.expression,
        })
    }
}

impl OrientationExpression {
    /// Unfold the tree into denominator products, one vector per additive term.
    ///
    /// The products exclude [`Self::weight`], which multiplies each of them.
    pub fn denominator_products(&self) -> Vec<Vec<SurfaceId>> {
        self.expression
            .root_to_leaf_paths()
            .into_iter()
            .map(|path| path.into_iter().copied().collect())
            .collect()
    }

    /// Lower this weighted tree to a sum of products of inverse surfaces.
    pub fn to_atom(&self) -> Atom {
        self.weight.to_atom() * self.expression.to_atom_inverse()
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
    use bincode::config;
    use bincode::serde as bincode_serde;

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
            weight: OrientationWeight::UNIT,
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

    /// Unit-weight, weighted and emptied expressions, each followed by a
    /// second weighted orientation so a misread field shifts later bytes.
    fn binary_controls() -> Vec<CffExpression> {
        let unit = expression_with_raised_surfaces();
        let mut weighted = unit.clone();
        weighted.orientations_mut()[0].weight = OrientationWeight::from_halvings(1);
        let mut emptied = weighted.clone();
        emptied.orientations_mut()[0]
            .expression
            .retain_branches_with_value_count(&SurfaceId::Unit, 1);
        assert!(emptied.orientations()[0].expression.is_empty());
        [unit, weighted, emptied]
            .into_iter()
            .map(|mut expression| {
                let mut trailing = expression.orientations()[0].clone();
                trailing.id = OrientationId(1);
                trailing.weight = OrientationWeight::from_halvings(2);
                expression.orientations.push(trailing);
                expression
            })
            .collect()
    }

    #[test]
    fn orientations_round_trip_through_derived_and_serde_bincode() {
        for expression in binary_controls() {
            for orientation in expression.orientations() {
                let bytes = bincode::encode_to_vec(orientation, config::standard()).unwrap();
                let (decoded, read): (OrientationExpression, usize) =
                    bincode::decode_from_slice(&bytes, config::standard()).unwrap();
                assert_eq!((&decoded, read), (orientation, bytes.len()));

                let bytes = bincode_serde::encode_to_vec(orientation, config::standard()).unwrap();
                let (decoded, read): (OrientationExpression, usize) =
                    bincode_serde::decode_from_slice(&bytes, config::standard()).unwrap();
                assert_eq!((&decoded, read), (orientation, bytes.len()));
            }
        }
    }

    #[test]
    fn expressions_round_trip_through_bincode() {
        for expression in binary_controls() {
            let bytes = bincode::encode_to_vec(&expression, config::standard()).unwrap();
            let (decoded, read): (CffExpression, usize) =
                bincode::decode_from_slice(&bytes, config::standard()).unwrap();
            assert_eq!((decoded, read), (expression, bytes.len()));
        }
    }

    #[test]
    fn json_omits_only_unit_weights() {
        for expression in binary_controls() {
            let json = serde_json::to_value(&expression).unwrap();
            for (orientation, encoded) in expression
                .orientations()
                .iter()
                .zip(json["orientations"].as_array().unwrap())
            {
                assert_eq!(
                    encoded.get("weight").is_none(),
                    orientation.weight.is_unit()
                );
            }
            assert_eq!(
                serde_json::from_value::<CffExpression>(json).unwrap(),
                expression
            );
        }
    }

    #[test]
    fn serialized_orientations_without_a_weight_keep_unit_weight() {
        let mut expression = expression_with_raised_surfaces();
        let legacy = serde_json::to_value(&expression).unwrap();
        assert!(legacy["orientations"][0].get("weight").is_none());
        assert_eq!(
            serde_json::from_value::<CffExpression>(legacy).unwrap(),
            expression
        );

        expression.orientations_mut()[0].weight = OrientationWeight::from_halvings(1);
        let weighted = serde_json::to_value(&expression).unwrap();
        assert_eq!(
            serde_json::from_value::<CffExpression>(weighted).unwrap(),
            expression
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
