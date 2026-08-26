use bincode::{Decode, Encode};
use linnet::half_edge::involution::EdgeIndex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use symbolica::{atom::Atom, function};

use crate::symbols::S;

#[derive(
    Debug,
    Clone,
    Copy,
    Default,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Encode,
    Decode,
    JsonSchema,
)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum MediumMode {
    #[default]
    Vacuum,
    ThermodynamicEquilibrium,
    ZeroTemperatureEquilibrium,
}

impl MediumMode {
    pub fn is_finite_temperature(&self) -> bool {
        matches!(self, Self::ThermodynamicEquilibrium)
    }
}

#[derive(
    Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Encode, Decode,
)]
pub struct ThermalDistributionFactor {
    pub edge_id: EdgeIndex,
    pub sign: i32,
    pub derivative_order: usize,
}

impl ThermalDistributionFactor {
    pub fn to_atom(self, is_finite_temperature: bool) -> Atom {
        function!(
            S.thermal_distribution,
            self.edge_id.0 as i64,
            self.derivative_order as i64,
            i64::from(is_finite_temperature),
            self.sign
        )
    }
}

#[derive(
    Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Encode, Decode,
)]
pub struct ThermalNumerator {
    pub positive_energies: Vec<EdgeIndex>,
    pub negative_energies: Vec<EdgeIndex>,
}

impl ThermalNumerator {
    pub(crate) fn from_edge_lists_canonicalized(
        mut positive_energies: Vec<EdgeIndex>,
        mut negative_energies: Vec<EdgeIndex>,
    ) -> (Self, i32) {
        positive_energies.sort();
        negative_energies.sort();
        let swap = match negative_energies.len().cmp(&positive_energies.len()) {
            std::cmp::Ordering::Greater => true,
            std::cmp::Ordering::Less => false,
            std::cmp::Ordering::Equal => positive_energies > negative_energies,
        };
        if swap {
            std::mem::swap(&mut positive_energies, &mut negative_energies);
        }
        (
            Self {
                positive_energies,
                negative_energies,
            },
            if swap { -1 } else { 1 },
        )
    }

    pub(crate) fn is_trivial(&self) -> bool {
        self.positive_energies.len() == 1 && self.negative_energies.is_empty()
    }

    pub fn to_atom(&self, is_finite_temperature: bool) -> Atom {
        let product = |sign: i32| {
            self.positive_energies
                .iter()
                .map(|edge| (*edge, sign))
                .chain(self.negative_energies.iter().map(|edge| (*edge, -sign)))
                .fold(Atom::num(1), |acc, (edge_id, sign)| {
                    acc * ThermalDistributionFactor {
                        edge_id,
                        sign,
                        derivative_order: 0,
                    }
                    .to_atom(is_finite_temperature)
                })
        };
        product(1) - product(-1)
    }
}

/// Products produced by thermal edge contractions stay separate from rational
/// CFF coefficients, so variant fusion cannot discard distribution functions.
#[derive(
    Debug,
    Clone,
    Default,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Encode,
    Decode,
)]
pub struct ThermalWeight {
    pub medium_mode: MediumMode,
    pub numerators: Vec<ThermalNumerator>,
    pub distributions: Vec<ThermalDistributionFactor>,
}

impl ThermalWeight {
    pub fn to_atom(&self) -> Atom {
        let finite = self.medium_mode.is_finite_temperature();
        self.numerators
            .iter()
            .map(|n| n.to_atom(finite))
            .chain(
                self.distributions
                    .iter()
                    .map(|factor| factor.to_atom(finite)),
            )
            .fold(Atom::num(1), |acc, factor| acc * factor)
    }

    pub(crate) fn canonicalize(&mut self) {
        self.numerators.sort();
        self.distributions.sort();
    }

    pub fn remap_internal_edges(&mut self, edge_map: &std::collections::BTreeMap<usize, usize>) {
        let remap = |edge: &mut EdgeIndex| {
            edge.0 = edge_map.get(&edge.0).copied().unwrap_or(edge.0);
        };
        for numerator in &mut self.numerators {
            for edge in numerator
                .positive_energies
                .iter_mut()
                .chain(&mut numerator.negative_energies)
            {
                remap(edge);
            }
            numerator.positive_energies.sort();
            numerator.negative_energies.sort();
        }
        for factor in &mut self.distributions {
            remap(&mut factor.edge_id);
        }
        self.canonicalize();
    }

    pub(crate) fn product(&self, rhs: &Self) -> Self {
        let mut result = self.clone();
        if result.medium_mode == MediumMode::Vacuum {
            result.medium_mode = rhs.medium_mode;
        }
        result.numerators.extend(rhs.numerators.iter().cloned());
        result
            .distributions
            .extend(rhs.distributions.iter().copied());
        result.canonicalize();
        result
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use symbolica::atom::AtomCore;

    use super::*;
    use crate::{Generate3DExpressionOptions, expression::AllOrientations, generate_3d_expression};

    #[test]
    fn thermal_numerator_keeps_sign_and_distribution_structure() {
        let (numerator, sign) = ThermalNumerator::from_edge_lists_canonicalized(
            vec![EdgeIndex(1)],
            vec![EdgeIndex(4), EdgeIndex(2)],
        );
        assert_eq!(sign, -1);
        assert_eq!(
            numerator.positive_energies,
            vec![EdgeIndex(2), EdgeIndex(4)]
        );
        assert_eq!(numerator.negative_energies, vec![EdgeIndex(1)]);
        for finite in [false, true] {
            let factor = |edge, sign| {
                ThermalDistributionFactor {
                    edge_id: EdgeIndex(edge),
                    sign,
                    derivative_order: 0,
                }
                .to_atom(finite)
            };
            assert_eq!(
                numerator.to_atom(finite),
                factor(2, 1) * factor(4, 1) * factor(1, -1)
                    - factor(2, -1) * factor(4, -1) * factor(1, 1)
            );
        }
    }

    #[test]
    fn thermal_weights_survive_fusion_remapping_and_serialization() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let mut expression = generate_3d_expression(
            &parsed,
            &Generate3DExpressionOptions {
                medium_mode: MediumMode::ThermodynamicEquilibrium,
                ..Default::default()
            },
        )
        .unwrap()
        .expression;
        assert_eq!(expression.orientations.len().0, 16);
        let before = expression.to_atom(AllOrientations);
        expression = expression.fuse_compatible_variants();
        assert_eq!(before, expression.to_atom(AllOrientations));
        let map = crate::EnergyEdgeIndexMap {
            internal: (0..4).map(|edge| (edge, edge + 10)).collect(),
            external: BTreeMap::new(),
            orientation_edge_count: 14,
        };
        expression = expression.remap_energy_edge_indices(&map);
        for orientation in &expression.orientations {
            for variant in &orientation.variants {
                for numerator in &variant.thermal_weight.numerators {
                    assert!(
                        numerator
                            .positive_energies
                            .iter()
                            .chain(&numerator.negative_energies)
                            .all(|edge| edge.0 >= 10)
                    );
                }
                assert!(
                    variant
                        .thermal_weight
                        .distributions
                        .iter()
                        .all(|factor| factor.edge_id.0 >= 10)
                );
            }
        }
        #[cfg(feature = "eval")]
        {
            let json = serde_json::to_string(&expression).unwrap();
            let decoded: crate::ThreeDExpression<crate::OrientationID> =
                serde_json::from_str(&json).unwrap();
            assert_eq!(
                expression.to_atom(AllOrientations),
                decoded.to_atom(AllOrientations)
            );
        }
    }

    #[test]
    fn thermal_box_vacuum_limit_agrees_with_vacuum_cff() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let vacuum = generate_3d_expression(&parsed, &Generate3DExpressionOptions::default())
            .unwrap()
            .expression;
        let vacuum_atom = vacuum
            .surfaces
            .substitute_energies(&vacuum.to_atom(AllOrientations), &[]);
        for medium_mode in [
            MediumMode::ThermodynamicEquilibrium,
            MediumMode::ZeroTemperatureEquilibrium,
        ] {
            let thermal = generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    medium_mode,
                    ..Default::default()
                },
            )
            .unwrap()
            .expression;
            let mut atom = thermal
                .surfaces
                .substitute_energies(&thermal.to_atom(AllOrientations), &[]);
            for edge in 0..4 {
                for sign in [-1, 1] {
                    for derivative_order in 0..4 {
                        atom = atom
                            .replace(
                                ThermalDistributionFactor {
                                    edge_id: EdgeIndex(edge),
                                    sign,
                                    derivative_order,
                                }
                                .to_atom(medium_mode.is_finite_temperature()),
                            )
                            .with(Atom::num(i64::from(sign == 1 && derivative_order == 0)));
                    }
                }
            }
            assert_eq!((atom - &vacuum_atom).expand(), Atom::Zero);
        }
    }

    #[test]
    fn thermal_tadpole_retains_both_orientation_weights() {
        let mut parsed = crate::graph_io::test_graphs::box_graph();
        parsed.internal_edges.truncate(1);
        parsed.internal_edges[0].head = 0;
        parsed.external_edges.clear();
        parsed.external_names.clear();
        parsed.internal_edges[0]
            .signature
            .external_signature
            .clear();
        parsed.node_name_to_internal = BTreeMap::from([("v0".to_string(), 0)]);
        for medium_mode in [
            MediumMode::ThermodynamicEquilibrium,
            MediumMode::ZeroTemperatureEquilibrium,
        ] {
            let expression = generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    medium_mode,
                    ..Default::default()
                },
            )
            .unwrap()
            .expression;
            assert_eq!(expression.orientations.len().0, 2);
            let factors = expression
                .orientations
                .iter()
                .flat_map(|orientation| &orientation.variants)
                .flat_map(|variant| &variant.thermal_weight.distributions)
                .copied()
                .collect::<Vec<_>>();
            assert_eq!(
                factors,
                vec![
                    ThermalDistributionFactor {
                        edge_id: EdgeIndex(0),
                        sign: 1,
                        derivative_order: 0
                    },
                    ThermalDistributionFactor {
                        edge_id: EdgeIndex(0),
                        sign: -1,
                        derivative_order: 0
                    },
                ]
            );
        }
    }
}
