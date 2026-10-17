use crate::cff::VertexSet;
use crate::graph::LoopMomentumBasis;
use crate::momentum::sample::{ExternalFourMomenta, ExternalThreeMomenta, LoopMomenta};
use crate::utils::{
    F, FloatLike, cut_energy, external_energy_atom_from_index,
    newton_solver::{NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity},
    ose_atom_from_index,
};
use bincode_trait_derive::{Decode, Encode};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use rug::Rational;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use symbolica::atom::Atom;
use tracing::warn;
use typed_index_collections::TiVec;

use super::esurface::Esurface;
use super::esurface::ExternalShift;
use super::radial::EnergySurfaceRay;
pub use super::surface::HsurfaceID;

pub type HsurfaceCollection = TiVec<HsurfaceID, Hsurface>;
pub type HsurfaceCache<T> = TiVec<HsurfaceID, T>;

#[derive(Serialize, Deserialize, Debug, Clone, Encode, Decode)]
/// H-surface of the supergraph, is most likely E-surface of the amplitude, kind of badly named.
pub struct Hsurface {
    pub positive_energies: Vec<EdgeIndex>,
    pub negative_energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
}

impl PartialEq for Hsurface {
    fn eq(&self, other: &Self) -> bool {
        self.positive_energies == other.positive_energies
            && self.negative_energies == other.negative_energies
            && self.external_shift == other.external_shift
    }
}

impl Eq for Hsurface {}

impl Hsurface {
    /// Prepare the retained signed equation and enumerate its simple positive roots.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn positive_radius_roots<T: FloatLike>(
        &self,
        loop_momenta: &LoopMomenta<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> eyre::Result<(EnergySurfaceRay<T>, Vec<NewtonIterationResult<T>>)> {
        let center = LoopMomenta::from_iter(
            loop_momenta
                .iter()
                .map(|momentum| momentum.map_ref(&|value| value.zero())),
        );
        let ray = EnergySurfaceRay::routed(
            self.positive_energies
                .iter()
                .map(|&edge| (edge, 1))
                .chain(self.negative_energies.iter().map(|&edge| (edge, -1))),
            &self.external_shift,
            loop_momenta,
            &center,
            external_momenta,
            masses,
            lmb,
        );
        let roots = ray.positive_radius_roots(e_cm, diagnostics, identity)?;
        Ok((ray, roots))
    }

    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        let (symbolic_positive_energies, symbolic_negative_energies) =
            [&self.positive_energies, &self.negative_energies]
                .iter()
                .map(|energies| {
                    energies
                        .iter()
                        .map(|i| {
                            if cut_edges.contains(i) {
                                cut_energy(*i)
                            } else {
                                ose_atom_from_index(*i)
                            }
                        })
                        .collect_vec()
                })
                .collect_tuple()
                .unwrap_or_else(|| unreachable!());

        let symbolic_shift = self
            .external_shift
            .iter()
            .fold(Atom::new(), |sum, (i, sign)| {
                Atom::num(*sign) * external_energy_atom_from_index(*i) + &sum
            });

        let symbolic_sum_positive_energies = symbolic_positive_energies
            .iter()
            .fold(Atom::new(), |sum, e| sum + e);

        let symbolic_sum_negative_energies = symbolic_negative_energies
            .iter()
            .fold(Atom::new(), |sum, e| sum + e);

        symbolic_sum_positive_energies - &symbolic_sum_negative_energies + &symbolic_shift
    }

    /// Whether the surface vanishes for every loop momentum, so that a denominator
    /// built from it divides by zero at every sample. On-shell energies of different
    /// spatial routings or squared masses are independent functions of the loop momenta,
    /// so this happens exactly when the positive and negative energies pair up into
    /// edges with both mass and spatial momentum equal up to sign, and the external
    /// energy shift is zero. All comparisons are exact: a near coincidence stays
    /// finite in some precision and is left to the stability checks.
    pub(crate) fn vanishes_identically(
        &self,
        lmb: &LoopMomentumBasis,
        masses: &EdgeVec<Option<Complex<F<f64>>>>,
        externals: &ExternalFourMomenta<F<f64>>,
    ) -> eyre::Result<bool> {
        let spatial_externals: ExternalThreeMomenta<F<f64>> =
            externals.iter().map(|momentum| momentum.spatial).collect();
        let mut unpaired = self.negative_energies.clone();
        for &positive in &self.positive_energies {
            let mut partner = None;
            for (position, &negative) in unpaired.iter().enumerate() {
                // An edge without a mass value never pairs.
                if masses[positive]
                    .is_some_and(|m| masses[negative].is_some_and(|n| m == n || m == -n))
                    && lmb.edge_signatures[positive].spatial_equality_up_to_sign(
                        &lmb.edge_signatures[negative],
                        &spatial_externals,
                    )?
                {
                    partner = Some(position);
                    break;
                }
            }
            let Some(position) = partner else {
                return Ok(false);
            };
            unpaired.swap_remove(position);
        }
        if !unpaired.is_empty() {
            return Ok(false);
        }
        let mut shift = Rational::new();
        for (edge, sign) in &self.external_shift {
            let signature = lmb.edge_signatures[*edge].external.to_momtrop_format();
            for (coefficient, momentum) in signature.into_iter().zip(externals) {
                if coefficient != 0 {
                    let energy =
                        Rational::from_f64(momentum.temporal.value.0).ok_or_else(|| {
                            eyre::eyre!("nonfinite external energy {}", momentum.temporal.value)
                        })?;
                    shift += energy * (sign * coefficient as i64);
                }
            }
        }
        Ok(shift == 0)
    }

    #[allow(dead_code)]
    pub(crate) fn equality_under_energy_conservation(
        &self,
        other: &Esurface,
        constraints: &[&Esurface],
    ) -> Option<bool> {
        if !self.external_shift.is_empty() {
            warn!("this is not handled yet");
            return None;
        }
        constraints
            .iter()
            .find(|esurface| {
                self.negative_energies
                    .iter()
                    .all(|index| esurface.energies.contains(index))
            })
            .map(|constraint| {
                let energies_to_be_added = constraint
                    .energies
                    .iter()
                    .filter(|index| !self.negative_energies.contains(index));

                let mut new_positive_energies = self.positive_energies.clone();
                new_positive_energies.extend(energies_to_be_added);
                new_positive_energies.sort();

                let new_esurface = Esurface {
                    energies: new_positive_energies,
                    external_shift: constraint.external_shift.clone(),
                    vertex_set: VertexSet::dummy(),
                };

                other == &new_esurface
            })
    }

    pub fn equality_by_try_convert(&self, other: &Esurface) -> bool {
        let negative_as_external_shift = self
            .negative_energies
            .iter()
            .map(|e| (*e, -1))
            .collect_vec();

        let self_as_esurface = Esurface {
            energies: self.positive_energies.clone(),
            external_shift: negative_as_external_shift,
            vertex_set: VertexSet::dummy(),
        };

        self_as_esurface == *other
    }
}

#[cfg(test)]
mod tests {

    use linnet::half_edge::involution::EdgeIndex;
    use symbolica::atom::{Atom, AtomCore};

    use symbolica::parse;

    use crate::cff::{VertexSet, esurface::Esurface};

    use super::Hsurface;

    type RadialTestEnergy = (i64, [f64; 3], [f64; 3], f64);

    fn radial_roots(
        energies: &[RadialTestEnergy],
        incoming_energy: f64,
    ) -> eyre::Result<Vec<crate::utils::newton_solver::NewtonIterationResult<f64>>> {
        use crate::{
            graph::LoopMomentumBasis,
            momentum::{FourMomentum, ThreeMomentum, signature::LoopExtSignature},
            utils::{
                F,
                newton_solver::{RadialRootDiagnostics, RadialRootIdentity},
                test_utils::dummy_hedge_graph,
            },
        };
        use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};

        let count = energies.len();
        let graph = dummy_hedge_graph(2 * count + 1);
        let signatures = (0..2 * count + 1).map(|edge| {
            let mut internal = vec![0_isize; count];
            let mut external = vec![0_isize; count + 1];
            if edge < count {
                internal[edge] = 1;
                external[edge] = 1;
            } else {
                external[edge - count] = 1;
            }
            LoopExtSignature::from((internal, external))
        });
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: (0..count).map(EdgeIndex).collect(),
            ext_edges: (count..2 * count + 1).map(EdgeIndex).collect(),
            edge_signatures: graph.new_edgevec_from_iter(signatures).unwrap(),
        };
        let loop_momenta = energies
            .iter()
            .map(|(_, v, _, _)| ThreeMomentum::new(F(v[0]), F(v[1]), F(v[2])))
            .collect();
        let external_momenta = energies
            .iter()
            .map(|(_, _, b, _)| FourMomentum::from_args(F(0.0), F(b[0]), F(b[1]), F(b[2])))
            .chain(std::iter::once(FourMomentum::from_args(
                F(incoming_energy),
                F(0.0),
                F(0.0),
                F(0.0),
            )))
            .collect();
        let masses = graph
            .new_edgevec_from_iter(
                (0..2 * count + 1).map(|edge| F(energies.get(edge).map_or(0.0, |energy| energy.3))),
            )
            .unwrap();
        let surface = Hsurface {
            positive_energies: energies
                .iter()
                .enumerate()
                .flat_map(|(edge, (sign, _, _, _))| {
                    std::iter::repeat_n(EdgeIndex(edge), (*sign).max(0) as usize)
                })
                .collect(),
            negative_energies: energies
                .iter()
                .enumerate()
                .flat_map(|(edge, (sign, _, _, _))| {
                    std::iter::repeat_n(EdgeIndex(edge), (-sign).max(0) as usize)
                })
                .collect(),
            external_shift: vec![(EdgeIndex(2 * count), -1)],
            vertex_set: VertexSet::dummy(),
        };
        surface
            .positive_radius_roots(
                &loop_momenta,
                &external_momenta,
                &masses,
                &lmb,
                &F(2.0),
                &mut RadialRootDiagnostics::default(),
                &RadialRootIdentity::new("signed radial test".to_string()),
            )
            .map(|(_, roots)| roots)
    }

    #[test]
    fn signed_radial_sunset_roots_at_rest_and_with_spatial_external_momentum() {
        let origin = [0.0; 3];
        for sign in [-1, 1] {
            for pz in [0.0, 0.5] {
                let energies = [
                    (1, [1.0, 0.0, 0.0], origin, 0.0),
                    (1, [0.0, 1.0, 0.0], origin, 0.0),
                    (sign, [-1.0, -1.0, 0.0], [0.0, 0.0, pz], 0.0),
                ];
                let roots = radial_roots(&energies, 2.0).unwrap();
                assert_eq!(roots.len(), 1, "sign={sign}, pz={pz}");
                let expected = if pz == 0.0 {
                    2.0 / (2.0 + sign as f64 * 2.0_f64.sqrt())
                } else {
                    (8.0 - sign as f64 * 34.0_f64.sqrt()) / 4.0
                };
                assert!((roots[0].solution.0 - expected).abs() < 1.0e-12);
                assert!(roots[0].derivative_at_solution.0 > 0.0);
            }
        }
    }

    #[test]
    fn signed_radial_enumerates_both_crossings_and_retains_signed_jacobians() {
        let energies = [
            (2, [1.0, 0.0, 0.0], [0.0; 3], 1.0),
            (-1, [1.9, 0.0, 0.0], [0.0; 3], 0.0),
        ];
        let roots = radial_roots(&energies, 1.0).unwrap();
        assert_eq!(roots.len(), 2);
        let discriminant = 3.8_f64.powi(2) - 4.0 * 0.39 * 3.0;
        for (index, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let expected = (3.8 + sign * discriminant.sqrt()) / 0.78;
            assert!((roots[index].solution.0 - expected).abs() < 1.0e-11);
            assert!(roots[index].derivative_at_solution.0 * sign > 0.0);
        }
        let negated =
            energies.map(|(sign, velocity, offset, mass)| (-sign, velocity, offset, mass));
        let reversed = radial_roots(&negated, -1.0).unwrap();
        for (root, reverse) in roots.iter().zip(&reversed) {
            assert!((root.solution.0 - reverse.solution.0).abs() < 1.0e-12);
            assert!(
                (root.derivative_at_solution.0 + reverse.derivative_at_solution.0).abs() < 1.0e-12
            );
        }
    }

    #[test]
    fn signed_radial_proves_closed_channels_and_rejects_the_squared_opposite_sheet() {
        let origin = [0.0; 3];
        // A negative discriminant proves the complete absence of a real root.
        assert!(
            radial_roots(
                &[
                    (2, [1.0, 0.0, 0.0], origin, 1.0),
                    (-1, [1.9, 0.0, 0.0], origin, 0.0),
                ],
                0.0
            )
            .unwrap()
            .is_empty()
        );
        // sqrt(r²+1)+2 has a positive root after squaring, on the wrong sheet.
        assert!(
            radial_roots(&[(1, [1.0, 0.0, 0.0], origin, 1.0)], -2.0)
                .unwrap()
                .is_empty()
        );
        // An affine equation with only a negative root is also a closed channel.
        assert!(
            radial_roots(&[(-1, [1.0, 0.0, 0.0], origin, 0.0)], 1.0)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn signed_radial_distinguishes_excluded_endpoints_from_nonisolated_cuts() {
        let origin = [0.0; 3];
        assert!(
            radial_roots(&[(1, [1.0, 0.0, 0.0], origin, 0.0)], 0.0)
                .unwrap()
                .is_empty()
        );
        // sqrt(r²+1)-1 only touches zero at the excluded endpoint.
        assert!(
            radial_roots(&[(1, [1.0, 0.0, 0.0], origin, 1.0)], 1.0)
                .unwrap()
                .is_empty()
        );
        let tangent =
            radial_roots(&[(1, [1.0, 0.0, 0.0], [-2.0, 0.0, 0.0], 1.0)], 1.0).unwrap_err();
        assert!(tangent.to_string().contains("tangent"));
        let cusp = radial_roots(&[(1, [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.0)], 0.0).unwrap_err();
        assert!(cusp.to_string().contains("tangent"));
        let identically_zero = radial_roots(
            &[
                (1, [1.0, 0.0, 0.0], origin, 0.0),
                (-1, [1.0, 0.0, 0.0], origin, 0.0),
            ],
            0.0,
        )
        .unwrap_err();
        assert!(identically_zero.to_string().contains("identically zero"));
    }

    #[test]
    fn signed_radial_unsupported_and_nonfinite_rays_are_errors_not_closed_channels() {
        let origin = [0.0; 3];
        let unsupported = radial_roots(
            &[
                (1, [1.0, 0.0, 0.0], origin, 1.0),
                (-1, [2.0, 0.0, 0.0], origin, 1.0),
            ],
            1.0,
        )
        .unwrap_err();
        assert!(unsupported.to_string().contains("unsupported ray"));
        assert!(radial_roots(&[(1, [f64::NAN, 0.0, 0.0], origin, 0.0)], 1.0).is_err());
    }

    #[test]
    fn test_equality_under_energy_conservation() {
        let constraint = Esurface {
            energies: vec![EdgeIndex::from(1), EdgeIndex::from(2)],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        };

        let hsurface = Hsurface {
            positive_energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            negative_energies: vec![EdgeIndex::from(2)],
            external_shift: vec![],
            vertex_set: VertexSet::dummy(),
        };

        let other = Esurface {
            energies: vec![EdgeIndex::from(1), EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        };

        assert!(
            hsurface
                .equality_under_energy_conservation(&other, &[&constraint])
                .unwrap()
        );
    }

    #[test]
    fn surfaces_vanish_identically_only_for_paired_energies_without_energy_shift() {
        use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
        use spenso::algebra::complex::Complex;

        use crate::{
            graph::LoopMomentumBasis,
            momentum::{
                FourMomentum,
                sample::{ExternalFourMomenta, ExternalIndex},
                signature::LoopExtSignature,
            },
            utils::{F, test_utils::dummy_hedge_graph},
        };

        // Edges 0 and 4 carry k, edge 1 carries k + p, and edges 2 and 3 are the
        // incoming and outgoing external p.
        let graph = dummy_hedge_graph(5);
        let edge_signatures = graph
            .new_edgevec_from_iter(
                [
                    (vec![1], vec![0, 0]),
                    (vec![1], vec![1, 0]),
                    (vec![0], vec![1, 0]),
                    (vec![0], vec![0, 1]),
                    (vec![1], vec![0, 0]),
                ]
                .map(LoopExtSignature::from),
            )
            .unwrap();
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![EdgeIndex(0)].into(),
            ext_edges: vec![EdgeIndex(2), EdgeIndex(3)].into(),
            edge_signatures,
        };
        let surface = |positive: &[usize], negative: &[usize], external_shift| Hsurface {
            positive_energies: positive.iter().copied().map(EdgeIndex).collect(),
            negative_energies: negative.iter().copied().map(EdgeIndex).collect(),
            external_shift,
            vertex_set: VertexSet::dummy(),
        };
        let masses = |values: [Option<f64>; 5]| {
            graph
                .new_edgevec_from_iter(values.map(|mass| mass.map(|m| Complex::new_re(F(m)))))
                .unwrap()
        };
        let externals = |energy: f64, py: f64| -> ExternalFourMomenta<F<f64>> {
            let p = FourMomentum::from_args(F(energy), F(0.0), F(py), F(0.0));
            [p, p].into_iter().collect()
        };
        let equal = masses([Some(1.0), Some(1.0), None, None, Some(1.0)]);
        let rest = externals(0.0, 0.0);

        // A divided difference over distinct mass parameters with equal values.
        let cycle = surface(&[4], &[0], vec![]);
        assert!(cycle.vanishes_identically(&lmb, &equal, &rest).unwrap());
        let split = masses([Some(1.0), Some(1.0), None, None, Some(2.0)]);
        assert!(!cycle.vanishes_identically(&lmb, &split, &rest).unwrap());
        for (mass, vanishes) in [(-1.0, true), ((-1.0_f64).next_up(), false), (-2.0, false)] {
            let signed = masses([Some(1.0), Some(1.0), None, None, Some(mass)]);
            assert_eq!(
                cycle.vanishes_identically(&lmb, &signed, &rest).unwrap(),
                vanishes,
                "mass {mass} paired with 1"
            );
        }
        // Reversing the spatial routing preserves the on-shell energy as well.
        let mut reversed_lmb = lmb.clone();
        reversed_lmb.edge_signatures[EdgeIndex(4)] = LoopExtSignature::from((vec![-1], vec![0, 0]));
        let opposite = masses([Some(1.0), Some(1.0), None, None, Some(-1.0)]);
        assert!(
            cycle
                .vanishes_identically(&reversed_lmb, &opposite, &rest)
                .unwrap()
        );
        // Complex masses agree up to sign, but conjugates need not have equal squares.
        for (mass, vanishes) in [
            (Complex::new(F(1.0), F(2.0)), true),
            (Complex::new(F(-1.0), F(-2.0)), true),
            (Complex::new(F(1.0), F(-2.0)), false),
        ] {
            let mut complex = equal.clone();
            complex[EdgeIndex(0)] = Some(Complex::new(F(1.0), F(2.0)));
            complex[EdgeIndex(4)] = Some(mass);
            assert_eq!(
                cycle.vanishes_identically(&lmb, &complex, &rest).unwrap(),
                vanishes,
                "mass {mass} paired with 1 + 2i"
            );
        }
        let unknown = masses([Some(1.0), Some(1.0), None, None, None]);
        assert!(!cycle.vanishes_identically(&lmb, &unknown, &rest).unwrap());
        let unpaired = surface(&[0, 1], &[4], vec![]);
        assert!(!unpaired.vanishes_identically(&lmb, &equal, &rest).unwrap());

        // The bubble E(k + p) - E(k) - p0 vanishes only at p = 0. A spatial momentum
        // splits the energies, and an energy alone leaves the constant -p0.
        let bubble = surface(&[1], &[0], vec![(EdgeIndex(2), -1)]);
        assert!(bubble.vanishes_identically(&lmb, &equal, &rest).unwrap());
        for moving in [externals(0.0, 0.5), externals(0.5, 0.0)] {
            assert!(!bubble.vanishes_identically(&lmb, &equal, &moving).unwrap());
        }
        let heavy = masses([Some(1.0), Some(2.0), None, None, Some(1.0)]);
        assert!(!bubble.vanishes_identically(&lmb, &heavy, &rest).unwrap());
        let mut nonfinite = rest.clone();
        nonfinite[ExternalIndex(0)].temporal.value = F(f64::NAN);
        assert!(
            bubble
                .vanishes_identically(&lmb, &equal, &nonfinite)
                .is_err()
        );
    }

    #[test]
    fn cache_equality_includes_the_external_shift_but_not_the_shore() {
        let base = Hsurface {
            positive_energies: vec![EdgeIndex::from(0)],
            negative_energies: vec![EdgeIndex::from(1)],
            external_shift: vec![(EdgeIndex::from(2), 1)],
            vertex_set: VertexSet::from_usize(0),
        };
        let different_shore = Hsurface {
            vertex_set: VertexSet::from_usize(1),
            ..base.clone()
        };
        let different_shift = Hsurface {
            external_shift: vec![(EdgeIndex::from(3), 1)],
            ..base.clone()
        };

        assert_eq!(base, different_shore);
        assert_ne!(base, different_shift);
        assert_ne!(
            base.to_atom(&[]).to_canonical_string(),
            different_shift.to_atom(&[]).to_canonical_string()
        );
    }

    #[test]
    fn test_to_atom() {
        let external_shift = vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)];
        let h_surface = Hsurface {
            positive_energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
            negative_energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift,
            vertex_set: VertexSet::dummy(),
        };

        let h_surface_atom = h_surface.to_atom(&[]);
        let expected_atom = parse!(
            "OSE(0) + OSE(1) - OSE(2) - OSE(3) - Q(4, spenso::cind(0)) + Q(5, spenso::cind(0))"
        );
        let diff = h_surface_atom - &expected_atom;
        let diff = diff.expand();
        assert_eq!(diff, Atom::new());
    }
}
