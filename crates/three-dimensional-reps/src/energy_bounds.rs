use itertools::Itertools;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{MomentumSignature, utils::rank_i64};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnergyDivergenceReport {
    pub edge_degree_bounds: Vec<usize>,
    /// Coordinate-axis scaling certificates, omitting loop energies on which
    /// no denominator depends.
    pub loops: Vec<EnergyDirectionReport>,
    /// Certificates for the scaling subspaces of the denominator arrangement.
    pub directions: Vec<EnergyDirectionReport>,
    /// Whether every coordinate-axis scaling limit is certified convergent.
    pub coordinate_convergent: bool,
    /// Whether every arrangement scaling limit is certified convergent.
    pub directional_convergent: bool,
    /// Whether all scaling limits considered by this report are certified.
    pub convergent: bool,
}

/// Conservative power counting on one large-energy scaling subspace.
///
/// `convergent == false` means that bound-only power counting is
/// inconclusive; it does not establish divergence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnergyDirectionReport {
    pub loop_id: Option<usize>,
    /// Denominator edge IDs of occurrences which vanish on the
    /// scaling subspace. Repeated IDs retain denominator multiplicity.
    pub zero_edges: Vec<usize>,
    /// Denominator edge IDs of occurrences which scale on the
    /// subspace. Repeated IDs retain denominator multiplicity.
    pub active_edges: Vec<usize>,
    /// Dimension of an arrangement scaling subspace; coordinate reports use
    /// `None` and are one-dimensional.
    pub nullity: Option<usize>,
    pub numerator_degree_bound: usize,
    pub denominator_degree: usize,
    /// Numerator minus denominator degree, before the radial measure of the
    /// scaling subspace is included.
    pub divergence_degree: isize,
    /// A one-way convergence certificate; `false` means inconclusive.
    pub convergent: bool,
}

#[derive(Debug, Error)]
pub enum EnergyBoundsError {
    #[error("energy-degree bound edge index {0} is outside 0..{1}")]
    EdgeOutOfRange(usize, usize),
    #[error("energy-denominator edge index {0} is outside 0..{1}")]
    DenominatorEdgeOutOfRange(usize, usize),
    #[error(
        "positive energy-degree bound for edge {0} depends on a loop-energy direction outside the active denominator span"
    )]
    NumeratorOutsideDenominatorSpan(usize),
    #[error("auto numerator generation requires at least one external momentum symbol")]
    AutoNumeratorNeedsExternal,
    #[error("energy degree arithmetic exceeds the supported integer range")]
    DegreeOverflow,
    #[error(
        "cannot enumerate an energy-denominator arrangement with {n_edges} edges using a {mask_bits}-bit subset mask"
    )]
    ArrangementTooLarge { n_edges: usize, mask_bits: u32 },
}

pub type Result<T> = std::result::Result<T, EnergyBoundsError>;

pub fn normalize_energy_degree_bounds(
    bounds: &[(usize, usize)],
    n_internal: usize,
) -> Result<Vec<usize>> {
    let mut out = vec![0; n_internal];
    for (edge, degree) in bounds {
        if *edge >= n_internal {
            return Err(EnergyBoundsError::EdgeOutOfRange(
                *edge,
                n_internal.saturating_sub(1),
            ));
        }
        out[*edge] = *degree;
    }
    Ok(out)
}

/// Certifies convergence from edge-indexed numerator degree bounds where
/// conservative power counting suffices.
///
/// `denominator_edge_ids` is the multiset of denominator occurrences, while
/// `bounds` belongs to numerator-energy rows. The two sets are kept
/// separate because cuts and contractions can remove a denominator without
/// removing its physical energy from the numerator, and propagator powers can
/// repeat a denominator without repeating numerator ownership. The bounds
/// discard coefficients and correlations between edge energies. Consequently,
/// this can certify convergence, while failure to certify must remain
/// advisory.
pub fn energy_divergence_report(
    signatures: &[MomentumSignature],
    denominator_edge_ids: &[usize],
    bounds: &[(usize, usize)],
) -> Result<EnergyDivergenceReport> {
    let bounds = normalize_energy_degree_bounds(bounds, signatures.len())?;
    for edge_id in denominator_edge_ids {
        if *edge_id >= signatures.len() {
            return Err(EnergyBoundsError::DenominatorEdgeOutOfRange(
                *edge_id,
                signatures.len().saturating_sub(1),
            ));
        }
    }
    let denominator_rows = denominator_edge_ids
        .iter()
        .map(|edge_id| {
            signatures[*edge_id]
                .loop_signature
                .iter()
                .map(|value| i64::from(*value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let denominator_rank = rank_i64(&denominator_rows);
    for (edge_id, degree) in bounds.iter().enumerate() {
        if *degree == 0 {
            continue;
        }
        let row = signatures[edge_id]
            .loop_signature
            .iter()
            .map(|value| i64::from(*value))
            .collect::<Vec<_>>();
        if row.iter().all(|value| *value == 0) {
            continue;
        }
        let trial_rank = rank_i64(
            &denominator_rows
                .iter()
                .cloned()
                .chain(std::iter::once(row))
                .collect::<Vec<_>>(),
        );
        if trial_rank != denominator_rank {
            return Err(EnergyBoundsError::NumeratorOutsideDenominatorSpan(edge_id));
        }
    }
    let n_loops = signatures
        .first()
        .map(|signature| signature.loop_signature.len())
        .unwrap_or(0);

    let loops = (0..n_loops)
        .filter_map(|loop_id| {
            let active_edges = denominator_edge_ids
                .iter()
                .filter_map(|edge_id| {
                    (signatures[*edge_id].loop_signature[loop_id] != 0).then_some(*edge_id)
                })
                .collect::<Vec<_>>();
            if active_edges.is_empty() {
                return None;
            }
            let active_numerator_edges = signatures
                .iter()
                .enumerate()
                .filter_map(|(edge_id, signature)| {
                    (bounds[edge_id] != 0 && signature.loop_signature[loop_id] != 0)
                        .then_some(edge_id)
                })
                .collect::<Vec<_>>();
            Some(direction_report(
                Some(loop_id),
                Vec::new(),
                active_edges,
                &active_numerator_edges,
                None,
                &bounds,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let directions = arrangement_direction_reports(
        signatures,
        denominator_edge_ids,
        &denominator_rows,
        denominator_rank,
        &bounds,
    )?;
    let coordinate_convergent = loops.iter().all(|item| item.convergent);
    let directional_convergent = directions.iter().all(|item| item.convergent);

    Ok(EnergyDivergenceReport {
        edge_degree_bounds: bounds,
        loops,
        directions,
        coordinate_convergent,
        directional_convergent,
        convergent: coordinate_convergent && directional_convergent,
    })
}

pub fn auto_numerator_expr_for_bounds(
    external_count: usize,
    bounds: &[(usize, usize)],
    n_internal: usize,
) -> Result<String> {
    let bounds = normalize_energy_degree_bounds(bounds, n_internal)?;
    if bounds.iter().any(|degree| *degree != 0) && external_count == 0 {
        return Err(EnergyBoundsError::AutoNumeratorNeedsExternal);
    }

    let mut factors = Vec::new();
    for (edge_id, degree) in bounds.iter().enumerate() {
        for power_index in 0..*degree {
            let external_id = (edge_id + power_index) % external_count;
            factors.push(format!("dot(edges[{edge_id}], ext[{external_id}])"));
        }
    }
    Ok(if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join(" * ")
    })
}

fn arrangement_direction_reports(
    signatures: &[MomentumSignature],
    denominator_edge_ids: &[usize],
    rows: &[Vec<i64>],
    denominator_rank: usize,
    bounds: &[usize],
) -> Result<Vec<EnergyDirectionReport>> {
    if denominator_rank == 0 {
        return Ok(Vec::new());
    }
    let n_edges = rows.len();
    let arrangement_count = u32::try_from(n_edges)
        .ok()
        .and_then(|shift| 1usize.checked_shl(shift))
        .ok_or(EnergyBoundsError::ArrangementTooLarge {
            n_edges,
            mask_bits: usize::BITS,
        })?;
    let mut seen_active = std::collections::BTreeSet::new();
    let mut out = Vec::new();

    for mask in 0..arrangement_count {
        let zero_rows = (0..n_edges)
            .filter(|edge_id| mask & (1usize << edge_id) != 0)
            .map(|edge_id| rows[edge_id].clone())
            .collect::<Vec<_>>();
        let zero_rank = rank_i64(&zero_rows);
        let closure = (0..n_edges)
            .filter(|edge_id| {
                let trial = zero_rows
                    .iter()
                    .cloned()
                    .chain(std::iter::once(rows[*edge_id].clone()))
                    .collect::<Vec<_>>();
                rank_i64(&trial) == zero_rank
            })
            .collect::<Vec<_>>();
        let closure_rank = rank_i64(
            &closure
                .iter()
                .map(|edge_id| rows[*edge_id].clone())
                .collect_vec(),
        );
        if closure_rank >= denominator_rank {
            continue;
        }
        let active_occurrences = (0..n_edges)
            .filter(|edge_id| !closure.contains(edge_id))
            .collect::<Vec<_>>();
        if !seen_active.insert(active_occurrences.clone()) {
            continue;
        }
        let active_numerator_edges = signatures
            .iter()
            .enumerate()
            .filter_map(|(edge_id, signature)| {
                if bounds[edge_id] == 0 {
                    return None;
                }
                let trial = zero_rows
                    .iter()
                    .cloned()
                    .chain(std::iter::once(
                        signature
                            .loop_signature
                            .iter()
                            .map(|value| i64::from(*value))
                            .collect::<Vec<_>>(),
                    ))
                    .collect::<Vec<_>>();
                (rank_i64(&trial) != closure_rank).then_some(edge_id)
            })
            .collect::<Vec<_>>();
        out.push(direction_report(
            None,
            closure
                .into_iter()
                .map(|occurrence| denominator_edge_ids[occurrence])
                .collect(),
            active_occurrences
                .into_iter()
                .map(|occurrence| denominator_edge_ids[occurrence])
                .collect(),
            &active_numerator_edges,
            Some(denominator_rank - closure_rank),
            bounds,
        )?);
    }

    out.sort_by_key(|item| {
        (
            i128::try_from(item.divergence_degree).expect("isize always fits in i128")
                + i128::try_from(item.nullity.unwrap_or(1)).expect("usize fits in i128"),
            item.active_edges.len(),
            item.active_edges.clone(),
        )
    });
    out.reverse();
    Ok(out)
}

fn direction_report(
    loop_id: Option<usize>,
    zero_edges: Vec<usize>,
    active_edges: Vec<usize>,
    active_numerator_edges: &[usize],
    nullity: Option<usize>,
    bounds: &[usize],
) -> Result<EnergyDirectionReport> {
    let numerator_degree_bound = active_numerator_edges
        .iter()
        .try_fold(0usize, |degree, edge_id| {
            degree.checked_add(bounds[*edge_id])
        })
        .ok_or(EnergyBoundsError::DegreeOverflow)?;
    let denominator_degree = active_edges
        .len()
        .checked_mul(2)
        .ok_or(EnergyBoundsError::DegreeOverflow)?;
    let divergence_degree = isize::try_from(numerator_degree_bound)
        .ok()
        .and_then(|numerator_degree| {
            isize::try_from(denominator_degree)
                .ok()
                .and_then(|denominator_degree| numerator_degree.checked_sub(denominator_degree))
        })
        .ok_or(EnergyBoundsError::DegreeOverflow)?;
    let convergent = numerator_degree_bound
        .checked_add(nullity.unwrap_or(1))
        .is_some_and(|degree| degree < denominator_degree);
    Ok(EnergyDirectionReport {
        loop_id,
        zero_edges,
        active_edges,
        nullity,
        numerator_degree_bound,
        denominator_degree,
        divergence_degree,
        convergent,
    })
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn energy_uv_check_scans_noncoordinate_loop_directions() {
        let report = energy_divergence_report(
            &[
                MomentumSignature {
                    loop_signature: vec![1, 0],
                    external_signature: vec![],
                },
                MomentumSignature {
                    loop_signature: vec![0, 1],
                    external_signature: vec![],
                },
                MomentumSignature {
                    loop_signature: vec![1, 1],
                    external_signature: vec![],
                },
            ],
            &[0, 1, 2],
            &[(0, 2), (1, 2)],
        )
        .unwrap();

        assert!(report.coordinate_convergent);
        assert!(!report.directional_convergent);
        assert!(!report.convergent);
    }

    #[test]
    fn arrangement_certificate_uses_nullity_as_radial_dimension() {
        let report = energy_divergence_report(
            &[
                MomentumSignature {
                    loop_signature: vec![1, 0],
                    external_signature: vec![],
                },
                MomentumSignature {
                    loop_signature: vec![0, 1],
                    external_signature: vec![],
                },
            ],
            &[0, 1],
            &[(0, 1), (1, 1)],
        )
        .unwrap();
        let joint_scaling = report
            .directions
            .iter()
            .find(|direction| direction.nullity == Some(2))
            .unwrap();

        assert_eq!(joint_scaling.numerator_degree_bound, 2);
        assert_eq!(joint_scaling.denominator_degree, 4);
        assert_eq!(joint_scaling.divergence_degree, -2);
        assert!(!joint_scaling.convergent);
    }

    #[test]
    fn convergence_certificate_rejects_unrepresentable_degree() {
        let result = energy_divergence_report(
            &[MomentumSignature {
                loop_signature: vec![1],
                external_signature: vec![],
            }],
            &[0],
            &[(0, isize::MAX as usize + 1)],
        );

        assert!(matches!(result, Err(EnergyBoundsError::DegreeOverflow)));
    }

    #[test]
    fn convergence_certificate_rejects_oversized_arrangement_mask() {
        let signatures = (0..usize::BITS)
            .map(|_| MomentumSignature {
                loop_signature: vec![1],
                external_signature: vec![],
            })
            .collect::<Vec<_>>();
        let denominator_edge_ids = (0..signatures.len()).collect::<Vec<_>>();
        let result = energy_divergence_report(&signatures, &denominator_edge_ids, &[]);

        assert!(matches!(
            result,
            Err(EnergyBoundsError::ArrangementTooLarge {
                n_edges,
                mask_bits
            }) if n_edges == usize::BITS as usize && mask_bits == usize::BITS
        ));
    }

    #[test]
    fn normal_box_bounds_are_convergent_for_affine_numerators() {
        let signatures = crate::graph_io::test_graphs::box_graph()
            .internal_edges
            .into_iter()
            .map(|edge| edge.signature)
            .collect::<Vec<_>>();
        let denominator_edge_ids = (0..signatures.len()).collect::<Vec<_>>();
        let report =
            energy_divergence_report(&signatures, &denominator_edge_ids, &[(0, 1), (1, 1)])
                .unwrap();
        assert!(report.convergent);
    }

    #[test]
    fn certificate_separates_denominator_occurrences_from_numerator_edges() {
        let signatures = [
            MomentumSignature {
                loop_signature: vec![1, 0],
                external_signature: vec![],
            },
            MomentumSignature {
                loop_signature: vec![1, 0],
                external_signature: vec![],
            },
            MomentumSignature {
                loop_signature: vec![0, 0],
                external_signature: vec![],
            },
        ];
        let report = energy_divergence_report(&signatures, &[0, 0], &[(1, 2), (2, 100)]).unwrap();

        assert_eq!(report.loops.len(), 1);
        assert_eq!(report.loops[0].loop_id, Some(0));
        assert_eq!(report.loops[0].active_edges, vec![0, 0]);
        assert_eq!(report.loops[0].numerator_degree_bound, 2);
        assert_eq!(report.loops[0].denominator_degree, 4);
        assert_eq!(report.directions[0].nullity, Some(1));
    }

    #[test]
    fn certificate_rejects_bounded_numerator_outside_denominator_span() {
        let signatures = [
            MomentumSignature {
                loop_signature: vec![1, 0],
                external_signature: vec![],
            },
            MomentumSignature {
                loop_signature: vec![0, 1],
                external_signature: vec![],
            },
        ];

        assert!(matches!(
            energy_divergence_report(&signatures, &[0], &[(1, 1)]),
            Err(EnergyBoundsError::NumeratorOutsideDenominatorSpan(1))
        ));
    }

    #[test]
    fn auto_numerator_expr_saturates_sparse_energy_bounds() {
        let expr = auto_numerator_expr_for_bounds(3, &[(0, 1), (1, 1), (3, 4)], 6).unwrap();
        assert_eq!(
            expr,
            "dot(edges[0], ext[0]) * dot(edges[1], ext[1]) * dot(edges[3], ext[0]) * dot(edges[3], ext[1]) * dot(edges[3], ext[2]) * dot(edges[3], ext[0])"
        );
    }
}
