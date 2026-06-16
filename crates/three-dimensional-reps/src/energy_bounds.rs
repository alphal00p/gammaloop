use itertools::Itertools;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{MomentumSignature, utils::rank_i64};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnergyDivergenceReport {
    pub edge_degree_bounds: Vec<usize>,
    /// Coordinate-axis scaling certificates, one per loop energy.
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
    pub zero_edges: Vec<usize>,
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
/// The bounds discard coefficients and correlations between edge energies.
/// Consequently, this can certify convergence, while failure to certify must
/// remain advisory.
pub fn energy_divergence_report(
    signatures: &[MomentumSignature],
    bounds: &[(usize, usize)],
) -> Result<EnergyDivergenceReport> {
    let bounds = normalize_energy_degree_bounds(bounds, signatures.len())?;
    let n_loops = signatures
        .first()
        .map(|signature| signature.loop_signature.len())
        .unwrap_or(0);

    let loops = (0..n_loops)
        .map(|loop_id| {
            let active_edges = signatures
                .iter()
                .enumerate()
                .filter_map(|(edge_id, signature)| {
                    (signature.loop_signature[loop_id] != 0).then_some(edge_id)
                })
                .collect::<Vec<_>>();
            direction_report(Some(loop_id), Vec::new(), active_edges, None, &bounds)
        })
        .collect::<Result<Vec<_>>>()?;
    let directions = arrangement_direction_reports(signatures, &bounds)?;
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
    bounds: &[usize],
) -> Result<Vec<EnergyDirectionReport>> {
    if signatures.is_empty() {
        return Ok(Vec::new());
    }
    let rows = signatures
        .iter()
        .map(|signature| {
            signature
                .loop_signature
                .iter()
                .map(|value| i64::from(*value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let n_loops = rows[0].len();
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
        if closure_rank >= n_loops {
            continue;
        }
        let active_edges = (0..n_edges)
            .filter(|edge_id| !closure.contains(edge_id))
            .collect::<Vec<_>>();
        if !seen_active.insert(active_edges.clone()) {
            continue;
        }
        out.push(direction_report(
            None,
            closure,
            active_edges,
            Some(n_loops - closure_rank),
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
    nullity: Option<usize>,
    bounds: &[usize],
) -> Result<EnergyDirectionReport> {
    let numerator_degree_bound = active_edges
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
        let result = energy_divergence_report(&signatures, &[]);

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
        let report = energy_divergence_report(&signatures, &[(0, 1), (1, 1)]).unwrap();
        assert!(report.convergent);
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
