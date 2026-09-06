use thiserror::Error;

#[derive(Debug, Error)]
pub enum EnergyBoundsError {
    #[error("energy-degree bound edge index {0} is outside 0..{1}")]
    EdgeOutOfRange(usize, usize),
    #[error("auto numerator generation requires at least one external momentum symbol")]
    AutoNumeratorNeedsExternal,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_numerator_expr_saturates_sparse_energy_bounds() {
        let expr = auto_numerator_expr_for_bounds(3, &[(0, 1), (1, 1), (3, 4)], 6).unwrap();
        assert_eq!(
            expr,
            "dot(edges[0], ext[0]) * dot(edges[1], ext[1]) * dot(edges[3], ext[0]) * dot(edges[3], ext[1]) * dot(edges[3], ext[2]) * dot(edges[3], ext[0])"
        );
    }
}
