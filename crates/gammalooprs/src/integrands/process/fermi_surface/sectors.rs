use std::collections::BTreeMap;

use color_eyre::eyre::{Result, ensure, eyre};
use symbolica::atom::{Atom, AtomCore, AtomView};
use three_dimensional_reps::ThermalDistributionFactor;

use crate::{
    graph::Graph, integrands::process::param_builder::ThermalDistributionCall,
    momentum::SignOrZero, numerator::symbolica_ext::NumeratorAtomExt, utils::GS, uv::uv_graph::UVE,
};

use super::FermiSurfaceProduct;

/// A coefficient and its independent, oriented zero-temperature delta support.
/// Numerator functions remain opaque; their existing definitions accompany the
/// coefficient into evaluator construction.
pub(crate) struct FermiSurfaceSector {
    pub(crate) product: FermiSurfaceProduct,
    pub(crate) orientations: Vec<i32>,
    pub(crate) coefficient: Atom,
}

impl FermiSurfaceSector {
    fn calls(expression: &Atom) -> Result<Vec<(Atom, ThermalDistributionCall)>> {
        let mut calls = Vec::new();
        let mut error = None;
        expression.visitor(&mut |atom| {
            if error.is_some() {
                return false;
            }
            if let AtomView::Fun(function) = atom
                && function.get_symbol() == GS.thermal_weight_wrapper
                && function.get_nargs() != 1
            {
                error = Some(eyre!(
                    "Thermal weight wrapper must have one argument, got {atom}"
                ));
                return false;
            }
            if let AtomView::Fun(function) = atom
                && function.get_symbol() == GS.thermal_distribution
            {
                let call = ThermalDistributionCall::try_from(atom).and_then(|call| {
                    ensure!(
                        call.temperature_flag.is_zero() || call.temperature_flag.is_one(),
                        "Thermal distribution temperature flag must be zero or one, got {}",
                        call.temperature_flag
                    );
                    Ok(call)
                });
                match call {
                    Ok(call) => {
                        if !calls
                            .iter()
                            .any(|(key, _): &(Atom, _)| key.as_view() == atom)
                        {
                            calls.push((atom.to_owned(), call));
                        }
                    }
                    Err(parse_error) => error = Some(parse_error),
                }
            }
            true
        });
        if let Some(error) = error {
            return Err(error);
        }
        Ok(calls)
    }

    /// Detect distribution derivatives before resolving physical orientations.
    pub(crate) fn is_present(expression: &Atom) -> Result<bool> {
        Ok(Self::calls(expression)?
            .iter()
            .any(|(_, call)| call.temperature_flag.is_zero() && call.derivative_order > 0))
    }

    /// Separate regular terms and Fermi supports after residue and physical
    /// orientation selection. The fifth N argument must be a definite sign;
    /// distinct signs in an explicit sum remain distinct supports.
    pub(crate) fn extract(graph: &Graph, expression: &Atom) -> Result<(Atom, Vec<Self>)> {
        let calls = Self::calls(expression)?;
        let mut expression = expression.unwrap_function(GS.thermal_weight_wrapper);
        let mut active = Vec::new();
        for (key, call) in calls {
            if !call.temperature_flag.is_zero() || call.derivative_order == 0 {
                continue;
            }
            let (_, _, edge) = graph
                .iter_edges()
                .find(|(_, edge, _)| *edge == call.edge)
                .ok_or_else(|| eyre!("Unknown thermal distribution edge {}", call.edge))?;
            if edge.data.chemical_potential_atom().is_none() {
                // At zero temperature an edge without a chemical potential has
                // constant weight, so every positive-order derivative vanishes.
                expression = expression.replace(key.to_pattern()).with(Atom::Zero);
                continue;
            }
            let sign = i64::try_from(call.thermal_sign.as_view())
                .map_err(|_| eyre!("Fermi-surface thermal sign must be resolved in {key}"))?;
            let orientation = i64::try_from(call.orientation_sign.as_view()).map_err(|_| {
                eyre!("Fermi-surface chemical orientation must be resolved in {key}")
            })?;
            ensure!(
                matches!(sign, -1 | 1) && matches!(orientation, -1 | 1),
                "Fermi-surface thermal sign and chemical orientation must be +1 or -1 in {key}"
            );
            active.push((
                key,
                ThermalDistributionFactor {
                    edge_id: call.edge,
                    sign: sign as i32,
                    derivative_order: call.derivative_order,
                },
                orientation as i32,
            ));
        }
        if active.is_empty() || expression.is_zero() {
            return Ok((expression, Vec::new()));
        }
        let contains_active = |expression: AtomView<'_>| {
            let mut found = false;
            expression.visitor(&mut |part| {
                found |= active.iter().any(|(key, _, _)| key.as_view() == part);
                !found
            });
            found
        };
        let mut invalid = None;
        expression.visitor(&mut |part| {
            if invalid.is_some() || active.iter().any(|(key, _, _)| key.as_view() == part) {
                return false;
            }
            if matches!(part, AtomView::Fun(_) | AtomView::Pow(_)) && contains_active(part) {
                invalid = Some(part.to_owned());
                return false;
            }
            true
        });
        ensure!(
            invalid.is_none(),
            "Fermi distributions must occur polynomially outside functions, inverses and powers: {}",
            invalid.unwrap_or(Atom::Zero)
        );

        let keys = active
            .iter()
            .map(|(key, _, _)| key.clone())
            .collect::<Vec<_>>();
        let mut bulk = Atom::Zero;
        let mut grouped = BTreeMap::<Vec<(ThermalDistributionFactor, i32)>, Atom>::new();
        // Collect only tagged distributions: graph numerator bodies and their
        // factorized symbolic references never enter polynomial expansion.
        for (key, coefficient) in expression.coefficient_list::<u32>(&keys) {
            ensure!(
                !contains_active(coefficient.as_view()),
                "Fermi distribution remains hidden in coefficient {coefficient}"
            );
            if key.is_one() {
                bulk += coefficient;
                continue;
            }
            let mut support = active
                .iter()
                .filter(|(candidate, _, _)| {
                    let mut found = false;
                    key.visitor(&mut |part| {
                        found |= candidate.as_view() == part;
                        !found
                    });
                    found
                })
                .map(|(_, factor, orientation)| (*factor, *orientation))
                .collect::<Vec<_>>();
            support.sort();
            *grouped.entry(support).or_insert(Atom::Zero) += coefficient;
        }
        let mut sectors = Vec::with_capacity(grouped.len());
        for (support, coefficient) in grouped {
            if coefficient.is_zero() {
                continue;
            }
            let (factors, orientations): (Vec<_>, Vec<_>) = support.into_iter().unzip();
            let product = FermiSurfaceProduct::new(graph, &factors)?;
            ensure!(
                !coefficient.contains_symbol(GS.heaviside),
                "An explicit step in a Fermi coefficient has no certified independent routing"
            );
            for (_, step) in Self::calls(&coefficient)? {
                if !step.temperature_flag.is_zero() || step.derivative_order != 0 {
                    continue;
                }
                let (_, _, edge) = graph
                    .iter_edges()
                    .find(|(_, edge, _)| *edge == step.edge)
                    .ok_or_else(|| eyre!("Unknown thermal distribution edge {}", step.edge))?;
                if edge.data.chemical_potential_atom().is_none() {
                    continue;
                }
                let routing = product
                    .lmb()
                    .edge_signatures
                    .get(step.edge)
                    .ok_or_else(|| {
                        eyre!("Missing routing for thermal step on edge {}", step.edge)
                    })?;
                ensure!(
                    routing
                        .internal
                        .iter()
                        .take(factors.len())
                        .all(|sign| *sign == SignOrZero::Zero),
                    "Thermal step on edge {} depends on an active Fermi momentum; its boundary intersection is unsupported in production localization",
                    step.edge
                );
            }
            sectors.push(Self {
                product,
                orientations,
                coefficient,
            });
        }
        Ok((bulk, sectors))
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{function, symbol};

    use super::*;
    use crate::integrands::process::fermi_surface::routing;

    #[test]
    fn fermi_sector_extraction_preserves_factorized_coefficients_and_regular_terms() -> Result<()> {
        let graph = routing::test_graph()?;
        let x = Atom::var(symbol!("fermi_sector_x"));
        let y = Atom::var(symbol!("fermi_sector_y"));
        let numerator = function!(symbol!("fermi_sector_numerator"), &x, &y);
        let coefficient = numerator * (&x + Atom::num(1)) * (&y + Atom::num(2));
        let independent_step = GS.thermal_distribution(3, 0, 0, -1, -1);
        let bulk = Atom::num(7) * GS.thermal_distribution(1, 3, 1, 1, 1);
        for derivative_order in [2, 3] {
            let delta = GS.thermal_distribution(1, derivative_order as i64, 0, -1, 1);
            let expression = &bulk
                + function!(GS.thermal_weight_wrapper, &delta) * &coefficient * &independent_step;
            let (actual_bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
            assert_eq!(actual_bulk, bulk);
            assert_eq!(sectors.len(), 1);
            assert_eq!(sectors[0].coefficient, &coefficient * &independent_step);
            assert_eq!(sectors[0].orientations, [1]);
            assert_eq!(
                sectors[0].product.factors(),
                &[ThermalDistributionFactor {
                    edge_id: linnet::half_edge::involution::EdgeIndex(1),
                    sign: -1,
                    derivative_order,
                }]
            );
            assert_eq!(sectors[0].product.derivative_orders, [derivative_order - 1]);
            assert_eq!(sectors[0].product.shape.len(), derivative_order);
            assert_eq!(
                actual_bulk + &sectors[0].coefficient * delta,
                expression.unwrap_function(GS.thermal_weight_wrapper)
            );
        }
        Ok(())
    }

    #[test]
    fn fermi_sectors_keep_opposite_chemical_supports_and_independent_products() -> Result<()> {
        let graph = routing::test_graph()?;
        let plus = GS.thermal_distribution(1, 1, 0, 1, 1);
        let minus = GS.thermal_distribution(1, 1, 0, 1, -1);
        let second = GS.thermal_distribution(3, 2, 0, -1, -1);
        let expression =
            Atom::num(2) * &plus + Atom::num(3) * &minus + Atom::num(5) * &plus * &second;
        let (bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
        assert!(bulk.is_zero());
        assert_eq!(sectors.len(), 3);
        let reconstructed = sectors.iter().fold(bulk, |sum, sector| {
            sum + sector
                .product
                .factors()
                .iter()
                .zip(&sector.orientations)
                .fold(
                    sector.coefficient.clone(),
                    |coefficient, (factor, orientation)| {
                        coefficient
                            * GS.thermal_distribution(
                                factor.edge_id.0 as i64,
                                factor.derivative_order as i64,
                                0,
                                factor.sign,
                                *orientation,
                            )
                    },
                )
        });
        assert_eq!(reconstructed, expression);
        assert!(sectors.iter().any(|sector| sector.orientations == [1, -1]));
        Ok(())
    }

    #[test]
    fn fermi_sector_detector_accepts_unresolved_orientation_and_ignores_regular_weights()
    -> Result<()> {
        let graph = routing::test_graph()?;
        let unresolved = GS.thermal_distribution(
            1,
            1,
            0,
            1,
            GS.sign(linnet::half_edge::involution::EdgeIndex(1)),
        );
        assert!(FermiSurfaceSector::is_present(&unresolved)?);
        assert!(FermiSurfaceSector::extract(&graph, &unresolved).is_err());
        for expression in [
            GS.thermal_distribution(1, 0, 0, 1, 1),
            GS.thermal_distribution(1, 2, 1, 1, 1),
        ] {
            assert!(!FermiSurfaceSector::is_present(&expression)?);
            let (bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
            assert_eq!(bulk, expression);
            assert!(sectors.is_empty());
        }
        for derivative_order in [1, 3] {
            let (bulk, sectors) = FermiSurfaceSector::extract(
                &graph,
                &GS.thermal_distribution(4, derivative_order, 0, 1, 1),
            )?;
            assert!(bulk.is_zero());
            assert!(sectors.is_empty());
        }
        Ok(())
    }

    #[test]
    fn fermi_sector_extraction_rejects_undefined_products_and_uncertified_steps() -> Result<()> {
        let graph = routing::test_graph()?;
        let delta = GS.thermal_distribution(1, 1, 0, 1, 1);
        for malformed in [
            delta.clone().pow(2),
            delta.clone().pow(-1),
            function!(symbol!("fermi_sector_hidden"), &delta),
            &delta * GS.thermal_distribution(0, 1, 0, 1, 1),
            &delta * GS.thermal_distribution(1, 2, 0, 1, 1),
            &delta * GS.thermal_distribution(1, 0, 0, 1, 1),
            &delta * GS.thermal_distribution(0, 0, 0, 1, 1),
            &delta * GS.heaviside(Atom::var(symbol!("fermi_sector_hidden_step"))),
            GS.thermal_distribution(99, 1, 0, 1, 1),
            GS.thermal_distribution(1, 1, 0, 0, 1),
            GS.thermal_distribution(1, 1, 0, 1, 0),
            GS.thermal_distribution(1, -1, 0, 1, 1),
            GS.thermal_distribution(1, 1, 2, 1, 1),
            function!(GS.thermal_distribution, 1, 1, 0, 1),
        ] {
            assert!(
                FermiSurfaceSector::extract(&graph, &malformed).is_err(),
                "accepted {malformed}"
            );
        }
        Ok(())
    }
}
