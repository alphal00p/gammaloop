use std::collections::BTreeMap;

use color_eyre::eyre::{Result, ensure, eyre};
use linnet::half_edge::involution::EdgeIndex;
use symbolica::atom::{Atom, AtomCore, Indeterminate, Symbol};
use three_dimensional_reps::ThermalDistributionFactor;

use crate::{graph::Graph, utils::GS, uv::uv_graph::UVE};

use super::FermiSurfaceProduct;

/// One smooth coefficient times explicit thermal distribution factors.
///
/// Order zero remains a step; order r>0 represents delta^(r-1). Differentiating
/// the factors symbolically exposes boundary and intersection terms before any
/// smooth numerical callback is built. This does not assemble production sectors.
#[derive(Clone, Debug)]
pub struct ThermalBoundaryTerm {
    coefficient: Atom,
    factors: Vec<ThermalDistributionFactor>,
}

impl ThermalBoundaryTerm {
    pub fn new(
        graph: &Graph,
        coefficient: Atom,
        factors: &[ThermalDistributionFactor],
    ) -> Result<Self> {
        ensure!(
            !coefficient.contains_symbol(GS.thermal_distribution)
                && !coefficient.contains_symbol(GS.heaviside),
            "Thermal boundary coefficient must keep distributions in its explicit factor list"
        );
        let mut factors = factors.to_vec();
        factors.sort();
        let mut occurrences = BTreeMap::<EdgeIndex, (usize, bool)>::new();
        for factor in &factors {
            let (count, active) = occurrences.entry(factor.edge_id).or_default();
            *count += 1;
            *active |= factor.derivative_order > 0;
            ensure!(
                *count == 1 || !*active,
                "Coincident thermal factors on edge {} require a regulated product; a step cannot multiply a delta on the same edge",
                factor.edge_id
            );
            // Reuse the localizer's graph, particle, sign and routing checks.
            // An undifferentiated step is checked as a potential single shell.
            if factor.derivative_order == 0 {
                FermiSurfaceProduct::new(
                    graph,
                    &[ThermalDistributionFactor {
                        derivative_order: 1,
                        ..*factor
                    }],
                )?;
            }
        }
        let active = factors
            .iter()
            .filter(|factor| factor.derivative_order > 0)
            .copied()
            .collect::<Vec<_>>();
        if !active.is_empty() {
            FermiSurfaceProduct::new(graph, &active)?;
            // Raised-edge groups already certify full routing up to sign and
            // equal masses, including every fixed external shift. The thermal
            // sign does not resolve sigma, so equal chemical potentials up to
            // sign can still put a step and a delta on the same physical shell.
            for group in graph.get_raised_edge_groups() {
                for delta in active
                    .iter()
                    .filter(|factor| group.contains(&factor.edge_id))
                {
                    for step in factors.iter().filter(|factor| {
                        factor.derivative_order == 0 && group.contains(&factor.edge_id)
                    }) {
                        let delta_mu = graph[delta.edge_id].chemical_potential_atom();
                        let step_mu = graph[step.edge_id].chemical_potential_atom();
                        ensure!(
                            delta_mu != step_mu && delta_mu != step_mu.map(|mu| -mu),
                            "Potentially coincident thermal shells on edges {} and {} require resolved orientations or a regulated product; a step cannot multiply a delta on the same physical shell",
                            step.edge_id,
                            delta.edge_id
                        );
                    }
                }
            }
        }
        Ok(Self {
            coefficient,
            factors,
        })
    }

    pub fn coefficient(&self) -> &Atom {
        &self.coefficient
    }

    pub fn factors(&self) -> &[ThermalDistributionFactor] {
        &self.factors
    }

    /// Differentiate the smooth coefficient and each thermal factor by Leibniz.
    ///
    /// The coefficient and supplied positive energies must already be explicit
    /// expressions in `variable`; chemical potentials, orientation and masses
    /// are held fixed. Each factor contributes E'(x) N^(r+1)(E(x)), without an
    /// extra thermal sign. Repeated calls also differentiate the generated E'
    /// coefficients and therefore include nonlinear chain-rule terms.
    ///
    /// When applying the normal derivative of a localized delta, pass only its
    /// complete coefficient (including the host profile, volume and Jacobian)
    /// and residual thermal factors here. The host delta is not part of that
    /// coefficient. Before localizing a generated intersection, validate
    /// the union of its active factors and the host with `FermiSurfaceProduct::new`.
    /// Remaining order-zero factors stay symbolic until their normal dependence
    /// has been resolved; they must not be differentiated in a numeric callback.
    pub fn differentiate(
        &self,
        graph: &Graph,
        variable: Symbol,
        energies: &BTreeMap<EdgeIndex, Atom>,
    ) -> Result<Vec<Self>> {
        if self.coefficient.is_zero() {
            return Ok(Vec::new());
        }
        let variable = Indeterminate::from(variable);
        let mut terms =
            BTreeMap::from([(self.factors.clone(), self.coefficient.derivative(&variable))]);
        for (index, factor) in self.factors.iter().enumerate() {
            let energy = energies.get(&factor.edge_id).ok_or_else(|| {
                eyre!(
                    "Missing explicit energy for thermal boundary edge {}",
                    factor.edge_id
                )
            })?;
            let energy_derivative = energy.derivative(&variable);
            if energy_derivative.is_zero() {
                continue;
            }
            let mut factors = self.factors.clone();
            factors[index].derivative_order = factor
                .derivative_order
                .checked_add(1)
                .ok_or_else(|| eyre!("Thermal boundary derivative order overflow"))?;
            factors.sort();
            let coefficient = &self.coefficient * energy_derivative;
            let total = terms.entry(factors).or_insert(Atom::Zero);
            *total = &*total + coefficient;
        }
        terms
            .into_iter()
            .filter(|(_, coefficient)| !coefficient.is_zero())
            .map(|(factors, coefficient)| Self::new(graph, coefficient, &factors))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{domains::dual::HyperDual, prelude::Real, symbol};

    use super::*;
    use crate::{
        dot,
        graph::parse::IntoGraph,
        integrands::process::fermi_surface::routing,
        momentum::{ThreeMomentum, sample::LoopMomenta},
        utils::{F, hyperdual_utils::new_constant},
    };
    use spenso::algebra::complex::Complex;

    #[test]
    fn thermal_boundary_nonlinear_second_derivative_keeps_step_and_delta_sectors() -> Result<()> {
        let graph = routing::test_graph()?;
        let variable = symbol!("fermi_boundary_nonlinear_x");
        let x = Atom::var(variable);
        let energies = BTreeMap::from([(EdgeIndex(1), x.clone().pow(2) + &x + Atom::num(7))]);
        for sign in [-1, 0, 1] {
            let term = ThermalBoundaryTerm::new(
                &graph,
                x.clone().pow(2),
                &[ThermalDistributionFactor {
                    edge_id: EdgeIndex(1),
                    sign,
                    derivative_order: 0,
                }],
            )?;
            let mut coefficients = BTreeMap::<usize, Atom>::new();
            for first in term.differentiate(&graph, variable, &energies)? {
                for second in first.differentiate(&graph, variable, &energies)? {
                    assert_eq!(second.factors()[0].sign, sign);
                    let total = coefficients
                        .entry(second.factors()[0].derivative_order)
                        .or_insert(Atom::Zero);
                    *total = &*total + second.coefficient();
                }
            }
            // For c=x^2 and E=x^2+x+7:
            // (c N(E))''=c''N+(2c'E'+cE'')N'+c(E')^2 N''.
            let expected = [
                Atom::num(2),
                Atom::num(10) * x.clone().pow(2) + Atom::num(4) * &x,
                x.clone().pow(2) * (Atom::num(2) * &x + Atom::num(1)).pow(2),
            ];
            assert_eq!(coefficients.len(), expected.len());
            for (order, expected) in expected.into_iter().enumerate() {
                assert_eq!(
                    (&coefficients[&order] - expected).expand(),
                    Atom::Zero,
                    "order={order}, thermal sign={sign}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn thermal_boundary_mixed_derivative_generates_a_localizable_intersection() -> Result<()> {
        let graph = routing::test_graph()?;
        let x_symbol = symbol!("fermi_boundary_mixed_x");
        let y_symbol = symbol!("fermi_boundary_mixed_y");
        let x = Atom::var(x_symbol);
        let y = Atom::var(y_symbol);
        let energies = BTreeMap::from([
            (EdgeIndex(1), &x + &y),
            (EdgeIndex(3), &x + Atom::num(2) * &y),
        ]);
        let term = ThermalBoundaryTerm::new(
            &graph,
            Atom::num(1),
            &[1, 3].map(|edge| ThermalDistributionFactor {
                edge_id: EdgeIndex(edge),
                sign: -1,
                derivative_order: 0,
            }),
        )?;
        let mut coefficients = BTreeMap::<Vec<ThermalDistributionFactor>, Atom>::new();
        for first in term.differentiate(&graph, x_symbol, &energies)? {
            for second in first.differentiate(&graph, y_symbol, &energies)? {
                let total = coefficients
                    .entry(second.factors().to_vec())
                    .or_insert(Atom::Zero);
                *total = &*total + second.coefficient();
            }
        }
        assert_eq!(coefficients.len(), 3);
        for (orders, expected) in [([2, 0], 1), ([1, 1], 3), ([0, 2], 2)] {
            let factors = [1, 3]
                .into_iter()
                .zip(orders)
                .map(|(edge, derivative_order)| ThermalDistributionFactor {
                    edge_id: EdgeIndex(edge),
                    sign: -1,
                    derivative_order,
                })
                .collect::<Vec<_>>();
            assert_eq!(coefficients[&factors], Atom::num(expected));
        }
        let (factors, coefficient) = coefficients
            .into_iter()
            .find(|(factors, _)| factors.iter().all(|factor| factor.derivative_order == 1))
            .unwrap();
        let product = FermiSurfaceProduct::new(&graph, &factors)?;
        let coefficient = coefficient
            .evaluator::<Atom>(&[])
            .build()?
            .map_coeff(&|value| value.re.to_f64())
            .evaluate_single(&[]);
        let actual = product.localize(
            &LoopMomenta(vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0)); 2]),
            &[F(0.0), F(0.0)],
            &[F(1.0), F(2.0)],
            |scale| Ok((-scale.clone()).exp()),
            |momenta| {
                let value = new_constant(&momenta.0[0].px, &F(coefficient));
                Ok(HyperDual::from_values(
                    product.shape.clone(),
                    value.values.into_iter().map(Complex::new_re).collect(),
                ))
            },
        )?;
        // Each generated ordinary delta gives E^3 exp(-E) for |k|=1.
        assert!((actual.re.0 - 24.0 * (-3.0_f64).exp()).abs() < 3e-15);
        assert_eq!(actual.im, F(0.0));
        Ok(())
    }

    #[test]
    fn thermal_boundary_rejects_missing_energies_and_invalid_promoted_products() -> Result<()> {
        let graph = routing::test_graph()?;
        let variable = symbol!("fermi_boundary_validation_x");
        let step = ThermalDistributionFactor {
            edge_id: EdgeIndex(1),
            sign: 1,
            derivative_order: 0,
        };
        let single = ThermalBoundaryTerm::new(&graph, Atom::num(1), &[step])?;
        assert!(
            single
                .differentiate(&graph, variable, &BTreeMap::new())
                .is_err()
        );
        let energies = BTreeMap::from([(EdgeIndex(1), Atom::var(variable))]);
        let repeated = ThermalBoundaryTerm::new(&graph, Atom::num(1), &[step, step])?;
        let error = repeated
            .differentiate(&graph, variable, &energies)
            .unwrap_err();
        assert!(error.to_string().contains("Coincident thermal factors"));
        assert!(
            ThermalBoundaryTerm::new(
                &graph,
                Atom::num(1),
                &[
                    step,
                    ThermalDistributionFactor {
                        derivative_order: 1,
                        ..step
                    }
                ],
            )
            .is_err()
        );
        let dependent = ThermalBoundaryTerm::new(
            &graph,
            Atom::num(1),
            &[
                ThermalDistributionFactor {
                    edge_id: EdgeIndex(0),
                    derivative_order: 1,
                    ..step
                },
                step,
            ],
        )?;
        let energies = BTreeMap::from([
            (EdgeIndex(0), Atom::var(variable)),
            (EdgeIndex(1), Atom::var(variable)),
        ]);
        let error = dependent
            .differentiate(&graph, variable, &energies)
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("independent loop-momentum routes")
        );
        assert!(ThermalBoundaryTerm::new(&graph, GS.heaviside(Atom::var(variable)), &[]).is_err());
        assert!(
            ThermalBoundaryTerm::new(
                &graph,
                Atom::num(1),
                &[ThermalDistributionFactor { sign: 2, ..step }]
            )
            .is_err()
        );
        let constant = ThermalBoundaryTerm::new(&graph, Atom::num(1), &[])?;
        assert!(
            constant
                .differentiate(&graph, variable, &BTreeMap::new())?
                .is_empty()
        );
        Ok(())
    }

    #[test]
    fn thermal_boundary_rejects_coincident_shells_on_distinct_chain_edges() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let graph: Graph = dot!(digraph coincident_fermi_chain {
            node [num=1]
            edge [num=1 particle="d"]
            A -> B [id=0]
            B -> A [id=1]
        })?;
        let variable = symbol!("fermi_boundary_coincident_x");
        for sign in [-1, 0, 1] {
            let error = ThermalBoundaryTerm::new(
                &graph,
                Atom::var(variable),
                &[
                    ThermalDistributionFactor {
                        edge_id: EdgeIndex(0),
                        sign,
                        derivative_order: 0,
                    },
                    ThermalDistributionFactor {
                        edge_id: EdgeIndex(1),
                        sign: 1,
                        derivative_order: 1,
                    },
                ],
            )
            .unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains("Potentially coincident thermal shells")
            );
        }
        Ok(())
    }
}
