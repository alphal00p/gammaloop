//! Localization of zero-temperature distribution derivatives on independent cyclic chains.
//!
//! Each factor `N^(r)(E)` becomes `delta^(r-1)(E - sigma*mu)`. The CFF coefficient
//! already owns the cyclic-chain sign and factorial. This module supplies the
//! distribution action and spatial measure, and evaluates compiled smooth
//! coefficients in the amplitude's existing momentum and parameter layout.
//!
//! In a basis containing one edge from each chain, insert a normalized positive
//! scale integral and dilate that edge momentum, `q = t*k`. With shell offset
//! `u = E - sigma*mu`, the inverse map is analytic:
//! `t(u) = sqrt((sigma*mu + u)^2 - m^2) / |k|`. Its complete weight is
//! `h(t) * t^3 / g'(t) = h(t) * E * (E^2 - m^2) / |k|^4`.
//! Evaluating this weight and the coefficient on energy-offset jets implements
//! the same change-of-variable derivative as a raised-cut residue, including
//! every derivative of the profile, Jacobian and coefficient.

use bincode_trait_derive::{Decode, Encode};
use color_eyre::eyre::{Result, eyre};
use spenso::algebra::complex::Complex;
use symbolica::{
    domains::{
        dual::{DualNumberStructure, HyperDual},
        float::SingleFloat,
    },
    prelude::Real,
};
use three_dimensional_reps::thermal::ThermalDistributionFactor;

use crate::{
    graph::LoopMomentumBasis,
    integrands::process::sampling::maps::SamplingEvaluationError,
    momentum::sample::LoopMomenta,
    utils::{F, FloatLike, Length, hyperdual_utils::new_constant},
};

mod boundary;
mod evaluation;
mod routing;
mod sectors;

pub use boundary::ThermalBoundaryTerm;
pub(crate) use evaluation::FermiSurfaceEvaluator;
pub(crate) use sectors::FermiSurfaceSector;

/// A product of Fermi distributions whose cyclic-chain momenta are independent.
///
/// The first loop coordinates of `lmb()` follow `factors()` in order. External
/// shifts are retained in the basis routing, so the selected coordinates are
/// the actual edge momenta about which the Fermi spheres are centred.
#[derive(Clone, Debug, Encode, Decode)]
pub struct FermiSurfaceProduct {
    factors: Vec<ThermalDistributionFactor>,
    lmb: LoopMomentumBasis,
    derivative_orders: Vec<usize>,
    shape: Vec<Vec<usize>>,
}

impl FermiSurfaceProduct {
    /// Localize a smooth coefficient while preserving the spatial integration dimension.
    ///
    /// `loop_momenta` is in `lmb()`. Masses and oriented chemical potentials
    /// (`sigma*mu`, not their absolute values) follow `factors()` in order;
    /// only the squared mass enters the shell. `profile(t)` must integrate to
    /// one over `t > 0`.
    /// Both callbacks must return jets with the supplied shape. The coefficient
    /// includes all routed numerator/denominator dependence; the localizer owns
    /// the scale profiles, spatial dilation factors and delta Jacobians.
    ///
    /// A future conditional Cutkosky localization can be evaluated inside the
    /// coefficient callback. Its root and Jacobian must retain their dependence
    /// on all shell-offset jets before the Fermi derivatives are extracted.
    /// No delta or discontinuous step function may be numerically differentiated
    /// inside this smooth callback. Use [`ThermalBoundaryTerm`] to expand such
    /// derivatives symbolically first, and validate any generated intersection
    /// together with its host surface using `FermiSurfaceProduct::new`.
    pub fn localize<T: FloatLike>(
        &self,
        loop_momenta: &LoopMomenta<F<T>>,
        masses: &[F<T>],
        oriented_chemical_potentials: &[F<T>],
        mut profile: impl FnMut(&HyperDual<F<T>>) -> Result<HyperDual<F<T>>>,
        coefficient: impl FnOnce(&LoopMomenta<HyperDual<F<T>>>) -> Result<HyperDual<Complex<F<T>>>>,
    ) -> Result<Complex<F<T>>> {
        if loop_momenta.len() != self.lmb.loop_edges.len()
            || masses.len() != self.factors.len()
            || oriented_chemical_potentials.len() != self.factors.len()
        {
            return Err(eyre!(
                "Fermi localization requires {} loop momenta and {} masses/chemical potentials",
                self.lmb.loop_edges.len(),
                self.factors.len()
            ));
        }
        let representative = &oriented_chemical_potentials[0];
        for (index, (mass, nu)) in masses.iter().zip(oriented_chemical_potentials).enumerate() {
            if !mass.is_finite() || !nu.is_finite() {
                return Err(eyre!(
                    "Fermi surface on edge {} requires a finite mass and finite oriented chemical potential",
                    self.factors[index].edge_id
                ));
            }
        }
        let masses = masses.iter().map(|mass| mass.abs()).collect::<Vec<_>>();
        if loop_momenta.iter().any(|momentum| {
            [&momentum.px, &momentum.py, &momentum.pz]
                .into_iter()
                .any(|value| !value.is_finite())
        }) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "Fermi localization",
                detail: "Fermi localization requires finite loop momenta".to_owned(),
            }
            .into());
        }
        // Empty support is a physical zero, independent of the smooth coefficient.
        // A zero chemical potential leaves the occupation constant on E > 0,
        // also for a massless fermion, exactly as for a particle without one.
        if masses
            .iter()
            .zip(oriented_chemical_potentials)
            .any(|(mass, nu)| nu < mass || nu.is_zero())
        {
            return Ok(Complex::new_re(representative.zero()));
        }
        let mut template = HyperDual::from_values(
            self.shape.clone(),
            vec![representative.zero(); self.shape.len()],
        );
        template.values[0] = representative.one();
        let mut momenta = LoopMomenta::from_iter(
            loop_momenta
                .iter()
                .map(|momentum| momentum.map_ref(&|value| new_constant(&template, value))),
        );
        let mut measures = Vec::with_capacity(self.factors.len());

        for (axis, ((factor, mass), nu)) in self
            .factors
            .iter()
            .zip(&masses)
            .zip(oriented_chemical_potentials)
            .enumerate()
        {
            if nu == mass {
                return Err(eyre!(
                    "Degenerate Fermi surface on edge {} at sigma*mu = mass; the regular-shell localization does not define the onset",
                    factor.edge_id
                ));
            }
            let momentum = &loop_momenta[axis.into()];
            let norm_squared = momentum.norm_squared();
            if norm_squared <= representative.zero() || !norm_squared.is_finite() {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "Fermi localization",
                    detail: format!(
                        "Fermi surface on edge {} requires a finite nonzero radial sampling momentum",
                        factor.edge_id
                    ),
                }
                .into());
            }
            let mut energy = new_constant(&template, nu);
            if self.derivative_orders[axis] > 0 {
                let mut unit = vec![0; self.factors.len()];
                unit[axis] = 1;
                let index = self
                    .shape
                    .iter()
                    .position(|orders| *orders == unit)
                    .expect("mixed derivative shape contains each active first derivative");
                energy.values[index] = representative.one();
            }
            // Take square roots before multiplying, preserving representable
            // radii even when their squares underflow. The massless radius is
            // exactly E; using it directly also avoids cancelling higher jets.
            let radius = if *mass == mass.zero() {
                energy.clone()
            } else {
                let mass = new_constant(&template, mass);
                (energy.clone() - &mass).sqrt() * (energy.clone() + &mass).sqrt()
            };
            let scale = radius / new_constant(&template, &norm_squared.sqrt());
            if scale.values[0] <= representative.zero()
                || scale.values.iter().any(|value| !value.is_finite())
            {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "Fermi localization",
                    detail: format!("Nonfinite Fermi root jet on edge {}", factor.edge_id),
                }
                .into());
            }
            let density = profile(&scale)?;
            if density.get_shape() != template.get_shape()
                || density.values.len() != self.shape.len()
            {
                return Err(eyre!(
                    "Fermi scale profile returned an incompatible derivative shape"
                ));
            }
            momenta[axis.into()] =
                momentum.map_ref(&|value| new_constant(&template, value) * &scale);
            measures.push((density, energy, scale, norm_squared));
        }

        let value = coefficient(&momenta)?;
        if value.get_shape() != template.get_shape() || value.values.len() != self.shape.len() {
            return Err(eyre!(
                "Fermi coefficient returned an incompatible derivative shape"
            ));
        }
        // E*t*t/|k|^2 equals t^3/g'(t) without forming |k|^4. Keep its factors
        // separate until the coefficient is available: either multiplying the
        // small profile first or the large geometry first can lose a finite
        // fully weighted result. Logs below only choose the arithmetic order.
        let mut factors = Vec::with_capacity(5 * measures.len());
        for (density, energy, scale, norm_squared) in measures {
            for (factor, divide) in [
                (density, false),
                (energy, false),
                (scale.clone(), false),
                (scale, false),
                (new_constant(&template, &norm_squared), true),
            ] {
                if factor.values.iter().any(|value| !value.is_finite()) {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "Fermi localization",
                        detail: "Nonfinite Fermi measure jet".to_owned(),
                    }
                    .into());
                }
                let magnitude = factor
                    .values
                    .iter()
                    .map(|value| value.abs())
                    .max_by(|left, right| left.partial_cmp(right).unwrap())
                    .unwrap();
                let score = if divide {
                    -magnitude.log()
                } else {
                    magnitude.log()
                };
                factors.push((factor, divide, score));
            }
        }
        let index = self
            .shape
            .iter()
            .position(|orders| orders == &self.derivative_orders)
            .expect("mixed derivative shape contains the requested derivative");
        let mut components = [representative.zero(), representative.zero()];
        for (component, imaginary) in components.iter_mut().zip([false, true]) {
            let mut weighted = HyperDual::from_values(
                self.shape.clone(),
                value
                    .values
                    .iter()
                    .map(|value| {
                        if imaginary {
                            value.im.clone()
                        } else {
                            value.re.clone()
                        }
                    })
                    .collect(),
            );
            if weighted.values.iter().any(|value| !value.is_finite()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "Fermi localization",
                    detail: "Nonfinite Fermi coefficient jet".to_owned(),
                }
                .into());
            }
            // Balance each observable component independently. A large imaginary
            // coefficient must not determine the scaling of a small real one.
            let mut remaining = (0..factors.len()).collect::<Vec<_>>();
            remaining
                .sort_by(|left, right| factors[*left].2.partial_cmp(&factors[*right].2).unwrap());
            while !remaining.is_empty() {
                let magnitude = weighted
                    .values
                    .iter()
                    .map(|value| value.abs())
                    .max_by(|left, right| left.partial_cmp(right).unwrap())
                    .unwrap();
                let next = if magnitude >= representative.one() {
                    remaining.remove(0)
                } else {
                    remaining.pop().unwrap()
                };
                let (factor, divide, _) = &factors[next];
                if *divide {
                    for value in &mut weighted.values {
                        *value /= &factor.values[0];
                    }
                } else {
                    weighted *= factor;
                }
                if weighted.values.iter().any(|value| !value.is_finite()) {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "Fermi localization",
                        detail: "Nonfinite localized Fermi contribution".to_owned(),
                    }
                    .into());
                }
            }
            *component = weighted.values[index].clone();
            // Convert Taylor coefficients to ordinary distribution derivatives,
            // avoiding an independently overflowing integer or float factorial.
            for order in &self.derivative_orders {
                for factor in 1..=*order {
                    *component *= -representative.from_usize(factor);
                }
            }
        }
        let [real, imaginary] = components;
        let result = Complex::new(real, imaginary);
        if !result.re.is_finite() || !result.im.is_finite() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "Fermi localization",
                detail: "Nonfinite localized Fermi contribution".to_owned(),
            }
            .into());
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod residue_tests;

#[cfg(test)]
mod subtraction_tests;
