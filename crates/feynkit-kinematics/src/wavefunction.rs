//! Numerical external states in GammaLoop's MadGraph phase convention.
//!
//! These are four-dimensional external vectors and chiral-basis spinors. The
//! dimension of internal Lorentz/Dirac algebra is a separate caller choice.
use std::{fmt, str::FromStr};

use numerica::domains::float::Complex;
use thiserror::Error;

use crate::{FourMomentum, Helicity, KinematicScalar, Sign};

#[cfg(test)]
mod tests;

/// The existing scalar, vector and fermion external-state conventions.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum WavefunctionKind {
    Scalar,
    Epsilon,
    EpsilonBar,
    U,
    UBar,
    V,
    VBar,
}

impl WavefunctionKind {
    pub const fn bar(self) -> Self {
        match self {
            Self::Scalar => Self::Scalar,
            Self::Epsilon => Self::EpsilonBar,
            Self::EpsilonBar => Self::Epsilon,
            Self::U => Self::UBar,
            Self::UBar => Self::U,
            Self::V => Self::VBar,
            Self::VBar => Self::V,
        }
    }

    pub const fn is_spinor(self) -> bool {
        matches!(self, Self::U | Self::UBar | Self::V | Self::VBar)
    }

    pub const fn component_count(self) -> usize {
        if matches!(self, Self::Scalar) { 1 } else { 4 }
    }
}

impl fmt::Display for WavefunctionKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Scalar => "scalar",
            Self::Epsilon => "epsilon",
            Self::EpsilonBar => "epsilon_bar",
            Self::U => "u",
            Self::UBar => "u_bar",
            Self::V => "v",
            Self::VBar => "v_bar",
        })
    }
}

impl FromStr for WavefunctionKind {
    type Err = WavefunctionError;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "scalar" => Ok(Self::Scalar),
            "epsilon" => Ok(Self::Epsilon),
            "epsilon_bar" => Ok(Self::EpsilonBar),
            "u" => Ok(Self::U),
            "u_bar" => Ok(Self::UBar),
            "v" => Ok(Self::V),
            "v_bar" => Ok(Self::VBar),
            _ => Err(WavefunctionError::UnknownKind(value.into())),
        }
    }
}

#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum WavefunctionError {
    #[error("unknown wavefunction kind '{0}'")]
    UnknownKind(String),
    #[error("{kind} requires {expected} components, received {actual}")]
    ComponentCount {
        kind: WavefunctionKind,
        expected: usize,
        actual: usize,
    },
    #[error("external momentum or wavefunction components are not finite")]
    NonFinite,
    #[error("{0} does not admit the requested helicity")]
    Helicity(WavefunctionKind),
    #[error(
        "longitudinal polarization requires positive mass squared and nonzero spatial momentum"
    )]
    UndefinedLongitudinal,
}

/// Numerical components with their vector/spinor adjoint convention.
#[derive(Clone, Debug, PartialEq)]
pub struct Wavefunction<T: KinematicScalar> {
    kind: WavefunctionKind,
    components: Vec<Complex<T>>,
}

impl<T: KinematicScalar> Wavefunction<T> {
    /// Import a scalar factor or already constructed vector/spinor state.
    pub fn from_components(
        kind: WavefunctionKind,
        components: Vec<Complex<T>>,
    ) -> Result<Self, WavefunctionError> {
        if components.len() != kind.component_count() {
            return Err(WavefunctionError::ComponentCount {
                kind,
                expected: kind.component_count(),
                actual: components.len(),
            });
        }
        if components
            .iter()
            .any(|c| !c.re.is_finite() || !c.im.is_finite())
        {
            return Err(WavefunctionError::NonFinite);
        }
        Ok(Self { kind, components })
    }

    pub fn kind(&self) -> WavefunctionKind {
        self.kind
    }
    pub fn components(&self) -> &[Complex<T>] {
        &self.components
    }
    pub fn into_components(self) -> Vec<Complex<T>> {
        self.components
    }

    /// Complex conjugation for vectors/scalars; the chiral Dirac adjoint for spinors.
    pub fn bar(&self) -> Self {
        let mut components = self
            .components
            .iter()
            .map(Complex::conj)
            .collect::<Vec<_>>();
        if self.kind.is_spinor() {
            components.swap(0, 2);
            components.swap(1, 3);
        }
        Self {
            kind: self.kind.bar(),
            components,
        }
    }
}

impl<T: KinematicScalar> FourMomentum<T> {
    /// GammaLoop's two real transverse basis vectors, in `(E,x,y,z)` order.
    ///
    /// The beam-axis branches and phases are inherited from its HELAS/MadGraph
    /// implementation. This basis also fixes the convention at zero momentum.
    pub fn transverse_basis(&self) -> [[T; 4]; 2] {
        let pt = self.pt();
        let p = self.spatial.norm();
        let zero = pt.zero();
        let (a, b) = if pt == zero {
            (
                [zero.clone(), pt.one(), zero.clone(), zero.clone()],
                [
                    zero.clone(),
                    zero.clone(),
                    if self.spatial.pz > zero {
                        pt.one()
                    } else {
                        -pt.one()
                    },
                    zero,
                ],
            )
        } else {
            (
                [
                    zero.clone(),
                    self.spatial.px.clone() * self.spatial.pz.clone() / (pt.clone() * p.clone()),
                    self.spatial.py.clone() * self.spatial.pz.clone() / (pt.clone() * p.clone()),
                    -(pt.clone() / p),
                ],
                [
                    zero.clone(),
                    -(self.spatial.py.clone() / pt.clone()),
                    self.spatial.px.clone() / pt,
                    zero,
                ],
            )
        };
        [a, b]
    }

    /// The two-component helicity spinor used by both `u` and `v`.
    pub fn helicity_spinor(&self, helicity: Sign) -> [Complex<T>; 2] {
        let p = self.spatial.norm();
        let zero = self.temporal.value.zero();
        if self.spatial.pz <= zero && self.spatial.py == zero && self.spatial.px == zero {
            let z = Complex::new(zero.clone(), zero);
            let one = z.one();
            return match helicity {
                Sign::Positive => [z, -one],
                Sign::Negative => [one, z],
            };
        }
        let (upper, factor) = if self.spatial.pz < zero {
            // Rationalize p+pz = pt^2/(p-pz), retaining the off-axis azimuthal phase.
            let pt = self.pt();
            let p_minus_pz = p.clone() - &self.spatial.pz;
            let lower = (p_minus_pz.clone() / (p.from_usize(2) * &p)).sqrt();
            (pt.clone() / p_minus_pz * &lower, lower / pt)
        } else {
            let p_plus_pz = p.clone() + &self.spatial.pz;
            let factor = (p.from_usize(2) * p * &p_plus_pz).sqrt().inv();
            (factor.clone() * p_plus_pz, factor)
        };
        let mut xi = [
            Complex::new(upper, zero),
            Complex::new(
                self.spatial.px.clone() * &factor,
                self.spatial.py.clone() * factor,
            ),
        ];
        if helicity == Sign::Negative {
            xi.swap(0, 1);
            xi[0].re = -xi[0].re.clone();
        }
        xi
    }

    /// GammaLoop's principal square root of `E ± |p|`.
    pub fn spinor_energy_factor(&self, sign: Sign) -> Complex<T> {
        let p = self.spatial.norm();
        let value = self.temporal.value.clone() + if sign == Sign::Positive { p } else { -p };
        if value > value.zero() {
            Complex::new(value.sqrt(), value.zero())
        } else {
            Complex::new(value.zero(), (-value).sqrt())
        }
    }

    /// Construct one fixed external state, without averaging or summing helicities.
    ///
    /// Vectors use signature `+---`. Fermions use GammaLoop's chiral gamma-matrix
    /// basis and MadGraph phases, including the negative-z axis limit. All states
    /// are four-dimensional; keep the internal regulator dimension separately.
    /// Longitudinal states with zero mass or zero spatial momentum are undefined
    /// in the inherited helicity convention and return an error.
    pub fn wavefunction(
        &self,
        kind: WavefunctionKind,
        helicity: Helicity,
    ) -> Result<Wavefunction<T>, WavefunctionError> {
        if [
            &self.temporal.value,
            &self.spatial.px,
            &self.spatial.py,
            &self.spatial.pz,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(WavefunctionError::NonFinite);
        }
        if matches!(
            kind,
            WavefunctionKind::EpsilonBar | WavefunctionKind::UBar | WavefunctionKind::VBar
        ) {
            return Ok(self.wavefunction(kind.bar(), helicity)?.bar());
        }
        let zero = self.temporal.value.zero();
        let components = match kind {
            WavefunctionKind::Scalar => {
                if !helicity.is_zero() {
                    return Err(WavefunctionError::Helicity(kind));
                }
                vec![Complex::new(self.temporal.value.one(), zero)]
            }
            WavefunctionKind::Epsilon if helicity.is_zero() => {
                let mass2 = self.mass_squared();
                let p = self.spatial.norm();
                if mass2 <= zero || p == zero {
                    return Err(WavefunctionError::UndefinedLongitudinal);
                }
                let mass = mass2.sqrt();
                let scale = self.temporal.value.clone() / (mass.clone() * p.clone());
                [
                    p / mass,
                    self.spatial.px.clone() * scale.clone(),
                    self.spatial.py.clone() * scale.clone(),
                    self.spatial.pz.clone() * scale,
                ]
                .into_iter()
                .map(|x| Complex::new(x, zero.clone()))
                .collect()
            }
            WavefunctionKind::Epsilon => {
                let [a, b] = self.transverse_basis();
                let scale = self.temporal.value.from_usize(2).sqrt().inv();
                a.into_iter()
                    .zip(b)
                    .map(|(a, b)| {
                        Complex::new(
                            (if helicity.is_plus() { -a } else { a }) * scale.clone(),
                            -b * scale.clone(),
                        )
                    })
                    .collect()
            }
            WavefunctionKind::U | WavefunctionKind::V => {
                let sign = if helicity.is_plus() {
                    Sign::Positive
                } else if helicity.is_minus() {
                    Sign::Negative
                } else {
                    return Err(WavefunctionError::Helicity(kind));
                };
                let xi = self.helicity_spinor(if kind == WavefunctionKind::U {
                    sign
                } else {
                    -sign
                });
                let (lower, upper) = if kind == WavefunctionKind::U {
                    (
                        self.spinor_energy_factor(-sign),
                        self.spinor_energy_factor(sign),
                    )
                } else {
                    let a = self.spinor_energy_factor(sign);
                    let b = self.spinor_energy_factor(-sign);
                    if sign == Sign::Positive {
                        (-a, b)
                    } else {
                        (a, -b)
                    }
                };
                vec![
                    lower.clone() * xi[0].clone(),
                    lower * xi[1].clone(),
                    upper.clone() * xi[0].clone(),
                    upper * xi[1].clone(),
                ]
            }
            _ => unreachable!("adjoints returned above"),
        };
        Wavefunction::from_components(kind, components)
    }
}

#[cfg(test)]
mod helicity_spinor_tests {
    use numerica::domains::float::Float;

    use super::*;

    fn check_near_negative_axis_spinors<T: KinematicScalar>(offsets: &[T], tolerance: &T) {
        let one = tolerance.one();
        let zero = tolerance.zero();
        let tolerance_squared = tolerance.clone() * tolerance;
        for (x, y) in [(3, 4), (-3, 4), (3, -4), (-3, -4), (0, 5), (5, 0)] {
            for helicity in [Sign::Positive, Sign::Negative] {
                let mut previous: Option<([Complex<T>; 2], T)> = None;
                for offset in offsets {
                    let p = FourMomentum::from_args(
                        one.from_usize(2),
                        offset.clone() * one.from_i64(x),
                        offset.clone() * one.from_i64(y),
                        -one.clone(),
                    );
                    let xi = p.helicity_spinor(helicity);
                    let norm = xi[0].norm_squared() + xi[1].norm_squared();
                    assert!((norm - &one).pow(2) < tolerance_squared, "{p:?}: {xi:?}");

                    // Independent half-angle reference, including the azimuthal phase.
                    let half_angle =
                        (offset.clone() * one.from_usize(5)).atan2(&one) / one.from_usize(2);
                    let small = half_angle.sin();
                    let large = half_angle.cos();
                    let real = one.from_i64(x) / one.from_usize(5) * &large;
                    let imag = one.from_i64(y) / one.from_usize(5) * large;
                    let target = match helicity {
                        Sign::Positive => {
                            [Complex::new(small, zero.clone()), Complex::new(real, imag)]
                        }
                        Sign::Negative => {
                            [Complex::new(-real, imag), Complex::new(small, zero.clone())]
                        }
                    };
                    let eigenvalue = match helicity {
                        Sign::Positive => p.spatial.norm(),
                        Sign::Negative => -p.spatial.norm(),
                    };
                    let sigma_xi = [
                        &xi[0] * &p.spatial.pz
                            + Complex::new(p.spatial.px.clone(), -p.spatial.py.clone()) * &xi[1],
                        Complex::new(p.spatial.px.clone(), p.spatial.py.clone()) * &xi[0]
                            - &xi[1] * &p.spatial.pz,
                    ];
                    for (i, expected) in target.iter().enumerate() {
                        // Relative component checks retain sensitivity to the small component.
                        assert!(
                            (&xi[i] - expected).norm_squared()
                                < tolerance_squared.clone() * expected.norm_squared(),
                            "{p:?}, {helicity:?}: {xi:?} != {target:?}"
                        );
                        assert!(
                            (&sigma_xi[i] - &xi[i] * &eigenvalue).norm_squared()
                                < tolerance_squared.clone()
                                    * xi[i].norm_squared()
                                    * eigenvalue.clone()
                                    * &eigenvalue,
                            "{p:?}, {helicity:?}: sigma.p xi != h |p| xi"
                        );
                    }
                    if let Some((previous_xi, previous_offset)) = previous {
                        let bound = previous_offset * one.from_usize(5);
                        for (current, previous) in xi.iter().zip(previous_xi) {
                            assert!((current - previous).norm_squared() < bound.clone() * &bound);
                        }
                    }
                    previous = Some((xi, offset.clone()));
                }
            }
        }
    }

    #[test]
    fn near_negative_axis_spinors_are_normalized_eigenvectors_with_continuous_phases() {
        let offsets = [1e-3, 1e-6, 1e-7, 1e-8, 1e-10, 1e-20];
        check_near_negative_axis_spinors(&offsets, &2e-14);
    }

    #[test]
    fn rounded_negative_axis_keeps_nonzero_transverse_components() {
        let rounded_axis = FourMomentum::from_args(2., 3e-20, 4e-20, -1.);
        assert_eq!(rounded_axis.spatial.pz, -rounded_axis.spatial.norm());
        check_near_negative_axis_spinors(&[1e-20], &2e-14);
    }

    #[test]
    fn helicity_spinor_phases_are_continuous_between_hemispheres() {
        let equator = FourMomentum::from_args(6., 3., 4., 0.);
        for helicity in [Sign::Positive, Sign::Negative] {
            let target = equator.helicity_spinor(helicity);
            for pz in [-1e-12, 1e-12] {
                let p = FourMomentum::from_args(6., 3., 4., pz);
                for (a, b) in p.helicity_spinor(helicity).into_iter().zip(target) {
                    assert!((a - b).norm_squared() < 1e-24);
                }
            }
        }
    }

    #[test]
    fn near_negative_axis_spinors_preserve_mpfr_precision() {
        let offsets = ["1e-3", "1e-20", "1e-38", "1e-40", "1e-60", "1e-100"]
            .map(|value| Float::parse(value, Some(256)).unwrap());
        let tolerance = Float::parse("1e-60", Some(256)).unwrap();
        check_near_negative_axis_spinors(&offsets, &tolerance);
    }

    #[test]
    fn exact_axis_and_rest_spinor_phases_are_preserved() {
        for pz in [-2., 0.] {
            let p = FourMomentum::from_args(4., 0., 0., pz);
            let zero = Complex::new(0., 0.);
            let one = Complex::new(1., 0.);
            assert_eq!(p.helicity_spinor(Sign::Positive), [zero, -one]);
            assert_eq!(p.helicity_spinor(Sign::Negative), [one, zero]);
        }
    }
}
