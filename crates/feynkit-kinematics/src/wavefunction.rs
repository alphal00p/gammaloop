//! Numerical external states in GammaLoop's MadGraph phase convention.
//!
//! These are four-dimensional external vectors and chiral-basis spinors. The
//! dimension of internal Lorentz/Dirac algebra is a separate caller choice.
use std::array;
use std::{fmt, str::FromStr};

use numerica::domains::float::Complex;
use thiserror::Error;

use crate::{Axis, FourMomentum, Helicity, KinematicScalar, Rotation, Sign, ThreeMomentum};

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
    #[error("longitudinal polarization requires positive mass squared")]
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

    /// The state under an active spatial rotation of its event.
    ///
    /// Vectors rotate their spatial part. Both chiral blocks of a spinor
    /// transform with the rotation's SU(2) matrix, and adjoints with its
    /// complex conjugate.
    ///
    /// Returns an error if the rotated components are not finite.
    pub fn rotate(&self, rotation: &Rotation<T>) -> Result<Self, WavefunctionError> {
        let mut components = self.components.clone();
        match self.kind {
            WavefunctionKind::Scalar => {}
            WavefunctionKind::Epsilon | WavefunctionKind::EpsilonBar => {
                let part = |part: fn(&Complex<T>) -> T| {
                    rotation.rotate_three(&ThreeMomentum::new(
                        part(&components[1]),
                        part(&components[2]),
                        part(&components[3]),
                    ))
                };
                let (re, im) = (part(|c| c.re.clone()), part(|c| c.im.clone()));
                components[1] = Complex::new(re.px, im.px);
                components[2] = Complex::new(re.py, im.py);
                components[3] = Complex::new(re.pz, im.pz);
            }
            kind => {
                let mut matrix = spin_matrix(rotation, &components[0].re);
                if !matches!(kind, WavefunctionKind::U | WavefunctionKind::V) {
                    matrix = matrix.map(|row| row.map(|entry| entry.conj()));
                }
                for block in components.as_chunks_mut::<2>().0 {
                    let [a, b] = [block[0].clone(), block[1].clone()];
                    block[0] = matrix[0][0].clone() * a.clone() + matrix[0][1].clone() * b.clone();
                    block[1] = matrix[1][0].clone() * a + matrix[1][1].clone() * b;
                }
            }
        }
        Self::from_components(self.kind, components)
    }
}

type SpinMatrix<T> = [[Complex<T>; 2]; 2];

/// The SU(2) matrix `cos(t/2) - i sin(t/2) n.sigma` of a rotation by `t` about `n`.
fn axis_spin_matrix<T: KinematicScalar>(axis: Axis, cos: T, sin: T) -> SpinMatrix<T> {
    let zero = cos.zero();
    let c = |re: &T, im: &T| Complex::new(re.clone(), im.clone());
    let minus_sin = -sin.clone();
    match axis {
        Axis::X => [
            [c(&cos, &zero), c(&zero, &minus_sin)],
            [c(&zero, &minus_sin), c(&cos, &zero)],
        ],
        Axis::Y => [
            [c(&cos, &zero), c(&minus_sin, &zero)],
            [c(&sin, &zero), c(&cos, &zero)],
        ],
        Axis::Z => [
            [c(&cos, &minus_sin), c(&zero, &zero)],
            [c(&zero, &zero), c(&cos, &sin)],
        ],
    }
}

fn spin_matrix_product<T: KinematicScalar>(a: &SpinMatrix<T>, b: &SpinMatrix<T>) -> SpinMatrix<T> {
    array::from_fn(|i| {
        array::from_fn(|j| a[i][0].clone() * b[0][j].clone() + a[i][1].clone() * b[1][j].clone())
    })
}

/// The SU(2) matrix of a rotation, `Rz(gamma) Ry(beta) Rx(alpha)` for Euler
/// angles; `unit` supplies the scalar precision.
fn spin_matrix<T: KinematicScalar>(rotation: &Rotation<T>, unit: &T) -> SpinMatrix<T> {
    let about = |angle: &T, axis| {
        let half = angle.clone() / angle.from_usize(2);
        axis_spin_matrix(axis, half.cos(), half.sin())
    };
    match rotation {
        Rotation::Identity => axis_spin_matrix(Axis::Z, unit.one(), unit.zero()),
        Rotation::QuarterTurn { axis } => {
            let half = unit.from_usize(2).sqrt().inv();
            axis_spin_matrix(*axis, half.clone(), half)
        }
        Rotation::Euler { alpha, beta, gamma } => spin_matrix_product(
            &spin_matrix_product(&about(gamma, Axis::Z), &about(beta, Axis::Y)),
            &about(alpha, Axis::X),
        ),
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
        if (self.spatial.pz < zero && self.spatial.py == zero && self.spatial.px == zero)
            || self.spatial.pz == -p.clone()
        {
            let z = Complex::new(zero.clone(), zero);
            let one = z.one();
            return match helicity {
                Sign::Positive => [z, -one],
                Sign::Negative => [one, z],
            };
        }
        let factor = (p.from_usize(2) * p.clone() * (p.clone() + self.spatial.pz.clone()))
            .sqrt()
            .inv();
        let mut xi = [
            Complex::new(factor.clone() * (p + self.spatial.pz.clone()), zero),
            Complex::new(
                self.spatial.px.clone() * factor.clone(),
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
    /// Longitudinal states require positive mass squared. At zero spatial
    /// momentum the longitudinal state is HELAS's rest-frame `(0,0,0,1)`, the
    /// `|p| -> 0` limit of the moving state for momentum along `+z`.
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
                if mass2 <= zero {
                    return Err(WavefunctionError::UndefinedLongitudinal);
                }
                // Normalize the spatial direction before squaring components.
                // A tiny nonzero momentum still fixes a longitudinal axis.
                let spatial_scale = self
                    .spatial
                    .components()
                    .into_iter()
                    .map(|component| component.norm())
                    .fold(zero.clone(), |largest, component| {
                        if component > largest {
                            component
                        } else {
                            largest
                        }
                    });
                let components = if spatial_scale == zero {
                    [zero.clone(), zero.clone(), zero.clone(), zero.one()]
                } else {
                    let scaled = self
                        .spatial
                        .map_ref(|component| component.clone() / spatial_scale.clone());
                    let scaled_norm = scaled.norm();
                    let mass = mass2.sqrt();
                    let energy_over_mass = self.temporal.value.clone() / mass.clone();
                    [
                        spatial_scale * scaled_norm.clone() / mass,
                        scaled.px / scaled_norm.clone() * energy_over_mass.clone(),
                        scaled.py / scaled_norm.clone() * energy_over_mass.clone(),
                        scaled.pz / scaled_norm * energy_over_mass,
                    ]
                };
                components
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

    /// The external state of this momentum after active rotations of its event.
    ///
    /// `rotations` were applied to the event in order and `self` is the rotated
    /// momentum. The result is the state of the unrotated momentum transformed
    /// by each rotation, so a state at rest keeps its rotated quantization axis.
    /// For a moving state it differs from [`Self::wavefunction`] of the rotated
    /// momentum only by a phase.
    pub fn wavefunction_after_rotations(
        &self,
        kind: WavefunctionKind,
        helicity: Helicity,
        rotations: &[Rotation<T>],
    ) -> Result<Wavefunction<T>, WavefunctionError> {
        let reference = rotations
            .iter()
            .rev()
            .fold(self.clone(), |p, rotation| rotation.inverse_rotate_four(&p));
        rotations.iter().try_fold(
            reference.wavefunction(kind, helicity)?,
            |state, rotation| state.rotate(rotation),
        )
    }
}
