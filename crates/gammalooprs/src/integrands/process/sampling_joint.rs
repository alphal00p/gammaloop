//! A compact full-circle chart for two energy sums with one shared energy.
//!
//! The disk certificate concerns the exact represented prepared parameters.
//! A physical binder must additionally preserve the chosen branch/radius across
//! precision rescue and account for uncertainty in its host preparation.

use std::sync::{Arc, Mutex};

use color_eyre::eyre::{Result, WrapErr, eyre};
use rug::{
    Float,
    float::{Round, Special},
};
use symbolica::{
    atom::{Atom, AtomCore},
    prelude::{Real, SingleFloat},
    try_parse,
};

use crate::utils::{F, FloatLike};

use super::{
    sampling_context::{SamplingMapContext, SamplingProposalDecision},
    sampling_evaluator::SamplingExpressionEvaluator,
    sampling_maps::{
        SamplingEvaluationError, SamplingJacobian, SamplingMapComponent, SamplingMapContract,
        SamplingMapEvaluation, SamplingSupport, SurfaceRadialMap,
    },
};

/// Parameters of E0(x)+E1(x+a)-C1 and E0(x)+E2(x+b)-C2.
/// The graph matcher establishes the shared signed routing and masses.
#[derive(Clone, Debug)]
pub struct SharedEnergyJointGeometry<T: FloatLike = f64> {
    pub shifts: [[T; 3]; 2],
    pub masses: [T; 3],
    pub energy_sums: [T; 2],
}

pub type SharedEnergyJointGeometryEvaluator<T> =
    Arc<dyn Fn(&[T]) -> Result<SharedEnergyJointGeometry<T>> + Send + Sync>;

/// One compact conditional chart, or the normalized ordinary chart when its
/// complement-only domain policy declines focusing. Its static support remains
/// restricted, so the resolved catalogue must also retain a full-support map.
pub struct SharedEnergyJointMap<T: FloatLike = f64> {
    geometry: SharedEnergyJointGeometryEvaluator<T>,
    context_dimension: usize,
    max_radius: T,
    normal_scale: T,
    fallback: SurfaceRadialMap<T>,
    program: Mutex<SamplingExpressionEvaluator>,
}

impl<T: FloatLike> Clone for SharedEnergyJointMap<T> {
    fn clone(&self) -> Self {
        Self {
            geometry: self.geometry.clone(),
            context_dimension: self.context_dimension,
            max_radius: self.max_radius.clone(),
            normal_scale: self.normal_scale.clone(),
            fallback: self.fallback.clone(),
            program: Mutex::new(
                self.program
                    .lock()
                    .expect("joint evaluator poisoned")
                    .clone(),
            ),
        }
    }
}

impl<T: FloatLike> std::fmt::Debug for SharedEnergyJointMap<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SharedEnergyJointMap")
            .field("context_dimension", &self.context_dimension)
            .field("max_radius", &self.max_radius)
            .field("normal_scale", &self.normal_scale)
            .finish_non_exhaustive()
    }
}

// Each endpoint of the conservative box range is itself enclosed. This keeps
// interval dependency/overestimation distinct from MPFR rounding uncertainty:
// only a certified sign of the mathematical lower endpoint changes the disk.
// This arithmetic is deliberately private to the nine rational predicates.
#[derive(Clone, Debug)]
struct JointRange {
    lower: [Float; 2],
    upper: [Float; 2],
}

impl JointRange {
    const PRECISION: u32 = 2048;
    const PRECISIONS: [u32; 5] = [128, 256, 512, 1024, Self::PRECISION];

    fn native<T: FloatLike>(value: &T, precision: u32) -> Self {
        let (lo, hi) = value.mpfr_enclosure(precision);
        Self {
            lower: [lo.clone(), hi.clone()],
            upper: [lo, hi],
        }
    }

    fn integer(value: i32, precision: u32) -> Self {
        let value = Float::with_val(precision, value);
        Self {
            lower: [value.clone(), value.clone()],
            upper: [value.clone(), value],
        }
    }

    fn bounds(lower: Self, upper: Self) -> Self {
        Self {
            lower: lower.lower,
            upper: upper.upper,
        }
    }

    fn add_bound(a: &[Float; 2], b: &[Float; 2]) -> [Float; 2] {
        [
            Float::with_val_round(a[0].prec(), &a[0] + &b[0], Round::Down).0,
            Float::with_val_round(a[0].prec(), &a[1] + &b[1], Round::Up).0,
        ]
    }

    fn mul_bound(a: &[Float; 2], b: &[Float; 2]) -> [Float; 2] {
        let products = a.iter().flat_map(|x| b.iter().map(move |y| (x, y)));
        let mut lower = Float::with_val(a[0].prec(), Special::Infinity);
        let mut upper = -lower.clone();
        for (x, y) in products {
            lower = lower.min(&Float::with_val_round(a[0].prec(), x * y, Round::Down).0);
            upper = upper.max(&Float::with_val_round(a[0].prec(), x * y, Round::Up).0);
        }
        [lower, upper]
    }

    fn add(&self, rhs: &Self) -> Self {
        Self {
            lower: Self::add_bound(&self.lower, &rhs.lower),
            upper: Self::add_bound(&self.upper, &rhs.upper),
        }
    }

    fn neg(&self) -> Self {
        Self {
            lower: [-self.upper[1].clone(), -self.upper[0].clone()],
            upper: [-self.lower[1].clone(), -self.lower[0].clone()],
        }
    }

    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    fn mul(&self, rhs: &Self) -> Self {
        let products = [&self.lower, &self.upper].into_iter().flat_map(|a| {
            [&rhs.lower, &rhs.upper]
                .into_iter()
                .map(move |b| Self::mul_bound(a, b))
        });
        let infinity = Float::with_val(self.lower[0].prec(), Special::Infinity);
        let mut lower = [infinity.clone(), infinity.clone()];
        let mut upper = [-infinity.clone(), -infinity];
        for product in products {
            for index in 0..2 {
                lower[index] = lower[index].clone().min(&product[index]);
                upper[index] = upper[index].clone().max(&product[index]);
            }
        }
        Self { lower, upper }
    }

    fn square(&self) -> Self {
        // The dependency-aware lower square avoids artificially negative
        // lower bounds when a normal coordinate crosses zero.
        let zero = Float::with_val(self.lower[0].prec(), 0);
        let a = Self::mul_bound(&self.lower, &self.lower);
        let b = Self::mul_bound(&self.upper, &self.upper);
        let upper = [a[0].clone().max(&b[0]), a[1].clone().max(&b[1])];
        let lower = if self.lower[1] <= 0 && self.upper[0] >= 0 {
            [zero.clone(), zero]
        } else if self.lower[0] > 0 || self.upper[1] < 0 {
            [a[0].clone().min(&b[0]), a[1].clone().min(&b[1])]
        } else {
            // Only arithmetic uncertainty remains about crossing zero.
            [zero, a[1].clone().min(&b[1])]
        };
        Self { lower, upper }
    }

    fn reciprocal_positive(&self) -> Result<Self> {
        if self.lower[0] <= 0 {
            return Err(SamplingEvaluationError::UncertainGeometry {
                detail: "joint certificate denominator is not certified positive".into(),
            }
            .into());
        }
        let reciprocal = |bound: &[Float; 2]| {
            [
                Float::with_val_round(self.lower[0].prec(), 1 / &bound[1], Round::Down).0,
                Float::with_val_round(self.lower[0].prec(), 1 / &bound[0], Round::Up).0,
            ]
        };
        Ok(Self {
            lower: reciprocal(&self.upper),
            upper: reciprocal(&self.lower),
        })
    }

    /// Sign of the mathematical conservative lower endpoint. A negative
    /// endpoint declines a box; it does not prove physical absence.
    fn lower_positive(&self) -> Result<bool> {
        if self.lower[0] > 0 {
            return Ok(true);
        }
        if self.lower[1] <= 0 {
            return Ok(false);
        }
        Err(SamplingEvaluationError::UncertainGeometry {
            detail: "joint certificate lower-endpoint sign is unresolved by directed arithmetic"
                .into(),
        }
        .into())
    }

    fn sqrt_point(&self) -> Result<Self> {
        if self.lower[0] < 0 {
            return Err(SamplingEvaluationError::UncertainGeometry {
                detail: "joint support energy square is not certified nonnegative".into(),
            }
            .into());
        }
        let mut lo = self.lower[0].clone();
        let mut hi = self.upper[1].clone();
        lo.sqrt_round(Round::Down);
        hi.sqrt_round(Round::Up);
        Ok(Self {
            lower: [lo.clone(), hi.clone()],
            upper: [lo, hi],
        })
    }
}

impl<T: FloatLike> SharedEnergyJointGeometry<T> {
    fn validate(&self) -> Result<()> {
        if self
            .shifts
            .iter()
            .flatten()
            .chain(&self.masses)
            .chain(&self.energy_sums)
            .any(|x| !x.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint geometry",
                detail: "non-finite prepared parameter".into(),
            }
            .into());
        }
        if self.masses.iter().any(|m| m < &m.zero()) {
            return Err(eyre!("joint energy masses must be nonnegative"));
        }
        Ok(())
    }

    /// Prepare radius-independent directed ranges for one native geometry.
    /// Only this call's radius trials share them; each new conditional point
    /// prepares its own ranges. Keep the original arithmetic association and
    /// lower-endpoint predicate order when evaluating each candidate disk.
    fn prepare_certificate(
        &self,
        scale: &T,
        precision: u32,
    ) -> Result<impl Fn(&T) -> Result<bool> + use<T>> {
        let a = self.shifts[0]
            .each_ref()
            .map(|value| JointRange::native(value, precision));
        let b = self.shifts[1]
            .each_ref()
            .map(|value| JointRange::native(value, precision));
        let m = self
            .masses
            .each_ref()
            .map(|value| JointRange::native(value, precision));
        let dot = |a: &[JointRange; 3], b: &[JointRange; 3]| {
            a.iter()
                .zip(b)
                .fold(JointRange::integer(0, precision), |sum, (a, b)| {
                    sum.add(&a.mul(b))
                })
        };
        let aa = dot(&a, &a);
        let ab = dot(&a, &b);
        let bb = dot(&b, &b);
        let gram = aa.mul(&bb).sub(&ab.square());
        // A nonpositive Gram declines every candidate before checking any
        // denominator. The caller still performs all native radius halvings
        // and their underflow checks before recording the ordinary decision.
        let fixed = if gram.lower_positive()? {
            let inverse_scale = JointRange::native(scale, precision).reciprocal_positive()?;
            let energy_sums = self
                .energy_sums
                .each_ref()
                .map(|value| JointRange::native(value, precision));
            let inverse_two = JointRange::integer(2, precision).reciprocal_positive()?;
            let mass_squares = m.each_ref().map(JointRange::square);
            let inverse_gram = gram.reciprocal_positive()?;
            Some((
                energy_sums,
                mass_squares,
                [inverse_scale, inverse_two, inverse_gram],
                JointRange::integer(1, precision),
            ))
        } else {
            None
        };
        Ok(move |radius: &T| {
            let Some((energy_sums, mass_squares, [inverse_scale, inverse_two, inverse_gram], one)) =
                &fixed
            else {
                return Ok(false);
            };
            let r = JointRange::native(radius, precision);
            let z = r.mul(inverse_scale);
            let s = energy_sums[0].add(&JointRange::bounds(r.neg(), r));
            let t = energy_sums[1].add(&JointRange::bounds(z.neg(), z));
            let da = s
                .square()
                .add(&mass_squares[0])
                .sub(&mass_squares[1])
                .sub(&aa)
                .mul(inverse_two);
            let db = t
                .square()
                .add(&mass_squares[0])
                .sub(&mass_squares[2])
                .sub(&bb)
                .mul(inverse_two);
            let d = [
                bb.mul(&da).sub(&ab.mul(&db)).mul(inverse_gram),
                aa.mul(&db).sub(&ab.mul(&da)).mul(inverse_gram),
            ];
            let e = [
                ab.mul(&t).sub(&bb.mul(&s)).mul(inverse_gram),
                ab.mul(&s).sub(&aa.mul(&t)).mul(inverse_gram),
            ];
            let plane_dot = |x: &[JointRange; 2], y: &[JointRange; 2]| {
                aa.mul(&x[0])
                    .mul(&y[0])
                    .add(&ab.mul(&x[0].mul(&y[1]).add(&x[1].mul(&y[0]))))
                    .add(&bb.mul(&x[1]).mul(&y[1]))
            };
            let k = plane_dot(&e, &e).sub(one);
            let b = plane_dot(&d, &e);
            let c = mass_squares[0].add(&plane_dot(&d, &d));
            let discriminant = b.square().sub(&k.mul(&c));
            let limits = [
                b.neg().sub(&k.mul(&m[0])),
                k.mul(&s.sub(&m[1])).add(&b),
                k.mul(&t.sub(&m[2])).add(&b),
            ];
            for predicate in [k, discriminant.clone()].into_iter().chain(
                limits
                    .into_iter()
                    .flat_map(|l| [l.clone(), l.square().sub(&discriminant)]),
            ) {
                if !predicate.lower_positive()? {
                    return Ok(false);
                }
            }
            Ok(true)
        })
    }

    fn residuals(&self, point: &[T]) -> ([F<T>; 2], [F<T>; 3]) {
        let zero = F(self.masses[0].zero());
        let energies = std::array::from_fn(|index| {
            let square = point.iter().enumerate().fold(
                F(self.masses[index].clone()).square(),
                |sum, (component, x)| {
                    let shifted = F(x.clone())
                        + if index == 0 {
                            zero.clone()
                        } else {
                            F(self.shifts[index - 1][component].clone())
                        };
                    sum + shifted.square()
                },
            );
            square.sqrt()
        });
        (
            [
                &energies[0] + &energies[1] - F(self.energy_sums[0].clone()),
                &energies[0] + &energies[2] - F(self.energy_sums[1].clone()),
            ],
            energies,
        )
    }

    fn support_relation(&self, point: &[T], radius: &T, scale: &T) -> Result<std::cmp::Ordering> {
        // Refine only arithmetic at the identical point and disk. An unresolved
        // boundary cannot decline support or change the native source law.
        let evaluate = |precision| -> Result<Option<std::cmp::Ordering>> {
            let energies = (0..3)
                .map(|index| {
                    let square = point.iter().enumerate().fold(
                        JointRange::native(&self.masses[index], precision).square(),
                        |sum, (component, x)| {
                            let x = JointRange::native(x, precision);
                            let shifted = if index == 0 {
                                x
                            } else {
                                x.add(&JointRange::native(
                                    &self.shifts[index - 1][component],
                                    precision,
                                ))
                            };
                            sum.add(&shifted.square())
                        },
                    );
                    square.sqrt_point()
                })
                .collect::<Result<Vec<_>>>()?;
            let h = energies[0]
                .add(&energies[1])
                .sub(&JointRange::native(&self.energy_sums[0], precision));
            let z = energies[0]
                .add(&energies[2])
                .sub(&JointRange::native(&self.energy_sums[1], precision));
            let difference = h
                .square()
                .add(&z.mul(&JointRange::native(scale, precision)).square())
                .sub(&JointRange::native(radius, precision).square());
            if difference.lower[0] > 0 {
                return Ok(Some(std::cmp::Ordering::Greater));
            }
            if difference.upper[1] < 0 {
                return Ok(Some(std::cmp::Ordering::Less));
            }
            Ok(None)
        };
        for precision in JointRange::PRECISIONS {
            match evaluate(precision) {
                Ok(Some(relation)) => return Ok(relation),
                Ok(None) => {}
                Err(error)
                    if precision != JointRange::PRECISION
                        && matches!(
                            error.downcast_ref::<SamplingEvaluationError>(),
                            Some(SamplingEvaluationError::UncertainGeometry { .. })
                        ) => {}
                Err(error) => return Err(error),
            }
        }
        Err(SamplingEvaluationError::UncertainGeometry {
            detail: "joint inverse lies on an uncertified normal-disk boundary".into(),
        }
        .into())
    }
}

impl<T: FloatLike> SharedEnergyJointMap<T> {
    pub fn new(
        geometry: SharedEnergyJointGeometryEvaluator<T>,
        context_dimension: usize,
        max_radius: T,
        normal_scale: T,
        beta: f64,
        program: SamplingExpressionEvaluator,
    ) -> Result<Self> {
        if [&max_radius, &normal_scale]
            .iter()
            .any(|x| !x.is_finite() || **x <= x.zero())
        {
            return Err(eyre!(
                "joint normal radius and scale must be positive and finite"
            ));
        }
        if program.parameter_count() != 17
            || program.output_count() != 3
            || program.derivative_parameters() != [0, 1, 2]
        {
            return Err(eyre!(
                "joint chart requires its three-output, seventeen-parameter eager program with active columns [0, 1, 2]"
            ));
        }
        let fallback = SurfaceRadialMap::absent(3, vec![max_radius.zero(); 3], beta, 1.0)?;
        Ok(Self {
            geometry,
            context_dimension,
            max_radius,
            normal_scale,
            fallback,
            program: Mutex::new(program),
        })
    }

    /// Compile the neutral map once; native bindings and worker clones reuse
    /// these programs. Prepared inputs occupy inactive derivative columns.
    pub fn compile_program() -> Result<SamplingExpressionEvaluator> {
        Self::compile_program_with_derivatives(&[0, 1, 2])
    }

    fn compile_program_with_derivatives(
        derivative_parameters: &[usize],
    ) -> Result<SamplingExpressionEvaluator> {
        let names = [
            "jr", "jt", "jp", "rho", "alpha", "tau", "ax", "ay", "az", "bx", "by", "bz", "m0",
            "m1", "m2", "c1", "c2",
        ];
        let p = names
            .iter()
            .map(|name| try_parse!(*name).map_err(|e| eyre!(e)))
            .collect::<Result<Vec<_>>>()?;
        let a = &p[6..9];
        let b = &p[9..12];
        let dot = |a: &[Atom], b: &[Atom]| {
            a.iter()
                .zip(b)
                .fold(Atom::num(0), |sum, (a, b)| sum + a * b)
        };
        let anorm = dot(a, a).sqrt();
        let ea = a.iter().map(|x| x / &anorm).collect::<Vec<_>>();
        let along = dot(b, &ea);
        let bp = b
            .iter()
            .zip(&ea)
            .map(|(b, ea)| b - &along * ea)
            .collect::<Vec<_>>();
        let bnorm = dot(&bp, &bp).sqrt();
        let eb = bp.iter().map(|x| x / &bnorm).collect::<Vec<_>>();
        let normal = [
            &ea[1] * &eb[2] - &ea[2] * &eb[1],
            &ea[2] * &eb[0] - &ea[0] * &eb[2],
            &ea[0] * &eb[1] - &ea[1] * &eb[0],
        ];
        let r = &p[3] * &p[0];
        let theta = &p[5] * &p[1];
        let phi = &p[5] * &p[2];
        let s = &p[15] + &r * theta.cos();
        let t = &p[16] + &r * theta.sin() / &p[4];
        let d1 = (&s * &s + &p[12] * &p[12] - &p[13] * &p[13] - &anorm * &anorm)
            / (Atom::num(2) * &anorm);
        let e1 = -&s / &anorm;
        let d2 = ((&t * &t + &p[12] * &p[12] - &p[14] * &p[14] - dot(b, b)) / Atom::num(2)
            - &along * &d1)
            / &bnorm;
        let e2 = (-&t - &along * &e1) / &bnorm;
        let k = &e1 * &e1 + &e2 * &e2 - Atom::num(1);
        let big_b = &d1 * &e1 + &d2 * &e2;
        let c = &p[12] * &p[12] + &d1 * &d1 + &d2 * &d2;
        let delta = (&big_b * &big_b - &k * c).sqrt() / &k;
        let u = -&big_b / &k + &delta * phi.cos();
        let w = k.sqrt() * delta * phi.sin();
        let outputs =
            (0..3).map(|i| (&d1 + &e1 * &u) * &ea[i] + (&d2 + &e2 * &u) * &eb[i] + &w * &normal[i]);
        SamplingExpressionEvaluator::new(outputs, p.clone(), derivative_parameters)
    }

    fn prepare(
        &self,
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<(SharedEnergyJointGeometry<T>, Option<T>)> {
        if context.previous.len() != self.context_dimension {
            return Err(eyre!(
                "joint chart context has dimension {}, expected {}",
                context.previous.len(),
                self.context_dimension
            ));
        }
        if context.previous.iter().any(|x| !x.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint context",
                detail: "non-finite previous coordinates".into(),
            }
            .into());
        }
        let geometry = (self.geometry)(context.previous)?;
        geometry.validate()?;
        let mut radius = F(self.max_radius.clone());
        let two = radius.from_usize(2);
        let mut compact = None;
        {
            let mut precisions = JointRange::PRECISIONS.into_iter().peekable();
            let mut certificate = geometry.prepare_certificate(
                &self.normal_scale,
                *precisions
                    .peek()
                    .expect("joint precision ladder is nonempty"),
            );
            // A deterministic complement-only policy, not a claim of absence.
            // Refine arithmetic ambiguity at this same native radius before
            // halving; exhausted geometric boxes use the normalized ordinary map.
            // Fixed terms live only in this preparation. Once refinement is
            // needed, later radius trials reuse that more precise certificate.
            for dyadic_exponent in 0..32 {
                let admitted = loop {
                    let precision = *precisions.peek().expect("joint precision ceiling retained");
                    match certificate {
                        Ok(prepared) => {
                            let result = prepared(&radius.0);
                            certificate = Ok(prepared);
                            match result {
                                Err(error)
                                    if precision != JointRange::PRECISION
                                        && matches!(
                                            error.downcast_ref::<SamplingEvaluationError>(),
                                            Some(SamplingEvaluationError::UncertainGeometry { .. })
                                        ) => {}
                                result => break result?,
                            }
                        }
                        Err(error)
                            if precision != JointRange::PRECISION
                                && matches!(
                                    error.downcast_ref::<SamplingEvaluationError>(),
                                    Some(SamplingEvaluationError::UncertainGeometry { .. })
                                ) => {}
                        Err(error) => return Err(error),
                    }
                    precisions.next();
                    certificate = geometry.prepare_certificate(
                        &self.normal_scale,
                        *precisions.peek().expect("joint precision ceiling retained"),
                    );
                };
                if admitted {
                    compact = Some((dyadic_exponent, radius.0));
                    break;
                }
                radius /= &two;
                if !radius.is_finite() || radius <= radius.zero() {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "joint disk",
                        detail: "dyadic radius underflow".into(),
                    }
                    .into());
                }
            }
        }
        let decision = compact.as_ref().map_or(
            SamplingProposalDecision::Ordinary,
            |(dyadic_exponent, _)| SamplingProposalDecision::Compact {
                dyadic_exponent: *dyadic_exponent,
            },
        );
        // Revalidate the complete native choice, not merely whether a stored
        // smaller disk would also be valid. Different valid radii are different laws.
        context.validate_decision(decision).wrap_err_with(|| format!(
            "joint native decision {decision:?}, precision {:?}, max radius {:?}, normal scale {:?}, represented geometry {geometry:?}",
            T::sampling_precision(), self.max_radius, self.normal_scale,
        ))?;
        Ok((geometry, compact.map(|(_, radius)| radius)))
    }

    fn parameters(
        &self,
        geometry: &SharedEnergyJointGeometry<T>,
        radius: &T,
        coordinates: &[T],
    ) -> Vec<T> {
        coordinates
            .iter()
            .cloned()
            .chain([radius.clone(), self.normal_scale.clone(), radius.TAU()])
            .chain(geometry.shifts.iter().flatten().cloned())
            .chain(geometry.masses.iter().cloned())
            .chain(geometry.energy_sums.iter().cloned())
            .collect()
    }
}

// Native inverse geometry. The production forward determinant comes from the
// Symbolica map above; the inverse density uses the original supplied point's
// energies and residual radius, rather than a reconstructed forward point.
struct JointCircle<T: FloatLike> {
    axes: [[F<T>; 3]; 3],
    k: F<T>,
    center: F<T>,
    delta: F<T>,
    cross_norm: F<T>,
}

impl<T: FloatLike> JointCircle<T> {
    fn new(geometry: &SharedEnergyJointGeometry<T>, residuals: &[F<T>; 2]) -> Result<Self> {
        let a = geometry.shifts[0].each_ref().map(|x| F(x.clone()));
        let b = geometry.shifts[1].each_ref().map(|x| F(x.clone()));
        let zero = a[0].zero();
        let dot = |a: &[F<T>; 3], b: &[F<T>; 3]| {
            a.iter()
                .zip(b)
                .fold(zero.clone(), |sum, (a, b)| sum + a * b)
        };
        let anorm = dot(&a, &a).sqrt();
        let ea = a.each_ref().map(|x| x / &anorm);
        let along = dot(&b, &ea);
        let bp = std::array::from_fn(|i| &b[i] - &along * &ea[i]);
        let bnorm = dot(&bp, &bp).sqrt();
        let eb = bp.each_ref().map(|x| x / &bnorm);
        let n = [
            &ea[1] * &eb[2] - &ea[2] * &eb[1],
            &ea[2] * &eb[0] - &ea[0] * &eb[2],
            &ea[0] * &eb[1] - &ea[1] * &eb[0],
        ];
        let s = F(geometry.energy_sums[0].clone()) + &residuals[0];
        let t = F(geometry.energy_sums[1].clone()) + &residuals[1];
        let m = geometry.masses.each_ref().map(|x| F(x.clone()).square());
        let two = s.from_usize(2);
        let d1 = (s.square() + &m[0] - &m[1] - anorm.square()) / (&two * &anorm);
        let e1 = -&s / &anorm;
        let d2 = ((t.square() + &m[0] - &m[2] - dot(&b, &b)) / &two - &along * &d1) / &bnorm;
        let e2 = (-&t - &along * &e1) / &bnorm;
        let k = e1.square() + e2.square() - s.one();
        let big_b = &d1 * &e1 + &d2 * &e2;
        let discriminant = big_b.square() - &k * (&m[0] + d1.square() + d2.square());
        let center = -big_b / &k;
        let delta = discriminant.sqrt() / &k;
        let cross_norm = anorm * bnorm;
        if [&k, &center, &delta, &cross_norm]
            .into_iter()
            .any(|x| !x.is_finite())
            || k <= zero
            || delta <= zero
            || cross_norm <= zero
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint circle",
                detail: "certified circle is not representable in the active native precision"
                    .into(),
            }
            .into());
        }
        Ok(Self {
            axes: [ea, eb, n],
            k,
            center,
            delta,
            cross_norm,
        })
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for SharedEnergyJointMap<T> {
    fn dimensions(&self) -> usize {
        3
    }
    fn output_dimensions(&self) -> usize {
        3
    }
    fn name(&self) -> &'static str {
        "shared_energy_joint"
    }
    fn contract(&self) -> SamplingMapContract {
        SamplingMapContract {
            support: SamplingSupport::Restricted,
            requires_context: self.context_dimension > 0,
            requires_proposal_policy: true,
            jacobian: SamplingJacobian::ExactForward,
        }
    }

    fn forward(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        if coordinates.len() != 3
            || coordinates.iter().enumerate().any(|(i, x)| {
                !x.is_finite() || x < &x.zero() || x >= &x.one() || (i == 0 && x == &x.zero())
            })
        {
            return Err(eyre!(
                "joint chart requires radial coordinate in (0,1) and half-open angular coordinates in [0,1)"
            ));
        }
        let (geometry, radius) = self.prepare(context)?;
        let Some(radius) = radius else {
            // The full-circle seams are regular in the focused chart. The
            // ordinary polar fallback has a different, measure-zero angular
            // boundary there; preserve that law rather than remapping it.
            if coordinates[1..].iter().any(|x| x == &x.zero()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "joint ordinary fallback angular boundary",
                    detail: "exact angular zero lies on the ordinary chart boundary".into(),
                }
                .into());
            }
            let mut output = SamplingMapComponent::forward(
                &self.fallback,
                coordinates,
                &mut context.reborrow(&[], None),
            )?;
            output
                .diagnostics
                .push("joint: normalized ordinary fallback; disk policy declined".into());
            return Ok(output);
        };
        let output = self
            .program
            .lock()
            .map_err(|_| eyre!("joint evaluator poisoned"))?
            .evaluate_with_real_jacobian(
                &self.parameters(&geometry, &radius, coordinates),
                Some(&[0, 1, 2]),
            )?;
        let jacobian = F(output.absolute_determinant());
        let (residuals, _) = geometry.residuals(&output.values);
        let r = F(radius.clone()) * F(coordinates[0].clone());
        let angle = F(radius.TAU()) * F(coordinates[1].clone());
        let residual = (&residuals[0] - &r * angle.cos())
            .abs()
            .max((residuals[1].clone() - r * angle.sin() / F(self.normal_scale.clone())).abs());
        SamplingMapEvaluation {
            coordinates: coordinates.to_vec(),
            point: output.values,
            inverse_jacobian: jacobian.inv().0,
            jacobian: jacobian.0,
            residual: residual.0,
            support: SamplingSupport::Restricted,
            diagnostics: vec![format!("joint: certified normal disk rho={radius}")],
        }
        .validate("joint forward")
    }

    fn inverse(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        self.inverse_evaluation(point, context, true)
    }

    fn inverse_density(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<T>> {
        self.inverse_evaluation(point, context, false)
            .map(|evaluation| evaluation.map(|evaluation| evaluation.inverse_jacobian))
    }
}

impl<T: FloatLike> SharedEnergyJointMap<T> {
    fn inverse_evaluation(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
        reconstruct: bool,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        if point.len() != 3 {
            return Err(eyre!("joint chart point must have three components"));
        }
        if point.iter().any(|x| !x.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint inverse point",
                detail: "non-finite supplied point".into(),
            }
            .into());
        }
        let (geometry, radius) = self.prepare(context)?;
        let Some(radius) = radius else {
            let mut output = SamplingMapComponent::inverse(
                &self.fallback,
                point,
                &mut context.reborrow(&[], None),
            )?;
            if let Some(output) = &mut output {
                output
                    .diagnostics
                    .push("joint: normalized ordinary fallback; disk policy declined".into());
            }
            return Ok(output);
        };
        if geometry.support_relation(point, &radius, &self.normal_scale)?
            == std::cmp::Ordering::Greater
        {
            return Ok(None);
        }
        let (residuals, energies) = geometry.residuals(point);
        let scaled_z = F(self.normal_scale.clone()) * &residuals[1];
        let r = (residuals[0].square() + scaled_z.square()).sqrt();
        if !r.is_finite() || r <= r.zero() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint inverse normal radius",
                detail: "R=0 is singular; a rounded zero requires native retry".into(),
            }
            .into());
        }
        let circle = JointCircle::new(&geometry, &residuals)?;
        let w = point
            .iter()
            .zip(&circle.axes[2])
            .fold(r.zero(), |sum, (x, n)| sum + F(x.clone()) * n);
        let tau = r.TAU();
        let wrap = |angle: F<T>| {
            if angle < angle.zero() {
                angle + &tau
            } else {
                angle
            }
        };
        let theta = wrap(scaled_z.atan2(&residuals[0]));
        let phi = wrap(
            (w / (circle.k.sqrt() * &circle.delta))
                .atan2(&((&energies[0] - &circle.center) / &circle.delta)),
        );
        let coordinates = vec![
            (r.clone() / F(radius.clone())).0,
            (theta / &tau).0,
            (phi / &tau).0,
        ];
        if coordinates[0] <= radius.zero()
            || coordinates
                .iter()
                .any(|x| !x.is_finite() || x < &x.zero() || x >= &x.one())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint inverse coordinates",
                detail: "derived radial or half-open angular coordinate rounded outside its chart"
                    .into(),
            }
            .into());
        }
        let jacobian =
            tau.square() * F(radius.clone()) * &r * &energies[0] * &energies[1] * &energies[2]
                / (F(self.normal_scale.clone()) * &circle.cross_norm * circle.k.sqrt());
        // A density-only query discards this private diagnostic field; full
        // inverses always fill it from the independent forward reconstruction.
        let mut residual = r.zero();
        if reconstruct {
            // Reconstruction is diagnostic only: the density above is evaluated
            // at the original supplied point, including its original normal radius.
            let reconstructed = self
                .program
                .lock()
                .map_err(|_| eyre!("joint evaluator poisoned"))?
                .evaluate(&self.parameters(&geometry, &radius, &coordinates))?;
            for (actual, expected) in reconstructed.iter().zip(point) {
                if !actual.re.is_finite() || actual.im != actual.im.zero() {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "joint inverse reconstruction",
                        detail: "recovered cube does not produce finite real Cartesian coordinates"
                            .into(),
                    }
                    .into());
                }
                residual = residual.max((&actual.re - F(expected.clone())).abs());
            }
        }
        let evaluation = SamplingMapEvaluation {
            coordinates,
            point: point.to_vec(),
            inverse_jacobian: jacobian.inv().0,
            jacobian: jacobian.0,
            residual: residual.0,
            support: SamplingSupport::Restricted,
            diagnostics: vec![format!(
                "joint: certified normal disk rho={radius}; density at supplied point"
            )],
        }
        .validate("joint inverse")?;
        Ok(Some(evaluation))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{ArbPrec, QuadFloat};

    fn geometry<T: FloatLike>() -> SharedEnergyJointGeometry<T> {
        let zero = F::<T>::default();
        let n = |value| zero.from_usize(value).0;
        SharedEnergyJointGeometry {
            shifts: [[n(3), n(0), n(0)], [n(1), n(4), n(0)]],
            masses: [n(2), n(3), n(0)],
            energy_sums: [n(10), n(11)],
        }
    }

    fn map<T: FloatLike>(
        geometry: SharedEnergyJointGeometry<T>,
        program: SamplingExpressionEvaluator,
    ) -> SharedEnergyJointMap<T> {
        SharedEnergyJointMap::new(
            Arc::new(move |_| Ok(geometry.clone())),
            0,
            (F::<T>::default().one() / F::<T>::default().from_usize(8)).0,
            (F::<T>::default().from_usize(3) / F::<T>::default().from_usize(2)).0,
            4.0,
            program,
        )
        .unwrap()
    }

    #[test]
    fn joint_sampling_certificate_encloses_unequal_mass_nonorthogonal_geometry() {
        let geometry = geometry::<f64>();
        // Independently reduced exact-rational oracle at h=z=0:
        // Gram=144,k=1985/144,D=7829/36, all six energy-sign bounds >0.
        let certificate = geometry.prepare_certificate(&1.5, 2048).unwrap();
        assert!(certificate(&0.0).unwrap());
        assert!(certificate(&0.0001).unwrap());
        assert!(!certificate(&100.0).unwrap());
        let mut collinear = geometry.clone();
        collinear.shifts[1] = [6.0, 0.0, 0.0];
        assert!(!collinear.prepare_certificate(&1.5, 2048).unwrap()(&0.01).unwrap());

        // The outer endpoint enclosure separates a failing conservative box
        // from an arithmetic sign that cannot be decided at this precision.
        let uncertain = JointRange {
            lower: [Float::with_val(2048, -1), Float::with_val(2048, 1)],
            upper: [Float::with_val(2048, 2), Float::with_val(2048, 2)],
        };
        assert!(
            uncertain
                .lower_positive()
                .unwrap_err()
                .downcast_ref::<SamplingEvaluationError>()
                .is_some()
        );
        let negative =
            JointRange::bounds(JointRange::integer(-1, 2048), JointRange::integer(1, 2048));
        assert!(!negative.lower_positive().unwrap());
    }

    #[test]
    fn joint_sampling_prepared_certificate_reuses_only_fixed_geometry() {
        fn check<T: FloatLike>() {
            let geometry = geometry::<T>();
            let one = F::<T>::default().one();
            let scale = one.from_usize(3) / one.from_usize(2);
            let certificate = geometry.prepare_certificate(&scale.0, 2048).unwrap();
            let mut accepted = 0;
            let mut declined = 0;
            // Revisit large and small radii in both orders, including the
            // admitted 1/10000 disk and declined radius 100 from the rational
            // fixture above. A prepared closure must agree with fresh
            // preparation and cannot retain a previous candidate's ranges or
            // outcome.
            for (numerator, denominator) in [
                (100, 1),
                (1, 1),
                (1, 10000),
                (0, 1),
                (1, 10000),
                (1, 1),
                (100, 1),
            ] {
                let radius = one.from_usize(numerator) / one.from_usize(denominator);
                let result = certificate(&radius.0).unwrap();
                assert_eq!(
                    result,
                    geometry.prepare_certificate(&scale.0, 2048).unwrap()(&radius.0).unwrap()
                );
                // Lower working precision certifies the same conservative
                // endpoint signs; it must not define a different box policy.
                for precision in JointRange::PRECISIONS {
                    assert_eq!(
                        result,
                        geometry.prepare_certificate(&scale.0, precision).unwrap()(&radius.0)
                            .unwrap(),
                        "working precision {precision}, radius {radius}"
                    );
                }
                if result {
                    accepted += 1;
                } else {
                    declined += 1;
                }
            }
            assert!(accepted > 0 && declined > 0);
            let mut changed = geometry.clone();
            changed.energy_sums[0] = one.0.clone();
            let changed_certificate = changed.prepare_certificate(&scale.0, 2048).unwrap();
            // E0+E1 >= m0+m1=5 independently excludes this changed C1=1
            // zero-radius disk; the original rational fixture remains regular.
            assert!(!changed_certificate(&one.zero().0).unwrap());
            assert!(certificate(&one.zero().0).unwrap());
            assert!(geometry.prepare_certificate(&one.zero().0, 2048).is_err());
            changed.shifts[1] = [one.from_usize(6).0, one.zero().0, one.zero().0];
            // False Gram still precedes the scale denominator check.
            assert!(!changed.prepare_certificate(&one.zero().0, 2048).unwrap()(&one.0).unwrap());
        }
        check::<f64>();
        check::<QuadFloat>();
        check::<ArbPrec>();
    }

    #[test]
    fn joint_sampling_certificate_refines_before_changing_the_native_radius() {
        crate::initialisation::test_initialise().unwrap();
        fn check<T: FloatLike>(program: SamplingExpressionEvaluator) {
            let one = F::<T>::default().one();
            let two = one.from_usize(2);
            let delta = two.powi(-100);
            let radius = two.powi(-250);
            let geometry = SharedEnergyJointGeometry {
                shifts: [
                    [one.from_usize(4).0, one.zero().0, one.zero().0],
                    [one.from_usize(4).0, delta.0, one.zero().0],
                ],
                masses: [one.zero().0, one.zero().0, one.zero().0],
                energy_sums: [one.from_usize(8).0, one.from_usize(8).0],
            };
            // The exact Gram is 16*2^-200 > 0; at 128 bits the subtraction
            // cannot certify its sign. For h=z=0 the limiting k=3 and D=36
            // are positive, and the fixed2048 owner admits this tiny disk.
            let error = match geometry.prepare_certificate(&one.0, 128) {
                Ok(_) => panic!("the unresolved low-precision Gram was accepted"),
                Err(error) => error,
            };
            assert!(matches!(
                error.downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { .. })
            ));
            assert!(geometry.prepare_certificate(&one.0, 2048).unwrap()(&radius.0).unwrap());
            let mut chart = map(geometry, program);
            chart.normal_scale = one.0;
            chart.max_radius = radius.0.clone();
            let (_, selected) = chart
                .prepare(&mut SamplingMapContext::detached(&[]))
                .unwrap();
            assert_eq!(
                selected,
                Some(radius.0),
                "uncertainty must not halve the disk"
            );
        }
        let program = SharedEnergyJointMap::<f64>::compile_program().unwrap();
        check::<f64>(program.clone());
        check::<QuadFloat>(program.clone());
        check::<ArbPrec>(program);
    }

    #[test]
    fn joint_sampling_support_refines_separated_native_bits_without_moving_the_boundary() {
        use symbolica::domains::float::DoubleFloat;

        fn check<T: FloatLike>(energy_sum: T) {
            let one = F::<T>::default().one();
            let delta = one.from_usize(2).powi(-400);
            let geometry = SharedEnergyJointGeometry {
                shifts: [
                    [one.from_usize(4).0, one.zero().0, one.zero().0],
                    [one.zero().0, one.from_usize(4).0, one.zero().0],
                ],
                masses: [one.zero().0, one.zero().0, one.zero().0],
                energy_sums: [energy_sum, one.from_usize(8).0],
            };
            let point = [one.zero().0, one.zero().0, one.from_usize(3).0];
            // Exact energies are (3,5,5), hence H=-2^-400 and Z=0. The
            // original low limb is unresolved at 128 bits and must survive
            // refinement from the native input, including a separated Quad limb.
            let low_h =
                JointRange::integer(8, 128).sub(&JointRange::native(&geometry.energy_sums[0], 128));
            assert!(low_h.lower[0] < 0 && low_h.lower[1] == 0);
            for (factor, expected) in [
                (-1, std::cmp::Ordering::Greater),
                (1, std::cmp::Ordering::Less),
            ] {
                let radius = &delta * one.from_usize(2).powi(factor);
                assert_eq!(
                    geometry
                        .support_relation(&point, &radius.0, &one.0)
                        .unwrap(),
                    expected
                );
            }
            let error = geometry
                .support_relation(&point, &delta.0, &one.0)
                .unwrap_err();
            assert!(matches!(
                error.downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { detail })
                    if detail.contains("normal-disk boundary")
            ));
        }
        check(QuadFloat::from(DoubleFloat::from_compensated_sum(
            8.0,
            2.0_f64.powi(-400),
        )));
        let one = F::<ArbPrec>::default().one();
        check((one.from_usize(8) + one.from_usize(2).powi(-400)).0);
    }

    #[test]
    fn joint_sampling_full_circle_eager_jacobian_matches_cartesian_difference() {
        crate::initialisation::test_initialise().unwrap();
        let map = map(
            geometry::<f64>(),
            SharedEnergyJointMap::<f64>::compile_program().unwrap(),
        );
        for phi in [0.0, 0.13, 0.5, 0.83] {
            let cube = [0.41, 0.27, phi];
            let output = map
                .forward(&cube, &mut SamplingMapContext::detached(&[]))
                .unwrap();
            let inverse = map
                .inverse(&output.point, &mut SamplingMapContext::detached(&[]))
                .unwrap()
                .unwrap();
            assert!((output.jacobian * inverse.inverse_jacobian - 1.0).abs() < 2e-9);
            for (a, b) in cube.iter().zip(&inverse.coordinates) {
                let difference = (a - b).abs();
                assert!(difference.min((1.0 - difference).abs()) < 2e-9);
            }
            // Three-dimensional Cartesian FD includes both square-root circle
            // branches and their finite joining points. Angular input is periodic.
            let step = 1e-5;
            let mut matrix = [[0.0; 3]; 3];
            for column in 0..3 {
                let mut plus = cube;
                let mut minus = cube;
                plus[column] += step;
                minus[column] -= step;
                if column == 2 {
                    plus[column] = plus[column].rem_euclid(1.0);
                    minus[column] = minus[column].rem_euclid(1.0);
                }
                let plus = map
                    .forward(&plus, &mut SamplingMapContext::detached(&[]))
                    .unwrap();
                let minus = map
                    .forward(&minus, &mut SamplingMapContext::detached(&[]))
                    .unwrap();
                for (row, values) in matrix.iter_mut().enumerate() {
                    values[column] = (plus.point[row] - minus.point[row]) / (2.0 * step);
                }
            }
            let determinant = matrix[0][0]
                * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
                - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
                + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
            assert!(
                (determinant.abs() / output.jacobian - 1.0).abs() < 3e-6,
                "phi={phi}: FD={determinant}, eager={}",
                output.jacobian
            );
        }
        assert!(
            map.inverse(
                &[100.0, 100.0, 100.0],
                &mut SamplingMapContext::detached(&[])
            )
            .unwrap()
            .is_none()
        );
    }

    #[test]
    fn joint_sampling_active_program_matches_full_dual_at_arb_precision() {
        crate::initialisation::test_initialise().unwrap();
        let mut active = SharedEnergyJointMap::<ArbPrec>::compile_program().unwrap();
        let mut full = SharedEnergyJointMap::<ArbPrec>::compile_program_with_derivatives(
            &(0..17).collect::<Vec<_>>(),
        )
        .unwrap();
        let one = F::<ArbPrec>::default().one();
        let mut geometry = geometry::<ArbPrec>();
        geometry.energy_sums[0] =
            (F(geometry.energy_sums[0].clone()) + &one / one.from_usize(10).powi(25)).0;
        let chart = map(geometry, active.clone());
        let (geometry, radius) = chart
            .prepare(&mut SamplingMapContext::detached(&[]))
            .unwrap();
        let radius = radius.unwrap();
        assert_eq!(geometry.masses[2], one.zero().0);
        for numerator in [0, 13, 50, 83] {
            let cube = [41, 27, numerator].map(|n| (one.from_usize(n) / one.from_usize(100)).0);
            let parameters = chart.parameters(&geometry, &radius, &cube);
            let scalar = active.evaluate(&parameters).unwrap();
            let reduced = active
                .evaluate_with_real_jacobian(&parameters, None)
                .unwrap();
            let previous = full
                .evaluate_with_real_jacobian(&parameters, Some(&[0, 1, 2]))
                .unwrap();
            // The old seventeen-direction program is compiled from the same
            // factory. Require bit-identical primals, active derivatives and
            // determinant at identical Arb inputs, including a massless edge
            // and both circle seams. Scalar reconstruction skips vectorization
            // and its instruction regrouping: the observed difference is one
            // Arb ulp, so only this separately compiled primal has a small
            // native-rounding bound. Physical and J*q tolerances are unchanged.
            assert_eq!(reduced, previous, "circle fraction {numerator}/100");
            for (value, previous) in scalar.iter().zip(previous.values) {
                let previous = F(previous);
                let difference = (&value.re - &previous).abs();
                let bound = one.epsilon() * previous.abs().max(one.clone()) * one.from_usize(8);
                assert!(
                    difference <= bound,
                    "circle fraction {numerator}/100: scalar difference {difference}, native bound {bound}"
                );
                assert_eq!(value.im, one.zero());
            }
        }
    }

    #[test]
    fn joint_sampling_native_bindings_keep_sub_double_geometry_and_worker_buffers() {
        crate::initialisation::test_initialise().unwrap();
        fn check<T: FloatLike>(program: SamplingExpressionEvaluator) {
            let one = F::<T>::default().one();
            let delta = one.clone() / one.from_usize(10).powi(25);
            let original = geometry::<T>();
            let mut changed = original.clone();
            changed.energy_sums[0] = (F(changed.energy_sums[0].clone()) + &delta).0;
            let original = map(original, program.clone());
            let changed = map(changed, program);
            let cube = [0.41, 0.27, 0.13].map(T::from_f64_exact_binary);
            let base = original
                .forward(&cube, &mut SamplingMapContext::detached(&[]))
                .unwrap();
            let output = changed
                .forward(&cube, &mut SamplingMapContext::detached(&[]))
                .unwrap();
            let inverse = changed
                .inverse(&output.point, &mut SamplingMapContext::detached(&[]))
                .unwrap()
                .unwrap();
            let error = (F(output.jacobian.clone()) * F(inverse.inverse_jacobian) - &one).abs();
            assert!(
                error < one.clone() / one.from_usize(10).powi(27),
                "native J*q error {error}"
            );
            assert!(
                output.point.iter().zip(&base.point).any(|(a, b)| a != b),
                "native geometry perturbation was narrowed away"
            );
            // A clone must own independent mutable evaluator buffers. Holding
            // the parent lock cannot prevent a worker from evaluating.
            let worker = changed.clone();
            let _parent_lock = changed.program.lock().unwrap();
            let threaded = std::thread::spawn(move || {
                worker
                    .forward(&cube, &mut SamplingMapContext::detached(&[]))
                    .unwrap()
            })
            .join()
            .unwrap();
            assert_eq!(threaded.point, output.point);
        }
        let program = SharedEnergyJointMap::<f64>::compile_program().unwrap();
        check::<QuadFloat>(program.clone());
        check::<ArbPrec>(program);
    }

    #[test]
    fn joint_sampling_policy_rejects_another_valid_radius_and_ordinary_branch() {
        use crate::integrands::process::{SamplingChannelId, SamplingChannelRuntimeContexts};
        crate::initialisation::test_initialise().unwrap();
        let program = SharedEnergyJointMap::<f64>::compile_program().unwrap();
        let mut canonical = map(geometry::<ArbPrec>(), program.clone());
        canonical.max_radius = F::<ArbPrec>::default().from_usize(100).0;
        let mut metadata = crate::integrands::evaluation::EvaluationMetaData::new_empty();
        metadata.sampling_proposal_policies.begin_collection();
        let (_, radius) = {
            let mut rows =
                SamplingChannelRuntimeContexts::for_draw(1, 0, SamplingChannelId(0), &mut metadata);
            canonical
                .prepare(&mut rows.for_channel(SamplingChannelId(0)).unwrap())
                .unwrap()
        };
        let radius = radius.unwrap();
        assert!(radius < canonical.max_radius);
        metadata.sampling_proposal_policies.seal();
        let mut enlarged = geometry::<ArbPrec>();
        let two = F::<ArbPrec>::default().from_usize(2);
        for value in enlarged
            .shifts
            .iter_mut()
            .flatten()
            .chain(&mut enlarged.masses)
            .chain(&mut enlarged.energy_sums)
        {
            *value = (F(value.clone()) * &two).0;
        }
        // The retained radius remains certified in the enlarged geometry.
        // Its own policy chooses twice that radius, so accepting validity alone
        // would silently replace the normalized law for this original draw.
        assert!(
            enlarged
                .prepare_certificate(&canonical.normal_scale, 2048)
                .unwrap()(&radius)
            .unwrap()
        );
        let mut changed = map(enlarged, program.clone());
        changed.max_radius = canonical.max_radius.clone();
        let (_, changed_radius) = changed
            .prepare(&mut SamplingMapContext::detached(&[]))
            .unwrap();
        assert_eq!(changed_radius.unwrap(), (F(radius) * &two).0);
        let mut ordinary_geometry = geometry::<ArbPrec>();
        ordinary_geometry.energy_sums = [two.one().0.clone(), two.one().0];
        let ordinary = map(ordinary_geometry, program);
        for changed in [&changed, &ordinary] {
            let mut rows =
                SamplingChannelRuntimeContexts::for_draw(1, 0, SamplingChannelId(0), &mut metadata);
            let error = changed
                .prepare(&mut rows.for_channel(SamplingChannelId(0)).unwrap())
                .unwrap_err();
            assert!(
                error.downcast_ref::<SamplingEvaluationError>().is_some(),
                "{error:?}"
            );
        }
        let mut rows =
            SamplingChannelRuntimeContexts::for_draw(1, 0, SamplingChannelId(0), &mut metadata);
        canonical
            .prepare(&mut rows.for_channel(SamplingChannelId(0)).unwrap())
            .unwrap();
        assert_eq!(metadata.sampling_proposal_policies.len(), 1);
    }

    #[test]
    fn joint_sampling_declined_domain_reuses_normalized_ordinary_map() {
        crate::initialisation::test_initialise().unwrap();
        let mut geometry = geometry::<f64>();
        geometry.energy_sums = [1.0, 1.0];
        let map = map(
            geometry,
            SharedEnergyJointMap::<f64>::compile_program().unwrap(),
        );
        let cube = [0.41, 0.27, 0.13];
        let output = map
            .forward(&cube, &mut SamplingMapContext::detached(&[]))
            .unwrap();
        let ordinary = SamplingMapComponent::forward(
            &map.fallback,
            &cube,
            &mut SamplingMapContext::detached(&[]),
        )
        .unwrap();
        assert_eq!(output.point, ordinary.point);
        assert_eq!(output.jacobian, ordinary.jacobian);
        assert!(output.diagnostics.iter().any(|d| d.contains("fallback")));
        let inverse = map
            .inverse(&output.point, &mut SamplingMapContext::detached(&[]))
            .unwrap()
            .unwrap();
        assert_eq!(
            map.inverse_density(&output.point, &mut SamplingMapContext::detached(&[]))
                .unwrap(),
            Some(inverse.inverse_jacobian),
        );
        assert!((output.jacobian * inverse.inverse_jacobian - 1.0).abs() < 1e-12);
        assert_eq!(map.contract().support, SamplingSupport::Restricted);
        let boundary = map
            .forward(&[0.41, 0.0, 0.13], &mut SamplingMapContext::detached(&[]))
            .unwrap_err();
        assert!(boundary.downcast_ref::<SamplingEvaluationError>().is_some());
        assert!(
            boundary
                .to_string()
                .contains("ordinary fallback angular boundary")
        );

        // A cached false Gram must not bypass the original native halving
        // guard, even though every candidate would decline focusing.
        let mut tiny = map.clone();
        let mut collinear = (tiny.geometry)(&[]).unwrap();
        collinear.shifts[1] = [6.0, 0.0, 0.0];
        tiny.geometry = Arc::new(move |_| Ok(collinear.clone()));
        tiny.max_radius = f64::from_bits(1);
        let error = tiny
            .prepare(&mut SamplingMapContext::detached(&[]))
            .unwrap_err();
        assert!(matches!(
            error.downcast_ref::<SamplingEvaluationError>(),
            Some(SamplingEvaluationError::Unrepresentable { operation: "joint disk", detail })
                if detail.contains("underflow")
        ));
    }

    #[test]
    fn joint_sampling_exact_normal_origin_is_an_error_not_outside_support() {
        crate::initialisation::test_initialise().unwrap();
        let geometry = SharedEnergyJointGeometry {
            shifts: [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0]],
            masses: [0.0, 0.0, 0.0],
            energy_sums: [8.0, 8.0],
        };
        let map = map(
            geometry,
            SharedEnergyJointMap::<f64>::compile_program().unwrap(),
        );
        // Original unsquared energies are exactly (3,5,5), hence h=z=0.
        let error = map
            .inverse(&[0.0, 0.0, 3.0], &mut SamplingMapContext::detached(&[]))
            .unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());
        assert!(error.to_string().contains("R=0"), "{error}");
        let density_error = map
            .inverse_density(&[0.0, 0.0, 3.0], &mut SamplingMapContext::detached(&[]))
            .unwrap_err();
        assert!(
            density_error
                .downcast_ref::<SamplingEvaluationError>()
                .is_some()
        );
        assert_eq!(density_error.to_string(), error.to_string());
    }
}
