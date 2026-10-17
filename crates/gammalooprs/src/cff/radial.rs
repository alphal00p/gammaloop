//! Retained affine energy equations shared by LU roots, jets and host certificates.

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use rug::{Float, float::Round};
use symbolica::domains::{
    dual::HyperDual,
    float::{FloatLike as SymFloatLike, Real},
};

use super::esurface::{Esurface, ExternalShift};
use crate::{
    graph::LoopMomentumBasis,
    integrands::process::sampling::maps::SamplingEvaluationError,
    momentum::{
        Rotatable, Rotation, ThreeMomentum,
        sample::{ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopMomenta},
    },
    utils::{
        ArbPrec, F, FloatLike, compute_loop_part, compute_t_part_of_shift_part,
        hyperdual_utils::new_constant,
        newton_solver::{
            NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity,
            SafeguardedNewtonError,
        },
    },
};

/// One ordered occurrence of an on-shell energy in a represented radial equation.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct RadialEnergy<T: FloatLike> {
    pub(super) edge: EdgeIndex,
    pub(super) coefficient: i64,
    pub(super) velocity: ThreeMomentum<F<T>>,
    pub(super) offset: ThreeMomentum<F<T>>,
    pub(super) mass: F<T>,
}
// Ephemeral endpoint tuples for the same affine energy fold, with no stored ray.
type EnclosedRayEnergy = (i64, [[Float; 2]; 3], [[Float; 2]; 3], [Float; 2]);

/// One represented affine energy equation, retaining its ordered edge occurrences.
/// This is native geometry, not a cross-precision or selected-host transport record.
#[derive(Debug, Clone)]
pub(crate) struct EnergySurfaceRay<T: FloatLike> {
    // A Vec deliberately preserves repeated energy occurrences; an edge-keyed
    // map would not. Algebraic reduction is only a derived root-isolation view.
    pub(super) energies: Vec<RadialEnergy<T>>,
    pub(super) shift: F<T>,
}

impl EnergySurfaceRay<ArbPrec> {
    /// Retain the canonical ordered equation while materializing its native
    /// coefficients. Physical adoption still authenticates and checks the ray.
    pub(crate) fn materialize<T: FloatLike>(&self) -> Result<EnergySurfaceRay<T>> {
        let spatial = |p: &ThreeMomentum<F<ArbPrec>>| -> Result<ThreeMomentum<F<T>>> {
            Ok(ThreeMomentum::new(
                F::from_arb(&p.px.0)?,
                F::from_arb(&p.py.0)?,
                F::from_arb(&p.pz.0)?,
            ))
        };
        Ok(EnergySurfaceRay {
            energies: self
                .energies
                .iter()
                .map(|energy| {
                    Ok(RadialEnergy {
                        edge: energy.edge,
                        coefficient: energy.coefficient,
                        velocity: spatial(&energy.velocity)?,
                        offset: spatial(&energy.offset)?,
                        mass: F::from_arb(&energy.mass.0)?,
                    })
                })
                .collect::<Result<_>>()?,
            shift: F::from_arb(&self.shift.0)?,
        })
    }
}

impl<T: FloatLike> EnergySurfaceRay<T> {
    /// Route once, retaining the equation shared by localization and residue jets.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn routed(
        energies: impl Iterator<Item = (EdgeIndex, i64)>,
        external_shift: &ExternalShift,
        velocity: &LoopMomenta<F<T>>,
        center: &LoopMomenta<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Self {
        let spatial: ExternalThreeMomenta<F<T>> = external_momenta
            .iter()
            .map(|momentum| momentum.spatial.clone())
            .collect();
        Self {
            energies: energies
                .map(|(edge, coefficient)| {
                    let signature = &lmb.edge_signatures[edge];
                    RadialEnergy {
                        edge,
                        coefficient,
                        velocity: compute_loop_part(&signature.internal, velocity),
                        offset: signature.compute_momentum(center, &spatial),
                        mass: masses[edge].clone(),
                    }
                })
                .collect(),
            shift: external_shift
                .iter()
                .map(|(edge, coefficient)| {
                    external_momenta[ExternalIndex(0)]
                        .temporal
                        .value
                        .from_i64(*coefficient)
                        * compute_t_part_of_shift_part(
                            &lmb.edge_signatures[*edge].external,
                            external_momenta,
                        )
                })
                .reduce(|sum, term| sum + term)
                .unwrap_or_else(|| external_momenta[ExternalIndex(0)].temporal.value.zero()),
        }
    }

    /// Embed the represented equation exactly, preserving each ordered energy
    /// occurrence. No routing, coefficient reconstruction or root solve occurs.
    pub(crate) fn to_arb_exact(&self) -> Result<EnergySurfaceRay<ArbPrec>> {
        let spatial = |p: &ThreeMomentum<F<T>>| -> Result<ThreeMomentum<F<ArbPrec>>> {
            Ok(ThreeMomentum::new(
                p.px.to_arb_exact()?,
                p.py.to_arb_exact()?,
                p.pz.to_arb_exact()?,
            ))
        };
        Ok(EnergySurfaceRay {
            energies: self
                .energies
                .iter()
                .map(|energy| {
                    Ok(RadialEnergy {
                        edge: energy.edge,
                        coefficient: energy.coefficient,
                        velocity: spatial(&energy.velocity)?,
                        offset: spatial(&energy.offset)?,
                        mass: energy.mass.to_arb_exact()?,
                    })
                })
                .collect::<Result<_>>()?,
            shift: self.shift.to_arb_exact()?,
        })
    }

    const ENCLOSURE_PRECISION: u32 = 2048;

    // These directed operations serve the represented and original-routed
    // affine-energy checks; they do not introduce a second certificate engine.
    fn enclosed_add(a: &[Float; 2], b: &[Float; 2]) -> [Float; 2] {
        [
            Float::with_val_round(Self::ENCLOSURE_PRECISION, &a[0] + &b[0], Round::Down).0,
            Float::with_val_round(Self::ENCLOSURE_PRECISION, &a[1] + &b[1], Round::Up).0,
        ]
    }
    fn enclosed_mul(a: &[Float; 2], b: &[Float; 2]) -> [Float; 2] {
        let lower = a
            .iter()
            .flat_map(|x| {
                b.iter().map(move |y| {
                    Float::with_val_round(Self::ENCLOSURE_PRECISION, x * y, Round::Down).0
                })
            })
            .reduce(|a, b| a.min(&b))
            .unwrap();
        let upper = a
            .iter()
            .flat_map(|x| {
                b.iter().map(move |y| {
                    Float::with_val_round(Self::ENCLOSURE_PRECISION, x * y, Round::Up).0
                })
            })
            .reduce(|a, b| a.max(&b))
            .unwrap();
        [lower, upper]
    }
    fn enclosed_square(a: &[Float; 2]) -> [Float; 2] {
        let upper = a[0].clone().abs().max(&a[1].clone().abs());
        let lower = if a[0] <= 0 && a[1] >= 0 {
            Float::with_val(Self::ENCLOSURE_PRECISION, 0)
        } else {
            a[0].clone().abs().min(&a[1].clone().abs())
        };
        [
            Float::with_val_round(Self::ENCLOSURE_PRECISION, &lower * &lower, Round::Down).0,
            Float::with_val_round(Self::ENCLOSURE_PRECISION, &upper * &upper, Round::Up).0,
        ]
    }
    fn enclosed_divide(a: &[Float; 2], b: &[Float; 2]) -> [Float; 2] {
        debug_assert!(b[0] > 0);
        let reciprocal = [
            Float::with_val_round(Self::ENCLOSURE_PRECISION, b[1].recip_ref(), Round::Down).0,
            Float::with_val_round(Self::ENCLOSURE_PRECISION, b[0].recip_ref(), Round::Up).0,
        ];
        Self::enclosed_mul(a, &reciprocal)
    }

    pub(super) fn enclosed_native(value: &F<T>) -> [Float; 2] {
        let (lo, hi) = value.0.mpfr_enclosure(Self::ENCLOSURE_PRECISION);
        [lo, hi]
    }

    pub(super) fn enclosed_sum(terms: impl Iterator<Item = (i64, [Float; 2])>) -> [Float; 2] {
        let zero = Float::with_val(Self::ENCLOSURE_PRECISION, 0);
        terms.fold([zero.clone(), zero], |sum, (coefficient, value)| {
            let coefficient = Float::with_val(Self::ENCLOSURE_PRECISION, coefficient);
            Self::enclosed_add(
                &sum,
                &Self::enclosed_mul(&[coefficient.clone(), coefficient], &value),
            )
        })
    }

    pub(super) fn evaluate_enclosed(
        energies: impl Iterator<Item = EnclosedRayEnergy>,
        shift: [Float; 2],
        radius: &[Float; 2],
        with_derivative: bool,
    ) -> Result<([Float; 2], [Float; 2])> {
        let zero = Float::with_val(Self::ENCLOSURE_PRECISION, 0);
        let point = |x: &Float| [x.clone(), x.clone()];
        let uncertain = |detail: &str| SamplingEvaluationError::UncertainGeometry {
            detail: format!("radial energy enclosure: {detail}"),
        };
        let mut value = shift;
        let mut slope = point(&zero);
        for (coefficient, velocity, offset, mass) in energies {
            // Exact zero velocity makes the energy constant, even at E=0.
            // Never infer this from a norm that may have underflowed.
            let constant = velocity.iter().flatten().all(|value| value == &zero);
            let mut norm = Self::enclosed_square(&mass);
            let mut dot = point(&zero);
            for (v, b) in velocity.into_iter().zip(offset) {
                let q = Self::enclosed_add(&Self::enclosed_mul(&v, radius), &b);
                norm = Self::enclosed_add(&norm, &Self::enclosed_square(&q));
                dot = Self::enclosed_add(&dot, &Self::enclosed_mul(&v, &q));
            }
            if norm.iter().chain(&dot).any(|x| !x.is_finite()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "radial energy enclosure",
                    detail: "nonfinite directed energy or numerator".into(),
                }
                .into());
            }
            let energy = [
                Float::with_val_round(Self::ENCLOSURE_PRECISION, norm[0].sqrt_ref(), Round::Down).0,
                Float::with_val_round(Self::ENCLOSURE_PRECISION, norm[1].sqrt_ref(), Round::Up).0,
            ];
            let coefficient = Float::with_val(Self::ENCLOSURE_PRECISION, coefficient);
            let coefficient = point(&coefficient);
            value = Self::enclosed_add(&value, &Self::enclosed_mul(&coefficient, &energy));
            if with_derivative && !constant {
                if energy[0] <= 0 {
                    return Err(uncertain(
                        "energy lower bound touches zero on the fixed candidate interval",
                    )
                    .into());
                }
                slope = Self::enclosed_add(
                    &slope,
                    &Self::enclosed_mul(&coefficient, &Self::enclosed_divide(&dot, &energy)),
                );
            }
        }
        if value.iter().chain(&slope).any(|x| !x.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "radial energy enclosure",
                detail: "nonfinite directed residual or derivative".into(),
            }
            .into());
        }
        Ok((value, slope))
    }

    fn evaluate_interval(
        &self,
        radius: &[Float; 2],
        with_derivative: bool,
    ) -> Result<([Float; 2], [Float; 2])> {
        let native = Self::enclosed_native;
        Self::evaluate_enclosed(
            self.energies.iter().map(|energy| {
                (
                    energy.coefficient,
                    [
                        &energy.velocity.px,
                        &energy.velocity.py,
                        &energy.velocity.pz,
                    ]
                    .map(native),
                    [&energy.offset.px, &energy.offset.py, &energy.offset.pz].map(native),
                    native(&energy.mass),
                )
            }),
            native(&self.shift),
            radius,
            with_derivative,
        )
    }

    /// Compare original isotropic normal equations on the actual completed
    /// canonical/native LU points. Half the relative budget controls their
    /// displacement; half controls the native point's host constraint defect.
    /// This does not certify arbitrary multiplier functions or raised jets.
    pub(crate) fn verify_normal_alignment(
        canonical_normals: [[Float; 2]; 2],
        native_normals: [[Float; 2]; 2],
        native_host_residual: [Float; 2],
        tolerance: f64,
    ) -> Result<()> {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(eyre!(
                "hosted normal alignment requires a finite positive budget"
            ));
        }
        let bounds = canonical_normals
            .iter()
            .chain(&native_normals)
            .chain(std::iter::once(&native_host_residual));
        for bound in bounds {
            if bound.iter().any(|x| !x.is_finite()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "hosted normal alignment",
                    detail: "nonfinite directed normal or host residual".into(),
                }
                .into());
            }
            if bound[0] > bound[1] {
                return Err(eyre!("hosted normal alignment requires ordered bounds"));
            }
        }
        let radius_squared = Self::enclosed_add(
            &Self::enclosed_square(&canonical_normals[0]),
            &Self::enclosed_square(&canonical_normals[1]),
        );
        let radius_lower = Float::with_val_round(
            Self::ENCLOSURE_PRECISION,
            radius_squared[0].sqrt_ref(),
            Round::Down,
        )
        .0;
        if radius_lower <= 0 {
            return Err(SamplingEvaluationError::UncertainGeometry {
                detail: "hosted normal alignment has no strictly positive canonical radius bound"
                    .into(),
            }
            .into());
        }
        let difference: [_; 2] = std::array::from_fn(|axis| {
            let canonical = &canonical_normals[axis];
            Self::enclosed_add(
                &native_normals[axis],
                &[-canonical[1].clone(), -canonical[0].clone()],
            )
        });
        let error_squared = Self::enclosed_add(
            &Self::enclosed_square(&difference[0]),
            &Self::enclosed_square(&difference[1]),
        );
        let error_upper = Float::with_val_round(
            Self::ENCLOSURE_PRECISION,
            error_squared[1].sqrt_ref(),
            Round::Up,
        )
        .0;
        let budget = Float::with_val(Self::ENCLOSURE_PRECISION, tolerance);
        let half_budget =
            Float::with_val_round(Self::ENCLOSURE_PRECISION, &budget / 2, Round::Down).0;
        let allowed = Float::with_val_round(
            Self::ENCLOSURE_PRECISION,
            &half_budget * &radius_lower,
            Round::Down,
        )
        .0;
        let host_upper = native_host_residual[0]
            .clone()
            .abs()
            .max(&native_host_residual[1].clone().abs());
        if [&radius_lower, &error_upper, &host_upper, &allowed]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "hosted normal alignment",
                detail: "nonfinite directed comparison".into(),
            }
            .into());
        }
        if error_upper > allowed || host_upper > allowed {
            return Err(SamplingEvaluationError::UncertainGeometry {
                detail: "original normal displacement or completed host defect exceeds its half sampling budget".into(),
            }.into());
        }
        Ok(())
    }

    #[inline]
    pub(crate) fn evaluate(&self, radius: &F<T>) -> (F<T>, F<T>) {
        let zero = radius.zero();
        let (derivative, energy_sum) = self
            .energies
            .iter()
            .map(
                |RadialEnergy {
                     coefficient,
                     velocity,
                     offset,
                     mass,
                     ..
                 }| {
                    let momentum = velocity * radius + offset;
                    let energy = (momentum.norm_squared() + mass * mass).sqrt();
                    // At a massless endpoint E(r)=r|v| the radial right derivative is
                    // |v|, whereas the two-sided formula q.v/E would evaluate 0/0.
                    // Share this convention with sampling charts, using exact source
                    // zeros: a computed E=0 can instead be numerical underflow.
                    let derivative =
                        if velocity.px == zero && velocity.py == zero && velocity.pz == zero {
                            zero.clone()
                        } else if radius == &zero
                            && mass == &zero
                            && momentum.px == zero
                            && momentum.py == zero
                            && momentum.pz == zero
                        {
                            velocity.norm_squared().sqrt()
                        } else {
                            momentum * velocity / &energy
                        };
                    if *coefficient == 1 {
                        (derivative, energy)
                    } else {
                        let coefficient = radius.from_i64(*coefficient);
                        (derivative * &coefficient, energy * coefficient)
                    }
                },
            )
            .fold(
                (zero.clone(), zero.clone()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );
        (energy_sum + &self.shift, derivative)
    }

    /// The positive-E LU seed policy, also used by native partial host rays.
    pub(crate) fn solve_lu_cut(
        &self,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> std::result::Result<NewtonIterationResult<T>, SafeguardedNewtonError<T>> {
        assert!(
            self.energies.iter().all(|energy| energy.coefficient == 1),
            "the positive LU seed policy requires an E-surface ray"
        );
        let guess = Esurface::radius_guess_from_terms(
            &self.shift,
            self.energies.iter().map(|energy| {
                (
                    energy.velocity.norm_squared(),
                    energy.velocity.clone() * &energy.offset,
                )
            }),
        );
        crate::debug_tags!(#integration, #cut, #solver;
            radial_root = %identity,
            initial_guess = %guess,
            residual_tolerance = %(e_cm * guess.epsilon() * guess.from_i64((4 * self.energies.len() + 1) as i64)),
            "LU radial root setup"
        );
        // The energy equation is a sum of on-shell square roots.  Include a
        // small forward-error budget for the momentum norm, square root, and
        // accumulation of each edge instead of treating the whole sum as one
        // elementary operation.  This remains a caller-specific tolerance;
        // generic safeguarded-Newton users retain their original tolerance.
        self.refine_root(&guess.zero(), &guess, 1, 64, e_cm, diagnostics, identity)
    }

    #[allow(clippy::too_many_arguments)]
    fn refine_root(
        &self,
        lower: &F<T>,
        upper: &F<T>,
        crossing: i64,
        max_bracket_expansions: usize,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> std::result::Result<NewtonIterationResult<T>, SafeguardedNewtonError<T>> {
        let mut result = diagnostics.solve(
            identity,
            lower,
            upper,
            |radius| {
                let (value, derivative) = self.evaluate(radius);
                if crossing > 0 {
                    (value, derivative)
                } else {
                    (-value, -derivative)
                }
            },
            &upper.from_i64((4 * self.energies.len() + 1) as i64),
            2000,
            max_bracket_expansions,
            e_cm,
        )?;
        if crossing < 0 {
            result.derivative_at_solution = -result.derivative_at_solution;
            result.error_of_function = -result.error_of_function;
        }
        Ok(result)
    }

    /// Certify this retained candidate against both represented rays. This is
    /// a fixed directed check, not another root solve. The exponent bounds all
    /// pulled-back host-null dimensions; it need not be the actual map dimension.
    /// Higher raised jets and complete conditional densities are not enclosed.
    pub(crate) fn verify_lu_candidate(
        &self,
        completed: &Self,
        solution: &NewtonIterationResult<T>,
        dimension_bound: usize,
        tolerance: &F<T>,
    ) -> Result<()> {
        if [self, completed]
            .iter()
            .any(|ray| ray.energies.iter().any(|energy| energy.coefficient != 1))
        {
            return Err(eyre!("LU host adoption requires positive E-surface rays"));
        }
        let precision = Self::ENCLOSURE_PRECISION;
        let uncertain = |detail: &str| SamplingEvaluationError::UncertainGeometry {
            detail: format!("LU host adoption: {detail}"),
        };
        if self
            .energies
            .iter()
            .map(|entry| (entry.edge, entry.coefficient))
            .ne(completed
                .energies
                .iter()
                .map(|entry| (entry.edge, entry.coefficient)))
        {
            return Err(eyre!(
                "LU host adoption has incompatible ordered energy occurrences"
            ));
        }
        for ray in [self, completed] {
            if !ray.shift.0.is_finite()
                || ray.energies.iter().any(|energy| {
                    let (v, b, m) = (&energy.velocity, &energy.offset, &energy.mass);
                    [&v.px, &v.py, &v.pz, &b.px, &b.py, &b.pz, m]
                        .iter()
                        .any(|x| !x.0.is_finite())
                })
            {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "LU host adoption",
                    detail: "nonfinite represented ray coefficients".into(),
                }
                .into());
            }
        }
        let t = &solution.solution;
        let derivative = &solution.derivative_at_solution;
        if [t, derivative, tolerance].iter().any(|x| !x.0.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "LU host adoption",
                detail: "nonfinite candidate, derivative or budget".into(),
            }
            .into());
        }
        if t <= &t.zero() || derivative <= &t.zero() || tolerance <= &t.zero() {
            return Err(uncertain("candidate, derivative and budget must be positive").into());
        }
        let native = Self::enclosed_native;
        let point = |x: &Float| [x.clone(), x.clone()];
        let zero = Float::with_val(precision, 0);
        let one = Float::with_val(precision, 1);
        let delta = native(tolerance)[0]
            .clone()
            .min(&(Float::with_val(precision, 1) >> 2));
        let divisor = Float::with_val(precision, 4 * (dimension_bound + 1));
        let r = Float::with_val_round(precision, &delta / &divisor, Round::Down).0;
        let t0 = native(t);
        let lo_factor = Float::with_val_round(precision, &one - &r, Round::Down).0;
        let hi_factor = Float::with_val_round(precision, &one + &r, Round::Up).0;
        let interval = [
            Float::with_val_round(precision, &t0[0] * &lo_factor, Round::Down).0,
            Float::with_val_round(precision, &t0[1] * &hi_factor, Round::Up).0,
        ];
        if interval
            .iter()
            .chain(&t0)
            .chain([&delta, &r])
            .any(|x| !x.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "LU host certificate",
                detail: "nonfinite directed candidate interval".into(),
            }
            .into());
        }
        if interval[0] <= 0 || interval[0] >= interval[1] {
            return Err(uncertain("fixed candidate interval is not positive and distinct").into());
        }
        let allowed = [
            Float::with_val_round(precision, &one - &delta, Round::Up).0,
            Float::with_val_round(precision, &one + &delta, Round::Down).0,
        ];
        let residual_budget = Float::with_val_round(precision, &t0[0] * &r, Round::Down).0;
        for ray in [self, completed] {
            if ray.evaluate_interval(&point(&zero), false)?.0[1] >= 0
                || ray.evaluate_interval(&point(&interval[0]), false)?.0[1] >= 0
                || ray.evaluate_interval(&point(&interval[1]), false)?.0[0] <= 0
            {
                return Err(uncertain("origin/interior/exterior signs are not certified on the fixed candidate interval").into());
            }
            let (_, slope) = ray.evaluate_interval(&interval, true)?;
            if slope[0] <= 0 {
                return Err(uncertain(
                    "radial derivative is not certified positive throughout the interval",
                )
                .into());
            }
            let residual = ray.evaluate_interval(&t0, false)?.0;
            let residual = residual[0].clone().abs().max(&residual[1].clone().abs());
            if Float::with_val_round(precision, &residual / &slope[0], Round::Up).0
                > residual_budget
            {
                return Err(uncertain(
                    "candidate error exceeds its allocated relative root budget",
                )
                .into());
            }
            let ratio = Self::enclosed_divide(&slope, &native(derivative));
            if ratio.iter().any(|x| !x.is_finite())
                || ratio[0] < allowed[0]
                || ratio[1] > allowed[1]
            {
                return Err(
                    uncertain("first derivative variation exceeds the sampling budget").into(),
                );
            }
        }
        // Bound implemented/exact and its reciprocal explicitly, so the
        // statement does not depend on the chosen relative-error denominator.
        for ratio in [
            Self::enclosed_divide(&interval, &t0),
            Self::enclosed_divide(&t0, &interval),
        ] {
            let mut determinant_ratio = point(&one);
            for _ in 0..dimension_bound {
                determinant_ratio = Self::enclosed_mul(&determinant_ratio, &ratio);
            }
            if determinant_ratio.iter().any(|x| !x.is_finite())
                || determinant_ratio[0] < allowed[0]
                || determinant_ratio[1] > allowed[1]
            {
                return Err(uncertain(
                    "worst-case LU inverse-volume variation exceeds the sampling budget",
                )
                .into());
            }
        }
        Ok(())
    }

    /// Enumerate all simple positive roots when the ray reduces to
    /// `A*r + B + s*sqrt(a*r² + 2*b*r + c)`. This includes a massless sunset
    /// with its external momentum routed through one cut edge, also away from rest.
    /// Identical affine energies are combined before counting nonhomogeneous terms.
    /// Other rays and unresolved tangencies fail explicitly: a finite scan cannot
    /// establish that every physical root was found.
    pub(crate) fn positive_radius_roots(
        &self,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> eyre::Result<Vec<NewtonIterationResult<T>>> {
        let zero = e_cm.zero();
        let uncertain = |detail: &str| eyre::eyre!("{identity}: signed radial equation: {detail}");
        if !e_cm.0.is_finite() || e_cm <= &zero {
            return Err(uncertain("the residual scale must be finite and positive"));
        }
        // Collapse only an isolation view. The represented ray, native fold and
        // residue jets retain every signed occurrence and its physical edge.
        let mut terms: Vec<(i64, &RadialEnergy<T>)> = Vec::new();
        for energy in &self.energies {
            let (velocity, offset, mass) = (&energy.velocity, &energy.offset, &energy.mass);
            if [
                &velocity.px,
                &velocity.py,
                &velocity.pz,
                &offset.px,
                &offset.py,
                &offset.pz,
                mass,
            ]
            .into_iter()
            .any(|value| !value.0.is_finite())
            {
                return Err(uncertain("nonfinite routed energy"));
            }
            if let Some((multiplicity, _)) = terms.iter_mut().find(|(_, previous)| {
                (mass == &previous.mass || mass == &-previous.mass.clone())
                    && ((velocity == &previous.velocity && offset == &previous.offset)
                        || (velocity == &-previous.velocity.clone()
                            && offset == &-previous.offset.clone()))
            }) {
                *multiplicity += energy.coefficient;
            } else {
                terms.push((energy.coefficient, energy));
            }
        }
        terms.retain(|(sign, _)| *sign != 0);

        // Share the existing directed arithmetic. The classification below is
        // about the represented ray, and does not infer a missing root from a
        // failed Newton iteration or a floating-point discriminant alone.
        let add = Self::enclosed_add;
        let mul = Self::enclosed_mul;
        let square = Self::enclosed_square;
        let native = Self::enclosed_native;
        let neg = |x: &[Float; 2]| [-x[1].clone(), -x[0].clone()];
        let scale = |x: &[Float; 2], factor: i64| mul(x, &native(&zero.from_i64(factor)));
        let sqrt = |x: &[Float; 2]| -> eyre::Result<[Float; 2]> {
            if x[0] < 0 || x.iter().any(|value| !value.is_finite()) {
                return Err(uncertain("an energy or discriminant is not sign-certified"));
            }
            Ok([
                Float::with_val_round(x[0].prec(), x[0].sqrt_ref(), Round::Down).0,
                Float::with_val_round(x[1].prec(), x[1].sqrt_ref(), Round::Up).0,
            ])
        };
        let sign = |x: &[Float; 2]| -> eyre::Result<i64> {
            if x.iter().any(|value| !value.is_finite()) {
                Err(uncertain("nonfinite coefficient enclosure"))
            } else if x[0] > 0 {
                Ok(1)
            } else if x[1] < 0 {
                Ok(-1)
            } else if x[0] == 0 && x[1] == 0 {
                Ok(0)
            } else {
                Err(uncertain(
                    "a coefficient sign is unresolved at the current precision",
                ))
            }
        };
        let divide = |x: &[Float; 2], y: &[Float; 2]| -> eyre::Result<[Float; 2]> {
            match sign(y)? {
                1 => Ok(Self::enclosed_divide(x, y)),
                -1 => Ok(Self::enclosed_divide(&neg(x), &neg(y))),
                _ => Err(uncertain("a degenerate equation cannot be divided")),
            }
        };
        let coefficients = |v: &ThreeMomentum<F<T>>, b: &ThreeMomentum<F<T>>, mass: &F<T>| {
            let mut a = native(&zero);
            let mut dot = native(&zero);
            let mut c = square(&native(mass));
            for (velocity, offset) in [&v.px, &v.py, &v.pz].into_iter().zip([&b.px, &b.py, &b.pz]) {
                let velocity = native(velocity);
                let offset = native(offset);
                a = add(&a, &square(&velocity));
                dot = add(&dot, &mul(&velocity, &offset));
                c = add(&c, &square(&offset));
            }
            (a, dot, c)
        };
        let is_zero_vector = |v: &ThreeMomentum<F<T>>| v.px == zero && v.py == zero && v.pz == zero;
        let mut slope = native(&zero);
        if !self.shift.0.is_finite() {
            return Err(uncertain("nonfinite external energy shift"));
        }
        let mut shift = native(&self.shift);
        let mut curved = None;
        for (index, (multiplicity, energy)) in terms.iter().enumerate() {
            let (velocity, offset, mass) = (&energy.velocity, &energy.offset, &energy.mass);
            let (a, _, c) = coefficients(velocity, offset, mass);
            if is_zero_vector(velocity) {
                shift = add(&shift, &scale(&sqrt(&c)?, *multiplicity));
            } else if is_zero_vector(offset) && mass == &zero {
                slope = add(&slope, &scale(&sqrt(&a)?, *multiplicity));
            } else if curved.replace(index).is_some() {
                return Err(uncertain(
                    "unsupported ray with more than one distinct nonhomogeneous energy",
                ));
            }
        }

        let linear_roots = |a: &[Float; 2], b: &[Float; 2]| -> eyre::Result<Vec<[Float; 2]>> {
            if sign(a)? == 0 {
                return if sign(b)? == 0 {
                    Err(uncertain(
                        "the radial equation is identically zero or algebraically degenerate",
                    ))
                } else {
                    Ok(Vec::new())
                };
            }
            Ok(vec![divide(&neg(b), a)?])
        };
        let mut candidates = if let Some(index) = curved {
            let (multiplicity, energy) = &terms[index];
            let (velocity, offset, mass) = (&energy.velocity, &energy.offset, &energy.mass);
            let signed_slope = scale(&slope, *multiplicity);
            let signed_shift = scale(&shift, *multiplicity);
            if signed_slope[0] >= 0
                && signed_shift[0] >= 0
                && (signed_slope[0] > 0 || signed_shift[0] > 0)
            {
                // The unsquared equation cannot vanish on the positive ray:
                // its affine part has the same strict sign as the radical.
                return Ok(Vec::new());
            }
            let (a, b, c) = coefficients(velocity, offset, mass);
            let multiplicity_squared = multiplicity * multiplicity;
            let quadratic = add(&square(&slope), &neg(&scale(&a, multiplicity_squared)));
            let linear = scale(
                &add(&mul(&slope, &shift), &neg(&scale(&b, multiplicity_squared))),
                2,
            );
            let constant = add(&square(&shift), &neg(&scale(&c, multiplicity_squared)));
            if sign(&quadratic)? == 0 {
                linear_roots(&linear, &constant)?
            } else if sign(&constant)? == 0 {
                // The other polynomial root is exactly the excluded endpoint r=0.
                linear_roots(&quadratic, &linear)?
            } else {
                let discriminant = add(
                    &square(&linear),
                    &neg(&scale(&mul(&quadratic, &constant), 4)),
                );
                match sign(&discriminant)? {
                    -1 => Vec::new(),
                    0 => {
                        let root = divide(&neg(&linear), &scale(&quadratic, 2))?;
                        if root[1] <= 0 {
                            Vec::new()
                        } else {
                            return Err(uncertain(
                                "a positive tangent or unresolved repeated root is not a simple cut",
                            ));
                        }
                    }
                    _ => {
                        let radical = sqrt(&discriminant)?;
                        [
                            add(&neg(&linear), &neg(&radical)),
                            add(&neg(&linear), &radical),
                        ]
                        .iter()
                        .map(|numerator| divide(numerator, &scale(&quadratic, 2)))
                        .collect::<eyre::Result<Vec<_>>>()?
                    }
                }
            }
        } else {
            linear_roots(&slope, &shift)?
        };
        candidates.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());
        if candidates.windows(2).any(|pair| pair[0][1] >= pair[1][0]) {
            return Err(uncertain("distinct polynomial roots are not separated"));
        }

        let materialize = |interval: &[Float; 2]| -> eyre::Result<F<T>> {
            let midpoint = Float::with_val(interval[0].prec(), &interval[0] + &interval[1]) / 2;
            F::from_arb(&ArbPrec::from(midpoint))
        };
        let mut roots = Vec::new();
        for (index, candidate) in candidates.iter().enumerate() {
            if candidate[1] <= 0 {
                continue;
            }
            if candidate[0] <= 0 {
                return Err(uncertain(
                    "a root cannot be separated from the radial endpoint",
                ));
            }
            if let Some(curved_index) = curved {
                // Squaring adds the opposite square-root sheet. Reject it by
                // the unsquared sign constraint, never by a small residual test.
                let constraint =
                    scale(&add(&mul(&slope, candidate), &shift), terms[curved_index].0);
                match sign(&constraint)? {
                    1 => continue,
                    -1 => (),
                    _ => return Err(uncertain("a root lies at a vanishing cut energy")),
                }
            }
            let crossing = sign(&self.evaluate_interval(candidate, true)?.1)?;
            if crossing == 0 {
                return Err(uncertain("a positive tangent root is not a simple cut"));
            }
            let radius = materialize(candidate)?;
            let four = zero.from_i64(4);
            let mut width = &radius / &four;
            if let Some(previous) = index.checked_sub(1).map(|i| &candidates[i]) {
                let separation = (&radius - materialize(previous)?) / &four;
                if separation < width {
                    width = separation;
                }
            }
            if let Some(next) = candidates.get(index + 1) {
                let separation = (materialize(next)? - &radius) / &four;
                if separation < width {
                    width = separation;
                }
            }
            let lower = &radius - &width;
            let upper = &radius + &width;
            if lower <= zero || lower >= radius || upper <= radius {
                return Err(uncertain("the positive root bracket is not representable"));
            }
            if sign(&self.evaluate_interval(&native(&lower), false)?.0)? != -crossing
                || sign(&self.evaluate_interval(&native(&upper), false)?.0)? != crossing
            {
                return Err(uncertain("the isolated root bracket is not sign-certified"));
            }
            let branch_identity =
                RadialRootIdentity::new(format!("{identity} signed root {index}"));
            let result = self
                .refine_root(
                    &lower,
                    &upper,
                    crossing,
                    0,
                    e_cm,
                    diagnostics,
                    &branch_identity,
                )
                .map_err(|error| {
                    uncertain(&format!("isolated root refinement failed: {error:?}"))
                })?;
            roots.push(result);
        }
        Ok(roots)
    }

    /// Original eta jet at the same represented coefficients used by the root.
    /// The caller retains the existing dual shape and factorial convention.
    #[inline]
    pub(crate) fn evaluate_dual(&self, radius: &HyperDual<F<T>>) -> HyperDual<F<T>> {
        let zero = self.shift.zero();
        let energy_sum = self
            .energies
            .iter()
            .map(
                |RadialEnergy {
                     coefficient,
                     velocity,
                     offset,
                     mass,
                     ..
                 }| {
                    let energy =
                        if velocity.px == zero && velocity.py == zero && velocity.pz == zero {
                            // A constant zero energy has a constant jet; asking the
                            // dual square root to differentiate it would produce 0/0.
                            new_constant(radius, &(offset.norm_squared() + mass * mass).sqrt())
                        } else {
                            let momentum = velocity.map_ref(&|v| new_constant(radius, v) * radius)
                                + offset.map_ref(&|b| new_constant(radius, b));
                            (momentum.norm_squared() + mass * mass).sqrt()
                        };
                    if *coefficient == 1 {
                        energy
                    } else {
                        let coefficient = self.shift.from_i64(*coefficient);
                        energy * new_constant(radius, &coefficient)
                    }
                },
            )
            .reduce(|sum, energy| sum + energy)
            .unwrap_or_else(|| radius.zero());
        energy_sum + new_constant(radius, &self.shift)
    }
}

impl<T: FloatLike> Rotatable for EnergySurfaceRay<T> {
    fn rotate(&self, rotation: &Rotation) -> Self {
        Self {
            energies: self
                .energies
                .iter()
                .map(|energy| RadialEnergy {
                    edge: energy.edge,
                    coefficient: energy.coefficient,
                    velocity: energy.velocity.rotate(rotation),
                    offset: energy.offset.rotate(rotation),
                    mass: energy.mass.clone(),
                })
                .collect(),
            shift: self.shift.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::QuadFloat;
    use itertools::Itertools;

    #[test]
    fn canonical_ray_embedding_preserves_ordered_occurrences_and_native_coefficients() {
        use crate::utils::SamplingFloat;

        fn check<T: FloatLike>(marker: F<T>) {
            let one = marker.one();
            let velocity = ThreeMomentum::new(marker.clone(), -marker.clone(), one.zero());
            let offset = ThreeMomentum::new(
                one.from_usize(2),
                &marker / one.from_usize(3),
                F(T::from_f64_exact_binary(-0.0)),
            );
            let repeated = RadialEnergy {
                edge: EdgeIndex(4),
                coefficient: 1,
                velocity: velocity.clone(),
                offset: offset.clone(),
                mass: one.clone(),
            };
            let ray = EnergySurfaceRay {
                energies: vec![
                    repeated.clone(),
                    RadialEnergy {
                        edge: EdgeIndex(7),
                        coefficient: 1,
                        offset,
                        velocity,
                        mass: marker.clone(),
                    },
                    repeated,
                ],
                shift: -marker.clone(),
            };
            let canonical = ray.to_arb_exact().unwrap();
            assert_eq!(
                canonical
                    .energies
                    .iter()
                    .map(|entry| entry.edge)
                    .collect_vec(),
                [EdgeIndex(4), EdgeIndex(7), EdgeIndex(4)]
            );
            assert_eq!(canonical.energies[0], canonical.energies[2]);
            assert_eq!(
                canonical.energies[0].velocity.px,
                marker.to_arb_exact().unwrap()
            );
            // Preserve the native scalar's zero convention; a Quad sum need
            // not retain the sign of the binary64 token used to construct it.
            assert_eq!(
                canonical.energies[0]
                    .offset
                    .pz
                    .0
                    .mpfr_enclosure(1000)
                    .0
                    .is_sign_negative(),
                ray.energies[0].offset.pz.0.into_f64().is_sign_negative()
            );
            let roundtrip = canonical.materialize::<T>().unwrap();
            assert_eq!(roundtrip.energies, ray.energies);
            assert_eq!(roundtrip.shift, ray.shift);
            for field in 0..4 {
                let mut invalid = ray.clone();
                let nonfinite = F(T::from_f64_exact_binary(f64::INFINITY));
                match field {
                    0 => invalid.energies[0].velocity.px = nonfinite,
                    1 => invalid.energies[0].offset.py = nonfinite,
                    2 => invalid.energies[0].mass = nonfinite,
                    _ => invalid.shift = nonfinite,
                }
                assert!(invalid.to_arb_exact().is_err());
            }
        }
        let one = F::<QuadFloat>::default().one();
        check(one + one.from_i64(2).powi(-80));
        let one = F::<SamplingFloat>::default().one();
        check(&one + one.from_i64(2).powi(-180) + one.from_i64(2).powi(-240));
    }

    #[test]
    fn hosted_normal_alignment_enforces_directed_half_budgets() {
        use rug::Float;
        let point = |value: i32| {
            let x = Float::with_val(2048, value);
            [x.clone(), x]
        };
        let canonical = [point(3), point(4)];
        let half_budget: Float = Float::with_val(2048, 5) / 16;
        let shifted: Float = Float::with_val(2048, 3) + &half_budget;
        let native = [[shifted.clone(), shifted], point(4)];
        let host = [-half_budget.clone(), -half_budget.clone()];
        // R=5 and epsilon=1/8 allocate exactly 5/16 to each check.
        EnergySurfaceRay::<f64>::verify_normal_alignment(
            canonical.clone(),
            native.clone(),
            host.clone(),
            0.125,
        )
        .unwrap();
        let beyond: Float = Float::with_val(2048, 1) >> 1000;
        let mut bad_normal = native;
        bad_normal[0][1] += &beyond;
        let mut bad_host = host;
        bad_host[0] -= beyond;
        for (native, host) in [(bad_normal, point(0)), (canonical.clone(), bad_host)] {
            assert!(matches!(
                EnergySurfaceRay::<f64>::verify_normal_alignment(
                    canonical.clone(),
                    native,
                    host,
                    0.125
                )
                .unwrap_err()
                .downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { .. })
            ));
        }
        for canonical in [
            [point(0), point(0)],
            [
                [Float::with_val(2048, -1), Float::with_val(2048, 1)],
                point(0),
            ],
        ] {
            assert!(matches!(
                EnergySurfaceRay::<f64>::verify_normal_alignment(
                    canonical.clone(),
                    canonical,
                    point(0),
                    0.125
                )
                .unwrap_err()
                .downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { .. })
            ));
        }
        assert!(
            EnergySurfaceRay::<f64>::verify_normal_alignment(
                canonical.clone(),
                canonical.clone(),
                point(0),
                0.0
            )
            .is_err()
        );
        assert!(matches!(
            EnergySurfaceRay::<f64>::verify_normal_alignment(
                canonical.clone(),
                canonical,
                std::array::from_fn(|_| Float::with_val(2048, rug::float::Special::Infinity)),
                0.125
            )
            .unwrap_err()
            .downcast_ref::<SamplingEvaluationError>(),
            Some(SamplingEvaluationError::Unrepresentable { .. })
        ));
    }

    #[test]
    fn retained_lu_candidate_requires_directed_native_agreement() {
        use crate::{
            integrands::process::sampling::maps::SamplingEvaluationError,
            momentum::{Rotatable, Rotation, RotationMethod},
        };
        fn check<T: FloatLike>() {
            let one = F::<T>::default().one();
            let f = |n| one.from_usize(n);
            let vector = |x| ThreeMomentum::new(x, one.zero(), one.zero());
            // eta(t)=sqrt((3t+1)^2+9)-6, whose positive root is
            // (sqrt(27)-1)/3. The independent radical is not a Newton result.
            let ray = EnergySurfaceRay {
                energies: vec![RadialEnergy {
                    edge: EdgeIndex(1),
                    coefficient: 1,
                    velocity: vector(f(3)),
                    offset: vector(one.clone()),
                    mass: f(3),
                }],
                shift: -f(6),
            };
            let solution = ray
                .solve_lu_cut(
                    &f(6),
                    &mut RadialRootDiagnostics::default(),
                    &RadialRootIdentity::new("directed native host fixture".into()),
                )
                .unwrap();
            let expected = (f(27).sqrt() - &one) / f(3);
            assert!((&solution.solution - expected).abs() < one.epsilon() * f(100));
            let budget = &one / f(1_000_000);
            ray.verify_lu_candidate(&ray, &solution, 6, &budget)
                .unwrap();
            let rotation = Rotation::new(RotationMethod::EulerAngles(0.31, -0.17, 0.23));
            ray.rotate(&rotation)
                .verify_lu_candidate(&ray, &solution, 6, &budget)
                .unwrap();

            // The certificate encloses represented coefficients directly. A
            // shallow positive slope and tiny nonzero masses remain resolvable
            // even where native squaring would underflow; no Newton claim is
            // made for that underflowing native evaluator in this direct check.
            let tiny = &one / f(10).powi(200);
            let mut shallow = ray.clone();
            for energy in &mut shallow.energies {
                energy.velocity = &energy.velocity * &tiny;
                energy.offset = &energy.offset * &tiny;
                energy.mass *= &tiny;
            }
            shallow.shift *= &tiny;
            let mut shallow_solution = solution.clone();
            shallow_solution.derivative_at_solution *= &tiny;
            shallow_solution.error_of_function *= &tiny;
            shallow
                .verify_lu_candidate(&shallow, &shallow_solution, 6, &budget)
                .unwrap();
            let mut tiny_energy = ray.clone();
            tiny_energy.energies.push(RadialEnergy {
                edge: EdgeIndex(2),
                coefficient: 1,
                velocity: vector(one.zero()),
                offset: vector(one.zero()),
                mass: tiny,
            });
            tiny_energy
                .verify_lu_candidate(&tiny_energy, &solution, 6, &budget)
                .unwrap();
            tiny_energy.energies[1].mass = one.zero();
            tiny_energy
                .verify_lu_candidate(&tiny_energy, &solution, 6, &budget)
                .unwrap();
            // A varying energy that vanishes at the candidate is different:
            // its cusp must still fail the derivative certificate.
            tiny_energy.energies[1].velocity = vector(one.clone());
            tiny_energy.energies[1].offset = vector(-solution.solution.clone());
            let error = tiny_energy
                .verify_lu_candidate(&tiny_energy, &solution, 6, &budget)
                .unwrap_err();
            assert!(
                matches!(
                    error.downcast_ref::<SamplingEvaluationError>(),
                    Some(SamplingEvaluationError::UncertainGeometry { .. })
                ),
                "{error:?}"
            );
            assert!(
                error
                    .to_string()
                    .contains("energy lower bound touches zero"),
                "{error:?}"
            );

            let mut wrong = ray.clone();
            wrong.shift -= &one / f(100);
            let error = ray
                .verify_lu_candidate(&wrong, &solution, 6, &budget)
                .unwrap_err();
            assert!(
                matches!(
                    error.downcast_ref::<SamplingEvaluationError>(),
                    Some(SamplingEvaluationError::UncertainGeometry { .. })
                ),
                "{error:?}"
            );
            let mut wrong_derivative = solution.clone();
            wrong_derivative.derivative_at_solution *= f(11) / f(10);
            assert!(
                ray.verify_lu_candidate(&ray, &wrong_derivative, 6, &budget)
                    .is_err()
            );
            wrong = ray.clone();
            wrong.energies[0].edge = EdgeIndex(2);
            assert!(
                ray.verify_lu_candidate(&wrong, &solution, 6, &budget)
                    .unwrap_err()
                    .to_string()
                    .contains("ordered energy occurrences")
            );
            wrong = ray.clone();
            wrong.energies[0].mass = one.zero() / one.zero();
            assert!(matches!(
                ray.verify_lu_candidate(&wrong, &solution, 6, &budget)
                    .unwrap_err()
                    .downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::Unrepresentable { .. })
            ));

            // One rounded residual value cannot certify this much smaller
            // budget. Native Quad/Arb resolve the same irrational root, while
            // binary64 must ask for original-source precision rescue.
            let tight = &one / f(10).powi(25);
            let tight_result = ray.verify_lu_candidate(&ray, &solution, 6, &tight);
            let large = f(10).powi(20);
            let cancellation = (&large + &one) - &large;
            let mut completed = ray.clone();
            completed.energies[0].velocity.px += &cancellation - &one;
            let cancellation_result = ray.verify_lu_candidate(&completed, &solution, 6, &budget);
            if T::sampling_precision() == crate::utils::SamplingPrecision::Double {
                assert_eq!(cancellation, one.zero());
                for result in [tight_result, cancellation_result] {
                    assert!(matches!(
                        result
                            .unwrap_err()
                            .downcast_ref::<SamplingEvaluationError>(),
                        Some(SamplingEvaluationError::UncertainGeometry { .. })
                    ));
                }
            } else {
                assert_eq!(cancellation, one);
                tight_result.unwrap();
                cancellation_result.unwrap();
            }
        }
        check::<f64>();
        check::<QuadFloat>();
        check::<ArbPrec>();
    }

    #[test]
    fn signed_ray_uses_one_equation_for_values_enclosures_and_jets() {
        use crate::utils::hyperdual_utils::{extract_t_derivatives, simple_n_deriv_shape};

        fn check<T: FloatLike>() {
            let one = F::<T>::default().one();
            let f = |n| one.from_i64(n);
            let vector = |x| ThreeMomentum::new(x, one.zero(), one.zero());
            let repeated = RadialEnergy {
                edge: EdgeIndex(4),
                coefficient: 1,
                velocity: vector(f(3)),
                offset: vector(one.clone()),
                mass: f(3),
            };
            // eta(r) = 2 sqrt((3r+1)^2+9) - sqrt((4r-1)^2+16) - 2.
            let ray = EnergySurfaceRay {
                energies: vec![
                    repeated.clone(),
                    RadialEnergy {
                        edge: EdgeIndex(7),
                        coefficient: -1,
                        velocity: vector(f(4)),
                        offset: vector(-one.clone()),
                        mass: f(4),
                    },
                    repeated,
                ],
                shift: -f(2),
            };
            let expected = [f(3), f(12) / f(5), -f(94) / f(125), f(3384) / f(3125)];
            let (value, derivative) = ray.evaluate(&one);
            assert_eq!(value, expected[0]);
            assert!((derivative - &expected[1]).abs() < one.epsilon() * f(64));
            let interval = EnergySurfaceRay::<T>::enclosed_native(&one);
            let (value, derivative) = ray.evaluate_interval(&interval, true).unwrap();
            assert_eq!(value[0], 3);
            assert_eq!(value[1], 3);
            let exact_derivative = Float::with_val(2048, rug::Rational::from((12, 5)));
            assert!(derivative[0] <= exact_derivative && exact_derivative <= derivative[1]);
            for order in 1..=3 {
                let radius =
                    HyperDual::<F<T>>::new(simple_n_deriv_shape(order)).variable(0, one.clone());
                for (actual, expected) in extract_t_derivatives(ray.evaluate_dual(&radius))
                    .iter()
                    .zip(&expected)
                {
                    assert!(
                        (actual - expected).abs()
                            < one.epsilon() * f(256) * (one.clone() + expected.abs())
                    );
                }
            }
            let canonical = ray.to_arb_exact().unwrap();
            assert_eq!(
                canonical
                    .energies
                    .iter()
                    .map(|energy| (energy.edge, energy.coefficient))
                    .collect_vec(),
                [(EdgeIndex(4), 1), (EdgeIndex(7), -1), (EdgeIndex(4), 1)]
            );
            let materialized = canonical.materialize::<T>().unwrap();
            assert_eq!(materialized.energies, ray.energies);
            assert_eq!(materialized.shift, ray.shift);
        }
        check::<f64>();
        check::<QuadFloat>();
        check::<ArbPrec>();
    }

    #[test]
    fn signed_root_isolation_retains_occurrences_and_shared_jet_slopes() {
        use crate::utils::hyperdual_utils::{extract_t_derivatives, simple_n_deriv_shape};
        let vector = |x| ThreeMomentum::new(F(x), F(0.0), F(0.0));
        let repeated = RadialEnergy {
            edge: EdgeIndex(4),
            coefficient: 1,
            velocity: vector(1.0),
            offset: vector(0.0),
            mass: F(1.0),
        };
        let ray = EnergySurfaceRay {
            energies: vec![
                repeated.clone(),
                RadialEnergy {
                    edge: EdgeIndex(7),
                    coefficient: -1,
                    velocity: vector(1.9),
                    offset: vector(0.0),
                    mass: F(0.0),
                },
                repeated,
            ],
            shift: F(-1.0),
        };
        let retained = ray.energies.clone();
        let roots = ray
            .positive_radius_roots(
                &F(2.0),
                &mut RadialRootDiagnostics::default(),
                &RadialRootIdentity::new("shared signed ray".into()),
            )
            .unwrap();
        assert_eq!(roots.len(), 2);
        assert!(roots[0].derivative_at_solution < F(0.0));
        assert!(roots[1].derivative_at_solution > F(0.0));
        assert_eq!(ray.energies, retained);
        for root in &roots {
            let (residual, slope) = ray.evaluate(&root.solution);
            assert!(residual.abs() < F(1.0e-13));
            assert_eq!(slope, root.derivative_at_solution);
            let radius =
                HyperDual::<F<f64>>::new(simple_n_deriv_shape(1)).variable(0, root.solution);
            let jet = extract_t_derivatives(ray.evaluate_dual(&radius));
            assert!((jet[1] - slope).abs() < F(1.0e-14));
        }
        assert!(
            ray.verify_lu_candidate(&ray, &roots[0], 6, &F(1.0e-6))
                .unwrap_err()
                .to_string()
                .contains("positive E-surface rays")
        );
    }

    #[test]
    fn constant_zero_energy_preserves_signed_roots_and_jets() {
        use crate::utils::hyperdual_utils::{extract_t_derivatives, simple_n_deriv_shape};
        let vector = |x| ThreeMomentum::new(F(x), F(0.0), F(0.0));
        let ray = EnergySurfaceRay {
            energies: vec![RadialEnergy {
                edge: EdgeIndex(0),
                coefficient: -1,
                velocity: vector(1.0),
                offset: vector(0.0),
                mass: F(0.0),
            }],
            shift: F(2.0),
        };
        let constant = RadialEnergy {
            edge: EdgeIndex(1),
            coefficient: 1,
            velocity: vector(0.0),
            offset: vector(0.0),
            mass: F(0.0),
        };
        let mut with_constant = ray.clone();
        with_constant.energies.push(constant.clone());
        let identity = RadialRootIdentity::new("constant zero energy".into());
        let roots = ray
            .positive_radius_roots(&F(2.0), &mut RadialRootDiagnostics::default(), &identity)
            .unwrap();
        let with_roots = with_constant
            .positive_radius_roots(&F(2.0), &mut RadialRootDiagnostics::default(), &identity)
            .unwrap();
        assert_eq!(roots.len(), 1);
        assert_eq!(with_roots.len(), 1);
        assert_eq!(roots[0].solution, F(2.0));
        assert_eq!(with_roots[0].solution, roots[0].solution);
        assert_eq!(
            with_roots[0].derivative_at_solution,
            roots[0].derivative_at_solution
        );
        assert_eq!(with_constant.energies.len(), 2);
        assert_eq!(with_constant.evaluate(&F(2.0)), ray.evaluate(&F(2.0)));
        let interval = EnergySurfaceRay::<f64>::enclosed_native(&F(2.0));
        assert_eq!(
            with_constant.evaluate_interval(&interval, true).unwrap(),
            ray.evaluate_interval(&interval, true).unwrap()
        );
        let radius = HyperDual::<F<f64>>::new(simple_n_deriv_shape(3)).variable(0, F(2.0));
        assert_eq!(
            extract_t_derivatives(with_constant.evaluate_dual(&radius)),
            extract_t_derivatives(ray.evaluate_dual(&radius))
        );

        // A nonzero source velocity remains nonconstant even when its native
        // squared norm underflows. Its meaningful slope must not become zero.
        let mut underflowing = constant;
        underflowing.velocity = vector(1.0e-200);
        let underflowing = EnergySurfaceRay {
            energies: vec![underflowing],
            shift: F(0.0),
        };
        assert!(underflowing.evaluate(&F(1.0)).1.is_nan());
        let interval = EnergySurfaceRay::<f64>::enclosed_native(&F(1.0));
        let (_, slope) = underflowing.evaluate_interval(&interval, true).unwrap();
        assert!(slope[0] > 0);
    }
}
