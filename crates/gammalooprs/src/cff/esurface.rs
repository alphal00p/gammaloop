use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::HedgeGraph;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Flow, HedgePair};
use linnet::half_edge::subgraph::OrientedCut;
use ref_ops::RefNeg;
use rug::{Float, float::Round};
use serde::{Deserialize, Serialize};

use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::dual::HyperDual;
use symbolica::domains::float::{FloatLike as SymFloatLike, Real};
use symbolica::id::Replacement;
use symbolica::{function, parse};
use tracing::debug;
use typed_index_collections::TiVec;

use crate::cff::VertexSet;

use crate::cff::expression::{
    CFFExpression, OrientationID, RaisedEsurfaceDataView, RaisedEsurfaceGroupView,
};
pub use crate::cff::surface::EsurfaceID;
use crate::graph::{Graph, GraphGroupPosition, LmbIndex, LoopMomentumBasis};
use crate::{GammaLoopContext, define_index};

use crate::integrands::process::sampling_maps::SamplingEvaluationError;
use crate::integrands::process::{GenericEvaluator, ImplicitSurfaceRadialMap};
use crate::momentum::sample::{
    ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex, LoopMomenta, SubspaceData,
};
use crate::momentum::{Rotatable, Rotation, ThreeMomentum};
use crate::processes::CrossSectionCut;
use crate::utils::hyperdual_utils::new_constant;
use crate::utils::newton_solver::{
    NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity, SafeguardedNewtonError,
};
use crate::utils::{
    ArbPrec, DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, ESURFACE_SHIFT_THRESHOLD, F, FloatLike, GS,
    Length, compute_loop_part, compute_loop_part_subspace, compute_shift_part,
    compute_shift_part_subspace, compute_t_part_of_shift_part, cut_energy,
    external_energy_atom_from_index, ose_atom_from_index,
};
use crate::uv::uv_graph::UVE;
use color_eyre::Result;

use super::generation::ShiftRewrite;

/// Core esurface struct
#[derive(Serialize, Deserialize, Debug, Clone, bincode::Encode, bincode::Decode)]
pub struct Esurface {
    pub energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
    //#[bincode(with_serde)]
    //pub subspace_graph: InternalSubGraph,
}

// Edge identity, radial velocity, constant spatial offset, and mass.
type EsurfaceRayEnergy<T> = (EdgeIndex, ThreeMomentum<F<T>>, ThreeMomentum<F<T>>, F<T>);
// Ephemeral endpoint tuples for the same affine energy fold, with no stored ray.
type EnclosedRayEnergy = ([[Float; 2]; 3], [[Float; 2]; 3], [Float; 2]);

/// One represented affine energy equation, retaining its ordered edge occurrences.
/// This is native geometry, not a cross-precision or selected-host transport record.
#[derive(Debug, Clone)]
pub(crate) struct EsurfaceRay<T: FloatLike> {
    // A Vec deliberately preserves repeated energy occurrences; an edge-keyed
    // map would not. The tuple's fields are documented on EsurfaceRayEnergy.
    energies: Vec<EsurfaceRayEnergy<T>>,
    shift: F<T>,
}

impl EsurfaceRay<ArbPrec> {
    /// Retain the canonical ordered equation while materializing its native
    /// coefficients. Physical adoption still authenticates and checks the ray.
    pub(crate) fn materialize<T: FloatLike>(&self) -> Result<EsurfaceRay<T>> {
        let spatial = |p: &ThreeMomentum<F<ArbPrec>>| -> Result<ThreeMomentum<F<T>>> {
            Ok(ThreeMomentum::new(
                F::from_arb(&p.px.0)?,
                F::from_arb(&p.py.0)?,
                F::from_arb(&p.pz.0)?,
            ))
        };
        Ok(EsurfaceRay {
            energies: self
                .energies
                .iter()
                .map(|(edge, velocity, offset, mass)| {
                    Ok((
                        *edge,
                        spatial(velocity)?,
                        spatial(offset)?,
                        F::from_arb(&mass.0)?,
                    ))
                })
                .collect::<Result<_>>()?,
            shift: F::from_arb(&self.shift.0)?,
        })
    }
}

impl<T: FloatLike> EsurfaceRay<T> {
    /// Embed the represented equation exactly, preserving each ordered energy
    /// occurrence. No routing, coefficient reconstruction or root solve occurs.
    pub(crate) fn to_arb_exact(&self) -> Result<EsurfaceRay<ArbPrec>> {
        let spatial = |p: &ThreeMomentum<F<T>>| -> Result<ThreeMomentum<F<ArbPrec>>> {
            Ok(ThreeMomentum::new(
                p.px.to_arb_exact()?,
                p.py.to_arb_exact()?,
                p.pz.to_arb_exact()?,
            ))
        };
        Ok(EsurfaceRay {
            energies: self
                .energies
                .iter()
                .map(|(edge, velocity, offset, mass)| {
                    Ok((
                        *edge,
                        spatial(velocity)?,
                        spatial(offset)?,
                        mass.to_arb_exact()?,
                    ))
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

    fn enclosed_native(value: &F<T>) -> [Float; 2] {
        let (lo, hi) = value.0.mpfr_enclosure(Self::ENCLOSURE_PRECISION);
        [lo, hi]
    }

    fn enclosed_sum(terms: impl Iterator<Item = (i64, [Float; 2])>) -> [Float; 2] {
        let zero = Float::with_val(Self::ENCLOSURE_PRECISION, 0);
        terms.fold([zero.clone(), zero], |sum, (coefficient, value)| {
            let coefficient = Float::with_val(Self::ENCLOSURE_PRECISION, coefficient);
            Self::enclosed_add(
                &sum,
                &Self::enclosed_mul(&[coefficient.clone(), coefficient], &value),
            )
        })
    }

    fn evaluate_enclosed(
        energies: impl Iterator<Item = EnclosedRayEnergy>,
        shift: [Float; 2],
        radius: &[Float; 2],
        with_derivative: bool,
    ) -> Result<([Float; 2], [Float; 2])> {
        let zero = Float::with_val(Self::ENCLOSURE_PRECISION, 0);
        let point = |x: &Float| [x.clone(), x.clone()];
        let uncertain = |detail: &str| SamplingEvaluationError::UncertainGeometry {
            detail: format!("LU host adoption: {detail}"),
        };
        let mut value = shift;
        let mut slope = point(&zero);
        for (velocity, offset, mass) in energies {
            let mut norm = Self::enclosed_square(&mass);
            let mut dot = point(&zero);
            for (v, b) in velocity.into_iter().zip(offset) {
                let q = Self::enclosed_add(&Self::enclosed_mul(&v, radius), &b);
                norm = Self::enclosed_add(&norm, &Self::enclosed_square(&q));
                dot = Self::enclosed_add(&dot, &Self::enclosed_mul(&v, &q));
            }
            if norm.iter().chain(&dot).any(|x| !x.is_finite()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "LU host certificate",
                    detail: "nonfinite directed energy or numerator".into(),
                }
                .into());
            }
            let energy = [
                Float::with_val_round(Self::ENCLOSURE_PRECISION, norm[0].sqrt_ref(), Round::Down).0,
                Float::with_val_round(Self::ENCLOSURE_PRECISION, norm[1].sqrt_ref(), Round::Up).0,
            ];
            value = Self::enclosed_add(&value, &energy);
            if with_derivative {
                if energy[0] <= 0 {
                    return Err(uncertain(
                        "energy lower bound touches zero on the fixed candidate interval",
                    )
                    .into());
                }
                slope = Self::enclosed_add(&slope, &Self::enclosed_divide(&dot, &energy));
            }
        }
        if value.iter().chain(&slope).any(|x| !x.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "LU host certificate",
                detail: "nonfinite directed residual or derivative".into(),
            }
            .into());
        }
        Ok((value, slope))
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
            .map(|(_, velocity, offset, mass)| {
                let momentum = velocity * radius + offset;
                let energy = (momentum.norm_squared() + mass * mass).sqrt();
                // At a massless endpoint E(r)=r|v| the radial right derivative is
                // |v|, whereas the two-sided formula q.v/E would evaluate 0/0.
                // Share this convention with sampling charts, using exact source
                // zeros: a computed E=0 can instead be numerical underflow.
                let derivative = if radius == &zero
                    && mass == &zero
                    && momentum.px == zero
                    && momentum.py == zero
                    && momentum.pz == zero
                {
                    velocity.norm_squared().sqrt()
                } else {
                    momentum * velocity / &energy
                };
                (derivative, energy)
            })
            .fold(
                (zero.clone(), zero.clone()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );
        (energy_sum + &self.shift, derivative)
    }

    /// The sole LU seed/Newton policy, also used by native partial host rays.
    pub(crate) fn solve_lu_cut(
        &self,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> std::result::Result<NewtonIterationResult<T>, SafeguardedNewtonError<T>> {
        let guess = Esurface::radius_guess_from_terms(
            &self.shift,
            self.energies
                .iter()
                .map(|(_, v, b, _)| (v.norm_squared(), v.clone() * b)),
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
        let residual_tolerance = guess.from_i64((4 * self.energies.len() + 1) as i64);
        diagnostics.solve(
            identity,
            &guess.zero(),
            &guess,
            |t| self.evaluate(t),
            &residual_tolerance,
            2000,
            64,
            e_cm,
        )
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
        let precision = Self::ENCLOSURE_PRECISION;
        let uncertain = |detail: &str| SamplingEvaluationError::UncertainGeometry {
            detail: format!("LU host adoption: {detail}"),
        };
        if self
            .energies
            .iter()
            .map(|entry| entry.0)
            .ne(completed.energies.iter().map(|entry| entry.0))
        {
            return Err(eyre!(
                "LU host adoption has incompatible ordered energy occurrences"
            ));
        }
        for ray in [self, completed] {
            if !ray.shift.0.is_finite()
                || ray.energies.iter().any(|(_, v, b, m)| {
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
        let evaluate = |ray: &Self, radius: &[Float; 2], with_derivative: bool| {
            Self::evaluate_enclosed(
                ray.energies.iter().map(|(_, velocity, offset, mass)| {
                    (
                        [&velocity.px, &velocity.py, &velocity.pz].map(native),
                        [&offset.px, &offset.py, &offset.pz].map(native),
                        native(mass),
                    )
                }),
                native(&ray.shift),
                radius,
                with_derivative,
            )
        };
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
            if evaluate(ray, &point(&zero), false)?.0[1] >= 0
                || evaluate(ray, &point(&interval[0]), false)?.0[1] >= 0
                || evaluate(ray, &point(&interval[1]), false)?.0[0] <= 0
            {
                return Err(uncertain("origin/interior/exterior signs are not certified on the fixed candidate interval").into());
            }
            let (_, slope) = evaluate(ray, &interval, true)?;
            if slope[0] <= 0 {
                return Err(uncertain(
                    "radial derivative is not certified positive throughout the interval",
                )
                .into());
            }
            let residual = evaluate(ray, &t0, false)?.0;
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

    /// Original eta jet at the same represented coefficients used by the root.
    /// The caller retains the existing dual shape and factorial convention.
    #[inline]
    pub(crate) fn evaluate_dual(&self, radius: &HyperDual<F<T>>) -> HyperDual<F<T>> {
        let energy_sum = self
            .energies
            .iter()
            .map(|(_, velocity, offset, mass)| {
                let momentum = velocity.map_ref(&|v| new_constant(radius, v) * radius)
                    + offset.map_ref(&|b| new_constant(radius, b));
                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|sum, energy| sum + energy)
            .unwrap_or_else(|| radius.zero());
        energy_sum + new_constant(radius, &self.shift)
    }
}

impl<T: FloatLike> Rotatable for EsurfaceRay<T> {
    fn rotate(&self, rotation: &Rotation) -> Self {
        Self {
            energies: self
                .energies
                .iter()
                .map(|(edge, v, b, mass)| {
                    (*edge, v.rotate(rotation), b.rotate(rotation), mass.clone())
                })
                .collect(),
            shift: self.shift.clone(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NonExistingEsurfaceReason {
    NoExternalShift,
    ShiftNotNegative,
    NoRadialDependence,
    NoRealZero,
}

#[derive(Debug, Clone)]
pub(crate) enum EsurfaceExistence<T: FloatLike> {
    NonExisting {
        normalized_margin: Option<F<T>>,
        reason: NonExistingEsurfaceReason,
    },
    Pinched {
        normalized_margin: F<T>,
    },
    Existing {
        normalized_margin: F<T>,
    },
}

/// Pointwise existence classification exposed without the internal diagnostic margin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EsurfaceExistenceStatus {
    NonExisting,
    Pinched,
    Existing,
}

impl<T: FloatLike> EsurfaceExistence<T> {
    fn status(&self) -> EsurfaceExistenceStatus {
        match self {
            Self::NonExisting { .. } => EsurfaceExistenceStatus::NonExisting,
            Self::Pinched { .. } => EsurfaceExistenceStatus::Pinched,
            Self::Existing { .. } => EsurfaceExistenceStatus::Existing,
        }
    }

    pub(crate) fn is_existing(&self) -> bool {
        matches!(self, Self::Existing { .. })
    }

    pub(crate) fn normalized_margin(&self) -> Option<&F<T>> {
        match self {
            Self::NonExisting {
                normalized_margin, ..
            } => normalized_margin.as_ref(),
            Self::Pinched { normalized_margin } | Self::Existing { normalized_margin } => {
                Some(normalized_margin)
            }
        }
    }

    pub(crate) fn label(&self) -> &'static str {
        match self {
            Self::NonExisting { .. } => "non_existing",
            Self::Pinched { .. } => "pinched",
            Self::Existing { .. } => "existing",
        }
    }

    pub(crate) fn non_existing_reason(&self) -> Option<NonExistingEsurfaceReason> {
        match self {
            Self::NonExisting { reason, .. } => Some(*reason),
            Self::Pinched { .. } | Self::Existing { .. } => None,
        }
    }
}

pub(crate) fn esurface_value_is_strictly_inside<T: FloatLike>(value: &F<T>, e_cm: &F<T>) -> bool {
    let interior_tolerance = value.epsilon() * value.from_i64(8) * e_cm;
    !value.is_nan() && !value.is_infinite() && value < &(-interior_tolerance)
}

impl PartialEq for Esurface {
    fn eq(&self, other: &Self) -> bool {
        self.energies == other.energies && self.external_shift == other.external_shift
    }
}

impl Eq for Esurface {}

impl Esurface {
    /// Match two energy sums on one active spatial momentum. Their shared
    /// energy is identified by a single sign of the full internal routing,
    /// exact external spatial shift and mass expression, never by edge identity
    /// or a rounded-momentum tolerance.
    /// The returned common route has active coefficient +1: the kernel uses
    /// x=L+c0, so the existing affine embedding must return L=x-c0.
    // Physical binders must transport one prerequisite-only disk policy
    // through native precision retries alongside this reconstructed geometry.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sampling_joint_geometry_in_subspace<T: FloatLike>(
        &self,
        other: &Self,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        complement: &[LoopIndex],
    ) -> Result<(
        crate::integrands::process::sampling_joint::SharedEnergyJointGeometryEvaluator<T>,
        crate::momentum::signature::LoopExtSignature,
    )> {
        use crate::integrands::process::{
            SharedEnergyJointGeometry, sampling_maps::SamplingEvaluationError,
        };
        use crate::momentum::{SignOrZero, signature::LoopExtSignature};
        use std::sync::Arc;

        let lmb = subspace.get_lmb(all_lmbs);
        let active = subspace.iter_lmb_indices().collect_vec();
        let [active] = active.as_slice() else {
            return Err(eyre!(
                "joint energy sampling requires exactly one active loop cycle"
            ));
        };
        let active = *active;
        let mut covered = complement.iter().copied().chain([active]).collect_vec();
        covered.sort();
        if covered.iter().any(|index| index.0 >= lmb.loop_edges.len())
            || covered.windows(2).any(|pair| pair[0] == pair[1])
        {
            return Err(eyre!(
                "joint energy sampling has repeated or invalid active/prior loop indices"
            ));
        }
        if external_momenta.is_empty() || external_momenta.len() != lmb.ext_edges.len() {
            return Err(eyre!(
                "joint energy sampling requires all {} external ports, received {}",
                lmb.ext_edges.len(),
                external_momenta.len()
            ));
        }
        let surfaces = [self, other];
        for &surface in &surfaces {
            for &edge in &surface.energies {
                let signature = &lmb.edge_signatures[edge];
                for index in (0..lmb.loop_edges.len())
                    .map(LoopIndex)
                    .filter(|index| !covered.contains(index))
                {
                    if signature.internal[index] != SignOrZero::Zero {
                        return Err(eyre!(
                            "joint energy surface {:?} depends on unsampled parent edge {}",
                            surface.energies,
                            lmb.loop_edges[index]
                        ));
                    }
                }
                if !masses[edge].0.is_finite() {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "joint energy mass",
                        detail: format!("non-finite mass on edge {edge}"),
                    }
                    .into());
                }
                if masses[edge] < masses[edge].zero() {
                    return Err(eyre!(
                        "joint energy edge {edge} requires a nonnegative mass"
                    ));
                }
            }
        }
        let varying = surfaces.map(|surface| {
            surface
                .energies
                .iter()
                .copied()
                .filter(|edge| lmb.edge_signatures[*edge].internal[active] != SignOrZero::Zero)
                .sorted_by_key(|edge| edge.0)
                .collect_vec()
        });
        if varying.iter().any(|edges| edges.len() != 2) {
            return Err(eyre!(
                "joint energy sampling requires exactly two varying energy occurrences per equation, got {:?}",
                varying
            ));
        }
        for (surface, varying) in surfaces.iter().zip(&varying) {
            let contained = subspace.contains(&surface.energies, graph).collect_vec();
            if varying.iter().any(|edge| !contained.contains(edge)) {
                return Err(eyre!(
                    "joint varying energy routing is inconsistent with its selected cycle subgraph"
                ));
            }
        }
        // Keep the original global equations. Replacing this shift by a
        // cut-eliminated identity would move the target by a finite LU residual.
        let shifts =
            surfaces.map(|surface| surface.compute_shift_part_from_momenta(external_momenta, lmb));
        let externals: ExternalThreeMomenta<F<T>> =
            external_momenta.iter().map(|p| p.spatial.clone()).collect();
        if shifts
            .iter()
            .chain(externals.iter().flat_map(|p| [&p.px, &p.py, &p.pz]))
            .any(|value| !value.0.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "joint external data",
                detail: "non-finite external spatial vector or original energy shift".into(),
            }
            .into());
        }
        let mut shared = Vec::new();
        for (left, &a) in varying[0].iter().enumerate() {
            for (right, &b) in varying[1].iter().enumerate() {
                if graph[a].mass_atom() == graph[b].mass_atom()
                    && lmb.edge_signatures[a]
                        .spatial_equality_up_to_sign(&lmb.edge_signatures[b], &externals)?
                {
                    shared.push((left, right));
                }
            }
        }
        let [(left, right)] = shared.as_slice() else {
            return Err(eyre!(
                "joint energy sampling requires exactly one common full signed routing in the supplied spatial frame and mass expression; found {} candidates for {:?}",
                shared.len(),
                varying
            ));
        };
        let edges = [
            varying[0][*left],
            varying[0][1 - *left],
            varying[1][1 - *right],
        ];
        if masses[edges[0]] != masses[varying[1][*right]] {
            return Err(eyre!(
                "joint shared mass expression has inconsistent native values on edges {} and {}",
                edges[0],
                varying[1][*right]
            ));
        }
        let canonical = |edge| {
            let signature = &lmb.edge_signatures[edge];
            if signature.internal[active] == SignOrZero::Minus {
                LoopExtSignature {
                    internal: signature.internal.iter().map(|sign| -*sign).collect(),
                    external: signature.external.iter().map(|sign| -*sign).collect(),
                }
            } else {
                signature.clone()
            }
        };
        let routes = edges.map(canonical);
        let common = routes[0].clone();
        let fixed = surfaces.map(|surface| {
            surface
                .energies
                .iter()
                .copied()
                .filter(|edge| lmb.edge_signatures[*edge].internal[active] == SignOrZero::Zero)
                .collect_vec()
        });
        let lmb = lmb.clone();
        let masses = masses.clone();
        let complement = complement.to_vec();
        let geometry = Arc::new(move |context: &[T]| {
            if context.len() != 3 * complement.len() {
                return Err(eyre!(
                    "joint energy context has {} components, expected {}",
                    context.len(),
                    3 * complement.len()
                ));
            }
            if context.iter().any(|x| !x.is_finite()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "joint energy context",
                    detail: "non-finite declared prerequisite".into(),
                }
                .into());
            }
            let zero = masses[edges[0]].zero();
            let mut loops = LoopMomenta::from_iter(
                (0..lmb.loop_edges.len())
                    .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            for (&index, point) in complement.iter().zip(context.chunks_exact(3)) {
                loops[index] = ThreeMomentum::new(
                    F(point[0].clone()),
                    F(point[1].clone()),
                    F(point[2].clone()),
                );
            }
            let offsets: [ThreeMomentum<F<T>>; 3] = routes
                .each_ref()
                .map(|route| route.compute_momentum(&loops, &externals));
            let differences = [&offsets[1] - &offsets[0], &offsets[2] - &offsets[0]];
            let sums = fixed.iter().enumerate().map(|(index, fixed)| {
                let fixed_sum = fixed.iter().try_fold(zero.clone(), |sum, &edge| -> Result<F<T>> {
                    let momentum: ThreeMomentum<F<T>> = lmb.edge_signatures[edge].compute_momentum(&loops, &externals);
                    let energy = (momentum.norm_squared()+masses[edge].square()).sqrt();
                    if energy == zero && (masses[edge] != zero || momentum.px != zero || momentum.py != zero || momentum.pz != zero) {
                        return Err(SamplingEvaluationError::Unrepresentable { operation: "joint fixed energy", detail: format!("nonzero mass or momentum on edge {edge} produced zero energy") }.into());
                    }
                    Ok(sum+energy)
                })?;
                Ok((-(&shifts[index]+fixed_sum)).0)
            }).collect::<Result<Vec<T>>>()?;
            let energy_sums = [sums[0].clone(), sums[1].clone()];
            let geometry = SharedEnergyJointGeometry {
                shifts: differences.map(|v| [v.px.0, v.py.0, v.pz.0]),
                masses: edges.map(|edge| masses[edge].0.clone()),
                energy_sums,
            };
            if geometry
                .shifts
                .iter()
                .flatten()
                .chain(&geometry.energy_sums)
                .any(|x| !x.is_finite())
            {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "joint prepared energies",
                    detail: "non-finite routed offset or fixed energy sum".into(),
                }
                .into());
            }
            Ok(geometry)
        });
        Ok((geometry, common))
    }

    /// Compile an exact radial chart around this graph-routed energy surface.
    /// The callback retains the complete parent frame, masses and external
    /// data in the evaluation precision; its ray root is the same cut equation
    /// used by LU, without promoting already rounded lower-precision vectors.
    pub(crate) fn sampling_radial_map<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        masses: &EdgeVec<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        beta: f64,
        power: f64,
    ) -> Result<ImplicitSurfaceRadialMap<T>> {
        let dimension = 3 * lmb.loop_edges.len();
        let surface = self.clone();
        let lmb = lmb.clone();
        let masses = masses.clone();
        let external_momenta = external_momenta.clone();
        let zero = F::<T>::from_f64(0.0);
        let center = LoopMomenta::from_iter(
            (0..dimension / 3)
                .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
        );
        let evaluator = std::sync::Arc::new(move |direction: &[T], radius: T| {
            let radius = F(radius);
            let unit_loops = LoopMomenta::from_iter(direction.chunks_exact(3).map(|components| {
                ThreeMomentum::new(
                    F(components[0].clone()),
                    F(components[1].clone()),
                    F(components[2].clone()),
                )
            }));
            let (value, derivative) = surface.compute_self_and_r_derivative(
                &radius,
                &unit_loops,
                &center,
                &external_momenta,
                &masses,
                &lmb,
            );
            Ok((value.0, derivative.0))
        });
        ImplicitSurfaceRadialMap::new(dimension, vec![zero.0; dimension], beta, power, evaluator)
    }

    /// Bind a routed energy surface on active coordinates, conditional on an
    /// already sampled complement in the same complete parent LMB. Later
    /// coordinates may be omitted only when every energy row vanishes on them.
    /// Preparation
    /// never consumes active cube coordinates, so its derivatives occupy only
    /// the off-diagonal block of an ordered composition's Jacobian. Active
    /// outputs follow `subspace.iter_lmb_indices()` (parent-index order), while
    /// complement inputs follow the explicit `complement` order. A compiler
    /// must permute these blocks before attaching user-ordered edge metadata.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sampling_radial_map_in_subspace<T: FloatLike>(
        &self,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        complement: &[LoopIndex],
        settings: &crate::settings::RuntimeSettings,
        beta: f64,
        power: f64,
    ) -> Result<ImplicitSurfaceRadialMap<T>> {
        use crate::integrands::process::PreparedSurfaceStatus;
        use crate::integrands::process::sampling_maps::SamplingEvaluationError;
        use crate::momentum::SignOrZero;
        use crate::subtraction::overlap_subspace::{OverlapInput, find_center};
        use std::sync::Arc;
        use three_dimensional_reps::utils::rank_i64;

        if !settings.kinematics.e_cm.is_finite() || settings.kinematics.e_cm <= 0.0 {
            return Err(eyre!("sampling fiber requires finite positive e_cm"));
        }
        let lmb = subspace.get_lmb(all_lmbs);
        let active = subspace.iter_lmb_indices().collect_vec();
        let mut covered = active.iter().chain(complement).copied().collect_vec();
        covered.sort();
        if active.is_empty()
            || covered.windows(2).any(|pair| pair[0] == pair[1])
            || covered.iter().any(|index| index.0 >= lmb.loop_edges.len())
            || external_momenta.len() != lmb.ext_edges.len()
            || external_momenta.is_empty()
        {
            return Err(eyre!(
                "sampling fiber requires disjoint active/complement slots in its parent frame and complete external momenta: active {active:?}, complement {complement:?}, parent {:?}, external count {} (expected {})",
                lmb.loop_edges,
                external_momenta.len(),
                lmb.ext_edges.len(),
            ));
        }
        // Later blocks may be absent from this prefix only when their exact
        // routed rows vanish. Their zero placeholders then represent proven
        // spectators, never unknown active coordinates of the energy equation.
        for index in (0..lmb.loop_edges.len())
            .map(LoopIndex)
            .filter(|index| !covered.contains(index))
        {
            if self
                .energies
                .iter()
                .any(|edge| lmb.edge_signatures[*edge].internal[index] != SignOrZero::Zero)
            {
                return Err(eyre!(
                    "sampling fiber depends on unsampled parent edge {}",
                    lmb.loop_edges[index]
                ));
            }
        }
        for momentum in external_momenta {
            if !momentum.temporal.value.0.is_finite()
                || [
                    &momentum.spatial.px,
                    &momentum.spatial.py,
                    &momentum.spatial.pz,
                ]
                .iter()
                .any(|value| !value.0.is_finite())
            {
                return Err(eyre!("sampling fiber external momenta must be finite"));
            }
        }
        let rows = self
            .energies
            .iter()
            .map(|edge| {
                active
                    .iter()
                    .map(|&index| match lmb.edge_signatures[*edge].internal[index] {
                        SignOrZero::Minus => -1_i64,
                        SignOrZero::Zero => 0,
                        SignOrZero::Plus => 1,
                    })
                    .collect_vec()
            })
            .collect_vec();
        let rank = rank_i64(&rows);
        if rank != active.len() {
            return Err(eyre!(
                "sampling fiber {:?} has active routing rank {rank}, expected {}; spectator directions require a complement block",
                self.energies,
                active.len()
            ));
        }
        let varying = self
            .energies
            .iter()
            .zip(&rows)
            .filter_map(|(&edge, row)| row.iter().any(|&sign| sign != 0).then_some(edge))
            .collect_vec();
        let subspace_energies = subspace.contains(&self.energies, graph).collect_vec();
        if varying.iter().any(|edge| !subspace_energies.contains(edge)) {
            return Err(eyre!(
                "sampling fiber routing is inconsistent with its selected cycle subgraph"
            ));
        }
        for &edge in &self.energies {
            let mass = &masses[edge];
            if !mass.0.is_finite() || mass < &mass.zero() {
                return Err(eyre!(
                    "sampling fiber edge {edge} requires a finite nonnegative mass"
                ));
            }
        }
        let surface = Arc::new(self.clone());
        let lmbs = Arc::new(all_lmbs.clone());
        let graph = Arc::new(graph.clone());
        let subspace = Arc::new(subspace.clone());
        let masses = Arc::new(masses.clone());
        let externals = Arc::new(external_momenta.clone());
        let settings = Arc::new(settings.clone());
        let complement = Arc::new(complement.to_vec());
        let active = Arc::new(active);
        let zero = external_momenta[ExternalIndex(0)].temporal.value.zero();
        let dimension = active.len() * 3;
        let freeze_context = complement.is_empty();
        let n_loops = lmb.loop_edges.len();
        // One shared routed-ray evaluator serves both full and conditional
        // charts. Complement entries have zero radial velocity.
        let evaluator = {
            let (surface, lmbs, subspace, masses, externals, complement, active) = (
                surface.clone(),
                lmbs.clone(),
                subspace.clone(),
                masses.clone(),
                externals.clone(),
                complement.clone(),
                active.clone(),
            );
            Arc::new(
                move |direction: &[T], radius: T, center: &[T], context: &[T]| {
                    let zero = F(radius.zero());
                    let mut loops = LoopMomenta::from_iter(
                        (0..n_loops)
                            .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
                    );
                    let mut velocity = loops.clone();
                    for (&index, components) in complement.iter().zip(context.chunks_exact(3)) {
                        loops[index] = ThreeMomentum::new(
                            F(components[0].clone()),
                            F(components[1].clone()),
                            F(components[2].clone()),
                        );
                    }
                    for ((&index, components), unit) in active
                        .iter()
                        .zip(center.chunks_exact(3))
                        .zip(direction.chunks_exact(3))
                    {
                        loops[index] = ThreeMomentum::new(
                            F(components[0].clone()),
                            F(components[1].clone()),
                            F(components[2].clone()),
                        );
                        velocity[index] = ThreeMomentum::new(
                            F(unit[0].clone()),
                            F(unit[1].clone()),
                            F(unit[2].clone()),
                        );
                    }
                    let (value, derivative) = surface.compute_self_and_r_derivative(
                        &F(radius),
                        &velocity,
                        &loops,
                        &externals,
                        &masses,
                        subspace.get_lmb(&lmbs),
                    );
                    Ok((value.0, derivative.0))
                },
            )
        };
        let preparer = Arc::new(move |context: &[T]| {
            let zero = externals[ExternalIndex(0)].temporal.value.zero();
            if context.len() != complement.len() * 3
                || context.iter().any(|value| !value.is_finite())
            {
                return Err(eyre!(
                    "sampling fiber expected {} finite complement components, got {}",
                    complement.len() * 3,
                    context.len()
                ));
            }
            let lmb = subspace.get_lmb(&lmbs);
            let spatial: ExternalThreeMomenta<F<T>> =
                externals.iter().map(|p| p.spatial.clone()).collect();
            let mut center = LoopMomenta::from_iter(
                (0..n_loops).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            for (&index, components) in complement.iter().zip(context.chunks_exact(3)) {
                center[index] = ThreeMomentum::new(
                    F(components[0].clone()),
                    F(components[1].clone()),
                    F(components[2].clone()),
                );
            }
            let zero_velocity = LoopMomenta::from_iter(
                (0..n_loops).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
            );
            let mut constant = surface.compute_shift_part_from_momenta(&externals, lmb);
            // Include the inputs to routed differences, not just their small
            // residuals. Large common spatial shifts may otherwise round an
            // actually nonzero separation to zero before minimum testing.
            let mut scale = center
                .iter()
                .flat_map(|p| [&p.px, &p.py, &p.pz])
                .chain(externals.iter().flat_map(|p| {
                    [
                        &p.temporal.value,
                        &p.spatial.px,
                        &p.spatial.py,
                        &p.spatial.pz,
                    ]
                }))
                .fold(constant.abs(), |sum, component| sum + component.abs());
            for &edge in &surface.energies {
                if !varying.contains(&edge) {
                    let momentum = lmb.edge_signatures[edge].compute_momentum(&center, &spatial);
                    let energy = (momentum.norm_squared() + masses[edge].square()).sqrt();
                    constant += &energy;
                    scale += energy;
                }
            }
            let minimum = if active.len() == 1 && varying.len() == 2 {
                // For q_i = s_i k + b_i, s_i = +/-1, the exact convex minimum
                // is sqrt((b1-s1*s2*b2)^2 + (m1+m2)^2). With zero total mass
                // the minimizing set is a segment; its midpoint is a valid
                // proposal center, not a claim that the pinch is isolated.
                let [first, second] = [varying[0], varying[1]];
                let index = active[0];
                let sign = |edge| match lmb.edge_signatures[edge].internal[index] {
                    SignOrZero::Plus => zero.one(),
                    SignOrZero::Minus => -zero.one(),
                    SignOrZero::Zero => unreachable!("two varying rank-one energies"),
                };
                let first_shift = lmb.edge_signatures[first].compute_momentum(&center, &spatial);
                let second_shift = lmb.edge_signatures[second].compute_momentum(&center, &spatial);
                let separation = &first_shift - &(second_shift * (sign(first) * sign(second)));
                let total_mass = &masses[first] + &masses[second];
                let fraction = if total_mass > zero {
                    &masses[first] / &total_mass
                } else {
                    zero.one() / zero.from_i64(2)
                };
                center[index] = (&separation * fraction - first_shift) * sign(first);
                let energy = (separation.norm_squared() + total_mass.square()).sqrt();
                scale += &energy;
                Some(energy + &constant)
            } else {
                None
            };
            let tolerance = zero.epsilon()
                * zero.from_usize(64 * (surface.energies.len() + 1))
                * scale.max(zero.one());
            let normalized = |value: &F<T>| -> Result<T> {
                let margin = -value / F::<T>::from_f64(settings.kinematics.e_cm);
                if !margin.0.is_finite() {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "normalized sampling surface margin",
                        detail: "finite surface value and energy scale produce a nonfinite margin"
                            .to_owned(),
                    }
                    .into());
                }
                Ok(margin.0)
            };
            let status = if let Some(minimum) = minimum {
                if !minimum.0.is_finite() || minimum.abs() <= tolerance {
                    return Err(SamplingEvaluationError::UncertainGeometry { detail: format!("two-energy fiber minimum is not sign-certified: minimum={minimum}, tolerance={tolerance}, edges={:?}", surface.energies) }.into());
                }
                if minimum > zero {
                    PreparedSurfaceStatus::absent(format!(
                        "positive exact two-energy minimum {minimum}"
                    ))?
                } else {
                    PreparedSurfaceStatus::existing(Some(normalized(&minimum)?))?
                }
            } else {
                // The sum of masses is a necessary lower bound for every
                // routed topology. A positive bound certifies absence; a
                // failed approximate optimizer by itself never does.
                let lower = varying
                    .iter()
                    .fold(constant, |sum, &edge| sum + &masses[edge]);
                if lower > tolerance {
                    PreparedSurfaceStatus::absent(format!("positive energy lower bound {lower}"))?
                } else {
                    let value = surface
                        .compute_self_and_r_derivative(
                            &zero,
                            &zero_velocity,
                            &center,
                            &externals,
                            &masses,
                            lmb,
                        )
                        .0;
                    if !esurface_value_is_strictly_inside(
                        &value,
                        &F::<T>::from_f64(settings.kinematics.e_cm),
                    ) {
                        let thresholds = EsurfaceCollection::from(vec![surface.as_ref().clone()]);
                        let existing = ExistingThresholds::from(vec![EsurfaceID(0)]);
                        let input = OverlapInput {
                            graph: &graph,
                            settings: &settings,
                            subspace: &subspace,
                            threshold_subspaces: None,
                            lmbs: &lmbs,
                            thresholds: &thresholds,
                            edge_masses: masses.iter().map(|(_, mass)| mass.into_ff64()).collect(),
                            surface_kinematics: None,
                        };
                        let loops: LoopMomenta<F<f64>> =
                            center.iter().map(ThreeMomentum::to_f64).collect();
                        let external: ExternalFourMomenta<F<f64>> =
                            externals.iter().map(|p| p.to_f64()).collect();
                        if loops
                            .iter()
                            .flat_map(|p| [&p.px, &p.py, &p.pz])
                            .any(|v| !v.0.is_finite())
                            || external
                                .iter()
                                .flat_map(|p| {
                                    [
                                        &p.temporal.value,
                                        &p.spatial.px,
                                        &p.spatial.py,
                                        &p.spatial.pz,
                                    ]
                                })
                                .any(|v| !v.0.is_finite())
                            || input
                                .edge_masses
                                .iter()
                                .any(|(_, mass)| !mass.0.is_finite())
                        {
                            return Err(SamplingEvaluationError::UncertainGeometry {
                                detail: "native fiber data exceed the finite range of the approximate f64 SOCP center seed".to_owned(),
                            }.into());
                        }
                        let candidate = find_center(&input, &[ExistingEsurfaceId(0)], &existing, &loops, &external, false)
                            .map_err(|error| SamplingEvaluationError::UncertainGeometry { detail: format!("fiber center seed failed without an absence certificate: {error}") })?
                            .ok_or_else(|| SamplingEvaluationError::UncertainGeometry { detail: "fiber SOCP reported infeasibility without a native absence certificate".to_owned() })?;
                        for &index in active.iter() {
                            center[index] = ThreeMomentum::from_ff64(candidate[index]);
                        }
                    }
                    let value = surface
                        .compute_self_and_r_derivative(
                            &zero,
                            &zero_velocity,
                            &center,
                            &externals,
                            &masses,
                            lmb,
                        )
                        .0;
                    if !esurface_value_is_strictly_inside(
                        &value,
                        &F::<T>::from_f64(settings.kinematics.e_cm),
                    ) {
                        return Err(SamplingEvaluationError::UncertainGeometry {
                            detail: format!(
                                "fiber center is not natively certified inside: value={value}"
                            ),
                        }
                        .into());
                    }
                    PreparedSurfaceStatus::existing(Some(normalized(&value)?))?
                }
            };
            if status.is_existing() {
                // The analytic minimum does not excuse a rounded minimizing
                // point: the actual bound center must also be strictly inside.
                let value = surface
                    .compute_self_and_r_derivative(
                        &zero,
                        &zero_velocity,
                        &center,
                        &externals,
                        &masses,
                        lmb,
                    )
                    .0;
                if !esurface_value_is_strictly_inside(
                    &value,
                    &F::<T>::from_f64(settings.kinematics.e_cm),
                ) || value >= -&tolerance
                {
                    return Err(SamplingEvaluationError::UncertainGeometry {
                        detail: format!(
                            "prepared fiber center is not natively certified inside: value={value}, tolerance={tolerance}"
                        ),
                    }
                    .into());
                }
            }
            let center = active
                .iter()
                .flat_map(|&index| {
                    [
                        center[index].px.0.clone(),
                        center[index].py.0.clone(),
                        center[index].pz.0.clone(),
                    ]
                })
                .collect();
            Ok((center, status))
        });
        ImplicitSurfaceRadialMap::new(
            dimension,
            vec![zero.0; dimension],
            beta,
            power,
            Arc::new(|_, _| Err(eyre!("sampling fiber requires prepared complement data"))),
        )
        .and_then(|map| {
            let map = map
                .with_context_evaluator(evaluator)
                .with_context_preparer(preparer);
            if freeze_context {
                map.freeze_context(&[])
            } else {
                Ok(map)
            }
        })
    }

    pub(crate) fn has_radial_dependence_in_subspace(
        &self,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
    ) -> bool {
        let lmb = subspace.get_lmb(all_lmbs);
        subspace.contains(&self.energies, graph).any(|index| {
            subspace
                .project_loop_signature(&lmb.edge_signatures[index].internal)
                .any(|sign| sign.is_sign())
        })
    }

    pub(crate) fn external_shift_is_strictly_negative_for_positive_energies(
        &self,
        incoming_edges: &[EdgeIndex],
        outgoing_edges: &[EdgeIndex],
    ) -> bool {
        if incoming_edges.is_empty()
            || outgoing_edges.is_empty()
            || incoming_edges
                .iter()
                .any(|edge| outgoing_edges.contains(edge))
            || self
                .external_shift
                .iter()
                .any(|(edge, _)| !incoming_edges.contains(edge) && !outgoing_edges.contains(edge))
        {
            return false;
        }

        let coefficient = |edge: &EdgeIndex| {
            self.external_shift
                .iter()
                .filter(|(shift_edge, _)| shift_edge == edge)
                .map(|(_, coefficient)| i128::from(*coefficient))
                .sum::<i128>()
        };

        // Energy conservation makes external-energy coefficient vectors `c` and
        // `c + lambda * sigma` equivalent, where sigma is +1 for incoming and -1
        // for outgoing momenta. The shift is strictly negative for all positive
        // external energies if one representative is component-wise non-positive
        // and not identically zero.
        let lambda_lower_bound = outgoing_edges
            .iter()
            .map(&coefficient)
            .max()
            .expect("outgoing external edges were checked to be non-empty");
        let lambda_upper_bound = incoming_edges
            .iter()
            .map(|edge| -coefficient(edge))
            .min()
            .expect("incoming external edges were checked to be non-empty");

        if lambda_lower_bound > lambda_upper_bound {
            return false;
        }

        let lambda = lambda_lower_bound;
        let adjusted_coefficients = incoming_edges
            .iter()
            .map(|edge| coefficient(edge) + lambda)
            .chain(outgoing_edges.iter().map(|edge| coefficient(edge) - lambda));
        let mut has_strictly_negative_coefficient = false;
        for adjusted_coefficient in adjusted_coefficients {
            if adjusted_coefficient > 0 {
                return false;
            }
            has_strictly_negative_coefficient |= adjusted_coefficient < 0;
        }

        has_strictly_negative_coefficient
    }

    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        let symbolic_energies = self
            .energies
            .iter()
            .map(|i| {
                if cut_edges.contains(i) {
                    cut_energy(*i)
                } else {
                    ose_atom_from_index(*i)
                }
            })
            .collect_vec();

        let symbolic_shift = self
            .external_shift
            .iter()
            .fold(Atom::new(), |sum, (i, sign)| {
                external_energy_atom_from_index(*i) * &Atom::num(*sign) + &sum
            });

        let builder_atom = Atom::new();
        let energy_sum = symbolic_energies
            .iter()
            .fold(builder_atom, |acc, energy| acc + energy);

        energy_sum + &symbolic_shift
    }

    #[inline]
    pub(crate) fn compute_from_dual_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        dual_loop_moms: &LoopMomenta<HyperDual<F<T>>>,
        dual_external_moms: &ExternalFourMomenta<HyperDual<F<T>>>,
    ) -> HyperDual<F<T>> {
        let spatial_part_of_externals = dual_external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                let momentum = signature
                    .try_compute_momentum(&dual_loop_moms.0, &spatial_part_of_externals.raw)
                    .unwrap_or_else(|| unreachable!());
                let mass = &real_mass_vector[*index];

                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| dual_loop_moms[LoopIndex(0)].px.zero());

        let shift_part = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                let sign = energy_sum.values[0].from_i64(*sign);
                new_constant(&energy_sum, &sign)
                    * external_signature
                        .try_apply(&dual_external_moms.raw)
                        .map(|mom| mom.temporal.value)
                        .unwrap_or_else(|| energy_sum.zero())
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| energy_sum.zero());

        energy_sum + shift_part
    }

    /// Compute the value of the esurface from the momenta, needed to check if an arbitrary point
    /// is inside the esurface
    #[inline]
    pub(crate) fn compute_from_momenta<T: FloatLike>(
        &self,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) -> F<T> {
        let spatial_part_of_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let energy_sum = self
            .energies
            .iter()
            .map(|index| {
                let signature = &lmb.edge_signatures[*index];
                let momentum = signature.compute_momentum(loop_moms, &spatial_part_of_externals);
                let mass = &real_mass_vector[*index];

                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| loop_moms[LoopIndex(0)].px.zero());

        let shift_part = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                energy_sum.from_i64(*sign)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| energy_sum.zero());

        energy_sum + shift_part
    }

    fn classify_invariant_margin<T: FloatLike>(
        shift_part: &F<T>,
        invariant_margin: F<T>,
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistence<T> {
        // /!\ In alphaLoop this boundary needed special care so that `mul_unit` could not
        // create a discontinuous non-pinched-to-pinched transition for a massless 2->2
        // E-surface sandwich. Such a surface is `Existing` at non-collinear kinematics, where
        // its invariant margin is positive, and `Pinched` only at the collinear null boundary.
        // Here the normalized tolerance band around that boundary is classified as `Pinched`,
        // so perturbations within the band cannot toggle threshold subtraction. GammaLoop's
        // statuses are exclusive and only `Existing` surfaces are subtraction targets.
        // Keep the decision in the original, dimensionful scale. Besides avoiding an
        // unnecessary division, this preserves the previous existence boundary exactly.
        // Deserialization rejects invalid configured tolerances. Keep this defensive fallback for
        // settings mutated programmatically (for example through a language binding), so an
        // invalid sign can never invert the existence test.
        let normalized_margin_tolerance =
            if normalized_margin_tolerance.is_nan() || normalized_margin_tolerance.is_infinite() {
                F::from_f64(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD)
            } else {
                normalized_margin_tolerance.abs()
            };
        let invariant_tolerance = &normalized_margin_tolerance * e_cm * e_cm;
        let normalized_margin = &invariant_margin / (e_cm * e_cm);
        let shift_tolerance = F::from_f64(ESURFACE_SHIFT_THRESHOLD) * e_cm;

        if invariant_margin.abs() <= invariant_tolerance && shift_part <= &shift_tolerance {
            return EsurfaceExistence::Pinched { normalized_margin };
        }

        if shift_part >= &(-shift_tolerance) {
            return EsurfaceExistence::NonExisting {
                normalized_margin: Some(normalized_margin),
                reason: NonExistingEsurfaceReason::ShiftNotNegative,
            };
        }

        if invariant_margin > invariant_tolerance {
            EsurfaceExistence::Existing { normalized_margin }
        } else if invariant_margin >= -invariant_tolerance {
            EsurfaceExistence::Pinched { normalized_margin }
        } else {
            EsurfaceExistence::NonExisting {
                normalized_margin: Some(normalized_margin),
                reason: NonExistingEsurfaceReason::NoRealZero,
            }
        }
    }

    #[inline]
    #[allow(clippy::too_many_arguments)]
    /// Classify a surface in an active loop-momentum subspace. Pinched and
    /// non-existing surfaces remain distinct and are never subtraction targets.
    pub(crate) fn classify_existence_subspace<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        real_mass_vector: &EdgeVec<F<T>>,
        reversed_edges: &[EdgeIndex],
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistence<T> {
        //todo!("refactor for subspaces");
        if self.external_shift.is_empty() {
            debug!("esurface has no external shift, cannot exist");
            return EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoExternalShift,
            };
        }

        let shift_part = self.compute_shift_part_from_momenta_in_subspace(
            loop_moms,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            real_mass_vector,
        );

        let subspace_energy_indices = subspace.contains(&self.energies, graph).collect_vec();

        if !self.has_radial_dependence_in_subspace(subspace, all_lmbs, graph) {
            debug!(
                "esurface has no radial energy in this subspace, cannot bound a threshold region"
            );
            return EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoRadialDependence,
            };
        }

        let lmb = subspace.get_lmb(all_lmbs);
        let mass_sum: F<T> = subspace_energy_indices
            .iter()
            .map(|&index| &real_mass_vector[index])
            .fold(F::from_f64(0.0), |acc, x| acc + x);

        let zero_vector = ThreeMomentum::new(e_cm.zero(), e_cm.zero(), e_cm.zero());

        let graph_vector = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                compute_shift_part(external_signature, external_moms).spatial
                    * F::from_f64(*sign as f64)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| zero_vector.clone());

        let other_part = subspace
            .does_not_contain(&self.energies, graph)
            .map(|index| {
                let signature = &lmb.edge_signatures[index];
                let sign = if reversed_edges.contains(&index) {
                    -F::from_f64(1.0)
                } else {
                    F::from_f64(1.0)
                };

                signature.compute_momentum(
                    loop_moms,
                    &external_moms
                        .iter()
                        .map(|mom| mom.spatial.clone())
                        .collect::<TiVec<ExternalIndex, _>>(),
                ) * sign
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| zero_vector.clone());

        let shift_vector_sq = (&graph_vector + &other_part).norm_squared();
        let invariant_margin = &shift_part * &shift_part - &shift_vector_sq - &mass_sum * &mass_sum;
        let classification = Self::classify_invariant_margin(
            &shift_part,
            invariant_margin,
            e_cm,
            normalized_margin_tolerance,
        );

        if !classification.is_existing() {
            debug!(
                "subspace esurface classified as {}: shift_part^2: {}, shift_vector_sq: {}, mass_sum^2: {}, normalized_margin: {:?}",
                classification.label(),
                &shift_part * &shift_part,
                shift_vector_sq,
                &mass_sum * &mass_sum,
                classification.normalized_margin(),
            );
        }

        classification
    }

    #[inline]
    /// Classify a full-space surface. Pinched and non-existing surfaces remain
    /// distinct and are never subtraction targets.
    pub(crate) fn classify_existence<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistence<T> {
        //todo!("refactor for subspaces");
        if self.external_shift.is_empty() {
            return EsurfaceExistence::NonExisting {
                normalized_margin: None,
                reason: NonExistingEsurfaceReason::NoExternalShift,
            };
        }

        let shift_part = self.compute_shift_part_from_momenta(external_moms, lmb);
        let mass_sum: F<T> = self
            .energies
            .iter()
            .map(|index| &real_mass_vector[*index])
            .fold(F::from_f64(0.0), |acc, x| acc + x);

        let zero_vector = ThreeMomentum::new(e_cm.zero(), e_cm.zero(), e_cm.zero());

        let shift_vector = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                compute_shift_part(external_signature, external_moms).spatial
                    * F::from_f64(*sign as f64)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| zero_vector.clone());

        let shift_vector_sq = shift_vector.norm_squared();
        let invariant_margin = &shift_part * &shift_part - shift_vector_sq - &mass_sum * &mass_sum;

        Self::classify_invariant_margin(
            &shift_part,
            invariant_margin,
            e_cm,
            normalized_margin_tolerance,
        )
    }

    /// Classify a full-space surface while keeping normalized margins and rejection reasons
    /// internal to the subtraction implementation.
    pub fn existence_status<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
        real_mass_vector: &EdgeVec<F<T>>,
        e_cm: &F<T>,
        normalized_margin_tolerance: &F<T>,
    ) -> EsurfaceExistenceStatus {
        self.classify_existence(
            external_moms,
            lmb,
            real_mass_vector,
            e_cm,
            normalized_margin_tolerance,
        )
        .status()
    }

    /// Only compute the shift part, useful for center finding.
    pub(crate) fn compute_shift_part_from_momenta_in_subspace<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
    ) -> F<T> {
        let lmb = subspace.get_lmb(all_lmbs);

        let full_external_shift = self
            .external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                external_moms[ExternalIndex(0)]
                    .temporal
                    .value
                    .from_i64(*sign)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| external_moms[ExternalIndex(0)].temporal.value.zero());

        let spatial_externals = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect::<TiVec<ExternalIndex, _>>();

        let remaining_shift = subspace
            .does_not_contain(&self.energies, graph)
            .map(|index| {
                let signature = &lmb.edge_signatures[index];
                //panic!("signature: {:?}", signature);
                let momentum = signature.compute_momentum(loop_moms, &spatial_externals);
                let mass = &masses[index];

                (momentum.norm_squared() + mass * mass).sqrt()
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| full_external_shift.zero());

        full_external_shift + remaining_shift
    }

    /// Only compute the shift part, useful for center finding.
    pub(crate) fn compute_shift_part_from_momenta<T: FloatLike>(
        &self,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> F<T> {
        self.external_shift
            .iter()
            .map(|(index, sign)| {
                let external_signature = &lmb.edge_signatures[*index].external;
                external_moms[ExternalIndex(0)]
                    .temporal
                    .value
                    .from_i64(*sign)
                    * compute_t_part_of_shift_part(external_signature, external_moms)
            })
            .reduce(|acc, x| acc + x)
            .unwrap_or_else(|| external_moms[ExternalIndex(0)].temporal.value.zero())
    }

    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute_self_and_r_derivative_subspace<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops_in_subspace: &LoopMomenta<F<T>>,
        center_in_subspace: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        real_mass_vector: &EdgeVec<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect();

        let loops: LoopMomenta<F<T>> = shifted_unit_loops_in_subspace
            .iter_enumerated()
            .map(|(loop_index, shifted_unit_momenta)| {
                if subspace.contains_loop_index(loop_index) {
                    shifted_unit_momenta * radius + &center_in_subspace[loop_index]
                } else {
                    shifted_unit_momenta.clone()
                }
            })
            .collect();

        let shift = self.compute_shift_part_from_momenta_in_subspace(
            shifted_unit_loops_in_subspace,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            real_mass_vector,
        );

        let lmb = subspace.get_lmb(all_lmbs);
        let (derivative, energy_sum) = subspace
            .contains(&self.energies, graph)
            .map(|index| {
                let signature = &lmb.edge_signatures[index];

                let momentum = signature.compute_momentum(&loops, &spatial_part_of_externals);
                let unit_loop_part = compute_loop_part_subspace(
                    &signature.internal,
                    shifted_unit_loops_in_subspace,
                    subspace,
                );

                let energy = (momentum.norm_squared()
                    + &real_mass_vector[index] * &real_mass_vector[index])
                    .sqrt();

                let numerator = momentum * &unit_loop_part;

                (numerator / &energy, energy)
            })
            .fold(
                (radius.zero(), radius.zero()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );

        (energy_sum + shift, derivative)
    }

    /// Shared routed ray evaluation for physical roots and full-space/proper-fiber charts.
    #[inline]
    pub(crate) fn compute_self_and_r_derivative<T: FloatLike>(
        &self,
        radius: &F<T>,
        shifted_unit_loops: &LoopMomenta<F<T>>,
        center: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        real_mass_vector: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let spatial_part_of_externals: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|mom| mom.spatial.clone())
            .collect();

        let loops: LoopMomenta<F<T>> = shifted_unit_loops
            .iter()
            .zip(center.iter())
            .map(|(momentum, center)| momentum * radius + center)
            .collect();

        let shift = self.compute_shift_part_from_momenta(external_moms, lmb);

        let (derivative, energy_sum) = self
            .energies
            .iter()
            .map(|&index| {
                let signature = &lmb.edge_signatures[index];

                let momentum = signature.compute_momentum(&loops, &spatial_part_of_externals);
                let unit_loop_part = compute_loop_part(&signature.internal, shifted_unit_loops);

                let energy = (momentum.norm_squared()
                    + &real_mass_vector[index] * &real_mass_vector[index])
                    .sqrt();

                // At a massless endpoint E(r)=r|v| the radial right derivative is
                // |v|, whereas the two-sided formula q.v/E would evaluate 0/0.
                // Share this convention with sampling charts, using exact source
                // zeros: a computed E=0 can instead be numerical underflow.
                let zero = radius.zero();
                let derivative = if radius == &zero
                    && real_mass_vector[index] == zero
                    && momentum.px == zero
                    && momentum.py == zero
                    && momentum.pz == zero
                {
                    unit_loop_part.norm_squared().sqrt()
                } else {
                    momentum * &unit_loop_part / &energy
                };

                (derivative, energy)
            })
            .fold(
                (radius.zero(), radius.zero()),
                |(der_sum, en_sum), (der, en)| (der_sum + der, en_sum + en),
            );

        (energy_sum + shift, derivative)
    }

    /// Enclose the original graph equation at `radius * loop_momenta`.
    /// Route integer coefficients before rounding energy sums; importing an
    /// already represented `routed_ray` would omit this routing error. True
    /// external shifts and every repeated energy occurrence remain unchanged.
    pub(crate) fn evaluate_routed_enclosed<T: FloatLike>(
        &self,
        radius: &F<T>,
        loop_momenta: &LoopMomenta<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> Result<[Float; 2]> {
        use linnet::num_traits::SignOrZero;
        if loop_momenta.len() != lmb.loop_edges.len()
            || external_momenta.len() != lmb.ext_edges.len()
        {
            return Err(eyre!(
                "routed E-surface enclosure requires the complete parent and external frame"
            ));
        }
        if !radius.0.is_finite()
            || loop_momenta
                .iter()
                .flat_map(|p| [&p.px, &p.py, &p.pz])
                .chain(external_momenta.iter().flat_map(|p| {
                    [
                        &p.temporal.value,
                        &p.spatial.px,
                        &p.spatial.py,
                        &p.spatial.pz,
                    ]
                }))
                .any(|x| !x.0.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "routed E-surface enclosure",
                detail: "nonfinite original radius or momentum component".into(),
            }
            .into());
        }
        for edge in self
            .energies
            .iter()
            .chain(self.external_shift.iter().map(|(edge, _)| edge))
        {
            let signature = lmb.edge_signatures.get(*edge).ok_or_else(|| {
                eyre!("routed E-surface enclosure has no signature for edge {edge}")
            })?;
            // Empty signatures are exact zero routes; nonempty signatures must
            // describe the full frame, never a silently truncated prefix.
            if (!signature.internal.is_empty() && signature.internal.len() != loop_momenta.len())
                || (!signature.external.is_empty()
                    && signature.external.len() != external_momenta.len())
            {
                return Err(eyre!(
                    "routed E-surface enclosure has an incomplete signature for edge {edge}"
                ));
            }
        }
        let coefficient = |sign: &SignOrZero| match sign {
            SignOrZero::Plus => 1,
            SignOrZero::Minus => -1,
            SignOrZero::Zero => 0,
        };
        let spatial = |p: &ThreeMomentum<F<T>>, axis| {
            EsurfaceRay::<T>::enclosed_native([&p.px, &p.py, &p.pz][axis])
        };
        let mut energies = Vec::with_capacity(self.energies.len());
        for &edge in &self.energies {
            let mass = masses
                .get(edge)
                .ok_or_else(|| eyre!("routed E-surface enclosure has no mass for edge {edge}"))?;
            if !mass.0.is_finite() {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "routed E-surface enclosure",
                    detail: format!("nonfinite original mass on edge {edge}"),
                }
                .into());
            }
            let signature = &lmb.edge_signatures[edge];
            let velocity = std::array::from_fn(|axis| {
                EsurfaceRay::<T>::enclosed_sum(
                    signature
                        .internal
                        .iter()
                        .zip(loop_momenta.iter())
                        .map(|(sign, p)| (coefficient(sign), spatial(p, axis))),
                )
            });
            let offset = std::array::from_fn(|axis| {
                EsurfaceRay::<T>::enclosed_sum(
                    signature
                        .external
                        .iter()
                        .zip(external_momenta.iter())
                        .map(|(sign, p)| (coefficient(sign), spatial(&p.spatial, axis))),
                )
            });
            energies.push((velocity, offset, EsurfaceRay::<T>::enclosed_native(mass)));
        }
        let shift =
            EsurfaceRay::<T>::enclosed_sum(self.external_shift.iter().map(|(edge, scale)| {
                let temporal = EsurfaceRay::<T>::enclosed_sum(
                    lmb.edge_signatures[*edge]
                        .external
                        .iter()
                        .zip(external_momenta.iter())
                        .map(|(sign, p)| {
                            (
                                coefficient(sign),
                                EsurfaceRay::<T>::enclosed_native(&p.temporal.value),
                            )
                        }),
                );
                (*scale, temporal)
            }));
        Ok(EsurfaceRay::<T>::evaluate_enclosed(
            energies.into_iter(),
            shift,
            &EsurfaceRay::<T>::enclosed_native(radius),
            false,
        )?
        .0)
    }

    /// Explicit prepared-ray boundary for LU root and jet evaluation. Routing
    /// before radial scaling changes finite-precision association. Existing
    /// non-LU amplitude/static/fiber evaluation keeps scale-then-route semantics.
    /// The center is supplied explicitly; physical LU uses generation zero.
    pub(crate) fn routed_ray<T: FloatLike>(
        &self,
        velocity: &LoopMomenta<F<T>>,
        center: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> EsurfaceRay<T> {
        let spatial: ExternalThreeMomenta<F<T>> = external_moms
            .iter()
            .map(|momentum| momentum.spatial.clone())
            .collect();
        EsurfaceRay {
            energies: self
                .energies
                .iter()
                .map(|&edge| {
                    let signature = &lmb.edge_signatures[edge];
                    (
                        edge,
                        compute_loop_part(&signature.internal, velocity),
                        signature.compute_momentum(center, &spatial),
                        masses[edge].clone(),
                    )
                })
                .collect(),
            shift: self.compute_shift_part_from_momenta(external_moms, lmb),
        }
    }

    /// Solve the physical LU scaling about zero generation momentum. Sampling
    /// hosts use the same equation and policy; the caller owns the identity and
    /// precision history, so sharing this entry does not imply shared root data.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_lu_cut<T: FloatLike>(
        &self,
        loop_momenta: &LoopMomenta<F<T>>,
        external_momenta: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        lmb: &LoopMomentumBasis,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> std::result::Result<(EsurfaceRay<T>, NewtonIterationResult<T>), SafeguardedNewtonError<T>>
    {
        let zero = loop_momenta[LoopIndex(0)].px.zero();
        let center = LoopMomenta::from_iter(
            (0..loop_momenta.len())
                .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
        );
        let ray = self.routed_ray(loop_momenta, &center, external_momenta, masses, lmb);
        let solution = ray.solve_lu_cut(e_cm, diagnostics, identity)?;
        Ok((ray, solution))
    }

    // #[inline]
    /// the "loops_unit_in_subspace" means that the loop momenta that are part of the subspace are jointly normalized to unit length
    /// A ray with no radial dependence has no isolated root and returns zero guesses.
    pub(crate) fn get_radius_guess_subspace<T: FloatLike>(
        &self,
        loops_unit_in_subspace: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        subspace: &SubspaceData,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        graph: &Graph,
        masses: &EdgeVec<F<T>>,
    ) -> (F<T>, F<T>) {
        let const_builder = &loops_unit_in_subspace[LoopIndex(0)].px;

        let esurface_shift = self.compute_shift_part_from_momenta_in_subspace(
            loops_unit_in_subspace,
            external_moms,
            subspace,
            all_lmbs,
            graph,
            masses,
        );

        debug!("shift part for radius guess: {}", esurface_shift);

        let mut radius_guess = const_builder.zero();
        let mut denominator = const_builder.zero();

        let lmb = subspace.get_lmb(all_lmbs);

        debug!("unit loops in subspace: {:?}", loops_unit_in_subspace);
        let external_3_momenta = external_moms.iter().map(|x| x.spatial.clone()).collect();

        //println!("got to energy loop");
        for energy in subspace.contains(&self.energies, graph) {
            debug!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[energy];
            //println!("signature {:?}", signature);

            let unit_loop_part =
                compute_loop_part_subspace(&signature.internal, loops_unit_in_subspace, subspace);
            //println!("computed_loop_part {:?}", unit_loop_part);

            let three_shift = compute_shift_part_subspace(
                &signature.internal,
                &signature.external,
                loops_unit_in_subspace,
                &external_3_momenta,
                subspace,
            );
            //./bprintln!("computed_shift {:?}", shift);

            let norm_unit_loop_part_squared = unit_loop_part.norm_squared();
            // An energy with zero radial momentum is constant along this ray. It contributes
            // no large-radius growth or directional shift, even when its mass is nonzero.
            if norm_unit_loop_part_squared == const_builder.zero() {
                continue;
            }
            let loop_dot_shift = &unit_loop_part * three_shift;

            debug!(
                "norm unit loop part squared: {}",
                norm_unit_loop_part_squared
            );

            radius_guess += loop_dot_shift.abs() / &norm_unit_loop_part_squared;
            debug!("current radius guess: {}", radius_guess);
            denominator += norm_unit_loop_part_squared.sqrt();
            debug!("current denominator: {}", denominator);
        }

        if denominator != const_builder.zero() {
            radius_guess += esurface_shift.abs() / denominator;
        }
        debug!("final radius guess: {}", radius_guess);
        let negative_radius = radius_guess.ref_neg();
        (radius_guess, negative_radius)
    }

    /// A ray with no radial dependence has no isolated root and returns zero guesses.
    pub(crate) fn get_radius_guess<T: FloatLike>(
        &self,
        unit_loops: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
        lmb: &LoopMomentumBasis,
    ) -> (F<T>, F<T>) {
        let shift = self.compute_shift_part_from_momenta(external_moms, lmb);
        //println!("got to energy loop");
        let terms = self.energies.iter().map(|&energy| {
            //println!("computing contribution for energy {:?}", energy);
            let signature = &lmb.edge_signatures[energy];
            //println!("signature {:?}", signature);
            let unit_loop_part = compute_loop_part(&signature.internal, unit_loops);
            //println!("computed_loop_part {:?}", unit_loop_part);
            let shift = compute_shift_part(&signature.external, external_moms).spatial;
            //./bprintln!("computed_shift {:?}", shift);
            (unit_loop_part.norm_squared(), &unit_loop_part * shift)
        });
        // Existing CT callers intentionally estimate with true external offsets,
        // even when their later root evaluation uses a nonzero overlap center.
        let guess = Self::radius_guess_from_terms(&shift, terms);
        let negative = guess.ref_neg();
        (guess, negative)
    }

    fn radius_guess_from_terms<T: FloatLike>(
        shift: &F<T>,
        terms: impl Iterator<Item = (F<T>, F<T>)>,
    ) -> F<T> {
        let zero = shift.zero();
        let mut guess = zero.clone();
        let mut denominator = zero.clone();
        // Each pair is (|v|^2, v.b); LU supplies its represented ray, while
        // non-LU callers preserve the existing external-only seed convention.
        for (norm_squared, dot_offset) in terms {
            // Constant energies have no directional contribution to the radius estimate.
            if norm_squared == zero {
                continue;
            }
            guess += dot_offset.abs() / &norm_squared;
            denominator += norm_squared.sqrt();
        }
        if denominator != zero {
            guess += shift.abs() / denominator;
        }
        guess
    }

    pub(crate) fn canonicalize_shift(&mut self, shift_rewrite: &ShiftRewrite) {
        if let Some(dep_mom_pos) = self
            .external_shift
            .iter()
            .position(|(index, _)| *index == shift_rewrite.dependent_momentum)
        {
            let (_, dep_mom_sign) = self.external_shift.remove(dep_mom_pos);

            let external_shift = shift_rewrite
                .dependent_momentum_expr
                .iter()
                .map(|(index, sign)| (*index, dep_mom_sign * sign))
                .collect();

            self.external_shift = add_external_shifts(&self.external_shift, &external_shift);
        }
    }

    pub(crate) fn new_from_cut_left<E, V, H>(
        graph: &HedgeGraph<E, V, H>,
        cut: &CrossSectionCut,
        initial_state_cut: Option<&OrientedCut>,
    ) -> Self {
        let edges = graph
            .iter_edges_of(&cut.cut)
            .map(|(_, id, _)| id)
            .sorted()
            .collect();

        let external_shift = if let Some(is_cut) = initial_state_cut {
            graph
                .iter_edges_of(is_cut)
                .map(|(_, edge_index, __)| (edge_index, -1))
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .collect()
        } else {
            graph
                .iter_edges_of(&cut.left)
                .filter_map(|(hedge_pair, edge_index, _)| match hedge_pair {
                    HedgePair::Unpaired { flow, .. } => match flow {
                        Flow::Sink => Some((edge_index, -1)),
                        Flow::Source => Some((edge_index, 1)),
                    },
                    _ => None,
                })
                .sorted_by(|a, b| a.0.cmp(&b.0))
                .collect()
        };

        let vertex_set = graph
            .iter_nodes_of(&cut.left)
            .map(|(node_id, _, _)| VertexSet::from_usize(node_id.into()))
            .reduce(|acc, v| acc.join(&v))
            .unwrap();

        Self {
            energies: edges,
            external_shift,
            vertex_set,
            //subspace_graph: graph.full_graph(),
        }
    }

    pub(crate) fn lmb_atom(&self, graph: &Graph, lmb_reps: &[Replacement]) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom();
                let emr_symbols = (0..3)
                    .map(|i| function!(GS.emr_mom, usize::from(*index), i + 1))
                    .collect_vec();

                (&emr_symbols[0] * &emr_symbols[0]
                    + &emr_symbols[1] * &emr_symbols[1]
                    + &emr_symbols[2] * &emr_symbols[2]
                    + &mass_symbol * &mass_symbol)
                    .sqrt()
            })
            .chain(self.external_shift.iter().map(|(index, sign)| {
                function!(GS.emr_mom, usize::from(*index), 0) * Atom::num(*sign)
            }))
            .reduce(|sum, atom| sum + atom)
            .unwrap_or_else(Atom::new)
            .replace_multiple(lmb_reps)
            .replace(parse!("ZERO"))
            .with(Atom::new())
            .expand() // ensure canonical form
    }

    // more readable version for debugging, because it doesn't write out components
    pub(crate) fn lmb_atom_simplified(&self, graph: &Graph, lmb_reps: &[Replacement]) -> Atom {
        self.energies
            .iter()
            .map(|index| {
                let mass_symbol = graph.underlying[*index].mass_atom();
                let emr_symbol = function!(GS.emr_mom, usize::from(*index));

                (&emr_symbol * &emr_symbol + &mass_symbol * &mass_symbol).sqrt()
            })
            .chain(self.external_shift.iter().map(|(index, sign)| {
                function!(GS.emr_mom, usize::from(*index), 0) * Atom::num(*sign)
            }))
            .reduce(|sum, atom| sum + atom)
            .unwrap_or_else(Atom::new)
            .replace_multiple(lmb_reps)
            .replace(parse!("ZERO"))
            .with(Atom::new())
            .expand() // ensure canonical form
    }
}

define_index! {pub struct GroupEsurfaceId;}

pub type EsurfaceCollection = TiVec<EsurfaceID, Esurface>;

pub type EsurfaceCache<T> = TiVec<EsurfaceID, T>;

/// Container for esurfaces that exist at a given point in the phase space
pub type ExistingEsurfaces = TiVec<ExistingEsurfaceId, GroupEsurfaceId>;
pub type ExistingThresholds = TiVec<ExistingEsurfaceId, EsurfaceID>;

pub(crate) fn get_representative<T: Copy>(
    esurface_map: &TiVec<GraphGroupPosition, Option<T>>,
) -> Result<(GraphGroupPosition, T)> {
    for (group_pos, esurface_option) in esurface_map.iter_enumerated() {
        if let Some(esurface_id) = esurface_option {
            return Ok((group_pos, *esurface_id));
        }
    }

    Err(eyre!(
        "No representative esurface found, esurface map corrupted"
    ))
}

/// Index in the list of all existing esurfaces, essentially a pointer to a pointer to an esurface
#[derive(
    Debug,
    From,
    Into,
    Copy,
    Clone,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
    Encode,
    Decode,
)]
pub struct ExistingEsurfaceId(usize);

impl Display for ExistingEsurfaceId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ExistingEsurfaceID({})", self.0)
    }
}

pub type ExternalShift = Vec<(EdgeIndex, i64)>;

/// add two external shifts, eliminates zero signs and sorts
pub(crate) fn add_external_shifts(lhs: &ExternalShift, rhs: &ExternalShift) -> ExternalShift {
    let mut res = lhs.clone();

    for rhs_element in rhs.iter() {
        if let Some(lhs_element) = res
            .iter_mut()
            .find(|lhs_element| rhs_element.0 == lhs_element.0)
        {
            lhs_element.1 += rhs_element.1;
        } else {
            res.push(*rhs_element)
        }
    }

    res.retain(|(_index, sign)| *sign != 0);
    res.sort_by_key(|(index, _)| *index);
    res
}

define_index!(
    pub struct RaisedEsurfaceId;
);

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct RaisedEsurfaceData {
    pub raised_groups: TiVec<RaisedEsurfaceId, RaisedEsurfaceGroup>,
    pub pass_two_evaluator: Option<Vec<GenericEvaluator>>,
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct RaisedEsurfaceGroup {
    pub esurface_ids: Vec<EsurfaceID>,
    pub max_occurence: usize,
}

impl RaisedEsurfaceGroupView for RaisedEsurfaceGroup {
    fn esurface_ids(&self) -> &[EsurfaceID] {
        &self.esurface_ids
    }

    fn max_occurrence(&self) -> usize {
        self.max_occurence
    }
}

impl RaisedEsurfaceDataView for RaisedEsurfaceData {
    fn for_each_raised_group(&self, mut f: impl FnMut(&dyn RaisedEsurfaceGroupView)) {
        for group in &self.raised_groups {
            f(group);
        }
    }
}

impl Graph {
    pub(crate) fn normalize_esurface_with_raised_edge_groups(
        esurface: &Esurface,
        raised_edges: &[Vec<EdgeIndex>],
    ) -> Esurface {
        let mut normalized = esurface.clone();
        for energy in normalized.energies.iter_mut() {
            if let Some(representative) = raised_edges
                .iter()
                .find(|group| group.contains(energy))
                .and_then(|group| group.first())
            {
                *energy = *representative;
            }
        }
        normalized.energies.sort();
        normalized
    }

    pub(crate) fn determine_raised_esurfaces_from_expression(
        &self,
        expr: &CFFExpression<OrientationID>,
    ) -> RaisedEsurfaceData {
        let raised_edges = self.get_raised_edge_groups();
        let normalized_cut_esurfaces = expr
            .surfaces
            .esurface_cache
            .iter()
            .map(|esurface| {
                Self::normalize_esurface_with_raised_edge_groups(esurface, &raised_edges)
            })
            .collect::<TiVec<EsurfaceID, _>>();

        let mut raised_groups = TiVec::<RaisedEsurfaceId, RaisedEsurfaceGroup>::new();

        for (esurface_id, normalized_cut_esurface) in normalized_cut_esurfaces.iter_enumerated() {
            let raised_esurface_group_id = raised_groups.iter_enumerated().find_map(
                |(raised_esurface_group_id, esurface_group)| {
                    if esurface_group
                        .esurface_ids
                        .iter()
                        .all(|esurface_id_in_group| {
                            normalized_cut_esurfaces[*esurface_id_in_group].energies
                                == normalized_cut_esurface.energies
                                && normalized_cut_esurfaces[*esurface_id_in_group].external_shift
                                    == normalized_cut_esurface.external_shift
                        })
                    {
                        Some(raised_esurface_group_id)
                    } else {
                        None
                    }
                },
            );

            if let Some(found_group_id) = raised_esurface_group_id {
                raised_groups[found_group_id].esurface_ids.push(esurface_id);
            } else {
                raised_groups.push(RaisedEsurfaceGroup {
                    esurface_ids: vec![esurface_id],
                    max_occurence: 0,
                });
            }
        }

        let mut result = RaisedEsurfaceData {
            raised_groups,
            pass_two_evaluator: None,
        };

        let mut expression_copy = expr.clone();
        expression_copy.normalize_wrt_all_raisings(&result);

        for cut_group in result.raised_groups.iter_mut() {
            let representative_esurface_id = cut_group.esurface_ids[0];

            let max_occurence_for_this_id = expression_copy
                .orientations
                .iter()
                .map(|orientation_expression| {
                    orientation_expression.max_effective_denominator_value_count_on_branch(
                        &crate::cff::surface::HybridSurfaceID::Esurface(representative_esurface_id),
                    )
                })
                .max()
                .unwrap_or(0);

            cut_group.max_occurence = max_occurence_for_this_id;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use linnet::half_edge::HedgeGraph;
    use linnet::half_edge::builder::HedgeGraphBuilder;
    use linnet::half_edge::involution::{EdgeIndex, Flow, Orientation};
    use linnet::half_edge::nodestore::NodeStorageVec;
    use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
    use symbolica::atom::{Atom, AtomCore};
    use symbolica::parse;

    use crate::cff::VertexSet;
    use crate::graph::{Graph, LmbIndex, LoopMomentumBasis, parse::from_dot::IntoGraph};
    use crate::initialisation::test_initialise;
    use crate::integrands::process::SamplingMapAffine;
    use crate::momentum::{
        FourMomentum, ThreeMomentum,
        sample::{ExternalFourMomenta, LoopMomenta, SubspaceData},
        signature::LoopExtSignature,
    };
    use crate::processes::CrossSectionCut;
    use crate::{
        cff::{esurface::Esurface, generation::ShiftRewrite},
        dot,
        utils::{
            ArbPrec, DEFAULT_ESURFACE_EXISTENCE_THRESHOLD, ESURFACE_SHIFT_THRESHOLD, F, FloatLike,
            QuadFloat,
            newton_solver::{
                RadialRootDiagnostics, RadialRootIdentity, SafeguardedNewtonError,
                safeguarded_newton_iteration_and_derivative,
            },
            test_utils::dummy_hedge_graph,
        },
    };
    use typed_index_collections::ti_vec;

    use super::{EsurfaceExistence, add_external_shifts};

    #[test]
    fn canonical_ray_embedding_preserves_ordered_occurrences_and_native_coefficients() {
        use super::EsurfaceRay;
        use crate::utils::SamplingFloat;

        fn check<T: FloatLike>(marker: F<T>) {
            let one = marker.one();
            let velocity = ThreeMomentum::new(marker.clone(), -marker.clone(), one.zero());
            let offset = ThreeMomentum::new(
                one.from_usize(2),
                &marker / one.from_usize(3),
                F(T::from_f64_exact_binary(-0.0)),
            );
            let repeated = (EdgeIndex(4), velocity.clone(), offset.clone(), one.clone());
            let ray = EsurfaceRay {
                energies: vec![
                    repeated.clone(),
                    (EdgeIndex(7), offset, velocity, marker.clone()),
                    repeated,
                ],
                shift: -marker.clone(),
            };
            let canonical = ray.to_arb_exact().unwrap();
            assert_eq!(
                canonical.energies.iter().map(|entry| entry.0).collect_vec(),
                [EdgeIndex(4), EdgeIndex(7), EdgeIndex(4)]
            );
            assert_eq!(canonical.energies[0], canonical.energies[2]);
            assert_eq!(canonical.energies[0].1.px, marker.to_arb_exact().unwrap());
            // Preserve the native scalar's zero convention; a Quad sum need
            // not retain the sign of the binary64 token used to construct it.
            assert_eq!(
                canonical.energies[0]
                    .2
                    .pz
                    .0
                    .mpfr_enclosure(1000)
                    .0
                    .is_sign_negative(),
                ray.energies[0].2.pz.0.into_f64().is_sign_negative()
            );
            let roundtrip = canonical.materialize::<T>().unwrap();
            assert_eq!(roundtrip.energies, ray.energies);
            assert_eq!(roundtrip.shift, ray.shift);
            for field in 0..4 {
                let mut invalid = ray.clone();
                let nonfinite = F(T::from_f64_exact_binary(f64::INFINITY));
                match field {
                    0 => invalid.energies[0].1.px = nonfinite,
                    1 => invalid.energies[0].2.py = nonfinite,
                    2 => invalid.energies[0].3 = nonfinite,
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
    fn joint_sampling_matches_distinct_reversed_edges_and_unequal_masses() {
        test_initialise().unwrap();
        // Splitting the kite's common line through E gives distinct edge IDs
        // with opposite complete affine routes, without altering its two cycles.
        let graph: Graph = dot!(digraph joint_serial_energy {
            ext_in [style=invis]
            ext_out [style=invis]
            node [num=1]
            edge [num=1 mass=1]
            ext_in -> A:0 [id=0 mass=0]
            C:1 -> ext_out [id=1 mass=0]
            A -> B [id=2 mass=2]
            B -> C [id=3]
            C -> D [id=4 lmb_id=0]
            D -> A [id=5 mass=3]
            B -> E [id=6 lmb_id=1]
            D -> E [id=7]
        })
        .unwrap();
        let lmbs = ti_vec![graph.loop_momentum_basis.clone()];
        let id = LmbIndex::from(0);
        let lmb = &lmbs[id];
        let subspace = SubspaceData::new_from_parent_basis_edges(
            &[EdgeIndex(6)],
            &graph.full_filter(),
            id,
            &graph,
            &lmbs,
        )
        .unwrap();
        let active = subspace.iter_lmb_indices().next().unwrap();
        let prior = lmb
            .loop_edges
            .iter_enumerated()
            .find_map(|(index, edge)| (*edge == EdgeIndex(4)).then_some(index))
            .unwrap();
        let left = Esurface {
            energies: vec![EdgeIndex(2), EdgeIndex(4), EdgeIndex(6)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let right = Esurface {
            energies: vec![EdgeIndex(3), EdgeIndex(5), EdgeIndex(7)],
            ..left.clone()
        };
        let masses = graph
            .underlying
            .new_edgevec_from_iter([0., 0., 2., 1., 1., 3., 1., 1.].map(F))
            .unwrap();
        let externals = ExternalFourMomenta::from_iter(
            [FourMomentum::from_args(F(26_f64.sqrt()), F(0.), F(0.), F(1.)); 2],
        );
        let (prepare, common) = left
            .sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &graph,
                &masses,
                &externals,
                &[prior],
            )
            .unwrap();
        let geometry = prepare(&[1., 0., -0.5]).unwrap();
        assert_eq!(geometry.masses, [1., 2., 3.]);
        assert_eq!(geometry.shifts, [[1., 0., 0.5], [1., 0., -0.5]]);
        assert_ne!(
            lmb.edge_signatures[EdgeIndex(6)],
            lmb.edge_signatures[EdgeIndex(7)]
        );
        assert!(lmb.edges_are_raised(EdgeIndex(6), EdgeIndex(7)));
        assert_eq!(common, lmb.edge_signatures[EdgeIndex(6)]);
        let mut loops = LoopMomenta::from_iter([ThreeMomentum::new(F(0.), F(0.), F(0.)); 2]);
        loops[prior] = ThreeMomentum::new(F(1.), F(0.), F(-0.5));
        for components in [[0.5, -1., 1.], [-0.7, 0.3, -0.4]] {
            let x = ThreeMomentum::new(F(components[0]), F(components[1]), F(components[2]));
            loops[active] = x;
            let common_energy = (x.norm_squared() + F(1.)).sqrt();
            for (i, surface) in [&left, &right].into_iter().enumerate() {
                let [a, b, c] = geometry.shifts[i];
                let partner = x + ThreeMomentum::new(F(a), F(b), F(c));
                let expected = common_energy
                    + (partner.norm_squared() + F(geometry.masses[i + 1]).square()).sqrt()
                    - F(geometry.energy_sums[i]);
                assert!(
                    (surface.compute_from_momenta(lmb, &masses, &loops, &externals) - expected)
                        .abs()
                        < F(1e-12)
                );
            }
        }
        // Fixed energy multiplicities remain part of the original equation.
        let mut repeated = left.clone();
        repeated.energies.push(EdgeIndex(4));
        let (prepare_repeated, _) = repeated
            .sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &graph,
                &masses,
                &externals,
                &[prior],
            )
            .unwrap();
        assert!(
            (prepare_repeated(&[1., 0., -0.5]).unwrap().energy_sums[0]
                - (geometry.energy_sums[0] - 1.5))
                .abs()
                < 1e-12
        );

        let mut inconsistent_masses = masses.clone();
        inconsistent_masses[EdgeIndex(7)] = F(2.);
        assert!(
            left.sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &graph,
                &inconsistent_masses,
                &externals,
                &[prior],
            )
            .err()
            .unwrap()
            .to_string()
            .contains("inconsistent native values")
        );
        // Changing only the symbolic mass distinguishes formal identity from
        // accidentally equal numeric values in a supplied mass cache.
        let mut different_mass_graph = graph.clone();
        different_mass_graph.underlying[EdgeIndex(7)].particle =
            crate::graph::edge::PossibleParticle::JustMass { expr: parse!("2") };
        assert!(
            left.sampling_joint_geometry_in_subspace(
                &right,
                &subspace,
                &lmbs,
                &different_mass_graph,
                &masses,
                &externals,
                &[prior],
            )
            .err()
            .unwrap()
            .to_string()
            .contains("found 0 candidates")
        );
    }

    #[test]
    fn joint_sampling_matches_exact_external_cm_relation_only() {
        test_initialise().unwrap();
        // The two extra external ports between the distinct common-energy
        // edges sum to zero spatially in CM, but carry nonzero total energy.
        let graph: Graph = dot!(digraph joint_external_energy {
            ext_in [style=invis]
            ext_out [style=invis]
            beam_a [style=invis]
            beam_b [style=invis]
            node [num=1]
            edge [num=1 mass=1]
            ext_in -> A:0 [id=0 mass=0]
            C:1 -> ext_out [id=1 mass=0]
            A -> B [id=2 mass=2]
            B -> C [id=3]
            C -> D [id=4 lmb_id=0]
            D -> A [id=5 mass=3]
            B -> E [id=6 lmb_id=1]
            D -> E [id=7]
            beam_a -> E [id=8 mass=0]
            beam_b -> E [id=9 mass=0]
        })
        .unwrap();
        fn check<T: FloatLike>(graph: &Graph) {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let lmbs = ti_vec![graph.loop_momentum_basis.clone()];
            let id = LmbIndex::from(0);
            let lmb = &lmbs[id];
            let subspace = SubspaceData::new_from_parent_basis_edges(
                &[EdgeIndex(6)],
                &graph.full_filter(),
                id,
                graph,
                &lmbs,
            )
            .unwrap();
            let active = subspace.iter_lmb_indices().next().unwrap();
            let prior = lmb
                .loop_edges
                .iter_enumerated()
                .find_map(|(index, edge)| (*edge == EdgeIndex(4)).then_some(index))
                .unwrap();
            let left = Esurface {
                energies: vec![EdgeIndex(2), EdgeIndex(4), EdgeIndex(6)],
                external_shift: vec![(EdgeIndex(0), -1)],
                vertex_set: VertexSet::dummy(),
            };
            let right = Esurface {
                energies: vec![EdgeIndex(3), EdgeIndex(5), EdgeIndex(7)],
                ..left.clone()
            };
            let masses = graph
                .underlying
                .new_edgevec_from_iter(
                    [0, 0, 2, 1, 1, 3, 1, 1, 0, 0].map(|mass| one.from_usize(mass)),
                )
                .unwrap();
            assert!(
                !lmb.edge_signatures[EdgeIndex(6)]
                    .equality_up_to_sign(&lmb.edge_signatures[EdgeIndex(7)])
            );
            assert!(!lmb.edges_are_raised(EdgeIndex(6), EdgeIndex(7)));
            for boosted in [false, true] {
                let externals = ExternalFourMomenta::from_iter(lmb.ext_edges.iter().map(|edge| {
                    let (energy, z) = match edge.0 {
                        0 => (20, 0),
                        1 => (26, 0),
                        8 => (3, 2),
                        9 => (3, -2),
                        _ => unreachable!(),
                    };
                    let (energy, z) = (one.from_i64(energy), one.from_i64(z));
                    // Exact rational Lorentz boost: gamma=5/4, gamma*v=3/4.
                    let (energy, z) = if boosted {
                        (
                            (one.from_i64(5) * &energy + one.from_i64(3) * &z) / one.from_i64(4),
                            (one.from_i64(3) * energy + one.from_i64(5) * z) / one.from_i64(4),
                        )
                    } else {
                        (energy, z)
                    };
                    FourMomentum::from_args(energy, zero.clone(), zero.clone(), z)
                }));
                let matched = left.sampling_joint_geometry_in_subspace(
                    &right,
                    &subspace,
                    &lmbs,
                    graph,
                    &masses,
                    &externals,
                    &[prior],
                );
                if boosted {
                    assert!(
                        matched
                            .err()
                            .unwrap()
                            .to_string()
                            .contains("found 0 candidates")
                    );
                    continue;
                }
                let (prepare, common) = matched.unwrap();
                let prior_values = [one.clone(), zero.clone(), -(&one / one.from_i64(2))];
                let geometry = prepare(&prior_values.clone().map(|x| x.0)).unwrap();
                assert_eq!(
                    geometry.masses,
                    [one.0.clone(), one.from_i64(2).0, one.from_i64(3).0]
                );
                let mut loops = LoopMomenta::from_iter(
                    (0..2).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
                );
                loops[prior] = ThreeMomentum::new(
                    prior_values[0].clone(),
                    prior_values[1].clone(),
                    prior_values[2].clone(),
                );
                let spatial: crate::momentum::sample::ExternalThreeMomenta<F<T>> =
                    externals.iter().map(|p| p.spatial.clone()).collect();
                for values in [[1, -2, 3], [-3, 1, -1]] {
                    let values = values.map(|value| one.from_i64(value) / one.from_i64(4));
                    loops[active] =
                        ThreeMomentum::new(values[0].clone(), values[1].clone(), values[2].clone());
                    let x: ThreeMomentum<F<T>> = common.compute_momentum(&loops, &spatial);
                    let e0 = (x.norm_squared() + F(geometry.masses[0].clone()).square()).sqrt();
                    for (index, surface) in [&left, &right].into_iter().enumerate() {
                        let shift = geometry.shifts[index].each_ref().map(|x| F(x.clone()));
                        let partner = &x
                            + &ThreeMomentum::new(
                                shift[0].clone(),
                                shift[1].clone(),
                                shift[2].clone(),
                            );
                        let expected = &e0
                            + (partner.norm_squared()
                                + F(geometry.masses[index + 1].clone()).square())
                            .sqrt()
                            - F(geometry.energy_sums[index].clone());
                        let original =
                            surface.compute_from_momenta(lmb, &masses, &loops, &externals);
                        assert!(
                            (original - expected).abs() < one.epsilon().sqrt() * one.from_i64(100)
                        );
                    }
                }
                // External temporal shifts and every fixed occurrence survive
                // the spatial alias; no on-cut energy substitution is made.
                let fixed: ThreeMomentum<F<T>> =
                    lmb.edge_signatures[EdgeIndex(4)].compute_momentum(&loops, &spatial);
                let fixed_energy = (fixed.norm_squared() + masses[EdgeIndex(4)].square()).sqrt();
                assert!(
                    (F(geometry.energy_sums[0].clone())
                        + fixed_energy
                        + left.compute_shift_part_from_momenta(&externals, lmb))
                    .abs()
                        < one.epsilon().sqrt() * one.from_i64(100)
                );
            }
        }
        check::<f64>(&graph);
        check::<QuadFloat>(&graph);
        check::<ArbPrec>(&graph);
    }

    #[test]
    fn routed_enclosure_preserves_original_terms_and_native_rescale() {
        use super::{EsurfaceRay, SamplingEvaluationError};
        use crate::momentum::sample::LoopIndex;
        use rug::Float;

        fn check<T: FloatLike>() {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let graph = dummy_hedge_graph(11);
            let mut routes = (0..5)
                .map(|i| {
                    let mut internal = vec![0; 5];
                    internal[i] = 1;
                    LoopExtSignature::from((internal, vec![0; 2]))
                })
                .collect::<Vec<_>>();
            routes.extend([
                LoopExtSignature::from((vec![1, 1, 1, -1, -1], vec![1, 0])),
                LoopExtSignature::from((vec![0; 5], vec![0, 1])),
                LoopExtSignature::from((vec![0; 5], vec![0, 0])),
                LoopExtSignature::from((vec![0; 5], vec![1, -1])),
                LoopExtSignature::from((vec![0; 5], vec![1, 0])),
                LoopExtSignature::from((vec![0; 5], vec![0, 1])),
            ]);
            let lmb = LoopMomentumBasis {
                tree: SuBitGraph::empty(0),
                loop_edges: (0..5).map(EdgeIndex).collect(),
                ext_edges: vec![EdgeIndex(9), EdgeIndex(10)].into(),
                edge_signatures: graph.new_edgevec_from_iter(routes).unwrap(),
            };
            let big = one.clone() / one.epsilon().square();
            let middle = big.sqrt();
            let vector = |x| ThreeMomentum::new(x, zero.clone(), zero.clone());
            // The complete signed route is exactly one. Three separated
            // scales also expose cancellation in a two-limb Quad sum.
            let loops = LoopMomenta::from_iter(
                [big.clone(), middle.clone(), one.clone(), middle, big].map(&vector),
            );
            let externals = ExternalFourMomenta::from_iter([
                FourMomentum::from_args(
                    one.from_i64(7),
                    one.from_i64(3),
                    zero.clone(),
                    zero.clone(),
                ),
                FourMomentum::from_args(
                    one.from_i64(11),
                    one.from_i64(4),
                    zero.clone(),
                    zero.clone(),
                ),
            ]);
            let mut masses = graph
                .new_edgevec_from_iter(vec![one.from_i64(3); 11])
                .unwrap();
            masses[EdgeIndex(7)] = zero.clone();
            let surface = Esurface {
                energies: vec![EdgeIndex(5), EdgeIndex(5), EdgeIndex(6), EdgeIndex(7)],
                external_shift: vec![(EdgeIndex(8), 3), (EdgeIndex(6), -1)],
                vertex_set: VertexSet::dummy(),
            };
            // Two sqrt(4^2+3^2), one fixed sqrt(4^2+3^2), one exact
            // massless zero, and shift 3*(7-11)-11 give exactly -8.
            let bounds = surface
                .evaluate_routed_enclosed(&one, &loops, &externals, &masses, &lmb)
                .unwrap();
            assert_eq!(
                bounds,
                [Float::with_val(2048, -8), Float::with_val(2048, -8)]
            );
            let center = LoopMomenta::from_iter((0..5).map(|_| vector(zero.clone())));
            let represented = surface.routed_ray(&loops, &center, &externals, &masses, &lmb);
            assert_eq!(represented.energies[0].1.px, zero);
            let old_value = represented.evaluate(&one).0;
            assert!(
                old_value < one.from_i64(-9),
                "rounded routing must fail the exact -8 oracle"
            );

            // Include the physical scalar multiplication itself: evaluating
            // exact tau*K is different from evaluating the rounded rescale.
            let tau = &one + one.epsilon();
            let mut loops = center;
            loops[LoopIndex(0)] = vector(tau.clone());
            let scalar = Esurface {
                energies: vec![EdgeIndex(0)],
                external_shift: vec![],
                vertex_set: VertexSet::dummy(),
            };
            masses[EdgeIndex(0)] = zero.clone();
            let mathematical = scalar
                .evaluate_routed_enclosed(&tau, &loops, &externals, &masses, &lmb)
                .unwrap();
            let completed = loops.rescale(&tau, None);
            let actual = scalar
                .evaluate_routed_enclosed(&one, &completed, &externals, &masses, &lmb)
                .unwrap();
            assert!(
                mathematical[0] > actual[1],
                "the omitted epsilon^2 term must be resolved"
            );
            let actual_native = EsurfaceRay::<T>::enclosed_native(&completed[LoopIndex(0)].px);
            assert!(actual[0] <= actual_native[0] && actual[1] >= actual_native[1]);

            let missing = LoopMomenta::from_iter(loops.iter().take(4).cloned());
            assert!(
                scalar
                    .evaluate_routed_enclosed(&one, &missing, &externals, &masses, &lmb)
                    .is_err()
            );
            loops[LoopIndex(0)].px = zero.clone() / zero.clone();
            assert!(matches!(
                scalar
                    .evaluate_routed_enclosed(&one, &loops, &externals, &masses, &lmb)
                    .unwrap_err()
                    .downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::Unrepresentable { .. })
            ));
        }
        check::<f64>();
        check::<QuadFloat>();
        check::<ArbPrec>();
    }

    #[test]
    fn hosted_normal_alignment_enforces_directed_half_budgets() {
        use super::{EsurfaceRay, SamplingEvaluationError};
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
        EsurfaceRay::<f64>::verify_normal_alignment(
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
                EsurfaceRay::<f64>::verify_normal_alignment(canonical.clone(), native, host, 0.125)
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
                EsurfaceRay::<f64>::verify_normal_alignment(
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
            EsurfaceRay::<f64>::verify_normal_alignment(
                canonical.clone(),
                canonical.clone(),
                point(0),
                0.0
            )
            .is_err()
        );
        assert!(matches!(
            EsurfaceRay::<f64>::verify_normal_alignment(
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
        use super::EsurfaceRay;
        use crate::{
            integrands::process::sampling_maps::SamplingEvaluationError,
            momentum::{Rotatable, Rotation, RotationMethod},
        };
        fn check<T: FloatLike>() {
            let one = F::<T>::default().one();
            let f = |n| one.from_usize(n);
            let vector = |x| ThreeMomentum::new(x, one.zero(), one.zero());
            // eta(t)=sqrt((3t+1)^2+9)-6, whose positive root is
            // (sqrt(27)-1)/3. The independent radical is not a Newton result.
            let ray = EsurfaceRay {
                energies: vec![(EdgeIndex(1), vector(f(3)), vector(one.clone()), f(3))],
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
            for (_, v, b, mass) in &mut shallow.energies {
                *v = &*v * &tiny;
                *b = &*b * &tiny;
                *mass *= &tiny;
            }
            shallow.shift *= &tiny;
            let mut shallow_solution = solution.clone();
            shallow_solution.derivative_at_solution *= &tiny;
            shallow_solution.error_of_function *= &tiny;
            shallow
                .verify_lu_candidate(&shallow, &shallow_solution, 6, &budget)
                .unwrap();
            let mut tiny_energy = ray.clone();
            tiny_energy
                .energies
                .push((EdgeIndex(2), vector(one.zero()), vector(one.zero()), tiny));
            tiny_energy
                .verify_lu_candidate(&tiny_energy, &solution, 6, &budget)
                .unwrap();
            tiny_energy.energies[1].3 = one.zero();
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
            wrong.energies[0].0 = EdgeIndex(2);
            assert!(
                ray.verify_lu_candidate(&wrong, &solution, 6, &budget)
                    .unwrap_err()
                    .to_string()
                    .contains("ordered energy occurrences")
            );
            wrong = ray.clone();
            wrong.energies[0].3 = one.zero() / one.zero();
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
            completed.energies[0].1.px += &cancellation - &one;
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
    fn radial_guesses_and_lu_roots_handle_constant_massive_energies() {
        test_initialise().unwrap();
        // A ttH cut with back-to-back tops and a Higgs at rest has the analytic root
        // sqrt(((Q-mH)/2)^2-mt^2)/|p|. The Higgs momentum vanishes by cancellation of
        // two nonzero basis momenta, not because its graph signature is constant.
        let graph: Graph = dot!(digraph massive_energy_at_rest {
            ext [style=invis]
            node [num=1]
            edge [num=1 mass=0]
            ext -> a:0 [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2 lmb_id=1]
            a -> b [id=3]
            b:1 -> ext [id=4]
        })
        .unwrap();
        let lmbs = ti_vec![graph.loop_momentum_basis.clone()];
        let lmb = &lmbs[LmbIndex::from(0)];
        let subspace = SubspaceData::new_from_parent_basis_edges(
            &[EdgeIndex(1), EdgeIndex(2)],
            &graph.underlying.full_filter(),
            LmbIndex::from(0),
            &graph,
            &lmbs,
        )
        .unwrap();
        let surface = Esurface {
            energies: vec![EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = graph
            .underlying
            .new_edgevec_from_iter([F(0.0), F(173.0), F(173.0), F(125.0), F(0.0)])
            .unwrap();
        // Supply both external ports as future-directed momenta; the graph's incoming and
        // outgoing edge flows provide their relative sign in momentum conservation.
        let externals = ExternalFourMomenta::from_iter(
            [FourMomentum::from_args(F(1000.0), F(0.0), F(0.0), F(0.0)); 2],
        );
        let center = LoopMomenta::from_iter([ThreeMomentum::new(F(0.0), F(0.0), F(0.0)); 2]);

        for momentum in [0.0, 100.0] {
            let loops = LoopMomenta::from_iter([
                ThreeMomentum::new(F(momentum), F(0.0), F(0.0)),
                ThreeMomentum::new(F(-momentum), F(0.0), F(0.0)),
            ]);
            let guesses = [
                surface.get_radius_guess(&loops, &externals, lmb),
                surface.get_radius_guess_subspace(
                    &loops, &externals, &subspace, &lmbs, &graph, &masses,
                ),
            ];
            for (positive, negative) in guesses {
                assert!(positive.0.is_finite());
                assert_eq!(negative, -positive);
                let result = safeguarded_newton_iteration_and_derivative(
                    &F(0.0),
                    &positive,
                    |radius| {
                        surface.compute_self_and_r_derivative(
                            radius, &loops, &center, &externals, &masses, lmb,
                        )
                    },
                    &F(1.0),
                    2000,
                    64,
                    &F(1000.0),
                );
                if momentum == 0.0 {
                    assert_eq!(positive, F(0.0));
                    // The constant surface is negative everywhere: no finite radial root
                    // exists. Reject the ray instead of propagating a NaN or inventing a root.
                    assert!(matches!(
                        result,
                        Err(SafeguardedNewtonError::InvalidOutside { .. })
                    ));
                } else {
                    let result = result.unwrap();
                    let top_energy = (1000.0_f64 - 125.0) / 2.0;
                    let expected = (top_energy.powi(2) - 173.0_f64.powi(2)).sqrt() / momentum;
                    assert!((result.solution.0 - expected).abs() < 1.0e-14);
                    assert!(result.derivative_at_solution.0 > 0.0);
                    assert!(result.error_of_function.0.abs() < 1.0e-12);
                }
            }
        }

        fn check_native<T: FloatLike>(graph: &Graph, surface: &Esurface) {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let f = |value| one.from_i64(value);
            let lmb = &graph.loop_momentum_basis;
            let mut masses = graph
                .underlying
                .new_edgevec_from_iter([0, 173, 173, 125, 0].map(f))
                .unwrap();
            let e_cm = f(1000);
            let externals = ExternalFourMomenta::from_iter((0..2).map(|_| {
                FourMomentum::from_args(e_cm.clone(), zero.clone(), zero.clone(), zero.clone())
            }));
            let zero_vector = ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone());
            let center = LoopMomenta::from_iter([zero_vector.clone(), zero_vector.clone()]);
            let perturbation = &one / f(10).powi(25);
            let p = f(100) + &perturbation;
            let loops = LoopMomenta::from_iter([
                ThreeMomentum::new(p.clone(), zero.clone(), zero.clone()),
                ThreeMomentum::new(-&p, zero.clone(), zero.clone()),
            ]);
            let identity = RadialRootIdentity::new("native shared LU host".into());
            let mut diagnostics = RadialRootDiagnostics::default();
            let mut prior_policy = diagnostics.clone();
            let (ray, result) = surface
                .solve_lu_cut(
                    &loops,
                    &externals,
                    &masses,
                    lmb,
                    &e_cm,
                    &mut diagnostics,
                    &identity,
                )
                .unwrap();
            // The extracted entry must preserve the existing physical policy on
            // identical routed inputs and histories, including diagnostic data.
            let guess = surface.get_radius_guess(&loops, &externals, lmb).0;
            assert_eq!(
                Esurface::radius_guess_from_terms(
                    &ray.shift,
                    ray.energies
                        .iter()
                        .map(|(_, v, b, _)| (v.norm_squared(), v.clone() * b)),
                ),
                guess
            );
            let prior = prior_policy
                .solve(
                    &identity,
                    &zero,
                    &guess,
                    |t| {
                        surface.compute_self_and_r_derivative(
                            t, &loops, &center, &externals, &masses, lmb,
                        )
                    },
                    &one,
                    2000,
                    64,
                    &e_cm,
                )
                .unwrap();
            assert_eq!(result.solution, prior.solution);
            assert_eq!(result.derivative_at_solution, prior.derivative_at_solution);
            assert_eq!(result.error_of_function, prior.error_of_function);
            assert_eq!(result.num_iterations_used, prior.num_iterations_used);
            let top_energy = (&e_cm - f(125)) / f(2);
            let radial_momentum = (top_energy.square() - f(173).square()).sqrt();
            assert!((&result.solution - &radial_momentum / &p).abs() < one.epsilon() * f(128));
            if perturbation > one.epsilon() * f(1000) {
                assert!(
                    (&result.solution - radial_momentum / f(100)).abs() > &one / f(10).powi(28)
                );
            }

            // A retained ray preserves energy multiplicity and a nonzero affine
            // center. Independently eta(t)=2*sqrt((3*t+1)^2+9)+5-e_cm.
            use crate::utils::hyperdual_utils::{
                extract_t_derivatives, new_constant, simple_n_deriv_shape,
            };
            use symbolica::domains::dual::HyperDual;
            let mut jet_surface = surface.clone();
            jet_surface.energies = vec![EdgeIndex(1), EdgeIndex(1), EdgeIndex(3)];
            let mut jet_masses = masses.clone();
            jet_masses[EdgeIndex(1)] = f(3);
            jet_masses[EdgeIndex(3)] = f(4);
            let jet_velocity = LoopMomenta::from_iter([
                ThreeMomentum::new(f(3), zero.clone(), zero.clone()),
                ThreeMomentum::new(-f(3), zero.clone(), zero.clone()),
            ]);
            let jet_center = LoopMomenta::from_iter([
                ThreeMomentum::new(one.clone(), zero.clone(), zero.clone()),
                ThreeMomentum::new(-&one, zero.clone(), zero.clone()),
            ]);
            let jet_externals =
                ExternalFourMomenta::from_iter((0..2).map(|_| {
                    FourMomentum::from_args(e_cm.clone(), zero.clone(), f(3), zero.clone())
                }));
            let jet_ray = jet_surface.routed_ray(
                &jet_velocity,
                &jet_center,
                &jet_externals,
                &jet_masses,
                lmb,
            );
            assert_eq!(
                jet_ray.energies.iter().map(|term| term.0).collect_vec(),
                jet_surface.energies
            );
            assert_eq!(jet_ray.energies[2].1.norm_squared(), zero);
            assert_eq!(jet_ray.energies[2].2.norm_squared(), f(9));
            // A CT seed still ignores the center, whereas the explicit LU ray
            // seed consumes its actual affine offsets. Here the difference is2/3.
            let old_seed = jet_surface
                .get_radius_guess(&jet_velocity, &jet_externals, lmb)
                .0;
            assert_eq!(old_seed, &e_cm / f(6));
            let ray_seed = Esurface::radius_guess_from_terms(
                &jet_ray.shift,
                jet_ray
                    .energies
                    .iter()
                    .map(|(_, v, b, _)| (v.norm_squared(), v.clone() * b)),
            );
            assert!((ray_seed - old_seed - f(2) / f(3)).abs() < one.epsilon() * f(1024));
            let expected = [
                f(15) - &e_cm,
                f(24) / f(5),
                f(162) / f(125),
                -f(5832) / f(3125),
            ];
            let scalar = jet_ray.evaluate(&one);
            assert_eq!(scalar.0, expected[0]);
            assert!((&scalar.1 - &expected[1]).abs() < one.epsilon() * f(128));
            for order in 1..=3 {
                let t =
                    HyperDual::<F<T>>::new(simple_n_deriv_shape(order)).variable(0, one.clone());
                let actual = extract_t_derivatives(jet_ray.evaluate_dual(&t));
                // Independent legacy path: first scale complete loop jets, then
                // route original Esurface energies and true external ports.
                let dual_loops = LoopMomenta::from_iter(
                    jet_velocity.iter().zip(jet_center.iter()).map(|(v, b)| {
                        v.map_ref(&|x| new_constant(&t, x) * &t)
                            + b.map_ref(&|x| new_constant(&t, x))
                    }),
                );
                let dual_externals = jet_externals
                    .iter()
                    .map(|p| p.map_ref(&|x| new_constant(&t, x)))
                    .collect();
                let original = extract_t_derivatives(jet_surface.compute_from_dual_momenta(
                    lmb,
                    &jet_masses,
                    &dual_loops,
                    &dual_externals,
                ));
                for ((actual, original), expected) in actual.iter().zip(&original).zip(&expected) {
                    let tolerance = one.epsilon() * f(2048) * (one.clone() + expected.abs());
                    assert!((actual - expected).abs() < tolerance);
                    assert!((actual - original).abs() < tolerance);
                }
            }

            // Reassociation is deliberately confined to explicit LU preparation.
            // With exact binary64 t=0.1 and k=[10^16,-10^16+2], the old
            // scale-then-route momentum is 1/8; route-then-scale gives 2*t.
            // This is a negative equivalence control, not a tolerance failure.
            let mut cancellation = surface.clone();
            cancellation.energies = vec![EdgeIndex(3)];
            cancellation.external_shift.clear();
            let large = f(10).powi(16);
            let velocity = LoopMomenta::from_iter([
                ThreeMomentum::new(large.clone(), zero.clone(), zero.clone()),
                ThreeMomentum::new(-large + f(2), zero.clone(), zero.clone()),
            ]);
            let t = F(T::from_f64_exact_binary(0.1));
            for (offset, mass) in [(0, 0), (3, 0), (3, 4)] {
                let center = LoopMomenta::from_iter([
                    ThreeMomentum::new(f(offset), zero.clone(), zero.clone()),
                    zero_vector.clone(),
                ]);
                let mut case_masses = masses.clone();
                case_masses[EdgeIndex(3)] = f(mass);
                let represented =
                    cancellation.routed_ray(&velocity, &center, &externals, &case_masses, lmb);
                let prepared = represented.evaluate(&t);
                let original = cancellation.compute_self_and_r_derivative(
                    &t,
                    &velocity,
                    &center,
                    &externals,
                    &case_masses,
                    lmb,
                );
                let prepared_p = f(2) * &t + f(offset);
                let prepared_e = (prepared_p.square() + f(mass).square()).sqrt();
                let tolerance = one.epsilon() * f(128);
                assert!((&prepared.0 - &prepared_e).abs() < tolerance);
                assert!((&prepared.1 - f(2) * &prepared_p / &prepared_e).abs() < tolerance);
                if one.epsilon() > &one / f(10).powi(20) {
                    let original_p = &one / f(8) + f(offset);
                    let original_e = (original_p.square() + f(mass).square()).sqrt();
                    assert!((&original.0 - &original_e).abs() < tolerance);
                    assert!((&original.1 - f(2) * original_p / original_e).abs() < tolerance);
                    assert!((&prepared.0 - &original.0).abs() > &one / f(100));
                } else {
                    // Quad/Arb have enough mantissa for these exact represented
                    // products; they recover the same point by either ordering.
                    assert!((&original.0 - &prepared.0).abs() < tolerance);
                    assert!((&original.1 - &prepared.1).abs() < tolerance);
                }
            }

            // At the origin combine a massless moving energy with a shifted
            // massive energy E=sqrt(4^2+3^2). Their right derivative is exactly5.
            let mut endpoint = surface.clone();
            endpoint.energies = vec![EdgeIndex(1), EdgeIndex(2)];
            masses[EdgeIndex(1)] = zero.clone();
            masses[EdgeIndex(2)] = f(3);
            let velocity = LoopMomenta::from_iter([
                ThreeMomentum::new(f(3), f(4), zero.clone()),
                zero_vector.clone(),
            ]);
            let mut shifted_center = LoopMomenta::from_iter([
                zero_vector.clone(),
                ThreeMomentum::new(f(4), zero.clone(), zero.clone()),
            ]);
            let (value, derivative) = endpoint.compute_self_and_r_derivative(
                &zero,
                &velocity,
                &shifted_center,
                &externals,
                &masses,
                lmb,
            );
            assert_eq!(value, f(5) - &e_cm);
            assert_eq!(derivative, f(5));
            let prepared_endpoint = endpoint
                .routed_ray(&velocity, &shifted_center, &externals, &masses, lmb)
                .evaluate(&zero);
            assert_eq!(prepared_endpoint, (f(5) - &e_cm, f(5)));
            let stationary = LoopMomenta::from_iter([zero_vector.clone(), zero_vector]);
            assert_eq!(
                endpoint
                    .compute_self_and_r_derivative(
                        &zero,
                        &stationary,
                        &shifted_center,
                        &externals,
                        &masses,
                        lmb,
                    )
                    .1,
                zero
            );
            assert!(matches!(
                endpoint.solve_lu_cut(
                    &stationary,
                    &externals,
                    &masses,
                    lmb,
                    &e_cm,
                    &mut diagnostics,
                    &identity,
                ),
                Err(SafeguardedNewtonError::InvalidOutside { .. })
            ));

            // Squaring a nonzero source can underflow in Double/Quad. It must
            // never activate the exact massless/zero-momentum endpoint rule.
            let tiny = &one / f(10).powi(200);
            masses[EdgeIndex(1)] = tiny.clone();
            let derivative = endpoint
                .compute_self_and_r_derivative(
                    &zero,
                    &velocity,
                    &shifted_center,
                    &externals,
                    &masses,
                    lmb,
                )
                .1;
            let prepared_derivative = endpoint
                .routed_ray(&velocity, &shifted_center, &externals, &masses, lmb)
                .evaluate(&zero)
                .1;
            for derivative in [derivative, prepared_derivative] {
                if tiny.square() == zero {
                    assert!(derivative.is_nan() || derivative.is_infinite());
                } else {
                    assert_eq!(derivative, zero);
                }
            }
            masses[EdgeIndex(1)] = zero.clone();
            for (axis, expected) in [3, 4, 0].into_iter().enumerate() {
                shifted_center[crate::momentum::sample::LoopIndex(0)] = ThreeMomentum::new(
                    if axis == 0 {
                        tiny.clone()
                    } else {
                        zero.clone()
                    },
                    if axis == 1 {
                        tiny.clone()
                    } else {
                        zero.clone()
                    },
                    if axis == 2 {
                        tiny.clone()
                    } else {
                        zero.clone()
                    },
                );
                let derivative = endpoint
                    .compute_self_and_r_derivative(
                        &zero,
                        &velocity,
                        &shifted_center,
                        &externals,
                        &masses,
                        lmb,
                    )
                    .1;
                let prepared_derivative = endpoint
                    .routed_ray(&velocity, &shifted_center, &externals, &masses, lmb)
                    .evaluate(&zero)
                    .1;
                for derivative in [derivative, prepared_derivative] {
                    if tiny.square() == zero {
                        assert!(derivative.is_nan() || derivative.is_infinite());
                    } else {
                        assert!((derivative - f(expected)).abs() < one.epsilon() * f(128));
                    }
                }
            }
        }
        check_native::<f64>(&graph, &surface);
        check_native::<QuadFloat>(&graph, &surface);
        check_native::<ArbPrec>(&graph, &surface);
    }

    #[test]
    fn canonical_physical_lu_root_resolves_retained_bracket_endpoints() {
        use crate::graph::FeynmanGraph;

        test_initialise().unwrap();
        // Baseline seed20011, worker3/local239 (global5156), source path[0,5].
        // Cut group3 is physical Cut0, with its original ordered four energies.
        // These are the captured Arb generation momenta, without the radial
        // normalization used by the separate foreign-inverse fixture below.
        let point_tokens = [
            "-4877.19776949134799589974803772414639433359513358235797727662621078949562551592647720875226397943918395045206756176248678175688837208619050791653582445021695271333556149595551567501518142933108041137862389137548159866615474120165599255689473831877256382039353906809830046294831502842769433115811496491403",
            "-4788.44733559297503338424048189604502643310578502030319826258363405583106128848247171958913511942781252928660066253754883065799513281132098117063664934324987293459853837736080937683626624511055230677016357961716720986212711446609804288980413634695779067271806710579027258896440419592959858535084502130291",
            "3449.85408446744743054600317691517351068755983604965291557283080813569562933280322362036960861849284721431171841172810290687504905161462871600062758043010601204035005302723970251590824854187235628530596853916572662417827550093988644061689354024522549573776864100886485399680047404526188593294956650329940",
            "70.3244241810218751702648147574915682722146332904085090192206088762519914112850473659048343153323638641871729435254245633882918549420905038752476781927385371433230303313665642902999745355570884098566729530821574404942919072073238328487435044514999060063382762063688077567584466375120955594649302749638168",
            "-64.7328550651196388632156519795545666707240287311537237777794147531761978347671932403072110940638881039840834582472804137358931989273824856974894575353144549925155422489503341679707811903159498611335677205735964537298380313103709897338608109368125932232793979426015474979910153144543330056886932939528991",
            "-169.977204441411724276834141865099406683621147453741490083449379145954509733051600240055550798830052895001473550375457779362103819932642751869589067397325831680056922569024797824237561594653585076675793132823695793970660453958165008987223868047366419040597994417711631570323172275534429627197335733541788",
            "-5600.61087887521441317726380154507077804125909265810807210495844671671183617788774688843570145751202477847145016280355634030258604974973758865789530484255123560331590332195975778091351865880370425224431732203854524561456548571792312894398114116052933817481418085818812657379218499842947239271592857749704",
            "-4836.61374618501855280967146959499945066804232727342341015614259662078809271214898964547576447118364148478316361672773076137899768111124298404537445866753685338684537418353640716742530483139927834572180382685089716954278888134794244645196893419036539682711826176677044169785265849490484473737712006654560",
            "2145.99486696308129706260024015689439317496385917360897844622508358240397624209431777021609681539543126318894513773663359142436896535260327930425968202779550650784950814483566205233631415943233817778512294707166462039932281821355935955750715055287115993275252104683925707031965340569344306593459551977201",
            "-5604.88372374262804889294640877016787850008172170528245661491078660880869657816973625663719378778526660563427937564485272926173353092285235941079848809391075384808960803275406005472826675962967469084057925785326996659483564308288383384251734955274712244877309369072980850467041673500458058150368522257495",
            "-4334.45541003902767715452592577800941446238649756064773722205221184221360633788778921553682091482858995046598733903274399368838883790258941797495609397561874360921914142839313072883007139115101612942134481476731356575319298406454126144036222827728908340679754158679444863572311869978356712039675439410147",
            "2429.79731946564390222685915546149324045671034707482980645420787377878023354031617617014193097522475262733934128132427049128461231717208810751567187554211866137651923930428300030230687501307593042508234180738169107139713264377038173966078706808479870077056431717387171214082081679052824398764469662844001",
        ];
        let point = point_tokens
            .iter()
            .map(|token| {
                let value = token.parse::<ArbPrec>().unwrap();
                assert_eq!(value.to_string(), *token);
                F(value)
            })
            .collect_vec();
        let one = point[0].one();
        let zero = one.zero();
        let model = crate::utils::load_generic_model("sm");
        let graph: Graph = include_str!("../../../../tests/resources/graphs/GL638.dot")
            .into_graph(&model)
            .unwrap();
        let lmb = &graph.loop_momentum_basis;
        assert_eq!(
            lmb.loop_edges.raw,
            vec![EdgeIndex(3), EdgeIndex(4), EdgeIndex(7), EdgeIndex(10)]
        );
        let surface = Esurface {
            energies: vec![EdgeIndex(2), EdgeIndex(6), EdgeIndex(12), EdgeIndex(13)],
            external_shift: vec![(EdgeIndex(0), -1), (EdgeIndex(1), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = graph.get_real_mass_vector::<ArbPrec>(&model);
        for (edge, mass) in [(2, 125), (6, 173), (12, 173), (13, 0)] {
            assert_eq!(masses[EdgeIndex(edge)], one.from_i64(mass));
        }
        let externals = ExternalFourMomenta::from_iter([1, -1].map(|sign| {
            FourMomentum::from_args(
                one.from_i64(500),
                zero.clone(),
                zero.clone(),
                one.from_i64(500 * sign),
            )
        }));
        assert_eq!(externals.len(), lmb.ext_edges.len());
        let momentum = LoopMomenta::from_iter(
            point
                .chunks_exact(3)
                .map(|v| ThreeMomentum::new(v[0].clone(), v[1].clone(), v[2].clone())),
        );
        let center = LoopMomenta::from_iter(
            (0..4).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
        );
        let e_cm = one.from_i64(1000);
        let identity =
            RadialRootIdentity::new("canonical physical overlap graph 'GL638' cut group 3".into());
        let ray = surface.routed_ray(&momentum, &center, &externals, &masses, lmb);
        let guess = Esurface::radius_guess_from_terms(
            &ray.shift,
            ray.energies
                .iter()
                .map(|(_, v, b, _)| (v.norm_squared(), v.clone() * b)),
        );
        let strict = safeguarded_newton_iteration_and_derivative(
            &zero,
            &guess,
            |r| ray.evaluate(r),
            &one,
            2000,
            64,
            &e_cm,
        )
        .unwrap_err();
        assert_eq!(
            format!(
                "sampling root is not certified at the current precision: {identity}: {strict:?}"
            ),
            "sampling root is not certified at the current precision: canonical physical overlap graph 'GL638' cut group 3: DidNotConverge { result: NewtonIterationResult { solution: F(VarFloat { float: 4.50096019182270793871156294720615053060112718361364618170932397734514711080987392331940322755100771019674059655008057236752740928154857006112161128562114855544422368074316508870600192219593565649427544978773164779121432253631794038806540811255813327701451237828689832176990555566921873234752788019943785e-2 }), derivative_at_solution: F(VarFloat { float: 15947.7515847915197382229906009716238684820639052766863102939247189292647433796806064815357943825914286672680520425869860307247589322869726693916868899181132224865883496823643904736774712938449859097464691001739892743096176224973672876047140085992981227937651521443373989485856759341750304597183394266484 }), error_of_function: F(VarFloat { float: -2.86698583604188839625755508139156634506370492325388705163790645185321035051758274904404068484469074644867580906595607334158997766559958751695315195045047528315210956637693486164171287342869042640210070772612642749923485444545310860535350764889727519890646541984079814999901126931044646673003826046088028e-298 }), num_iterations_used: 10 }, lower_bound: F(VarFloat { float: 4.50096019182270793871156294720615053060112718361364618170932397734514711080987392331940322755100771019674059655008057236752740928154857006112161128562114855544422368074316508870600192219593565649427544978773164779121432253631794038806540811255813327701451237828689832176990555566921873234752788019943785e-2 }), upper_bound: F(VarFloat { float: 4.50096019182270793871156294720615053060112718361364618170932397734514711080987392331940322755100771019674059655008057236752740928154857006112161128562114855544422368074316508870600192219593565649427544978773164779121432253631794038806540811255813327701451237828689832176990555566921873234752788019943844e-2 }) }"
        );
        let (_, production) = surface
            .solve_lu_cut(
                &momentum,
                &externals,
                &masses,
                lmb,
                &e_cm,
                &mut RadialRootDiagnostics::default(),
                &identity,
            )
            .unwrap();
        // The production LU caller includes the existing energy-sum forward-error
        // budget. Preserve the retained strict fixture below at its unit budget.
        let production_budget =
            one.from_i64((4 * ray.energies.len() + 1) as i64) * one.epsilon() * &e_cm;
        let (production_value, production_derivative) = ray.evaluate(&production.solution);
        assert_eq!(production.error_of_function, production_value);
        assert_eq!(production.derivative_at_solution, production_derivative);
        assert!(production.error_of_function.abs() <= production_budget);
        let certified = RadialRootDiagnostics::default()
            .solve(
                &identity,
                &zero,
                &guess,
                |r| ray.evaluate(r),
                &one,
                2000,
                64,
                &e_cm,
            )
            .unwrap();
        let SafeguardedNewtonError::DidNotConverge {
            result,
            lower_bound,
            upper_bound,
        } = strict
        else {
            panic!("retained physical source did not reproduce an exhausted bracket");
        };
        assert_eq!(result.solution, lower_bound);
        assert_eq!(result.num_iterations_used, 10);
        assert!(upper_bound > lower_bound);
        let midpoint = (&lower_bound + &upper_bound) / one.from_i64(2);
        assert!(midpoint == lower_bound || midpoint == upper_bound);
        let width = &upper_bound - &lower_bound;
        let budget = one.epsilon() * &e_cm;
        assert_eq!(certified.solution, upper_bound);
        assert_eq!(certified.num_iterations_used, result.num_iterations_used);
        let (upper_value, upper_derivative) = ray.evaluate(&upper_bound);
        assert_eq!(certified.error_of_function, upper_value);
        assert_eq!(certified.derivative_at_solution, upper_derivative);
        // Strict native convergence still fails. Only the already represented
        // upper endpoint satisfies the unchanged bracket-resolution allowance.
        assert!(
            result.error_of_function.abs() > &budget + result.derivative_at_solution.abs() * &width
        );
        assert!(certified.error_of_function.abs() > budget);
        assert!(
            certified.error_of_function.abs()
                <= &budget + certified.derivative_at_solution.abs() * &width
        );
        crate::debug_tags!(#integration, #cut, #solver;
            stage = "canonical_physical_lu_root_setup",
            cut_group = 3, physical_cut = 0,
            initial_guess = %guess, original_residual_budget = %budget,
            bracket_width = %width, file.ray = ?ray,
            file.original_momenta = ?momentum, file.original_result = ?result,
            "actual captured physical ray and unchanged LU solver policy"
        );
        for (endpoint, radius) in [("lower", &lower_bound), ("upper", &upper_bound)] {
            let (value, derivative) = ray.evaluate(radius);
            let enclosed = surface
                .evaluate_routed_enclosed(radius, &momentum, &externals, &masses, lmb)
                .unwrap();
            assert!(!derivative.is_nan() && !derivative.is_infinite() && derivative > zero);
            assert!(enclosed.iter().all(|value| value.is_finite()));
            let resolution_budget = &budget + derivative.abs() * &width;
            let allowed = budget.0.mpfr_enclosure(2048).0;
            let native_strict = value.abs() <= budget;
            let native_resolution = value.abs() <= resolution_budget;
            let original_strict = enclosed[0] >= -allowed.clone() && enclosed[1] <= allowed;
            assert!(
                original_strict,
                "{endpoint} violates the original equation budget"
            );
            crate::debug_tags!(#integration, #cut, #solver;
                stage = "canonical_physical_lu_root_endpoint", endpoint,
                radius = %radius, native_residual = %value, derivative = %derivative,
                original_residual_lower = %enclosed[0], original_residual_upper = %enclosed[1],
                original_residual_budget = %budget, resolution_budget = %resolution_budget,
                native_strict, native_resolution, original_strict,
                "both represented endpoints satisfy the unchanged original equation budget"
            );
        }

        // This adversarial callback preserves the actual energy values but deliberately
        // supplies the wrong derivative only at the alternate endpoint. Its residual
        // fits the resolution allowance; its own four-probe slope check must reject it.
        let inconsistent_derivative = |radius: &F<ArbPrec>| {
            let (value, derivative) = ray.evaluate(radius);
            let derivative = if radius == &upper_bound {
                derivative * one.from_i64(8)
            } else {
                derivative
            };
            (value, derivative)
        };
        let invalid_strict = safeguarded_newton_iteration_and_derivative(
            &zero,
            &guess,
            inconsistent_derivative,
            &one,
            2000,
            64,
            &e_cm,
        )
        .unwrap_err();
        let SafeguardedNewtonError::DidNotConverge {
            result: invalid_result,
            lower_bound: invalid_lower,
            upper_bound: invalid_upper,
        } = &invalid_strict
        else {
            panic!("adversarial callback changed the structural failure");
        };
        assert!(invalid_result.num_iterations_used >= 4);
        assert_eq!(invalid_result.solution, lower_bound);
        assert_eq!(*invalid_lower, lower_bound);
        assert_eq!(*invalid_upper, upper_bound);
        let (alternate_value, wrong_derivative) = inconsistent_derivative(&upper_bound);
        assert!(alternate_value.abs() <= &budget + wrong_derivative * &width);
        let rejected = RadialRootDiagnostics::default()
            .solve(
                &RadialRootIdentity::new("adversarial alternate endpoint derivative".into()),
                &zero,
                &guess,
                inconsistent_derivative,
                &one,
                2000,
                64,
                &e_cm,
            )
            .unwrap_err();
        assert_eq!(format!("{rejected:?}"), format!("{invalid_strict:?}"));
    }

    #[test]
    fn phase_space_cut_root_replays_canonical_foreign_inverse_cancellation() {
        test_initialise().unwrap();
        // GL638 Halton3363: the hosted joint's ordinary fallback generated this
        // canonical Arb point; its foreign full-frame Cut1 inverse failed before
        // any reference/physical body. Retain the native decimal roundtrip, not
        // a binary64 reconstruction or an independently regenerated joint point.
        let point_tokens = [
            "1367030.43662744423077657886622585118732404500993473554913895566241589500405947461690714412624268513505631015396234810943269596104532748643105542702039862197867565164412838526283892048404678602045626445614508293241427506110755435656591802035619036221554197413841701638040025597247140012517189434530153459",
            "-999526.399472714004924745558487742174280959297534606425446121987615879279562127703414975258631734220750778229983972859653308430370782154297373722895399005861599636821957013566636477314723339751660664618818230868728997962470940365780821107652219318608507461393737546471739313569778449847787262829689519758",
            "2635.11804627574968679150746698546684926537860473844219904054202664480773875967585945065405808941879393444480927960481154917330356034222219635456028440278614622109927632916893794632929631394563782912542725019630912377572312073034188155233244619198043408564989307377224208055269346704158128964908317509261",
            "1365298.80606570193439178248903075440797721030001172865681315316962534564202453943423169414748443288269134401947007088168439700475928251218359436660404679423063978012317674014693339414507122222142523195226338848798195882508017496231651900921236382771006795390348616128129934316532814543030209744364019054",
            "-998820.613196962571373332987243603087494541655288997511144379181267051699516068044611894994190494668286551419090745020640141977947029490154255331734395335596048929725132710856235917568105358581132488124086315318516657933430936244305844004808412923698398375511933676277573855820903335530154557900187401664",
            "4047.33414152831088829375438466435223275171308987396130932582985387723022229862922976703274469623350667935654058159331604678380195069378771816145293238142522165286453283051283913992888932063257542647725648351805459013522266494043094042838131736080143438861266226902012570363568976019577762021664023114866",
            "83.1126856152680405091850152707032759635065978000760419277583832932540620572340228677113831897018166702219588733109015730699449632190446808998629297952866898031107347012837750524595236018672538656817674901345397345906420467235996584856394111008677691361716531755180039690828818424181738021720937612332473",
            "-271.487004006624936032383971262627992814142127178531082780852027992007247020700756462748911221240546859750792472835768342005119164451772938283340118890996475035825441694527852209889256539059182981432179743946981544937115591798425012272630650241061657324339012716835440520492257904951141777212048569398409",
            "172.995634312193812329535193383815352801038985837137105437759642412099825918135903635674987213191903064253661504641496794322344612352796885741009454426914017751085087430094408938579737245424065728815449873693652886682643449058779221848518639106186622340605844416298128898323757003199986295961514068931809",
            "1553.02030744971213963777112883628816116222224922513300779830478460724114360528704175205365061869941782426355915003669694077799625540057885933310993881768684532880891858480897867129984665510634951301156097819370899150217888483360555670282615505920358106216430177949458381682435533431403640861028867643634",
            "175.729482877721774403936070331243420758342427008390700023593859226796553067648366117048302143782877238854817038614857585311307016571639203050561588272527571998881087153974255763006207376888422035958124493436136546263533358525466966789793340591585338882735253021052931213562265925007864963381639808647379",
            "-987.094530658251250328865816150283711625817471698266080780545993723465686910923865574545177650018084715406989468479547700982468885609732012850096019949134333598256299704715871688857759497730140969322324491488236509562871514705347225367092074540791495561129260238451255593578254459645239533939527551314165",
        ];
        let point = point_tokens
            .iter()
            .map(|token| {
                let value = token.parse::<ArbPrec>().unwrap();
                assert_eq!(value.to_string(), *token);
                F(value)
            })
            .collect_vec();
        let one = point[0].one();
        let zero = one.zero();
        let f = |value| one.from_i64(value);
        let model = crate::utils::load_generic_model("sm");
        let graph: Graph = include_str!("../../../../tests/resources/graphs/GL638.dot")
            .into_graph(&model)
            .unwrap();
        let lmb = &graph.loop_momentum_basis;
        assert_eq!(
            lmb.loop_edges.raw,
            vec![EdgeIndex(3), EdgeIndex(4), EdgeIndex(7), EdgeIndex(10)]
        );
        for (edge, coefficients) in [(2, [-1, 1, 0, 0]), (6, [-1, 1, 0, 1]), (10, [0, 0, 0, 1])] {
            assert_eq!(
                lmb.edge_signatures[EdgeIndex(edge)]
                    .internal
                    .to_momtrop_format(),
                coefficients
            );
        }
        let surface = Esurface {
            energies: vec![EdgeIndex(2), EdgeIndex(6), EdgeIndex(10)],
            external_shift: vec![(EdgeIndex(0), -1), (EdgeIndex(1), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = graph
            .underlying
            .new_edgevec_from_iter(lmb.edge_signatures.iter().map(|(edge, _)| {
                f(match edge.0 {
                    2 => 125,
                    6 | 10 => 173,
                    _ => 0,
                })
            }))
            .unwrap();
        let externals = ExternalFourMomenta::from_iter([1, -1].map(|sign| {
            FourMomentum::from_args(f(500), zero.clone(), zero.clone(), f(500 * sign))
        }));
        assert_eq!(externals.len(), lmb.ext_edges.len());
        let center = LoopMomenta::from_iter(
            (0..4).map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
        );
        // Preserve the inverse owner's full12D sum/order and normalization.
        // Host-null p/s components must not be dropped before forming this ray.
        let radius = point
            .iter()
            .map(F::square)
            .fold(zero.clone(), |sum, x| sum + x)
            .sqrt();
        let direction = point.iter().map(|x| x / &radius).collect_vec();
        let velocity = LoopMomenta::from_iter(
            direction
                .chunks_exact(3)
                .map(|v| ThreeMomentum::new(v[0].clone(), v[1].clone(), v[2].clone())),
        );
        let original = |r: &F<ArbPrec>| {
            surface.compute_self_and_r_derivative(r, &velocity, &center, &externals, &masses, lmb)
        };
        let ray = surface.routed_ray(&velocity, &center, &externals, &masses, lmb);
        let routed = |r: &F<ArbPrec>| ray.evaluate(r);
        let scale = original(&zero).0.abs().max(one.clone());
        assert_eq!(scale, f(529));
        let guess = f(300);
        let tolerance = f(64);
        let identity = RadialRootIdentity::new("canonical foreign phase-space Cut1".into());
        let old_strict = safeguarded_newton_iteration_and_derivative(
            &zero, &guess, original, &tolerance, 2048, 96, &scale,
        );
        let old_checked = RadialRootDiagnostics::default().solve(
            &identity, &zero, &guess, original, &tolerance, 2048, 96, &scale,
        );
        let ray_strict = safeguarded_newton_iteration_and_derivative(
            &zero, &guess, routed, &tolerance, 2048, 96, &scale,
        );
        let ray_checked = RadialRootDiagnostics::default().solve(
            &identity, &zero, &guess, routed, &tolerance, 2048, 96, &scale,
        );
        crate::debug_tags!(#integration, #sampling, #solver;
            stage = "canonical_foreign_cut_root_comparison",
            scale = %scale, tolerance = %tolerance,
            old_strict = ?old_strict, old_checked = ?old_checked,
            ray_strict = ?ray_strict, ray_checked = ?ray_checked,
            "same native direction and root budget; only radial routing association differs"
        );
        let Err(SafeguardedNewtonError::DidNotConverge {
            result,
            lower_bound,
            upper_bound,
        }) = old_strict
        else {
            panic!("retained point did not reproduce the original failed callback");
        };
        assert_eq!(result.solution, lower_bound);
        assert_eq!(result.solution, F("397339.292129129457670657255502631235707430728003236007665188851741059988603886168137387085943384940050545361933742066907875904266411120132268928043162217340042566193200203939814362605526391210977153331875750728687297380447984755419924107738341193852462036361671032388546985919083686130672457601933412077".parse::<ArbPrec>().unwrap()));
        assert_eq!(result.derivative_at_solution, F("1.89172234332619088523547124801376213069340504387632086738281927967636989017896303335256085095364995340008004517730539670922548200475769065127437640250637540730631026978180801048759551386539060201137497670459856172911183019727298113312180354967670498695170211446748295063372406403797004835125320002494423e-3".parse::<ArbPrec>().unwrap()));
        assert_eq!(result.error_of_function, F("-1.14679433441675535850302203255662653802548196930155482065516258074128414020703309961761627393787629857947032362638242933663599106623983500678126078018019011326084382655077394465668514937147617056084028309045057099969394177818124344214140305955891007956258616793631925999960450772417858669201530418435211e-296".parse::<ArbPrec>().unwrap()));
        assert_eq!(result.num_iterations_used, 17);
        assert!(upper_bound > lower_bound);
        assert!(result.error_of_function.abs() > one.epsilon() * &tolerance * &scale);
        // The existing wrapper is tested independently: a tiny bracket alone
        // must not silently enlarge this callback's residual allowance.
        assert!(old_checked.is_err());
        let prepared = ray_strict.expect("pre-routed same-budget comparison did not converge");
        let checked = ray_checked.unwrap();
        assert_eq!(prepared.solution, checked.solution);
        assert!(prepared.derivative_at_solution > zero);
        assert!(prepared.error_of_function.abs() <= one.epsilon() * &tolerance * &scale);
        let residual = surface
            .evaluate_routed_enclosed(&prepared.solution, &velocity, &externals, &masses, lmb)
            .unwrap();
        let allowed = (one.epsilon() * &tolerance * &scale)
            .0
            .mpfr_enclosure(2048)
            .0;
        assert!(residual[0] >= -allowed.clone() && residual[1] <= allowed);
        crate::debug_tags!(#integration, #sampling, #solver;
            stage = "canonical_foreign_cut_original_residual",
            native_root = %prepared.solution,
            original_residual_lower = %residual[0], original_residual_upper = %residual[1],
            "directed original routing at the pre-routed callback root"
        );
    }

    #[test]
    fn implicit_sampling_map_uses_graph_esurface_root_and_roundtrips() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph implicit_sampling_esurface {
            ext [style=invis]
            node [num=1]
            edge [num=1 mass=0]
            ext -> a:0 [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2 lmb_id=1]
            a -> b [id=3]
            b:1 -> ext [id=4]
        })
        .unwrap();
        let lmb = graph.loop_momentum_basis.clone();
        let surface = Esurface {
            energies: vec![EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = graph
            .underlying
            .new_edgevec_from_iter([F(0.0), F(173.0), F(173.0), F(125.0), F(0.0)])
            .unwrap();
        let externals = ExternalFourMomenta::from_iter(
            [FourMomentum::from_args(F(1000.0), F(0.0), F(0.0), F(0.0)); 2],
        );
        let map = surface
            .sampling_radial_map(&lmb, &masses, &externals, 400.0, 2.0)
            .unwrap();
        // Exercise both sides of the cut shell with a sizeable branch
        // split; the outer compactification carries its (1-split) factor.
        for radial_coordinate in [0.1, 0.9] {
            let coordinates = [radial_coordinate, 0.27, 0.61, 0.39, 0.72, 0.58];
            let forward = map.forward(&coordinates).unwrap();
            let step = 1.0e-5;
            let mut derivative_matrix = vec![vec![0.0; 6]; 6];
            for axis in 0..6 {
                let mut plus = coordinates;
                let mut minus = coordinates;
                plus[axis] += step;
                minus[axis] -= step;
                let plus = map.forward(&plus).unwrap();
                let minus = map.forward(&minus).unwrap();
                for (component, row) in derivative_matrix.iter_mut().enumerate().take(6) {
                    row[axis] = (plus.point[component] - minus.point[component]) / (2.0 * step);
                }
            }
            let numerical_jacobian = SamplingMapAffine::new(derivative_matrix, vec![0.0; 6])
                .unwrap()
                .determinant();
            assert!((numerical_jacobian / forward.jacobian - 1.0).abs() < 1.0e-4);

            assert!(
                forward
                    .diagnostics
                    .iter()
                    .any(|diagnostic| diagnostic == "implicit_surface:regular_root")
            );
            assert!(forward.jacobian.is_finite() && forward.jacobian > 0.0);
            let inverse = map.inverse(&forward.point).unwrap();
            assert!(inverse.residual < 1.0e-9, "{}", inverse.residual);
            for (actual, expected) in inverse.coordinates.iter().zip(coordinates) {
                assert!((actual - expected).abs() < 1.0e-9);
            }
        }

        // The graph factory must bind source data in the evaluation precision.
        // These two energies are indistinguishable in f64 but define different
        // physical shells and therefore different forward maps at Quad precision.
        let one = F::<crate::utils::f128>::from_f64(1.0);
        let energy = one.from_i64(1000);
        let displacement = one.from_i64(10).powi(-20);
        assert_eq!(energy.into_f64(), (energy + displacement).into_f64());
        let native_masses = graph
            .underlying
            .new_edgevec_from_iter(
                masses
                    .iter()
                    .map(|(_, mass)| F::<crate::utils::f128>::from_ff64(*mass)),
            )
            .unwrap();
        let coordinates = [0.19, 0.27, 0.61, 0.39, 0.72, 0.58]
            .map(|value| F::<crate::utils::f128>::from_f64(value).0);
        let mut points = Vec::new();
        for energy in [energy, energy + displacement] {
            let externals = ExternalFourMomenta::from_iter(
                [FourMomentum::from_args(energy, one.zero(), one.zero(), one.zero()); 2],
            );
            let map = surface
                .sampling_radial_map(&lmb, &native_masses, &externals, 400.0, 2.0)
                .unwrap();
            let forward = map.forward(&coordinates).unwrap();
            let inverse = map.inverse(&forward.point).unwrap();
            for (actual, expected) in inverse.coordinates.iter().zip(coordinates) {
                assert!((F(*actual) - F(expected)).abs() < one.from_i64(10).powi(-25));
            }
            points.push(forward.point);
        }
        assert_ne!(points[0], points[1]);
        assert_eq!(
            points[0]
                .iter()
                .map(|value| F(*value).into_f64())
                .collect::<Vec<_>>(),
            points[1]
                .iter()
                .map(|value| F(*value).into_f64())
                .collect::<Vec<_>>(),
        );
    }

    #[test]
    fn classification_preserves_the_previous_existing_predicate() {
        let e_cm = F(10.0);
        let shift_tolerance = F(ESURFACE_SHIFT_THRESHOLD) * e_cm;
        let normalized_margin_tolerance = F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD);
        let invariant_tolerance = normalized_margin_tolerance * e_cm * e_cm;

        for shift_factor in [-2, -1, 0, 1, 2] {
            let shift_part = shift_tolerance * shift_tolerance.from_i64(shift_factor);
            for margin_factor in [-2, -1, 0, 1, 2] {
                let invariant_margin =
                    invariant_tolerance * invariant_tolerance.from_i64(margin_factor);
                let was_existing =
                    shift_part < -&shift_tolerance && invariant_margin > invariant_tolerance;
                let classification = Esurface::classify_invariant_margin(
                    &shift_part,
                    invariant_margin,
                    &e_cm,
                    &normalized_margin_tolerance,
                );

                assert_eq!(
                    classification.is_existing(),
                    was_existing,
                    "existence changed for shift factor {shift_factor} and margin factor {margin_factor}",
                );
            }
        }
    }

    #[test]
    fn invalid_programmatic_tolerances_cannot_invert_classification() {
        let e_cm = F(10.0);
        let shift_part = F(-1.0);
        let invariant_margin = F(0.0);

        for tolerance in [
            F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            F(-DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            F(f64::NAN),
            F(f64::INFINITY),
        ] {
            assert!(matches!(
                Esurface::classify_invariant_margin(
                    &shift_part,
                    invariant_margin,
                    &e_cm,
                    &tolerance,
                ),
                EsurfaceExistence::Pinched { .. }
            ));
        }
    }

    #[test]
    fn positive_external_energy_filter_uses_energy_conservation() {
        let incoming_edges = (0..4).map(EdgeIndex::from).collect_vec();
        let outgoing_edges = (4..9).map(EdgeIndex::from).collect_vec();
        let shift_is_negative = |external_shift: &[(usize, i64)]| {
            Esurface {
                energies: vec![],
                external_shift: external_shift
                    .iter()
                    .map(|(edge, coefficient)| (EdgeIndex::from(*edge), *coefficient))
                    .collect(),
                vertex_set: VertexSet::dummy(),
            }
            .external_shift_is_strictly_negative_for_positive_energies(
                &incoming_edges,
                &outgoing_edges,
            )
        };

        assert!(
            shift_is_negative(&[(0, -1), (1, -1)]),
            "a negative proper subset of incoming energies must be retained"
        );
        assert!(
            shift_is_negative(&[(4, -1)]),
            "a negative proper subset of outgoing energies must be retained"
        );
        assert!(
            shift_is_negative(&[(0, -1), (1, -1), (2, -1), (3, -1), (4, 1)]),
            "energy conservation turns minus all incoming plus one outgoing into minus the remaining outgoing energies"
        );
        assert!(
            !shift_is_negative(&[
                (0, -1),
                (1, -1),
                (2, -1),
                (3, -1),
                (4, 1),
                (5, 1),
                (6, 1),
                (7, 1),
                (8, 1),
            ]),
            "the energy-conservation identity is zero rather than strictly negative"
        );
        assert!(
            !shift_is_negative(&[(0, -1), (4, 1)]),
            "a sign-indefinite difference must not be accepted from positivity alone"
        );
        assert!(
            !shift_is_negative(&[(9, -1)]),
            "a shift outside the external-energy partition must be rejected conservatively"
        );
    }

    #[test]
    fn massless_two_to_two_surface_is_existing_away_from_collinear_pinch() {
        let dummy_graph = dummy_hedge_graph(4);
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![EdgeIndex::from(2)].into(),
            ext_edges: vec![].into(),
            edge_signatures: dummy_graph
                .new_edgevec_from_iter(vec![
                    LoopExtSignature::from((vec![0], vec![1, 0])),
                    LoopExtSignature::from((vec![0], vec![0, 1])),
                    LoopExtSignature::from((vec![1], vec![0, 0])),
                    LoopExtSignature::from((vec![-1], vec![-1, -1])),
                ])
                .unwrap(),
        };
        let esurface = Esurface {
            energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift: vec![(EdgeIndex::from(0), -1), (EdgeIndex::from(1), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = dummy_graph.new_edgevec_from_iter(vec![F(0.0); 4]).unwrap();
        // Either side of the 2->2 sandwich determines the same total four-momentum by
        // momentum conservation, so the incoming pair is sufficient for this classification.
        let classify_pair = |second_spatial_momentum: (f64, f64), threshold: f64| {
            let external_momenta = ExternalFourMomenta::from_iter([
                FourMomentum::from_args(F(5.0), F(5.0), F(0.0), F(0.0)),
                FourMomentum::from_args(
                    F(5.0),
                    F(second_spatial_momentum.0),
                    F(second_spatial_momentum.1),
                    F(0.0),
                ),
            ]);
            esurface.classify_existence(&external_momenta, &lmb, &masses, &F(10.0), &F(threshold))
        };

        assert!(matches!(
            classify_pair((0.0, 5.0), DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            EsurfaceExistence::Existing { .. }
        ));
        assert!(matches!(
            classify_pair((5.0, 0.0), DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            EsurfaceExistence::Pinched { .. }
        ));
        assert!(matches!(
            classify_pair((0.0, 5.0), 1.0),
            EsurfaceExistence::Pinched { .. }
        ));
    }

    #[test]
    fn classifies_existing_pinched_and_non_existing_surfaces() {
        let dummy_graph = dummy_hedge_graph(5);
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![EdgeIndex::from(2), EdgeIndex::from(3)].into(),
            ext_edges: vec![].into(),
            edge_signatures: dummy_graph
                .new_edgevec_from_iter(vec![
                    LoopExtSignature::from((vec![0, 0], vec![1])),
                    LoopExtSignature::from((vec![0, 0], vec![-1])),
                    LoopExtSignature::from((vec![1, 0], vec![0])),
                    LoopExtSignature::from((vec![0, 1], vec![0])),
                    LoopExtSignature::from((vec![1, 1], vec![-1])),
                ])
                .unwrap(),
        };
        let esurface = Esurface {
            energies: vec![EdgeIndex::from(2), EdgeIndex::from(3), EdgeIndex::from(4)],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        let masses = dummy_graph.new_edgevec_from_iter(vec![F(0.0); 5]).unwrap();

        let classification = |energy| {
            let external_momenta = ExternalFourMomenta::from_iter([FourMomentum::from_args(
                F(energy),
                F(-10.0),
                F(0.0),
                F(0.0),
            )]);
            esurface.classify_existence(
                &external_momenta,
                &lmb,
                &masses,
                &F(10.0),
                &F(DEFAULT_ESURFACE_EXISTENCE_THRESHOLD),
            )
        };

        assert!(matches!(
            classification(11.0),
            EsurfaceExistence::Existing { .. }
        ));
        assert!(matches!(
            classification(10.0),
            EsurfaceExistence::Pinched { .. }
        ));
        assert!(matches!(
            classification(10.0 + 1.0e-8),
            EsurfaceExistence::Pinched { .. }
        ));
        assert!(matches!(
            classification(9.0),
            EsurfaceExistence::NonExisting { .. }
        ));
    }

    #[test]
    fn test_esurface() {
        let dummy_graph = dummy_hedge_graph(5);

        let _energies_cache = dummy_graph
            .new_edgevec_from_iter([F(1.), F(2.), F(3.), F(4.), F(5.)])
            .unwrap();

        let energies = vec![EdgeIndex::from(0), EdgeIndex::from(1), EdgeIndex::from(2)];

        let external_shift = vec![(EdgeIndex::from(3), 1), (EdgeIndex::from(4), 1)];

        let mut esurface = Esurface {
            energies,
            external_shift,
            vertex_set: VertexSet::dummy(),
            //subspace_graph: dummy_graph.full_graph(),
        };

        let shift_rewrite = ShiftRewrite {
            dependent_momentum: EdgeIndex::from(4),
            dependent_momentum_expr: vec![
                (EdgeIndex::from(1), -1),
                (EdgeIndex::from(2), -1),
                (EdgeIndex::from(3), -1),
            ],
        };

        esurface.canonicalize_shift(&shift_rewrite);

        assert_eq!(
            esurface.external_shift,
            vec![(EdgeIndex::from(1), -1), (EdgeIndex::from(2), -1)]
        );

        let energies = vec![EdgeIndex::from(0), EdgeIndex::from(2)];

        let external_shift = vec![(EdgeIndex::from(1), -1)];

        let _esurface = Esurface {
            energies,
            external_shift,
            vertex_set: VertexSet::dummy(),
            //subspace_graph: dummy_graph.full_graph(),
        };
    }

    #[test]
    fn test_add_external_shifts() {
        let shift_1 = vec![
            (EdgeIndex::from(0), 1),
            (EdgeIndex::from(1), 1),
            (EdgeIndex::from(2), -1),
        ];
        let shift_2 = vec![(EdgeIndex::from(1), -1), (EdgeIndex::from(2), 1)];

        let add = add_external_shifts(&shift_1, &shift_2);

        assert_eq!(add, vec![(EdgeIndex::from(0), 1)]);

        let shift_3 = vec![(EdgeIndex::from(3), 1), (EdgeIndex::from(4), -1)];
        let shift_4 = vec![
            (EdgeIndex::from(0), 1),
            (EdgeIndex::from(1), 1),
            (EdgeIndex::from(2), 1),
            (EdgeIndex::from(4), 1),
        ];

        let add = add_external_shifts(&shift_3, &shift_4);

        assert_eq!(
            add,
            vec![
                (EdgeIndex::from(0), 1),
                (EdgeIndex::from(1), 1),
                (EdgeIndex::from(2), 1),
                (EdgeIndex::from(3), 1)
            ]
        );
    }

    #[test]
    fn test_esurface_equality() {
        let esurface_1 = Esurface {
            energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
            vertex_set: VertexSet::dummy(),
            //subspace_graph: unsafe { InternalSubGraph::new_unchecked(SuBitGraph::new()) },
        };

        let esurface_2 = Esurface {
            energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), 1), (EdgeIndex::from(1), 1)],
            vertex_set: VertexSet::dummy(),
            //subspace_graph: unsafe { InternalSubGraph::new_unchecked(SuBitGraph::new()) },
        };

        assert_eq!(esurface_1, esurface_2);
    }

    mod failing {
        use super::*;

        #[test]
        fn test_to_atom() {
            let external_shift = vec![(EdgeIndex::from(1), -1)];

            let esurface = Esurface {
                energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
                external_shift,
                vertex_set: VertexSet::dummy(),
                // subspace_graph: unsafe { InternalSubGraph::new_unchecked(SuBitGraph::new()) },
            };

            let esurface_atom = esurface.to_atom(&[]);
            let expected_atom = parse!("Q(2, cind(0)) + Q(3, cind(0)) - P(1, cind(0))");

            let diff = esurface_atom - &expected_atom;
            let diff = diff.expand();
            assert_eq!(diff, Atom::new());
        }

        #[test]
        fn test_from_cut_left_dt() {
            let mut hedge_graph_builder = HedgeGraphBuilder::new();
            let nodes = (0..4)
                .map(|_| hedge_graph_builder.add_node(()))
                .collect_vec();

            hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);

            hedge_graph_builder.add_external_edge(
                nodes[0],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );
            hedge_graph_builder.add_external_edge(
                nodes[3],
                (),
                Orientation::Undirected,
                Flow::Source,
            );

            let double_triangle: HedgeGraph<(), (), ()> =
                hedge_graph_builder.build::<NodeStorageVec<()>>();
            let node_0 = double_triangle.iter_crown(nodes[0]).into();
            let node_3 = double_triangle.iter_crown(nodes[3]).into();

            let cuts = double_triangle.all_cuts(node_0, node_3);

            let cross_section_cuts = cuts
                .into_iter()
                .map(|(node_l, cut, node_r)| CrossSectionCut {
                    cut,
                    left: node_l,
                    right: node_r,
                })
                .map(|cut| Esurface::new_from_cut_left(&double_triangle, &cut, None))
                .collect_vec();

            let expected_esurfaces = vec![
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(2), EdgeIndex::from(4)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(3), EdgeIndex::from(4)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(1), EdgeIndex::from(2), EdgeIndex::from(3)],
                    external_shift: vec![(EdgeIndex::from(5), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: double_triangle.full_graph(),
                },
            ];

            for expected_esurface in expected_esurfaces {
                assert!(cross_section_cuts.contains(&expected_esurface));
            }
        }

        #[test]
        fn test_from_cut_left_box() {
            let mut hedge_graph_builder = HedgeGraphBuilder::new();
            let nodes = (0..4)
                .map(|_| hedge_graph_builder.add_node(()))
                .collect_vec();

            hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);
            hedge_graph_builder.add_edge(nodes[3], nodes[0], (), Orientation::Undirected);

            hedge_graph_builder.add_external_edge(
                nodes[0],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );
            hedge_graph_builder.add_external_edge(
                nodes[1],
                (),
                Orientation::Undirected,
                Flow::Source,
            );
            hedge_graph_builder.add_external_edge(
                nodes[2],
                (),
                Orientation::Undirected,
                Flow::Source,
            );
            hedge_graph_builder.add_external_edge(
                nodes[3],
                (),
                Orientation::Undirected,
                Flow::Sink,
            );

            let box_graph: HedgeGraph<(), (), ()> =
                hedge_graph_builder.build::<NodeStorageVec<()>>();

            let node_0 = box_graph.iter_crown(nodes[0]).into();
            let node_2 = box_graph.iter_crown(nodes[2]).into();

            let cuts = box_graph.all_cuts(node_0, node_2);
            assert_eq!(cuts.len(), 4);

            let cross_section_cuts = cuts
                .into_iter()
                .map(|(node_l, cut, node_r)| CrossSectionCut {
                    cut,
                    left: node_l,
                    right: node_r,
                })
                .map(|cut| Esurface::new_from_cut_left(&box_graph, &cut, None))
                .collect_vec();

            let expected_esurfaces = vec![
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(3)],
                    external_shift: vec![(EdgeIndex::from(4), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(0), EdgeIndex::from(2)],
                    external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(7), -1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(1), EdgeIndex::from(3)],
                    external_shift: vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
                Esurface {
                    energies: vec![EdgeIndex::from(1), EdgeIndex::from(2)],
                    external_shift: vec![
                        (EdgeIndex::from(4), -1),
                        (EdgeIndex::from(5), 1),
                        (EdgeIndex::from(7), -1),
                    ],
                    vertex_set: VertexSet::dummy(),
                    //subspace_graph: box_graph.full_graph(),
                },
            ];

            for expected_esurface in expected_esurfaces {
                assert!(cross_section_cuts.contains(&expected_esurface));
            }
        }
    }
}
