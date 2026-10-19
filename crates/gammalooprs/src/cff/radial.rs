//! Retained affine energy equations shared by LU roots, jets and host certificates.

use std::{rc::Rc, time::Instant};

use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use rug::{
    Float,
    float::{Round, Special},
};
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
        ArbPrec, F, FloatLike, SamplingPrecision, compute_loop_part, compute_t_part_of_shift_part,
        hyperdual_utils::new_constant,
        newton_solver::{
            NewtonIterationResult, RadialRootDiagnostics, RadialRootIdentity,
            SafeguardedNewtonBracket, SafeguardedNewtonError,
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
type EnclosedRayEnergy<B> = (i64, [[B; 2]; 3], [[B; 2]; 3], [B; 2]);

/// One represented affine energy equation, retaining its ordered edge occurrences.
/// This is native geometry, not a cross-precision or selected-host transport record.
#[derive(Debug, Clone)]
pub(crate) struct EnergySurfaceRay<T: FloatLike> {
    // A Vec deliberately preserves repeated energy occurrences; an edge-keyed
    // map would not. Algebraic reduction is only a derived root-isolation view.
    pub(super) energies: Vec<RadialEnergy<T>>,
    pub(super) shift: F<T>,
}

/// The stage that decided one signed radial isolation call.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
enum RadialIsolationPath {
    /// An error preceded any decision: native checks, collapse or convex origin enclosure.
    #[default]
    Undecided,
    ConvexInterior,
    ConvexReceding,
    /// Every energy cancels in the collapsed isolation view.
    Constant,
    Isolated,
}

/// How directed origin enclosures at one precision decide the convex
/// positive-E shortcut.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConvexOrigin {
    /// The origin is interior: the unique outward crossing is refined.
    Interior,
    /// The origin is exterior with nonnegative slope: no positive root.
    Receding,
    /// These bounds exclude both shortcuts at every precision.
    NotApplicable,
    /// A shortcut threshold is unresolved at this precision, or its origin
    /// enclosure failed.
    Unresolved,
}

/// How the binary64 first pass of one isolation call ended.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
enum Binary64Outcome {
    /// The call was decided before isolation.
    #[default]
    NotAttempted,
    /// The floating-point environment lacks gradual underflow, so no
    /// binary64 tier ran.
    Unavailable,
    /// Isolation was unresolved at binary64 resolution. Only the summary
    /// event reads the detail.
    Uncertain(&'static str),
    /// A bracket did not materialize as increasing native radii.
    Unmaterialized,
    /// Neither binary64 nor reference bounds certify a bracket's signs, or
    /// the reference enclosure failed.
    Uncertified,
    /// A bracket fails the safeguarded solver's start conditions.
    Unrefinable,
    /// Every bracket was certified before its root was refined.
    Certified,
}

/// Work of one signed radial isolation call, reported only by debug events.
#[derive(Debug, Default)]
struct RadialIsolationTrace {
    path: RadialIsolationPath,
    /// Collapsed signed terms of the isolation view.
    terms: usize,
    /// The binary64 tier of the convex positive-E shortcut, when it ran;
    /// `binary64` records an unavailable environment.
    binary64_convex: Option<ConvexOrigin>,
    binary64_convex_evaluations: usize,
    /// Reference origin enclosures of the convex shortcut, which run unless
    /// the binary64 tier resolves it.
    convex_evaluations: usize,
    binary64: Binary64Outcome,
    binary64_work: RadialIsolationWork,
    /// Work of the reference pass, which runs unless binary64 certifies.
    work: RadialIsolationWork,
    /// Materialized brackets sign-certified by binary64 bounds, in either
    /// pass: the reference pass tries binary64 bounds first as well.
    binary64_certificates: usize,
    /// Materialized brackets sign-certified at the reference precision, in
    /// either pass, where binary64 bounds do not certify them.
    certificates: usize,
}

/// Work of one compactified isolation pass.
#[derive(Debug, Default)]
struct RadialIsolationWork {
    /// Subdivision cells, finite-upper probes and split retries (the budget counter).
    visited: usize,
    /// Directed compactified endpoint enclosures.
    endpoints: usize,
}

/// One bound of a directed affine-energy enclosure `[lower, upper]`. On
/// enclosures of finite natives, every operation encloses its exact real
/// result outward for all reals in its operands and never returns NaN: lower
/// bounds lie in [-inf, MAX] and upper bounds in [-MAX, +inf]. A finite point
/// enclosure, in particular [0, 0], is exact; isolation reads exact
/// cancellations and zeros from points. Nonfinite natives are outside this
/// contract: they propagate NaN or infinities, which the `is_finite_bound`
/// checks must reject before any sign is read.
pub(crate) trait DirectedBound: Clone + PartialOrd + PartialOrd<Self::Zero> {
    /// The exact zero in the form that bounds compare against most cheaply.
    type Zero;
    const ZERO: Self::Zero;
    /// Visited cells, finite-upper probes and split retries allowed in one
    /// compactified isolation pass.
    const CELL_BUDGET: usize;
    /// Whether the floating-point environment supports these rounding rules.
    /// Only the binary64 tiers depend on it, so callers choose them by
    /// `f64::available()`; the reference is always available.
    fn available() -> bool;
    fn zero() -> Self;
    fn one() -> Self;
    /// An outward enclosure of a finite native. A magnitude beyond this
    /// bound's range widens to an infinite outward bound. Only an exact zero
    /// factor absorbs it; otherwise it propagates to the finiteness checks.
    fn enclose_native<T: FloatLike>(value: &F<T>) -> [Self; 2];
    fn enclose_multiplicity(value: i64) -> [Self; 2];
    fn negated(&self) -> Self;
    fn abs_bound(&self) -> Self;
    // Given zeros of both signs, min/max may return either (MPFR: min -0,
    // max +0). Nothing observes it: comparisons treat +-0 as equal, and
    // reciprocals, the only operation whose nonzero result depends on the
    // sign of a zero, are taken only of a strictly positive divisor in
    // `enclosed_divide`.
    fn min_bound(self, other: &Self) -> Self;
    fn max_bound(self, other: &Self) -> Self;
    fn add_down(&self, other: &Self) -> Self;
    fn add_up(&self, other: &Self) -> Self;
    fn mul_down(&self, other: &Self) -> Self;
    fn mul_up(&self, other: &Self) -> Self;
    fn recip_down(&self) -> Self;
    fn recip_up(&self) -> Self;
    fn sqrt_down(&self) -> Self;
    fn sqrt_up(&self) -> Self;
    /// A representable point of `[self, other]` whenever `self + other` is
    /// finite, as for the compactified u in [0, 1].
    fn bisect(&self, other: &Self) -> Self;
    fn is_finite_bound(&self) -> bool;

    // These directed operations serve the represented and original-routed
    // affine-energy checks; they do not introduce a second certificate engine.
    // They are written once over these bounds, whose reference instance is
    // MPFR at ENCLOSURE_PRECISION. Signed radial isolation, its convex origin
    // shortcut and its bracket sign checks try binary64 bounds first; the
    // reference decides whatever they leave unresolved, and all other checks
    // use the reference only.
    fn enclosed_add(a: &[Self; 2], b: &[Self; 2]) -> [Self; 2] {
        [a[0].add_down(&b[0]), a[1].add_up(&b[1])]
    }
    /// Each directed product bounds its exact corner outward, and the operand
    /// signs select the extremal exact corner (one of two when both straddle
    /// zero), whose directed product therefore bounds the whole product range.
    /// Under monotone rounding (MPFR) it equals the all-corner extremum.
    fn enclosed_mul(a: &[Self; 2], b: &[Self; 2]) -> [Self; 2] {
        let straddles = |x: &[Self; 2]| x[0] < Self::ZERO && x[1] > Self::ZERO;
        // Products commute exactly, so `b` straddles zero only if both do.
        let (a, b) = if straddles(b) { (b, a) } else { (a, b) };
        if straddles(b) {
            return [
                a[0].mul_down(&b[1]).min_bound(&a[1].mul_down(&b[0])),
                a[0].mul_up(&b[0]).max_bound(&a[1].mul_up(&b[1])),
            ];
        }
        // A one-signed `b` makes the product monotone in `a`; the sign of the
        // selected bound of `a` then selects the bound of `b`.
        let (lower, upper) = if b[0] >= Self::ZERO {
            (&a[0], &a[1])
        } else {
            (&a[1], &a[0])
        };
        [
            lower.mul_down(if *lower >= Self::ZERO { &b[0] } else { &b[1] }),
            upper.mul_up(if *upper >= Self::ZERO { &b[1] } else { &b[0] }),
        ]
    }
    /// The upper bound squares the bound of larger magnitude; the lower bound
    /// squares the other one unless `a` straddles zero.
    fn enclosed_square(a: &[Self; 2]) -> [Self; 2] {
        if a[0] < Self::ZERO && a[1] > Self::ZERO {
            [
                Self::zero(),
                a[0].mul_up(&a[0]).max_bound(&a[1].mul_up(&a[1])),
            ]
        } else if a[1] <= Self::ZERO {
            [a[1].mul_down(&a[1]), a[0].mul_up(&a[0])]
        } else {
            [a[0].mul_down(&a[0]), a[1].mul_up(&a[1])]
        }
    }
    fn enclosed_divide(a: &[Self; 2], b: &[Self; 2]) -> [Self; 2] {
        debug_assert!(b[0] > Self::ZERO);
        Self::enclosed_mul(a, &[b[1].recip_down(), b[0].recip_up()])
    }
    fn enclosed_neg(a: &[Self; 2]) -> [Self; 2] {
        [a[1].negated(), a[0].negated()]
    }
    /// The enclosed value is a sum of squares, so a negative lower bound is
    /// only rounding and is clamped to zero.
    fn enclosed_sqrt(a: &[Self; 2]) -> [Self; 2] {
        let lower = if a[0] < Self::ZERO {
            Self::zero()
        } else {
            a[0].sqrt_down()
        };
        [lower, a[1].sqrt_up()]
    }
    /// Multiplicities +-1 scale without rounding.
    fn enclosed_scale(a: [Self; 2], coefficient: i64) -> [Self; 2] {
        match coefficient {
            1 => a,
            -1 => Self::enclosed_neg(&a),
            _ => Self::enclosed_mul(&a, &Self::enclose_multiplicity(coefficient)),
        }
    }
    fn enclosed_sum(terms: impl Iterator<Item = (i64, [Self; 2])>) -> [Self; 2] {
        terms.fold([Self::zero(), Self::zero()], |sum, (coefficient, value)| {
            Self::enclosed_add(&sum, &Self::enclosed_scale(value, coefficient))
        })
    }
}

/// The precision of the reference MPFR bounds.
pub(crate) const ENCLOSURE_PRECISION: u32 = 2048;

/// MPFR bounds at [`ENCLOSURE_PRECISION`], the reference certificate.
impl DirectedBound for Float {
    // A sign test; comparing with a binary64 zero would convert it first.
    type Zero = Special;
    const ZERO: Special = Special::Zero;
    const CELL_BUDGET: usize = 65536;

    fn available() -> bool {
        true
    }
    fn zero() -> Self {
        Self::with_val(ENCLOSURE_PRECISION, 0)
    }
    fn one() -> Self {
        Self::with_val(ENCLOSURE_PRECISION, 1)
    }
    fn enclose_native<T: FloatLike>(value: &F<T>) -> [Self; 2] {
        let (lo, hi) = value.0.mpfr_enclosure(ENCLOSURE_PRECISION);
        [lo, hi]
    }
    fn enclose_multiplicity(value: i64) -> [Self; 2] {
        let value = Self::with_val(ENCLOSURE_PRECISION, value);
        [value.clone(), value]
    }
    fn negated(&self) -> Self {
        -self.clone()
    }
    fn abs_bound(&self) -> Self {
        self.clone().abs()
    }
    fn min_bound(self, other: &Self) -> Self {
        self.min(other)
    }
    fn max_bound(self, other: &Self) -> Self {
        self.max(other)
    }
    fn add_down(&self, other: &Self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self + other, Round::Down).0
    }
    fn add_up(&self, other: &Self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self + other, Round::Up).0
    }
    fn mul_down(&self, other: &Self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self * other, Round::Down).0
    }
    fn mul_up(&self, other: &Self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self * other, Round::Up).0
    }
    fn recip_down(&self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self.recip_ref(), Round::Down).0
    }
    fn recip_up(&self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self.recip_ref(), Round::Up).0
    }
    fn sqrt_down(&self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self.sqrt_ref(), Round::Down).0
    }
    fn sqrt_up(&self) -> Self {
        Self::with_val_round(ENCLOSURE_PRECISION, self.sqrt_ref(), Round::Up).0
    }
    fn bisect(&self, other: &Self) -> Self {
        Self::with_val(ENCLOSURE_PRECISION, self + other) / 2
    }
    fn is_finite_bound(&self) -> bool {
        self.is_finite()
    }
}

/// Binary64 bounds, a first attempt at the same certificate. Each result is
/// the nearest one stepped one ulp outward. Nearest rounding errs by at most
/// half an ulp, also on the subnormal grid, so the outward neighbour bounds
/// the exact result in any rounding mode; a lower bound that overflows steps
/// back to MAX, which the exact result exceeds. A result stays unstepped only
/// where an operand rule makes it exact.
impl DirectedBound for f64 {
    type Zero = f64;
    const ZERO: f64 = 0.0;
    // Binary64 resolution is exhausted long before the reference budget.
    const CELL_BUDGET: usize = 256;

    /// The exact-operand rules assume gradual underflow. Flush-to-zero or
    /// denormals-are-zero, set by foreign code, would make a rounded
    /// subnormal look exact; comparing bits keeps DAZ from hiding the probe.
    fn available() -> bool {
        // A volatile load cannot be folded, so the sum is computed at run
        // time in the current floating-point environment.
        static TINY: u64 = 1;
        // SAFETY: a reference to a static is valid, aligned and initialized.
        let tiny = f64::from_bits(unsafe { std::ptr::read_volatile(&TINY) });
        (tiny + tiny).to_bits() == 2
    }
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn enclose_native<T: FloatLike>(value: &F<T>) -> [Self; 2] {
        // A nonfinite native has no enclosure. NaN reaches the finiteness
        // checks even beside an exact zero factor, which absorbs infinities.
        if !value.0.is_finite() {
            return [f64::NAN; 2];
        }
        // Every lane converts to the nearest binary64 value.
        let nearest = value.0.into_f64();
        if F(T::from_f64_exact_binary(nearest)) == *value {
            [nearest; 2]
        } else {
            Outward::around(nearest)
        }
    }
    fn enclose_multiplicity(value: i64) -> [Self; 2] {
        let nearest = value as f64;
        if nearest as i128 == i128::from(value) {
            [nearest; 2]
        } else {
            Outward::around(nearest)
        }
    }
    fn negated(&self) -> Self {
        -self
    }
    fn abs_bound(&self) -> Self {
        self.abs()
    }
    fn min_bound(self, other: &Self) -> Self {
        self.min(*other)
    }
    fn max_bound(self, other: &Self) -> Self {
        self.max(*other)
    }
    fn add_down(&self, other: &Self) -> Self {
        Outward::Down.sum(*self, *other)
    }
    fn add_up(&self, other: &Self) -> Self {
        Outward::Up.sum(*self, *other)
    }
    fn mul_down(&self, other: &Self) -> Self {
        Outward::Down.product(*self, *other)
    }
    fn mul_up(&self, other: &Self) -> Self {
        Outward::Up.product(*self, *other)
    }
    fn recip_down(&self) -> Self {
        Outward::Down.step(1.0 / self)
    }
    fn recip_up(&self) -> Self {
        Outward::Up.step(1.0 / self)
    }
    fn sqrt_down(&self) -> Self {
        Outward::Down.root(*self)
    }
    fn sqrt_up(&self) -> Self {
        Outward::Up.root(*self)
    }
    fn bisect(&self, other: &Self) -> Self {
        0.5 * (self + other)
    }
    fn is_finite_bound(&self) -> bool {
        self.is_finite()
    }
}

/// The outward direction of one binary64 bound.
#[derive(Clone, Copy)]
enum Outward {
    Down,
    Up,
}

impl Outward {
    fn step(self, nearest: f64) -> f64 {
        match self {
            Self::Down => nearest.next_down(),
            Self::Up => nearest.next_up(),
        }
    }
    /// Both outward neighbours of an inexact nearest value.
    fn around(nearest: f64) -> [f64; 2] {
        [Self::Down.step(nearest), Self::Up.step(nearest)]
    }
    /// x+0 and exact cancellation are exact: with gradual underflow a
    /// nonzero exact sum of binary64 values never rounds to zero.
    fn sum(self, a: f64, b: f64) -> f64 {
        let sum = a + b;
        if a == 0.0 || b == 0.0 || sum == 0.0 {
            sum
        } else {
            self.step(sum)
        }
    }
    /// Bounds enclose finite reals, so a zero factor gives an exact zero even
    /// beside an unbounded bound; factors +-1 are exact as well.
    fn product(self, a: f64, b: f64) -> f64 {
        if a == 0.0 || b == 0.0 {
            if a.is_nan() || b.is_nan() {
                f64::NAN
            } else {
                0.0
            }
        } else if a.abs() == 1.0 || b.abs() == 1.0 {
            a * b
        } else {
            self.step(a * b)
        }
    }
    /// The square root of zero is exact.
    fn root(self, x: f64) -> f64 {
        if x == 0.0 { 0.0 } else { self.step(x.sqrt()) }
    }
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

    pub(super) fn evaluate_enclosed<B: DirectedBound>(
        energies: impl Iterator<Item = EnclosedRayEnergy<B>>,
        shift: [B; 2],
        radius: &[B; 2],
        with_derivative: bool,
    ) -> Result<([B; 2], [B; 2])> {
        let zero = B::zero();
        let point = |x: &B| [x.clone(), x.clone()];
        let uncertain = |detail: &str| SamplingEvaluationError::UncertainGeometry {
            detail: format!("radial energy enclosure: {detail}"),
        };
        let mut value = shift;
        let mut slope = point(&zero);
        for (coefficient, velocity, offset, mass) in energies {
            // Exact zero velocity makes the energy constant, even at E=0.
            // Never infer this from a norm that may have underflowed.
            let constant = velocity.iter().flatten().all(|value| *value == B::ZERO);
            // Only the derivative needs the numerator q.v.
            let derivative = with_derivative && !constant;
            let mut norm = B::enclosed_square(&mass);
            let mut dot = point(&zero);
            for (v, b) in velocity.into_iter().zip(offset) {
                let q = B::enclosed_add(&B::enclosed_mul(&v, radius), &b);
                norm = B::enclosed_add(&norm, &B::enclosed_square(&q));
                if derivative {
                    dot = B::enclosed_add(&dot, &B::enclosed_mul(&v, &q));
                }
            }
            if norm.iter().chain(&dot).any(|x| !x.is_finite_bound()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "radial energy enclosure",
                    detail: "nonfinite directed energy or numerator".into(),
                }
                .into());
            }
            let energy = B::enclosed_sqrt(&norm);
            if derivative {
                if energy[0] <= B::ZERO {
                    return Err(uncertain(
                        "energy lower bound touches zero on the fixed candidate interval",
                    )
                    .into());
                }
                slope = B::enclosed_add(
                    &slope,
                    &B::enclosed_scale(B::enclosed_divide(&dot, &energy), coefficient),
                );
            }
            value = B::enclosed_add(&value, &B::enclosed_scale(energy, coefficient));
        }
        if value.iter().chain(&slope).any(|x| !x.is_finite_bound()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "radial energy enclosure",
                detail: "nonfinite directed residual or derivative".into(),
            }
            .into());
        }
        Ok((value, slope))
    }

    pub(crate) fn evaluate_interval<B: DirectedBound>(
        &self,
        radius: &[B; 2],
        with_derivative: bool,
    ) -> Result<([B; 2], [B; 2])> {
        let native = B::enclose_native;
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
        let radius_squared = Float::enclosed_add(
            &Float::enclosed_square(&canonical_normals[0]),
            &Float::enclosed_square(&canonical_normals[1]),
        );
        let radius_lower = Float::with_val_round(
            ENCLOSURE_PRECISION,
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
            Float::enclosed_add(
                &native_normals[axis],
                &[-canonical[1].clone(), -canonical[0].clone()],
            )
        });
        let error_squared = Float::enclosed_add(
            &Float::enclosed_square(&difference[0]),
            &Float::enclosed_square(&difference[1]),
        );
        let error_upper =
            Float::with_val_round(ENCLOSURE_PRECISION, error_squared[1].sqrt_ref(), Round::Up).0;
        let budget = Float::with_val(ENCLOSURE_PRECISION, tolerance);
        let half_budget = Float::with_val_round(ENCLOSURE_PRECISION, &budget / 2, Round::Down).0;
        let allowed = Float::with_val_round(
            ENCLOSURE_PRECISION,
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
            residual_tolerance = %(e_cm * guess.epsilon() * self.refinement_tolerance(&guess)),
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
            self.oriented_equation(crossing),
            &self.refinement_tolerance(upper),
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

    /// Whether [`Self::refine_root`] starts from this bracket without
    /// expanding it: the solver's own start check, which records no
    /// diagnostics.
    fn refinable(&self, lower: &F<T>, upper: &F<T>, crossing: i64, e_cm: &F<T>) -> bool {
        SafeguardedNewtonBracket::new(
            lower,
            upper,
            self.oriented_equation(crossing),
            &self.refinement_tolerance(upper),
            0,
            e_cm,
        )
        .is_ok()
    }

    /// The equation of [`Self::refine_root`]: eta oriented so that a
    /// `crossing` bracket ascends through its root.
    fn oriented_equation(&self, crossing: i64) -> impl Fn(&F<T>) -> (F<T>, F<T>) + '_ {
        move |radius| {
            let (value, derivative) = self.evaluate(radius);
            if crossing > 0 {
                (value, derivative)
            } else {
                (-value, -derivative)
            }
        }
    }

    /// The residual tolerance of [`Self::refine_root`], in units of the
    /// solver's epsilon times the energy scale, at the precision of `reference`.
    fn refinement_tolerance(&self, reference: &F<T>) -> F<T> {
        reference.from_i64((4 * self.energies.len() + 1) as i64)
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
        let precision = ENCLOSURE_PRECISION;
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
        let native = Float::enclose_native;
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
            let ratio = Float::enclosed_divide(&slope, &native(derivative));
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
            Float::enclosed_divide(&interval, &t0),
            Float::enclosed_divide(&t0, &interval),
        ] {
            let mut determinant_ratio = point(&one);
            for _ in 0..dimension_bound {
                determinant_ratio = Float::enclosed_mul(&determinant_ratio, &ratio);
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

    /// Enumerate simple positive roots of the represented affine-energy equation.
    /// Directed value and derivative enclosures on a compactified half-line
    /// certify every empty cell and isolate every simple crossing, first with
    /// binary64 bounds and, wherever those leave anything unresolved, again at
    /// ENCLOSURE_PRECISION. Convex positive-E rays retain their native interior
    /// shortcut. Unresolved tangencies, cusps and exhausted work or precision
    /// budgets are errors, never evidence that a channel is absent.
    pub(crate) fn positive_radius_roots(
        &self,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
    ) -> eyre::Result<Vec<NewtonIterationResult<T>>> {
        crate::debug_tags!(#integration, #cut, #solver, #inspect, #dump;
            radial_root = %identity,
            file.ray = ?self,
            "radial root isolation ray"
        );
        // Time the call whenever any debug event can be recorded, so that the
        // summary event below carries its time under every filter.
        let started = (tracing::Level::DEBUG <= tracing::level_filters::LevelFilter::current())
            .then(Instant::now);
        let mut trace = RadialIsolationTrace::default();
        let result = self.enumerate_positive_radius_roots(e_cm, diagnostics, identity, &mut trace);
        crate::debug_tags!(#integration, #cut, #solver, #summary;
            radial_root = %identity,
            lane = ?T::sampling_precision(),
            energies = self.energies.len(),
            terms = trace.terms,
            path = ?trace.path,
            binary64_convex = ?trace.binary64_convex,
            binary64_convex_evaluations = trace.binary64_convex_evaluations,
            convex_evaluations = trace.convex_evaluations,
            binary64 = ?trace.binary64,
            binary64_visited = trace.binary64_work.visited,
            binary64_endpoints = trace.binary64_work.endpoints,
            binary64_certificates = trace.binary64_certificates,
            visited = trace.work.visited,
            endpoints = trace.work.endpoints,
            certificates = trace.certificates,
            roots = result.as_ref().map_or(0, Vec::len),
            failed = result.is_err(),
            elapsed_us = started.map(|start| start.elapsed().as_secs_f64() * 1e6),
            "radial root isolation"
        );
        result
    }

    /// The enumeration behind [`Self::positive_radius_roots`]; `trace` holds
    /// the work done up to whichever exit is taken.
    fn enumerate_positive_radius_roots(
        &self,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
        trace: &mut RadialIsolationTrace,
    ) -> eyre::Result<Vec<NewtonIterationResult<T>>> {
        let zero = e_cm.zero();
        let uncertain = |detail: &str| eyre::eyre!("{identity}: signed radial equation: {detail}");
        let available = f64::available();
        if !available {
            trace.binary64 = Binary64Outcome::Unavailable;
        }
        if !e_cm.0.is_finite() || e_cm <= &zero {
            return Err(uncertain("the residual scale must be finite and positive"));
        }
        if !self.shift.0.is_finite() {
            return Err(uncertain("nonfinite external energy shift"));
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
                *multiplicity = multiplicity
                    .checked_add(energy.coefficient)
                    .ok_or_else(|| uncertain("energy multiplicity overflow"))?;
            } else {
                terms.push((energy.coefficient, energy));
            }
        }
        terms.retain(|(sign, _)| *sign != 0);
        trace.terms = terms.len();

        if self.energies.iter().all(|energy| energy.coefficient == 1)
            && self.energies.iter().any(|energy| {
                energy.velocity.px != zero
                    || energy.velocity.py != zero
                    || energy.velocity.pz != zero
            })
        {
            // Convex positive-E geometry gives a unique outward crossing from
            // an interior origin. Keep its inexpensive native LU policy for
            // arbitrary multiplicity; no subdivision is needed.
            // Binary64 bounds decide the shortcut first. Only an unresolved
            // or unavailable binary64 tier leaves the origin enclosures to the
            // reference precision, whose errors propagate as before.
            trace.binary64_convex = available.then(|| {
                self.convex_origin::<f64>(&mut trace.binary64_convex_evaluations)
                    .unwrap_or(ConvexOrigin::Unresolved)
            });
            let origin = match trace.binary64_convex {
                Some(resolved) if resolved != ConvexOrigin::Unresolved => resolved,
                _ => self.convex_origin::<Float>(&mut trace.convex_evaluations)?,
            };
            match origin {
                ConvexOrigin::Interior => {
                    trace.path = RadialIsolationPath::ConvexInterior;
                    return self
                        .solve_lu_cut(e_cm, diagnostics, identity)
                        .map(|root| vec![root])
                        .map_err(|error| {
                            uncertain(&format!("convex root refinement failed: {error:?}"))
                        });
                }
                ConvexOrigin::Receding => {
                    trace.path = RadialIsolationPath::ConvexReceding;
                    return Ok(Vec::new());
                }
                // Otherwise interval isolation decides the ray, as before.
                ConvexOrigin::NotApplicable | ConvexOrigin::Unresolved => {}
            }
        }

        if terms.is_empty() {
            trace.path = RadialIsolationPath::Constant;
            return if self.shift == zero {
                Err(uncertain("the radial equation is identically zero"))
            } else {
                Ok(Vec::new())
            };
        }
        trace.path = RadialIsolationPath::Isolated;
        let mut refine = |index: usize, lower: &F<T>, upper: &F<T>, crossing: i64| {
            self.refine_root(
                lower,
                upper,
                crossing,
                0,
                e_cm,
                diagnostics,
                &RadialRootIdentity::new(format!("{identity} signed root {index}")),
            )
            .map_err(|error| uncertain(&format!("isolated root refinement failed: {error:?}")))
        };
        // The binary64 pass certifies every bracket before refining any root.
        // Its brackets and shared radicals need not match the reference's, so
        // its roots stand even where the reference pass would have failed.
        if available {
            match self
                .compactified::<f64>(&terms)
                .isolate(&mut trace.binary64_work)
                .map_err(Binary64Outcome::Uncertain)
                .and_then(|cells| self.binary64_brackets(cells, e_cm, trace))
            {
                Ok(brackets) => {
                    trace.binary64 = Binary64Outcome::Certified;
                    return brackets
                        .iter()
                        .enumerate()
                        .map(|(index, (lower, upper, crossing))| {
                            refine(index, lower, upper, *crossing)
                        })
                        .collect();
                }
                Err(outcome) => trace.binary64 = outcome,
            }
        }
        // The reference pass refines each bracket as soon as it is certified.
        let brackets = self
            .compactified::<Float>(&terms)
            .isolate(&mut trace.work)
            .map_err(uncertain)?;
        brackets
            .into_iter()
            .enumerate()
            .map(|(index, (left, right, crossing))| {
                let lower = Self::materialized_radius(&left)?;
                let upper = Self::materialized_radius(&right)?;
                if lower >= upper {
                    return Err(uncertain("the isolated root bracket is not representable"));
                }
                if !self.bracket_signs_certified(&lower, &upper, crossing, available, trace)? {
                    return Err(uncertain(
                        "the materialized root bracket is not sign-certified",
                    ));
                }
                refine(index, &lower, &upper, crossing)
            })
            .collect()
    }

    /// The convex positive-E origin shortcut at the precision of `B`. Both
    /// predicates are monotone in the bounds, so a shortcut that `B` resolves
    /// is decided identically by any finer enclosure of the same expressions.
    fn convex_origin<B: DirectedBound>(&self, evaluations: &mut usize) -> Result<ConvexOrigin> {
        let zero = self.shift.zero();
        let origin = B::enclose_native(&zero);
        *evaluations += 1;
        let value = self.evaluate_interval(&origin, false)?.0;
        if value[1] < B::ZERO {
            return Ok(ConvexOrigin::Interior);
        }
        // An exterior origin with nonnegative slope cannot intersect a
        // convex E-surface in the positive direction. Require smooth
        // energies here; origin cusps remain with interval isolation.
        let smooth = self.energies.iter().all(|energy| {
            energy.mass != zero
                || energy.offset.px != zero
                || energy.offset.py != zero
                || energy.offset.pz != zero
        });
        // Whether every finer enclosure also rejects the receding shortcut,
        // and evaluates no slope that could fail.
        let receding_excluded = if value[0] > B::ZERO && smooth {
            *evaluations += 1;
            let slope = self.evaluate_interval(&origin, true)?.1;
            if slope[0] >= B::ZERO {
                return Ok(ConvexOrigin::Receding);
            }
            slope[1] < B::ZERO
        } else {
            !smooth || value[1] <= B::ZERO
        };
        Ok(if value[0] >= B::ZERO && receding_excluded {
            ConvexOrigin::NotApplicable
        } else {
            ConvexOrigin::Unresolved
        })
    }

    /// The brackets of the binary64 pass. Every isolated cell is materialized,
    /// sign-certified and checked against the solver's start conditions
    /// before any root is refined, so falling back to the reference records
    /// no diagnostics.
    fn binary64_brackets(
        &self,
        cells: Vec<Bracket<f64>>,
        e_cm: &F<T>,
        trace: &mut RadialIsolationTrace,
    ) -> std::result::Result<Vec<Bracket<F<T>>>, Binary64Outcome> {
        // Wherever both passes reach the same u their radii are bit-identical.
        let radius = |u: f64| Self::binary64_radius(u).map_err(|_| Binary64Outcome::Unmaterialized);
        cells
            .into_iter()
            .map(|(left, right, crossing)| {
                let (lower, upper) = (radius(left)?, radius(right)?);
                if lower >= upper {
                    return Err(Binary64Outcome::Unmaterialized);
                }
                // This pass runs only where binary64 bounds are available.
                if !matches!(
                    self.bracket_signs_certified(&lower, &upper, crossing, true, trace),
                    Ok(true)
                ) {
                    return Err(Binary64Outcome::Uncertified);
                }
                if !self.refinable(&lower, &upper, crossing, e_cm) {
                    return Err(Binary64Outcome::Unrefinable);
                }
                Ok((lower, upper, crossing))
            })
            .collect()
    }

    /// The native radius r=u/(1-u) of a compactified endpoint: the midpoint of
    /// its directed reference enclosure, rounded to ArbPrec and then to the
    /// lane.
    fn materialized_radius(u: &Float) -> eyre::Result<F<T>> {
        let point = [u.clone(), u.clone()];
        let one = [Float::one(), Float::one()];
        let denominator = Float::enclosed_add(&one, &Float::enclosed_neg(&point));
        let radius = Float::enclosed_divide(&point, &denominator);
        F::from_arb(&ArbPrec::from(radius[0].bisect(&radius[1])))
    }

    /// [`Self::materialized_radius`] of a binary64 endpoint, which lifts
    /// exactly. In the Double lane, wherever 1-u is a binary64 value, one
    /// division rounds to the same radius: a quotient of binary64 values lies
    /// at least 2^-108 (relative) from every rounding midpoint, beyond the
    /// 2^-999 by which the reference midpoint, rounded to 1000 bits, can err.
    fn binary64_radius(u: f64) -> eyre::Result<F<T>> {
        let denominator = 1.0 - u;
        // Sterbenz's lemma makes the second subtraction exact below u=1/2, and
        // from there on the first, after which the second returns u. So this
        // tests whether the first subtraction is exact.
        if T::sampling_precision() == SamplingPrecision::Double
            && (0.0..1.0).contains(&u)
            && 1.0 - denominator == u
        {
            return Ok(F(T::from_f64_exact_binary(u / denominator)));
        }
        Self::materialized_radius(&Float::with_val(ENCLOSURE_PRECISION, u))
    }

    /// Whether directed enclosures at the native radii certify the opposite
    /// signs of a `crossing` bracket: binary64 bounds first if `binary64`
    /// (they are available), then the reference, whose errors propagate.
    fn bracket_signs_certified(
        &self,
        lower: &F<T>,
        upper: &F<T>,
        crossing: i64,
        binary64: bool,
        trace: &mut RadialIsolationTrace,
    ) -> Result<bool> {
        if binary64
            && self
                .bracket_signs::<f64>(lower, upper, crossing)
                .is_ok_and(|certified| certified)
        {
            trace.binary64_certificates += 1;
            return Ok(true);
        }
        let certified = self.bracket_signs::<Float>(lower, upper, crossing)?;
        trace.certificates += usize::from(certified);
        Ok(certified)
    }

    fn bracket_signs<B: DirectedBound>(
        &self,
        lower: &F<T>,
        upper: &F<T>,
        crossing: i64,
    ) -> Result<bool> {
        let lower = self.evaluate_interval(&B::enclose_native(lower), false)?.0;
        let upper = self.evaluate_interval(&B::enclose_native(upper), false)?.0;
        Ok(if crossing > 0 {
            lower[1] < B::ZERO && upper[0] > B::ZERO
        } else {
            lower[0] > B::ZERO && upper[1] < B::ZERO
        })
    }

    /// The compactified isolation view of the collapsed signed terms.
    fn compactified<B: DirectedBound>(
        &self,
        terms: &[(i64, &RadialEnergy<T>)],
    ) -> CompactifiedRay<B> {
        let add = B::enclosed_add;
        let square = B::enclosed_square;
        let native = B::enclose_native;
        let neg = B::enclosed_neg;
        let zero = [B::zero(), B::zero()];
        let energies = terms
            .iter()
            .enumerate()
            .map(|(index, &(coefficient, energy))| {
                let far_class = terms[..index]
                    .iter()
                    .position(|(_, other)| {
                        energy.velocity == other.velocity
                            || energy.velocity == -other.velocity.clone()
                    })
                    .unwrap_or(index);
                let origin = [
                    native(&energy.offset.px),
                    native(&energy.offset.py),
                    native(&energy.offset.pz),
                    native(&energy.mass),
                ];
                let end = [
                    native(&energy.velocity.px),
                    native(&energy.velocity.py),
                    native(&energy.velocity.pz),
                    zero.clone(),
                ];
                let direction: [_; 4] =
                    std::array::from_fn(|axis| add(&end[axis], &neg(&origin[axis])));
                let speed = B::enclosed_sqrt(
                    &direction
                        .iter()
                        .fold(zero.clone(), |sum, value| add(&sum, &square(value))),
                );
                CompactifiedEnergy {
                    coefficient,
                    origin,
                    end,
                    direction,
                    speed,
                    far_class,
                }
            })
            .collect();
        CompactifiedRay {
            energies,
            shift: native(&self.shift),
        }
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

/// Compactify the complete positive ray: r=u/(1-u), 0<u<1.
/// Each (1-u)E(r) is the Euclidean norm of the affine four-vector
/// ((1-u)b+uv, (1-u)m). Its convexity bounds its derivative by the
/// one-sided endpoint derivatives, even when the signed sum is not convex.
struct CompactifiedRay<B> {
    energies: Vec<CompactifiedEnergy<B>>,
    shift: [B; 2],
}

/// One collapsed term, whose affine four-vector runs from `origin` (b, m) at
/// u=0 to `end` (v, 0) at u=1.
struct CompactifiedEnergy<B> {
    coefficient: i64,
    origin: [[B; 2]; 4],
    end: [[B; 2]; 4],
    direction: [[B; 2]; 4],
    speed: [B; 2],
    /// The first collapsed term whose native velocity is +-this one. Terms of
    /// one far class have the same true radical |v| at u=1.
    far_class: usize,
}

/// The compactified value at one endpoint and, per term, its scaled left and
/// right slopes and affine four-vector. Pending cells share it through `Rc`
/// with their neighbours.
type EndpointData<B> = ([B; 2], Vec<([B; 2], [B; 2], [[B; 2]; 4])>);

/// The lower and upper endpoints of a cell or radial bracket holding one
/// simple root, and the root's crossing direction.
type Bracket<X> = (X, X, i64);

impl<B: DirectedBound> CompactifiedRay<B> {
    fn endpoint(&self, u: &B) -> std::result::Result<EndpointData<B>, &'static str> {
        let add = B::enclosed_add;
        let mul = B::enclosed_mul;
        let neg = B::enclosed_neg;
        let scale = B::enclosed_scale;
        let square = B::enclosed_square;
        let sqrt = B::enclosed_sqrt;
        let point = |value: &B| [value.clone(), value.clone()];
        let zero = point(&B::zero());
        let one = B::one();
        let fraction = point(u);
        let remainder = add(&point(&one), &neg(&fraction));
        let mut value = mul(&self.shift, &remainder);
        let mut data = Vec::with_capacity(self.energies.len());
        // Exact equal squared norms at an endpoint share one radical.
        // In particular, massive energies with equal velocities can cancel
        // at u=1. Do not replace such exact cancellation by a rounding band
        // which would make the excluded infinite endpoint look like a root.
        // At u=1 each compactified energy is |v|, so one far velocity class
        // shares one radical even when its |v|^2 has no exact enclosure.
        let at_infinity = *u == one;
        let mut shared_norms: Vec<(usize, [B; 2], i64)> = Vec::new();
        for CompactifiedEnergy {
            coefficient,
            origin,
            end,
            direction,
            speed,
            far_class,
        } in &self.energies
        {
            let vector: [_; 4] = std::array::from_fn(|axis| {
                add(&mul(&origin[axis], &remainder), &mul(&end[axis], &fraction))
            });
            let norm = vector
                .iter()
                .fold(zero.clone(), |sum, component| add(&sum, &square(component)));
            let energy = sqrt(&norm);
            let exact = norm[0] == norm[1];
            if exact || at_infinity {
                if let Some((_, _, previous)) =
                    shared_norms.iter_mut().find(|(class, shared, _)| {
                        (exact && shared[0] == shared[1] && shared[0] == norm[0])
                            || (at_infinity && class == far_class)
                    })
                {
                    *previous = previous
                        .checked_add(*coefficient)
                        .ok_or("energy multiplicity overflow")?;
                } else {
                    shared_norms.push((*far_class, norm, *coefficient));
                }
            } else {
                value = add(&value, &scale(energy.clone(), *coefficient));
            }
            let (left, right) = if energy[0] > B::ZERO {
                let numerator = vector
                    .iter()
                    .zip(direction)
                    .fold(zero.clone(), |sum, (q, v)| add(&sum, &mul(q, v)));
                let derivative = scale(B::enclosed_divide(&numerator, &energy), *coefficient);
                (derivative.clone(), derivative)
            } else if energy[1] == B::ZERO {
                // At an exact zero of a norm the left/right derivatives
                // are -/+ the affine speed; a constant zero has both zero.
                (
                    scale(neg(speed), *coefficient),
                    scale(speed.clone(), *coefficient),
                )
            } else {
                // A rounded enclosure touching zero has no resolved side.
                let bound = scale([speed[1].negated(), speed[1].clone()], *coefficient);
                (bound.clone(), bound)
            };
            data.push((left, right, vector));
        }
        for (_, norm, coefficient) in shared_norms {
            if coefficient != 0 {
                value = add(&value, &scale(sqrt(&norm), coefficient));
            }
        }
        // An infinite native bound surfaces here at u=0 or u=1, which isolation
        // evaluates before deciding any cell; a far class with zero net
        // coefficient still exposes it through its slopes. Vectors and cell
        // derivative folds need no check: infinite outward bounds are valid,
        // and enclosures of finite natives are never NaN.
        if value
            .iter()
            .chain(
                data.iter()
                    .flat_map(|(left, right, _)| left.iter().chain(right)),
            )
            .any(|number| !number.is_finite_bound())
        {
            return Err("nonfinite compactified energy enclosure");
        }
        Ok((value, data))
    }

    /// Certify every cell of the compactified half-line as empty or as holding
    /// exactly one simple crossing; returns the crossing cells in radial order.
    fn isolate(
        &self,
        work: &mut RadialIsolationWork,
    ) -> std::result::Result<Vec<Bracket<B>>, &'static str> {
        let add = B::enclosed_add;
        let neg = B::enclosed_neg;
        let (zero, one) = (B::zero(), B::one());
        let mut endpoint = |u: &B| {
            work.endpoints += 1;
            self.endpoint(u)
        };
        let mut pending = vec![(
            zero.clone(),
            one.clone(),
            Rc::new(endpoint(&zero)?),
            Rc::new(endpoint(&one)?),
        )];
        let mut brackets = Vec::new();
        while let Some((left, right, left_data, right_data)) = pending.pop() {
            work.visited += 1;
            if work.visited > B::CELL_BUDGET {
                return Err(
                    "cannot isolate all simple roots within the subdivision budget; tangent or unresolved geometry",
                );
            }
            let mut derivative = neg(&self.shift);
            let mut analytic = true;
            for ((energy, (_, left_slope, left_vector)), (right_slope, _, right_vector)) in
                self.energies.iter().zip(&left_data.1).zip(&right_data.1)
            {
                let CompactifiedEnergy {
                    coefficient,
                    origin,
                    end,
                    ..
                } = energy;
                let bound = if *coefficient > 0 {
                    [left_slope[0].clone(), right_slope[1].clone()]
                } else {
                    [right_slope[0].clone(), left_slope[1].clone()]
                };
                derivative = add(&derivative, &bound);
                // A nonzero affine component throughout the open cell excludes
                // a cusp. Constant zero energies are harmless. Analyticity lets
                // a monotone, nonconstant cell exclude a zero endpoint without
                // mistaking a piecewise-linear zero plateau for an absent cut.
                analytic &= origin.iter().chain(end).flatten().all(|x| *x == B::ZERO)
                    || left_vector.iter().zip(right_vector).any(|(a, b)| {
                        (a[0] >= B::ZERO && b[0] >= B::ZERO && (a[0] > B::ZERO || b[0] > B::ZERO))
                            || (a[1] <= B::ZERO
                                && b[1] <= B::ZERO
                                && (a[1] < B::ZERO || b[1] < B::ZERO))
                    });
            }
            let increasing = derivative[0] >= B::ZERO;
            let decreasing = derivative[1] <= B::ZERO;
            let left_value = &left_data.0;
            let right_value = &right_data.0;
            if (increasing && (left_value[0] > B::ZERO || right_value[1] < B::ZERO))
                || (decreasing && (left_value[1] < B::ZERO || right_value[0] > B::ZERO))
                || (analytic
                    && ((increasing
                        && ((left_value[0] >= B::ZERO && right_value[0] > B::ZERO)
                            || (right_value[1] <= B::ZERO && left_value[1] < B::ZERO)))
                        || (decreasing
                            && ((left_value[1] <= B::ZERO && right_value[1] < B::ZERO)
                                || (right_value[0] >= B::ZERO && left_value[0] > B::ZERO)))))
            {
                continue;
            }
            if analytic {
                let crossing = if derivative[0] > B::ZERO
                    && left_value[1] < B::ZERO
                    && right_value[0] > B::ZERO
                {
                    Some(1)
                } else if derivative[1] < B::ZERO
                    && left_value[0] > B::ZERO
                    && right_value[1] < B::ZERO
                {
                    Some(-1)
                } else {
                    None
                };
                if let Some(crossing) = crossing {
                    let mut upper = right.clone();
                    if upper == one {
                        // This whole cell already has exactly one root. Find a
                        // finite upper endpoint while retaining its certified
                        // lower endpoint; contracting both sides toward u=1
                        // can round a large root's lower radius onto the root.
                        let mut probe = left.clone();
                        loop {
                            let next = probe.bisect(&right);
                            work.visited += 1;
                            if next == probe || next == right || work.visited > B::CELL_BUDGET {
                                return Err("cannot separate a finite root from infinity");
                            }
                            let value = endpoint(&next)?.0;
                            if (crossing > 0 && value[0] > B::ZERO)
                                || (crossing < 0 && value[1] < B::ZERO)
                            {
                                upper = next;
                                break;
                            }
                            probe = next;
                        }
                    }
                    brackets.push((left, upper, crossing));
                    continue;
                }
            }
            let midpoint = left.bisect(&right);
            if midpoint == left || midpoint == right {
                return Err(
                    "cannot separate a tangent, cusp, or unresolved root at enclosure precision",
                );
            }
            let midpoint_data = endpoint(&midpoint)?;
            if midpoint_data.0[0] == B::ZERO && midpoint_data.0[1] == B::ZERO {
                let slope = midpoint_data
                    .1
                    .iter()
                    .fold(neg(&self.shift), |sum, (_, right, _)| add(&sum, right));
                if (slope[0] <= B::ZERO && slope[1] >= B::ZERO)
                    || midpoint_data.1.iter().any(|(left, right, _)| left != right)
                {
                    return Err("a tangent or vanishing-energy root is not a simple cut");
                }
            }
            let reach_radius = midpoint
                .add_up(&left.negated())
                .max_bound(&right.add_up(&midpoint.negated()));
            let maximum_derivative = derivative[0]
                .abs_bound()
                .max_bound(&derivative[1].abs_bound());
            let reach = maximum_derivative.mul_up(&reach_radius);
            let enclosure = add(&midpoint_data.0, &[reach.negated(), reach]);
            if enclosure[0] > B::ZERO || enclosure[1] < B::ZERO {
                continue;
            }
            if derivative[0] == B::ZERO && derivative[1] == B::ZERO {
                return Err("the radial equation is identically zero on a cell");
            }
            // Keep subdivision points away from roots. Otherwise an exact
            // dyadic root would remain an uncertified endpoint of both children.
            // Every accepted interior endpoint has a directed nonzero sign.
            let mut split = midpoint;
            let mut split_data = midpoint_data;
            while split_data.0[0] <= B::ZERO && split_data.0[1] >= B::ZERO {
                work.visited += 1;
                let next = left.bisect(&split);
                if next == left || next == split || work.visited > B::CELL_BUDGET {
                    return Err(
                        "cannot separate a tangent or nonisolated root from subdivision endpoints",
                    );
                }
                split = next;
                split_data = endpoint(&split)?;
            }
            let split_data = Rc::new(split_data);
            pending.push((split.clone(), right, Rc::clone(&split_data), right_data));
            pending.push((left, split, left_data, split_data));
        }
        // Enumeration finishes before evaluating any residue. Root identities
        // follow physical radial order and remain stable across precision lanes.
        brackets.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        Ok(brackets)
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
    use rand::{Rng, SeedableRng, rngs::StdRng};

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
            let interval = Float::enclose_native(&one);
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

        fn check<T: FloatLike>() {
            let native = |x| F(T::from_f64_exact_binary(x));
            let vector = |x| ThreeMomentum::new(native(x), native(0.0), native(0.0));
            let (zero, e_cm) = (native(0.0), native(2.0));
            // Exercise homogeneous and massive absorbed energies on the same chart.
            for negative_mass in [0.0, 0.25] {
                let repeated = RadialEnergy {
                    edge: EdgeIndex(4),
                    coefficient: 1,
                    velocity: vector(1.0),
                    offset: vector(0.0),
                    mass: native(1.0),
                };
                let ray = EnergySurfaceRay {
                    energies: vec![
                        repeated.clone(),
                        RadialEnergy {
                            edge: EdgeIndex(7),
                            coefficient: -1,
                            velocity: vector(1.9),
                            offset: vector(0.0),
                            mass: native(negative_mass),
                        },
                        repeated,
                    ],
                    shift: native(-1.0),
                };
                let retained = ray.energies.clone();
                let (roots, trace) = traced_roots(
                    &ray,
                    &e_cm,
                    &mut RadialRootDiagnostics::default(),
                    "shared signed ray",
                );
                // Every lane isolates and certifies this ray in binary64 alone.
                assert_eq!(trace.binary64, Binary64Outcome::Certified, "{trace:?}");
                assert_eq!((trace.work.visited, trace.work.endpoints), (0, 0));
                let roots = roots.unwrap();
                assert_eq!(roots.len(), 2);
                assert!(roots[0].derivative_at_solution < zero);
                assert!(roots[1].derivative_at_solution > zero);
                assert_eq!(ray.energies, retained);
                for root in &roots {
                    let epsilon = root.solution.epsilon();
                    let (residual, slope) = ray.evaluate(&root.solution);
                    assert!(residual.abs() <= &epsilon * native(64.0) * &e_cm);
                    assert_eq!(slope, root.derivative_at_solution);
                    let radius = HyperDual::<F<T>>::new(simple_n_deriv_shape(1))
                        .variable(0, root.solution.clone());
                    let jet = extract_t_derivatives(ray.evaluate_dual(&radius));
                    assert!((&jet[1] - &slope).abs() <= epsilon * native(32.0));
                }
                assert!(
                    ray.verify_lu_candidate(&ray, &roots[0], 6, &native(1.0e-6))
                        .unwrap_err()
                        .to_string()
                        .contains("positive E-surface rays")
                );
            }
        }
        check::<f64>();
        check::<QuadFloat>();
        check::<ArbPrec>();
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
        let interval = Float::enclose_native(&F(2.0));
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
        let interval = Float::enclose_native(&F(1.0));
        let (_, slope) = underflowing.evaluate_interval(&interval, true).unwrap();
        assert!(slope[0] > 0);
    }

    #[test]
    fn velocity_classes_cancel_inexact_norms() {
        // `tiny` is so far below the unit component that |v|^2 has no exact
        // enclosure at ENCLOSURE_PRECISION.
        fn check<T: FloatLike>(tiny: F<T>) {
            let one = tiny.one();
            let f = |n| one.from_i64(n);
            let velocity = ThreeMomentum::new(one.clone(), tiny, one.zero());
            let energy = |edge, coefficient, mass| RadialEnergy {
                edge: EdgeIndex(edge),
                coefficient,
                velocity: velocity.clone(),
                offset: ThreeMomentum::new(one.zero(), one.zero(), one.zero()),
                mass,
            };
            // eta(r) = sqrt(|v r|^2 + 4) - sqrt(|v r|^2 + 1) - 1/2 decreases
            // through one root at |v| r = sqrt(105)/4 and tends to -1/2, while
            // both compactified energies tend to the same radical |v| at u=1.
            let ray = EnergySurfaceRay {
                energies: vec![energy(0, 1, f(2)), energy(1, -1, one.clone())],
                shift: -(&one / f(2)),
            };
            let terms = uncollapsed_terms(&ray);
            // Each radical alone is a rounding band at u=1. Two opposite bands
            // would leave the excluded infinite endpoint unresolved; the shared
            // far-class radical cancels exactly.
            let infinity = Float::one();
            let lone = ray
                .compactified::<Float>(&terms[..1])
                .endpoint(&infinity)
                .unwrap()
                .0;
            assert!(lone[0] < lone[1], "{lone:?}");
            let compactified = ray.compactified::<Float>(&terms);
            let far = compactified.endpoint(&infinity).unwrap().0;
            assert!(far.iter().all(Float::is_zero), "{far:?}");

            // The reference isolation itself resolves the cancelling far endpoint.
            let brackets = compactified
                .isolate(&mut RadialIsolationWork::default())
                .unwrap();
            assert_eq!(brackets.len(), 1);
            assert_eq!(brackets[0].2, -1);
            // The default path already resolves this far class in binary64,
            // where |v|^2 is a rounding band as well; the reference never runs.
            let (roots, trace) = traced_roots(
                &ray,
                &f(2),
                &mut RadialRootDiagnostics::default(),
                "inexact velocity class",
            );
            assert_eq!(trace.path, RadialIsolationPath::Isolated);
            assert_eq!(trace.binary64, Binary64Outcome::Certified, "{trace:?}");
            assert_eq!(trace.terms, 2);
            assert_eq!((trace.work.visited, trace.work.endpoints), (0, 0));
            assert_eq!((trace.binary64_certificates, trace.certificates), (1, 0));
            let roots = roots.unwrap();
            assert_eq!(roots.len(), 1);
            let expected = f(105).sqrt() / f(4);
            assert!((&roots[0].solution - expected).abs() < one.epsilon() * f(1024));
            assert!(roots[0].derivative_at_solution < one.zero());
        }
        // Components 2^1074 apart need 2149 bits for |v|^2 in every lane, and
        // 2^-1100 widens to [-2^-1074, 2^-1074] in binary64.
        let tiny = f64::from_bits(1);
        check(F(tiny));
        check(F(QuadFloat::from_f64_exact_binary(tiny)));
        let one = F::<ArbPrec>::default().one();
        check(one.from_i64(2).powi(-1100));
    }

    /// A test energy: coefficient, velocity, offset and mass.
    type TestEnergy = (i64, [f64; 3], [f64; 3], f64);

    /// A ray of energies whose binary64 coordinates embed exactly in lane T.
    fn exact_ray<T: FloatLike>(energies: &[TestEnergy], shift: f64) -> EnergySurfaceRay<T> {
        let native = |x: f64| F(T::from_f64_exact_binary(x));
        let vector = |[x, y, z]: [f64; 3]| ThreeMomentum::new(native(x), native(y), native(z));
        EnergySurfaceRay {
            energies: energies
                .iter()
                .enumerate()
                .map(
                    |(edge, &(coefficient, velocity, offset, mass))| RadialEnergy {
                        edge: EdgeIndex(edge),
                        coefficient,
                        velocity: vector(velocity),
                        offset: vector(offset),
                        mass: native(mass),
                    },
                )
                .collect(),
            shift: native(shift),
        }
    }

    const ORIGIN: [f64; 3] = [0.0; 3];

    fn along_x(value: f64) -> [f64; 3] {
        [value, 0.0, 0.0]
    }

    /// Every energy of a ray as its own isolation term.
    fn uncollapsed_terms<T: FloatLike>(ray: &EnergySurfaceRay<T>) -> Vec<(i64, &RadialEnergy<T>)> {
        ray.energies
            .iter()
            .map(|energy| (energy.coefficient, energy))
            .collect()
    }

    /// The cells of one isolation pass.
    type Isolation<B> = std::result::Result<Vec<Bracket<B>>, &'static str>;

    /// Both isolation passes over the uncollapsed terms of a ray.
    fn both_isolations<T: FloatLike>(
        ray: &EnergySurfaceRay<T>,
    ) -> (Isolation<f64>, Isolation<Float>) {
        let terms = uncollapsed_terms(ray);
        (
            ray.compactified::<f64>(&terms)
                .isolate(&mut RadialIsolationWork::default()),
            ray.compactified::<Float>(&terms)
                .isolate(&mut RadialIsolationWork::default()),
        )
    }

    fn traced_roots<T: FloatLike>(
        ray: &EnergySurfaceRay<T>,
        e_cm: &F<T>,
        diagnostics: &mut RadialRootDiagnostics,
        identity: &str,
    ) -> (
        eyre::Result<Vec<NewtonIterationResult<T>>>,
        RadialIsolationTrace,
    ) {
        let mut trace = RadialIsolationTrace::default();
        let roots = ray.enumerate_positive_radius_roots(
            e_cm,
            diagnostics,
            &RadialRootIdentity::new(identity.into()),
            &mut trace,
        );
        (roots, trace)
    }

    /// A reference bracket at the native radii the reference pass materializes.
    fn materialized((left, right, crossing): &Bracket<Float>) -> Bracket<F<f64>> {
        let radius = |u| EnergySurfaceRay::<f64>::materialized_radius(u).unwrap();
        (radius(left), radius(right), *crossing)
    }

    /// Whether NaN-free binary64 bounds contain a reference enclosure.
    fn contains(outer: [f64; 2], inner: [Float; 2]) -> bool {
        !outer[0].is_nan() && !outer[1].is_nan() && outer[0] <= inner[0] && inner[1] <= outer[1]
    }

    #[test]
    fn binary64_bounds_enclose_exact_results() {
        assert!(f64::available());
        // A directed reference enclosure of the exact result: containment
        // also makes every binary64 point exact.
        let check = |operation: &str, bound: [f64; 2], reference: [Float; 2]| {
            assert!(
                contains(bound, reference.clone()),
                "{operation}: {bound:?} does not contain {reference:?}"
            );
        };
        let lift = |x: f64| Float::with_val(53, x);
        let tiny = f64::from_bits(1);
        let mut values = vec![
            0.0,
            tiny,
            3.0 * tiny,
            f64::MIN_POSITIVE.next_down(),
            f64::MIN_POSITIVE,
            2.0_f64.powi(-537),
            1.0e-300,
            0.1,
            1.0 / 3.0,
            0.5,
            1.0_f64.next_down(),
            1.0,
            1.0_f64.next_up(),
            2.0,
            3.0,
            2.0_f64.powi(53) + 2.0,
            1.0e300,
            f64::MAX,
        ];
        values.extend(values.clone().into_iter().map(|x| -x));
        for &a in &values {
            for &b in &values {
                // Exact binary64 sums need at most 2098 bits.
                let sum = Float::with_val(2200, &lift(a) + &lift(b));
                check(
                    &format!("{a:e} + {b:e}"),
                    [a.add_down(&b), a.add_up(&b)],
                    [sum.clone(), sum],
                );
                let product = Float::with_val(2200, &lift(a) * &lift(b));
                check(
                    &format!("{a:e} * {b:e}"),
                    [a.mul_down(&b), a.mul_up(&b)],
                    [product.clone(), product],
                );
            }
            if a > 0.0 {
                let reference = [Round::Down, Round::Up]
                    .map(|round| Float::with_val_round(2048, lift(a).recip_ref(), round).0);
                check(
                    &format!("1 / {a:e}"),
                    [a.recip_down(), a.recip_up()],
                    reference,
                );
            }
            if a >= 0.0 {
                let reference = [Round::Down, Round::Up]
                    .map(|round| Float::with_val_round(2048, lift(a).sqrt_ref(), round).0);
                check(
                    &format!("sqrt {a:e}"),
                    [a.sqrt_down(), a.sqrt_up()],
                    reference,
                );
            }
        }
        // A zero factor beside an unbounded bound is an exact zero, while NaN
        // from a nonfinite native still propagates.
        for unbounded in [f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                [0.0_f64.mul_down(&unbounded), unbounded.mul_up(&-0.0)],
                [0.0; 2]
            );
        }
        assert!(0.0_f64.mul_down(&f64::NAN).is_nan());
        // Overflow stays outward, and lower bounds never reach +inf, so no
        // operation meets inf-inf.
        assert_eq!(
            [f64::MAX.add_down(&f64::MAX), f64::MAX.add_up(&f64::MAX)],
            [f64::MAX, f64::INFINITY]
        );
        assert_eq!(
            [
                f64::NEG_INFINITY.add_down(&f64::MAX),
                f64::INFINITY.add_up(&-f64::MAX)
            ],
            [f64::NEG_INFINITY, f64::INFINITY]
        );
        for (value, exact) in [
            (0, true),
            (-1, true),
            (7, true),
            (1 << 53, true),
            ((1 << 53) + 1, false),
            (-(1 << 53) - 1, false),
            (i64::MAX, false),
            (i64::MIN, true),
        ] {
            let bound = f64::enclose_multiplicity(value);
            let reference = Float::with_val(64, value);
            check(&format!("{value}"), bound, [reference.clone(), reference]);
            assert_eq!(bound[0] == bound[1], exact, "{value}: {bound:?}");
        }
        for (left, right) in [
            (0.0, 1.0),
            (0.5, 0.75),
            (1.0_f64.next_down(), 1.0),
            (0.0, tiny),
        ] {
            let middle = left.bisect(&right);
            assert!(left <= middle && middle <= right);
        }

        // Every binary64 interval operation contains its reference result, so
        // whatever binary64 bounds certify, the reference certifies as well.
        let points = [
            -3.0,
            -1.0,
            -0.1,
            0.0,
            1.0e-300,
            1.0 / 3.0,
            1.0,
            2.5,
            1.0e300,
        ];
        let intervals = points
            .iter()
            .flat_map(|&a| {
                points
                    .iter()
                    .filter(move |&&b| a <= b)
                    .map(move |&b| [a, b])
            })
            .collect_vec();
        let reference = |interval: &[f64; 2]| interval.map(|x| Float::with_val(2048, x));
        for a in &intervals {
            assert!(contains(
                f64::enclosed_square(a),
                Float::enclosed_square(&reference(a))
            ));
            if a[0] >= 0.0 {
                assert!(contains(
                    f64::enclosed_sqrt(a),
                    Float::enclosed_sqrt(&reference(a))
                ));
            }
            for b in &intervals {
                assert!(contains(
                    f64::enclosed_add(a, b),
                    Float::enclosed_add(&reference(a), &reference(b))
                ));
                assert!(contains(
                    f64::enclosed_mul(a, b),
                    Float::enclosed_mul(&reference(a), &reference(b))
                ));
                if b[0] > 0.0 {
                    assert!(contains(
                        f64::enclosed_divide(a, b),
                        Float::enclosed_divide(&reference(a), &reference(b))
                    ));
                }
            }
        }
    }

    #[test]
    fn sign_case_enclosure_product_matches_all_corners() {
        // Every ordered interval of these bounds: all nine sign classes of
        // two operands, points, and zeros of either sign.
        fn intervals<B: DirectedBound>(bounds: &[B]) -> Vec<[B; 2]> {
            bounds
                .iter()
                .flat_map(|a| {
                    bounds
                        .iter()
                        .filter(move |&b| a <= b)
                        .map(move |b| [a.clone(), b.clone()])
                })
                .collect()
        }
        // The all-corner products and magnitude squares that the sign cases
        // replace.
        fn corners<B: DirectedBound>(a: &[B; 2], b: &[B; 2]) -> [B; 2] {
            let lower = a
                .iter()
                .flat_map(|x| b.iter().map(move |y| x.mul_down(y)))
                .reduce(|a, b| a.min_bound(&b))
                .unwrap();
            let upper = a
                .iter()
                .flat_map(|x| b.iter().map(move |y| x.mul_up(y)))
                .reduce(|a, b| a.max_bound(&b))
                .unwrap();
            [lower, upper]
        }
        fn square<B: DirectedBound>(a: &[B; 2]) -> [B; 2] {
            let upper = a[0].abs_bound().max_bound(&a[1].abs_bound());
            let lower = if a[0] <= B::ZERO && a[1] >= B::ZERO {
                B::zero()
            } else {
                a[0].abs_bound().min_bound(&a[1].abs_bound())
            };
            [lower.mul_down(&lower), upper.mul_up(&upper)]
        }

        // Reference bounds agree bitwise, up to the sign of zero. Thirds and
        // sevenths make the directed products inexact.
        let precision = ENCLOSURE_PRECISION;
        let magnitudes = [
            Float::with_val(precision, 0),
            Float::with_val(precision, 1) / 7,
            Float::with_val(precision, 1) / 3,
            Float::with_val(precision, 1),
            Float::with_val(precision, 7) / 3,
        ];
        let bounds = magnitudes
            .iter()
            .flat_map(|x| [x.clone(), -x.clone()])
            .collect_vec();
        let reference = intervals(&bounds);
        for a in &reference {
            assert_eq!(Float::enclosed_square(a), square(a), "{a:?}^2");
            for b in &reference {
                assert_eq!(Float::enclosed_mul(a, b), corners(a, b), "{a:?} * {b:?}");
            }
        }

        // Binary64 bounds contain the exact products. Binary64 corners are
        // exact at the reference precision; MPFR makes a zero beside an
        // unbounded bound NaN, which the extrema skip, and a zero point
        // times any enclosure is exactly zero.
        let tiny = f64::from_bits(1);
        let mut bounds = vec![0.0, tiny, 1.0e-300, 1.0 / 3.0, 1.0, 1.5, 1.0e300, f64::MAX];
        bounds.extend(bounds.clone().into_iter().map(|x| -x));
        bounds.extend([f64::INFINITY, f64::NEG_INFINITY]);
        let lift = |x: &[f64; 2]| x.map(|x| Float::with_val(precision, x));
        let exact = |a: &[f64; 2], b: &[f64; 2]| {
            if a == &[0.0; 2] || b == &[0.0; 2] {
                lift(&[0.0; 2])
            } else {
                corners(&lift(a), &lift(b))
            }
        };
        // Lower bounds never reach +inf and upper bounds never reach -inf.
        let binary64 = intervals(&bounds)
            .into_iter()
            .filter(|[a, b]| *a < f64::INFINITY && *b > f64::NEG_INFINITY)
            .collect_vec();
        for a in &binary64 {
            assert!(
                contains(f64::enclosed_square(a), square(&lift(a))),
                "{a:?}^2"
            );
            for b in &binary64 {
                assert!(
                    contains(f64::enclosed_mul(a, b), exact(a, b)),
                    "{a:?} * {b:?}"
                );
            }
        }
    }

    #[test]
    fn binary64_native_enclosures_are_exact_or_outward() {
        fn check<T: FloatLike>(value: F<T>, point: bool) -> [f64; 2] {
            let bound = f64::enclose_native(&value);
            assert!(
                contains(bound, Float::enclose_native(&value)),
                "{value}: {bound:?}"
            );
            assert_eq!(bound[0] == bound[1], point, "{value}: {bound:?}");
            bound
        }
        let tiny = f64::from_bits(1);
        for value in [0.0, -0.0, tiny, 0.1, -f64::MAX] {
            check(F(value), true);
        }
        let one = F::<QuadFloat>::default().one();
        check(one.from_i64(3), true);
        // A second limb below binary64 resolution widens by one ulp.
        assert_eq!(
            check(one + one.from_i64(2).powi(-80), false),
            [1.0_f64.next_down(), 1.0_f64.next_up()]
        );
        let one = F::<ArbPrec>::default().one();
        check(&one / one.from_i64(3), false);
        // Magnitudes beyond binary64 range widen to the outward limits.
        assert_eq!(check(one.from_i64(2).powi(-1100), false), [-tiny, tiny]);
        let huge = one.from_i64(10).powi(400);
        assert_eq!(check(huge.clone(), false), [f64::MAX, f64::INFINITY]);
        assert_eq!(check(-huge, false), [f64::NEG_INFINITY, -f64::MAX]);
        // Nonfinite natives have no enclosure.
        for value in [f64::NAN, f64::INFINITY] {
            assert!(f64::enclose_native(&F(value)).iter().all(|x| x.is_nan()));
        }
    }

    #[test]
    fn binary64_radii_match_reference_materialization() {
        let lifted = |u: f64| Float::with_val(ENCLOSURE_PRECISION, u);
        // Dyadic u near 0, 1/2 and 1 at depths 1-53, where 1-u is a binary64
        // value, and at depths 54-70, where it is not.
        let mut corpus = vec![0.0];
        for depth in 1..=70 {
            let top = 2.0_f64.powi(depth.min(53));
            let numerators = [1.0, 3.0, top / 2.0 + 1.0, top - 3.0, top - 1.0];
            corpus.extend(
                numerators
                    .map(|numerator| numerator * 2.0_f64.powi(-depth))
                    .into_iter()
                    .filter(|u| (0.0..1.0).contains(u)),
            );
        }
        // Random 53-bit u in [0, 1), and in [0, 1/2) on a finer grid, where 1-u
        // is a binary64 value only for an even numerator.
        let mut rng = StdRng::seed_from_u64(0x9e37_79b9_7f4a_7c15);
        for _ in 0..2048 {
            let numerator = (rng.random::<u64>() >> 11) as f64;
            corpus.extend([numerator * 2.0_f64.powi(-53), numerator * 2.0_f64.powi(-54)]);
        }
        let mut divided = 0;
        for u in corpus {
            let reference = EnergySurfaceRay::<f64>::materialized_radius(&lifted(u)).unwrap();
            let radius = EnergySurfaceRay::<f64>::binary64_radius(u).unwrap();
            assert_eq!(radius.0.to_bits(), reference.0.to_bits(), "u = {u:e}");
            if 1.0 - (1.0 - u) == u {
                divided += 1;
                assert_eq!(
                    (u / (1.0 - u)).to_bits(),
                    reference.0.to_bits(),
                    "u = {u:e}"
                );
            }
            // Other lanes keep the reference rounding to their own precision.
            assert_eq!(
                EnergySurfaceRay::<QuadFloat>::binary64_radius(u).unwrap(),
                EnergySurfaceRay::<QuadFloat>::materialized_radius(&lifted(u)).unwrap(),
                "u = {u:e}"
            );
        }
        assert!(divided > 2048, "{divided}");
    }

    #[test]
    fn binary64_isolation_agrees_with_reference() {
        /// The roots of one lane, as binary64 radii with their crossing
        /// directions, counting isolated and binary64-certified calls.
        fn lane_roots<T: FloatLike>(
            energies: &[TestEnergy],
            shift: f64,
            [isolated, certified]: &mut [usize; 2],
        ) -> Option<Vec<(f64, bool)>> {
            let ray = exact_ray::<T>(energies, shift);
            let (roots, trace) = traced_roots(
                &ray,
                &F::<T>::default().one(),
                &mut RadialRootDiagnostics::default(),
                "binary64 agreement",
            );
            if trace.path == RadialIsolationPath::Isolated {
                *isolated += 1;
                if trace.binary64 == Binary64Outcome::Certified {
                    *certified += 1;
                    // Only the native solve can still fail, at the same
                    // residual floor as from the reference brackets. Its
                    // start conditions were checked before any solve.
                    if let Err(error) = &roots {
                        let error = error.to_string();
                        assert!(
                            error.contains("isolated root refinement failed")
                                && !error.contains("InvalidInside")
                                && !error.contains("InvalidOutside"),
                            "{ray:?}: {error}"
                        );
                    }
                }
            }
            let roots = roots.ok()?;
            Some(
                roots
                    .iter()
                    .map(|root| {
                        (
                            root.solution.into_f64(),
                            root.derivative_at_solution > root.solution.zero(),
                        )
                    })
                    .collect(),
            )
        }

        let mut rng = StdRng::seed_from_u64(0x9e37_79b9_7f4a_7c15);
        let mut draw = |scale: f64| scale * rng.random_range(-1.0..1.0);
        let add = |a: [f64; 3], b: [f64; 3]| [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
        let neg = |a: [f64; 3]| a.map(|value| -value);
        let mut thermal_a = Vec::new();
        let mut thermal_b = Vec::new();
        let mut two_loop = Vec::new();
        let mut signed = Vec::new();
        for _ in 0..200 {
            let v = [draw(1.0), draw(1.0), draw(1.0)];
            let p = [draw(1.0), draw(1.0), draw(1.0)];
            let (m1, m2, p0) = (draw(1.0).abs(), draw(1.0).abs(), draw(2.0));
            // E(k+p) - E(k) = p0 and E(k+p) - E(k-p) = p0.
            thermal_a.push((vec![(1, v, p, m1), (-1, v, ORIGIN, m2)], -p0));
            thermal_b.push((vec![(1, v, p, m1), (-1, v, neg(p), m1)], -p0));
            let w = [draw(1.0), draw(1.0), draw(1.0)];
            let c = [draw(1.0), draw(1.0), draw(1.0)];
            two_loop.push((
                vec![
                    (1, v, c, m1),
                    (1, w, p, m2),
                    (-1, add(v, w), add(add(c, p), along_x(p0)), m1 * m2),
                ],
                -p0.abs(),
            ));
            let count = 2 + (draw(1.0).abs() * 3.0) as usize;
            let energies = (0..count)
                .map(|index| {
                    (
                        if index == 0 || draw(1.0) < 0.0 { -1 } else { 1 },
                        [draw(1.0), draw(1.0), draw(1.0)],
                        [draw(1.0), draw(1.0), draw(1.0)],
                        draw(1.0).abs(),
                    )
                })
                .collect_vec();
            signed.push((energies, draw(2.0)));
        }
        // The signed hsurface fixtures, with each incoming energy as -shift.
        let fixtures = vec![
            (
                vec![
                    (1, along_x(1.0), ORIGIN, 0.0),
                    (1, [0.0, 1.0, 0.0], ORIGIN, 0.0),
                    (-1, [-1.0, -1.0, 0.0], [0.0, 0.0, 0.5], 0.0),
                ],
                -2.0,
            ),
            (
                vec![
                    (1, along_x(1.0), ORIGIN, 1.0),
                    (1, along_x(1.0), ORIGIN, 1.0),
                    (-1, along_x(1.9), ORIGIN, 0.0),
                ],
                -1.0,
            ),
            (
                vec![
                    (1, along_x(1.0), ORIGIN, 1.0),
                    (1, along_x(1.0), ORIGIN, 1.0),
                    (-1, along_x(1.9), ORIGIN, 0.0),
                ],
                0.0,
            ),
            (vec![(-1, along_x(1.0), ORIGIN, 0.0)], -1.0),
            (
                vec![
                    (1, along_x(1.0), ORIGIN, 1.0),
                    (-1, along_x(2.0), ORIGIN, 1.0),
                ],
                1.0,
            ),
            (
                vec![
                    (1, along_x(1.0), ORIGIN, 1.8),
                    (-1, along_x(1.0), ORIGIN, 0.3),
                ],
                -1.0,
            ),
            (
                vec![
                    (1, along_x(1.0), [0.0, 0.0, 2.0], 0.5),
                    (-1, along_x(1.0), ORIGIN, 0.5),
                ],
                -1.0,
            ),
            (
                vec![
                    (1, along_x(1.0), ORIGIN, 1.0),
                    (-1, along_x(2.0), ORIGIN, 2.0),
                    (1, along_x(3.0), ORIGIN, 3.0),
                ],
                -4.0,
            ),
            (
                vec![
                    (1, along_x(1.0), ORIGIN, 1.0),
                    (-1, along_x(1.0), ORIGIN, 2.0),
                ],
                -1.0,
            ),
        ];
        // The thermal and two-loop families must stay on the binary64 path in
        // every lane, and so must every fixture, including its equal-velocity
        // rays.
        for (family, minimum_rate, rays) in [
            ("E(k+p)-E(k)", 0.95, thermal_a),
            ("E(k+p)-E(k-p)", 0.95, thermal_b),
            ("two-loop", 0.95, two_loop),
            ("signed", 0.0, signed),
            ("hsurface fixtures", 1.0, fixtures),
        ] {
            let mut counts = [[0; 2]; 3];
            for (energies, shift) in &rays {
                let ray = exact_ray::<f64>(energies, *shift);
                let (fast, reference) = both_isolations(&ray);
                if let Ok(fast) = &fast {
                    let reference = reference
                        .as_ref()
                        .unwrap_or_else(|error| panic!("{family} {ray:?}: {error}"));
                    assert_eq!(fast.len(), reference.len(), "{family} {ray:?}");
                    for ((left, right, crossing), (lower, upper, reference_crossing)) in
                        fast.iter().zip(reference)
                    {
                        // Both brackets hold the same unique root.
                        assert_eq!(crossing, reference_crossing, "{family} {ray:?}");
                        assert!(*left < *upper && *lower < *right, "{family} {ray:?}");
                    }
                }
                let double = lane_roots::<f64>(energies, *shift, &mut counts[0]);
                // Wherever two lanes refine, they find the same roots.
                for (lane, roots) in [
                    (
                        "Quad",
                        lane_roots::<QuadFloat>(energies, *shift, &mut counts[1]),
                    ),
                    (
                        "Arb",
                        lane_roots::<ArbPrec>(energies, *shift, &mut counts[2]),
                    ),
                ] {
                    let (Some(double), Some(roots)) = (&double, roots) else {
                        continue;
                    };
                    assert_eq!(roots.len(), double.len(), "{family} {lane} {ray:?}");
                    for ((root, ascending), (expected, expected_ascending)) in
                        roots.into_iter().zip(double)
                    {
                        assert_eq!(ascending, *expected_ascending, "{family} {lane} {ray:?}");
                        assert!(
                            (root - expected).abs() <= 1.0e-9 * expected.max(1.0),
                            "{family} {lane} {ray:?}: {root} against {expected}"
                        );
                    }
                }
            }
            for (lane, [isolated, certified]) in ["Double", "Quad", "Arb"].into_iter().zip(counts) {
                assert!(
                    certified as f64 >= minimum_rate * isolated as f64,
                    "{family} {lane}: {certified}/{isolated} certified in binary64"
                );
            }
        }
    }

    #[test]
    fn forced_fallbacks_report_reference_outcomes() {
        let fallback = |trace: &RadialIsolationTrace| {
            assert_eq!(trace.path, RadialIsolationPath::Isolated, "{trace:?}");
            assert!(
                !matches!(
                    trace.binary64,
                    Binary64Outcome::Certified | Binary64Outcome::NotAttempted
                ),
                "{trace:?}"
            );
            assert!(trace.work.endpoints > 0, "{trace:?}");
        };
        let diagnostics = &mut RadialRootDiagnostics::default();

        // Squares of 2^-665-scaled coordinates underflow in binary64; the
        // reference keeps both roots of the unscaled ray.
        let two_roots = [
            (1, along_x(1.0), ORIGIN, 1.0),
            (1, along_x(1.0), ORIGIN, 1.0),
            (-1, along_x(1.9), ORIGIN, 0.0),
        ];
        let expected = exact_ray::<f64>(&two_roots, -1.0)
            .positive_radius_roots(
                &F(2.0),
                diagnostics,
                &RadialRootIdentity::new("unscaled".into()),
            )
            .unwrap();
        let scale = 2.0_f64.powi(-665);
        let scaled = two_roots.map(|(coefficient, velocity, offset, mass)| {
            (
                coefficient,
                velocity.map(|v| v * scale),
                offset,
                mass * scale,
            )
        });
        let (roots, trace) = traced_roots(
            &exact_ray::<ArbPrec>(&scaled, -scale),
            &F::<ArbPrec>::default().one(),
            diagnostics,
            "scaled",
        );
        fallback(&trace);
        let roots = roots.unwrap();
        assert_eq!(roots.len(), expected.len());
        for (root, expected) in roots.iter().zip(&expected) {
            assert!((root.solution.into_f64() / expected.solution.0 - 1.0).abs() < 1.0e-12);
        }

        // Velocities one ulp apart leave the far endpoint unresolved in
        // binary64; the reference proves the channel closed. Below binary64
        // resolution, a Quad velocity mismatch does the same.
        let near_equal = |velocity| {
            [
                (1, along_x(1.0), ORIGIN, 1.0),
                (-1, along_x(velocity), ORIGIN, 2.0),
            ]
        };
        let (roots, trace) = traced_roots(
            &exact_ray::<f64>(&near_equal(1.0_f64.next_up()), -0.5),
            &F(2.0),
            diagnostics,
            "near-equal velocities",
        );
        fallback(&trace);
        assert!(roots.unwrap().is_empty());
        let mut quad = exact_ray::<QuadFloat>(&near_equal(1.0), -0.5);
        let one = F::<QuadFloat>::default().one();
        quad.energies[1].velocity.px = one + one.from_i64(2).powi(-80);
        let (roots, trace) = traced_roots(&quad, &one.from_i64(2), diagnostics, "Quad velocities");
        fallback(&trace);
        assert!(roots.unwrap().is_empty());

        // Roots beyond binary64 compactified resolution. Binary64 isolation
        // fails before any bracket materializes, so each root is solved once.
        for velocity in [1.0e-17, 2.0_f64.powi(-200)] {
            let identity = format!("far root {velocity:e}");
            let (roots, trace) = traced_roots(
                &exact_ray::<f64>(&[(-1, along_x(velocity), ORIGIN, 0.0)], 1.0),
                &F(2.0),
                diagnostics,
                &identity,
            );
            fallback(&trace);
            assert!(
                matches!(trace.binary64, Binary64Outcome::Uncertain(_)),
                "{trace:?}"
            );
            let roots = roots.unwrap();
            assert_eq!(roots.len(), 1);
            assert!((roots[0].solution.0 * velocity - 1.0).abs() < 1.0e-12);
            assert_eq!(
                diagnostics.solve_counts(&RadialRootIdentity::new(format!(
                    "{identity} signed root 0"
                ))),
                (1, 1)
            );
        }
        // An Arb velocity beyond binary64 range encloses to [MAX, inf], which
        // reaches the finiteness check of the compactified endpoints.
        let one = F::<ArbPrec>::default().one();
        let huge = one.from_i64(10).powi(400);
        let mut ray = exact_ray::<ArbPrec>(&[(-1, along_x(1.0), ORIGIN, 0.0)], 1.0);
        ray.energies[0].velocity.px = huge.clone();
        let (roots, trace) = traced_roots(&ray, &one, diagnostics, "huge Arb velocity");
        fallback(&trace);
        assert_eq!(
            trace.binary64,
            Binary64Outcome::Uncertain("nonfinite compactified energy enclosure")
        );
        let roots = roots.unwrap();
        assert_eq!(roots.len(), 1);
        assert!((&roots[0].solution * &huge - &one).abs() < one.epsilon() * one.from_i64(1024));

        // A root at r=2^-1000 is certified in both passes, but its inside
        // value lies within the solver's residual margin: the reference
        // reports the refinement failure.
        let (roots, trace) = traced_roots(
            &exact_ray::<f64>(&[(-1, along_x(1.0), ORIGIN, 0.0)], 2.0_f64.powi(-1000)),
            &F(2.0),
            diagnostics,
            "boundary root",
        );
        fallback(&trace);
        assert_eq!(trace.binary64, Binary64Outcome::Unrefinable, "{trace:?}");
        let error = roots.unwrap_err().to_string();
        assert!(error.contains("isolated root refinement failed"), "{error}");

        // Tangent and cusp rays keep the reference diagnosis.
        for (energy, shift) in [
            ((1, along_x(1.0), along_x(-2.0), 1.0), -1.0),
            ((1, along_x(1.0), along_x(-1.0), 0.0), 0.0),
        ] {
            let (roots, trace) = traced_roots(
                &exact_ray::<f64>(&[energy], shift),
                &F(2.0),
                diagnostics,
                "tangent",
            );
            fallback(&trace);
            let error = roots.unwrap_err().to_string();
            assert!(error.contains("tangent"), "{error}");
        }
    }

    #[test]
    fn binary64_brackets_can_refine_where_reference_brackets_cannot() {
        // Subdivision trees need not nest: where a binary64 value band at a
        // split point straddles zero, the split moves, and so can the brackets.
        // Both outcomes below are sound; only the precision escalation differs.
        let certified_roots = |ray: &EnergySurfaceRay<f64>| {
            let (roots, trace) = traced_roots(
                ray,
                &F(1.0),
                &mut RadialRootDiagnostics::default(),
                "binary64 brackets",
            );
            assert_eq!(trace.binary64, Binary64Outcome::Certified, "{trace:?}");
            roots.unwrap()
        };

        // eta(r) = sqrt((r-2)^2 + 2^-52) - 1 has roots at 2 -+ sqrt(1-2^-52).
        // At u=1/2 (r=1) the reference resolves eta ~ 2^-53 > 0 and brackets
        // the first root by [1/2, 5/8], but the native value at r=1 rounds to
        // -0, inside the solver's residual margin: refining fails there, and
        // the caller escalates precision. The binary64 value band at u=1/2
        // straddles zero, so its split moves to u=1/4, and its brackets refine
        // in Double.
        let ray = exact_ray::<f64>(&[(1, along_x(1.0), along_x(-2.0), 2.0_f64.powi(-26))], -1.0);
        let (fast, reference) = both_isolations(&ray);
        assert_eq!(fast.unwrap(), [(0.25, 0.625, -1), (0.71875, 0.8125, 1)]);
        let reference = &reference.unwrap()[0];
        assert_eq!(
            (reference.0.to_f64(), reference.1.to_f64(), reference.2),
            (0.5, 0.625, -1)
        );
        let (lower, upper, crossing) = materialized(reference);
        let error = ray
            .refine_root(
                &lower,
                &upper,
                crossing,
                0,
                &F(1.0),
                &mut RadialRootDiagnostics::default(),
                &RadialRootIdentity::new("reference bracket".into()),
            )
            .unwrap_err();
        assert!(
            matches!(error, SafeguardedNewtonError::InvalidInside { .. }),
            "{error:?}"
        );
        let roots = certified_roots(&ray);
        assert_eq!(roots.len(), 2);
        for (root, expected) in roots.iter().zip([1.0, 3.0]) {
            assert!((root.solution.0 - expected).abs() < 1.0e-12, "{roots:?}");
        }

        // E(vr+b1) - E(vr+b2) with decimal momenta has E1^2-E2^2 = 3.48-0.36r,
        // so its root lies at r=29/3, i.e. u=29/32, up to input rounding. The
        // reference brackets it by [115/128, 29/32], whose materialized upper
        // radius is not sign-certified, so the reference pass fails and the
        // caller escalates precision. The binary64 bracket [231/256, 117/128]
        // is certified and refines in Double.
        let velocity = [0.5, 0.0, 0.3];
        let ray = exact_ray::<f64>(
            &[
                (1, velocity, [-1.8, 1.9, -2.6], 0.1),
                (-1, velocity, [-1.8, 1.7, -2.0], 0.1),
            ],
            0.0,
        );
        let (fast, reference) = both_isolations(&ray);
        assert_eq!(fast.unwrap(), [(0.90234375, 0.9140625, -1)]);
        let reference = &reference.unwrap()[0];
        assert_eq!(
            (reference.0.to_f64(), reference.1.to_f64(), reference.2),
            (0.8984375, 0.90625, -1)
        );
        let (lower, upper, crossing) = materialized(reference);
        let mut trace = RadialIsolationTrace::default();
        assert!(
            !ray.bracket_signs_certified(&lower, &upper, crossing, true, &mut trace)
                .unwrap()
        );
        let roots = certified_roots(&ray);
        assert_eq!(roots.len(), 1);
        assert!(
            (roots[0].solution.0 - 29.0 / 3.0).abs() < 1.0e-12,
            "{roots:?}"
        );
    }

    /// eta(r) = sqrt((256r-261)^2+144) - s with s just above 2729/13 has two
    /// roots, and its compactified slope at u=1/2 is s-2729/13 > 0, below
    /// binary64 resolution. The reference brackets the second root by
    /// [1/2, 3/4]; binary64 cannot certify the slope at u=1/2 and subdivides to
    /// [5/8, 3/4], whose inside value -43.8 lies closer to the root. The
    /// returned residual scale sets the solver's start margin to 47.5, which
    /// rejects only that binary64 bracket.
    fn shallow_slope_ray() -> (EnergySurfaceRay<f64>, F<f64>) {
        let shift = 2729.0_f64 / 13.0;
        let shift = if Float::with_val(64, shift) * 13 > 2729 {
            shift
        } else {
            shift.next_up()
        };
        let ray = exact_ray(&[(1, along_x(256.0), along_x(-261.0), 12.0)], -shift);
        // The margin is epsilon times the solver tolerance 5 times the scale.
        (ray, F(47.5 / (5.0 * f64::EPSILON)))
    }

    #[test]
    fn fallback_refines_each_root_once() {
        // The second binary64 bracket fails the solver start check after the
        // first was materialized and certified.
        let (ray, e_cm) = shallow_slope_ray();
        let (fast, reference) = both_isolations(&ray);
        assert_eq!(fast.unwrap(), [(0.0, 0.25, -1), (0.625, 0.75, 1)]);
        let reference = reference.unwrap();
        assert_eq!(
            reference
                .iter()
                .map(|(left, right, crossing)| (left.to_f64(), right.to_f64(), *crossing))
                .collect_vec(),
            [(0.0, 0.25, -1), (0.5, 0.75, 1)]
        );
        let mut diagnostics = RadialRootDiagnostics::default();
        let (roots, trace) = traced_roots(&ray, &e_cm, &mut diagnostics, "fallback fixture");
        assert_eq!(trace.binary64, Binary64Outcome::Unrefinable, "{trace:?}");
        // Binary64 bounds certify both brackets in each pass.
        assert_eq!((trace.binary64_certificates, trace.certificates), (4, 0));
        let roots = roots.unwrap();
        assert_eq!(roots.len(), 2);
        assert!(roots[0].derivative_at_solution < F(0.0));
        assert!(roots[1].derivative_at_solution > F(0.0));
        for (index, (root, bracket)) in roots.iter().zip(&reference).enumerate() {
            // One recorded observation and one current-pass solve per root.
            assert_eq!(
                diagnostics.solve_counts(&RadialRootIdentity::new(format!(
                    "fallback fixture signed root {index}"
                ))),
                (1, 1),
                "{diagnostics:?}"
            );
            // The reference pass result, refined from the reference brackets.
            let (lower, upper, crossing) = materialized(bracket);
            let expected = ray
                .refine_root(
                    &lower,
                    &upper,
                    crossing,
                    0,
                    &e_cm,
                    &mut RadialRootDiagnostics::default(),
                    &RadialRootIdentity::new("reference fixture".into()),
                )
                .unwrap();
            assert_eq!(root.solution, expected.solution);
        }
    }

    #[test]
    fn binary64_brackets_stop_at_the_first_failing_cell() {
        // A later cell that fails to materialize (a point) or to certify (an
        // inverted crossing) ends the binary64 pass after the first bracket
        // was certified. `binary64_brackets` takes no diagnostics, so the
        // fallback to the reference records none.
        let (ray, e_cm) = shallow_slope_ray();
        for (cell, outcome) in [
            ((0.625, 0.625, 1), Binary64Outcome::Unmaterialized),
            ((0.625, 0.75, -1), Binary64Outcome::Uncertified),
        ] {
            let mut trace = RadialIsolationTrace::default();
            let result = ray.binary64_brackets(vec![(0.0, 0.25, -1), cell], &e_cm, &mut trace);
            assert_eq!(result.unwrap_err(), outcome);
            assert_eq!((trace.binary64_certificates, trace.certificates), (1, 0));
        }
    }

    #[test]
    fn convex_origin_tiers_agree() {
        // Both tiers of one ray, with their origin enclosure counts.
        fn tiers<T: FloatLike>(
            ray: &EnergySurfaceRay<T>,
        ) -> (Result<ConvexOrigin>, ConvexOrigin, [usize; 2]) {
            let mut evaluations = [0; 2];
            let binary64 = ray.convex_origin::<f64>(&mut evaluations[0]);
            let reference = ray.convex_origin::<Float>(&mut evaluations[1]).unwrap();
            (binary64, reference, evaluations)
        }
        let proportional = (1..=7)
            .map(|scale| {
                let scale = f64::from(scale);
                (1, along_x(3.0 * scale), ORIGIN, 4.0 * scale)
            })
            .collect_vec();
        let sunset = [
            (1, [1.0, 0.0, 0.0], ORIGIN, 0.0),
            (1, [0.0, 1.0, 0.0], ORIGIN, 0.0),
            (1, [-1.0, -1.0, 0.0], [0.0, 0.0, 0.5], 0.0),
        ];
        let receding = [(1, along_x(1.0), ORIGIN, 1.0)];
        let exterior = [
            (1, along_x(1.0), along_x(-3.0), 1.0),
            (1, along_x(1.0), along_x(-3.0), 2.0),
        ];
        let massless = [(1, along_x(1.0), ORIGIN, 0.0)];
        for (energies, shift, decided, evaluations, roots) in [
            // Interior origins, including seven proportional radicals.
            (&proportional[..], -140.0, ConvexOrigin::Interior, 1, 1),
            (&sunset[..], -2.0, ConvexOrigin::Interior, 1, 1),
            // sqrt(r^2+1)+2 recedes from its exterior origin.
            (&receding[..], 2.0, ConvexOrigin::Receding, 2, 0),
            // An exterior origin with incoming slope, and a massless energy
            // through the origin, are left to interval isolation.
            (&exterior[..], -4.0, ConvexOrigin::NotApplicable, 2, 2),
            (&massless[..], 1.0, ConvexOrigin::NotApplicable, 1, 0),
        ] {
            let ray = exact_ray::<f64>(energies, shift);
            let (binary64, reference, counts) = tiers(&ray);
            assert_eq!((binary64.unwrap(), reference), (decided, decided));
            assert_eq!(counts, [evaluations; 2]);
            // The driver takes the binary64 decision without reference work.
            let diagnostics = &mut RadialRootDiagnostics::default();
            let (found, trace) = traced_roots(&ray, &F(2.0), diagnostics, "convex tier");
            assert_eq!(trace.binary64_convex, Some(decided));
            assert_eq!(
                (trace.binary64_convex_evaluations, trace.convex_evaluations),
                (evaluations, 0)
            );
            assert_eq!(found.unwrap().len(), roots);
        }

        // A mass whose square underflows in binary64 leaves no positive
        // binary64 energy bound for the slope; the reference recedes.
        let underflowing = exact_ray::<f64>(&[(1, along_x(1.0), ORIGIN, 1.0e-200)], 1.0);
        let (binary64, reference, counts) = tiers(&underflowing);
        let error = binary64.unwrap_err().to_string();
        assert!(error.contains("energy lower bound touches zero"), "{error}");
        assert_eq!((reference, counts), (ConvexOrigin::Receding, [2, 2]));
        let diagnostics = &mut RadialRootDiagnostics::default();
        let (found, trace) = traced_roots(&underflowing, &F(2.0), diagnostics, "underflow");
        assert!(found.unwrap().is_empty());
        assert_eq!(trace.path, RadialIsolationPath::ConvexReceding);
        assert_eq!(trace.binary64_convex, Some(ConvexOrigin::Unresolved));
        assert_eq!(
            (trace.binary64_convex_evaluations, trace.convex_evaluations),
            (2, 2)
        );

        // eta(r) = r - 2^-60 in Quad: the binary64 origin value straddles
        // zero, and the reference proves the origin interior.
        let mut quad = exact_ray::<QuadFloat>(&[(1, along_x(1.0), along_x(1.0), 0.0)], -1.0);
        let one = F::<QuadFloat>::default().one();
        let tiny = one.from_i64(2).powi(-60);
        quad.shift = -(one + tiny);
        let (binary64, reference, counts) = tiers(&quad);
        assert_eq!(
            (binary64.unwrap(), reference, counts),
            (ConvexOrigin::Unresolved, ConvexOrigin::Interior, [1, 1])
        );
        let diagnostics = &mut RadialRootDiagnostics::default();
        let (found, trace) = traced_roots(&quad, &one.from_i64(2), diagnostics, "Quad interior");
        assert_eq!(trace.path, RadialIsolationPath::ConvexInterior);
        assert_eq!(trace.binary64_convex, Some(ConvexOrigin::Unresolved));
        assert_eq!(
            (trace.binary64_convex_evaluations, trace.convex_evaluations),
            (1, 1)
        );
        let found = found.unwrap();
        assert_eq!(found.len(), 1);
        assert!((found[0].solution / tiny - one).abs() < one.from_i64(2).powi(-20));
    }
}
