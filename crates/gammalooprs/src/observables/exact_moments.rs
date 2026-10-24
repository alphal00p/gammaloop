//! Exact additive statistics of final, finite, represented binary64 samples.
//!
//! This does not change the integrand's precision or weighting boundary. Callers
//! must first finish each correlated sample, including cut zeros, and pass its
//! final represented value exactly once. These helpers do not by themselves
//! enable an exact consumer mode: every enclosing count, bin and capability must
//! also be checked. Legacy rounded moments cannot initialize this state.

use serde::{Deserialize, Deserializer, Serialize, Serializer};
use symbolica::domains::{
    float::{Float, Real, RealLike, RoundingDirection},
    integer::Integer,
    rational::Rational,
};
use thiserror::Error;

pub const MAX_SAMPLE_COUNT: u64 = i64::MAX as u64;
pub const MAX_CANONICAL_OBJECT_BYTES: usize = 12_288;
pub const MAX_EXACT_OBJECTS: usize = 16_384;
pub const MAX_JSON_FRAME_BYTES: usize = 16 * 1024 * 1024;
const PRECISIONS: [u32; 8] = [128, 256, 512, 1024, 2048, 4096, 8192, 16_384];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MomentOrder {
    Second,
    Fourth,
}

#[derive(Debug, Clone, Error, PartialEq, Eq)]
pub enum ExactMomentError {
    #[error("sample must be finite binary64")]
    NonFiniteSample,
    #[error("sample count exceeds the v1 i64 range")]
    CountOverflow,
    #[error("exact moment orders differ")]
    OrderMismatch,
    #[error("invalid exact moment payload: {0}")]
    InvalidPayload(String),
    #[error("exact moment resource limit exceeded: {0}")]
    ResourceLimit(&'static str),
    #[error("this estimator requires at least {0} samples")]
    InsufficientSamples(u64),
    #[error("negative square-root argument")]
    NegativeSquareRoot,
    #[error("nonzero result rounds below the binary64 reporting range")]
    ReportingUnderflow,
    #[error("result exceeds the finite binary64 reporting range")]
    ReportingOverflow,
    #[error("unable to certify final binary64 rounding within the precision cap")]
    RoundingNotCertified,
    #[error("unsupported exact represented-sample transformation: {0}")]
    Unsupported(&'static str),
}

type Result<T> = std::result::Result<T, ExactMomentError>;

/// Zero rounding is permitted only for non-authoritative raw diagnostics in a
/// mandatory exact mode. Primary estimates and uncertainties use `Strict`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportPolicy {
    Strict,
    ExactModeRawDiagnostic,
}

/// Fields are private so every state is constructed from actual samples or a
/// checked canonical payload. Rational's backend-specific serde is never used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactF64Moments {
    sample_count: u64,
    sum: Rational,
    sum_abs: Rational,
    sum_squared: Rational,
    higher: Option<(Rational, Rational)>,
}

impl ExactF64Moments {
    pub fn empty(order: MomentOrder) -> Self {
        Self {
            sample_count: 0,
            sum: Rational::zero(),
            sum_abs: Rational::zero(),
            sum_squared: Rational::zero(),
            higher: (order == MomentOrder::Fourth).then(|| (Rational::zero(), Rational::zero())),
        }
    }

    pub fn sample_count(&self) -> u64 {
        self.sample_count
    }

    pub fn order(&self) -> MomentOrder {
        if self.higher.is_some() {
            MomentOrder::Fourth
        } else {
            MomentOrder::Second
        }
    }

    pub fn sum(&self) -> &Rational {
        &self.sum
    }
    pub fn sum_abs(&self) -> &Rational {
        &self.sum_abs
    }
    pub fn sum_squared(&self) -> &Rational {
        &self.sum_squared
    }
    pub fn higher(&self) -> Option<(&Rational, &Rational)> {
        self.higher.as_ref().map(|(third, fourth)| (third, fourth))
    }

    /// Transactional append of one fully weighted, grouped sample, including 0.
    pub fn push(&mut self, sample: f64) -> Result<()> {
        let x = Rational::try_from(sample).map_err(|_| ExactMomentError::NonFiniteSample)?;
        let mut next = self.clone();
        next.sample_count = checked_count(self.sample_count, 1)?;
        next.sum += &x;
        next.sum_abs += x.abs();
        let square = &x * &x;
        next.sum_squared += &square;
        if let Some((third, fourth)) = &mut next.higher {
            *third += &square * &x;
            *fourth += &square * &square;
        }
        // Actual finite samples preserve the moment inequalities by construction.
        // Count bounds imply the fixed dyadic resource bounds, without a GCD or
        // decimal serialization on every sample.
        *self = next;
        Ok(())
    }

    /// Missing/cut bins contribute zero for each additional independent sample.
    pub fn extend_zeros(&mut self, count: u64) -> Result<()> {
        self.sample_count = checked_count(self.sample_count, count)?;
        Ok(())
    }

    pub fn merge(&mut self, other: &Self) -> Result<()> {
        if self.order() != other.order() {
            return Err(ExactMomentError::OrderMismatch);
        }
        let mut next = self.clone();
        next.sample_count = checked_count(self.sample_count, other.sample_count)?;
        next.sum += &other.sum;
        next.sum_abs += &other.sum_abs;
        next.sum_squared += &other.sum_squared;
        if let (Some((a3, a4)), Some((b3, b4))) = (&mut next.higher, &other.higher) {
            *a3 += b3;
            *a4 += b4;
        }
        next.validate()?;
        *self = next;
        Ok(())
    }

    /// Exact multiplication by an arbitrary f64 does not imply that each
    /// product is a represented f64; re-rounding cannot be recovered from sums.
    pub fn rescale(&mut self, _factor: f64) -> Result<()> {
        Err(ExactMomentError::Unsupported(
            "rescale requires a tagged scaled-basis domain or the original samples",
        ))
    }

    /// A spatial rebin needs cross-bin products, which marginal moments omit.
    pub fn merge_correlated_bin(&mut self, _other: &Self) -> Result<()> {
        Err(ExactMomentError::Unsupported(
            "spatial rebin requires cross-bin moments",
        ))
    }

    pub fn validate_count(&self, enclosing_count: u64) -> Result<()> {
        if self.sample_count != enclosing_count {
            return Err(invalid("sample_count differs from enclosing accumulator"));
        }
        Ok(())
    }

    /// Exact variance numerator N*S2-S1^2. No rounded mean is involved.
    pub fn variance_numerator(&self) -> Rational {
        &Rational::from(self.sample_count) * &self.sum_squared - &self.sum * &self.sum
    }

    pub fn absolute_variance_numerator(&self) -> Rational {
        &Rational::from(self.sample_count) * &self.sum_squared - &self.sum_abs * &self.sum_abs
    }

    pub fn mean(&self) -> Result<f64> {
        self.require_samples(1)?;
        certified_round_f64(
            &(&self.sum / &Rational::from(self.sample_count)),
            ReportPolicy::Strict,
        )
    }

    pub fn mean_abs(&self) -> Result<f64> {
        self.require_samples(1)?;
        certified_round_f64(
            &(&self.sum_abs / &Rational::from(self.sample_count)),
            ReportPolicy::Strict,
        )
    }

    pub fn population_variance(&self) -> Result<f64> {
        self.require_samples(1)?;
        certified_round_f64(
            &(self.variance_numerator() / Rational::from(self.sample_count).pow(2)),
            ReportPolicy::Strict,
        )
    }

    /// GammaBoard's population convention: SEM^2=(N*S2-S1^2)/N^3.
    pub fn population_sem(&self) -> Result<f64> {
        self.require_samples(2)?;
        certified_sqrt_f64(&(self.variance_numerator() / Rational::from(self.sample_count).pow(3)))
    }

    pub fn population_absolute_sem(&self) -> Result<f64> {
        self.require_samples(2)?;
        certified_sqrt_f64(
            &(self.absolute_variance_numerator() / Rational::from(self.sample_count).pow(3)),
        )
    }

    /// Histogram convention: SEM^2=(N*S2-S1^2)/(N^2*(N-1)).
    pub fn histogram_sem(&self) -> Result<f64> {
        self.require_samples(2)?;
        let denominator =
            Rational::from(self.sample_count).pow(2) * Rational::from(self.sample_count - 1);
        certified_sqrt_f64(&(self.variance_numerator() / denominator))
    }

    fn require_samples(&self, minimum: u64) -> Result<()> {
        if self.sample_count < minimum {
            Err(ExactMomentError::InsufficientSamples(minimum))
        } else {
            Ok(())
        }
    }

    /// Parse directly from original JSON bytes: a preceding `serde_json::Value`
    /// conversion would erase duplicate-key evidence. The outer protocol must
    /// separately enforce its aggregate frame/object limits.
    pub fn from_json(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_JSON_FRAME_BYTES {
            return Err(ExactMomentError::ResourceLimit("JSON frame"));
        }
        serde_json::from_slice(bytes).map_err(|error| invalid(error.to_string()))
    }

    pub fn to_canonical_json(&self) -> Result<Vec<u8>> {
        let bytes = serde_json::to_vec(self).map_err(|error| invalid(error.to_string()))?;
        if bytes.len() > MAX_CANONICAL_OBJECT_BYTES {
            return Err(ExactMomentError::ResourceLimit("canonical object"));
        }
        Ok(bytes)
    }

    fn validate(&self) -> Result<()> {
        if self.sample_count > MAX_SAMPLE_COUNT {
            return Err(ExactMomentError::CountOverflow);
        }
        let n = Rational::from(self.sample_count);
        let max = Rational::try_from(f64::MAX).expect("finite constant");
        for (value, order) in [(&self.sum, 1), (&self.sum_abs, 1), (&self.sum_squared, 2)]
            .into_iter()
            .chain(
                self.higher
                    .iter()
                    .flat_map(|(third, fourth)| [(third, 3), (fourth, 4)]),
            )
        {
            check_rational_domain(value, order)?;
            if value.abs() > &n * &max.pow(order as u64) {
                return Err(invalid(
                    "moment exceeds finite represented-input magnitude bound",
                ));
            }
        }
        if self.sample_count == 0 {
            // Magnitude bounds above force every sum to zero.
            return Ok(());
        }
        if self.sum_squared.is_negative() || self.sum_abs < self.sum.abs() {
            return Err(invalid(
                "negative second moment or inconsistent absolute sum",
            ));
        }
        let abs_square = self.sum_abs.pow(2);
        if self.sum_squared > abs_square || abs_square > &n * &self.sum_squared {
            return Err(invalid("absolute-sum/second-moment inequalities violated"));
        }
        if self.sample_count == 1
            && (self.sum_abs != self.sum.abs() || self.sum_squared != self.sum.pow(2))
        {
            return Err(invalid("one-sample moments are inconsistent"));
        }
        if let Some((third, fourth)) = &self.higher {
            if fourth.is_negative()
                || &n * fourth < self.sum_squared.pow(2)
                || &self.sum_squared * fourth < third.pow(2)
            {
                return Err(invalid("negative Hankel principal minor"));
            }
            let determinant = &n * &(&self.sum_squared * fourth - third.pow(2))
                - &self.sum * &(&self.sum * fourth - &self.sum_squared * third)
                + &self.sum_squared * &(&self.sum * third - self.sum_squared.pow(2));
            if determinant.is_negative() {
                return Err(invalid("negative Hankel determinant"));
            }
            if self.sample_count == 1 && (*third != self.sum.pow(3) || *fourth != self.sum.pow(4)) {
                return Err(invalid("one-sample higher moments are inconsistent"));
            }
        }
        Ok(())
    }
}

fn checked_count(left: u64, right: u64) -> Result<u64> {
    left.checked_add(right)
        .filter(|count| *count <= MAX_SAMPLE_COUNT)
        .ok_or(ExactMomentError::CountOverflow)
}

fn invalid(message: impl Into<String>) -> ExactMomentError {
    ExactMomentError::InvalidPayload(message.into())
}

// A finite x=k*2^-1074 has |k|<2^2098. At N<2^63, the order-r
// fixed-grid numerator needs at most 2098*r+63 bits. Reduction only decreases it.
fn check_rational_domain(value: &Rational, order: usize) -> Result<()> {
    if value.numerator_ref().significant_bits() > (2098 * order + 63) as u64 {
        return Err(ExactMomentError::ResourceLimit("moment numerator bits"));
    }
    let denominator = value.denominator_ref();
    let bits = denominator.significant_bits();
    if bits == 0
        || bits > (1074 * order + 1) as u64
        || *denominator != Integer::from(2).pow(bits - 1)
    {
        return Err(invalid(
            "moment denominator is not a bounded positive power of two",
        ));
    }
    Ok(())
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct RationalWire {
    numerator: String,
    denominator: String,
}

impl From<&Rational> for RationalWire {
    fn from(value: &Rational) -> Self {
        Self {
            numerator: value.numerator_ref().to_string(),
            denominator: value.denominator_ref().to_string(),
        }
    }
}

impl RationalWire {
    fn check_lexical(&self, order: usize) -> Result<()> {
        const NUM_DIGITS: [usize; 4] = [651, 1283, 1914, 2546];
        const DEN_DIGITS: [usize; 4] = [324, 647, 971, 1294];
        let digits = self.numerator.strip_prefix('-').unwrap_or(&self.numerator);
        if digits.len() > NUM_DIGITS[order - 1] || self.denominator.len() > DEN_DIGITS[order - 1] {
            return Err(ExactMomentError::ResourceLimit("moment decimal digits"));
        }
        if !canonical_unsigned(digits)
            || self.numerator == "-0"
            || !canonical_unsigned(&self.denominator)
            || self.denominator == "0"
        {
            return Err(invalid("noncanonical integer spelling"));
        }
        Ok(())
    }

    fn decode(self, order: usize) -> Result<Rational> {
        self.check_lexical(order)?;
        let numerator: Integer = self.numerator.parse().map_err(invalid)?;
        let denominator: Integer = self.denominator.parse().map_err(invalid)?;
        // Lexical checking has already established a strictly positive denominator.
        let value = Rational::new(numerator.clone(), denominator.clone());
        if value.numerator_ref() != &numerator || value.denominator_ref() != &denominator {
            return Err(invalid("rational is not reduced or zero is not 0/1"));
        }
        check_rational_domain(&value, order)?;
        Ok(value)
    }
}

fn canonical_unsigned(value: &str) -> bool {
    value == "0"
        || (!value.is_empty()
            && value.as_bytes()[0].is_ascii_digit()
            && value.as_bytes()[0] != b'0'
            && value.bytes().all(|byte| byte.is_ascii_digit()))
}

fn present_rational<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> std::result::Result<Option<RationalWire>, D::Error> {
    // Unlike Option<T>'s default decoder, present null is invalid.
    RationalWire::deserialize(deserializer).map(Some)
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MomentsWire {
    version: u8,
    sample_count: String,
    max_order: u8,
    sum: RationalWire,
    sum_abs: RationalWire,
    sum_squared: RationalWire,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "present_rational"
    )]
    sum_cubed: Option<RationalWire>,
    #[serde(
        default,
        skip_serializing_if = "Option::is_none",
        deserialize_with = "present_rational"
    )]
    sum_fourth: Option<RationalWire>,
}

impl From<&ExactF64Moments> for MomentsWire {
    fn from(value: &ExactF64Moments) -> Self {
        Self {
            version: 1,
            sample_count: value.sample_count.to_string(),
            max_order: if value.higher.is_some() { 4 } else { 2 },
            sum: (&value.sum).into(),
            sum_abs: (&value.sum_abs).into(),
            sum_squared: (&value.sum_squared).into(),
            sum_cubed: value.higher.as_ref().map(|(third, _)| third.into()),
            sum_fourth: value.higher.as_ref().map(|(_, fourth)| fourth.into()),
        }
    }
}

impl TryFrom<MomentsWire> for ExactF64Moments {
    type Error = ExactMomentError;
    fn try_from(wire: MomentsWire) -> Result<Self> {
        if wire.version != 1 || !matches!(wire.max_order, 2 | 4) {
            return Err(invalid("unsupported version or moment order"));
        }
        if wire.sample_count.len() > 19 || !canonical_unsigned(&wire.sample_count) {
            return Err(invalid("noncanonical or oversized sample_count"));
        }
        let sample_count = wire
            .sample_count
            .parse::<u64>()
            .map_err(|_| ExactMomentError::CountOverflow)?;
        if sample_count > MAX_SAMPLE_COUNT {
            return Err(ExactMomentError::CountOverflow);
        }
        if (wire.max_order == 4) != (wire.sum_cubed.is_some() && wire.sum_fourth.is_some())
            || (wire.max_order == 2 && (wire.sum_cubed.is_some() || wire.sum_fourth.is_some()))
        {
            return Err(invalid("higher moments do not match max_order"));
        }
        // Check ALL strings and total canonical bytes before any big-integer parse.
        for (value, order) in [(&wire.sum, 1), (&wire.sum_abs, 1), (&wire.sum_squared, 2)]
            .into_iter()
            .chain(wire.sum_cubed.iter().map(|value| (value, 3)))
            .chain(wire.sum_fourth.iter().map(|value| (value, 4)))
        {
            value.check_lexical(order)?;
        }
        if serde_json::to_vec(&wire)
            .map_err(|error| invalid(error.to_string()))?
            .len()
            > MAX_CANONICAL_OBJECT_BYTES
        {
            return Err(ExactMomentError::ResourceLimit("canonical object"));
        }
        let higher = match (wire.sum_cubed, wire.sum_fourth) {
            (Some(third), Some(fourth)) => Some((third.decode(3)?, fourth.decode(4)?)),
            (None, None) => None,
            _ => return Err(invalid("incomplete higher moments")),
        };
        let value = Self {
            sample_count,
            sum: wire.sum.decode(1)?,
            sum_abs: wire.sum_abs.decode(1)?,
            sum_squared: wire.sum_squared.decode(2)?,
            higher,
        };
        value.validate()?;
        Ok(value)
    }
}

impl Serialize for ExactF64Moments {
    fn serialize<S: Serializer>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error> {
        MomentsWire::from(self).serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for ExactF64Moments {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> std::result::Result<Self, D::Error> {
        MomentsWire::deserialize(deserializer)?
            .try_into()
            .map_err(serde::de::Error::custom)
    }
}

/// Correctly round a general rational, not only a dyadic stored moment. Final
/// means and variance denominators contain factors of N and can be non-dyadic.
pub fn certified_round_f64(value: &Rational, policy: ReportPolicy) -> Result<f64> {
    let magnitude = certified_positive(value.abs(), false, policy)?;
    // Unique exact zero is +0; a negative diagnostic underflow retains -0.
    Ok(if value.is_negative() {
        -magnitude
    } else {
        magnitude
    })
}

/// Take the square root in arbitrary precision before binary64 conversion,
/// then certify its final rounding cell with exact squared midpoint bounds.
pub fn certified_sqrt_f64(value: &Rational) -> Result<f64> {
    if value.is_negative() {
        return Err(ExactMomentError::NegativeSquareRoot);
    }
    certified_positive(value.clone(), true, ReportPolicy::Strict)
}

fn certified_positive(value: Rational, square_root: bool, policy: ReportPolicy) -> Result<f64> {
    if value.is_zero() {
        return Ok(0.0);
    }
    let threshold = |boundary: Rational| {
        if square_root {
            boundary.pow(2)
        } else {
            boundary
        }
    };
    let half_min =
        Rational::try_from(f64::from_bits(1)).expect("finite constant") / Rational::from(2);
    if value <= threshold(half_min) {
        return if policy == ReportPolicy::ExactModeRawDiagnostic {
            Ok(0.0)
        } else {
            Err(ExactMomentError::ReportingUnderflow)
        };
    }
    let max = Rational::try_from(f64::MAX).expect("finite constant");
    let conceptual_next = Rational::from(Integer::from(2).pow(1024));
    // At this midpoint the even conceptual next value overflows; MAX is odd.
    if value >= threshold((&max + &conceptual_next) / Rational::from(2)) {
        return Err(ExactMomentError::ReportingOverflow);
    }
    for precision in PRECISIONS {
        let candidate = Float::from_rational_round(&value, precision, RoundingDirection::Nearest);
        let candidate = if square_root {
            candidate.sqrt()
        } else {
            candidate
        }
        .to_f64();
        let candidate = if candidate.is_infinite() {
            f64::MAX
        } else if candidate == 0.0 {
            f64::from_bits(1)
        } else {
            candidate
        };
        // A rounded AP approximation may land on the other side of a binary64
        // tie. Certify adjacent cells as well; accepting one is still exact.
        for candidate in [candidate, candidate.next_down(), candidate.next_up()] {
            if candidate > 0.0
                && candidate.is_finite()
                && positive_cell_contains(&value, candidate, square_root)
            {
                return Ok(candidate);
            }
        }
    }
    Err(ExactMomentError::RoundingNotCertified)
}

fn positive_cell_contains(value: &Rational, candidate: f64, square_root: bool) -> bool {
    let center = Rational::try_from(candidate).expect("positive finite candidate");
    let previous =
        Rational::try_from(candidate.next_down()).expect("nonnegative finite predecessor");
    let next = if candidate == f64::MAX {
        Rational::from(Integer::from(2).pow(1024))
    } else {
        Rational::try_from(candidate.next_up()).expect("finite successor")
    };
    let mut lower = (&previous + &center) / Rational::from(2);
    let mut upper = (&center + &next) / Rational::from(2);
    if square_root {
        lower = lower.pow(2);
        upper = upper.pow(2);
    }
    if candidate.to_bits() & 1 == 0 {
        *value >= lower && *value <= upper
    } else {
        *value > lower && *value < upper
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{Value, json};

    fn moments(samples: &[f64], order: MomentOrder) -> ExactF64Moments {
        let mut result = ExactF64Moments::empty(order);
        for &sample in samples {
            result.push(sample).unwrap();
        }
        result
    }

    fn rational(value: f64) -> Rational {
        Rational::try_from(value).unwrap()
    }

    fn decode(value: Value) -> Result<ExactF64Moments> {
        ExactF64Moments::from_json(&serde_json::to_vec(&value).unwrap())
    }

    fn wire(value: &ExactF64Moments) -> Value {
        serde_json::to_value(value).unwrap()
    }

    #[test]
    fn high_offset_pair_preserves_exact_variance_and_named_sem_conventions() {
        let sample = moments(&[134_217_728.0, 134_217_729.0], MomentOrder::Second);
        assert_eq!(sample.sum(), &Rational::from(268_435_457));
        assert_eq!(sample.variance_numerator(), Rational::one());
        assert_eq!(sample.histogram_sem().unwrap(), 0.5);
        assert_eq!(sample.population_sem().unwrap(), (0.125_f64).sqrt());
        assert_eq!(sample.population_variance().unwrap(), 0.25);
        assert_eq!(sample.mean().unwrap(), 134_217_728.5);
        let expected = br#"{"version":1,"sample_count":"2","max_order":2,"sum":{"numerator":"268435457","denominator":"1"},"sum_abs":{"numerator":"268435457","denominator":"1"},"sum_squared":{"numerator":"36028797287399425","denominator":"1"}}"#;
        assert_eq!(sample.to_canonical_json().unwrap(), expected);
        assert_eq!(ExactF64Moments::from_json(expected).unwrap(), sample);
    }

    #[test]
    fn adjacent_midpoints_and_three_samples_are_partition_invariant() {
        let a = 1.0_f64;
        let b = a.next_up();
        let c = b.next_up();
        for values in [[a, b], [b, c]] {
            let full = moments(&values, MomentOrder::Fourth);
            assert_eq!(full.histogram_sem().unwrap(), (values[1] - values[0]) / 2.0);
            let mut partition = moments(&values[..1], MomentOrder::Fourth);
            partition
                .merge(&moments(&values[1..], MomentOrder::Fourth))
                .unwrap();
            assert_eq!(full, partition);
        }
        let full = moments(&[a, b, c], MomentOrder::Fourth);
        let h = rational(b - a);
        assert_eq!(full.variance_numerator(), Rational::from(6) * h.pow(2));
        let mut partition = moments(&[a, c], MomentOrder::Fourth);
        partition
            .merge(&moments(&[b], MomentOrder::Fourth))
            .unwrap();
        assert_eq!(full, partition);
        assert_eq!(
            full.histogram_sem().unwrap(),
            certified_sqrt_f64(&(h.pow(2) / Rational::from(3))).unwrap()
        );
    }

    #[test]
    fn all_partition_boundaries_and_reversed_order_have_identical_exact_sums() {
        let samples = [
            f64::from_bits(1),
            1.0,
            -1.0,
            1.0_f64.next_up(),
            -12.5,
            12.5,
            0.0,
            1e100,
            -1e100,
        ];
        for order in [MomentOrder::Second, MomentOrder::Fourth] {
            let full = moments(&samples, order);
            for split in 0..=samples.len() {
                let mut partition = moments(&samples[..split], order);
                partition.merge(&moments(&samples[split..], order)).unwrap();
                assert_eq!(partition, full);
                let bytes = partition.to_canonical_json().unwrap();
                assert_eq!(ExactF64Moments::from_json(&bytes).unwrap(), full);
            }
            let reverse: Vec<_> = samples.iter().copied().rev().collect();
            assert_eq!(moments(&reverse, order), full);
        }
    }

    #[test]
    fn absolute_sum_is_of_final_grouped_samples_and_cut_zeros_count() {
        let mut signed = moments(&[-3.0, 3.0], MomentOrder::Fourth);
        assert_eq!(signed.sum(), &Rational::zero());
        assert_eq!(signed.sum_abs(), &Rational::from(6));
        assert_eq!(signed.population_absolute_sem().unwrap(), 0.0);
        signed.extend_zeros(2).unwrap();
        assert_eq!(signed, moments(&[-3.0, 3.0, 0.0, 0.0], MomentOrder::Fourth));
        assert_eq!(signed.sample_count(), 4);
        assert_eq!(signed.mean_abs().unwrap(), 1.5);
        // The complete +3/-3 counterevent group supplies one final sample zero.
        let cancelled_group = moments(&[0.0], MomentOrder::Second);
        assert_eq!(cancelled_group.sum_abs(), &Rational::zero());
        assert_eq!(cancelled_group.sample_count(), 1);
    }

    #[test]
    fn empty_single_and_constant_states_do_not_claim_precision() {
        let empty = ExactF64Moments::empty(MomentOrder::Fourth);
        assert_eq!(empty.mean(), Err(ExactMomentError::InsufficientSamples(1)));
        assert_eq!(
            empty.histogram_sem(),
            Err(ExactMomentError::InsufficientSamples(2))
        );
        let single = moments(&[7.0], MomentOrder::Fourth);
        assert_eq!(
            single.population_sem(),
            Err(ExactMomentError::InsufficientSamples(2))
        );
        assert_eq!(single.population_variance().unwrap(), 0.0);
        let mut constant = moments(&[7.0; 9], MomentOrder::Fourth);
        assert_eq!(constant.mean().unwrap(), 7.0);
        assert_eq!(constant.histogram_sem().unwrap(), 0.0);
        let before = constant.clone();
        constant.merge(&empty).unwrap();
        assert_eq!(constant, before);
        let zeros = moments(&[0.0, -0.0], MomentOrder::Fourth);
        assert_eq!(zeros.mean().unwrap().to_bits(), 0);
        assert_eq!(zeros.histogram_sem().unwrap().to_bits(), 0);
    }

    #[test]
    fn exact_subnormal_square_survives_without_a_binary64_variance() {
        let tiny = f64::from_bits(1);
        let sample = moments(&[tiny, 3.0 * tiny], MomentOrder::Fourth);
        assert!(!sample.sum_squared().is_zero());
        assert_eq!(
            sample.population_variance(),
            Err(ExactMomentError::ReportingUnderflow)
        );
        assert_eq!(sample.histogram_sem().unwrap(), tiny);
        assert_eq!(sample.population_sem().unwrap(), tiny);
        assert_eq!(sample.mean().unwrap(), 2.0 * tiny);
        assert_eq!(
            certified_round_f64(sample.sum_squared(), ReportPolicy::ExactModeRawDiagnostic)
                .unwrap(),
            0.0
        );
        assert_eq!(
            ExactF64Moments::from_json(&sample.to_canonical_json().unwrap()).unwrap(),
            sample
        );
    }

    #[test]
    fn final_rounding_is_exact_at_signed_midpoints_and_odd_denominators() {
        let a = 1.0_f64;
        let b = a.next_up();
        let c = b.next_up();
        for (left, right, expected) in [(a, b, a), (b, c, c)] {
            let midpoint = (rational(left) + rational(right)) / Rational::from(2);
            assert_eq!(
                certified_round_f64(&midpoint, ReportPolicy::Strict).unwrap(),
                expected
            );
            assert_eq!(
                certified_round_f64(&midpoint.neg(), ReportPolicy::Strict).unwrap(),
                -expected
            );
            assert_eq!(certified_sqrt_f64(&midpoint.pow(2)).unwrap(), expected);
        }
        assert_eq!(
            certified_round_f64(&Rational::new(1, 3), ReportPolicy::Strict).unwrap(),
            1.0 / 3.0
        );
        assert_eq!(
            certified_sqrt_f64(&Rational::new(1, 3)).unwrap(),
            (1.0_f64 / 3.0).sqrt()
        );
        assert_eq!(
            certified_sqrt_f64(&Rational::from(-1)),
            Err(ExactMomentError::NegativeSquareRoot)
        );
    }

    #[test]
    fn zero_cell_is_checked_exactly_and_diagnostic_underflow_is_explicit() {
        let half_min = rational(f64::from_bits(1)) / Rational::from(2);
        let perturbation = Rational::new(Integer::one(), Integer::from(2).pow(2000));
        let below = &half_min - &perturbation;
        let above = &half_min + &perturbation;
        for value in [&below, &half_min] {
            assert_eq!(
                certified_round_f64(value, ReportPolicy::Strict),
                Err(ExactMomentError::ReportingUnderflow)
            );
            assert_eq!(
                certified_round_f64(value, ReportPolicy::ExactModeRawDiagnostic)
                    .unwrap()
                    .to_bits(),
                0
            );
            assert_eq!(
                certified_round_f64(&value.neg(), ReportPolicy::ExactModeRawDiagnostic)
                    .unwrap()
                    .to_bits(),
                (-0.0_f64).to_bits()
            );
            assert_eq!(
                certified_sqrt_f64(&value.pow(2)),
                Err(ExactMomentError::ReportingUnderflow)
            );
        }
        assert_eq!(
            certified_round_f64(&above, ReportPolicy::Strict).unwrap(),
            f64::from_bits(1)
        );
        assert_eq!(
            certified_sqrt_f64(&above.pow(2)).unwrap(),
            f64::from_bits(1)
        );
        for value in [0.0, -0.0] {
            assert_eq!(
                certified_round_f64(&rational(value), ReportPolicy::Strict)
                    .unwrap()
                    .to_bits(),
                0
            );
        }
    }

    #[test]
    fn finite_overflow_boundary_and_tie_are_not_saturated() {
        let max = rational(f64::MAX);
        let boundary = (&max + &Rational::from(Integer::from(2).pow(1024))) / Rational::from(2);
        let just_below = &boundary - &Rational::one();
        assert_eq!(
            certified_round_f64(&max, ReportPolicy::Strict).unwrap(),
            f64::MAX
        );
        assert_eq!(
            certified_round_f64(&just_below, ReportPolicy::Strict).unwrap(),
            f64::MAX
        );
        assert_eq!(certified_sqrt_f64(&just_below.pow(2)).unwrap(), f64::MAX);
        for value in [&boundary, &(&boundary + &Rational::one())] {
            assert_eq!(
                certified_round_f64(value, ReportPolicy::Strict),
                Err(ExactMomentError::ReportingOverflow)
            );
            assert_eq!(
                certified_round_f64(&value.neg(), ReportPolicy::Strict),
                Err(ExactMomentError::ReportingOverflow)
            );
            assert_eq!(
                certified_sqrt_f64(&value.pow(2)),
                Err(ExactMomentError::ReportingOverflow)
            );
        }
    }

    #[test]
    fn wire_rejects_duplicates_unknowns_nulls_and_wrong_order_fields() {
        let sample = moments(&[1.0, -1.0], MomentOrder::Second);
        let original = String::from_utf8(sample.to_canonical_json().unwrap()).unwrap();
        for duplicate in [
            original.replacen("\"version\":1", "\"version\":1,\"version\":1", 1),
            original.replacen(
                "\"numerator\":\"0\"",
                "\"numerator\":\"0\",\"numerator\":\"0\"",
                1,
            ),
        ] {
            assert!(ExactF64Moments::from_json(duplicate.as_bytes()).is_err());
        }
        for key in ["sum_cubed", "sum_fourth", "unrecognized"] {
            let mut bad = wire(&sample);
            bad[key] = Value::Null;
            assert!(decode(bad).is_err(), "accepted {key}");
        }
        for wrong in [json!(1.0), json!("1"), json!(null)] {
            let mut bad = wire(&sample);
            bad["version"] = wrong;
            assert!(decode(bad).is_err());
        }
        let fourth = moments(&[1.0], MomentOrder::Fourth);
        for key in ["sum_cubed", "sum_fourth"] {
            let mut bad = wire(&fourth);
            bad.as_object_mut().unwrap().remove(key);
            assert!(decode(bad).is_err());
        }
    }

    #[test]
    fn wire_checks_canonical_grammar_sizes_and_dyadic_domain() {
        let sample = moments(&[1.0], MomentOrder::Second);
        for text in ["", "+1", "-0", "01", " 1", "1 ", "1_0", "1.0", "1e0", "١"] {
            let mut bad = wire(&sample);
            bad["sum"]["numerator"] = json!(text);
            assert!(decode(bad).is_err(), "accepted {text:?}");
        }
        for text in ["0", "-1", "+1", "01", "3", "0.5"] {
            let mut bad = wire(&sample);
            bad["sum"]["denominator"] = json!(text);
            assert!(decode(bad).is_err(), "accepted denominator {text:?}");
        }
        let mut unreduced = wire(&sample);
        unreduced["sum"] = json!({"numerator":"2","denominator":"2"});
        assert!(decode(unreduced).is_err());
        let mut zero = wire(&moments(&[0.0], MomentOrder::Second));
        zero["sum"]["denominator"] = json!("2");
        assert!(decode(zero).is_err());
        let mut huge = wire(&sample);
        huge["sum"]["numerator"] = json!("1".repeat(652));
        assert!(decode(huge).is_err());
        let mut exponent = wire(&sample);
        exponent["sum"]["denominator"] = json!(Integer::from(2).pow(1075).to_string());
        assert!(decode(exponent).is_err());
        for text in [
            "-1",
            "+1",
            "01",
            "9223372036854775808",
            "99999999999999999999",
        ] {
            let mut bad = wire(&sample);
            bad["sample_count"] = json!(text);
            assert!(decode(bad).is_err());
        }
    }

    #[test]
    fn wire_rejects_inconsistent_empty_one_sample_absolute_and_hankel_moments() {
        let sample = moments(&[1.0], MomentOrder::Fourth);
        for (key, text) in [
            ("sum_abs", "0"),
            ("sum_squared", "2"),
            ("sum_cubed", "2"),
            ("sum_fourth", "2"),
        ] {
            let mut bad = wire(&sample);
            bad[key]["numerator"] = json!(text);
            assert!(decode(bad).is_err(), "accepted inconsistent {key}");
        }
        let mut empty = wire(&sample);
        empty["sample_count"] = json!("0");
        assert!(decode(empty).is_err());
        let signed = moments(&[-1.0, 1.0], MomentOrder::Fourth);
        let mut negative_minor = wire(&signed);
        negative_minor["sum_cubed"]["numerator"] = json!("3");
        assert!(decode(negative_minor).is_err());
        let mut impossible_abs = wire(&signed);
        impossible_abs["sum_abs"]["numerator"] = json!("1");
        assert!(decode(impossible_abs).is_err());
        let mut negative_square = wire(&signed);
        negative_square["sum_squared"]["numerator"] = json!("-2");
        assert!(decode(negative_square).is_err());
    }

    #[test]
    fn maximum_finite_samples_and_maximum_zero_count_round_trip_within_caps() {
        let sample = moments(
            &[f64::MAX, -f64::MAX, f64::from_bits(1)],
            MomentOrder::Fourth,
        );
        let bytes = sample.to_canonical_json().unwrap();
        assert!(bytes.len() <= MAX_CANONICAL_OBJECT_BYTES);
        assert_eq!(ExactF64Moments::from_json(&bytes).unwrap(), sample);
        // Exact storage is valid even when the legacy required raw-f64 fields
        // would reject the payload; no finite-raw schema claim is made here.
        assert_eq!(
            certified_round_f64(sample.sum_squared(), ReportPolicy::ExactModeRawDiagnostic),
            Err(ExactMomentError::ReportingOverflow)
        );
        let mut zeros = ExactF64Moments::empty(MomentOrder::Fourth);
        zeros.extend_zeros(MAX_SAMPLE_COUNT).unwrap();
        assert_eq!(
            ExactF64Moments::from_json(&zeros.to_canonical_json().unwrap()).unwrap(),
            zeros
        );
    }

    #[test]
    fn failures_are_transactional_and_unsupported_transforms_do_not_drop_authority() {
        let mut sample = moments(&[1.0], MomentOrder::Second);
        let before = sample.clone();
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(sample.push(value), Err(ExactMomentError::NonFiniteSample));
            assert_eq!(sample, before);
        }
        assert_eq!(
            sample.merge(&moments(&[2.0], MomentOrder::Fourth)),
            Err(ExactMomentError::OrderMismatch)
        );
        assert_eq!(sample, before);
        for factor in [0.0, 1.0, 2.0, f64::NAN] {
            assert!(matches!(
                sample.rescale(factor),
                Err(ExactMomentError::Unsupported(_))
            ));
        }
        assert!(matches!(
            sample.merge_correlated_bin(&before),
            Err(ExactMomentError::Unsupported(_))
        ));
        assert_eq!(sample, before);
        assert!(sample.validate_count(2).is_err());
        let mut full = ExactF64Moments::empty(MomentOrder::Second);
        full.extend_zeros(MAX_SAMPLE_COUNT).unwrap();
        let before_full = full.clone();
        assert_eq!(full.push(0.0), Err(ExactMomentError::CountOverflow));
        assert_eq!(
            full.extend_zeros(u64::MAX),
            Err(ExactMomentError::CountOverflow)
        );
        assert_eq!(full.merge(&before), Err(ExactMomentError::CountOverflow));
        assert_eq!(full, before_full);
    }
}
