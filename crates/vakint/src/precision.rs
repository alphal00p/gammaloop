//! Numerical Laurent coefficients with source-aware, componentwise accuracy.
//!
//! Working precision controls arithmetic. It never increases the accuracy of a
//! finite source. Exact algebra precedes this boundary, and uncertain zeros stay
//! distinguishable from algebraic zeros until checked symbolic conversion.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::ops::Div;

mod materialize;
pub(crate) use materialize::{materialize, snapshot_reference};
pub(crate) mod legacy;

#[cfg(test)]
mod tests;

use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::float::{Complex, Float, Real, RealBall, SingleFloat, TrackedComplex},
};

pub use symbolica::domains::float::{
    ErrorQuality, TrackedComponentPrecision as ComponentPrecision,
    TrackedPrecision as CoefficientPrecision,
};

use crate::{
    Vakint, VakintError, VakintSettings, get_integer_from_atom,
    inspection::FiniteSource,
    utils::vakint_macros::{vk_parse, vk_symbol},
};

/// Source capability, independent of cancellation in an individual evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AvailablePrecision {
    ArbitraryPrecision,
    FixedPrecision(u32),
}

/// A named contributor to numerical uncertainty.
#[derive(Clone, Debug)]
pub struct PrecisionSource {
    pub identity: String,
    pub master_id: Option<String>,
    pub epsilon_power: Option<i64>,
    pub available_decimal_precision: Option<u32>,
    pub source_sha256: String,
}

impl From<&FiniteSource> for PrecisionSource {
    fn from(source: &FiniteSource) -> Self {
        Self {
            identity: source.identity.clone(),
            master_id: source.master_id.clone(),
            epsilon_power: source.epsilon_power,
            available_decimal_precision: Some(source.significant_digits),
            source_sha256: source.source_sha256.clone(),
        }
    }
}

/// A source limit or a final loss of accuracy; retained as structured evidence.
#[derive(Clone, Debug)]
pub struct PrecisionWarning {
    pub epsilon_power: Option<i64>,
    pub requested_decimal_precision: u32,
    pub achieved_decimal_precision: Option<u32>,
    pub sources: Vec<PrecisionSource>,
    pub data_version: String,
    pub message: String,
}

/// One numerical coefficient. Access to the midpoint does not discard its error.
#[derive(Clone, Debug)]
pub struct NumericalCoefficient {
    value: Complex<Float>,
    precision: CoefficientPrecision,
    maximum_export_bits: u32,
}

impl NumericalCoefficient {
    pub(crate) fn from_tracked(value: TrackedComplex, requested: u32) -> Self {
        // Keep the evaluated center and its original enclosure together. A
        // checked symbolic export may round the center, but raw inspection and
        // comparisons must not lose digits merely because the error is larger.
        Self {
            value: Complex::new(value.ball.re.center.clone(), value.ball.im.center.clone()),
            precision: value.precision(),
            maximum_export_bits: decimal_bits(requested),
        }
    }

    /// Import declared numeric values. Explicit zero components are exact;
    /// nonzero components carry the precision already present in each Float.
    pub fn from_value(value: Complex<Float>, requested: u32) -> Self {
        let bits = decimal_bits(requested.saturating_add(32));
        Self::from_tracked(declared_value(value, bits), requested)
    }

    pub fn value(&self) -> &Complex<Float> {
        &self.value
    }
    pub fn precision(&self) -> &CoefficientPrecision {
        &self.precision
    }

    /// Lower a presentation limit while retaining the measured absolute error.
    pub fn capped_decimal_precision(&self, requested: u32) -> Self {
        let value = self.as_tracked();
        let maximum_bits = self.maximum_export_bits.min(decimal_bits(requested));
        let rounded = value.rounded(Some(maximum_bits)).unwrap_or(value);
        let mut result = Self::from_tracked(rounded, requested);
        result.maximum_export_bits = maximum_bits;
        result
    }

    pub(crate) fn as_tracked(&self) -> TrackedComplex {
        TrackedComplex::from_ball(
            Complex::new(
                RealBall::new(
                    self.value.re.clone(),
                    self.precision.real().absolute_error().clone(),
                ),
                RealBall::new(
                    self.value.im.clone(),
                    self.precision.imaginary().absolute_error().clone(),
                ),
            ),
            self.precision
                .real()
                .quality()
                .combine(self.precision.imaginary().quality()),
        )
    }

    /// Export only when both components have a representable accuracy.
    pub fn to_atom(&self) -> Result<Atom, VakintError> {
        self.check_export()?;
        let value = self
            .as_tracked()
            .to_float_with_max_precision(self.maximum_export_bits)
            .map_err(VakintError::NumericalPrecision)?;
        Ok(Atom::num(value))
    }

    /// Real constants and canonical master coefficients have an exact imaginary zero.
    pub fn to_real_float(&self) -> Result<Float, VakintError> {
        self.check_export()?;
        self.real_midpoint()?;
        self.as_tracked()
            .to_float_with_max_precision(self.maximum_export_bits)
            .map(|value| value.re)
            .map_err(VakintError::NumericalPrecision)
    }

    /// Inspect the midpoint of a real value, including an uncertain real zero.
    /// Its precision record remains authoritative; symbolic export is checked.
    pub fn real_midpoint(&self) -> Result<&Float, VakintError> {
        if !self.precision.imaginary().exact_zero() {
            return Err(VakintError::NumericalPrecision(
                "a real coefficient requires an exactly zero imaginary component".into(),
            ));
        }
        Ok(&self.value.re)
    }

    fn check_export(&self) -> Result<(), VakintError> {
        for (name, component) in [
            ("real", self.precision.real()),
            ("imaginary", self.precision.imaginary()),
        ] {
            if !component.exact_zero()
                && !component
                    .significant_digits()
                    .is_some_and(|digits| digits >= 1)
            {
                return Err(VakintError::NumericalPrecision(format!(
                    "{name} component has no reliable significant decimal digit; midpoint {}, absolute uncertainty {}",
                    if name == "real" {
                        &self.value.re
                    } else {
                        &self.value.im
                    },
                    component.absolute_error(),
                )));
            }
        }
        Ok(())
    }

    fn exact_zero(bits: u32) -> Self {
        Self::from_tracked(
            TrackedComplex::from_ball(
                Complex::new(
                    RealBall::exact(Float::new(bits)),
                    RealBall::exact(Float::new(bits)),
                ),
                ErrorQuality::Bounded,
            ),
            ((u64::from(bits) * 30103) / 100000) as u32,
        )
    }
}

/// An ordered Laurent series and the accuracy measured during its evaluation.
#[derive(Clone, Debug)]
pub struct NumericalEvaluationResult {
    coefficients: BTreeMap<i64, NumericalCoefficient>,
    requested_decimal_precision: u32,
    warnings: Vec<PrecisionWarning>,
}

impl Default for NumericalEvaluationResult {
    fn default() -> Self {
        Self::new(16)
    }
}

impl NumericalEvaluationResult {
    pub fn new(requested_decimal_precision: u32) -> Self {
        Self {
            coefficients: BTreeMap::new(),
            requested_decimal_precision,
            warnings: Vec::new(),
        }
    }

    pub fn coefficients(&self) -> &BTreeMap<i64, NumericalCoefficient> {
        &self.coefficients
    }
    pub fn requested_decimal_precision(&self) -> u32 {
        self.requested_decimal_precision
    }
    pub fn warnings(&self) -> &[PrecisionWarning] {
        &self.warnings
    }

    pub fn from_values(input: Vec<(i64, Complex<Float>)>, requested: u32) -> Self {
        let bits = decimal_bits(requested.saturating_add(32));
        let mut combined = BTreeMap::<i64, TrackedComplex>::new();
        let mut sources = BTreeMap::<i64, Vec<PrecisionSource>>::new();
        let mut warnings = Vec::new();
        for (index, (power, value)) in input.into_iter().enumerate() {
            let value = declared_value(value, bits);
            let precision = value.precision();
            let available = [precision.real(), precision.imaginary()]
                .into_iter()
                .filter(|part| !part.exact_zero())
                .map(|part| part.significant_digits().unwrap_or(0))
                .min();
            if let Some(available) = available {
                let source = PrecisionSource {
                    identity: format!("imported epsilon^{power} input {index}"),
                    master_id: None,
                    epsilon_power: Some(power),
                    available_decimal_precision: Some(available),
                    source_sha256: String::new(),
                };
                sources.entry(power).or_default().push(source.clone());
                if available < requested {
                    warnings.push(PrecisionWarning {
                        epsilon_power: Some(power),
                        requested_decimal_precision: requested,
                        achieved_decimal_precision: None,
                        sources: vec![source],
                        data_version: "explicit numerical inputs".into(),
                        message: "imported value cannot support the requested precision".into(),
                    });
                }
            }
            combined
                .entry(power)
                .and_modify(|sum| *sum = sum.add(&value))
                .or_insert(value);
        }
        let coefficients = combined.into_iter().map(|(power, value)| {
            let coefficient = NumericalCoefficient::from_tracked(value, requested);
            let achieved = [coefficient.precision().real(), coefficient.precision().imaginary()].into_iter()
                .filter(|part| !part.exact_zero())
                .map(|part| part.significant_digits().unwrap_or(0)).min();
            if achieved.is_some_and(|digits| digits < requested) {
                warnings.push(PrecisionWarning {
                    epsilon_power: Some(power), requested_decimal_precision: requested,
                    achieved_decimal_precision: achieved,
                    sources: sources.remove(&power).unwrap_or_default(),
                    data_version: "explicit numerical inputs".into(),
                    message: "input uncertainty, repeated-order cancellation and rounding limit this Laurent coefficient".into(),
                });
            }
            (power, coefficient)
        }).collect();
        Self {
            coefficients,
            requested_decimal_precision: requested,
            warnings,
        }
    }

    pub fn to_atom(&self, epsilon_symbol: Symbol) -> Result<Atom, VakintError> {
        let mut result = Atom::Zero;
        for (&power, coefficient) in &self.coefficients {
            let value = coefficient.to_atom().map_err(|error| {
                VakintError::NumericalPrecision(format!("coefficient epsilon^{power}: {error}"))
            })?;
            result += value * Atom::var(epsilon_symbol).pow(Atom::num(power));
        }
        Ok(result)
    }

    pub fn from_atom(
        input: AtomView,
        epsilon_symbol: Symbol,
        settings: &VakintSettings,
    ) -> Result<Self, VakintError> {
        evaluate_series(input, epsilon_symbol, &HashMap::new(), settings, None)
    }

    /// Import decimal values at the requested precision, retaining lower explicit
    /// backtick precision annotations. Constructor warnings are stored in the result.
    pub fn from_vec(
        input: Vec<(i64, (String, String))>,
        settings: &VakintSettings,
    ) -> Result<Self, VakintError> {
        if settings.run_time_decimal_precision == 0 {
            return Err(VakintError::NumericalPrecision(
                "requested decimal precision must be positive".into(),
            ));
        }
        // One input-rounding guard bit accounts for measuring the Float import
        // bound relative to the smallest possible true magnitude.
        let bits = decimal_bits(settings.run_time_decimal_precision).saturating_add(1);
        let parse = |text: &str| -> Result<Float, VakintError> {
            if text
                .split_once('`')
                .is_some_and(|(_, suffix)| !suffix.trim().is_empty())
            {
                // Symbolica interprets the suffix in decimal digits. An explicit
                // binary precision argument overrides it, so inspect the native
                // parsed precision before choosing the lower storage precision.
                let declared = Float::parse(text, None).map_err(VakintError::NumericalPrecision)?;
                if declared.prec() <= bits {
                    return Ok(declared);
                }
            }
            Float::parse(text, Some(bits)).map_err(VakintError::NumericalPrecision)
        };
        let values = input
            .into_iter()
            .map(|(power, (real, imaginary))| {
                Ok((power, Complex::new(parse(&real)?, parse(&imaginary)?)))
            })
            .collect::<Result<Vec<_>, VakintError>>()?;
        Ok(Self::from_values(
            values,
            settings.run_time_decimal_precision,
        ))
    }

    /// Midpoints with each component's achieved Float precision retained.
    pub fn get_epsilon_coefficients(&self) -> Vec<(i64, Complex<Float>)> {
        self.coefficients
            .iter()
            .map(|(&power, coefficient)| (power, coefficient.value.clone()))
            .collect()
    }

    pub fn get_epsilon_coefficient(&self, power: i64) -> Complex<Float> {
        self.coefficients.get(&power).map_or_else(
            || {
                Complex::new(
                    Float::new(decimal_bits(self.requested_decimal_precision)),
                    Float::new(decimal_bits(self.requested_decimal_precision)),
                )
            },
            |coefficient| coefficient.value.clone(),
        )
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.values().all(|coefficient| {
            coefficient.precision.real().exact_zero()
                && coefficient.precision.imaginary().exact_zero()
        })
    }

    pub(crate) fn align_orders(&mut self, other: &Self) {
        for &power in other.coefficients.keys() {
            self.coefficients.entry(power).or_insert_with(|| {
                NumericalCoefficient::exact_zero(decimal_bits(self.requested_decimal_precision))
            });
        }
    }

    /// Independent statistical errors combine in quadrature for each component.
    pub fn aggregate_errors(&self, other: &Self) -> Result<Self, VakintError> {
        let requested = self
            .requested_decimal_precision
            .min(other.requested_decimal_precision);
        let bits = decimal_bits(requested.saturating_add(32));
        let orders: BTreeSet<_> = self
            .coefficients
            .keys()
            .chain(other.coefficients.keys())
            .copied()
            .collect();
        let mut result = Self::new(requested);
        let names = [
            "aggregate_error_real_a",
            "aggregate_error_real_b",
            "aggregate_error_imaginary_a",
            "aggregate_error_imaginary_b",
        ];
        let formula = vk_parse!("sqrt(aggregate_error_real_a^2+aggregate_error_real_b^2)+𝑖*sqrt(aggregate_error_imaginary_a^2+aggregate_error_imaginary_b^2)")
            .map_err(|error| VakintError::EvaluationError(error.to_string()))?;
        for order in orders {
            let a = self.coefficients.get(&order).map_or_else(
                || Self::zero_tracked(bits),
                NumericalCoefficient::as_tracked,
            );
            let b = other.coefficients.get(&order).map_or_else(
                || Self::zero_tracked(bits),
                NumericalCoefficient::as_tracked,
            );
            let quality = a.error_quality.combine(b.error_quality);
            let mut inputs = HashMap::new();
            for (name, component) in names
                .iter()
                .zip([a.ball.re, b.ball.re, a.ball.im, b.ball.im])
            {
                inputs.insert(
                    Atom::var(vk_symbol!(*name)),
                    TrackedComplex::from_ball(
                        Complex::new(component, RealBall::exact(Float::new(bits))),
                        quality,
                    ),
                );
            }
            let coefficient = formula
                .evaluate_tracked(&inputs, bits)
                .map_err(|error| VakintError::NumericalPrecision(error.to_string()))?;
            result.coefficients.insert(
                order,
                NumericalCoefficient::from_tracked(coefficient, requested),
            );
        }
        result
            .warnings
            .extend(self.warnings.iter().chain(&other.warnings).cloned());
        Ok(result)
    }

    fn zero_tracked(bits: u32) -> TrackedComplex {
        TrackedComplex::from_ball(
            Complex::new(
                RealBall::exact(Float::new(bits)),
                RealBall::exact(Float::new(bits)),
            ),
            ErrorQuality::Bounded,
        )
    }
}

impl NumericalEvaluationResult {
    pub fn does_approx_match(
        &self,
        other: &NumericalEvaluationResult,
        error: Option<&NumericalEvaluationResult>,
        threshold: f64,
        max_pull: f64,
    ) -> (bool, String) {
        let mut powers = HashSet::new();
        for power in self.coefficients.keys() {
            powers.insert(*power);
        }
        for power in other.coefficients.keys() {
            powers.insert(*power);
        }
        let mut power_vec = powers.iter().cloned().collect::<Vec<_>>();
        power_vec.sort();

        for power in power_vec {
            let (self_val, other_val) = (
                self.get_epsilon_coefficient(power),
                other.get_epsilon_coefficient(power),
            );
            let comparisons = [
                (
                    "real",
                    self_val.re.clone(),
                    other_val.re.clone(),
                    error.map(|e| e.get_epsilon_coefficient(power).re.clone()),
                ),
                (
                    "imaginary",
                    self_val.im.clone(),
                    other_val.im.clone(),
                    error.map(|e| e.get_epsilon_coefficient(power).im.clone()),
                ),
            ];
            for (part, r, o, e) in comparisons.iter() {
                let f_max_pull = Float::with_val(r.prec(), max_pull);
                let f_threshold = Float::with_val(r.prec(), threshold);
                let delta = (r.clone() - o).norm();
                if let Some(err) = e
                    && !err.norm().is_zero()
                    && delta > f_max_pull * err.norm()
                {
                    return (
                        false,
                        format!(
                            "{} part of ε^{} coefficient does not match within max pull: {} != {} (pull = {})",
                            part,
                            power,
                            r,
                            o,
                            delta / err.norm()
                        ),
                    );
                }
                let scale = (r.norm() + o.norm()) / Float::with_val(r.prec(), 2.0);
                if scale.is_zero() {
                    if delta > f_threshold {
                        return (
                            false,
                            format!(
                                "{} part of ε^{} coefficient does not match within abs. error required: {} != {} (abs. error = {})",
                                part, power, r, o, delta
                            ),
                        );
                    }
                } else {
                    let rel_error = delta.div(scale);
                    if rel_error > f_threshold {
                        return (
                            false,
                            format!(
                                "{} part of ε^{} coefficient does not match within rel. error required: {} != {} (rel. error = {})",
                                part, power, r, o, rel_error
                            ),
                        );
                    }
                }
            }
        }
        (true, "matches!".into())
    }
}

pub(crate) fn decimal_bits(digits: u32) -> u32 {
    // Upper bound on log2(10), plus the leading-bit convention used by the
    // tracked relative error. No f64 conversion or decimal-digit truncation.
    ((u128::from(digits) * 3_321_928_094_887_362_348).div_ceil(1_000_000_000_000_000_000) + 1)
        .clamp(2, u128::from(u32::MAX)) as u32
}

/// A literal/parameter zero component declares an exact zero. Uncertain zeros
/// produced by arithmetic enter through TrackedComplex and never use this path.
fn declared_value(value: Complex<Float>, bits: u32) -> TrackedComplex {
    let real_zero = value.re.is_zero();
    let imaginary_zero = value.im.is_zero();
    let mut tracked = TrackedComplex::from_float(value, bits);
    if real_zero {
        tracked.ball.re = RealBall::exact(Float::new(bits));
    }
    if imaginary_zero {
        tracked.ball.im = RealBall::exact(Float::new(bits));
    }
    tracked
}

pub(crate) fn evaluate_coefficient(
    expression: AtomView,
    sources: &BTreeMap<Atom, FiniteSource>,
    parameters: &HashMap<Atom, Complex<Float>>,
    settings: &VakintSettings,
) -> Result<(NumericalCoefficient, Vec<PrecisionWarning>), VakintError> {
    evaluate_coefficient_with_working_precision(
        expression,
        sources,
        parameters,
        settings,
        settings.run_time_decimal_precision,
    )
}

/// [`evaluate_coefficient`] with the arithmetic carried at `working` decimal
/// digits (plus the usual guard) instead of the requested precision.
///
/// Warnings and the export limit still refer to the request. Callers lower
/// `working` only when finite sources cap the result below the request anyway.
pub(crate) fn evaluate_coefficient_with_working_precision(
    expression: AtomView,
    sources: &BTreeMap<Atom, FiniteSource>,
    parameters: &HashMap<Atom, Complex<Float>>,
    settings: &VakintSettings,
    working: u32,
) -> Result<(NumericalCoefficient, Vec<PrecisionWarning>), VakintError> {
    Vakint::initialize_vakint_symbols();
    for identity in parameters.keys() {
        crate::inspection::validate_parameter(identity.as_view(), None)?;
    }
    let requested = settings.run_time_decimal_precision;
    if requested == 0 || working == 0 {
        return Err(VakintError::NumericalPrecision(
            "requested decimal precision must be positive".into(),
        ));
    }
    let bits = decimal_bits(
        working
            .min(requested)
            .checked_add(32)
            .ok_or_else(|| VakintError::NumericalPrecision("working precision overflow".into()))?,
    );
    // Give finite literals identities before expansion can combine their
    // midpoints. Exact polynomial cancellation then remains safe and complete.
    let staged = legacy::stage_table(
        Atom::var(vk_symbol!("numerical_coefficient")).as_view(),
        expression,
        vk_symbol!(&settings.epsilon_symbol),
        "Input",
    )?;
    let normalized = Vakint::normalize_gamma_functions(staged.as_view());
    let normalized = crate::rustred_evaluation::combine_exact(crate::hpl::normalize_constants(
        normalized.as_view(),
    )?);
    let mut combined_sources = sources.clone();
    combined_sources.extend(legacy::collect_sources(normalized.as_view())?);
    let sources = &combined_sources;
    let mut inputs = HashMap::new();
    // Match the legacy constant map: this exact default precedes caller inputs.
    inputs.insert(
        Atom::var(crate::symbols::S.cmplx_i),
        TrackedComplex::from_exact_rational(&Complex::new(0.into(), 1.into()), bits),
    );
    let mut used_sources = Vec::new();
    let mut versions = BTreeSet::new();
    let mut warnings = Vec::new();
    for (symbol, source) in sources {
        if normalized
            .pattern_match(&symbol.to_pattern(), None, None)
            .next()
            .is_none()
        {
            continue;
        }
        let value = if let Some(value) = legacy::input_value(symbol.as_view())? {
            // A caller's Float is a binary source, not a decimal table entry.
            // In particular, a valid 1–4 bit cancellation residue has zero
            // complete decimal digits but still has a usable absolute error.
            declared_value(Complex::new(value, Float::new(bits)), bits)
        } else {
            let limit =
                TrackedComplex::from_source_decimal(&source.value, source.significant_digits, bits)
                    .map_err(|error| {
                        VakintError::NumericalPrecision(format!("{}: {error}", source.identity))
                    })?;
            if let Some(bound) = &source.absolute_error {
                let mut value = TrackedComplex::from_decimal_with_error(
                    &source.value,
                    bound,
                    bits,
                    if source.error_is_bounded {
                        ErrorQuality::Bounded
                    } else {
                        ErrorQuality::Estimated
                    },
                )
                .map_err(|error| {
                    VakintError::NumericalPrecision(format!("{}: {error}", source.identity))
                })?;
                // Stored guard digits and a tighter certificate do not override the
                // dataset's declared maximum validated significant digits. A larger
                // certified error remains authoritative in the other direction.
                if value.ball.re.radius < limit.ball.re.radius {
                    value.ball.re.radius = limit.ball.re.radius;
                }
                value
            } else {
                limit
            }
        };
        inputs.insert(symbol.clone(), value);
        let identity = PrecisionSource::from(source);
        used_sources.push(identity.clone());
        versions.insert(source.data_version.clone());
        if requested > source.significant_digits {
            warnings.push(PrecisionWarning {
                epsilon_power: None,
                requested_decimal_precision: requested,
                achieved_decimal_precision: None,
                sources: vec![identity],
                data_version: source.data_version.clone(),
                message: format!(
                    "source {} cannot support the requested precision",
                    source.identity
                ),
            });
        }
    }
    for (symbol, value) in parameters {
        if normalized
            .pattern_match(&symbol.to_pattern(), None, None)
            .next()
            .is_none()
        {
            continue;
        }
        let available = [&value.re, &value.im]
            .into_iter()
            .filter(|component| !component.is_zero())
            .map(|component| {
                ((u128::from(component.prec().saturating_sub(1)) * 301_029_995_663_981_195)
                    / 1_000_000_000_000_000_000) as u32
            })
            .min();
        let identity = PrecisionSource {
            identity: symbol.to_canonical_string(),
            master_id: None,
            epsilon_power: None,
            available_decimal_precision: available,
            source_sha256: String::new(),
        };
        if available.is_some() {
            used_sources.push(identity.clone());
            versions.insert("runtime inputs".into());
        }
        if available.is_some_and(|available| requested > available) {
            warnings.push(PrecisionWarning {
                epsilon_power: None,
                requested_decimal_precision: requested,
                achieved_decimal_precision: None,
                sources: vec![identity],
                data_version: "runtime inputs".into(),
                message: format!("parameter {symbol} is limited by its stored Float precision"),
            });
        }
        inputs.insert(symbol.clone(), declared_value(value.clone(), bits));
    }
    let mut projections = BTreeMap::new();
    let mut issue = None;
    normalized.visitor(&mut |part| match crate::hpl::projection_kernel(part) {
        Ok(Some(kernel)) => {
            projections.insert(part.to_owned(), kernel);
            false
        }
        Ok(None) => true,
        Err(error) => {
            issue = Some(error);
            false
        }
    });
    if let Some(error) = issue {
        return Err(error);
    }
    for (symbol, (kernel, imaginary)) in projections {
        let value = kernel
            .evaluate_tracked(&inputs, bits)
            .map_err(|error| VakintError::EvaluationError(error.to_string()))?;
        let component = if imaginary {
            value.ball.im
        } else {
            value.ball.re
        };
        inputs.insert(
            symbol,
            TrackedComplex::from_ball(
                Complex::new(component, RealBall::exact(Float::new(bits))),
                value.error_quality,
            ),
        );
    }
    let value = normalized
        .evaluate_tracked(&inputs, bits)
        .map_err(|error| {
            VakintError::EvaluationError(format!(
                "coefficient {}: {error}",
                normalized.to_canonical_string()
            ))
        })?;
    let coefficient = NumericalCoefficient::from_tracked(value, requested);
    let achieved = [
        coefficient.precision.real(),
        coefficient.precision.imaginary(),
    ]
    .into_iter()
    .filter(|part| !part.exact_zero())
    .map(|part| part.significant_digits().unwrap_or(0))
    .min();
    if achieved.is_some_and(|digits| digits < requested) {
        if used_sources.is_empty() {
            used_sources.push(PrecisionSource {
                identity: normalized.to_canonical_string(),
                master_id: None,
                epsilon_power: None,
                available_decimal_precision: None,
                source_sha256: String::new(),
            });
        }
        warnings.push(PrecisionWarning {
            epsilon_power: None,
            requested_decimal_precision: requested,
            achieved_decimal_precision: achieved,
            sources: used_sources,
            data_version: versions.into_iter().collect::<Vec<_>>().join(", "),
            message: "source uncertainty, cancellation and rounding limit this Laurent coefficient"
                .into(),
        });
    }
    Ok((coefficient, warnings))
}

/// Bound intermediate epsilon degrees before entering Symbolica's polynomial
/// collector. Sums may acquire a common denominator, so summing the children's
/// budgets also covers the numerator degree of a rational Laurent expression.
/// Function arguments are opaque coefficients, including source identity tags.
fn checked_laurent_degree(expression: AtomView, epsilon: Symbol) -> Result<u64, VakintError> {
    let overflow = || {
        VakintError::EvaluationError(format!(
            "epsilon degree budget exceeds {}; numerical Laurent collection requires intermediate degrees to fit a signed 32-bit exponent",
            i32::MAX
        ))
    };
    let degree = match expression {
        AtomView::Var(variable) => u64::from(variable.get_symbol() == epsilon),
        AtomView::Num(_) | AtomView::Fun(_) => 0,
        AtomView::Add(sum) => sum.iter().try_fold(0_u64, |degree, term| {
            degree
                .checked_add(checked_laurent_degree(term, epsilon)?)
                .ok_or_else(overflow)
        })?,
        AtomView::Mul(product) => product.iter().try_fold(0_u64, |degree, factor| {
            degree
                .checked_add(checked_laurent_degree(factor, epsilon)?)
                .ok_or_else(overflow)
        })?,
        AtomView::Pow(power) => {
            let degree = checked_laurent_degree(power.get_base(), epsilon)?;
            if degree == 0 {
                0
            } else if let Some(exponent) = get_integer_from_atom(power.get_exp()) {
                degree
                    .checked_mul(exponent.unsigned_abs())
                    .ok_or_else(overflow)?
            } else if matches!(power.get_exp(), AtomView::Num(number) if number.get_coeff_view().is_integer())
            {
                return Err(overflow());
            } else {
                return Err(VakintError::EvaluationError(format!(
                    "epsilon power is not an integer Laurent monomial: {expression}"
                )));
            }
        }
    };
    if degree > i32::MAX as u64 {
        Err(overflow())
    } else {
        Ok(degree)
    }
}

pub(crate) fn evaluate_series(
    expression: AtomView,
    epsilon: Symbol,
    parameters: &HashMap<Atom, Complex<Float>>,
    settings: &VakintSettings,
    catalog: Option<&crate::rustred_evaluation::masters::MasterCatalog>,
) -> Result<NumericalEvaluationResult, VakintError> {
    Vakint::initialize_vakint_symbols();
    for identity in parameters.keys() {
        crate::inspection::validate_parameter(identity.as_view(), catalog)?;
    }
    checked_laurent_degree(expression, epsilon)?;
    let expression = legacy::stage_table(
        Atom::var(vk_symbol!("numerical_expression")).as_view(),
        expression,
        epsilon,
        "Input",
    )?;
    let expression = if catalog.is_none() {
        crate::rustred_evaluation::combine_exact(expression)
    } else {
        expression
    };
    if catalog.is_none() {
        let mut references = BTreeSet::new();
        expression.visitor(&mut |part| {
            if materialize::snapshot_reference(part) {
                references.insert(part.to_canonical_string());
                return false;
            }
            !legacy::is_source(part)
        });
        if !references.is_empty() {
            return Err(VakintError::NumericalPrecision(format!(
                "unresolved snapshot-bound masters or coefficients: {}; use numerical evaluation on the Vakint instance that owns their loaded data",
                references.into_iter().collect::<Vec<_>>().join("; "),
            )));
        }
    }
    // Validate retained, uncomputed and unknown terminals at the registry
    // boundary, before any evaluator can descend into their identity tags.
    // This also cancels exact zero combinations before reporting limitations.
    let expression = if let Some(catalog) = catalog {
        catalog.resolve_numerical(expression, settings)?
    } else {
        expression
    };
    let mut result = NumericalEvaluationResult::new(settings.run_time_decimal_precision);
    checked_laurent_degree(expression.as_view(), epsilon)?;
    for (power, coefficient) in expression.coefficient_list::<i32>(&[Atom::var(epsilon)]) {
        let order = if power == Atom::num(1) {
            0
        } else if power == Atom::var(epsilon) {
            1
        } else {
            match power.as_view() {
                AtomView::Pow(power) if power.get_base() == Atom::var(epsilon).as_view() => {
                    get_integer_from_atom(power.get_exp())
                }
                _ => None,
            }
            .ok_or_else(|| {
                VakintError::EvaluationError(format!(
                    "epsilon power is not an integer monomial: {power}"
                ))
            })?
        };
        let coefficient = coefficient.to_owned();
        let (coefficient, sources, origins) = if let Some(catalog) = catalog {
            let prepared = catalog.prepare_evaluation(coefficient)?;
            (prepared.expression, prepared.sources, prepared.origins)
        } else {
            (coefficient, BTreeMap::new(), Vec::new())
        };
        let (value, mut warnings) =
            evaluate_coefficient(coefficient.as_view(), &sources, parameters, settings)?;
        for warning in &mut warnings {
            attach_origins(
                std::slice::from_mut(warning),
                &origins,
                Some(order),
                catalog.map_or("", |catalog| catalog.version.as_str()),
            );
            log::warn!("{warning}");
        }
        result.coefficients.insert(order, value);
        result.warnings.extend(warnings);
    }
    Ok(result)
}

pub(crate) fn attach_origins(
    warnings: &mut [PrecisionWarning],
    origins: &[PrecisionSource],
    epsilon_power: Option<i64>,
    data_version: &str,
) {
    for warning in warnings {
        warning.epsilon_power = epsilon_power;
        if !data_version.is_empty()
            && !warning
                .data_version
                .split(", ")
                .any(|version| version == data_version)
        {
            if !warning.data_version.is_empty() {
                warning.data_version.push_str(", ");
            }
            warning.data_version.push_str(data_version);
        }
        if warning.data_version.is_empty() {
            warning.data_version = "standalone expression".into();
        }
        {
            let mut seen = warning
                .sources
                .iter()
                .map(|source| {
                    (
                        source.identity.clone(),
                        source.master_id.clone(),
                        source.epsilon_power,
                        source.available_decimal_precision,
                        source.source_sha256.clone(),
                    )
                })
                .collect::<BTreeSet<_>>();
            for source in origins.iter().filter(|source| {
                warning.achieved_decimal_precision.is_some()
                    || source
                        .available_decimal_precision
                        .is_some_and(|cap| cap < warning.requested_decimal_precision)
            }) {
                if seen.insert((
                    source.identity.clone(),
                    source.master_id.clone(),
                    source.epsilon_power,
                    source.available_decimal_precision,
                    source.source_sha256.clone(),
                )) {
                    warning.sources.push(source.clone());
                }
            }
        }
    }
}

/// Build symbolic momentum contractions so their products and cancellations are
/// performed by the tracked evaluator, rather than rounded before it starts.
pub(crate) fn evaluation_inputs(
    expression: AtomView,
    epsilon: Symbol,
    real: &HashMap<String, Float, ahash::RandomState>,
    complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
    externals: Option<&HashMap<usize, crate::Momentum, ahash::RandomState>>,
) -> Result<(Atom, HashMap<Atom, Complex<Float>>), VakintError> {
    Vakint::initialize_vakint_symbols();
    checked_laurent_degree(expression, epsilon)?;
    // Dot conversion expands its input. Retain independent finite estimates
    // before that expansion, just as at the final numerical boundary. Input
    // staging traverses the complete expression and does not split epsilon orders.
    let expression = legacy::stage_table(
        Atom::var(vk_symbol!("evaluation_input")).as_view(),
        expression,
        vk_symbol!("ep"),
        "Input",
    )?;
    let mut expression = Vakint::convert_to_dot_notation(expression.as_view());
    let mut parameters = HashMap::new();
    for (name, value) in real {
        let symbol = vk_parse!(name).map_err(|error| {
            VakintError::EvaluationError(format!("invalid parameter {name}: {error}"))
        })?;
        parameters.insert(
            symbol,
            Complex::new(value.clone(), Float::new(value.prec())),
        );
    }
    for (name, value) in complex {
        let symbol = vk_parse!(name).map_err(|error| {
            VakintError::EvaluationError(format!("invalid parameter {name}: {error}"))
        })?;
        parameters.insert(symbol, value.clone());
    }
    let mut external_dots = BTreeMap::new();
    if let Some(externals) = externals {
        let mut momenta = BTreeMap::new();
        for (&id, momentum) in externals {
            let values = match momentum {
                crate::Momentum::Real(value) => [&value.0, &value.1, &value.2, &value.3]
                    .map(|value| Complex::new(value.clone(), Float::new(value.prec()))),
                crate::Momentum::Complex(value) => [
                    value.0.clone(),
                    value.1.clone(),
                    value.2.clone(),
                    value.3.clone(),
                ],
            };
            let components = values
                .into_iter()
                .enumerate()
                .map(|(component, value)| {
                    let key = external_component(id, component);
                    parameters.insert(key.clone(), value);
                    key
                })
                .collect::<Vec<_>>();
            momenta.insert(id, components);
        }
        for (&left_id, left) in &momenta {
            for (&right_id, right) in &momenta {
                if left_id > right_id {
                    continue;
                }
                let dot = symbolica::function!(
                    crate::symbols::S.dot,
                    symbolica::function!(crate::symbols::S.p, Atom::num(left_id)),
                    symbolica::function!(crate::symbols::S.p, Atom::num(right_id))
                );
                let value = left[0].clone() * right[0].clone()
                    - left[1].clone() * right[1].clone()
                    - left[2].clone() * right[2].clone()
                    - left[3].clone() * right[3].clone();
                external_dots.insert(dot, value);
            }
        }
    }
    let dimension = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    let mut issue = None;
    expression = expression.replace_map(|part, _, out| {
        // Input bindings own complete subtrees, before either dot notation is
        // resolved. Both notations reuse the same component source identities.
        if parameters.contains_key(&part.to_owned()) || issue.is_some() {
            **out = part.to_owned();
            return;
        }
        let AtomView::Fun(function) = part else { return; };
        if function.get_symbol() == crate::SPENSO_TAG.pure_scalar
            || function.get_symbol().get_evaluation_info().is_some_and(|info|
                info.get_tag_count() > 0 && info.get_tag_count() == function.get_nargs())
        {
            **out = part.to_owned();
            return;
        }
        if let Some(value) = external_dots.get(&part.to_owned()) {
            **out = value.clone();
            return;
        }
        if ![crate::symbols::S.dot, crate::SPENSO_TAG.dot].contains(&function.get_symbol())
            || !function.iter().any(|argument| {
                matches!(argument, AtomView::Fun(vector)
                    if vector.get_symbol() == crate::symbols::S.p
                        && vector.iter().last().is_some_and(crate::VakintTerm::is_compact_spenso_minkowski_representation))
            })
        {
            return;
        }
        let failure = |detail: String| VakintError::EvaluationError(format!(
            "native external-vector dot {}: {detail}", part.to_canonical_string()
        ));
        if function.get_nargs() != 2 {
            issue = Some(failure("expected two external vectors".into()));
            return;
        }
        let mut ids = Vec::with_capacity(2);
        for argument in function.iter() {
            let AtomView::Fun(vector) = argument else {
                issue = Some(failure("expected external p(id) vectors".into()));
                return;
            };
            if vector.get_symbol() != crate::symbols::S.p
                || !matches!(vector.get_nargs(), 1 | 2)
            {
                issue = Some(failure("expected external p(id) vectors".into()));
                return;
            }
            if vector.get_nargs() == 2 {
                let representation = vector.get(1);
                if !crate::VakintTerm::is_compact_spenso_minkowski_representation(representation) {
                    issue = Some(failure("expected a compact Minkowski representation".into()));
                    return;
                }
                let AtomView::Fun(representation) = representation else {
                    issue = Some(failure("expected a compact Minkowski representation".into()));
                    return;
                };
                let source_dimension = representation.get(0);
                if !(source_dimension.to_owned() - &dimension).expand().is_zero() {
                    issue = Some(failure(format!(
                        "source dimension {source_dimension} differs from requested dimension {dimension}"
                    )));
                    return;
                }
            }
            let Ok(id) = usize::try_from(vector.get(0)) else {
                issue = Some(failure(format!("invalid external momentum id {}", vector.get(0))));
                return;
            };
            ids.push(id);
        }
        let bare = symbolica::function!(
            crate::symbols::S.dot,
            symbolica::function!(crate::symbols::S.p, ids[0]),
            symbolica::function!(crate::symbols::S.p, ids[1])
        );
        if let Some(value) = external_dots.get(&bare) {
            **out = value.clone();
        } else {
            issue = Some(failure(format!(
                "missing external momenta for p({}) and p({})", ids[0], ids[1]
            )));
        }
    });
    if let Some(issue) = issue {
        return Err(issue);
    }
    Ok((expression, parameters))
}

/// A scalar identifier that cannot be mistaken for indexed-vector syntax.
pub(crate) fn external_component(momentum: usize, component: usize) -> Atom {
    Atom::var(symbolica::symbol!(format!(
        "vakint_external_components::p{momentum}_{component}"
    )))
}

pub(crate) fn evaluate_with_errors(
    expression: AtomView,
    parameters: &HashMap<Atom, Complex<Float>>,
    settings: &VakintSettings,
    catalog: Option<&crate::rustred_evaluation::masters::MasterCatalog>,
) -> Result<(NumericalEvaluationResult, Option<NumericalEvaluationResult>), VakintError> {
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    if expression
        .pattern_match(&crate::symbols::S.error_flag.to_pattern(), None, None)
        .next()
        .is_none()
    {
        return Ok((
            evaluate_series(expression, epsilon, parameters, settings, catalog)?,
            None,
        ));
    }
    let (error, central) = crate::utils::split_linear_atom(
        expression,
        Atom::var(crate::symbols::S.error_flag_symbol).as_view(),
    );
    let mut central = evaluate_series(central.as_view(), epsilon, parameters, settings, catalog)?;
    let mut error = evaluate_series(error.as_view(), epsilon, parameters, settings, catalog)?;
    central.align_orders(&error);
    error.align_orders(&central);
    Ok((central, Some(error)))
}
