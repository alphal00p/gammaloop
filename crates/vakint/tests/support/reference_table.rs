//! Readable result tables in the spirit of the original MATAD/FMFT tests.
//!
//! A [`LaurentTable`] prints one block per integral: what was computed, then one
//! line per Laurent coefficient with vakint's value and, below it, every
//! reference value with its relative difference and verdict, then a one-line
//! summary (`PASS`/`FAIL`). The relative difference of a coefficient is
//! `max(|Δ re|, |Δ im|) / max(|re|, |im|)` of the reference (the absolute
//! difference, printed as `Δabs`, when the reference coefficient is zero).
//! Every order of the requested window is checked; a coefficient a side does
//! not list is an exact zero, so a reference can state that an order vanishes
//! (write it as `"0"`) and a spurious vakint pole fails.

use std::ops::RangeInclusive;

use symbolica::domains::float::{Complex, Float, RoundingDirection};

/// Binary precision of the comparison arithmetic (about 616 decimal digits,
/// above every stored reference).
const WORKING_BITS: u32 = 2048;

/// One set of reference Laurent coefficients and the agreement it requires.
pub struct Reference {
    /// Who produced the values, e.g. `"FMFT"` or `"pySecDec (Monte Carlo)"`.
    pub source: String,
    /// Largest accepted relative difference of a coefficient.
    pub tolerance: f64,
    /// `(ε power, real part, imaginary part)` as decimal strings.
    pub coefficients: Vec<(i64, String, String)>,
}

impl Reference {
    pub fn new(source: &str, tolerance: f64, coefficients: &[(i64, &str, &str)]) -> Self {
        Self {
            source: source.to_owned(),
            tolerance,
            coefficients: coefficients
                .iter()
                .map(|(power, real, imaginary)| {
                    (*power, (*real).to_owned(), (*imaginary).to_owned())
                })
                .collect(),
        }
    }
}

/// The Laurent coefficients vakint computed, with a description of the case.
pub struct LaurentTable {
    /// Integral or topology name, printed as the block title.
    pub title: String,
    /// `(label, text)` lines describing input, point, scheme and engine.
    pub facts: Vec<(String, String)>,
    /// Label of the computed values, e.g. `"vakint"`.
    pub engine: String,
    /// Computed `(ε power, value)` pairs; vanishing orders may be left out.
    pub values: Vec<(i64, Complex<Float>)>,
    /// Laurent orders the computation was asked for. vakint's coefficient of a
    /// requested order it does not list is an exact zero; `None` means the
    /// orders either side lists.
    pub orders: Option<RangeInclusive<i64>>,
    /// Significant digits printed per value.
    pub shown_digits: usize,
}

fn working(value: &Float) -> Float {
    Float::from_rational_round(
        &value.to_rational(),
        WORKING_BITS,
        RoundingDirection::Nearest,
    )
}

fn parse(text: &str) -> Result<Float, String> {
    Float::parse(text.trim(), Some(WORKING_BITS)).map_err(|error| format!("{text:?}: {error}"))
}

fn absolute(value: Float) -> Float {
    if value.is_negative() { -value } else { value }
}

fn larger(left: Float, right: Float) -> Float {
    if left > right { left } else { right }
}

fn is_zero(value: &Float) -> bool {
    value.to_rational().is_zero()
}

fn exact_zero() -> Complex<Float> {
    Complex::new(
        parse("0").expect("decimal zero"),
        parse("0").expect("decimal zero"),
    )
}

/// `±d.ddd…e±x` with `digits` significant digits.
pub fn decimal(value: &Float, digits: usize) -> String {
    if is_zero(value) {
        return "0".into();
    }
    // The precision of a multiple-precision float's `e` format counts significant digits.
    let text = format!("{:.*e}", digits.max(1), absolute(value.clone()));
    format!("{}{text}", if value.is_negative() { "-" } else { "+" })
}

/// Significant decimal digits written in a decimal string (`"-6.0152e-5"` has 5).
fn written_digits(text: &str) -> usize {
    let mantissa = text.trim().split(['e', 'E']).next().unwrap_or("");
    let digits = mantissa
        .chars()
        .filter(char::is_ascii_digit)
        .skip_while(|digit| *digit == '0')
        .count();
    digits.max(1)
}

/// A complex value showing only the components that matter at `digits`.
pub fn complex(value: &Complex<Float>, digits: usize) -> String {
    let scale = larger(absolute(value.re.clone()), absolute(value.im.clone()));
    let threshold = scale * &parse(&format!("1e-{}", digits + 2)).expect("decimal threshold");
    let show_re = absolute(value.re.clone()) > threshold;
    let show_im = absolute(value.im.clone()) > threshold;
    match (show_re, show_im) {
        (false, false) => "0".into(),
        (true, false) => decimal(&value.re, digits),
        (false, true) => format!("{} i", decimal(&value.im, digits)),
        (true, true) => format!(
            "{} {} i",
            decimal(&value.re, digits),
            decimal(&value.im, digits)
        ),
    }
}

/// Relative difference of `value` to `reference` (absolute when the reference vanishes).
pub fn relative_difference(value: &Complex<Float>, reference: &Complex<Float>) -> Float {
    let (re, im) = (working(&value.re), working(&value.im));
    let (reference_re, reference_im) = (working(&reference.re), working(&reference.im));
    let difference = larger(absolute(re - &reference_re), absolute(im - &reference_im));
    let scale = larger(absolute(reference_re), absolute(reference_im));
    if is_zero(&scale) {
        difference
    } else {
        difference / &scale
    }
}

fn order(power: i64) -> String {
    format!("ε^{power}")
}

impl LaurentTable {
    /// Print the table and check it against every reference the way
    /// `NumericalEvaluationResult::does_approx_match` compares results: at every
    /// order of the requested window [`Self::orders`] and at every order that
    /// vakint or a reference lists, a coefficient a side does not list is an
    /// exact zero. vakint leaves vanishing orders out of its output, so a
    /// spurious vakint pole fails against the (listed or implied) zero of the
    /// reference, and a reference written as `"0"` states that the order
    /// vanishes. Only a reference order outside the window that vakint did not
    /// compute fails as `NOT COMPUTED`.
    pub fn check(&self, references: &[Reference]) -> Result<(), String> {
        let mut text = format!("\n━━ {} ━━\n", self.title);
        for (label, fact) in &self.facts {
            text.push_str(&format!("  {label:<9} {fact}\n"));
        }
        let parsed = references
            .iter()
            .map(|reference| {
                reference
                    .coefficients
                    .iter()
                    .map(|(power, re, im)| {
                        // Never print more digits than the reference was given with.
                        let written = written_digits(re).max(written_digits(im));
                        Ok((*power, Complex::new(parse(re)?, parse(im)?), written))
                    })
                    .collect::<Result<Vec<_>, String>>()
            })
            .collect::<Result<Vec<_>, String>>()?;
        let in_window = |power: i64| {
            self.orders
                .as_ref()
                .is_none_or(|orders| orders.contains(&power))
        };
        let mut powers = self
            .orders
            .clone()
            .into_iter()
            .flatten()
            .chain(self.values.iter().map(|(power, _)| *power))
            .chain(parsed.iter().flatten().map(|(power, _, _)| *power))
            .collect::<Vec<_>>();
        powers.sort_unstable();
        powers.dedup();
        let width = references
            .iter()
            .map(|reference| reference.source.chars().count())
            .chain([self.engine.chars().count()])
            .max()
            .unwrap_or(0);
        let column = self.shown_digits + 12;
        let mut failures = Vec::new();
        let mut checked = 0;
        for power in powers {
            let listed = self
                .values
                .iter()
                .find(|(computed, _)| *computed == power)
                .map(|(_, value)| value.clone());
            // A requested order vakint does not list vanishes.
            let value = listed.clone().or_else(|| in_window(power).then(exact_zero));
            let shown = match (&listed, &value) {
                (Some(listed), _) => complex(listed, self.shown_digits),
                (None, Some(_)) => "0 (no such term)".to_owned(),
                (None, None) => "(not computed)".to_owned(),
            };
            text.push_str(&format!(
                "  {:<6} {:<width$}  {shown}\n",
                order(power),
                self.engine
            ));
            for (reference, coefficients) in references.iter().zip(&parsed) {
                let listed_reference = coefficients.iter().find(|(at, _, _)| *at == power);
                let (expected, reference_text) = match listed_reference {
                    Some((_, expected, written)) => (
                        expected.clone(),
                        complex(expected, self.shown_digits.min(*written)),
                    ),
                    // An order a reference does not list vanishes for it.
                    None if value.is_some() => (exact_zero(), "0 (not listed)".to_owned()),
                    None => continue,
                };
                checked += 1;
                let zero_reference = is_zero(&expected.re) && is_zero(&expected.im);
                // The difference is absolute when the reference vanishes.
                let measure = if zero_reference { "Δabs" } else { "Δrel" };
                let (verdict, difference) = match &value {
                    Some(value) => {
                        let difference = relative_difference(value, &expected);
                        let tolerance = Float::with_val(64, reference.tolerance);
                        let ok = difference <= tolerance;
                        if !ok {
                            failures.push(format!(
                                "{} vs {}: {measure} {:.1e} > {:.0e}",
                                order(power),
                                reference.source,
                                difference,
                                reference.tolerance
                            ));
                        }
                        (
                            if ok { "ok" } else { "MISMATCH" },
                            if is_zero(&difference) {
                                "0".to_owned()
                            } else {
                                format!("{difference:.1e}")
                            },
                        )
                    }
                    None => {
                        failures.push(format!(
                            "{} vs {}: vakint did not compute this order (outside the requested orders)",
                            order(power),
                            reference.source
                        ));
                        ("NOT COMPUTED", "-".to_owned())
                    }
                };
                text.push_str(&format!(
                    "  {:<6} {:<width$}  {reference_text:<column$} {measure} {difference:<8} (tol {:.0e})  {verdict}\n",
                    "",
                    reference.source,
                    reference.tolerance
                ));
            }
        }
        let sources = references
            .iter()
            .map(|reference| reference.source.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        if failures.is_empty() {
            text.push_str(&format!(
                "  PASS     {checked} reference coefficients agree ({sources})\n"
            ));
            print!("{text}");
            Ok(())
        } else {
            text.push_str(&format!(
                "  FAIL     {} of {checked} reference coefficients disagree: {}\n",
                failures.len(),
                failures.join("; ")
            ));
            print!("{text}");
            Err(format!("{}: {}", self.title, failures.join("; ")))
        }
    }
}

/// One master product of a result whose masters stay symbolic.
pub struct ProductRow {
    /// Short product label, e.g. `M1^2*M7`.
    pub label: String,
    /// Leading Laurent terms of vakint's coefficient of the product.
    pub vakint: String,
    /// What the coefficient was compared with, or `None` without a reference.
    pub reference: Option<String>,
}

/// The numbers in a product label (`M1^2*M7` -> `[1, 2, 7]`), for natural order.
fn label_key(label: &str) -> Vec<u64> {
    label
        .split(|character: char| !character.is_ascii_digit())
        .filter_map(|number| number.parse().ok())
        .collect()
}

/// Print a master-product table: header facts, one block per product (in
/// natural label order), the master legend and a final verdict line.
pub fn print_master_products(
    title: &str,
    facts: &[(String, String)],
    rows: &[ProductRow],
    legend: &[String],
    verdict: &str,
) {
    let mut text = format!("\n━━ {title} ━━\n");
    for (label, fact) in facts {
        text.push_str(&format!("  {label:<9} {fact}\n"));
    }
    let width = rows
        .iter()
        .map(|row| row.label.chars().count())
        .max()
        .unwrap_or(0)
        .max(7);
    text.push_str(&format!(
        "  {:>3}  {:<width$}  coefficient at d = 4-2ε (leading Laurent terms)\n",
        "#", "product"
    ));
    let mut ordered = rows.iter().collect::<Vec<_>>();
    ordered.sort_by_key(|row| label_key(&row.label));
    for (index, row) in ordered.into_iter().enumerate() {
        text.push_str(&format!(
            "  {:>3}  {:<width$}  vakint     {}\n",
            index + 1,
            row.label,
            row.vakint
        ));
        if let Some(reference) = &row.reference {
            text.push_str(&format!(
                "  {:>3}  {:<width$}  reference  {reference}\n",
                "", ""
            ));
        }
    }
    if !legend.is_empty() {
        text.push_str("  masters\n");
        for line in legend {
            text.push_str(&format!("    {line}\n"));
        }
    }
    text.push_str(&format!("  {verdict}\n"));
    print!("{text}");
}

/// A complex value from decimal real and imaginary parts.
pub fn complex_from_decimal(real: &str, imaginary: &str) -> Result<Complex<Float>, String> {
    Ok(Complex::new(parse(real)?, parse(imaginary)?))
}
