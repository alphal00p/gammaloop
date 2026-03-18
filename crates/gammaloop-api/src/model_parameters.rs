use color_eyre::Result;
use eyre::{eyre, Context};
use gammalooprs::{
    model::{Model, ParameterNature, ParameterType},
    utils::F,
};
use spenso::algebra::complex::Complex;

pub(crate) const MODEL_COMPLEX_VALUE_FORMAT_HINT: &str =
    "expects complex like 1.0, -1.0e13, 1.0+2.0i, or 1.0-2.0j";
pub(crate) const MODEL_REAL_VALUE_FORMAT_HINT: &str = "expects real like 1.0, -1.0, or 1.0e13";

pub(crate) fn model_value_format_hint(parameter_type: Option<ParameterType>) -> &'static str {
    match parameter_type {
        Some(ParameterType::Real) => MODEL_REAL_VALUE_FORMAT_HINT,
        _ => MODEL_COMPLEX_VALUE_FORMAT_HINT,
    }
}

pub(crate) fn external_model_parameter_type(
    model: &Model,
    parameter_name: &str,
) -> Option<ParameterType> {
    model
        .get_parameter_opt(parameter_name)
        .filter(|parameter| parameter.nature == ParameterNature::External)
        .map(|parameter| parameter.parameter_type.clone())
}

pub(crate) fn parse_model_parameter_value(value: &str) -> Result<Complex<F<f64>>> {
    let value = value.trim();
    if value.is_empty() {
        return Err(eyre!("Model parameter value cannot be empty"));
    }
    if value.starts_with('{') || value.starts_with('[') {
        return Err(eyre!(
            "Legacy model-parameter syntax like '{{re:...,im:...}}' or '[re,im]' is deprecated; use {MODEL_COMPLEX_VALUE_FORMAT_HINT}"
        ));
    }

    let (re, im) = parse_complex_literal(value)?;
    Ok(Complex::new(F(re), F(im)))
}

fn parse_complex_literal(value: &str) -> Result<(f64, f64)> {
    let value = value.trim();
    if value.is_empty() {
        return Err(eyre!("Complex literal cannot be empty"));
    }

    if let Some(unit) = value.chars().last() {
        if matches!(unit, 'i' | 'j' | 'I' | 'J') {
            let body = &value[..value.len() - unit.len_utf8()];
            let body = body.trim();
            if body.is_empty() {
                return Err(eyre!(
                    "Missing coefficient before imaginary unit in '{value}'"
                ));
            }

            if let Some(index) = find_imaginary_separator(body) {
                let re = parse_float_literal(&body[..index], value)?;
                let im = parse_imaginary_literal(&body[index..], value)?;
                return Ok((re, im));
            }

            let im = parse_imaginary_literal(body, value)?;
            return Ok((0.0, im));
        }
    }

    Ok((parse_float_literal(value, value)?, 0.0))
}

fn find_imaginary_separator(value: &str) -> Option<usize> {
    let bytes = value.as_bytes();
    for index in (1..bytes.len()).rev() {
        let current = bytes[index] as char;
        if !matches!(current, '+' | '-') {
            continue;
        }
        let previous = bytes[index - 1] as char;
        if matches!(previous, 'e' | 'E') {
            continue;
        }
        return Some(index);
    }
    None
}

fn parse_float_literal(value: &str, full_value: &str) -> Result<f64> {
    value
        .trim()
        .parse::<f64>()
        .with_context(|| format!("Failed to parse real part of complex literal '{full_value}'"))
}

fn parse_imaginary_literal(value: &str, full_value: &str) -> Result<f64> {
    let value = value.trim();
    if value == "+" {
        return Ok(1.0);
    }
    if value == "-" {
        return Ok(-1.0);
    }
    value.parse::<f64>().with_context(|| {
        format!("Failed to parse imaginary part of complex literal '{full_value}'")
    })
}

pub(crate) fn validate_model_parameter_type(
    parameter_name: &str,
    parameter_type: ParameterType,
    value: &Complex<F<f64>>,
) -> Result<()> {
    if parameter_type == ParameterType::Real && value.im.0 != 0.0 {
        return Err(eyre!(
            "Model parameter '{parameter_name}' is Real and cannot be assigned an imaginary component; {}",
            MODEL_REAL_VALUE_FORMAT_HINT
        ));
    }
    Ok(())
}
