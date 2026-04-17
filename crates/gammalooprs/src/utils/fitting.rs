use crate::utils::{F, FloatLike};
use color_eyre::eyre::{Result, eyre};
use std::cmp::Ordering;
use symbolica::domains::float::Real;

pub fn constant_dropped_fit_points<T: FloatLike>(
    start: &F<T>,
    stop: &F<T>,
    num: usize,
) -> Result<Vec<F<T>>> {
    if num < 3 {
        return Err(eyre!(
            "need at least 3 points to generate a grid usable by the constant-dropped fit"
        ));
    }

    if !start.0.is_finite() || !stop.0.is_finite() || start <= &start.zero() || stop <= &stop.zero()
    {
        return Err(eyre!(
            "fit points must be generated from finite positive bounds"
        ));
    }

    if stop == start {
        return Err(eyre!(
            "fit point range must have distinct bounds for multiplicative spacing"
        ));
    }

    let is_ascending = stop > start;
    let log_start = start.log10();
    let log_stop = stop.log10();
    let step_count = start.from_usize(num - 1);
    let log_step = (&log_stop - &log_start) / &step_count;
    if log_step.is_nan() || log_step.is_infinite() || log_step.abs().less_than_epsilon() {
        return Err(eyre!(
            "fit point range must yield a finite non-zero log step"
        ));
    }

    let ten = F::from_f64(10.0);
    let mut points = Vec::with_capacity(num);
    points.push(start.clone());

    for index in 1..(num - 1) {
        let step_index = log_step.from_usize(index);
        let exponent = &log_start + &(&log_step * &step_index);
        let point = ten.powf(&exponent);
        if point.is_nan() || point.is_infinite() {
            return Err(eyre!("generated a non-finite fit point"));
        }
        let previous = points.last().expect("point list is non-empty");
        let advances = if is_ascending {
            &point > previous
        } else {
            &point < previous
        };
        if !advances {
            return Err(eyre!(
                "range and point count collapse at current precision; generated fit points are not strictly monotonic"
            ));
        }

        points.push(point);
    }

    let previous = points.last().expect("point list is non-empty");
    let advances = if is_ascending {
        stop > previous
    } else {
        stop < previous
    };
    if !advances {
        return Err(eyre!(
            "range and point count collapse at current precision; final bound does not preserve the requested monotonic direction"
        ));
    }

    points.push(stop.clone());
    Ok(points)
}

pub fn log_log_slope_constant_dropped<T: FloatLike>(x: &[F<T>], y: &[F<T>]) -> Result<SlopeFit<T>> {
    let samples = collect_valid_samples(x, y)?;
    if samples.len() < 3 {
        return Err(eyre!(
            "need at least 3 valid samples for constant-dropped slope fit"
        ));
    }

    let transformed_samples = collect_log_difference_samples(&samples)?;
    linear_regression(&transformed_samples)
}

pub struct SlopeFit<T: FloatLike> {
    pub slope: F<T>,
    pub r_suared: F<T>,
}

impl<T: FloatLike> SlopeFit<T> {
    pub fn r_squared(&self) -> &F<T> {
        &self.r_suared
    }
}

fn collect_valid_samples<T: FloatLike>(x: &[F<T>], y: &[F<T>]) -> Result<Vec<(F<T>, F<T>)>> {
    if x.len() != y.len() {
        return Err(eyre!(
            "x and y must have the same length, got {} and {}",
            x.len(),
            y.len()
        ));
    }

    let mut samples = x
        .iter()
        .zip(y)
        .filter_map(|(x_value, y_value)| {
            if !x_value.0.is_finite() || !y_value.0.is_finite() || x_value <= &x_value.zero() {
                return None;
            }

            Some((x_value.clone(), y_value.clone()))
        })
        .collect::<Vec<_>>();

    samples.sort_by(|(lhs, _), (rhs, _)| lhs.partial_cmp(rhs).unwrap_or(Ordering::Equal));

    let mut deduped = Vec::with_capacity(samples.len());
    for sample in samples {
        if deduped
            .last()
            .is_some_and(|(previous_x, _)| previous_x == &sample.0)
        {
            continue;
        }
        deduped.push(sample);
    }

    Ok(deduped)
}

fn collect_log_difference_samples<T: FloatLike>(
    samples: &[(F<T>, F<T>)],
) -> Result<Vec<(F<T>, F<T>)>> {
    let reference_log_x = samples[0].0.log10();
    let next_log_x = samples[1].0.log10();
    let reference_step = &next_log_x - &reference_log_x;
    if reference_step.is_nan()
        || reference_step.is_infinite()
        || reference_step.abs().less_than_epsilon()
    {
        return Err(eyre!("x grid must have a non-zero multiplicative spacing"));
    }

    let tolerance = {
        let scale = &reference_step.abs() + &reference_step.one();
        F::from_f64(1.0e-9) * scale
    };

    let mut transformed = Vec::with_capacity(samples.len().saturating_sub(1));
    for pair in samples.windows(2) {
        let left_log_x = pair[0].0.log10();
        let right_log_x = pair[1].0.log10();
        let step = &right_log_x - &left_log_x;
        if (&step - &reference_step).abs() > tolerance {
            return Err(eyre!(
                "x grid must have constant multiplicative spacing for constant-dropped log-log fitting"
            ));
        }

        let difference_magnitude = (&pair[1].1 - &pair[0].1).abs();
        if !difference_magnitude.positive()
            || difference_magnitude.is_nan()
            || difference_magnitude.is_infinite()
        {
            return Err(eyre!(
                "adjacent y differences must stay finite and non-zero after dropping the constant term"
            ));
        }

        let log_difference = difference_magnitude.log10();
        if log_difference.is_nan() || log_difference.is_infinite() {
            return Err(eyre!(
                "adjacent y differences must be positive in magnitude for the log transform"
            ));
        }

        transformed.push((left_log_x, log_difference));
    }

    if transformed.len() < 2 {
        return Err(eyre!(
            "need at least 2 adjacent log-difference samples for linear regression"
        ));
    }

    Ok(transformed)
}

fn linear_regression<T: FloatLike>(samples: &[(F<T>, F<T>)]) -> Result<SlopeFit<T>> {
    let zero = samples[0].0.zero();
    let one = zero.one();
    let n = zero.from_usize(samples.len());

    let mut sum_x = zero.clone();
    let mut sum_y = zero.clone();
    let mut sum_xy = zero.clone();
    let mut sum_x2 = zero.clone();
    for (x_value, y_value) in samples {
        sum_x += x_value.clone();
        sum_y += y_value.clone();
        sum_xy += x_value * y_value;
        sum_x2 += x_value * x_value;
    }

    let denominator = &(&n * &sum_x2) - &(&sum_x * &sum_x);
    if denominator.abs().less_than_epsilon() {
        return Err(eyre!(
            "log-log regression is singular for the transformed difference samples"
        ));
    }

    let slope = (&(&n * &sum_xy) - &(&sum_x * &sum_y)) / &denominator;
    if slope.is_nan() || slope.is_infinite() {
        return Err(eyre!("log-log regression produced a non-finite slope"));
    }

    let mean_y = &sum_y / &n;
    let mut ss_tot = zero.clone();
    let mut ss_res = zero;
    let intercept = (&sum_y - &(&slope * &sum_x)) / &n;
    for (x_value, y_value) in samples {
        let prediction = &(&slope * x_value) + &intercept;
        let centered = y_value - &mean_y;
        let residual = y_value - &prediction;

        ss_tot += &centered * &centered;
        ss_res += &residual * &residual;
    }

    let r_squared = if ss_tot.less_than_epsilon() {
        one
    } else {
        &one - &(&ss_res / &ss_tot)
    };

    if r_squared.is_nan() || r_squared.is_infinite() {
        Err(eyre!("log-log regression produced a non-finite r-squared"))
    } else {
        Ok(SlopeFit {
            slope,
            r_suared: r_squared,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ff64_values(values: &[f64]) -> Vec<F<f64>> {
        values.iter().copied().map(F::from_f64).collect()
    }

    #[test]
    fn constant_dropped_fit_recovers_power_law_slope() {
        let coefficient = 3.2_f64;
        let slope = -1.75_f64;
        let offset = 0.6_f64;
        let ratio: f64 = 1.6;

        let x = constant_dropped_fit_points(
            &F::from_f64(0.2_f64),
            &F::from_f64(0.2_f64 * ratio.powi(7)),
            8,
        )
        .expect("fit point generation should succeed");
        let x_values = x
            .iter()
            .map(|x_value| x_value.into_f64())
            .collect::<Vec<_>>();
        let y = x
            .iter()
            .map(|x_value| coefficient * x_value.into_f64().powf(slope) + offset)
            .collect::<Vec<_>>();

        let fit = log_log_slope_constant_dropped(&x, &ff64_values(&y))
            .expect("power-law fit should succeed");

        assert!((x_values[0] - 0.2_f64).abs() < 1.0e-12);
        assert!((x_values[7] - 0.2_f64 * ratio.powi(7)).abs() < 1.0e-12);
        assert!((fit.slope.into_f64() - slope).abs() < 1.0e-10);
        assert!(fit.r_squared().into_f64() > 0.999_999_999);
    }

    #[test]
    fn constant_dropped_fit_supports_descending_geometric_grid() {
        let coefficient = 3.2_f64;
        let slope = -1.75_f64;
        let offset = 0.6_f64;
        let ratio: f64 = 1.6;

        let x = constant_dropped_fit_points(
            &F::from_f64(0.2_f64 * ratio.powi(7)),
            &F::from_f64(0.2_f64),
            8,
        )
        .expect("descending fit point generation should succeed");
        let x_values = x
            .iter()
            .map(|x_value| x_value.into_f64())
            .collect::<Vec<_>>();
        let y = x
            .iter()
            .map(|x_value| coefficient * x_value.into_f64().powf(slope) + offset)
            .collect::<Vec<_>>();

        let fit = log_log_slope_constant_dropped(&x, &ff64_values(&y))
            .expect("power-law fit should succeed on a descending geometric grid");

        assert!((x_values[0] - 0.2_f64 * ratio.powi(7)).abs() < 1.0e-12);
        assert!((x_values[7] - 0.2_f64).abs() < 1.0e-12);
        assert!(x_values.windows(2).all(|window| window[0] > window[1]));
        assert!((fit.slope.into_f64() - slope).abs() < 1.0e-10);
        assert!(fit.r_squared().into_f64() > 0.999_999_999);
    }

    #[test]
    fn constant_dropped_fit_points_require_strict_positive_range() {
        let points = constant_dropped_fit_points(
            &F::<f64>::from_f64(1.0_f64),
            &F::<f64>::from_f64(1.0_f64),
            8,
        );

        assert!(points.is_err());
    }

    #[test]
    fn constant_dropped_fit_rejects_non_geometric_grid() {
        let coefficient = 3.2_f64;
        let slope = -1.75_f64;
        let offset = 0.6_f64;

        let x = [
            0.2_f64, 0.37_f64, 0.55_f64, 0.92_f64, 1.3_f64, 1.85_f64, 2.75_f64, 3.6_f64,
        ];
        let y = x
            .iter()
            .map(|x_value| coefficient * (*x_value).powf(slope) + offset)
            .collect::<Vec<_>>();

        let fit = log_log_slope_constant_dropped(&ff64_values(&x), &ff64_values(&y));

        assert!(fit.is_err());
    }
}
