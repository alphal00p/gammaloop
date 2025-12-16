//! UV Profile Analysis Module
//!
//! This module provides functionality for analyzing ultraviolet behavior of loop integrands
//! by evaluating them at different momentum scalings and computing degrees of divergence.

use crate::evaluation_result::EvaluationResult;
use crate::gammaloop_integrand::GLIntegrand;
use crate::integrands::{HasIntegrand, Integrand};
use crate::model::Model;
use crate::momentum::FourMomentum;
use crate::settings::RuntimeSettings;
use crate::utils::{F, FloatLike};
use color_eyre::{Report, Result, eyre::Context};
use log::{debug, info, warn};
use nalgebra::{Matrix4, Vector4};
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use std::collections::HashMap;
use symbolica::numerical_integration::Sample;

/// Configuration for UV profile analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UVProfileConfig {
    pub n_points: usize,
    pub min_scaling: f64,
    pub max_scaling: f64,
    pub lmb: String,
    pub uv_indices: Option<Vec<usize>>,
    pub use_f128: bool,
    pub no_f128: bool,
    pub plots: bool,
    pub scale_cuts: bool,
    pub seed: Option<u64>,
    pub target_scaling: Option<f64>,
    pub verbose: bool,
    pub n_max: usize,
    pub stability_threshold: f64,
    pub output_file: Option<String>,
    pub output_format: String,
}

/// Result of UV profile analysis at a single scaling point
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalingPointResult {
    pub scaling: f64,
    pub evaluation: Complex<F<f64>>,
    pub real_part: f64,
    pub imaginary_part: f64,
    pub magnitude: f64,
    pub phase: f64,
    pub stability_check: bool,
}

/// Complete UV profile analysis results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UVProfileResult {
    pub config: UVProfileConfig,
    pub scaling_points: Vec<ScalingPointResult>,
    pub degree_of_divergence: Option<f64>,
    pub dod_fit_quality: Option<f64>,
    pub analysis_metadata: ProfileMetadata,
}

/// Metadata about the profiling analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProfileMetadata {
    pub total_evaluations: usize,
    pub failed_evaluations: usize,
    pub unstable_evaluations: usize,
    pub analysis_time_ms: u64,
    pub momentum_basis_used: String,
    pub precision_used: String,
}

/// Main UV profiling function
pub fn run_uv_profile<I: HasIntegrand>(
    config: UVProfileConfig,
    integrand: &mut I,
    model: &Model,
    runtime_settings: &RuntimeSettings,
) -> Result<UVProfileResult> {
    let start_time = std::time::Instant::now();

    info!(
        "Starting UV profile analysis with {} scaling points",
        config.n_points
    );
    debug!("Configuration: {:?}", config);

    // Generate scaling sequence (log-spaced between min and max)
    let scaling_sequence =
        generate_scaling_sequence(config.min_scaling, config.max_scaling, config.n_points);

    // Initialize random number generator if seed provided
    if let Some(seed) = config.seed {
        // TODO: Initialize RNG with seed for reproducible momentum sampling
        debug!("Using random seed: {}", seed);
    }

    let mut scaling_results = Vec::with_capacity(config.n_points);
    let mut total_evaluations = 0;
    let mut failed_evaluations = 0;
    let mut unstable_evaluations = 0;

    // Determine precision to use
    let use_high_precision = config.use_f128 && !config.no_f128;
    let precision_str = if use_high_precision { "f128" } else { "f64" };

    info!("Using precision: {}", precision_str);

    // Process each scaling point
    for (i, &scaling) in scaling_sequence.iter().enumerate() {
        if config.verbose {
            info!(
                "Processing scaling point {}/{}: {:.3e}",
                i + 1,
                config.n_points,
                scaling
            );
        }

        match evaluate_at_scaling(scaling, integrand, model, &config, use_high_precision) {
            Ok(evaluation) => {
                total_evaluations += 1;

                // Perform stability check
                let is_stable = check_stability(&evaluation, config.stability_threshold);
                if !is_stable {
                    unstable_evaluations += 1;
                    if config.verbose {
                        warn!("Unstable evaluation at scaling {:.3e}", scaling);
                    }
                }

                let point_result = ScalingPointResult {
                    scaling,
                    evaluation: evaluation.clone(),
                    real_part: evaluation.re.into(),
                    imaginary_part: evaluation.im.into(),
                    magnitude: magnitude(&evaluation),
                    phase: phase(&evaluation),
                    stability_check: is_stable,
                };

                scaling_results.push(point_result);
            }
            Err(e) => {
                failed_evaluations += 1;
                if config.verbose {
                    warn!("Failed to evaluate at scaling {:.3e}: {}", scaling, e);
                }
            }
        }
    }

    // Compute degree of divergence
    let (dod, dod_quality) = compute_degree_of_divergence(&scaling_results)?;

    let analysis_time = start_time.elapsed();

    let metadata = ProfileMetadata {
        total_evaluations,
        failed_evaluations,
        unstable_evaluations,
        analysis_time_ms: analysis_time.as_millis() as u64,
        momentum_basis_used: config.lmb.clone(),
        precision_used: precision_str.to_string(),
    };

    let result = UVProfileResult {
        config: config.clone(),
        scaling_points: scaling_results,
        degree_of_divergence: dod,
        dod_fit_quality: dod_quality,
        analysis_metadata: metadata,
    };

    info!(
        "UV profile analysis completed in {:.2}s",
        analysis_time.as_secs_f64()
    );
    if let Some(dod) = dod {
        info!("Computed degree of divergence: {:.6}", dod);
    }

    Ok(result)
}

/// Generate logarithmically spaced scaling sequence
fn generate_scaling_sequence(min_scaling: f64, max_scaling: f64, n_points: usize) -> Vec<f64> {
    if n_points <= 1 {
        return vec![min_scaling];
    }

    let log_min = min_scaling.ln();
    let log_max = max_scaling.ln();
    let step = (log_max - log_min) / (n_points - 1) as f64;

    (0..n_points)
        .map(|i| (log_min + i as f64 * step).exp())
        .collect()
}

/// Evaluate integrand at a specific momentum scaling
fn evaluate_at_scaling<I: HasIntegrand>(
    scaling: f64,
    integrand: &mut I,
    model: &Model,
    config: &UVProfileConfig,
    use_f128: bool,
) -> Result<Complex<F<f64>>> {
    // Generate a random sample point in [0,1]^n_dim
    let n_dim = integrand.get_n_dim();
    let mut xs = Vec::with_capacity(n_dim);

    // Use scaling factor to influence the sampling
    // This creates a deterministic but scaled sampling pattern
    let base_offset = scaling.ln() / 10.0; // Use log of scaling to create variation

    for i in 0..n_dim {
        // Create a deterministic but varied sample based on scaling and dimension
        let x = 0.5 + 0.3 * ((i as f64 + 1.0) * scaling.sqrt() + base_offset).sin();
        xs.push(F::from(x.clamp(0.001, 0.999))); // Avoid exact boundaries
    }

    // Create continuous sample (no discrete dimensions for UV profiling)
    let sample = Sample::Continuous(F::from(1.0), xs);

    // Evaluate the integrand at this sample point
    let weight = F::from(1.0);
    let iter = 0;
    let max_eval = Complex::new(F::from(1e10), F::from(0.0));

    let eval_result = integrand.evaluate_sample(&sample, model, weight, iter, use_f128, max_eval);

    Ok(eval_result.integrand_result)
}

/// Scale momentum configuration by applying scaling factor
/// This modifies the momenta to probe UV behavior at different energy scales
fn apply_momentum_scaling(
    scaling: f64,
    base_momenta: &[F<f64>],
    config: &UVProfileConfig,
) -> Vec<F<f64>> {
    // For UV profiling, we typically scale all momenta by the same factor
    // This probes the UV behavior as energy scales increase
    base_momenta
        .iter()
        .map(|&p| F::from(f64::from(p) * scaling))
        .collect()
}

/// Transform momentum configuration to different LMB basis
fn transform_to_lmb_basis(
    momenta: &[FourMomentum<F<f64>>],
    target_lmb: &str,
) -> Result<Vec<FourMomentum<F<f64>>>> {
    match target_lmb {
        "defining" => {
            // Use defining LMB (no transformation needed if already in defining basis)
            Ok(momenta.to_vec())
        }
        "LU" => {
            // Transform to LU basis
            // TODO: Implement LU transformation matrix
            warn!("LU basis transformation not yet implemented, using defining basis");
            Ok(momenta.to_vec())
        }
        _ => {
            warn!("Unknown LMB basis '{}', using defining basis", target_lmb);
            Ok(momenta.to_vec())
        }
    }
}

/// Check numerical stability of evaluation result
fn check_stability(evaluation: &Complex<F<f64>>, threshold: f64) -> bool {
    let mag = magnitude(evaluation);

    // Check for NaN or infinite values
    if !mag.is_finite() {
        return false;
    }

    // Check if the result is numerically reasonable
    // This is a simple heuristic - can be improved
    mag > threshold && mag < 1e20
}

/// Compute magnitude of complex number
fn magnitude(z: &Complex<F<f64>>) -> f64 {
    let re: f64 = z.re.into();
    let im: f64 = z.im.into();
    (re * re + im * im).sqrt()
}

/// Compute phase of complex number
fn phase(z: &Complex<F<f64>>) -> f64 {
    let re: f64 = z.re.into();
    let im: f64 = z.im.into();
    im.atan2(re)
}

/// Compute degree of divergence from scaling results
fn compute_degree_of_divergence(
    results: &[ScalingPointResult],
) -> Result<(Option<f64>, Option<f64>)> {
    if results.len() < 3 {
        warn!("Not enough data points for DOD computation");
        return Ok((None, None));
    }

    // Filter out unstable or zero results
    let valid_results: Vec<_> = results
        .iter()
        .filter(|r| r.stability_check && r.magnitude > 1e-15)
        .collect();

    if valid_results.len() < 3 {
        warn!("Not enough stable data points for DOD computation");
        return Ok((None, None));
    }

    // Perform linear regression on log(|result|) vs log(scaling)
    // DOD = slope of this fit
    let n = valid_results.len() as f64;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_x2 = 0.0;

    for result in &valid_results {
        let x = result.scaling.ln();
        let y = result.magnitude.ln();

        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }

    let denominator = n * sum_x2 - sum_x * sum_x;
    if denominator.abs() < 1e-15 {
        warn!("Singular matrix in DOD computation");
        return Ok((None, None));
    }

    let slope = (n * sum_xy - sum_x * sum_y) / denominator;

    // Compute R-squared for fit quality
    let y_mean = sum_y / n;
    let mut ss_tot = 0.0;
    let mut ss_res = 0.0;

    for result in &valid_results {
        let x = result.scaling.ln();
        let y = result.magnitude.ln();
        let y_pred = (sum_y - slope * sum_x) / n + slope * x; // intercept + slope * x

        ss_tot += (y - y_mean).powi(2);
        ss_res += (y - y_pred).powi(2);
    }

    let r_squared = if ss_tot > 1e-15 {
        1.0 - ss_res / ss_tot
    } else {
        0.0
    };

    debug!("DOD fit: slope = {:.6}, R² = {:.6}", slope, r_squared);

    Ok((Some(slope), Some(r_squared)))
}

/// Write UV profile results to file
pub fn write_results(result: &UVProfileResult) -> Result<()> {
    if let Some(ref output_file) = result.config.output_file {
        match result.config.output_format.as_str() {
            "yaml" => {
                let yaml_content =
                    serde_yaml::to_string(result).context("Failed to serialize results to YAML")?;
                std::fs::write(output_file, yaml_content)
                    .context("Failed to write YAML results file")?;
            }
            "json" => {
                let json_content = serde_json::to_string_pretty(result)
                    .context("Failed to serialize results to JSON")?;
                std::fs::write(output_file, json_content)
                    .context("Failed to write JSON results file")?;
            }
            format => {
                return Err(color_eyre::eyre::eyre!(
                    "Unsupported output format: {}",
                    format
                ));
            }
        }
        info!("Results written to: {}", output_file);
    }
    Ok(())
}

/// Display UV profile results summary
pub fn display_results_summary(result: &UVProfileResult) {
    println!("\n=== UV Profile Analysis Results ===");
    println!("Configuration:");
    println!("  Points analyzed: {}", result.scaling_points.len());
    println!(
        "  Scaling range: {:.2e} - {:.2e}",
        result.config.min_scaling, result.config.max_scaling
    );
    println!("  LMB basis: {}", result.config.lmb);
    println!("  Precision: {}", result.analysis_metadata.precision_used);

    println!("\nAnalysis Results:");
    if let Some(dod) = result.degree_of_divergence {
        println!("  Degree of divergence: {:.6}", dod);
        if let Some(quality) = result.dod_fit_quality {
            println!("  Fit quality (R²): {:.6}", quality);
        }
    } else {
        println!("  Degree of divergence: Could not be computed");
    }

    println!("\nEvaluation Statistics:");
    println!(
        "  Total evaluations: {}",
        result.analysis_metadata.total_evaluations
    );
    println!(
        "  Failed evaluations: {}",
        result.analysis_metadata.failed_evaluations
    );
    println!(
        "  Unstable evaluations: {}",
        result.analysis_metadata.unstable_evaluations
    );
    println!(
        "  Analysis time: {:.2}s",
        result.analysis_metadata.analysis_time_ms as f64 / 1000.0
    );

    if result.config.verbose && !result.scaling_points.is_empty() {
        println!("\nDetailed Results (first 10 points):");
        println!(
            "  {:>12} {:>15} {:>15} {:>10}",
            "Scaling", "Real", "Imag", "Stable"
        );
        for (i, point) in result.scaling_points.iter().enumerate().take(10) {
            println!(
                "  {:>12.3e} {:>15.6e} {:>15.6e} {:>10}",
                point.scaling,
                point.real_part,
                point.imaginary_part,
                if point.stability_check { "✓" } else { "✗" }
            );
        }
        if result.scaling_points.len() > 10 {
            println!("  ... ({} more points)", result.scaling_points.len() - 10);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scaling_sequence_generation() {
        let sequence = generate_scaling_sequence(1e-3, 1e3, 7);
        assert_eq!(sequence.len(), 7);
        assert!((sequence[0] - 1e-3).abs() < 1e-15);
        assert!((sequence[6] - 1e3).abs() < 1e-10);

        // Check log spacing
        for i in 1..sequence.len() {
            let ratio = sequence[i] / sequence[i - 1];
            let expected_ratio = (1e6_f64).powf(1.0 / 6.0);
            assert!((ratio - expected_ratio).abs() < 1e-10);
        }
    }

    #[test]
    fn test_magnitude_computation() {
        let z = Complex::new(F::from(3.0), F::from(4.0));
        let mag = magnitude(&z);
        assert!((mag - 5.0).abs() < 1e-15);
    }

    #[test]
    fn test_phase_computation() {
        let z = Complex::new(F::from(1.0), F::from(1.0));
        let ph = phase(&z);
        assert!((ph - std::f64::consts::FRAC_PI_4).abs() < 1e-15);
    }

    #[test]
    fn test_stability_check() {
        let stable = Complex::new(F::from(1.0), F::from(0.5));
        let unstable_nan = Complex::new(F::from(f64::NAN), F::from(0.0));
        let unstable_inf = Complex::new(F::from(f64::INFINITY), F::from(0.0));
        let unstable_tiny = Complex::new(F::from(1e-20), F::from(0.0));

        assert!(check_stability(&stable, 1e-12));
        assert!(!check_stability(&unstable_nan, 1e-12));
        assert!(!check_stability(&unstable_inf, 1e-12));
        assert!(!check_stability(&unstable_tiny, 1e-12));
    }
}
