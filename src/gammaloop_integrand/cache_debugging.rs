//! Essential cache debugging utilities for GammaLoop
//!
//! This module provides simple debugging tools to verify cache behavior
//! and detect performance issues with polarization caching.

use crate::{gammaloop_integrand::GammaloopIntegrand, status_debug, status_info, status_warn};

/// Check if cache debug mode is enabled via environment variable
pub fn is_debug_cache_enabled() -> bool {
    std::env::var("GAMMALOOP_DEBUG_CACHE")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false)
}

/// Simple cache health check
pub fn quick_cache_check<I: GammaloopIntegrand>(integrand: &I, context: &str) -> bool {
    let validation = integrand.validate_cache_consistency();
    let stats = integrand.get_cache_stats();

    if !validation.is_valid {
        status_warn!(
            "❌ Cache health check failed at {}: {}",
            context,
            validation.diagnostics
        );
        return false;
    }

    if stats.efficiency_ratio < 0.3 {
        status_warn!(
            "⚠️ Very low cache efficiency at {}: {:.1}%",
            context,
            stats.efficiency_ratio * 100.0
        );
        return false;
    }

    if is_debug_cache_enabled() {
        status_info!(
            "✅ Cache health check passed at {}: {:.1}% efficiency",
            context,
            stats.efficiency_ratio * 100.0
        );
    }

    true
}

/// Validate cache state before critical operations
pub fn validate_before_operation<I: GammaloopIntegrand>(
    integrand: &I,
    operation: &str,
) -> Result<(), String> {
    let validation = integrand.validate_cache_consistency();
    if !validation.is_valid {
        return Err(format!(
            "Cache validation failed before {}: {}",
            operation, validation.diagnostics
        ));
    }
    Ok(())
}

/// Debug cache state with detailed output (only when debug enabled)
pub fn debug_cache_state<I: GammaloopIntegrand>(integrand: &I, context: &str) {
    if !is_debug_cache_enabled() {
        return;
    }

    let validation = integrand.validate_cache_consistency();
    let stats = integrand.get_cache_stats();

    status_debug!("🔍 CACHE STATE at {}: {}", context, validation);
    status_debug!("   Stats: {}", stats);

    if !validation.is_valid {
        status_warn!("   ❌ CACHE CORRUPTION DETECTED at {}", context);
    } else if validation.has_rotations {
        status_debug!(
            "   🔄 {} rotational variants",
            validation.current_external_cache_id - validation.base_external_cache_id
        );
    }

    if stats.efficiency_ratio < 0.3 {
        status_warn!(
            "   ⚠️ Very low efficiency: {:.1}%",
            stats.efficiency_ratio * 100.0
        );
    } else if stats.efficiency_ratio > 0.8 {
        status_debug!(
            "   ✅ Excellent efficiency: {:.1}%",
            stats.efficiency_ratio * 100.0
        );
    }
}

/// Monitor cache during Monte Carlo iterations
pub fn monte_carlo_cache_monitor<I: GammaloopIntegrand>(
    integrand: &I,
    iteration: usize,
    check_interval: usize,
) {
    if !is_debug_cache_enabled() || !iteration.is_multiple_of(check_interval) {
        return;
    }

    let stats = integrand.get_cache_stats();
    status_info!(
        "Iteration {}: Cache efficiency {:.1}%, {} base configs, {} rotations",
        iteration,
        stats.efficiency_ratio * 100.0,
        stats.base_configurations,
        stats.rotational_variants
    );
}

/// Warn if cache efficiency falls below threshold
pub fn warn_cache_efficiency<I: GammaloopIntegrand>(
    integrand: &I,
    min_efficiency: f64,
    context: String,
) {
    let stats = integrand.get_cache_stats();
    if stats.efficiency_ratio < min_efficiency {
        status_warn!(
            "⚠️ Cache efficiency {:.1}% below threshold {:.1}% at: {}",
            stats.efficiency_ratio * 100.0,
            min_efficiency * 100.0,
            context
        );
    }
}

/// Example usage for development
pub fn setup_debug_environment() {
    unsafe {
        std::env::set_var("GAMMALOOP_DEBUG_CACHE", "1");
    }
    status_info!("🔧 Cache debug mode enabled");
    status_info!("   Use debug_cache!(integrand, \"context\") for state checks");
    status_info!("   Use monitor_cache!(integrand, condition, \"context\") for monitoring");
    status_info!(
        "   Use warn_cache_efficiency!(integrand, 0.5, \"context\") for efficiency warnings"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_debug_mode_detection() {
        // Test environment variable detection
        unsafe {
            std::env::remove_var("GAMMALOOP_DEBUG_CACHE");
            assert!(!is_debug_cache_enabled());

            std::env::set_var("GAMMALOOP_DEBUG_CACHE", "1");
            assert!(is_debug_cache_enabled());

            std::env::set_var("GAMMALOOP_DEBUG_CACHE", "true");
            assert!(is_debug_cache_enabled());

            std::env::set_var("GAMMALOOP_DEBUG_CACHE", "0");
            assert!(!is_debug_cache_enabled());
        }
    }
}
