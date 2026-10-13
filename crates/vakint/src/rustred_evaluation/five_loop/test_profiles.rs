//! Explicit test-only resource policy; never a production default or scope change.
use super::VakintSettings;

pub(super) const MATERIALIZATION_PROFILE: &str = "materialize-512m";

pub(super) fn apply_materialization_profile(settings: &mut VakintSettings) {
    let exact_algebra = settings
        .rustred_application_options
        .trace_limits
        .expansion
        .exact_algebra;
    let limits = &mut settings.rustred_application_options.trace_limits;
    limits.max_transport_operations = 512_000_000;
    limits.max_unique_nodes = 8_000_000;
    limits.max_transport_endpoints = 8_000_000;
    limits.expansion.max_native_polynomial_operations = 512_000_000;
    limits.expansion.max_native_exponent_entries = 128_000_000;
    limits.expansion.max_endpoint_power_entries = 128_000_000;
    limits.expansion.max_native_polynomial_terms = 8_000_000;
    limits.expansion.max_endpoints = 8_000_000;
    limits.expansion.max_retained_endpoint_key_bytes = 2_000_000_000;
    limits.expansion.max_retained_coefficient_clone_owned_bytes = 2_000_000_000;
    assert_eq!(limits.expansion.exact_algebra, exact_algebra);
}

#[test]
fn materialization_profile_preserves_other_policies() {
    let mut settings = VakintSettings::default();
    settings.rustred_scalar_numerator_limits.max_input_nodes = 17;
    settings
        .rustred_momentum_match_limits
        .max_verified_witnesses = 19;
    let original = settings.clone();
    apply_materialization_profile(&mut settings);
    let mut expected = original.rustred_application_options;
    expected.trace_limits = settings.rustred_application_options.trace_limits;
    assert_eq!(settings.rustred_application_options, expected);
    assert_eq!(
        settings.rustred_scalar_numerator_limits,
        original.rustred_scalar_numerator_limits
    );
    assert_eq!(
        settings.rustred_momentum_match_limits,
        original.rustred_momentum_match_limits
    );
    assert_eq!(
        settings.rustred_artifact_path,
        original.rustred_artifact_path
    );
    assert_eq!(
        settings.rustred_product_limits,
        original.rustred_product_limits
    );
}
