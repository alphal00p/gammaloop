//! Engine settings and external data shared by the RustRed integration tests.
//!
//! Include it with `#[path = "support/rustred.rs"] mod rustred_support;`.

use std::path::PathBuf;

use vakint::{
    EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions, RustRedMasterMode, VakintSettings,
};

// The one inclusion of `external_data.rs` in a test binary (`test_utils`
// reuses it): clippy rejects a file loaded as a module twice.
#[cfg(all(unix, not(target_arch = "wasm32")))]
#[path = "external_data.rs"]
pub mod external_data_fixture;

/// Settings that evaluate with RustRed alone, everything else at its default.
///
/// `form_guard` is the FORM executable path. Pass a path that does not exist
/// and names the test: a FORM call then fails loudly instead of silently
/// falling back to a FORM-based backend.
pub fn rustred_only(master_mode: RustRedMasterMode, form_guard: &str) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode },
        )]),
        form_exe_path: form_guard.into(),
        ..VakintSettings::default()
    }
}

/// The vakint-external-data checkout selected with `VAKINT_TEST_EXTERNAL_DATA`.
pub fn external_data_checkout() -> Option<PathBuf> {
    std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(PathBuf::from)
}

/// Load options for the external-data release that the vendored manifest
/// pins: read from the `VAKINT_TEST_EXTERNAL_DATA` checkout when it is set,
/// otherwise installed (or reused) from the published release.
#[cfg(all(unix, not(target_arch = "wasm32")))]
pub fn pinned_release() -> vakint::OfficialExternalDataOptions {
    vakint::OfficialExternalDataOptions {
        version: external_data_fixture::required_version(),
        local_artifacts_path: external_data_checkout(),
        ..Default::default()
    }
}
