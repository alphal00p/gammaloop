//! Settings shared by the p-integral integration tests (`docs/p_integrals.md`, section 7).
//!
//! Include it with `#[path = "support/p_integrals.rs"] mod p_integrals;`.

use vakint::{
    EvaluationMethod, EvaluationOrder, ForcerOptions, RustRedMasterMode, VakintDependency,
    VakintSettings,
};

/// FORM 5 for live Forcer runs: `VAKINT_TEST_FORM5`, otherwise the default
/// FORM (`FORM_PATH` or `form`) when it is FORM 5, as in the dev shell and CI.
pub fn forcer_form() -> Option<String> {
    std::env::var("VAKINT_TEST_FORM5").ok().or_else(|| {
        let form = std::env::var("FORM_PATH").unwrap_or_else(|_| "form".into());
        VakintDependency::Forcer
            .version(&form, "python3")
            .is_ok()
            .then_some(form)
    })
}

/// Settings that evaluate p-integrals with Forcer alone, in `master_mode`.
pub fn forcer_settings(form: &str, master_mode: RustRedMasterMode) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
            master_mode,
        })]),
        form_exe_path: form.into(),
        n_cores: 1,
        ..VakintSettings::default()
    }
}
