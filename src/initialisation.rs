use color_eyre::{config::HookBuilder, Result};
use symbolica::activate_oem_license;

use crate::utils::GS;

static INITIALISED: std::sync::Once = std::sync::Once::new();

pub(crate) fn initialise() -> Result<()> {
    INITIALISED.call_once(|| {
        let (panic, eyre) = HookBuilder::default()
            .capture_span_trace_by_default(cfg!(debug_assertions))
            .into_hooks();
        panic.install();
        eyre.install().unwrap();

        if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
            activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
        };

        let _ = GS.delta_vec;
        crate::set_interrupt_handler();
        crate::initialize_reps();
    });
    Ok(())
}

pub fn test_initialise() -> Result<()> {
    use crate::cli::tracing::init_test_tracing;

    init_test_tracing();
    crate::initialize_reps();

    Ok(())
}

pub fn bench_initialise() -> Result<()> {
    use crate::cli::tracing::init_bench_tracing;

    init_bench_tracing();
    crate::initialize_reps();

    Ok(())
}
