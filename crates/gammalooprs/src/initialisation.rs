use color_eyre::{Result, config::HookBuilder};
use spenso::symbolica_init::in_symbolica_initializer;
use symbolica::{activate_oem_license, initialize};

use crate::numerator::ufo::UFO;
use crate::utils::{GS, init_vakint};
static INITIALISED: std::sync::Once = std::sync::Once::new();

initialize!(|| {
    in_symbolica_initializer(|| {
        let _ = GS.force_in_initializer();
        let _ = UFO.force_in_initializer();
    });
});

pub fn initialise() -> Result<()> {
    INITIALISED.call_once(|| {
        if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
            activate_oem_license!("SYMBOLICA_OEM_KEY_ba2512eb");
        };

        let (panic, eyre) = HookBuilder::default()
            .capture_span_trace_by_default(cfg!(debug_assertions))
            .into_hooks();
        // println!("Installing panic and eyre hooks");
        panic.install();
        // Tests and embedded entry points may have already installed an eyre hook.
        let _ = eyre.install();

        // println!("Setting up interrupt handler");
        crate::set_interrupt_handler();
        // println!("Initialize_reps");
    });
    // println!("Initializing Vakint");
    init_vakint()?;
    Ok(())
}

pub fn test_initialise() -> Result<()> {
    use crate::utils::tracing::init_test_tracing;

    init_test_tracing();
    initialise()?;
    init_vakint()?;

    Ok(())
}

pub fn bench_initialise() -> Result<()> {
    use crate::utils::tracing::init_bench_tracing;

    init_bench_tracing();
    initialise()?;
    init_vakint()?;

    Ok(())
}
