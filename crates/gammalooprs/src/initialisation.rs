use color_eyre::{Result, config::HookBuilder};
use spenso::network::library::function_lib::INBUILTS;
use spenso::network::parsing::SPENSO_TAG;

use crate::utils::init_vakint;
use crate::{model::UFOSymbol, numerator::ufo::UFO, utils::GS};
static INITIALISED: std::sync::Once = std::sync::Once::new();

pub fn initialise() -> Result<()> {
    INITIALISED.call_once(|| {
        // Here it would not match the crate name
        // if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
        //     activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
        // };

        let (panic, eyre) = HookBuilder::default()
            .capture_span_trace_by_default(cfg!(debug_assertions))
            .into_hooks();
        // println!("Installing panic and eyre hooks");
        panic.install();
        // Tests and embedded entry points may have already installed an eyre hook.
        let _ = eyre.install();

        // println!("Initializing symbols");
        let _ = GS.delta_vec;
        let _ = INBUILTS.conj;
        let _ = SPENSO_TAG.tag;
        let _ = UFO.complexconjugate;

        // let _ = Symbol::id();
        let _ = UFOSymbol::zero();

        // println!("Setting up interrupt handler");
        crate::set_interrupt_handler();
        // println!("Initialize_reps");
        crate::initialize_reps();
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
