use color_eyre::{config::HookBuilder, Result};

use crate::{model::UFOSymbol, utils::GS};

static INITIALISED: std::sync::Once = std::sync::Once::new();

pub fn initialise() -> Result<()> {
    INITIALISED.call_once(|| {
        let (panic, eyre) = HookBuilder::default()
            .capture_span_trace_by_default(cfg!(debug_assertions))
            .into_hooks();
        panic.install();
        eyre.install().unwrap();

        let _ = GS.delta_vec;
        let _ = UFOSymbol::zero();

        crate::set_interrupt_handler();
        crate::initialize_reps();
    });
    Ok(())
}

pub fn test_initialise() -> Result<()> {
    use crate::utils::tracing::init_test_tracing;

    init_test_tracing();
    initialise()?;

    Ok(())
}

pub fn bench_initialise() -> Result<()> {
    use crate::utils::tracing::init_bench_tracing;

    init_bench_tracing();
    initialise()?;

    Ok(())
}
