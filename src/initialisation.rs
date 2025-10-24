use color_eyre::{config::HookBuilder, Result};
use symbolica::activate_oem_license;
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor, MATADOptions,
    PySecDecOptions, Vakint, VakintSettings,
};

use crate::utils::VAKINT;
use crate::{model::UFOSymbol, numerator::ufo::UFO, settings::GlobalSettings, utils::GS};
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
        panic.install();
        eyre.install().unwrap();

        let _ = GS.delta_vec;
        let _ = UFO.complexconjugate;
        let _ = UFOSymbol::zero();

        crate::set_interrupt_handler();
        crate::initialize_reps();
    });
    Ok(())
}

pub fn initialise_with_settings(settings: Option<&GlobalSettings>) -> Result<()> {
    let global_settings = if let Some(s) = settings {
        s
    } else {
        &GlobalSettings::default()
    };
    let vakint_settings = &global_settings.generation.uv.vakint;

    let mut vakint_evaluation_methods = EvaluationOrder::empty();
    for method in &vakint_settings.evaluation_methods {
        match method.as_str() {
            "alphaloop" => vakint_evaluation_methods
                .0
                .push(EvaluationMethod::AlphaLoop),
            "matad" => vakint_evaluation_methods
                .0
                .push(EvaluationMethod::MATAD(MATADOptions::default())),
            "fmft" => vakint_evaluation_methods
                .0
                .push(EvaluationMethod::FMFT(FMFTOptions::default())),
            "pysecdec" => vakint_evaluation_methods
                .0
                .push(EvaluationMethod::PySecDec(PySecDecOptions::default())),
            other => {
                return Err(color_eyre::eyre::eyre!(
                    "Unknown Vakint evaluation method: {}",
                    other
                ))
            }
        }
    }

    *VAKINT.write().unwrap() = Some(Vakint::new(Some(VakintSettings {
        allow_unknown_integrals: false,
        evaluation_order: vakint_evaluation_methods,
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        run_time_decimal_precision: vakint_settings.run_time_decimal_precision as u32,
        number_of_terms_in_epsilon_expansion: 5,
        temporary_directory: vakint_settings.temporary_directory.clone(),
        mu_r_sq_symbol: GS.mu_r_sq.get_name().to_string(),
        form_exe_path: vakint_settings.form_exe_path.clone(),
        python_exe_path: vakint_settings.python_exe_path.clone(),
        clean_tmp_dir: vakint_settings.clean_tmp_dir,
        ..VakintSettings::default()
    }))?);

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
