//! Explicit release installation or local-checkout loading, followed by a
//! numerical 1L evaluation. Requesting N always loads 1L through N.
//! cargo run --release -p vakint --example official_external_data -- MAX_LOOPS [PATH] [VERSION]

#[cfg(all(unix, not(target_arch = "wasm32")))]
fn main() -> Result<(), Box<dyn std::error::Error>> {
    use vakint::{
        EvaluationOrder, LoopNormalizationFactor, OfficialExternalDataOptions, Vakint,
        VakintSettings, vakint_parse,
    };

    let mut args = std::env::args().skip(1);
    let maximum: usize = args.next().ok_or("missing maximum loop count")?.parse()?;
    let path = args.next().map(Into::into);
    let version = args.next().unwrap_or_else(|| "latest".into());
    if args.next().is_some() {
        return Err("usage: official_external_data MAX_LOOPS [PATH] [VERSION]".into());
    }
    let mut engine = Vakint::new()?;
    let report = engine.load_official_external_data(
        maximum,
        &OfficialExternalDataOptions {
            version,
            local_artifacts_path: path,
            ..Default::default()
        },
    )?;
    println!("{}", serde_json::to_string_pretty(&report)?);
    let settings = VakintSettings {
        evaluation_order: EvaluationOrder::rustred_only(),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        form_exe_path: "/official-data-example-does-not-use-form".into(),
        ..Default::default()
    };
    let input = vakint_parse!("topo(I1L(1,1))")?;
    println!("{}", engine.evaluate_integral(&settings, input.as_view())?);
    Ok(())
}

#[cfg(not(all(unix, not(target_arch = "wasm32"))))]
fn main() {
    eprintln!(
        "The official download transport currently requires Linux or macOS; local package loading remains available."
    );
}
