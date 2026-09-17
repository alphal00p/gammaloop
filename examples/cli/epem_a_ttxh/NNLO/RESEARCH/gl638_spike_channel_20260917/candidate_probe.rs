// Read-only raw-map diagnostic at the retained exact canonical GL638 spike.
use color_eyre::{Result, eyre::{ensure, eyre}};
use gammaloop_api::{StateLoadOption, state::ProcessRef};
use gammalooprs::{initialisation::initialise, integrands::process::{GraphTerm, ProcessIntegrand, ProcessIntegrandImpl, SamplingChannelBridge, SamplingChannelId}, settings::RuntimeSettings, utils::{ArbPrec, FloatLike, SamplingFloat, SamplingPrecision, f128}};
use gammalooprs::utils::serde_utils::SmartSerde;
use serde_json::{Value, json};
use std::{fs, path::Path};
fn source<T: FloatLike>(
    bridge: &SamplingChannelBridge<T>,
    cube: &[f64],
) -> Result<(Vec<ArbPrec>, f64)> {
    let draw = bridge.forward(
        SamplingChannelId(3),
        &cube
            .iter()
            .map(|x| T::from_f64_exact_binary(*x))
            .collect::<Vec<_>>(),
    )?;
    let raw = draw
        .raw_coordinates
        .iter()
        .map(|x| {
            let (lo, hi) = x.mpfr_enclosure(1000);
            ensure!(
                lo == hi,
                "source not exactly representable in diagnostic Arb lane"
            );
            Ok(ArbPrec::from(lo))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((raw, draw.selected_factor()?.into_f64()))
}

fn main() -> Result<()> {
    initialise()?;
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let card = Path::new(&args[0]);
    let output = Path::new(&args[1]);
    let base = Path::new("target/gl638_spike_inverse_diagnostics");
    let retained: Value = serde_json::from_str(&fs::read_to_string(base.join("inverse_replay.json"))?)?;
    let cube = vec![0.65326967347543075, 0.99071624997667884, 0.34598647260057075, 0.73874307614468926, 0.90145414024659398, 0.087144592922320935, 0.58803498486189254, 0.099566804759442573, 0.8048060051539222, 0.67392249894159018, 0.66470428853177543, 0.060665899033751615];
    let mut loaded = StateLoadOption::read_only(Path::new("examples/cli/epem_a_ttxh/NNLO/gammaloop_state_GL638_optimized_lmbs")).load()?;
    loaded.state.activate_loaded_integrand_backends(false)?;
    let pid = loaded.state.resolve_process_ref(Some(&ProcessRef::Name("epem_a_tth".into())))?;
    let model = &loaded.state.model;
    let runtime = loaded.state.process_list.get_integrand_mut(pid, "NNLO")?;
    *runtime.get_mut_settings() = RuntimeSettings::from_file(base.join("frozen_optimized/integrands/epem_a_tth@NNLO/settings.toml"), "retained source")?;
    runtime.warm_up(model)?;
    let raw = {
        let ProcessIntegrand::CrossSection(process) = runtime else { return Err(eyre!("cross section required")); };
        let setup = process.data.graph_terms[0].sampling_setup();
        match process.sampling_source_policy()?.0 {
            SamplingPrecision::Quad => source(setup.sampling_bridge::<f128>()?, &cube)?,
            SamplingPrecision::Fixed256 => source(setup.sampling_bridge::<SamplingFloat>()?, &cube)?,
            SamplingPrecision::Arb => source(setup.sampling_bridge::<ArbPrec>()?, &cube)?,
            SamplingPrecision::Double => unreachable!(),
        }.0
    };
    ensure!(json!(raw.iter().map(ToString::to_string).collect::<Vec<_>>()) == retained["common_raw_coordinates_decimal"], "canonical point changed");
    *runtime.get_mut_settings() = RuntimeSettings::from_file(card, "candidate raw map probe")?;
    runtime.warm_up(model)?;
    let ProcessIntegrand::CrossSection(process) = runtime else { unreachable!() };
    let policy = process.sampling_source_policy()?;
    process.prepare_sampling_precision::<ArbPrec>()?;
    let bridge = process.data.graph_terms[0].sampling_setup().sampling_bridge::<ArbPrec>()?;
    let mut report = json!({"runtime_card":card,"scope":"raw maps only; no trained grid weighting", "canonical_point_matches_retained_trace":true, "source_policy":format!("{policy:?}"),"channels":[]});
    let mut total = 0.0;
    for (id, channel) in bridge.channels().iter().enumerate() {
        eprintln!("inverse channel {id} {}", channel.name);
        let row = match bridge.inverse(SamplingChannelId(id), &raw) {
            Ok(Some(inverse)) => {
                let density = inverse.map.inverse_jacobian.into_f64();
                total += density;
                let coordinates = inverse.map.coordinates.iter().map(FloatLike::into_f64).collect::<Vec<_>>();
                let roundtrip = bridge.forward(SamplingChannelId(id), &coordinates.iter().map(|x| ArbPrec::from_f64_exact_binary(*x)).collect::<Vec<_>>());
                let roundtrip = match roundtrip {
                    Ok(forward) => json!({"max_abs_raw_error_GeV":forward.raw_coordinates.iter().zip(&raw).map(|(a,b)|(a.into_f64()-b.into_f64()).abs()).fold(0.0, f64::max), "native_diagnostics":forward.map.diagnostics}),
                    Err(error) => json!({"error":format!("{error:#}")}),
                };
                json!({"id":id,"name":channel.name,"support":true,"density":density,"partition":inverse.partition.weight(id).unwrap().into_f64(),"J_times_partition":inverse.selected_factor()?.into_f64(),"coordinates":coordinates,"roundtrip":roundtrip,"diagnostics":inverse.map.diagnostics})
            },
            Ok(None) => json!({"id":id,"name":channel.name,"support":false}),
            Err(error) => json!({"id":id,"name":channel.name,"error":format!("{error:#}")}),
        };
        report["channels"].as_array_mut().unwrap().push(row);
        fs::write(output, serde_json::to_vec_pretty(&report)?)?;
    }
    report["raw_density_sum"] = json!(total);
    report["raw_density_ratio_to_optimized"] = json!(total/retained["modes"][0]["raw_map_density_sum"].as_f64().unwrap());
    report["raw_density_ratio_to_augmented"] = json!(total/retained["modes"][1]["raw_map_density_sum"].as_f64().unwrap());
    fs::write(output, serde_json::to_vec_pretty(&report)?)?;
    Ok(())
}
