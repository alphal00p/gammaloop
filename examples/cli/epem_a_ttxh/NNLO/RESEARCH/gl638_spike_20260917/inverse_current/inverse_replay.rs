// Read-only diagnostic client: existing maps, production evaluation, and Havana
// Grid::probe own all numerical operations. No live state is saved or trained.
use color_eyre::{
    Result,
    eyre::{ensure, eyre},
};
use gammaloop_api::{StateLoadOption, state::ProcessRef};
use gammalooprs::{
    initialisation::initialise,
    integrands::{
        evaluation::{GenericEvaluationResult, PreciseEvaluationResult},
        process::{
            GraphTerm, MomentumSpaceEvaluationInput, ProcessIntegrand, ProcessIntegrandImpl,
            SamplingChannelBridge, SamplingChannelId,
        },
    },
    momentum::ThreeMomentum,
    settings::RuntimeSettings,
    utils::{
        ArbPrec, F, FloatLike, SamplingFloat, SamplingPrecision, f128, serde_utils::SmartSerde,
    },
};
use serde_json::{Value, json};
use spenso::algebra::complex::Complex;
use std::{fs, path::Path};
use symbolica::numerical_integration::{Grid, Probe, Sample};

fn native<T: FloatLike>(r: GenericEvaluationResult<T>, precision: &str) -> Value {
    let remaining = r
        .parameterization_jacobian
        .clone()
        .unwrap_or_else(|| F(r.integrator_weight.0.one()))
        * &r.integrator_weight;
    json!({"precision":precision, "integrand_result":[r.integrand_result.re.to_string(),r.integrand_result.im.to_string()],
        "complete_estimator":[(&r.integrand_result.re * &remaining).0.into_f64(),(&r.integrand_result.im * &remaining).0.into_f64()],
        "parameterization_jacobian":r.parameterization_jacobian.map(|x|x.to_string()),
        "integrator_weight":r.integrator_weight.to_string(),"metadata_debug":format!("{:?}",r.evaluation_metadata),"metadata":r.evaluation_metadata})
}
fn evaluation(r: Result<PreciseEvaluationResult>) -> Value {
    match r {
        Ok(PreciseEvaluationResult::Double(r)) => native(r, "Double"),
        Ok(PreciseEvaluationResult::Quad(r)) => native(r, "Quad"),
        Ok(PreciseEvaluationResult::Arb(r)) => native(r, "Arb"),
        Err(e) => json!({"error":format!("{e:#}")}),
    }
}
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
    let base = Path::new(
        args.first()
            .map(String::as_str)
            .unwrap_or("target/gl638_spike_inverse_diagnostics"),
    );
    let nnlo = Path::new(
        args.get(1)
            .map(String::as_str)
            .unwrap_or("examples/cli/epem_a_ttxh/NNLO"),
    );
    let retained: Value =
        serde_json::from_str(&fs::read_to_string(base.join("frozen_inputs.json"))?)?;
    let cube = vec![
        0.65326967347543075,
        0.99071624997667884,
        0.34598647260057075,
        0.73874307614468926,
        0.90145414024659398,
        0.087144592922320935,
        0.58803498486189254,
        0.099566804759442573,
        0.8048060051539222,
        0.67392249894159018,
        0.66470428853177543,
        0.060665899033751615,
    ];
    let mut common_raw = Vec::<ArbPrec>::new();
    let mut physical = [0.0; 2];
    let mut report = json!({"scope":"same canonical physical point; production raw-map partition retained, frozen next-iteration adaptive densities probed", "modes":[]});
    for (mode, state_name) in [
        ("optimized", "gammaloop_state_GL638_optimized_lmbs"),
        (
            "augmented",
            "gammaloop_state_GL638_advanced_sampling_optimized_lmbs",
        ),
    ] {
        eprintln!("loading {mode}");
        let frozen = base.join(format!("frozen_{mode}"));
        let grid_export: Value =
            serde_json::from_str(&fs::read_to_string(base.join(format!("{mode}_grid.json")))?)?;
        let grid: Grid<F<f64>> = serde_json::from_value(grid_export["grid"].clone())?;
        let mut loaded = StateLoadOption::read_only(nnlo.join(state_name)).load()?;
        loaded.state.activate_loaded_integrand_backends(false)?;
        let pid = loaded
            .state
            .resolve_process_ref(Some(&ProcessRef::Name("epem_a_tth".into())))?;
        let model = &loaded.state.model;
        let runtime = loaded.state.process_list.get_integrand_mut(pid, "NNLO")?;
        *runtime.get_mut_settings() = RuntimeSettings::from_file(
            frozen.join("integrands/epem_a_tth@NNLO/settings.toml"),
            "retained frozen proposal",
        )?;
        runtime.warm_up(model)?;
        let ProcessIntegrand::CrossSection(process) = runtime else {
            return Err(eyre!("cross section required"));
        };
        let policy = process.sampling_source_policy()?;
        let source_factor = if mode == "optimized" {
            let setup = process.data.graph_terms[0].sampling_setup();
            let (raw, factor) = match policy.0 {
                SamplingPrecision::Quad => source(setup.sampling_bridge::<f128>()?, &cube)?,
                SamplingPrecision::Fixed256 => {
                    source(setup.sampling_bridge::<SamplingFloat>()?, &cube)?
                }
                SamplingPrecision::Arb => source(setup.sampling_bridge::<ArbPrec>()?, &cube)?,
                SamplingPrecision::Double => unreachable!(),
            };
            common_raw = raw;
            factor
        } else {
            0.0
        };
        process.prepare_sampling_precision::<ArbPrec>()?;
        let bridge = process.data.graph_terms[0]
            .sampling_setup()
            .sampling_bridge::<ArbPrec>()?
            .clone();
        let mut result = json!({"mode":mode,"completed_iteration":retained[mode]["completed_iteration"],"completed_points":retained[mode]["completed_points"],"source_policy":format!("{:?}",policy),"channel_count":bridge.channels().len(),"channels":[]});
        if mode == "augmented" {
            result["joint_channel_compiled_debug"] = json!(format!("{:?}", bridge.channels()[12]));
        }
        if mode == "optimized" {
            let sample = Sample::Uniform(F(1.0), vec![0, 3], cube.iter().map(|x| F(*x)).collect());
            let original = evaluation(runtime.evaluate_sample_precise(
                &sample,
                model,
                F(1.0),
                false,
                Complex::new_zero(),
            ));
            ensure!(
                original["error"].is_null() && original["metadata"]["is_nan"] == false,
                "source evaluation failed: {original}"
            );
            for k in 0..2 {
                physical[k] = original["complete_estimator"][k].as_f64().unwrap() / source_factor;
            }
            result["original_source_replay"] = original;
            result["original_source_J_times_partition"] = json!(source_factor);
            result["retained_negative_max_sample"] =
                retained[mode]["retained_negative_max_sample"].clone();
            let raw_strings = common_raw
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>();
            let trace: Value =
                serde_json::from_str(&fs::read_to_string(base.join("raw_native_trace.json"))?)?;
            let expected = trace["raw"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|x| x.as_str().unwrap().to_owned())
                })
                .collect::<Vec<_>>();
            ensure!(
                raw_strings == expected,
                "canonical raw point differs from retained native trace"
            );
            report["common_raw_coordinates_decimal"] = json!(raw_strings);
            report["physical_integrand_inferred_from_exact_original_draw"] = json!(physical);
        }
        let raw_input = MomentumSpaceEvaluationInput {
            loop_momenta: common_raw
                .chunks(3)
                .map(|x| ThreeMomentum {
                    px: F(x[0].into_f64()),
                    py: F(x[1].into_f64()),
                    pz: F(x[2].into_f64()),
                })
                .collect(),
            integrator_weight: F(1.0),
            graph_id: Some(0),
            group_id: None,
            orientation: None,
            channel_id: None,
        };
        result["physical_integrand_control_binary64_raw_input"] =
            evaluation(runtime.evaluate_momentum_configuration_precise(model, &raw_input, false));
        let mut raw_density_sum = 0.0;
        let mut adaptive_density_sum = 0.0;
        for i in 0..bridge.channels().len() {
            eprintln!("{mode} inverse channel {i}");
            let name = bridge.channels()[i].name.clone();
            let inverse = match bridge.inverse(SamplingChannelId(i), &common_raw) {
                Ok(Some(v)) => v,
                Ok(None) => {
                    result["channels"]
                        .as_array_mut()
                        .unwrap()
                        .push(json!({"id":i,"name":name,"support":false}));
                    continue;
                }
                Err(e) => {
                    result["channels"]
                        .as_array_mut()
                        .unwrap()
                        .push(json!({"id":i,"name":name,"inverse_error":format!("{e:#}")}));
                    continue;
                }
            };
            let coords = inverse
                .map
                .coordinates
                .iter()
                .map(FloatLike::into_f64)
                .collect::<Vec<_>>();
            let outer = grid
                .probe(&Probe::Discrete(
                    vec![0, i],
                    coords.iter().map(|x| Some(F(*x))).collect(),
                ))
                .map_err(|e| eyre!(e))?
                .0;
            let discrete = grid
                .probe(&Probe::Discrete(vec![0, i], vec![]))
                .map_err(|e| eyre!(e))?
                .0;
            let q = inverse.map.inverse_jacobian.into_f64();
            let factor = inverse.selected_factor()?.into_f64();
            raw_density_sum += q;
            adaptive_density_sum += q / outer;
            let sample =
                Sample::Uniform(F(outer), vec![0, i], coords.iter().map(|x| F(*x)).collect());
            let actual = evaluation(runtime.evaluate_sample_precise(
                &sample,
                model,
                F(outer),
                false,
                Complex::new_zero(),
            ));
            let forward = bridge.forward(
                SamplingChannelId(i),
                &coords
                    .iter()
                    .map(|x| ArbPrec::from_f64_exact_binary(*x))
                    .collect::<Vec<_>>(),
            );
            let roundtrip = forward.map(|v| {
                v.raw_coordinates
                    .iter()
                    .zip(&common_raw)
                    .map(|(a, b)| (a.into_f64() - b.into_f64()).abs())
                    .fold(0.0, f64::max)
            });
            result["channels"].as_array_mut().unwrap().push(json!({"id":i,"name":name,"support":true,"inverse_coordinates":coords,"inverse_coordinates_decimal":inverse.map.coordinates.iter().map(ToString::to_string).collect::<Vec<_>>(),"raw_map_density":q,"map_jacobian":inverse.map.jacobian.into_f64(),"raw_partition_alpha":inverse.partition.weight(i).unwrap().into_f64(),"J_times_alpha":factor,"discrete_probability":1.0/discrete,"continuous_cube_density":discrete/outer,"full_outer_inverse_adaptive_density":outer,"adaptive_raw_density_contribution":q/outer,"predicted_conditional_estimator_at_exact_common_raw":[physical[0]*factor*outer,physical[1]*factor*outer],"actual_production_evaluation_at_rounded_inverse_cube":actual,"rounded_cube_raw_roundtrip_max_abs_GeV":roundtrip.ok(),"map_residual":inverse.map.residual.to_string(),"map_diagnostics":inverse.map.diagnostics}));
            fs::write(
                base.join(format!("{mode}_partial.json")),
                serde_json::to_vec_pretty(&result)?,
            )?;
        }
        result["raw_map_density_sum"] = json!(raw_density_sum);
        result["adaptive_mixture_raw_density"] = json!(adaptive_density_sum);
        result["hypothetical_full_adaptive_mixture_estimator"] = json!([
            physical[0] / adaptive_density_sum,
            physical[1] / adaptive_density_sum
        ]);
        result["hypothetical_mixture_scope"] = json!(
            "Diagnostic only: production partition is raw-map density normalized before learned adaptive channel and cube factors, so production conditional estimates can differ by channel."
        );
        report["modes"].as_array_mut().unwrap().push(result);
        fs::write(
            base.join("inverse_replay.json"),
            serde_json::to_vec_pretty(&report)?,
        )?;
    }
    Ok(())
}
