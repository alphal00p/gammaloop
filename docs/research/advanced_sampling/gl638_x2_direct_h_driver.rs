// Standalone diagnostic: uses existing public production/acceptance owners.
use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use figment::{providers::Serialized, Figment};
use gammaloop_api::{commands::set::SetArgs, StateLoadOption};
use gammalooprs::{
    graph::GroupId,
    initialisation::initialise,
    integrands::process::{
        cross_section::CrossSectionGraphTerm,
        sampling_maps::{SamplingEvaluationError, SamplingMapAffine},
        GaussianReferenceFunction, GraphTerm, ProcessIntegrand, SamplingChannelId,
    },
    settings::RuntimeSettings,
    utils::{ArbPrec, FloatLike, QuadFloat, F},
    DependentMomentaConstructor,
};
use serde_json::{json, Value};
use std::{env, fs, time::Instant};

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
fn energy(v: &[f64], mass: f64) -> f64 {
    (norm(v).powi(2) + mass * mass).sqrt()
}
fn scaled(v: &[f64], scale: f64) -> Vec<f64> {
    v.iter().map(|x| x * scale).collect()
}
fn sum(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(a, b)| a + b).collect()
}

// Independent GL638 oracle from the recorded cut-1 equations, never used to
// parameterize or evaluate the production proposal.
fn geometry(raw: &[f64], radial_coordinate: Option<f64>, power: f64, beta: f64) -> Value {
    let p = &raw[..3];
    let a = raw[3..6]
        .iter()
        .zip(p)
        .map(|(k4, p)| k4 - p)
        .collect::<Vec<_>>();
    let t = &raw[9..12];
    let cut = |scale| {
        energy(&scaled(&a, scale), 125.0)
            + energy(&scaled(&sum(&a, t), scale), 173.0)
            + energy(&scaled(t, scale), 173.0)
            - 1000.0
    };
    let (mut lo, mut hi) = (0.0, 1.0);
    while cut(hi) < 0.0 {
        hi *= 2.0;
    }
    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        if cut(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let tau = (lo + hi) / 2.0;
    let (a, p, t) = (scaled(&a, tau), scaled(p, tau), scaled(t, tau));
    let b = energy(&sum(&a, &t), 173.0) + energy(&t, 173.0);
    let h = energy(&sum(&a, &p), 173.0) + energy(&p, 173.0) - b;
    let z = energy(&p, 173.0)
        + energy(&t, 173.0)
        + norm(&p.iter().zip(&t).map(|(p, t)| p - t).collect::<Vec<_>>())
        - 1000.0;
    let displacement = p
        .iter()
        .zip(&a)
        .map(|(p, a)| p + a / 2.0)
        .collect::<Vec<_>>();
    let radius = norm(&displacement);
    let n = scaled(&displacement, 1.0 / radius);
    let an = a.iter().zip(&n).map(|(a, n)| a * n).sum::<f64>();
    let root = ((b * b - norm(&a).powi(2) - 4.0 * 173.0f64.powi(2))
        / (4.0 * (1.0 - (an / b).powi(2))))
    .sqrt();
    let predicted = radial_coordinate.map(|u| {
        let split = root / (root + beta);
        if u < split {
            root * (1.0 - (1.0 - u / split).powf(power))
        } else {
            root + beta * ((u - split) / (1.0 - u)).powf(power)
        }
    });
    json!({"tau":tau,"cut_residual":cut(tau),"H":h,"Z":z,"prepared_center":scaled(&a,-0.5),"prepared_radius":radius,"analytic_root":root,"predicted_radius":predicted,"radius_relative_error":predicted.map(|r|(radius/r-1.0).abs())})
}

fn native_rays<T: FloatLike>(
    term: &CrossSectionGraphTerm,
    settings: &RuntimeSettings,
    externals: &[[f64; 4]],
    points: &[Value],
) -> Result<Vec<Value>> {
    let native_externals = externals
        .iter()
        .map(|p| p.map(|v| F::<T>::from_f64(v).0))
        .collect::<Vec<_>>();
    let parameterization = settings.sampling.get_parameterization_settings().unwrap();
    let bridge =
        term.compile_sampling_bridge(&parameterization, settings, &native_externals, None)?;
    points
        .iter()
        .map(|point| {
            let evaluated = (|| -> Result<Value> {
                let tokens = point["raw_point_tokens"].as_array().unwrap();
                // These tokens denote archived binary64 raw points. Preserve
                // those inputs while rebuilding all derived arithmetic natively.
                let raw = tokens
                    .iter()
                    .map(|v| v.as_str().unwrap().parse::<f64>())
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                let native = raw
                    .iter()
                    .map(|v| T::from_f64_exact_binary(*v))
                    .collect::<Vec<_>>();
                let inverse = bridge.inverse(SamplingChannelId(0), &native)?;
                let forward = bridge.forward(SamplingChannelId(0), &inverse.map.coordinates)?;
                let residual = forward
                    .raw_coordinates
                    .iter()
                    .zip(&native)
                    .map(|(a, b)| (F(a.clone()) - F(b.clone())).0.into_f64().abs())
                    .fold(0.0, f64::max);
                ensure!(residual < 1.0e-7, "{} roundtrip residual {residual}", point["name"]);
                Ok(json!({
                    "name": point["name"], "success": true,
                    "jacobian": inverse.map.jacobian.to_string(),
                    "density": inverse.map.inverse_jacobian.to_string(),
                    "log_density": inverse.partition.log_denominator.to_string(),
                    "roundtrip_residual_GeV": residual,
                    "inverse_coordinates": inverse.map.coordinates.iter().map(ToString::to_string).collect::<Vec<_>>(),
                    "geometry": geometry(&raw, None, parameterization.power, settings.kinematics.e_cm * parameterization.b)
                }))
            })();
            match evaluated {
                Ok(value) => Ok(value),
                Err(error) => Ok(json!({
                    "name": point["name"], "success": false,
                    "retryable": error.downcast_ref::<SamplingEvaluationError>().is_some(),
                    "error": format!("{error:?}")
                })),
            }
        })
        .collect()
}

fn main() -> Result<()> {
    initialise()?;
    let args = env::args().skip(1).collect::<Vec<_>>();
    ensure!(
        args.len() == 6,
        "state overlay rays points output_prefix expected_orientations"
    );
    let expected_orientations: usize = args[5].parse()?;
    let started = Instant::now();
    let mut loaded = StateLoadOption::read_only(&args[0]).load()?;
    ensure!(loaded.is_read_only_state(), "state must be read-only");
    loaded.state.activate_loaded_integrand_backends(false)?;
    let model = &loaded.state.model;
    let integrand = loaded.state.process_list.get_integrand_mut(0, "NNLO")?;
    let input = SetArgs::String {
        string: fs::read_to_string(&args[1])?,
    };
    let settings = input
        .merge_figment(Figment::from(Serialized::defaults(
            integrand.get_settings(),
        )))?
        .extract()?;
    *integrand.get_mut_settings() = settings;
    std::thread::scope(|scope| -> Result<()> {
        std::thread::Builder::new()
            .stack_size(128 * 1024 * 1024)
            .spawn_scoped(scope, || integrand.warm_up(model))?
            .join()
            .map_err(|_| eyre!("warmup thread panicked"))?
    })?;
    let warmup_seconds = started.elapsed().as_secs_f64();
    fs::write(
        format!("{}.settings.toml", args[4]),
        toml::to_string_pretty(integrand.get_settings())?,
    )?;
    let exposed_orientation_selector_count = integrand.group_orientation_count(GroupId(0)).unwrap();
    let settings = integrand.get_settings().clone();
    let parameterization = settings.sampling.get_parameterization_settings().unwrap();
    let ids = integrand.group_sampling_channel_ids(GroupId(0), &parameterization)?;
    ensure!(
        ids == vec![SamplingChannelId(0)],
        "gate requires exactly one canonical channel"
    );
    let external = settings
        .kinematics
        .externals
        .get_dependent_externals::<f64>(DependentMomentaConstructor::CrossSection)?;
    let externals = external
        .iter()
        .map(|p| {
            [
                p.temporal.value.0,
                p.spatial.px.0,
                p.spatial.py.0,
                p.spatial.pz.0,
            ]
        })
        .collect::<Vec<_>>();
    let ProcessIntegrand::CrossSection(prepared) = &*integrand else {
        return Err(eyre!("expected cross section"));
    };
    ensure!(
        prepared.data.graph_terms.len() == 1,
        "gate requires exactly one graph term"
    );
    let term = &prepared.data.graph_terms[0];
    ensure!(
        term.selected_production_orientation_keys().len() == expected_orientations,
        "expected {expected_orientations} production orientation keys, got {}; exposed selector count is {exposed_orientation_selector_count}",
        term.selected_production_orientation_keys().len()
    );
    ensure!(term.graph.name == "GL638", "expected GL638");
    let generation_parent = term
        .graph
        .loop_momentum_basis
        .loop_edges
        .iter()
        .map(|edge| edge.0)
        .collect::<Vec<_>>();
    ensure!(
        generation_parent == [3, 4, 7, 10],
        "archived rays require generation parent [3,4,7,10], got {generation_parent:?}"
    );
    let cuts = term
        .cut_esurface
        .iter_enumerated()
        .map(|(id, surface)| {
            json!({
                "id": id.0,
                "energy_edges": surface.energies.iter().map(|edge| edge.0).collect::<Vec<_>>(),
                "external_shift": format!("{:?}", surface.external_shift)
            })
        })
        .collect::<Vec<_>>();
    let bridge = term.compile_sampling_bridge(&parameterization, &settings, &externals, None)?;
    let cube = [
        0.23, 0.37, 0.41, 0.29, 0.53, 0.67, 0.43, 0.61, 0.71, 0.31, 0.59, 0.79,
    ];
    let point = bridge.forward(SamplingChannelId(0), &cube)?;
    let mut finite_differences = Vec::new();
    for step in [1.0e-5, 2.0e-6] {
        let mut jacobian = vec![vec![0.0; 12]; 12];
        for axis in 0..12 {
            let mut plus = cube;
            let mut minus = cube;
            plus[axis] += step;
            minus[axis] -= step;
            let plus = bridge.forward(SamplingChannelId(0), &plus)?;
            let minus = bridge.forward(SamplingChannelId(0), &minus)?;
            for (row, values) in jacobian.iter_mut().enumerate() {
                values[axis] =
                    (plus.raw_coordinates[row] - minus.raw_coordinates[row]) / (2.0 * step);
            }
        }
        let off_diagonal = jacobian[..3]
            .iter()
            .flat_map(|row| row[..9].iter())
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        ensure!(
            off_diagonal > 1.0e-3,
            "conditional off-diagonal derivatives must be nontrivial"
        );
        let determinant = SamplingMapAffine::new(jacobian, vec![0.0; 12])?.determinant();
        let relative_error = (determinant / point.map.jacobian - 1.0).abs();
        ensure!(
            relative_error < 3.0e-5,
            "full12D determinant mismatch {relative_error}"
        );
        finite_differences.push(json!({"step":step,"determinant":determinant,"analytic_map_J":point.map.jacobian,"relative_error":relative_error,"max_cut_to_active_derivative":off_diagonal}));
    }
    let analytic = geometry(
        &point.raw_coordinates,
        Some(cube[9]),
        parameterization.power,
        settings.kinematics.e_cm * parameterization.b,
    );
    ensure!(
        analytic["radius_relative_error"].as_f64().unwrap() < 1.0e-9,
        "graph chart disagrees with independent H ellipsoid: {analytic}"
    );
    let rays: Value = serde_json::from_str(&fs::read_to_string(&args[2])?)?;
    let points = rays["points"].as_array().unwrap();
    let f64_rays = native_rays::<f64>(term, &settings, &externals, points)?;
    let quad_rays = native_rays::<QuadFloat>(term, &settings, &externals, points)?;
    let arb_rays = native_rays::<ArbPrec>(term, &settings, &externals, points)?;
    ensure!(
        arb_rays.iter().all(|v| v["success"] == true),
        "final native inverse failures: {arb_rays:?}"
    );
    for (source, result) in points.iter().zip(&arb_rays) {
        ensure!(
            (source["H_GeV"].as_f64().unwrap() - result["geometry"]["H"].as_f64().unwrap()).abs()
                < 1.0e-8,
            "H oracle disagrees"
        );
        ensure!(
            (source["Z_GeV"].as_f64().unwrap() - result["geometry"]["Z"].as_f64().unwrap()).abs()
                < 1.0e-8,
            "Z oracle disagrees"
        );
    }
    let geometry_report = json!({"scope":"graph-backed direct H; no physical integration; H/Z values are independent f64 oracle, not direct native Esurface outputs","orientation_scope":"full visible orientation catalogue loaded; reference substitution does not evaluate physical orientation amplitudes","exposed_orientation_selector_count":exposed_orientation_selector_count,"expected_orientations":expected_orientations,"generation_parent":generation_parent,"cut_esurface_registry":cuts,"orientation_keys":term.selected_production_orientation_keys(),"warmup_seconds":warmup_seconds,"forward_point":point.raw_coordinates,"full_12D_finite_differences":finite_differences,"forward_J":point.map.jacobian,"analytic":analytic,"f64_rays":f64_rays,"quad_rays":quad_rays,"arb_rays":arb_rays});
    fs::write(
        format!("{}.geometry.json", args[4]),
        serde_json::to_string_pretty(&geometry_report)?,
    )?;
    let n: usize = args[3].parse()?;
    if n == 0 {
        println!("geometry-only determinant/native-inverse gate passed");
        return Ok(());
    }
    let coordinates = (1..=n)
        .map(|index| {
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
                .map(|base| {
                    let (mut index, mut fraction, mut value) = (index, 1.0, 0.0);
                    while index > 0 {
                        fraction /= base as f64;
                        value += (index % base) as f64 * fraction;
                        index /= base;
                    }
                    value
                })
                .to_vec()
        })
        .collect::<Vec<_>>();
    let reference = GaussianReferenceFunction::new(300.0, [30.0, -20.0, 10.0].repeat(4))?;
    let report =
        integrand.evaluate_reference_discrete_coordinates(&[0, 0], &coordinates, &reference)?;
    let passed = report.finite_sample_count == n
        && (report.normalization - 1.0).abs() < 0.06
        && (report.second_moment / report.expected_second_moment - 1.0).abs() < 0.08;
    let result = json!({"passed":passed,"scope":"all visible orientations summed, one graph, one channel; complete normalized reference estimator", "production_orientation_count": expected_orientations,"exposed_orientation_selector_count":exposed_orientation_selector_count, "dispersion_note":"Per-draw standard-error arithmetic on deterministic Halton values; not randomized-QMC uncertainty","points":n,"finite_points":report.finite_sample_count,"normalization":report.normalization,"normalization_dispersion":report.normalization_stderr,"raw_second_moment":report.second_moment,"expected_raw_second_moment":report.expected_second_moment,"second_moment_dispersion":report.second_moment_stderr,"elapsed_seconds":started.elapsed().as_secs_f64()});
    fs::write(
        format!("{}.reference.json", args[4]),
        serde_json::to_string_pretty(&result)?,
    )?;
    println!("{result}");
    ensure!(passed, "normalized reference acceptance failed");
    Ok(())
}
