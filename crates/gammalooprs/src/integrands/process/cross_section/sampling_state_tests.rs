//! Opt-in boundary gate for a generated state, using its actual production draw.
//!
//! Set `GL_SAMPLING_STATE_REQUEST` to a JSON object with `state`, `process`,
//! `integrand`, `settings` (a complete RuntimeSettings JSON file), `samples`
//! (a JSON array of complete Sample values, at most 32), and a new `report` path.
//! Also supply `graph`, `orientations`, `cuts`, `threshold_variants`,
//! `sampling_precision` (Debug spelling), `relative_tolerance` (number), and
//! `absolute_tolerance` (decimal string, in the returned physical units).
//! `minimum_hosts` and `minimum_normal_certificates` are integer arrays with
//! one entry per sample: declare which cases must exercise hosted/joint checks.
//! The caller authenticates the state/build/input hashes before and after this
//! test. This gate neither saves the state nor regenerates a failed draw.

use super::{EsurfaceRay, SamplingCatalogueEntry};
use crate::{
    GammaLoopContextContainer,
    graph::FeynmanGraph,
    initialisation::initialise,
    integrands::{
        evaluation::EvaluationMetaData,
        process::{
            EvaluationSource, EvaluationTarget, GraphTerm, ProcessIntegrand, ProcessIntegrandImpl,
            evaluate_all_rotations, evaluate_single,
        },
    },
    model::{InputParamCard, Model},
    momentum::{Rotation, RotationMethod, sample::Subspace},
    processes::ProcessList,
    utils::{ArbPrec, F, f128, newton_solver::RadialRootIdentity},
};
use color_eyre::eyre::{Result, WrapErr, ensure, eyre};
use serde_json::{Value, json};
use std::{fs, path::Path};
use symbolica::{domains::float::SingleFloat, numerical_integration::Sample};

// Complete estimators and cut totals use their own component magnitudes. Only
// their decomposition may use the sum of absolute contributions, separately
// for real and imaginary parts; a large imaginary part cannot mask a real error.
fn agrees(
    left: &F<ArbPrec>,
    right: &F<ArbPrec>,
    absolute: &F<ArbPrec>,
    relative: &F<ArbPrec>,
    component_scale: Option<&F<ArbPrec>>,
) -> bool {
    let scale = component_scale
        .cloned()
        .unwrap_or_else(|| left.abs().max(right.abs()));
    left.0.is_finite()
        && right.0.is_finite()
        && scale.0.is_finite()
        && scale >= scale.zero()
        && (left - right).abs() <= absolute + relative * scale
}

#[test]
fn saved_state_component_comparison_rejects_perturbations() -> Result<()> {
    let number = |value: &str| value.parse::<ArbPrec>().map(F);
    let zero = number("0")?;
    let relative = number("1e-10")?;
    let real_scale = number("2e-9")?;
    let imaginary_scale = number("1")?;
    let compare = |left, right, scale| agrees(left, right, &zero, &relative, scale);
    let noise = number("1e-308")?;
    assert!(!compare(&zero, &noise, None));
    assert!(compare(&zero, &noise, Some(&real_scale)));

    // These opposite term errors cancel in the cut total. They would be hidden
    // by an imaginary scale but must fail against the real contributions.
    let terms = [number("1e-9")?, number("-1e-9")?];
    let shifted = [&terms[0] + number("1e-18")?, &terms[1] - number("1e-18")?];
    for (term, shifted) in terms.iter().zip(&shifted) {
        assert!(!compare(term, shifted, Some(&real_scale)));
        assert!(compare(term, shifted, Some(&imaginary_scale)));
    }
    let original_sum = &terms[0] + &terms[1];
    let shifted_sum = &shifted[0] + &shifted[1];
    assert!(compare(&original_sum, &shifted_sum, None));

    // A tiny cut total still has to pass its own strict comparison, even when
    // the same difference is negligible at the decomposition's real scale.
    let total = number("1e-12")?;
    let changed_total = &total + number("1e-21")?;
    assert!(compare(&total, &changed_total, Some(&real_scale)));
    assert!(!compare(&total, &changed_total, None));
    Ok(())
}

#[test]
#[ignore = "requires an explicitly authenticated generated state and bounded Sample deck"]
fn native_sampling_saved_state_boundary() -> Result<()> {
    let request: Value =
        serde_json::from_reader(fs::File::open(std::env::var("GL_SAMPLING_STATE_REQUEST")?)?)?;
    let text = |key: &str| {
        request[key]
            .as_str()
            .ok_or_else(|| eyre!("missing string request field {key}"))
    };
    let count = |key: &str| {
        request[key]
            .as_u64()
            .map(|value| value as usize)
            .ok_or_else(|| eyre!("missing integer request field {key}"))
    };
    let state = Path::new(text("state")?);
    let report_path = Path::new(text("report")?);
    ensure!(!report_path.exists(), "refusing to overwrite gate evidence");
    ensure!(
        !report_path
            .parent()
            .unwrap()
            .canonicalize()?
            .starts_with(state.canonicalize()?),
        "the report must be outside the saved state"
    );
    let samples: Vec<Sample<F<f64>>> = serde_json::from_reader(fs::File::open(text("samples")?)?)?;
    ensure!(
        (1..=32).contains(&samples.len()),
        "request 1..=32 complete samples"
    );
    let minimum_hosts: Vec<usize> = serde_json::from_value(request["minimum_hosts"].clone())?;
    let minimum_normal_certificates: Vec<usize> =
        serde_json::from_value(request["minimum_normal_certificates"].clone())?;
    ensure!(
        minimum_hosts.len() == samples.len() && minimum_normal_certificates.len() == samples.len(),
        "minimum_hosts and minimum_normal_certificates need one entry per sample"
    );
    let relative = request["relative_tolerance"]
        .as_f64()
        .ok_or_else(|| eyre!("missing relative_tolerance"))?;
    ensure!(
        relative.is_finite() && relative > 0.0,
        "invalid relative tolerance"
    );
    let relative = F::<ArbPrec>::from_f64(relative);
    let absolute = F(text("absolute_tolerance")?.parse::<ArbPrec>()?);
    ensure!(
        absolute.0.is_finite() && absolute >= absolute.zero(),
        "invalid absolute tolerance"
    );
    let agrees = |left: &F<ArbPrec>, right: &F<ArbPrec>, scale: Option<&F<ArbPrec>>| {
        agrees(left, right, &absolute, &relative, scale)
    };

    initialise()?;
    // Preserve State::load's owner ordering, especially model parameter values
    // and custom UFO symbol printers before the Symbolica ID remapping.
    let mut model = Model::from_file(state.join("model.json"))?;
    model.apply_param_card(&InputParamCard::from_file(
        state.join("model_parameters.json"),
    )?)?;
    let state_map = symbolica::state::State::import(
        &mut fs::File::open(state.join("symbolica_state.bin"))?,
        None,
    )?;
    let mut processes = ProcessList::load(
        state,
        GammaLoopContextContainer {
            state_map: &state_map,
            model: &model,
        },
    )?;
    processes.activate_loaded_integrand_backends(false)?;
    let process_id = processes
        .processes
        .iter()
        .position(|process| process.definition.folder_name == text("process").unwrap())
        .ok_or_else(|| eyre!("requested process is absent"))?;
    let integrand = processes.get_integrand_mut(process_id, text("integrand")?)?;
    *integrand.get_mut_settings() = serde_json::from_reader(fs::File::open(text("settings")?)?)?;
    ensure!(
        integrand.get_settings().general.generate_events
            && integrand
                .get_settings()
                .general
                .store_additional_weights_in_event,
        "enable events and additional weights to compare original cuts and CTs"
    );
    integrand.warm_up(&model)?;
    let ProcessIntegrand::CrossSection(integrand) = integrand else {
        return Err(eyre!("this host gate requires a cross section"));
    };
    ensure!(
        integrand
            .get_rotations()
            .any(|rotation| !rotation.is_identity()),
        "configure a nonidentity physical stability probe"
    );
    ensure!(
        integrand.data.graph_terms.len() == 1 && integrand.data.graph_group_structure.len() == 1,
        "this bounded gate expects one graph and group"
    );
    let graph = integrand.get_graph(0);
    ensure!(graph.graph.name == text("graph")?, "graph identity changed");
    ensure!(
        graph.production_orientation_keys().len() == count("orientations")?
            && graph.selected_production_orientation_keys().len() == count("orientations")?,
        "orientation inventory changed"
    );
    ensure!(
        graph.cuts.len() == count("cuts")? && graph.cut_esurface.len() == count("cuts")?,
        "cut inventory changed"
    );
    ensure!(
        graph
            .threshold_counterterm_metadata()
            .ok_or_else(|| eyre!("missing threshold metadata"))?
            .variants
            .len()
            == count("threshold_variants")?,
        "threshold inventory changed"
    );
    let source_policy = *graph
        .sampling_setup()
        .sampling_source
        .as_ref()
        .ok_or_else(|| eyre!("missing source policy"))?;
    ensure!(
        format!("{:?}", source_policy.0) == text("sampling_precision")?,
        "unexpected source precision: {source_policy:?}"
    );
    let mut report =
        json!({"request": request, "source_policy": format!("{source_policy:?}"), "cases": []});
    let identity = Rotation::new(RotationMethod::Identity);
    let mut failed = false;
    for (index, sample) in samples.iter().enumerate() {
        let mut record = json!({"index": index, "sample": sample, "hosts": [], "lanes": []});
        let mut certified_hosts = 0;
        let mut normal_certificates = 0;
        let result = (|| -> Result<()> {
            let original = EvaluationSource::XSpace(sample);
            let mut metadata = EvaluationMetaData::new_empty();
            let mut anchor = original.prepare_draw(integrand, &mut metadata)?.unwrap();
            record["canonical_source"] = json!(format!("{anchor:?}"));
            let source_time = metadata.canonical_sampling_preparation_time;
            let policies = metadata.sampling_proposal_policies.clone();
            // Re-solve the ORIGINAL full graph equation at the exact promoted
            // point. A same-cube Arb forward would test a different source.
            for (_, rows) in &anchor.groups {
                for row in rows {
                    let graph = integrand.get_graph(row.graph_id);
                    let point = &row.sample;
                    let masses = graph.graph.get_real_mass_vector::<ArbPrec>(&model);
                    for host in &row.prepared_lu_hosts {
                        let cut = &graph.cut_esurface[host.plan.representative_cut_id];
                        let (ray, solved) = cut
                            .solve_lu_cut(
                                point.loop_moms(),
                                point.external_moms(),
                                &masses,
                                &graph.graph.loop_momentum_basis,
                                &F::from_f64(integrand.settings.kinematics.e_cm),
                                &mut EvaluationMetaData::new_empty().radial_root_diagnostics,
                                &RadialRootIdentity::new(format!(
                                    "state gate case {index}, {:?}",
                                    host.source
                                )),
                            )
                            .map_err(|error| eyre!("independent original host solve: {error:?}"))?;
                        record["hosts"]
                            .as_array_mut()
                            .unwrap()
                            .push(json!({"source": format!("{host:?}"),
                            "independent_root": format!("{solved:?}"), "normals": []}));
                        let host_record =
                            record["hosts"].as_array_mut().unwrap().last_mut().unwrap();
                        host.ray.verify_lu_candidate(
                            &ray,
                            &host.solution,
                            3 * host.plan.parent_lmb.len(),
                            &F::from_f64(source_policy.1),
                        )?;
                        certified_hosts += 1;
                        let independent_loops =
                            point.loop_moms().rescale(&solved.solution, Subspace::None);
                        let source_loops = point
                            .loop_moms()
                            .rescale(&host.solution.solution, Subspace::None);
                        let Some(SamplingCatalogueEntry::Named(channel)) =
                            graph.sampling_setup().sampling_catalogue.as_ref().and_then(
                                |catalogue| {
                                    catalogue
                                        .entries
                                        .get(host.source.generating_channel.index())
                                },
                            )
                        else {
                            return Err(eyre!("retained host lacks its named source channel"));
                        };
                        for block in &channel.blocks {
                            let sets = block.target.energy_edge_sets();
                            let [left, right] = sets.as_slice() else {
                                continue;
                            };
                            let Some(host_edges) = block.target.host_cut() else {
                                continue;
                            };
                            let group = &graph.cut_group_data.cut_groups[host.plan.cut_group_id];
                            let Some(cut_id) = group.cuts.iter().find(|id| {
                                let mut edges = graph.cut_esurface[**id]
                                    .energies
                                    .iter()
                                    .map(|edge| edge.0)
                                    .collect::<Vec<_>>();
                                edges.sort_unstable();
                                (channel.definition.on_cut.is_empty()
                                    || channel.definition.on_cut.contains(&id.0))
                                    && edges == host_edges
                            }) else {
                                continue;
                            };
                            let targets = [
                                graph.sampling_target_surface(
                                    *cut_id,
                                    block.target.cut_side(),
                                    left,
                                )?,
                                graph.sampling_target_surface(
                                    *cut_id,
                                    block.target.cut_side(),
                                    right,
                                )?,
                            ];
                            let original_normals = targets.map(|target| {
                                target.evaluate_routed_enclosed(
                                    &point.one(),
                                    &independent_loops,
                                    point.external_moms(),
                                    &masses,
                                    &graph.graph.loop_momentum_basis,
                                )
                            });
                            let [h, z] = original_normals;
                            let original_normals = [h?, z?];
                            let source_normals = targets.map(|target| {
                                target.evaluate_routed_enclosed(
                                    &point.one(),
                                    &source_loops,
                                    point.external_moms(),
                                    &masses,
                                    &graph.graph.loop_momentum_basis,
                                )
                            });
                            let [h, z] = source_normals;
                            let source_normals = [h?, z?];
                            let host_residual = cut.evaluate_routed_enclosed(
                                &point.one(),
                                &source_loops,
                                point.external_moms(),
                                &masses,
                                &graph.graph.loop_momentum_basis,
                            )?;
                            host_record["normals"].as_array_mut().unwrap().push(json!({"target": format!("{:?}", block.target),
                                "independent": format!("{original_normals:?}"), "source": format!("{source_normals:?}"),
                                "source_host_residual": format!("{host_residual:?}")}));
                            // The existing certificate uses the actual H/Z
                            // normal radius, including when it is very small.
                            EsurfaceRay::<ArbPrec>::verify_normal_alignment(
                                original_normals,
                                source_normals,
                                host_residual,
                                source_policy.1,
                            )?;
                            normal_certificates += 1;
                        }
                    }
                }
            }
            ensure!(
                certified_hosts >= minimum_hosts[index]
                    && normal_certificates >= minimum_normal_certificates[index],
                "case {index} certified {certified_hosts} hosts and {normal_certificates} normal pairs; expected at least {} hosts and {} normal pairs",
                minimum_hosts[index],
                minimum_normal_certificates[index]
            );
            anchor.prepare_physical_overlaps(integrand, &model, &mut metadata)?;
            let frozen = format!("{anchor:?}");
            let mut independent = anchor.clone();
            for (_, rows) in &mut independent.groups {
                for row in rows {
                    row.prepared_lu_hosts.clear();
                    row.physical_overlaps = None;
                }
            }
            independent.prepare_physical_overlaps(
                integrand,
                &model,
                &mut EvaluationMetaData::new_empty(),
            )?;
            let adopted = evaluate_single(
                integrand,
                EvaluationTarget::Physical(&model),
                &anchor,
                &identity,
                &mut EvaluationMetaData::new_empty(),
                Some(&anchor),
            )?;
            let solved = evaluate_single(
                integrand,
                EvaluationTarget::Physical(&model),
                &independent,
                &identity,
                &mut EvaluationMetaData::new_empty(),
                Some(&independent),
            )?;
            record["adopted_physics"] = json!(format!("{adopted:?}"));
            record["independent_physics"] = json!(format!("{solved:?}"));
            // Map J/partition factors are already inside these results. The
            // identical complete Sample outer weight cancels in this comparison.
            let mut pairs = vec![(
                &adopted.integrand_result,
                &solved.integrand_result,
                None,
                "complete estimator".to_owned(),
            )];
            ensure!(
                adopted.event_groups.len() == solved.event_groups.len(),
                "event group count changed"
            );
            for (left, right) in adopted.event_groups.iter().zip(solved.event_groups.iter()) {
                ensure!(left.len() == right.len(), "cut event count changed");
                for (left, right) in left.iter().zip(right.iter()) {
                    ensure!(
                        format!("{:?}", left.cut_info) == format!("{:?}", right.cut_info),
                        "cut event identity changed"
                    );
                    let cut = format!("{:?}", left.cut_info);
                    pairs.push((&left.weight, &right.weight, None, cut.clone()));
                    let left = left.additional_weights.threshold_counterterms.as_ref();
                    let right = right.additional_weights.threshold_counterterms.as_ref();
                    ensure!(left.is_some() == right.is_some(), "CT presence changed");
                    if let (Some(left), Some(right)) = (left, right) {
                        let scales = [false, true].map(|imaginary| {
                            let left_sum = std::iter::once(&left.original)
                                .chain(left.components.iter().map(|component| &component.weighted))
                                .map(|value| {
                                    if imaginary {
                                        value.im.abs()
                                    } else {
                                        value.re.abs()
                                    }
                                })
                                .fold(left.original.re.zero(), |sum, value| sum + value);
                            let right_sum = std::iter::once(&right.original)
                                .chain(right.components.iter().map(|component| &component.weighted))
                                .map(|value| {
                                    if imaginary {
                                        value.im.abs()
                                    } else {
                                        value.re.abs()
                                    }
                                })
                                .fold(right.original.re.zero(), |sum, value| sum + value);
                            left_sum.max(right_sum)
                        });
                        pairs.push((
                            &left.original,
                            &right.original,
                            Some(scales.clone()),
                            format!("{cut} Original"),
                        ));
                        ensure!(
                            left.components.len() == right.components.len(),
                            "CT count changed"
                        );
                        for (left, right) in left.components.iter().zip(&right.components) {
                            ensure!(
                                left.component_id == right.component_id
                                    && left.occurrence == right.occurrence,
                                "CT occurrence changed"
                            );
                            pairs.push((
                                &left.weighted,
                                &right.weighted,
                                Some(scales.clone()),
                                format!(
                                    "{cut} {:?} occurrence {:?}",
                                    left.component_id, left.occurrence
                                ),
                            ));
                        }
                    }
                }
            }
            let mut all_agree = true;
            record["physics_comparisons"] = json!([]);
            for (left, right, scales, label) in pairs {
                for (index, (component, left, right)) in [
                    ("real", &left.re, &right.re),
                    ("imaginary", &left.im, &right.im),
                ]
                .into_iter()
                .enumerate()
                {
                    let selected_scale = scales.as_ref().map(|scales| &scales[index]);
                    let scale = selected_scale
                        .cloned()
                        .unwrap_or_else(|| left.abs().max(right.abs()));
                    let passed = agrees(left, right, selected_scale);
                    all_agree &= passed;
                    record["physics_comparisons"].as_array_mut().unwrap().push(json!({
                        "label": label, "component": component,
                        "criterion": if scales.is_some() { "sum_absolute_same_component" } else { "own_component_relative" },
                        "left": left.to_string(), "right": right.to_string(),
                        "scale": scale.to_string(), "absolute_error": (left - right).abs().to_string(),
                        "absolute_budget": (&absolute + &relative * scale).to_string(), "passed": passed,
                    }));
                }
            }
            ensure!(
                all_agree,
                "adopted/original Arb physics mismatch; see physics_comparisons"
            );
            let prepared = EvaluationSource::Prepared {
                sample: &anchor,
                original: &original,
            };
            ensure!(
                prepared.prepare_draw(integrand, &mut metadata)?.is_none(),
                "prepared draw was remapped"
            );
            macro_rules! lane {
                ($numeric:ty, $required:expr) => {{
                    let result = (|| -> Result<_> {
                        let native = prepared.build_gamma_sample::<$numeric, _>(integrand, &mut metadata)?;
                        let rotated = native.rotate(&Rotation::new(RotationMethod::Pi2Z), 0, 0);
                        for ((_, original_rows), (_, rotated_rows)) in anchor.groups.iter().zip(&rotated.groups) {
                            for (original, rotated) in original_rows.iter().zip(rotated_rows) {
                                ensure!(std::sync::Arc::ptr_eq(original.physical_overlaps.as_ref().unwrap(), rotated.physical_overlaps.as_ref().unwrap()),
                                    "rotation replaced canonical overlap authority");
                            }
                        }
                        evaluate_all_rotations(integrand, EvaluationTarget::Physical(&model), &native,
                            &mut EvaluationMetaData::new_empty(), true, Some(&anchor))
                    })();
                    record["lanes"].as_array_mut().unwrap().push(json!({"precision": stringify!($numeric),
                        "result": format!("{result:?}"), "required": $required}));
                    if $required {
                        let (results, _, _) = result.wrap_err("required final Arb lane failed")?;
                        for result in results {
                            ensure!(agrees(&result.integrand_result.re.to_arb_exact()?, &adopted.integrand_result.re, None)
                                && agrees(&result.integrand_result.im.to_arb_exact()?, &adopted.integrand_result.im, None),
                                "final Arb materialization/probe changed the complete estimator");
                        }
                    }
                    ensure!(format!("{anchor:?}") == frozen, "physical lane mutated canonical source");
                    ensure!(metadata.sampling_proposal_policies == policies, "physical lane changed source decisions");
                    ensure!(metadata.canonical_sampling_preparation_time == source_time, "physical lane remapped the source");
                }};
            }
            lane!(f64, false);
            lane!(f128, false);
            lane!(ArbPrec, true);
            Ok(())
        })();
        record["certified_hosts"] = json!(certified_hosts);
        record["normal_certificates"] = json!(normal_certificates);
        record["status"] = json!(if result.is_ok() { "passed" } else { "failed" });
        if let Err(error) = result {
            failed = true;
            record["error"] = json!(format!("{error:?}"));
        }
        report["cases"].as_array_mut().unwrap().push(record);
        fs::write(report_path, serde_json::to_vec_pretty(&report)?)?;
    }
    ensure!(
        !failed,
        "saved-state boundary failures retained in {}",
        report_path.display()
    );
    Ok(())
}
