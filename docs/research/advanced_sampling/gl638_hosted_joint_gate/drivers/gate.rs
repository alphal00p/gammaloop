// Off-source GL638 gate client, migrated from the archived X2 replay owners.
// No sampler, retry loop or integration engine lives here: Grid::sample and
// ProcessIntegrand's precise source/physical APIs own every numerical operation.
// Compile only against root's frozen optimized libraries. Never save the state.
use color_eyre::{
    Result,
    eyre::{ensure, eyre},
};
use figment::{Figment, providers::Serialized};
use gammaloop_api::{
    CLISettings, StateLoadOption,
    commands::{
        integrate::{Integrate, RendererOption, ShowPhaseOption},
        set::SetArgs,
    },
    state::ProcessRef,
};
use gammalooprs::{
    graph::GroupId,
    initialisation::initialise,
    integrands::{
        HasIntegrand,
        evaluation::{GenericEvaluationResult, PreciseEvaluationResult, StabilityStatus},
        process::{
            GaussianReferenceFunction, GraphTerm, MomentumSpaceEvaluationInput, ProcessIntegrand,
            SamplingChannelId,
        },
    },
    model::Model,
    momentum::ThreeMomentum,
    utils::{ArbPrec, F, FloatLike},
};
use serde_json::{Value, json};
use spenso::algebra::complex::Complex;
use std::{collections::BTreeMap, env, fs, path::Path, process::Command, time::Instant};
use symbolica::numerical_integration::{MonteCarloRng, Sample};

// File checks run only outside numerical timings. The saved generated payload,
// metadata, settings and model must match the reviewed complete 35-file snapshot.
fn state_hashes(root: &Path) -> Result<BTreeMap<String, String>> {
    let mut directories = vec![root.to_path_buf()];
    let mut paths = Vec::new();
    while let Some(directory) = directories.pop() {
        for entry in fs::read_dir(directory)? {
            let path = entry?.path();
            if path.is_dir() {
                directories.push(path);
            } else if path.is_file() {
                paths.push(path);
            }
        }
    }
    paths.sort();
    let mut hashes = BTreeMap::new();
    for path in paths {
        let output = Command::new("sha256sum").arg(&path).output()?;
        ensure!(
            output.status.success(),
            "sha256sum failed for {}",
            path.display()
        );
        let text = String::from_utf8(output.stdout)?;
        let hash = text
            .split_whitespace()
            .next()
            .ok_or_else(|| eyre!("empty sha256sum output"))?;
        hashes.insert(
            path.strip_prefix(root)?.to_string_lossy().into_owned(),
            hash.to_owned(),
        );
    }
    Ok(hashes)
}

fn inventory(integrand: &ProcessIntegrand, manifest: &Value, mode: &str) -> Result<Value> {
    let ProcessIntegrand::CrossSection(process) = integrand else {
        return Err(eyre!("cross section required"));
    };
    ensure!(
        process.data.graph_terms.len() == 1 && process.data.graph_group_structure.len() == 1,
        "exactly one graph and graph group required"
    );
    let term = &process.data.graph_terms[0];
    ensure!(
        term.graph.name == manifest["graph_name"].as_str().unwrap(),
        "wrong graph name"
    );
    let keys = term.selected_production_orientation_keys();
    ensure!(
        keys.len() == 936 && term.production_orientation_keys().len() == 936,
        "full936 generated/selected orientations required, got {}/{}",
        term.production_orientation_keys().len(),
        keys.len()
    );
    ensure!(
        term.cuts.len() == 6 && term.cut_esurface.len() == 6,
        "all six physical cuts required"
    );
    let generation_parent = term
        .graph
        .loop_momentum_basis
        .loop_edges
        .iter()
        .map(|e| e.0)
        .collect::<Vec<_>>();
    ensure!(
        json!(generation_parent) == manifest["generation_parent_lmb"],
        "generation parent mismatch"
    );
    let cuts = term
        .cut_esurface
        .iter()
        .enumerate()
        .map(|(id, surface)| {
            let mut edges = surface.energies.iter().map(|e| e.0).collect::<Vec<_>>();
            edges.sort();
            json!({"id":id,"edges":edges,"original_equation":format!("{surface:?}")})
        })
        .collect::<Vec<_>>();
    ensure!(
        cuts[manifest["host_cut"]["id"].as_u64().unwrap() as usize]["edges"]
            == manifest["host_cut"]["edges"],
        "host edge identity mismatch"
    );
    let parameterization = integrand
        .get_settings()
        .sampling
        .get_parameterization_settings()
        .ok_or_else(|| eyre!("parameterization missing"))?;
    let ids = integrand.group_sampling_channel_ids(GroupId(0), &parameterization)?;
    ensure!(
        ids.len() == manifest["expected_channels"][mode].as_u64().unwrap() as usize,
        "unexpected canonical channel count in {mode}: {ids:?}"
    );
    let channels = ids
        .iter()
        .map(|id| {
            Ok(json!({"id":id.0,"label":term.sampling_channel_label(*id,&parameterization)?,
                "is_lmb":term.sampling_setup().sampling_channel_is_lmb(*id,&term.graph.name,&parameterization)?}))
        })
        .collect::<Result<Vec<_>>>()?;
    ensure!(
        integrand.group_orientation_count(GroupId(0)) == Some(1),
        "expected one exposed summed-orientation selector"
    );
    ensure!(
        term.threshold_counterterm_metadata().is_some(),
        "threshold metadata registry missing"
    );
    Ok(
        json!({"graph":term.graph.name,"generation_parent_lmb":generation_parent,
        "production_orientation_keys":keys,"exposed_orientation_selectors":1,"cuts":cuts,"channels":channels,
        "threshold_counterterm_metadata":term.threshold_counterterm_metadata(),
        "uv_preservation_evidence":"complete generated-state hashes identical; no regeneration, filter or physical-setting override"}),
    )
}

// Preserve native decimal totals/events. Formatting and JSON allocation are
// deliberately after the timed precise call; no map factors are applied again.
fn native_record<T: FloatLike>(r: GenericEvaluationResult<T>, precision: &str) -> Value {
    let final_unstable = r
        .evaluation_metadata
        .stability_results
        .last()
        .is_some_and(|r| matches!(r.status, StabilityStatus::Unstable(_)));
    let finite = !final_unstable
        && !r.evaluation_metadata.is_nan
        && r.integrand_result.re.0.is_finite()
        && r.integrand_result.im.0.is_finite();
    let mut events = Vec::new();
    let mut finite_events = true;
    for (group_id, group) in r.event_groups.iter().enumerate() {
        for event in group.iter() {
            finite_events &= event.weight.re.0.is_finite() && event.weight.im.0.is_finite();
            events.push(json!({"group":group_id,"cut_info":event.cut_info,
                "weight":{"re":event.weight.re.to_string(),"im":event.weight.im.to_string()},
                "additional_weights":event.additional_weights.weights.iter().map(|(key,value)|json!({"key":format!("{key:?}"),"re":value.re.to_string(),"im":value.im.to_string()})).collect::<Vec<_>>()}));
        }
    }
    json!({"valid":finite && finite_events,"precision":precision,
        "integrand_result":{"re":r.integrand_result.re.to_string(),"im":r.integrand_result.im.to_string()},
        "parameterization_jacobian":r.parameterization_jacobian.map(|v|v.to_string()),
        "integrator_weight":r.integrator_weight.to_string(),"events":events,
        "sampling_seconds":r.evaluation_metadata.parameterization_time.as_secs_f64(),
        "physical_seconds":r.evaluation_metadata.integrand_evaluation_time.as_secs_f64(),
        "evaluation_metadata":r.evaluation_metadata})
}
fn record(result: Result<PreciseEvaluationResult>, wall: f64) -> Value {
    let evaluation = match result {
        Ok(PreciseEvaluationResult::Double(r)) => native_record(r, "Double"),
        Ok(PreciseEvaluationResult::Quad(r)) => native_record(r, "Quad"),
        Ok(PreciseEvaluationResult::Arb(r)) => native_record(r, "Arb"),
        Err(error) => {
            json!({"valid":false,"error":format!("{error:#}"),"sampling_seconds":null,"physical_seconds":null})
        }
    };
    let p = evaluation["physical_seconds"].as_f64();
    json!({"wall_seconds":wall,"sampling_upper_seconds":p.map(|p|(wall-p).max(0.0)),"evaluation":evaluation})
}

fn replay(
    integrand: &mut ProcessIntegrand,
    model: &Model,
    requests: &[Value],
    info: &Value,
    partition_diagnostic: bool,
) -> Result<Vec<Value>> {
    let mut rows = Vec::new();
    for case in requests {
        ensure!(
            case["threshold_ct_enabled"] == true && case["momentum_space"] == true,
            "only original CT-on momentum-space cases supported"
        );
        let point = case["point_tokens"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().parse::<f64>())
            .collect::<std::result::Result<Vec<_>, _>>()?;
        ensure!(
            point.len() == 12 && point.iter().all(|v| v.is_finite()),
            "four finite original loop momenta required"
        );
        let channel_id = if let Some(name) = case["channel_name"].as_str() {
            let matches = info["channels"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|entry| entry["label"] == name)
                .collect::<Vec<_>>();
            ensure!(
                matches.len() == 1,
                "expected unique named canonical channel {name}"
            );
            Some(SamplingChannelId(
                matches[0]["id"].as_u64().unwrap() as usize
            ))
        } else {
            None
        };
        let input = MomentumSpaceEvaluationInput {
            loop_momenta: point
                .chunks_exact(3)
                .map(|v| ThreeMomentum::from([F(v[0]), F(v[1]), F(v[2])]))
                .collect(),
            integrator_weight: F(1.0),
            graph_id: channel_id.is_none().then_some(0),
            group_id: channel_id.map(|_| GroupId(0)),
            orientation: None,
            channel_id,
        };
        let at = Instant::now();
        let result = integrand.evaluate_momentum_configuration_precise(
            model,
            &input,
            case["forced_arb"].as_bool().unwrap(),
        );
        let wall = at.elapsed().as_secs_f64();
        let mut row = record(result, wall);
        row["case"] = case.clone();
        if partition_diagnostic {
            if let Some(channel_id) = channel_id {
                // Independent raw partition oracle at the original binary64
                // point, in canonical Arb1000. This diagnostic is excluded from
                // W/S/P and never runs in hard timing or worker warmup.
                let canonical = point
                    .iter()
                    .map(|v| ArbPrec::from_f64_exact_binary(*v))
                    .collect::<Vec<_>>();
                let ProcessIntegrand::CrossSection(process) = &*integrand else {
                    return Err(eyre!("cross section required"));
                };
                let inverse = process.data.graph_terms[0]
                    .sampling_setup()
                    .sampling_bridge::<ArbPrec>()
                    .and_then(|bridge| bridge.inverse(channel_id, &canonical));
                row["canonical_raw_partition"] = match inverse {
                    Ok(Some(mapped)) => match mapped.partition.weight(channel_id.0) {
                        Some(weight) => {
                            json!({"status":"inside_support","weight":F(weight).to_string(),"precision":"Arb1000", "source":"exact binary64 promotion; independent canonical inverse; raw input uses w only, never inverse Jacobian"})
                        }
                        None => {
                            json!({"status":"error","error":"selected partition weight missing"})
                        }
                    },
                    Ok(None) => {
                        json!({"status":"outside_support","weight":null,"scope":"the selected chart does not cover this raw point; not a physical-law failure and no replacement chart is selected"})
                    }
                    Err(error) => json!({"status":"error","error":format!("{error:#}")}),
                };
                row["canonical_raw_partition"]["channel_id"] = json!(channel_id.0);
            }
        }
        row["point_definition"] =
            json!("original binary64 tokens, exactly promoted by existing source owner");
        row["factor_scope"] = json!(if channel_id.is_some() {
            "partition-weighted F; do not compare directly with bare F"
        } else {
            "bare F; no sampling partition/Jacobian"
        });
        let mut cuts = row["evaluation"]["events"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|event| event["cut_info"]["cut_id"].as_u64())
            .collect::<Vec<_>>();
        cuts.sort();
        row["six_distinct_cut_events"] = json!(cuts == vec![0, 1, 2, 3, 4, 5]);
        rows.push(row);
    }
    Ok(rows)
}

// Deterministic acceptance reuses the same production map and per-draw sum
// owner as the archived X2 client. Eight finite points are smoke, not an
// acceptance claim; only a finite 8192/32768 report can pass the fixed bounds.
fn reference_stage(
    integrand: &mut ProcessIntegrand,
    model: &Model,
    manifest: &Value,
    full: bool,
) -> Result<Value> {
    let settings = integrand.get_settings().clone();
    *integrand.get_mut_settings() = SetArgs::String {
        string: "sampling.sampling_channels = \"summed\"".into(),
    }
    .merge_figment(Figment::from(Serialized::defaults(&settings)))?
    .extract()?;
    let reference = GaussianReferenceFunction::new(
        manifest["reference"]["width_GeV"].as_f64().unwrap(),
        manifest["reference"]["center_GeV"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect(),
    )?;
    let counts = if full { vec![8, 8192, 32768] } else { vec![8] };
    let coordinates = (1..=*counts.last().unwrap())
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
    let result = std::thread::scope(|scope| -> Result<Value> {
        std::thread::Builder::new().stack_size(128*1024*1024).spawn_scoped(scope,||->Result<Value>{
            let at=Instant::now();integrand.warm_up(model)?;let warmup_seconds=at.elapsed().as_secs_f64();
            let mut rows=Vec::new();let mut passed=false;
            for count in counts {
                let at=Instant::now();
                // Graph0 is the sole discrete choice. The existing summed
                // owner combines every channel before its squared statistic.
                let r=match integrand.evaluate_reference_discrete_coordinates(&[0],&coordinates[..count],&reference) {
                    Ok(r)=>r,
                    Err(error)=>{
                        rows.push(json!({"points":count,"finite":false,"error":format!("{error:#}"),"elapsed_seconds":at.elapsed().as_secs_f64(),"coordinate_sequence":"first N twelve-dimensional Halton points from fixed primes recorded in source"}));
                        return Ok(json!({"full_acceptance_requested":full,"passed":false,"finite":false,"warmup_seconds":warmup_seconds,"rows":rows}));
                    }
                };
                if r.finite_sample_count!=count {
                    rows.push(json!({"points":count,"finite_points":r.finite_sample_count,"finite":false,"elapsed_seconds":at.elapsed().as_secs_f64()}));
                    return Ok(json!({"full_acceptance_requested":full,"passed":false,"finite":false,"warmup_seconds":warmup_seconds,"rows":rows}));
                }
                let within=(r.normalization-1.0).abs()<0.06&&(r.second_moment/r.expected_second_moment-1.0).abs()<0.08;
                passed=count>=8192&&within;
                rows.push(json!({"points":count,"finite_points":r.finite_sample_count,"normalization":r.normalization,
                    "normalization_dispersion":r.normalization_stderr,"raw_second_moment":r.second_moment,
                    "expected_raw_second_moment":r.expected_second_moment,"second_moment_dispersion":r.second_moment_stderr,
                    "elapsed_seconds":at.elapsed().as_secs_f64(),"acceptance_passed":passed,"smoke_only":count==8}));
                if passed {break;}
            }
            Ok(json!({"full_acceptance_requested":full,"passed":passed,"finite":true,"warmup_seconds":warmup_seconds,"rows":rows,
                "dispersion_scope":"IID-formula diagnostic on deterministic Halton points, not a QMC confidence interval"}))
        })?.join().map_err(|_|eyre!("reference worker panicked"))?
    });
    *integrand.get_mut_settings() = settings;
    result
}

// Workers clone the existing physical owner; its caches, grid and precise call
// remain authoritative. Setup/lazy warmup is reported separately, not in W/S/P.
fn time_sources(
    integrand: &ProcessIntegrand,
    model: &Model,
    workers: usize,
    draws: usize,
    repetitions: usize,
    seed: u64,
    hard_requests: &[Value],
    info: &Value,
) -> Result<Vec<Value>> {
    // Finish every worker's setup before any measured calls. Two bounded
    // scoped phases avoid a barrier deadlock on warmup/spawn/panic errors.
    let prepared = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for worker_id in 0..workers {
            let mut worker = integrand.clone();
            handles.push(std::thread::Builder::new().stack_size(128*1024*1024).spawn_scoped(scope,move ||->Result<_>{
                let at=Instant::now();let warmed=worker.warm_up(model);
                let mut warm_rows=Vec::new();let mut samples=Vec::new();let mut warm_hard=Vec::new();
                if warmed.is_ok() {
                    let mut grid=worker.create_grid();let mut rng=MonteCarloRng::new(seed,worker_id);
                    for index in 0..draws+2 {
                        let mut sample=Sample::new();grid.sample(&mut rng,&mut sample);
                        if index<2 {
                            let at=Instant::now();
                            let result=worker.evaluate_sample_precise(&sample,model,sample.get_weight(),false,Complex::new_zero());
                            let wall=at.elapsed().as_secs_f64();warm_rows.push(json!({"sample":sample,"result":record(result,wall)}));
                        } else {samples.push(sample);}
                    }
                    // The same cloned integrand also warms each raw/inverse
                    // hard path. Events remain disabled as in production cards.
                    warm_hard=replay(&mut worker,model,hard_requests,info,false)?;
                }
                if warm_rows.iter().any(|r:&Value|r["result"]["evaluation"]["valid"]!=true)
                    || warm_hard.iter().any(|r|r["evaluation"]["valid"]!=true) {samples.clear();}
                let setup=json!({"worker":worker_id,"setup_seconds":at.elapsed().as_secs_f64(),
                    "warmup_error":warmed.err().map(|e|format!("{e:#}")),"warmup_draws":warm_rows,"warmup_hard":warm_hard});
                Ok((worker_id,worker,samples,setup))
            })?);
        }
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .map_err(|_| eyre!("timing warmup worker panicked"))?
            })
            .collect::<Result<Vec<_>>>()
    })?;
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (worker_id, mut worker, samples, mut setup) in prepared {
            handles.push(std::thread::Builder::new().stack_size(128*1024*1024).spawn_scoped(scope,move ||->Result<Value>{
                let mut rows=Vec::new();let mut maximum:Option<(f64,usize)>=None;
                for repetition in 0..repetitions {
                    for (draw,sample) in samples.iter().enumerate() {
                        let at=Instant::now();
                        let result=worker.evaluate_sample_precise(sample,model,sample.get_weight(),false,Complex::new_zero());
                        let wall=at.elapsed().as_secs_f64();let r=record(result,wall);
                        if repetition==0 && r["evaluation"]["valid"]==true {
                            let value=&r["evaluation"]["integrand_result"];
                            let score=(value["re"].as_str().unwrap().parse::<f64>()?.abs()+value["im"].as_str().unwrap().parse::<f64>()?.abs())*sample.get_weight().0.abs();
                            if maximum.is_none_or(|(old,_)|score>old) {maximum=Some((score,draw));}
                        }
                        rows.push(json!({"worker":worker_id,"repetition":repetition,"draw":draw,"sample":sample,"result":r}));
                    }
                }
                let mut hard_rows=Vec::new();let mut maximum_rows=Vec::new();
                if !samples.is_empty() {
                    for repetition in 0..repetitions {
                        for r in replay(&mut worker,model,hard_requests,info,false)? {
                            hard_rows.push(json!({"worker":worker_id,"repetition":repetition,"case":r["case"],"result":r}));
                        }
                        if let Some((score,draw))=maximum {
                            let sample=&samples[draw];let at=Instant::now();
                            let result=worker.evaluate_sample_precise(sample,model,sample.get_weight(),false,Complex::new_zero());
                            let wall=at.elapsed().as_secs_f64();
                            maximum_rows.push(json!({"worker":worker_id,"repetition":repetition,"draw":draw,"sample":sample,
                                "reported_f64_absolute_score":score,"result":record(result,wall)}));
                        }
                    }
                }
                setup["rows"]=json!(rows);setup["hard_rows"]=json!(hard_rows);setup["maximum_rows"]=json!(maximum_rows);
                Ok(setup)
            })?);
        }
        handles
            .into_iter()
            .map(|h| h.join().map_err(|_| eyre!("timing worker panicked"))?)
            .collect::<Result<Vec<_>>>()
    })
}

fn timing_summary(workers: &[Value], field: &str, selected_hard: Option<bool>) -> Value {
    let rows = workers
        .iter()
        .flat_map(|w| w[field].as_array().into_iter().flatten())
        .filter(|row| {
            selected_hard.is_none_or(|selected| row["case"]["channel_name"].is_string() == selected)
        })
        .collect::<Vec<_>>();
    let mut wall = 0.0;
    let mut s = 0.0;
    let mut p = 0.0;
    let mut upper = 0.0;
    let mut ratios = Vec::new();
    let mut invalid = 0;
    for row in &rows {
        let r = &row["result"];
        wall += r["wall_seconds"].as_f64().unwrap();
        if r["evaluation"]["valid"] != true {
            invalid += 1;
        }
        if let (Some(si), Some(pi), Some(ui)) = (
            r["evaluation"]["sampling_seconds"].as_f64(),
            r["evaluation"]["physical_seconds"].as_f64(),
            r["sampling_upper_seconds"].as_f64(),
        ) {
            s += si;
            p += pi;
            upper += ui;
            if pi > 0.0 {
                ratios.push(ui / pi);
            }
        }
    }
    ratios.sort_by(f64::total_cmp);
    let percentile = |q: f64| {
        if ratios.is_empty() {
            None
        } else {
            Some(ratios[((ratios.len() - 1) as f64 * q).ceil() as usize])
        }
    };
    let mut distributions = serde_json::Map::new();
    for metric in [
        "wall_seconds",
        "sampling_seconds",
        "physical_seconds",
        "sampling_upper_seconds",
    ] {
        let mut values = rows
            .iter()
            .filter_map(|row| {
                let r = &row["result"];
                r[metric]
                    .as_f64()
                    .or_else(|| r["evaluation"][metric].as_f64())
            })
            .collect::<Vec<_>>();
        values.sort_by(f64::total_cmp);
        let at = |q: f64| {
            if values.is_empty() {
                None
            } else {
                Some(values[((values.len() - 1) as f64 * q).ceil() as usize])
            }
        };
        distributions.insert(metric.into(), json!({"count":values.len(),"p50":at(0.5),"p90":at(0.9),"p99":at(0.99),"max":values.last()}));
    }
    let warmup_failed = workers.iter().any(|w| {
        !w["warmup_error"].is_null()
            || w["warmup_draws"]
                .as_array()
                .unwrap()
                .iter()
                .any(|r| r["result"]["evaluation"]["valid"] != true)
            || w["warmup_hard"]
                .as_array()
                .unwrap()
                .iter()
                .any(|r| r["evaluation"]["valid"] != true)
    });
    json!({"measured_calls":rows.len(),"invalid_calls":invalid,"warmup_failed":warmup_failed,
        "distributions_seconds":distributions,"sum_wall_seconds":wall,"sum_sampling_seconds":s,"sum_physical_seconds":p,"sum_sampling_upper_seconds":upper,
        "sampling_over_physical":if p>0.0{Some(s/p)}else{None},"conservative_over_physical":if p>0.0{Some(upper/p)}else{None},
        "conservative_ratio_p50":percentile(0.5),"conservative_ratio_p90":percentile(0.9),"conservative_ratio_p99":percentile(0.99),"conservative_ratio_max":ratios.last(),
        "representative_cap_passed":!warmup_failed&&invalid==0&&p>0.0&&upper/p<=0.10,
        "scope":"aggregate elapsed precise-call times across workers, not CPU or batch-wall time; excludes Grid/RNG generation; max_eval=0 matches fresh iteration1, not trained iterations"})
}

fn main() -> Result<()> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    ensure!(
        (4..=8).contains(&args.len()),
        "usage: gate MANIFEST OUTPUT_DIR MODES(comma-separated or all) STAGE(inventory|smoke|reference|replay|timing|preflight|pilot|diagnostic-pilot) [WORKERS=1] [DRAWS=16] [REPETITIONS=3] [SEED=1337]"
    );
    let manifest: Value = serde_json::from_str(&fs::read_to_string(&args[0])?)?;
    let output = Path::new(&args[1]);
    ensure!(
        !output.join("summary.json").exists(),
        "refusing to overwrite completed evidence"
    );
    fs::create_dir_all(output)?;
    let modes = if args[2] == "all" {
        manifest["proposal_order"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s.as_str().unwrap().to_owned())
            .collect::<Vec<_>>()
    } else {
        args[2].split(',').map(str::to_owned).collect()
    };
    ensure!(
        [
            "inventory",
            "smoke",
            "reference",
            "replay",
            "timing",
            "preflight",
            "all",
            "pilot",
            "diagnostic-pilot"
        ]
        .contains(&args[3].as_str()),
        "unknown stage"
    );
    let stage = &args[3];
    let with_reference = matches!(
        stage.as_str(),
        "smoke" | "reference" | "preflight" | "all" | "pilot" | "diagnostic-pilot"
    );
    let full_reference = with_reference && stage != "smoke";
    let with_replay = matches!(
        stage.as_str(),
        "smoke" | "replay" | "preflight" | "all" | "pilot" | "diagnostic-pilot"
    );
    let with_timing = matches!(
        stage.as_str(),
        "smoke" | "timing" | "preflight" | "all" | "pilot" | "diagnostic-pilot"
    );
    let with_pilots = matches!(stage.as_str(), "pilot" | "diagnostic-pilot");
    let diagnostic_pilots = stage == "diagnostic-pilot";
    let workers = args
        .get(4)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(1usize);
    let draws = args
        .get(5)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(16usize);
    let repetitions = args
        .get(6)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(3usize);
    let seed = args
        .get(7)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(1337u64);
    ensure!(
        (1..=20).contains(&workers) && draws > 0 && repetitions > 0,
        "workers1..20 and positive draws/repetitions required"
    );
    ensure!(
        !with_pilots || workers == 20,
        "pilot preflight requires the planned20-worker timing; use smoke for bounded development"
    );
    let state = Path::new(manifest["state"].as_str().unwrap());
    ensure!(
        !output.canonicalize()?.starts_with(state.canonicalize()?),
        "output must be outside saved state"
    );
    let expected: BTreeMap<String, String> = serde_json::from_str(&fs::read_to_string(
        manifest["saved_state_hashes"].as_str().unwrap(),
    )?)?;
    let before = state_hashes(state)?;
    fs::write(
        output.join("state_before.json"),
        serde_json::to_string_pretty(&before)?,
    )?;
    ensure!(
        before == expected,
        "saved state differs from complete reviewed snapshot"
    );
    let mut report = json!({"status":"started","mode":modes,"stage":stage,"workers":workers,"draws_per_worker":draws,"repetitions":repetitions,"seed":seed,"reports":[]});
    let execution = (|| -> Result<()> {
        initialise()?;
        let at = Instant::now();
        let mut loaded = StateLoadOption::read_only(state).load()?;
        ensure!(loaded.is_read_only_state(), "read-only load required");
        loaded.state.activate_loaded_integrand_backends(false)?;
        report["cold_load_seconds"] = json!(at.elapsed().as_secs_f64());
        let process_id = loaded.state.resolve_process_ref(Some(&ProcessRef::Name(
            manifest["process_name"].as_str().unwrap().to_owned(),
        )))?;
        let integrand_name = manifest["integrand_name"].as_str().unwrap();
        let original = loaded
            .state
            .process_list
            .get_integrand_mut(process_id, integrand_name)?
            .get_settings()
            .clone();
        let packet = Path::new(&args[0])
            .parent()
            .ok_or_else(|| eyre!("manifest parent missing"))?;
        let mut requests = Vec::new();
        for name in ["raw_smoke_requests.json", "selected_host_requests.json"] {
            let source: Value = serde_json::from_str(&fs::read_to_string(packet.join(name))?)?;
            requests.extend(source["cases"].as_array().unwrap().iter().cloned());
        }
        // Every requested mode completes preflight before any integration starts.
        // The same loaded State remains alive while model/integrand borrows are
        // released at the end of each mode and before Integrate::run.
        for mode in &modes {
            let model = &loaded.state.model;
            let integrand = loaded
                .state
                .process_list
                .get_integrand_mut(process_id, integrand_name)?;
            let card = fs::read_to_string(
                manifest["cards"][mode]
                    .as_str()
                    .ok_or_else(|| eyre!("unknown mode {mode}"))?,
            )?;
            *integrand.get_mut_settings() = SetArgs::String { string: card }
                .merge_figment(Figment::from(Serialized::defaults(&original)))?
                .extract()?;
            let production_settings = integrand.get_settings().clone();
            let at = Instant::now();
            std::thread::scope(|scope| -> Result<()> {
                std::thread::Builder::new()
                    .stack_size(128 * 1024 * 1024)
                    .spawn_scoped(scope, || integrand.warm_up(model))?
                    .join()
                    .map_err(|_| eyre!("warmup panic"))?
            })?;
            let warmup_seconds = at.elapsed().as_secs_f64();
            let info = inventory(integrand, &manifest, mode)?;
            let mut mode_report =
                json!({"mode":mode,"warmup_seconds":warmup_seconds,"inventory":info});
            fs::write(
                output.join(format!("{mode}.settings.toml")),
                toml::to_string_pretty(&production_settings)?,
            )?;
            let mode_result = (|| -> Result<()> {
                if with_reference {
                    let reference =
                        match reference_stage(integrand, model, &manifest, full_reference) {
                            Ok(reference) => reference,
                            Err(error) => {
                                json!({"passed":false,"finite":false,"error":format!("{error:#}")})
                            }
                        };
                    let passed = reference["passed"] == true;
                    let finite = reference["finite"] == true;
                    mode_report["reference"] = reference;
                    ensure!(
                        finite,
                        "{mode}: numerical reference failure retained with completed/failed batch details"
                    );
                    ensure!(
                        !full_reference || passed || diagnostic_pilots,
                        "{mode}: finite reference missed acceptance bounds; dependent stages not started"
                    );
                }
                if with_replay {
                    // Events are a diagnostic retention overlay only. Restore
                    // the exact production card before timing; no event I/O is timed.
                    *integrand.get_mut_settings() = SetArgs::String {
                        string: fs::read_to_string(packet.join("cards/events_overlay.toml"))?,
                    }
                    .merge_figment(Figment::from(Serialized::defaults(&production_settings)))?
                    .extract()?;
                    let mut cases = requests
                        .iter()
                        .filter(|r| r["mode"] == mode.as_str())
                        .cloned()
                        .collect::<Vec<_>>();
                    // Direct-H uses the same stored bare controls; only its
                    // sampling card changes, never the physical input/metadata.
                    if cases.is_empty() {
                        cases = requests
                            .iter()
                            .filter(|r| r["mode"] == "optimized_lmb")
                            .cloned()
                            .map(|mut r| {
                                r["mode"] = json!(mode);
                                r
                            })
                            .collect();
                    }
                    let rows = std::thread::scope(|scope| -> Result<_> {
                        std::thread::Builder::new()
                            .stack_size(128 * 1024 * 1024)
                            .spawn_scoped(scope, || {
                                integrand.warm_up(model)?;
                                replay(integrand, model, &cases, &info, true)
                            })?
                            .join()
                            .map_err(|_| eyre!("replay worker panicked"))?
                    })?;
                    mode_report["replay_case_count"] = json!(rows.len());
                    mode_report["replays"] = json!(rows);
                    let passed = !rows.is_empty()
                        && rows.iter().all(|r| {
                            r["evaluation"]["valid"] == true && r["six_distinct_cut_events"] == true
                        });
                    mode_report["replay_passed"] = json!(passed);
                    // Exhausted/invalid points remain in the report. No second
                    // retry engine or changed source is introduced to hide them.
                    ensure!(
                        passed,
                        "{mode}: raw/selected replay failed; dependent timing not started"
                    );
                    *integrand.get_mut_settings() = production_settings.clone();
                }
                if with_timing {
                    let mut hard = Vec::new();
                    let labels = info["channels"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .filter_map(|c| c["label"].as_str())
                        .collect::<Vec<_>>();
                    let primary = match mode.as_str() {
                        "joint_hz_plus_lmb" => "direct_joint_HZ",
                        "direct_h_p2" => "direct_H_p2",
                        _ => labels[0],
                    };
                    ensure!(
                        labels.contains(&primary),
                        "requested hard channel {primary} missing"
                    );
                    let full_lmb = info["channels"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .find(|channel| channel["is_lmb"] == true)
                        .and_then(|channel| channel["label"].as_str());
                    for case in requests
                        .iter()
                        .filter(|r| r["mode"] == "optimized_lmb" && r["forced_arb"] == false)
                    {
                        for selected in [false, true] {
                            let mut r = case.clone();
                            r["mode"] = json!(mode);
                            let selected_label = if mode == "joint_hz_plus_lmb"
                                && case["source_case"].as_str().unwrap().starts_with("soft22")
                            {
                                full_lmb.ok_or_else(|| {
                                    eyre!(
                                        "soft22 control requires explicit full-support LMB sibling"
                                    )
                                })?
                            } else {
                                primary
                            };
                            r["channel_name"] = if selected {
                                json!(selected_label)
                            } else {
                                Value::Null
                            };
                            r["selection_scope"] = json!(if mode == "joint_hz_plus_lmb"
                                && selected_label != primary
                            {
                                "explicit full-support soft control; no compact-joint coverage claim"
                            } else {
                                "explicit named chart; support must be established by the actual inverse"
                            });
                            r["name"] = json!(format!(
                                "{}_hard_{}",
                                case["source_case"].as_str().unwrap(),
                                if selected { "selected" } else { "bare" }
                            ));
                            hard.push(r);
                        }
                    }
                    ensure!(
                        hard.len() == 6,
                        "three stored hard points in both bare/selected modes required"
                    );
                    let at = Instant::now();
                    let rows = time_sources(
                        integrand,
                        model,
                        workers,
                        draws,
                        repetitions,
                        seed,
                        &hard,
                        &info,
                    )?;
                    mode_report["timing_batch_wall_seconds"] = json!(at.elapsed().as_secs_f64());
                    let representative = timing_summary(&rows, "rows", None);
                    let hard_selected = timing_summary(&rows, "hard_rows", Some(true));
                    let hard_bare = timing_summary(&rows, "hard_rows", Some(false));
                    let maxima = timing_summary(&rows, "maximum_rows", None);
                    let summaries = [&representative, &hard_selected, &hard_bare, &maxima];
                    let finite = summaries.iter().all(|s| {
                        s["invalid_calls"] == 0
                            && s["warmup_failed"] == false
                            && s["measured_calls"].as_u64().unwrap_or(0) > 0
                    });
                    let cap = workers == 20
                        && [&representative, &hard_selected, &maxima]
                            .iter()
                            .all(|s| s["representative_cap_passed"] == true);
                    mode_report["timing_summary"] = representative;
                    mode_report["hard_selected_timing"] = hard_selected;
                    mode_report["hard_bare_timing"] = hard_bare;
                    mode_report["representative_maximum_timing"] = maxima;
                    mode_report["timing_finite"] = json!(finite);
                    mode_report["timing_cap_passed"] = json!(cap);
                    mode_report["timing_workers"] = json!(rows);
                    ensure!(
                        finite,
                        "{mode}: exhausted/invalid timing result retained; no pilot can hide numerical failures"
                    );
                }
                Ok(())
            })();
            *integrand.get_mut_settings() = production_settings;
            let accepted = mode_report["reference"]["passed"] == true
                && mode_report["replay_passed"] == true
                && mode_report["timing_cap_passed"] == true;
            mode_report["preflight_accepted"] = json!(accepted);
            mode_report["error"] = mode_result
                .as_ref()
                .err()
                .map(|e| json!(format!("{e:#}")))
                .unwrap_or(Value::Null);
            fs::write(
                output.join(format!("{mode}.json")),
                serde_json::to_string_pretty(&mode_report)?,
            )?;
            report["reports"].as_array_mut().unwrap().push(mode_report);
            mode_result?;
        }
        if with_replay {
            // Decimal-only postprocessing of retained current-version values;
            // this script never loads a GammaLoop state or evaluates physics.
            let audit = Command::new("python")
                .arg(packet.join("audit_physical.py"))
                .arg(&args[0])
                .arg(output)
                .output()?;
            let audit_path = output.join("physics_audit.json");
            report["physics_audit"] = if audit_path.exists() {
                serde_json::from_str(&fs::read_to_string(&audit_path)?)?
            } else {
                json!({"error":"audit produced no result","stdout":String::from_utf8_lossy(&audit.stdout),"stderr":String::from_utf8_lossy(&audit.stderr)})
            };
            ensure!(
                audit.status.success(),
                "current native physical-accuracy checks failed; see physics_audit.json"
            );
        }
        let preflight_accepted = report["reports"]
            .as_array()
            .unwrap()
            .iter()
            .all(|r| r["preflight_accepted"] == true);
        report["preflight_accepted"] = json!(preflight_accepted);
        if with_pilots {
            ensure!(
                preflight_accepted || diagnostic_pilots,
                "pilot preflight missed acceptance; use explicit diagnostic-pilot stage only for bounded investigation"
            );
            let mut cli = CLISettings::default();
            cli.session.read_only_state = true;
            cli.state.folder = state.to_path_buf();
            let mut pilots = Vec::new();
            for pilot_seed in manifest["paired_seeds"].as_array().unwrap() {
                for mode in &modes {
                    let workspace = output
                        .join("pilot")
                        .join(format!("{mode}_seed{}", pilot_seed.as_u64().unwrap()))
                        .join("workspace");
                    ensure!(
                        !workspace.exists(),
                        "refusing to resume/overwrite pilot {}",
                        workspace.display()
                    );
                    {
                        let integrand = loaded
                            .state
                            .process_list
                            .get_integrand_mut(process_id, integrand_name)?;
                        *integrand.get_mut_settings() = SetArgs::String {
                            string: fs::read_to_string(manifest["cards"][mode].as_str().unwrap())?,
                        }
                        .merge_figment(Figment::from(Serialized::defaults(&original)))?
                        .extract()?;
                        let settings = integrand.get_mut_settings();
                        settings.integrator.seed = pilot_seed.as_u64().unwrap();
                        ensure!(
                            settings.integrator.n_start == 2048
                                && settings.integrator.n_max == 2048
                                && settings.integrator.n_increase == 0,
                            "only planned fresh2048-point single-iteration pilots supported"
                        );
                    }
                    let command = Integrate {
                        process: vec![ProcessRef::Name(
                            manifest["process_name"].as_str().unwrap().to_owned(),
                        )],
                        integrand_name: vec![integrand_name.to_owned()],
                        n_cores: Some(20),
                        workspace_path: Some(workspace.clone()),
                        batch_size: Some(256),
                        batch_timing: 0.0,
                        renderer: RendererOption::Tabled,
                        show_phase: ShowPhaseOption::Both,
                        show_max_weight_info: true,
                        show_max_weight_info_for_discrete_bins: true,
                        write_results_for_each_iteration: true,
                        no_stream_updates: true,
                        no_stream_iterations: true,
                        ..Integrate::default()
                    };
                    let at = Instant::now();
                    // Existing Integrate owns slot warming, worker dispatch,
                    // accumulation, absolute estimates, maxima and checkpoints.
                    let result = std::thread::scope(|scope| -> Result<_> {
                        std::thread::Builder::new()
                            .stack_size(128 * 1024 * 1024)
                            .spawn_scoped(scope, || command.run(&mut loaded.state, &cli))?
                            .join()
                            .map_err(|_| eyre!("pilot worker panicked"))?
                    });
                    let wall = at.elapsed().as_secs_f64();
                    let (record, finite) = match result {
                        Ok(value) => {
                            let slot = value
                                .result
                                .single_slot()
                                .ok_or_else(|| eyre!("single pilot slot required"))?;
                            let scalars = [
                                slot.integral.result.re.0,
                                slot.integral.result.im.0,
                                slot.integral.error.re.0,
                                slot.integral.error.im.0,
                                slot.absolute.integral.result.re.0,
                                slot.absolute.integral.result.im.0,
                                slot.absolute.integral.error.re.0,
                                slot.absolute.integral.error.im.0,
                            ];
                            let finite = slot.integral.neval == 2048
                                && slot.integration_statistics.nan_or_unstable_percentage == 0.0
                                && scalars.iter().all(|v| v.is_finite());
                            (json!({"output":value,"finite_complete":finite}), finite)
                        }
                        Err(error) => (
                            json!({"error":format!("{error:#}"),"finite_complete":false}),
                            false,
                        ),
                    };
                    let row = json!({"mode":mode,"seed":pilot_seed,"workspace":workspace,"elapsed_seconds":wall,"result":record,
                        "diagnostic_only":diagnostic_pilots,"preflight_accepted":preflight_accepted,
                        "scope":"fresh first iteration; equal seeds pair schedules, not physical points across maps; no second integration/statistics engine"});
                    pilots.push(row);
                    report["pilots"] = json!(pilots);
                    fs::write(
                        output.join("pilots.json"),
                        serde_json::to_string_pretty(&pilots)?,
                    )?;
                    ensure!(
                        finite,
                        "pilot numerical/completion failure retained; remaining pilots not started"
                    );
                    let integrand = loaded
                        .state
                        .process_list
                        .get_integrand_mut(process_id, integrand_name)?;
                    inventory(integrand, &manifest, mode)?;
                }
            }
        }

        Ok(())
    })();
    // Check persistence even on a propagated numerical/inventory failure.
    let after = state_hashes(state)?;
    fs::write(
        output.join("state_after.json"),
        serde_json::to_string_pretty(&after)?,
    )?;
    report["saved_state_unchanged"] = json!(before == after);
    report["error"] = execution
        .as_ref()
        .err()
        .map(|e| json!(format!("{e:#}")))
        .unwrap_or(Value::Null);
    let invalid = report["reports"].as_array().unwrap().iter().any(|r| {
        r.get("replay_passed") == Some(&Value::Bool(false))
            || r.get("timing_finite") == Some(&Value::Bool(false))
    });
    report["status"] = json!(if execution.is_ok() && before == after && !invalid {
        if diagnostic_pilots {
            "completed_diagnostic_pilot"
        } else {
            "completed_requested_stages"
        }
    } else {
        "failed"
    });
    report["scope"] = json!({"finite_smoke_is_not_normalization":true,"runtime_cap":"representative, stored hard selected/bare and replayed representative maxima, each reported separately",
        "precision_scope":"fresh first iteration max_eval=0; trained-iteration rescue history not certified",
        "maps_scope":"completed-point hosted joint; CT-star density cure remains pending","gain_claim":false});
    fs::write(
        output.join("summary.json"),
        serde_json::to_string_pretty(&report)?,
    )?;
    execution?;
    ensure!(before == after, "saved state changed");
    ensure!(!invalid, "bounded evaluations failed; see retained reports");
    Ok(())
}
