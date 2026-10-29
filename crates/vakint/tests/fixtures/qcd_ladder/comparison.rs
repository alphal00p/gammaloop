//! Reproducible physical-input comparisons, separate from scalar preparation.

use super::super::vakint::{
    AvailablePrecision, ConstantSymbol, EvaluationMethod, EvaluationOrder, LoopNormalizationFactor,
    MasterSymbol, OfficialExternalDataOptions, RustRedEvaluationOptions, RustRedMasterMode,
    RustRedMonitorOptions, RustRedMonitorScope, RustRedReductionStrategy, Vakint, VakintMaster,
    VakintSettings, vakint_parse, vakint_symbol,
};
use super::{LadderFixture, Projection, Result};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    path::{Path, PathBuf},
    time::Instant,
};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::{integer::Integer, rational::Rational};

pub struct RunOptions {
    pub strategy: RustRedReductionStrategy,
    pub artifacts: Option<PathBuf>,
    pub monitor: RustRedMonitorOptions,
    pub numerical_reference: bool,
}

pub struct CompletedRun {
    pub value: Atom,
    pub report: serde_json::Value,
}

/// Capture the retained call only. Symbolic expansion and numerical references
/// happen after this closure has returned and its instrumentation has closed.
fn capture_retained_call<T, E>(
    evaluate: impl FnOnce() -> std::result::Result<T, E>,
) -> (std::result::Result<T, E>, serde_json::Value) {
    #[cfg(feature = "experimental-rustred")]
    {
        let capture =
            super::super::vakint::rustred_evaluation::experimental::capture_reduction_timing(
                evaluate,
            );
        let timing = capture.timing;
        (
            capture.result,
            serde_json::json!({
                "reduction_math_seconds":timing.reduction_math.map(|duration| duration.as_secs_f64()),
                "reduction_math_timing":{
                    "scope_id":timing.scope_id,
                    "complete":timing.complete,
                    "gross_seconds":timing.gross.as_secs_f64(),
                    "measured_owner_preparation_seconds":timing.measured_owner_preparation.as_secs_f64(),
                    "interval_count":timing.interval_count,
                    "incompleteness_reason":timing.incompleteness_reason,
                },
            }),
        )
    }
    #[cfg(not(feature = "experimental-rustred"))]
    {
        (
            evaluate(),
            serde_json::json!({
                "reduction_math_seconds":null,
                "reduction_math_timing":{
                    "scope_id":"vakint.rustred.reduction-math.v1",
                    "complete":false,
                    "gross_seconds":null,
                    "measured_owner_preparation_seconds":null,
                    "interval_count":0,
                    "incompleteness_reason":"build lacks the experimental-rustred timing capture feature",
                },
            }),
        )
    }
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn write_identity(directory: &Path, name: &str, text: &str) -> Result<serde_json::Value> {
    std::fs::write(directory.join(name), text)?;
    Ok(serde_json::json!({"file":name,"bytes":text.len(),"sha256":sha256(text.as_bytes())}))
}

fn process_memory() -> serde_json::Value {
    let status = std::fs::read_to_string("/proc/self/status").ok();
    let bytes = |name: &str| {
        status.as_ref().and_then(|text| {
            text.lines().find_map(|line| {
                let number = line
                    .strip_prefix(name)?
                    .split_whitespace()
                    .next()?
                    .parse::<u64>()
                    .ok()?;
                number.checked_mul(1024)
            })
        })
    };
    serde_json::json!({"rss_bytes":bytes("VmRSS:"),"process_peak_rss_bytes":bytes("VmHWM:"),
        "scope":"Linux whole-process RSS; peak includes preparation/loading and earlier stages"})
}

fn finite_cap(availability: AvailablePrecision, cap: &mut Option<u32>) {
    if let AvailablePrecision::FixedPrecision(value) = availability {
        *cap = Some(cap.map_or(value, |old| old.min(value)));
    }
}

/// Inspect only identities which survive exact symbolic master expansion.
/// The inspection target is deliberately low: source capability, rather than
/// the accuracy of these inspection values, selects the subsequent evaluation.
fn reference_precision(engine: &Vakint, expression: &Atom) -> Result<(u32, serde_json::Value)> {
    let registered: BTreeSet<_> = engine
        .constant_symbols()?
        .into_iter()
        .map(|symbol| symbol.as_atom().clone())
        .collect();
    let mut identities = BTreeSet::new();
    expression.visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            if function.get_symbol() == vakint_symbol!("RustRedFiniteLiteral") {
                if let Some(coefficient) = function.iter().next() {
                    identities.insert(coefficient.to_owned());
                }
                return false;
            }
            if function.get_symbol() == vakint_symbol!("RustRedCoefficient")
                || function.get_symbol() == vakint_symbol!("RustRedConstant")
            {
                identities.insert(part.to_owned());
                return false;
            }
        }
        if registered.contains(&part.to_owned()) {
            identities.insert(part.to_owned());
            return false;
        }
        true
    });
    let mut masters = BTreeMap::<Atom, VakintMaster>::new();
    let mut cap = None;
    let mut sources = Vec::new();
    for identity in identities {
        let (available, version, source_sha256) = if let AtomView::Fun(function) =
            identity.as_view()
            && function.get_symbol() == vakint_symbol!("RustRedCoefficient")
        {
            let arguments = function.iter().collect::<Vec<_>>();
            let [_, master, order] = arguments.as_slice() else {
                return Err("malformed coefficient identity in physical reference".into());
            };
            let master = (*master).to_owned();
            if let std::collections::btree_map::Entry::Vacant(entry) = masters.entry(master.clone())
            {
                entry.insert(engine.get_master(&MasterSymbol::try_from(master.clone())?, 16)?);
            }
            let power = i64::try_from(*order)?;
            let coefficient = masters[&master]
                .coefficients()
                .get(&power)
                .ok_or("symbolic reference uses a missing master coefficient")?;
            (
                coefficient.available_precision(),
                coefficient.provenance().version().to_owned(),
                coefficient.provenance().source_sha256().to_owned(),
            )
        } else {
            let constant = engine.get_constant(&ConstantSymbol::try_from(identity.clone())?, 16)?;
            (
                constant.available_precision(),
                constant.provenance().version().to_owned(),
                constant.provenance().source_sha256().to_owned(),
            )
        };
        finite_cap(available, &mut cap);
        sources.push(serde_json::json!({"identity":identity.to_canonical_string(),"available_precision":format!("{available:?}"),"data_version":version,"source_sha256":source_sha256}));
    }
    let requested = cap.unwrap_or(300);
    Ok((
        requested,
        serde_json::json!({"requested_decimal_precision":requested,
        "selection":if cap.is_some() {"minimum validated fixed precision among surviving dependencies"} else {"all surviving dependencies arbitrary precision; explicit 300-digit comparison target, not a source maximum"},
        "surviving_dependencies":sources}),
    ))
}

pub fn run(
    prepared: &LadderFixture,
    directory: &Path,
    options: &RunOptions,
) -> Result<CompletedRun> {
    if prepared.loop_count == 5 && options.strategy != RustRedReductionStrategy::Ordered {
        return Err("five-loop physical benchmarks permit only ordered reduction; use the preserved historical DFS result for comparison".into());
    }
    if prepared.loop_count == 5 && options.monitor.enabled() {
        RustRedMonitorScope::run(&options.monitor, |scope| {
            run_scoped(prepared, directory, options, Some(scope))
        })
    } else {
        run_scoped(prepared, directory, options, None)
    }
}

#[cfg(test)]
#[test]
fn five_loop_benchmarks_reject_both_depth_first_strategies_before_loading() {
    let prepared = LadderFixture {
        loop_count: 5,
        input: Atom::one(),
        preparation_seconds: 0.0,
        scalar_terms: 1,
        contributions: 1,
        nonzero_conditions: Vec::new(),
        projection: Projection::ZeroMomentumMetric,
        timings: super::PreparationTimings::default(),
        prepared_input: None,
    };
    let directory = tempfile::tempdir().unwrap();
    // An invalid local manifest prevents network access even if this guard
    // regresses and the loader is accidentally reached.
    std::fs::write(directory.path().join("manifest.json"), "{}").unwrap();
    let output = directory.path().join("must-not-create");
    for strategy in [
        RustRedReductionStrategy::DepthFirst,
        RustRedReductionStrategy::DepthFirstWithoutMemoization,
    ] {
        let error = run(
            &prepared,
            &output,
            &RunOptions {
                strategy,
                artifacts: Some(directory.path().to_owned()),
                monitor: RustRedMonitorOptions::default(),
                numerical_reference: false,
            },
        )
        .err()
        .expect("a five-loop DFS benchmark must fail before any work");
        assert!(error.to_string().contains("permit only ordered reduction"));
        assert!(!output.exists());
    }
}

fn run_scoped(
    prepared: &LadderFixture,
    directory: &Path,
    options: &RunOptions,
    mut scope: Option<&mut RustRedMonitorScope<'_, '_>>,
) -> Result<CompletedRun> {
    std::fs::create_dir_all(directory)?;
    let total = Instant::now();
    let loading = Instant::now();
    if let Some(scope) = scope.as_deref_mut() {
        scope.phase("data_loading");
    }
    let mut engine = Vakint::new()?;
    let data = if prepared.loop_count >= 4 || options.artifacts.is_some() {
        let loaded = engine.load_official_external_data(
            prepared.loop_count,
            &OfficialExternalDataOptions {
                local_artifacts_path: options.artifacts.clone(),
                load_prepared_caches: false,
                ..Default::default()
            },
        )?;
        serde_json::json!({"version":loaded.version,"manifest_sha256":sha256(&std::fs::read(loaded.path.join("manifest.json"))?),"path":loaded.path})
    } else {
        let manifest = include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/manifest.json"));
        let metadata: serde_json::Value = serde_json::from_slice(manifest)?;
        serde_json::json!({"source":"vendored offline 1L–3L runtime registry", "version":metadata["version"], "manifest_sha256":sha256(manifest)})
    };
    let loading_seconds = loading.elapsed().as_secs_f64();
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        form_exe_path: "/ordered-rung-comparison-must-not-use-form".into(),
        use_dot_product_notation: true,
        number_of_terms_in_epsilon_expansion: prepared.loop_count as i64 + 1,
        ..Default::default()
    };
    settings.rustred_application_options.reduction_strategy = options.strategy;
    // Native owner-preparation diagnostics sum worker wall times. Keep one
    // worker so subtracting measured preparation from math spans is valid.
    settings.rustred_application_options.preparation_workers = 1;
    // Bind the fixture's exact native trace dimension identically for both
    // public evaluators; the raw prepared scalar payload remains unchanged.
    let native_d = super::dimension()?;
    let physical_d =
        Atom::num(4) - Atom::num(2) * Atom::var(vakint_symbol!(&settings.epsilon_symbol));
    let evaluation_input = prepared
        .input
        .replace(native_d.clone())
        .with(physical_d.clone());
    let reduction = Instant::now();
    let (retained_result, math_timing) = capture_retained_call(|| -> Result<_> {
        if prepared.loop_count == 5 {
            let output = if let Some(scope) = scope.as_deref_mut() {
                engine.evaluate_rustred_symbolic_in_scope(
                    &settings,
                    evaluation_input.as_view(),
                    scope,
                )?
            } else {
                engine.evaluate_rustred_symbolic_with_monitor(
                    &settings,
                    evaluation_input.as_view(),
                    &options.monitor,
                )?
            };
            Ok((output.value, Some(output.nonzero_conditions)))
        } else {
            Ok((
                engine.evaluate_with_monitor(
                    &settings,
                    evaluation_input.as_view(),
                    &options.monitor,
                )?,
                None,
            ))
        }
    });
    let (value, guards) = retained_result?;
    let reduction_seconds = reduction.elapsed().as_secs_f64();
    let retained_memory = process_memory();
    if let Some(scope) = scope.as_deref_mut() {
        scope.expression_phase("retained_expression_serialization", value.as_view());
    }
    let retained = write_identity(
        directory,
        "retained.symbolica",
        &value.to_canonical_string(),
    )?;
    let catalogue = engine
        .master_symbols()?
        .into_iter()
        .map(|master| master.as_atom().clone())
        .collect::<BTreeSet<_>>();
    if let Some(scope) = scope.as_deref_mut() {
        scope.expression_phase("catalogue_classification", value.as_view());
    }
    let mut masters = BTreeSet::new();
    value.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == vakint_symbol!("RustRedMaster")
        {
            masters.insert(part.to_owned());
            return false;
        }
        true
    });
    let unknown = masters
        .difference(&catalogue)
        .map(AtomCore::to_canonical_string)
        .collect::<Vec<_>>();
    let guards = match guards {
        Some(guards) => {
            let native_d = super::dimension()?;
            let physical_d =
                Atom::num(4) - Atom::num(2) * Atom::var(vakint_symbol!(&settings.epsilon_symbol));
            let all =
                guards
                    .into_iter()
                    .chain(prepared.nonzero_conditions.iter().map(|condition| {
                        condition.replace(native_d.clone()).with(physical_d.clone())
                    }))
                    .collect::<BTreeSet<_>>();
            if all.iter().any(Atom::is_zero) {
                return Err(
                    "physical reduction returned an identically zero nonzero condition".into(),
                );
            }
            let mut text = String::new();
            let mut progress = scope.as_deref_mut().map(|scope| {
                scope.work_phase("guard_serialization", all.len(), "guards serialized")
            });
            for (index, guard) in all.iter().enumerate() {
                if index > 0 {
                    text.push('\n');
                }
                text.push_str(&guard.to_canonical_string());
                if let Some(progress) = &mut progress {
                    progress.completed(index + 1);
                }
            }
            serde_json::json!({"payload":write_identity(directory,"nonzero-conditions.symbolica",&text)?,"count":all.len(),"scope":"sufficient conditions, including scalar preparation; exact exceptional-domain equality is not inferred from equal results"})
        }
        None => {
            serde_json::json!({"scope":"existing lower-loop public Atom API does not expose its internal guard set","preparation":prepared.nonzero_conditions.iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>()})
        }
    };
    let intrinsic_memoization = options.strategy == RustRedReductionStrategy::DepthFirst;
    let mut report = serde_json::json!({"schema":"vakint.ordered-rung-comparison.v1",
        "loop_count":prepared.loop_count,"gluon_rungs":prepared.loop_count-1,"projection":prepared.projection,
        "strategy":format!("{:?}",options.strategy),"persistent_cache":false,
        "intrinsic_memoization":intrinsic_memoization,
        "intrinsic_memoization_scope":"primary depth-first completed-result reuse; DFS-no-memo also bypasses dependency/product result memos, while Ordered checked opaque factorization may use a fresh isolated internal memo",
        "cache_policy":{"prepared_cache_loaded":false,"previous_run_cache_restored":false,
            "depth_first_intrinsic_memoization":intrinsic_memoization,
            "ordered_result_cache":"not supported"},
        "worker_policy":{"owner_preparation_workers":1,
            "scope":"single owner preparation worker required for native timing subtraction"},
        "input_sha256":sha256(prepared.input.to_canonical_string().as_bytes()),
        "evaluated_input_sha256":sha256(evaluation_input.to_canonical_string().as_bytes()),
        "dimension_binding":{"source":native_d.to_canonical_string(),"target":physical_d.to_canonical_string(),"scope":"exact fixture trace dimension; all public evaluations share this binding"},
        "preparation":{"reused_authenticated_bundle":prepared.prepared_input.is_some(),
            "historical_seconds":prepared.preparation_seconds,"timings":prepared.timings,
            "current_run_seconds":if prepared.prepared_input.is_some() {0.0} else {prepared.preparation_seconds},
            "timing_scope":"scalar preparation precedes retained evaluation; loaded-bundle historical work is not added to current timing"},
        "contributions":prepared.contributions,"data":data,"loading_seconds":loading_seconds,
        "retained_reduction_seconds":reduction_seconds,"retained":retained,"guards":guards,
        "reduction_math_seconds":math_timing["reduction_math_seconds"],
        "reduction_math_timing":math_timing["reduction_math_timing"],
        "retained_reduction_seconds_scope":"complete retained public call, including frontend and finalization; separate from reduction math",
        "reduction_completed":unknown.is_empty(),"catalogue_closure":{"closed":unknown.is_empty(),"master_count":masters.intersection(&catalogue).count(),"unknown_count":unknown.len(),"unknown":unknown},
        "memory_after_retained_reduction":retained_memory,"monitoring_enabled":options.monitor.enabled(),
        "normalization":"MSbar","numerical_reference":null});
    report["epsilon_symbol"] = settings.epsilon_symbol.clone().into();
    report["mu_r_sq_symbol"] = settings.mu_r_sq_symbol.clone().into();
    if options.numerical_reference {
        if prepared.loop_count > 4 {
            return Err("five-loop catalogue coefficients are uncomputed; numerical references are only available through four loops".into());
        }
        let reference = Instant::now();
        settings.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            })]);
        let symbolic = engine.evaluate(&settings, evaluation_input.as_view())?;
        let model = super::Model::from_json(include_str!("../qcd.json"))?;
        let symbolic = model.expand_couplings(&symbolic);
        // Exact reference-point substitutions preserve source accuracy and do
        // not introduce an f64 parameter cap. Coupling G=1 is stated explicitly.
        let symbolic = symbolic
            .replace(Atom::var(vakint_symbol!("mUV")))
            .with(Atom::one())
            .replace(Atom::var(vakint_symbol!(&settings.mu_r_sq_symbol)))
            .with(Atom::one())
            .replace(Atom::var(symbolica::symbol!("UFO::G")))
            .with(Atom::one());
        let (requested, selection) = reference_precision(&engine, &symbolic)?;
        settings.run_time_decimal_precision = requested;
        let (numerical, statistical) = engine.numerical_evaluation(
            &settings,
            symbolic.as_view(),
            &HashMap::default(),
            &HashMap::default(),
            None,
        )?;
        if statistical.is_some() {
            return Err(
                "deterministic rung reduction unexpectedly returned statistical uncertainty".into(),
            );
        }
        let coefficients = numerical.coefficients().iter().map(|(power, coefficient)| {
            let precision = coefficient.precision();
            serde_json::json!({"epsilon_power":power,"real":coefficient.value().re.to_string(),"imaginary":coefficient.value().im.to_string(),
                "real_exact_rational":coefficient.value().re.to_rational().to_string(),"imaginary_exact_rational":coefficient.value().im.to_rational().to_string(),
                "real_significant_digits":precision.real().significant_digits(),"imaginary_significant_digits":precision.imaginary().significant_digits(),
                "real_absolute_error":precision.real().absolute_error().to_string(),"imaginary_absolute_error":precision.imaginary().absolute_error().to_string(),
                "real_error_exact_rational":precision.real().absolute_error().to_rational().to_string(),"imaginary_error_exact_rational":precision.imaginary().absolute_error().to_rational().to_string(),
                "real_exact_zero":precision.real().exact_zero(),"imaginary_exact_zero":precision.imaginary().exact_zero(),
                "real_error_quality":format!("{:?}",precision.real().quality()),"imaginary_error_quality":format!("{:?}",precision.imaginary().quality())})
        }).collect::<Vec<_>>();
        report["numerical_reference"] = serde_json::json!({"selection":selection,"coefficients":coefficients,
            "point":{"mUV":1,"mu_r_squared":1,"UFO_G":1},"coupling_definitions":"exact expansion through the authenticated native QCD model (GC_11=i*G)","epsilon_range":"all available poles through epsilon^0; remainder O(epsilon)",
            "seconds":reference.elapsed().as_secs_f64(),"scope":"separate Symbolic-mode public evaluation plus tracked substitution; excluded from retained comparison time",
            "warnings":numerical.warnings().iter().map(ToString::to_string).collect::<Vec<_>>()});
    }
    report["completed_seconds"] = total.elapsed().as_secs_f64().into();
    report["completed_seconds_scope"] = "loading, retained evaluation, artifact serialization and optional numerical reference; scalar input preparation is recorded separately".into();
    if let Some(scope) = scope {
        scope.phase("report_serialization");
    }
    std::fs::write(
        directory.join("report.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    if !unknown.is_empty() {
        return Err(format!(
            "physical reference retained undeclared integrals: {}",
            unknown.join("; ")
        )
        .into());
    }
    Ok(CompletedRun { value, report })
}

fn rational_field(row: &serde_json::Value, field: &str) -> Result<Rational> {
    let text = row[field]
        .as_str()
        .ok_or_else(|| format!("missing exact reference field {field}"))?;
    let (numerator, denominator) = text.split_once('/').unwrap_or((text, "1"));
    let numerator: Integer = numerator
        .parse()
        .map_err(|error| format!("invalid reference numerator: {error}"))?;
    let denominator: Integer = denominator
        .parse()
        .map_err(|error| format!("invalid reference denominator: {error}"))?;
    if denominator == 0 {
        return Err("reference has a zero denominator".into());
    }
    Ok(Rational::from((numerator, denominator)))
}

/// Only exp(±γ)>0 for the registered Euler constant is flattened. Neither
/// stored artifacts nor arbitrary complex exponentials are rewritten.
pub fn canonical_retained(value: &Atom) -> Result<Atom> {
    let gamma = Atom::var(vakint_symbol!("EulerGamma"));
    let bases = [
        (vakint_parse!("exp(-EulerGamma)")?, -gamma.clone()),
        (vakint_parse!("exp(EulerGamma)")?, gamma),
    ];
    let mut result = value.clone();
    loop {
        let next = result.replace_map(|part, _, output| {
            if let AtomView::Pow(power) = part
                && let Some((_, logarithm)) = bases
                    .iter()
                    .find(|(base, _)| power.get_base() == base.as_view())
            {
                **output = Atom::var(symbolica::atom::Symbol::E)
                    .pow(logarithm.clone() * power.get_exp().to_owned());
            }
        });
        if next == result {
            return Ok(result);
        }
        result = next;
    }
}

/// Convert the preserved pySecDec result with the same public normalization
/// definitions as evaluation, before any Laurent truncation. The second Atom
/// is the exact applied factor for the conversion's provenance receipt.
pub fn convert_historical_retained(
    value: &Atom,
    report: &serde_json::Value,
) -> Result<(Atom, Atom)> {
    if report["normalization"] != "MSbar" {
        return Err("historical conversion requires an MSbar target report".into());
    }
    let loop_count = report["loop_count"]
        .as_u64()
        .ok_or("target report has no loop count")?;
    let mut settings = VakintSettings {
        epsilon_symbol: report["epsilon_symbol"]
            .as_str()
            .ok_or("target report has no epsilon symbol")?
            .to_owned(),
        mu_r_sq_symbol: report["mu_r_sq_symbol"]
            .as_str()
            .ok_or("target report has no renormalization-scale symbol")?
            .to_owned(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        ..Default::default()
    };
    let target = settings.get_integral_normalization_factor_atom()?;
    settings.integral_normalization_factor = LoopNormalizationFactor::pySecDec;
    let source = settings.get_integral_normalization_factor_atom()?;
    let factor = (target / source)
        .replace(Atom::var(vakint_symbol!("n_loops")))
        .with(Atom::num(loop_count))
        .replace(Atom::var(vakint_symbol!("log_mu_sq")))
        .with(
            symbolica::atom::Symbol::LOG.call(Atom::var(vakint_symbol!(&settings.mu_r_sq_symbol))),
        );
    let native_d = super::dimension()?;
    let physical_d =
        Atom::num(4) - Atom::num(2) * Atom::var(vakint_symbol!(&settings.epsilon_symbol));
    if report["dimension_binding"]["source"] != native_d.to_canonical_string()
        || report["dimension_binding"]["target"] != physical_d.to_canonical_string()
    {
        return Err("target report has an inconsistent physical dimension binding".into());
    }
    Ok((
        value.replace(native_d).with(physical_d) * factor.clone(),
        factor,
    ))
}

#[cfg(test)]
#[test]
fn comparison_power_identity_is_restricted_to_the_positive_euler_base() -> Result<()> {
    let nested = vakint_parse!("exp(-EulerGamma)^p")?;
    assert_eq!(
        canonical_retained(&nested)?,
        vakint_parse!("exp(-EulerGamma*p)")?
    );
    assert_eq!(
        canonical_retained(&vakint_parse!("exp(EulerGamma)^p")?)?,
        vakint_parse!("exp(EulerGamma*p)")?
    );
    let complex = vakint_parse!("exp(𝑖*𝜋)^p")?;
    assert_eq!(canonical_retained(&complex)?, complex);
    let unknown = vakint_parse!("exp(x)^p")?;
    assert_eq!(canonical_retained(&unknown)?, unknown);
    Ok(())
}

/// Compare serialized binary-rational centers and error radii exactly. Decimal
/// strings are display values only and never contribute formatting uncertainty.
pub fn compare_numerical_references(
    left: &serde_json::Value,
    right: &serde_json::Value,
) -> Result<()> {
    if left["point"] != right["point"] || left["epsilon_range"] != right["epsilon_range"] {
        return Err("physical numerical references use different points or epsilon ranges".into());
    }
    let rows = |reference: &serde_json::Value| -> Result<BTreeMap<i64, serde_json::Value>> {
        let coefficients = reference["coefficients"]
            .as_array()
            .ok_or("numerical reference has no coefficient inventory")?;
        let mut result = BTreeMap::new();
        for coefficient in coefficients {
            let order = coefficient["epsilon_power"]
                .as_i64()
                .ok_or("numerical reference has an invalid epsilon order")?;
            if result.insert(order, coefficient.clone()).is_some() {
                return Err("numerical reference repeats an epsilon order".into());
            }
        }
        Ok(result)
    };
    let left = rows(left)?;
    let right = rows(right)?;
    if left.keys().ne(right.keys()) {
        return Err("physical numerical references contain different Laurent orders".into());
    }
    for (order, left) in left {
        let right = &right[&order];
        for component in ["real", "imaginary"] {
            for row in [&left, right] {
                let exact_zero = row[format!("{component}_exact_zero")]
                    .as_bool()
                    .ok_or("missing component exact-zero status")?;
                if !exact_zero
                    && row[format!("{component}_significant_digits")]
                        .as_u64()
                        .unwrap_or(0)
                        < 64
                {
                    return Err(format!("physical reference at epsilon^{order} has fewer than 64 meaningful digits; inspect its source/rounding warnings").into());
                }
            }
            let l = rational_field(&left, &format!("{component}_exact_rational"))?;
            let r = rational_field(right, &format!("{component}_exact_rational"))?;
            let le = rational_field(&left, &format!("{component}_error_exact_rational"))?;
            let re = rational_field(right, &format!("{component}_error_exact_rational"))?;
            if le < 0 || re < 0 {
                return Err("physical reference has a negative error radius".into());
            }
            for (row, center, error) in [(&left, &l, &le), (right, &r, &re)] {
                if row[format!("{component}_exact_zero")] == true
                    && (!center.is_zero() || !error.is_zero())
                {
                    return Err(
                        "an exact-zero reference component must have zero center and error".into(),
                    );
                }
            }
            if (l - r).abs() > le + re {
                return Err(format!("physical numerical references disagree beyond tracked errors at epsilon^{order}").into());
            }
        }
    }
    Ok(())
}

pub fn compare(left: &CompletedRun, right: &CompletedRun) -> Result<()> {
    if left.report["input_sha256"] != right.report["input_sha256"] {
        return Err("physical runs have different scalar input identities".into());
    }
    for field in [
        "loop_count",
        "projection",
        "normalization",
        "data",
        "evaluated_input_sha256",
        "dimension_binding",
    ] {
        if left.report[field] != right.report[field] {
            return Err(format!("physical comparisons use different {field}").into());
        }
    }
    if !(canonical_retained(&left.value)? - canonical_retained(&right.value)?)
        .together()
        .expand()
        .is_zero()
    {
        return Err("physical reduction strategies differ exactly".into());
    }
    let left = &left.report["numerical_reference"];
    let right = &right.report["numerical_reference"];
    if !left.is_null() || !right.is_null() {
        compare_numerical_references(left, right)?;
    }
    Ok(())
}

/// Authenticate and compare separate completed processes without repeating IBPs.
pub fn compare_reports(
    dfs: &Path,
    dfs_no_memo: &Path,
    ordered: &Path,
    destination: &Path,
) -> Result<()> {
    let load = |directory: &Path, strategy: &str, memoization: bool| -> Result<CompletedRun> {
        let report: serde_json::Value =
            serde_json::from_slice(&std::fs::read(directory.join("report.json"))?)?;
        if report["schema"] != "vakint.ordered-rung-comparison.v1"
            || report["persistent_cache"] != false
            || report["reduction_completed"] != true
            || report["catalogue_closure"]["closed"] != true
            || report["catalogue_closure"]["unknown_count"] != 0
            || report["strategy"] != strategy
            || report["intrinsic_memoization"] != memoization
            || report["cache_policy"]["depth_first_intrinsic_memoization"] != memoization
            || report["cache_policy"]["prepared_cache_loaded"] != false
            || report["cache_policy"]["previous_run_cache_restored"] != false
        {
            return Err(
                format!("comparison requires a completed {strategy} report with the declared memoization policy and no loaded prepared or previous-run caches").into(),
            );
        }
        if report["retained"]["file"] != "retained.symbolica" {
            return Err("unexpected retained result filename".into());
        }
        let payload = std::fs::read(directory.join("retained.symbolica"))?;
        if report["retained"]["bytes"].as_u64() != Some(payload.len() as u64)
            || report["retained"]["sha256"].as_str() != Some(sha256(&payload).as_str())
        {
            return Err("retained physical result failed byte-length or SHA256 validation".into());
        }
        Vakint::initialize_vakint_symbols();
        let value = Atom::parse(&String::from_utf8(payload)?, "vakint", Default::default())?;
        Ok(CompletedRun { value, report })
    };
    let memoized = load(dfs, "DepthFirst", true)?;
    let nonmemoized = load(dfs_no_memo, "DepthFirstWithoutMemoization", false)?;
    let ordered_result = load(ordered, "Ordered", false)?;
    compare(&memoized, &nonmemoized)?;
    compare(&memoized, &ordered_result)?;
    let dfs_report = std::fs::canonicalize(dfs.join("report.json"))?;
    let dfs_no_memo_report = std::fs::canonicalize(dfs_no_memo.join("report.json"))?;
    let ordered_report = std::fs::canonicalize(ordered.join("report.json"))?;
    std::fs::write(
        destination,
        serde_json::to_vec_pretty(&serde_json::json!({
            "schema":"vakint.ordered-rung-parity.v2","exact_retained_equality":true,
            "numerical_enclosures_agree":if memoized.report["numerical_reference"].is_null() { None } else { Some(true) },
            "input_sha256":memoized.report["input_sha256"],
            "reports":{
                "dfs":{"path":dfs_report,"sha256":sha256(&std::fs::read(&dfs_report)?)},
                "dfs_no_memo":{"path":dfs_no_memo_report,"sha256":sha256(&std::fs::read(&dfs_no_memo_report)?)},
                "ordered":{"path":ordered_report,"sha256":sha256(&std::fs::read(&ordered_report)?)}},
            "retained":{"dfs":memoized.report["retained"]["sha256"],"dfs_no_memo":nonmemoized.report["retained"]["sha256"],"ordered":ordered_result.report["retained"]["sha256"]},
            "intrinsic_memoization":{"dfs":true,"dfs_no_memo":false,"ordered":false},
            "equality_scope":"exact retained equality after documented positive-real EulerGamma normalization; exact binary-rational tracked-error overlap when numerical references are present (rigorous bounds only when error quality is Bounded); no claim of identical exceptional-dimensional guard sets",
        "comparison_normalization":"(exp(±EulerGamma))^p = exp(±EulerGamma*p), using only the registered real Euler constant and its positive exponentials",
            "numerical_check":{"minimum_meaningful_digits":64,"point":{"mUV":1,"mu_r_squared":1,"UFO_G":1}}
        }))?,
    )?;
    Ok(())
}

pub fn prepared(
    directory: &Path,
    projection: Projection,
    loop_count: usize,
) -> Result<LadderFixture> {
    if LadderFixture::has_prepared_input(directory)? {
        LadderFixture::load_prepared_at_loop_count(directory, projection, loop_count)
    } else {
        let fixture = super::prepare_at_loop_count(directory, projection, loop_count)?;
        fixture.save_prepared(directory)?;
        Ok(fixture)
    }
}
