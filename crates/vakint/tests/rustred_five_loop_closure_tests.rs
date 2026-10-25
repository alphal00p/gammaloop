#![cfg(feature = "feynkit-ingress")]
//! Finite exact closure comparisons; no assertion of arbitrary-index coverage.

#[path = "support/external_data.rs"]
mod external_data_fixture;
mod five_loop_closure;

use std::{cell::Cell, collections::BTreeSet, path::PathBuf, time::Instant};

use five_loop_closure::{Pair, Result, cases};
use serde_json::json;
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, MasterSymbol,
    OfficialExternalDataOptions, RustRedEvaluationOptions, RustRedMasterMode,
    RustRedSymbolicEvaluation, Vakint, VakintSettings, vakint_symbol,
};

struct ClosureEngine {
    engine: Vakint,
    settings: VakintSettings,
    version: String,
    manifest_sha256: String,
    load_seconds: f64,
    completed_cases: Cell<usize>,
}

impl ClosureEngine {
    fn cache_statistics(&self) -> Result<serde_json::Value> {
        let saved = self.engine.rustred_application_cache_statistics()?;
        let factors = self.engine.rustred_application_factor_cache_statistics()?;
        Ok(json!({
            "saved_rule_applications": saved.rule_applications(),
            "saved_cache_hits": saved.cache_hits(),
            "saved_cached_integrals": saved.cached_integrals(),
            "saved_cached_coefficient_bytes": saved.cached_coefficient_bytes(),
            "cached_factors": factors.cached_factors,
            "cached_factor_coefficient_bytes": factors.coefficient_bytes,
        }))
    }
}

fn engine() -> Result<ClosureEngine> {
    let started = Instant::now();
    let mut engine = Vakint::new()?;
    let report = engine.load_official_external_data(
        5,
        &OfficialExternalDataOptions {
            version: external_data_fixture::required_version(),
            local_artifacts_path: std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(PathBuf::from),
            ..Default::default()
        },
    )?;
    // Five-loop vacuum master identities contain one family and fifteen indices.
    // Check the active campaign inventory, independently of which checkout was selected.
    let five_loop_masters = engine
        .master_symbols()?
        .into_iter()
        .filter(|master| {
            matches!(master.as_atom().as_view(), AtomView::Fun(function)
                if function.get_symbol() == vakint_symbol!("RustRedMaster")
                    && function.get_nargs() == 16)
        })
        .count();
    assert_eq!(
        five_loop_masters, 640,
        "closure tests require the authenticated 640-master campaign"
    );
    let settings = VakintSettings {
        epsilon_symbol: "eps".into(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::Custom("1".into()),
        form_exe_path: "/five-loop-closure-must-not-use-form".into(),
        verify_numerator_identification: true,
        ..Default::default()
    };
    // Reuse the exact options used by transactional data activation. Changing
    // preparation_workers here would unnecessarily load the application again.
    let manifest_sha256 = format!(
        "{:x}",
        Sha256::digest(std::fs::read(report.path.join("manifest.json"))?)
    );
    Ok(ClosureEngine {
        engine,
        settings,
        version: report.version,
        manifest_sha256,
        load_seconds: started.elapsed().as_secs_f64(),
        completed_cases: Cell::new(0),
    })
}

fn assert_catalogue(
    engine: &Vakint,
    settings: &VakintSettings,
    result: &RustRedSymbolicEvaluation,
) -> Result<usize> {
    log::info!(target: "vakint::closure", "phase=catalogue_visit_begin guards={}", result.nonzero_conditions.len());
    let catalogue = engine
        .master_symbols()?
        .into_iter()
        .collect::<BTreeSet<_>>();
    let mut pending = vec![result.value.as_view()];
    let mut masters = BTreeSet::new();
    while let Some(atom) = pending.pop() {
        match atom {
            AtomView::Fun(function) if function.get_symbol() == vakint_symbol!("RustRedMaster") => {
                let master = MasterSymbol::try_from(atom.to_owned())?;
                assert!(
                    catalogue.contains(&master),
                    "undeclared residual integral {master}"
                );
                masters.insert(master);
            }
            AtomView::Fun(function) => pending.extend(function.iter()),
            AtomView::Add(add) => pending.extend(add.iter()),
            AtomView::Mul(mul) => pending.extend(mul.iter()),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                pending.extend([base, exponent]);
            }
            _ => {}
        }
    }
    log::info!(target: "vakint::closure", "phase=guard_validation_begin masters={}", masters.len());
    let epsilon = vakint_symbol!(&settings.epsilon_symbol);
    for condition in &result.nonzero_conditions {
        // The native exact coefficient domain has already certified these
        // rational conditions as nonzero in Q(d), and Vakint checks their
        // dependence is only on d. The substitution d=4-2ε has inverse
        // ε=(4-d)/2, so every guard remains nonzero in Q(ε). A finite conjunction
        // therefore excludes only exceptional dimensions. Inspect every raw
        // returned guard without expanding their global LCM (140k conditions).
        assert!(
            !condition.is_zero(),
            "a reduction must not require zero to be nonzero"
        );
        assert!(
            condition
                .get_all_symbols(true)
                .iter()
                .all(|symbol| *symbol == epsilon),
            "unexpected dependency in a generic-d closure guard"
        );
    }
    log::info!(target: "vakint::closure", "phase=catalogue_visit_end");
    Ok(masters.len())
}

fn run_case(name: &str, context: &ClosureEngine) -> Result<()> {
    let ClosureEngine {
        engine,
        settings,
        version,
        manifest_sha256,
        load_seconds,
        ..
    } = context;
    log::info!(target: "vakint::closure", "case={name} phase=construct_fixture");
    let preparation = Instant::now();
    let case = cases()
        .into_iter()
        .find(|case| case.name == name)
        .ok_or("unknown closure case")?;
    let pair = Pair::new(case)?;
    let fixture_seconds = preparation.elapsed().as_secs_f64();
    let cache_before = context.cache_statistics()?;
    log::info!(target: "vakint::closure", "case={name} phase=fixture_ready fixture_seconds={fixture_seconds:.6}");
    log::info!(target: "vakint::closure", "case={name} phase=reduce_factored load_seconds={load_seconds:.6}");
    let started = Instant::now();
    let factored = engine.evaluate_rustred_symbolic(settings, pair.factored.as_view())?;
    let factored_seconds = started.elapsed().as_secs_f64();
    let cache_after_factored = context.cache_statistics()?;
    let second = Instant::now();
    log::info!(target: "vakint::closure", "case={name} phase=reduce_pinched factored_seconds={factored_seconds:.6}");
    let pinched = engine.evaluate_rustred_symbolic(settings, pair.pinched.as_view())?;
    let pinched_seconds = second.elapsed().as_secs_f64();
    let cache_after_pinched = context.cache_statistics()?;
    log::info!(target: "vakint::closure", "case={name} phase=compare_results pinched_seconds={pinched_seconds:.6}");
    assert!(
        !factored.value.is_zero(),
        "massive positive scalar pinch must not vanish"
    );
    assert!(
        (&factored.value - &pinched.value)
            .expand()
            .together()
            .is_zero(),
        "{name}: numerator/pinch mismatch"
    );
    let masters = assert_catalogue(engine, settings, &factored)?;
    assert_catalogue(engine, settings, &pinched)?;

    let expanded_started = Instant::now();
    log::info!(target: "vakint::closure", "case={name} phase=reduce_expanded");
    let expanded = engine.evaluate_rustred_symbolic(settings, pair.expanded.as_view())?;
    let expanded_seconds = expanded_started.elapsed().as_secs_f64();
    let cache_after_expanded = context.cache_statistics()?;
    assert_eq!(
        factored.value, expanded.value,
        "{name}: expansion changed exact reduction"
    );
    assert!(
        factored.nonzero_conditions == expanded.nonzero_conditions,
        "{name}: expansion changed guards ({} versus {} conditions)",
        factored.nonzero_conditions.len(),
        expanded.nonzero_conditions.len()
    );
    assert_eq!(factored.families, expanded.families);
    if pair.case.mass != 1 {
        let mut unit_case = pair.case.clone();
        unit_case.mass = 1;
        let unit_pair = Pair::new(unit_case)?;
        let unit = engine.evaluate_rustred_symbolic(settings, unit_pair.pinched.as_view())?;
        assert_mass_scaling(
            &factored.value,
            &unit.value,
            &Atom::num(pair.case.mass * pair.case.mass),
            pair.pinched_powers.iter().sum(),
            settings,
        );
    }
    let elapsed = started.elapsed().as_secs_f64();
    let cache_after = context.cache_statistics()?;
    let report = json!({
        "schema": "vakint.five-loop-closure.v1", "name": name, "data_version": version,
        "data_manifest_sha256": manifest_sha256,
        "parent_index": pair.case.parent, "mass": pair.case.mass, "rerouted": pair.case.reroute,
        "cancellations": pair.case.cancellations, "source_powers": pair.source_powers,
        "pinched_powers": pair.pinched_powers, "source_distinct_quadratics": pair.source_distinct_quadratics,
        "pinched_loop_blocks": pair.pinched_loop_blocks,
        "source_dot": pair.source_dot, "pinched_dot": pair.pinched_dot,
        "native_factored_numerator": pair.source_numerator.to_canonical_string(),
        "factored_input": pair.factored.to_canonical_string(), "expanded_input": pair.expanded.to_canonical_string(),
        "pinched_input": pair.pinched.to_canonical_string(), "exact_result": factored.value.to_canonical_string(),
        "left_guard_count": factored.nonzero_conditions.len(),
        "right_guard_count": pinched.nonzero_conditions.len(),
        "generic_d_guards_validated": true, "identical_route_guard_parity": true,
        "declared_retained_integrals": masters, "engine_prior_cases": context.completed_cases.get(), "fixture_seconds": fixture_seconds, "load_seconds": load_seconds,
        "cache_snapshots": { "before": cache_before, "after_factored": cache_after_factored,
            "after_pinched": cache_after_pinched, "after_expanded": cache_after_expanded, "after": cache_after },
        "factored_seconds": factored_seconds, "pinched_seconds": pinched_seconds,
        "expanded_seconds": expanded_seconds, "comparison_seconds": elapsed,
        "exact_equality": true, "all_five_loop_terminals_in_catalogue": true,
        "independent_unit_mass_scaling_checked": pair.case.mass != 1,
    });
    if let Some(directory) = std::env::var_os("VAKINT_FIVE_LOOP_CLOSURE_REPORTS") {
        let directory = PathBuf::from(directory);
        std::fs::create_dir_all(&directory)?;
        std::fs::write(
            directory.join(format!("{name}.json")),
            serde_json::to_vec_pretty(&report)?,
        )?;
    }
    println!(
        "case={name} exact_closure=passed data={version} comparison_seconds={elapsed:.6} load_seconds={load_seconds:.6}"
    );
    assert!(
        elapsed <= 240.0,
        "{name}: closure comparison exceeded the selected 240-second budget ({elapsed:.3}s)"
    );
    context
        .completed_cases
        .set(context.completed_cases.get() + 1);
    Ok(())
}

fn assert_mass_scaling(
    value: &Atom,
    unit: &Atom,
    mass_squared: &Atom,
    degree: i64,
    settings: &VakintSettings,
) {
    // Five unnormalized loop measures give (m²)^(5d/2-A), d=4-2ε.
    // All constructed numerator factors have already canceled their powers.
    let exponent =
        Atom::num(10 - degree) - Atom::num(5) * Atom::var(vakint_symbol!(&settings.epsilon_symbol));
    let expected = unit * mass_squared.clone().pow(exponent);
    assert!(
        (value - expected).expand().together().is_zero(),
        "incorrect five-loop mass homogeneity"
    );
}

#[test]
fn independent_native_pinch_fixtures_have_checked_mass_and_routing() -> Result<()> {
    let mut irreducible_parents = BTreeSet::new();
    let cases = cases();
    assert_eq!(cases.len(), 12);
    let mut names = BTreeSet::new();
    let mut classes = BTreeSet::new();
    for case in cases {
        assert!(names.insert(case.name), "duplicate closure case name");
        let pair = Pair::new(case)?;
        if pair.pinched_loop_blocks == [5] {
            irreducible_parents.insert(pair.case.parent);
        }
        assert!(!pair.source_numerator.is_one());
        assert!(pair.source_powers.iter().all(|power| *power > 0));
        assert!(pair.pinched_powers.iter().all(|power| *power > 0));
        assert!(pair.pinched_powers.iter().sum::<i64>() <= 10);
        assert_eq!(
            pair.source_powers.iter().sum::<i64>() - pair.pinched_powers.iter().sum::<i64>(),
            pair.case
                .cancellations
                .iter()
                .map(|(_, power)| power)
                .sum::<i64>(),
            "cancellation changed the total propagator degree"
        );
        classes.insert(if pair.case.cancellations.len() == 1 {
            "single_line"
        } else {
            "multiple_lines"
        });
        for &(slot, power) in pair.case.dots {
            if let Some((_, cancelled)) = pair.case.cancellations.iter().find(|(id, _)| *id == slot)
            {
                classes.insert(if *cancelled == power {
                    "complete_dot"
                } else {
                    "partial_dot"
                });
            }
        }
        if pair.case.mass != 1 {
            classes.insert("nonunit_mass");
        }
        if pair.case.reroute {
            classes.insert("rerouted");
        }
        if pair.case.parent.is_none() {
            assert_eq!(pair.source_powers.len(), 12);
            assert_eq!(pair.source_powers.iter().sum::<i64>(), 14);
            assert_eq!(pair.source_distinct_quadratics, 9);
            classes.insert("repeated_quadratics");
        }
        if pair.pinched_loop_blocks != [5] {
            classes.insert("factorized");
        }
    }
    assert_eq!(
        classes,
        BTreeSet::from([
            "complete_dot",
            "factorized",
            "multiple_lines",
            "nonunit_mass",
            "partial_dot",
            "repeated_quadratics",
            "rerouted",
            "single_line",
        ])
    );
    for parent in 0..4 {
        assert!(
            irreducible_parents.contains(&Some(parent)),
            "parent {parent} lacks a genuinely irreducible five-loop control"
        );
    }
    Ok(())
}

#[test_log::test]
fn retained_external_scalar_notation_preserves_full_five_loop_guards() -> Result<()> {
    let context = engine()?;
    let fixtures = five_loop_closure::factorized::fixtures(&context.engine)?;
    let pair = fixtures
        .iter()
        .find(|pair| pair.label == "single_positive_line")
        .ok_or("missing frozen factorized single-line fixture")?;
    let scalar =
        vakint_symbol!("dot").call((vakint_symbol!("p").call(1), vakint_symbol!("p").call(2)));
    let input = scalar * &pair.cancelled;
    let mut compact_settings = context.settings.clone();
    compact_settings.use_dot_product_notation = true;
    let compact = context
        .engine
        .evaluate_rustred_symbolic(&compact_settings, input.as_view())?;
    let mut indexed_settings = compact_settings.clone();
    indexed_settings.use_dot_product_notation = false;
    let indexed = context
        .engine
        .evaluate_rustred_symbolic(&indexed_settings, input.as_view())?;
    assert!(!compact.value.is_zero());
    assert!(compact.value.contains_symbol(vakint_symbol!("dot")));
    assert!(!indexed.value.contains_symbol(vakint_symbol!("dot")));
    let restored = Vakint::convert_to_dot_notation(indexed.value.as_view());
    assert!((restored - compact.value).expand().together().is_zero());
    assert_eq!(indexed.nonzero_conditions, compact.nonzero_conditions);
    assert_eq!(indexed.families, compact.families);
    Ok(())
}

macro_rules! closure_tests {
    ($($name:ident),+ $(,)?) => { $(
        #[test_log::test]
        #[ignore = "standalone cold pilot; ordinary grouped tests cover this same case"]
        fn $name() -> Result<()> { run_case(stringify!($name), &engine()?) }
    )+ };
}

closure_tests!(
    parent_32745_single,
    parent_31740_single,
    parent_30699_multiple,
    parent_30527_multiple,
    parent_32745_partial_dot,
    parent_32745_complete_dot,
    parent_31740_mixed_dots,
    parent_30699_nonunit_routing,
    ladder_left_partial_attachment,
    ladder_right_complete_attachment,
    ladder_attachment_and_rung,
    ladder_repeated_quadratics_nonunit,
);

fn run_factorized_case(name: &str, context: &ClosureEngine) -> Result<()> {
    let ClosureEngine {
        engine,
        settings,
        version,
        manifest_sha256,
        load_seconds,
        ..
    } = context;
    let setup = Instant::now();
    log::info!(target: "vakint::closure", "case={name} phase=construct_fixture load_seconds={load_seconds:.6}");
    let fixtures = five_loop_closure::factorized::fixtures(engine)?;
    let pair = fixtures
        .iter()
        .find(|pair| pair.label == name)
        .ok_or("unknown factorized fixture")?;
    let cache_before = context.cache_statistics()?;
    let started = Instant::now();
    log::info!(target: "vakint::closure", "case={name} phase=reduce_numerator setup_seconds={:.6}", setup.elapsed().as_secs_f64());
    let cancelled = engine.evaluate_rustred_symbolic(settings, pair.cancelled.as_view())?;
    let cache_after_numerator = context.cache_statistics()?;
    log::info!(target: "vakint::closure", "case={name} phase=reduce_pinched numerator_seconds={:.6}", started.elapsed().as_secs_f64());
    let pinched = engine.evaluate_rustred_symbolic(settings, pair.pinched.as_view())?;
    let cache_after_pinched = context.cache_statistics()?;
    assert!(!cancelled.value.is_zero());
    assert!(
        (&cancelled.value - &pinched.value)
            .expand()
            .together()
            .is_zero()
    );
    let masters = assert_catalogue(engine, settings, &cancelled)?;
    assert_catalogue(engine, settings, &pinched)?;
    if pair.mass_squared != Atom::one() {
        let unit_pair = fixtures
            .iter()
            .find(|pair| pair.mass_squared == Atom::one())
            .ok_or("missing independently constructed unit-mass control")?;
        let unit = engine.evaluate_rustred_symbolic(settings, unit_pair.pinched.as_view())?;
        assert_mass_scaling(
            &cancelled.value,
            &unit.value,
            &pair.mass_squared,
            pair.expected_key.powers().iter().sum(),
            settings,
        );
    }
    let elapsed = started.elapsed().as_secs_f64();
    let cache_after = context.cache_statistics()?;
    if let Some(directory) = std::env::var_os("VAKINT_FIVE_LOOP_CLOSURE_REPORTS") {
        let directory = PathBuf::from(directory);
        std::fs::create_dir_all(&directory)?;
        let report = json!({
            "schema": "vakint.five-loop-closure.v1", "name": name,
            "class": "factorized", "data_version": version, "load_seconds": load_seconds,
            "engine_prior_cases": context.completed_cases.get(),
            "cache_snapshots": { "before": cache_before, "after_numerator": cache_after_numerator,
                "after_pinched": cache_after_pinched, "after": cache_after },
            "data_manifest_sha256": manifest_sha256,
            "cancellations": pair.positive_cancellations,
            "mass_squared": pair.mass_squared.to_canonical_string(),
            "expanded_input": pair.cancelled.to_canonical_string(),
            "pinched_input": pair.pinched.to_canonical_string(),
            "pinched_dot": pair.pinched_dot,
            "exact_result": cancelled.value.to_canonical_string(),
            "left_guard_count": cancelled.nonzero_conditions.len(),
            "right_guard_count": pinched.nonzero_conditions.len(),
            "generic_d_guards_validated": true,
            "declared_retained_integrals": masters,
            "all_five_loop_terminals_in_catalogue": true,
            "comparison_seconds": elapsed, "exact_equality": true,
            "independent_unit_mass_scaling_checked": pair.mass_squared != Atom::one(),
        });
        std::fs::write(
            directory.join(format!("{name}.json")),
            serde_json::to_vec_pretty(&report)?,
        )?;
    }
    println!("case={name} exact_closure=passed data={version} comparison_seconds={elapsed:.6}");
    assert!(
        elapsed <= 240.0,
        "{name}: factorized pair exceeded the selected budget"
    );
    context
        .completed_cases
        .set(context.completed_cases.get() + 1);
    Ok(())
}

#[test_log::test]
#[ignore = "standalone cold pilot; ordinary grouped tests cover this same case"]
fn factorized_single_positive_line() -> Result<()> {
    run_factorized_case("single_positive_line", &engine()?)
}

#[test_log::test]
#[ignore = "standalone cold pilot; ordinary grouped tests cover this same case"]
fn factorized_multiple_lines_nonunit_mass() -> Result<()> {
    run_factorized_case("multiple_lines_nonunit_mass", &engine()?)
}

// Loading the authenticated 5L programs dominates these short comparisons.
// Four bounded groups reuse one immutable registry each, while every compared
// input still follows the complete public FeynKit/Vakint path independently.
#[test_log::test]
fn connected_parent_closure_group() -> Result<()> {
    let context = engine()?;
    for name in [
        "parent_32745_single",
        "parent_31740_single",
        "parent_30699_multiple",
        "parent_30527_multiple",
    ] {
        run_case(name, &context)?;
    }
    assert_eq!(context.completed_cases.get(), 4);
    Ok(())
}

#[test_log::test]
fn dotted_and_rerouted_closure_group() -> Result<()> {
    let context = engine()?;
    for name in [
        "parent_32745_partial_dot",
        "parent_32745_complete_dot",
        "parent_31740_mixed_dots",
        "parent_30699_nonunit_routing",
    ] {
        run_case(name, &context)?;
    }
    assert_eq!(context.completed_cases.get(), 4);
    Ok(())
}

#[test_log::test]
fn ladder_attachment_and_rung_closure_group() -> Result<()> {
    let context = engine()?;
    for name in [
        "ladder_left_partial_attachment",
        "ladder_right_complete_attachment",
        "ladder_attachment_and_rung",
        "ladder_repeated_quadratics_nonunit",
    ] {
        run_case(name, &context)?;
    }
    assert_eq!(context.completed_cases.get(), 4);
    Ok(())
}

#[test_log::test]
fn factorized_closure_group() -> Result<()> {
    let context = engine()?;
    for name in ["single_positive_line", "multiple_lines_nonunit_mass"] {
        run_factorized_case(name, &context)?;
    }
    assert_eq!(context.completed_cases.get(), 2);
    Ok(())
}
