//! Finite admitted physical-numerator controls, not an R2 rectangular census.
use std::{collections::BTreeSet, path::PathBuf, time::Instant};

use rustred::{
    family::IntegralKey,
    foundry::artifact::ClosedArtifact,
    reduction::{Reducer, ReductionLimits},
    scalar_numerator::FamilyScalarNumeratorService,
};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedSymbolicEvaluation, Vakint, VakintExpression, VakintSettings,
    diagram_integral::DiagramIntegralOptions, symbols::S, vakint_parse, vakint_symbol,
};

use crate::topologies::five_loop::data as five_loop_data;

#[path = "../../../examples/five_loop_cancellation_fixtures/mod.rs"]
pub(super) mod cancellation_fixtures;
#[path = "r2_fixtures.rs"]
mod fixtures;
use fixtures::{Fixture, Result};

fn lowering_keys(fixture: &Fixture) -> Result<()> {
    fixture.verify_native_family_map()?;
    let family = cancellation_fixtures::unit_core_family()?;
    let scalar = FamilyScalarNumeratorService::try_new(
        &family,
        S.dot,
        (1..=5).map(|axis| S.k.call(axis)).collect(),
        Default::default(),
    )?;
    let transported_input = fixture.transported_input();
    let expression = crate::VakintExpression::try_from(transported_input.as_view())?;
    let physical_slots = fixture
        .base
        .powers()
        .iter()
        .enumerate()
        .filter_map(|(slot, power)| (*power > 0).then_some(slot))
        .collect::<Vec<_>>();
    assert_eq!(
        physical_slots.len(),
        if fixture.name == "connected_irreducible_difference" {
            6
        } else {
            5
        }
    );
    for (edge, &slot) in physical_slots.iter().enumerate() {
        let properties = crate::get_prop_with_id(expression.0[0].integral.as_view(), edge)
            .ok_or("native DOT ingress lost a physical edge")?;
        let momentum = &properties[&vakint_symbol!("q_")];
        assert_eq!(properties[&vakint_symbol!("mUVsq_")], Atom::one());
        assert_eq!(
            crate::get_integer_from_atom(properties[&vakint_symbol!("pow_")].as_view()),
            Some(fixture.base.powers()[slot])
        );
        // Native lowering proves that numerator AND denominators underwent the
        // same exact map to the declared common-family coordinates.
        let denominator = S.dot.call((momentum, momentum)) - Atom::one();
        let mut unit = [0; 15];
        unit[slot] = 1;
        let coordinate =
            scalar.lower_weighted(&denominator, &IntegralKey::try_new(unit)?, &Atom::one())?;
        assert_eq!(coordinate.terms().len(), 1);
        assert!(
            coordinate.terms()[0]
                .integral()
                .powers()
                .iter()
                .all(|power| *power == 0)
        );
        assert_eq!(
            coordinate.terms()[0].coefficient().to_expression(),
            Atom::one()
        );
        assert_eq!(coordinate.terms()[0].scalar_spectator(), &Atom::one());
    }
    let expected = scalar.lower_weighted(&fixture.common_numerator, &fixture.base, &Atom::one())?;
    let transported = scalar.lower_weighted(
        &fixture.transported_numerator()?,
        &fixture.base,
        &Atom::one(),
    )?;
    assert_eq!(expected.terms(), transported.terms());
    assert_eq!(
        expected.nonzero_conditions(),
        transported.nonzero_conditions()
    );
    assert_eq!(
        expected
            .terms()
            .iter()
            .map(|term| term.integral().clone())
            .collect::<BTreeSet<_>>(),
        fixture.expected
    );
    assert!(!expected.is_zero());
    Ok(())
}

#[test]
fn r2_physical_fixture_native_lowering_and_dot_transport() -> Result<()> {
    Vakint::initialize_vakint_symbols();
    for fixture in [
        fixtures::homogeneous_factor()?,
        fixtures::cross_factor()?,
        fixtures::connected()?,
    ] {
        lowering_keys(&fixture)?;
        println!(
            "fixture={} native_dot_bytes={} scalar_degree={} momentum_degree={} expected_keys={}",
            fixture.name,
            fixture.dot.len(),
            fixture.scalar_degree,
            fixture.momentum_degree,
            fixture.expected.len()
        );
    }
    Ok(())
}

pub(super) fn engine() -> Result<(Vakint, VakintSettings)> {
    let package = PathBuf::from(
        std::env::var_os("VAKINT_SAVED_APPLICATION")
            .ok_or("set VAKINT_SAVED_APPLICATION to the externally exported R2/D10 publication")?,
    );
    let mut settings = VakintSettings {
        epsilon_symbol: "eps".into(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/r2-physical-control-must-not-use-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::Custom("1".into()),
        ..Default::default()
    };
    let started = Instant::now();
    let engine = Vakint::new()?.with_rustred_application(&mut settings, &package)?;
    println!(
        "r2_package_load_seconds={:.6}",
        started.elapsed().as_secs_f64()
    );
    Ok((engine, settings))
}

fn admitted(engine: &Vakint, settings: &VakintSettings, fixture: &Fixture) -> Result<()> {
    lowering_keys(fixture)?;
    let application = engine.saved_application(settings)?;
    let prepared = application
        .prepared
        .lock()
        .map_err(|_| "prepared application lock poisoned")?;
    let scalar = FamilyScalarNumeratorService::try_new(
        prepared.family(),
        S.dot,
        (1..=5).map(|axis| S.k.call(axis)).collect(),
        settings.rustred_scalar_numerator_limits,
    )?;
    let lowered = scalar.lower_weighted(&fixture.common_numerator, &fixture.base, &Atom::one())?;
    let keys = lowered
        .terms()
        .iter()
        .map(|term| term.integral().clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(keys, fixture.expected);
    for key in &keys {
        assert!(
            prepared.routed_entry_is_admitted(key)?,
            "unsupported planned endpoint: {:?}",
            key.powers()
        );
    }
    // Physical bases may have A=11/12 with no negative indices; that is NOT
    // their numerator's entry. Only the exactly lowered endpoints are admitted.
    if fixture.name != "connected_irreducible_difference" {
        assert!(!prepared.routed_entry_is_admitted(&fixture.base)?);
    }
    println!(
        "r2_fixture_admitted={} endpoint_count={} scalar_degree={} momentum_degree={} open_tensor_rank=0",
        fixture.name,
        keys.len(),
        fixture.scalar_degree,
        fixture.momentum_degree
    );
    Ok(())
}

fn evaluate(
    engine: &Vakint,
    settings: &VakintSettings,
    fixture: &Fixture,
) -> Result<RustRedSymbolicEvaluation> {
    admitted(engine, settings, fixture)?;
    let started = Instant::now();
    let first = engine.evaluate_rustred_symbolic(settings, fixture.input.as_view())?;
    println!(
        "case={} first_application_seconds={:.6}",
        fixture.name,
        started.elapsed().as_secs_f64()
    );
    let started = Instant::now();
    let warm = engine.evaluate_rustred_symbolic(settings, fixture.input.as_view())?;
    println!(
        "case={} warm_application_seconds={:.6}",
        fixture.name,
        started.elapsed().as_secs_f64()
    );
    assert_eq!(first.value, warm.value);
    assert_eq!(first.nonzero_conditions, warm.nonzero_conditions);
    assert_eq!(first.families, warm.families);
    assert!(!first.value.is_zero());
    assert!(
        first
            .nonzero_conditions
            .iter()
            .all(|condition| !condition.is_zero())
    );
    if !fixture.pieces.is_empty() {
        let mut reference = Atom::Zero;
        for (coefficient, piece) in &fixture.pieces {
            let reduced = engine.evaluate_rustred_symbolic(settings, piece.as_view())?;
            assert!(
                reduced
                    .nonzero_conditions
                    .iter()
                    .all(|condition| !condition.is_zero())
            );
            reference += Atom::num(*coefficient) * reduced.value;
        }
        assert!((&first.value - reference).expand().together().is_zero());
        println!("case={} exact_lowered_piece_parity=passed", fixture.name);
    }
    println!(
        "case={} coefficients={:?} witnesses={:?} factor_cache={:?}",
        fixture.name,
        engine.rustred_application_cache_statistics()?,
        engine.rustred_application_witness_cache_statistics()?,
        engine.rustred_application_factor_cache_statistics()?
    );
    Ok(first)
}

fn one_loop(power: i64, result: &RustRedSymbolicEvaluation) -> Result<Atom> {
    let artifact =
        ClosedArtifact::decode_durable(include_bytes!("../../../data/1L/ibps/rules.rrbin"))?;
    let label = result.families.iter().find_map(|(label,fingerprint)|
        (fingerprint == artifact.family_fingerprint()).then_some(label))
        .ok_or("factorized result has not been normalized to the independently shipped one-loop family")?;
    let (reduced, conditions, _) = Reducer::reduce_artifact_unit_mass_guarded(
        &artifact,
        &IntegralKey::try_new([power])?,
        ReductionLimits::default(),
    )?;
    let dimension = artifact.family().dimension().to_expression();
    let physical = vakint_parse!("4-2*eps")?;
    let mut value = Atom::Zero;
    for (key, coefficient) in reduced.terms() {
        value += coefficient
            .to_expression()
            .replace(dimension.to_pattern())
            .with(physical.to_pattern())
            * vakint_symbol!("RustRedMaster").call((label.clone(), Atom::num(key.powers()[0])));
    }
    // Keep the existing independent-factor guard convention: native provenance
    // can retain extra guards, but none of the factor's own guards may vanish.
    for condition in conditions {
        let condition = condition
            .to_expression()
            .replace(dimension.to_pattern())
            .with(physical.to_pattern());
        assert!(result.nonzero_conditions.contains(&condition));
    }
    Ok(value)
}

#[test]
#[ignore = "external R2/D10 package; symbolic homogeneous physical numerator"]
fn external_r2_homogeneous_factor_numerator_matches_independent_one_loop_products() -> Result<()> {
    let (engine, settings) = engine()?;
    let fixture = fixtures::homogeneous_factor()?;
    let result = evaluate(&engine, &settings, &fixture)?;
    let (j1, j2, j3) = (
        one_loop(1, &result)?,
        one_loop(2, &result)?,
        one_loop(3, &result)?,
    );
    let expected = (&j1 * &j3 - &j2 * &j2) * &j1 * &j2 * &j3;
    assert!((&result.value - expected).expand().together().is_zero());
    println!("homogeneous_factor_independent_parity=passed");
    Ok(())
}

#[test]
#[ignore = "external R2/D10 package; genuine cross-factor scalar degree2"]
fn external_r2_cross_factor_contraction_matches_independent_one_loop_products() -> Result<()> {
    let (engine, settings) = engine()?;
    let fixture = fixtures::cross_factor()?;
    let result = evaluate(&engine, &settings, &fixture)?;
    let (j1, j2, j3) = (
        one_loop(1, &result)?,
        one_loop(2, &result)?,
        one_loop(3, &result)?,
    );
    let d = vakint_parse!("4-2*eps")?;
    let expected = Atom::num(2) * (&j1 + &j2) * (&j2 + &j3) * &j1 * &j3 * &j3 / d;
    assert!((&result.value - expected).expand().together().is_zero());
    // 1/d belongs to this independent tensor derivation's domain; do not claim
    // that different exact derivations have identical conservative guard lists.
    println!("cross_factor_independent_parity=passed reference_condition=4-2*eps!=0");
    Ok(())
}

#[test]
#[ignore = "external R2/D10 package; potentially slow connected irreducible numerator"]
fn external_r2_connected_irreducible_numerator_matches_exact_lowered_pieces() -> Result<()> {
    let (engine, settings) = engine()?;
    let fixture = fixtures::connected()?;
    evaluate(&engine, &settings, &fixture)?;
    Ok(())
}

#[test]
#[ignore = "external R2/D10 package; six existing positive-line native DOT pinch controls"]
fn external_r2_existing_native_cancellation_controls() -> Result<()> {
    let (engine, settings) = engine()?;
    let fixtures = cancellation_fixtures::fixtures(&engine)?;
    cancellation_fixtures::verify_fixture_shapes(&fixtures)?;
    assert_eq!(fixtures.len(), 6);
    for fixture in fixtures {
        let application = engine.saved_application(&settings)?;
        assert!(
            application
                .prepared
                .lock()
                .map_err(|_| "prepared application lock poisoned")?
                .routed_entry_is_admitted(&fixture.expected_key)?
        );
        let started = Instant::now();
        let first = engine.evaluate_rustred_symbolic(&settings, fixture.cancelled.as_view())?;
        let pinched = engine.evaluate_rustred_symbolic(&settings, fixture.pinched.as_view())?;
        assert!(
            (&first.value - &pinched.value)
                .expand()
                .together()
                .is_zero()
        );
        assert_eq!(first.families, pinched.families);
        assert_eq!(first.nonzero_conditions, pinched.nonzero_conditions);
        let warm = engine.evaluate_rustred_symbolic(&settings, fixture.cancelled.as_view())?;
        assert_eq!(first.value, warm.value);
        assert_eq!(first.nonzero_conditions, warm.nonzero_conditions);
        assert_eq!(first.families, warm.families);
        println!(
            "case={} native_dot_cancellation=passed total_seconds={:.6}",
            fixture.label,
            started.elapsed().as_secs_f64()
        );
    }
    Ok(())
}
