//! Uncertified scalar inputs survive exact lowering without invoking a reducer.
use super::*;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedMasterMode, Vakint,
    vakint_parse,
};
use std::sync::Arc;

fn owned_family(source: &rustred::family::IntegralFamily) -> rustred::family::IntegralFamily {
    rustred::family::IntegralFamily::new(
        source.name(),
        source.loop_momenta().to_vec(),
        source.external_momenta().to_vec(),
        source.coefficient_context().clone(),
        source.dimension().clone(),
        source.denominators().to_vec(),
        source.external_gram().to_vec(),
        source.power_shifts().to_vec(),
    )
    .unwrap()
}

/// Produce a genuine source-replayed closed sunset with a degree-zero entry
/// certificate. This tests external bounded programs without forging metadata.
fn bounded_sunset() -> (Vakint, Arc<rustred::family::IntegralFamily>) {
    use rustred::foundry::artifact::SourcePortAudit;
    use rustred::sector::{Mask, zero};
    use rustred::solver::{SectorConfig, SectorSolveOptions, SectorSolver, SourceSystem};
    let mut engine = Vakint::new().unwrap();
    let source = super::super::artifact::ArtifactFamily::UnitMassVacuumK3
        .artifact()
        .unwrap();
    let family = owned_family(source.family());
    let input_family = Arc::new(owned_family(&family));
    let mut zeros = Vec::new();
    let mut nonzero = Vec::new();
    {
        let analyzer = zero::Analyzer::try_unrestricted(&family).unwrap();
        for mask in 0..8 {
            let sector = std::array::from_fn::<_, 3, _>(|axis| mask & (1 << axis) != 0);
            match analyzer.analyze(&Mask::try_new(sector).unwrap()).unwrap() {
                zero::Decision::ProvedZero(_) => zeros.push(sector),
                zero::Decision::Inconclusive(_) => nonzero.push(sector),
                zero::Decision::Excluded(_) => {
                    panic!("unrestricted fixture cannot exclude a sector")
                }
            }
        }
    }
    let zeros: Arc<[[bool; 3]]> = zeros.into();
    let system = SourceSystem::<3>::from_family(&family).unwrap();
    let sectors = nonzero
        .into_iter()
        .map(|sector| {
            let solution = SectorSolver::new(
                &system,
                sector,
                SectorConfig {
                    zero_sectors: zeros.clone(),
                    ..Default::default()
                },
            )
            .unwrap()
            .solve_sector(SectorSolveOptions::default())
            .unwrap();
            (sector, None, solution)
        })
        .collect::<Vec<_>>();
    let artifact = SourcePortAudit::try_new(&family, zeros)
        .unwrap()
        .install_complete_through_total_excess(family, sectors, 0)
        .unwrap();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("bounded-sunset.rrbin");
    std::fs::write(&path, artifact.encode_durable().unwrap()).unwrap();
    engine.load_rustred_artifacts(2, &path).unwrap();
    (engine, input_family)
}

fn settings(mode: RustRedMasterMode) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: mode },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        form_exe_path: "/partial-lower-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    }
}

#[test_log::test]
fn uncertified_lower_root_is_retained_exactly_and_rejected_only_by_numerical_materialization() {
    let engine = Vakint::new().unwrap();
    let data = engine.rustred_data().unwrap();
    let family = super::super::artifact::ArtifactFamily::UnitMassVacuumK3;
    let snapshot = data.artifact(family).unwrap();
    // The real sealed sunset needs one unit of i64 headroom at entry. This
    // request is structurally valid and cheap to lower, but has no certified
    // reduction. No synthetic certificate or giant reduction is constructed.
    let outside = snapshot.artifact().supported_root_power_bounds()[0]
        .upper()
        .checked_add(1)
        .unwrap();
    let powers = [outside, 1, 1];
    assert!(matches!(
        snapshot
            .program()
            .validate_entry(&IntegralKey::try_new(powers).unwrap()),
        Err(rustred::reduction::ReductionError::OutsideCertifiedRootDomain { .. })
    ));
    let input = vakint_parse!(&format!("topo(I2L(1,{outside},1,1))")).unwrap();
    let retained_settings = settings(RustRedMasterMode::Retain);
    let retained = engine
        .evaluate_integral(&retained_settings, input.as_view())
        .unwrap();
    let expected = data
        .masters
        .symbol(snapshot.artifact().family_fingerprint(), &powers);
    let unit_scale = retained
        .replace(Atom::var(vk_symbol!(&retained_settings.mu_r_sq_symbol)).to_pattern())
        .with(Atom::one());
    assert!(
        (unit_scale - expected.clone())
            .expand()
            .together()
            .is_zero()
    );

    let error = engine
        .evaluate_integral(&settings(RustRedMasterMode::Numerical), input.as_view())
        .unwrap_err()
        .to_string();
    assert!(error.contains("undeclared residual master"), "{error}");
    assert!(error.contains(&outside.to_string()), "{error}");
    assert!(
        error.contains("parent graph") && error.contains(".dot"),
        "{error}"
    );
    let error = engine
        .resolve_rustred_numerical(retained, &retained_settings)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("unresolved masters") && error.contains(&outside.to_string()),
        "{error}"
    );
}

#[test]
fn preserving_uncertified_lower_inputs_still_enforces_scalar_work_limits() {
    let engine = Vakint::new().unwrap();
    let mut settings = settings(RustRedMasterMode::Retain);
    settings.rustred_scalar_numerator_limits.max_input_nodes = 0;
    settings.rustred_scalar_numerator_limits =
        settings.rustred_scalar_numerator_limits.with_diagnostics();
    let matched = MatchedScalarFamily {
        family: super::super::artifact::ArtifactFamily::UnitMassVacuumK3,
        powers: vec![i64::MAX, 1, 1],
        loop_count: 2,
        mass_squared: Atom::one(),
        parent_coordinates: None,
    };
    let error = evaluate(
        &engine,
        matched,
        &settings,
        Atom::one().as_view(),
        &RustRedEvaluationOptions {
            master_mode: RustRedMasterMode::Retain,
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("limit"), "{error}");
}

#[test_log::test]
fn a_bounded_lower_program_preserves_outside_total_excess_roots() {
    let (engine, family) = bounded_sunset();
    let key = IntegralKey::try_new([2, 1, 1]).unwrap();
    let data = engine.rustred_data().unwrap();
    let snapshot = data
        .artifact(super::super::artifact::ArtifactFamily::UnitMassVacuumK3)
        .unwrap();
    assert!(matches!(
        snapshot.program().validate_entry(&key),
        Err(rustred::reduction::ReductionError::OutsideCertifiedTotalExcessDomain { .. })
    ));
    let config = settings(RustRedMasterMode::Retain);
    let actual = engine
        .evaluate_integral(
            &config,
            vakint_parse!("topo(I2L(1,2,1,1))").unwrap().as_view(),
        )
        .unwrap();
    let actual = actual
        .replace(Atom::var(vk_symbol!(&config.mu_r_sq_symbol)).to_pattern())
        .with(Atom::one());
    let expected = data.masters.symbol(family.fingerprint(), key.powers());
    assert!((actual - expected).expand().together().is_zero());
}

#[test]
fn transported_closed_factor_endpoints_keep_typed_scope_refusals() {
    use super::super::factor_provider::{FactorProviderError, ShippedFactorProvider};
    use rustred::reduction::terminal_normalization::{
        VacuumFactorProvider, VacuumProductReductionError,
    };
    let (engine, family) = bounded_sunset();
    let data = engine.rustred_data().unwrap();
    let snapshot = data
        .artifact(super::super::artifact::ArtifactFamily::UnitMassVacuumK3)
        .unwrap();
    snapshot
        .program()
        .validate_entry(&IntegralKey::try_new([1, 1, 0]).unwrap())
        .unwrap();
    let mut provider = ShippedFactorProvider::with_data(data);
    let target = IntegralKey::try_new([1, 1, -1]).unwrap();
    let error = provider
        .reduce(family.clone(), &target, Default::default())
        .unwrap_err();
    assert!(
        matches!(
            &error,
            FactorProviderError::Native(VacuumProductReductionError::UnresolvedFactor { .. })
        ),
        "{error}"
    );
    assert!(ShippedFactorProvider::error_is_missing_coverage(&error));
    let unsafe_error = provider
        .reduce(
            family,
            &IntegralKey::try_new([1]).unwrap(),
            Default::default(),
        )
        .unwrap_err();
    assert!(!ShippedFactorProvider::error_is_missing_coverage(
        &unsafe_error
    ));
}

#[test]
fn cross_integral_cancellation_precedes_master_coefficient_availability() {
    let engine = Vakint::new().unwrap();
    let mut config = settings(RustRedMasterMode::Numerical);
    config.epsilon_symbol = "ep".into();
    config.number_of_terms_in_epsilon_expansion = 8;
    // Each term alone needs coefficients above the vendored ep^4 data, but
    // the exact IBP identity makes the combined result zero at every order.
    let unavailable = engine
        .evaluate_integral(&config, vakint_parse!("topo(I1L(1,1))").unwrap().as_view())
        .unwrap_err()
        .to_string();
    assert!(
        unavailable.contains("unavailable epsilon coefficient"),
        "{unavailable}"
    );
    let input = vakint_parse!("topo(I1L(1,2))-(1-ep)*topo(I1L(1,1))").unwrap();
    assert!(
        engine
            .evaluate_integral(&config, input.as_view())
            .unwrap()
            .is_zero()
    );
}

#[test]
fn raw_scalar_normalization_keeps_dots_independently_of_public_notation() {
    let engine = Vakint::new().unwrap();
    let mut config = settings(RustRedMasterMode::Retain);
    config.epsilon_symbol = "ep".into();
    let integral = vakint_parse!("topo(I1L(1,2))").unwrap();
    let numerator = vakint_parse!("dot(p(1),p(2))").unwrap();
    let mut matched = engine
        .match_topologies_for_settings(&config, integral.as_view())
        .unwrap()
        .unwrap();
    matched.apply_replacement_rules().unwrap();
    config.use_dot_product_notation = false;
    let indexed = engine
        .rustred_reduce(&config, numerator.as_view(), &matched)
        .unwrap();
    config.use_dot_product_notation = true;
    let compact = engine
        .rustred_reduce(&config, numerator.as_view(), &matched)
        .unwrap();
    assert_eq!(indexed, compact);
    assert!(indexed.contains_symbol(S.dot));
    assert!(!indexed.contains_symbol(S.dot_dummy_ind));
}

#[test]
fn tadpole_recurrence_with_external_scalar_cancels_before_availability_in_both_notations() {
    let engine = Vakint::new().unwrap();
    let mut config = settings(RustRedMasterMode::Numerical);
    config.epsilon_symbol = "ep".into();
    config.number_of_terms_in_epsilon_expansion = 8;
    let input = vakint_parse!("dot(p(1),p(2))*(topo(I1L(1,2))-(1-ep)*topo(I1L(1,1)))").unwrap();
    for dots in [false, true] {
        config.use_dot_product_notation = dots;
        assert!(engine.evaluate(&config, input.as_view()).unwrap().is_zero());
        let tensor = engine.tensor_reduce(&config, input.as_view()).unwrap();
        assert!(
            engine
                .evaluate_integral(&config, tensor.as_view())
                .unwrap()
                .is_zero()
        );
    }
}

#[test]
fn cross_integral_cancellation_removes_uncomputed_residuals_before_numeric_checks() {
    let (engine, _) = bounded_sunset();
    let input = vakint_parse!("topo(I2L(1,2,1,1))-(dot(k(1),k(1))-1)*topo(I2L(1,3,1,1))").unwrap();
    // Both starting roots exceed the real degree-zero certificate. Exact
    // numerator lowering leaves the same residual on each side, which cancels.
    assert!(
        engine
            .evaluate_integral(&settings(RustRedMasterMode::Numerical), input.as_view())
            .unwrap()
            .is_zero()
    );
}

#[test]
fn missing_optional_four_loop_factor_data_is_typed_coverage_but_bad_input_is_not() {
    use super::super::factor_provider::{FactorProviderError, ShippedFactorProvider};
    use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
    use rustred::reduction::terminal_normalization::{
        VacuumFactorProvider, VacuumProductReductionError,
    };
    let engine = Vakint::new().unwrap();
    let family = Arc::new(
        Compiler::new(Limits::default())
            .unwrap()
            .compile_text(TextProject {
                name: Some("unloaded_four_loop_factor".into()),
                parameters: None,
                loop_momenta: (1..=4).map(|n| format!("k{n}")).collect(),
                external_momenta: vec![],
                dimension: "d".into(),
                propagators: [
                    "k1", "k2", "k3", "k4", "k1-k2", "k1-k3", "k1-k4", "k2-k3", "k2-k4", "k3-k4",
                ]
                .into_iter()
                .enumerate()
                .map(|(n, q)| TextPropagator {
                    id: format!("D{n}"),
                    expression: format!("({q})^2-1"),
                    target_power: 1,
                    power_shift: None,
                })
                .collect(),
                external_gram: vec![],
                numerator: None,
            })
            .unwrap()
            .into_lowered(LoweringLimits::default())
            .unwrap()
            .into_family(),
    );
    let mut provider = ShippedFactorProvider::with_data(engine.rustred_data().unwrap());
    let missing = provider
        .reduce(
            family.clone(),
            &IntegralKey::try_new([1; 10]).unwrap(),
            Default::default(),
        )
        .unwrap_err();
    assert!(matches!(
        missing,
        FactorProviderError::Native(VacuumProductReductionError::UnresolvedFactor { .. })
    ));
    let malformed = provider
        .reduce(
            family,
            &IntegralKey::try_new([1]).unwrap(),
            Default::default(),
        )
        .unwrap_err();
    assert!(!ShippedFactorProvider::error_is_missing_coverage(
        &malformed
    ));
}
