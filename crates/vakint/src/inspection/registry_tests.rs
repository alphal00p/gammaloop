use super::*;
use crate::precision::AvailablePrecision;

fn empty_catalog() -> MasterCatalog {
    crate::Vakint::initialize_vakint_symbols();
    let (mut catalog, _) = MasterCatalog::test_catalog(None, 20);
    catalog.masters.clear();
    catalog
}

fn finite_constant(catalog: &mut MasterCatalog, name: &str) -> Atom {
    let symbol = vk_parse!(name).unwrap();
    catalog.constants.insert(
        symbol.clone(),
        FiniteSource {
            value: "1.23456789012345678901234567890123456789".into(),
            significant_digits: 30,
            absolute_error: Some("1e-38".into()),
            error_is_bounded: true,
            identity: symbol.to_canonical_string(),
            master_id: None,
            epsilon_power: None,
            data_version: catalog.version.clone(),
            source_sha256: "test-certified-primitive".into(),
        },
    );
    let definition = definitions::external_fallback(&symbol, &serde_json::Value::Null, &[])
        .unwrap()
        .unwrap_or_else(|| {
            vk_parse!("vakint::Integral(vakint::t^2,vakint::Interval(vakint::t,0,1))").unwrap()
        });
    catalog
        .constant_definitions
        .insert(symbol.clone(), definition);
    symbol
}

#[test]
fn finite_primitives_share_snapshot_leaves_and_keep_formula_capabilities() {
    let mut catalog = empty_catalog();
    let period = finite_constant(&mut catalog, "K5Moment(1)");
    let coefficient = add_master(
        &mut catalog,
        1,
        "moment-master",
        "PR4",
        Some("2*K5Moment(1)"),
        20,
    );
    let prepared = catalog.prepare_evaluation(coefficient).unwrap();
    assert_eq!(
        prepared.available_precision(),
        AvailablePrecision::FixedPrecision(30)
    );
    assert_eq!(prepared.sources.len(), 1);
    let source = prepared.sources.keys().next().unwrap().clone();
    assert!(source.contains_symbol(vk_symbol!("RustRedConstant")));
    assert_eq!(
        catalog
            .prepare_evaluation(period.clone())
            .unwrap()
            .expression,
        source
    );
    let canceled = catalog.prepare_evaluation(period.clone() - period).unwrap();
    assert!(canceled.expression.is_zero());
    assert!(canceled.sources.is_empty());
    assert!(
        catalog
            .constant_symbols()
            .unwrap()
            .iter()
            .any(|entry| entry.as_atom() == &vk_parse!("K5Moment(1)").unwrap())
    );
    assert!(
        super::super::validate_parameter(
            vk_parse!("K5Moment(1)").unwrap().as_view(),
            Some(&catalog)
        )
        .is_err()
    );
}

#[test]
fn estimated_named_periods_share_formula_sources_without_becoming_exact() {
    let mut catalog = empty_catalog();
    let period = finite_constant(&mut catalog, "independent_period");
    catalog.constants.get_mut(&period).unwrap().error_is_bounded = false;
    let coefficient = add_master(
        &mut catalog,
        1,
        "period-dependent-master",
        "period-dependent-basis",
        Some("2*independent_period"),
        20,
    );
    let prepared = catalog.prepare_evaluation(coefficient.clone()).unwrap();
    assert_eq!(
        prepared.available_precision(),
        AvailablePrecision::FixedPrecision(30)
    );
    assert_eq!(prepared.sources.len(), 1);
    assert!(!prepared.sources.values().next().unwrap().error_is_bounded);
    let canceled = catalog
        .prepare_evaluation(coefficient - Atom::num(2) * period.clone())
        .unwrap();
    assert!(canceled.expression.is_zero());
    assert!(canceled.sources.is_empty());
    assert!(canceled.origins.is_empty());

    let catalog = Arc::new(catalog);
    let constant = catalog
        .get_constant(ConstantSymbol::try_from(period.clone()).unwrap(), 40)
        .unwrap();
    assert_eq!(
        constant.definition(),
        Some(&catalog.constant_definitions[&period])
    );
    assert!(constant.definition_unavailable_reason().is_none());
    assert!(constant.normalisation().starts_with("Plain real number"));
    assert_eq!(
        constant.available_precision(),
        AvailablePrecision::FixedPrecision(30)
    );
    assert_eq!(
        constant.precision().real().quality(),
        crate::ErrorQuality::Estimated
    );
    assert!(constant.warnings().iter().any(|warning| {
        warning.sources.iter().any(|source| {
            source.identity.contains("independent_period")
                && source.available_decimal_precision == Some(30)
        })
    }));
}

#[test]
fn finite_primitive_lookup_is_capped_and_rejects_another_snapshot() {
    let mut catalog = empty_catalog();
    let period = finite_constant(&mut catalog, "K5LogMoment(3)");
    let catalog = Arc::new(catalog);
    let record = catalog
        .get_constant(ConstantSymbol::try_from(period.clone()).unwrap(), 40)
        .unwrap();
    assert_eq!(
        record.available_precision(),
        AvailablePrecision::FixedPrecision(30)
    );
    assert!(record.warnings().iter().any(|warning| {
        warning
            .sources
            .iter()
            .any(|source| source.identity.contains("K5LogMoment"))
    }));
    let retained = catalog.prepare_evaluation(period).unwrap().expression;
    let mut other = empty_catalog();
    let _ = finite_constant(&mut other, "K5LogMoment(3)");
    // A genuine loader computes this tag from all authenticated capabilities
    // and payloads. The synthetic replacement changes it explicitly.
    let retained = retained
        .replace(catalog.coefficient_snapshot_tag().to_pattern())
        .with(Atom::var(vk_symbol!("different_constant_snapshot")));
    assert!(
        other
            .prepare_evaluation(retained.clone())
            .unwrap_err()
            .to_string()
            .contains("another data snapshot")
    );
    assert!(
        crate::NumericalEvaluationResult::from_atom(
            retained.as_view(),
            vk_symbol!("ep"),
            &VakintSettings::default()
        )
        .unwrap_err()
        .to_string()
        .contains("RustRedConstant")
    );
    assert!(
        catalog
            .get_constant(
                ConstantSymbol::try_from(vk_parse!("K5Moment(7)").unwrap()).unwrap(),
                16
            )
            .unwrap_err()
            .to_string()
            .contains("unknown RustRed constant")
    );
}

#[test]
fn finite_function_constants_materialize_inside_scalar_prefactors() {
    let mut catalog = empty_catalog();
    let _ = finite_constant(&mut catalog, "K5Moment(1)");
    let settings = VakintSettings {
        run_time_decimal_precision: 20,
        ..Default::default()
    };
    for input in [
        "K5Moment(1)",
        "exp(K5Moment(1))",
        "K5Moment(1)*unbound_spectator",
    ] {
        let expression = vk_parse!(input).unwrap();
        let evaluated = crate::precision::materialize(
            expression.as_view(),
            &settings,
            &HashMap::new(),
            Some(&catalog),
        )
        .unwrap();
        assert!(
            !evaluated.contains_symbol(vk_symbol!("K5Moment")),
            "{evaluated}"
        );
        assert!(
            !evaluated.contains_symbol(vk_symbol!("RustRedConstant")),
            "{evaluated}"
        );
        assert_eq!(
            evaluated.contains_symbol(vk_symbol!("unbound_spectator")),
            input.contains("unbound_spectator")
        );
    }
}

fn add_master(
    catalog: &mut MasterCatalog,
    index: i64,
    name: &str,
    basis: &str,
    exact: Option<&str>,
    digits: u32,
) -> Atom {
    let (other, _) = MasterCatalog::test_catalog(exact, digits);
    let mut master = other.masters.into_values().next().unwrap();
    master.identity = name.into();
    master.loop_count = 4;
    master.metadata = serde_json::json!({"exact_basis_expression": basis});
    let symbol =
        vk_symbol!("RustRedMaster").call((Atom::var(vk_symbol!("inspection_test_family")), index));
    catalog.masters.insert(symbol.clone(), master);
    catalog.coefficient_symbol(&symbol, 0)
}

#[test]
fn mixed_exact_and_finite_aliases_share_the_best_definition_in_either_order() {
    for exact_first in [false, true] {
        let mut catalog = empty_catalog();
        let a = add_master(
            &mut catalog,
            1,
            "A",
            "shared_basis",
            exact_first.then_some("5/2"),
            20,
        );
        let b = add_master(
            &mut catalog,
            2,
            "B",
            "shared_basis",
            (!exact_first).then_some("5/2"),
            50,
        );
        for symbol in [&a, &b] {
            let prepared = catalog.prepare_evaluation(symbol.clone()).unwrap();
            assert_eq!(prepared.expression, vk_parse!("5/2").unwrap());
            assert_eq!(
                prepared.available_precision(),
                AvailablePrecision::ArbitraryPrecision
            );
        }
        let prepared = catalog.prepare_evaluation(a - b).unwrap();
        assert!(prepared.expression.is_zero());
        assert!(prepared.sources.is_empty());
        assert!(prepared.origins.is_empty());
    }
}

#[test]
fn a_finite_exact_formula_does_not_downgrade_a_better_alias_payload() {
    let mut catalog = empty_catalog();
    let _ = add_master(&mut catalog, 1, "primitive", "PR4", None, 10);
    let exact = add_master(
        &mut catalog,
        2,
        "formula",
        "dependent_basis",
        Some("PR4fin+1"),
        20,
    );
    let numeric = add_master(
        &mut catalog,
        3,
        "better-payload",
        "dependent_basis",
        None,
        40,
    );
    for symbol in [&exact, &numeric] {
        let prepared = catalog.prepare_evaluation(symbol.clone()).unwrap();
        assert_eq!(
            prepared.available_precision(),
            AvailablePrecision::FixedPrecision(40)
        );
        assert_eq!(prepared.sources.len(), 1);
        assert_eq!(
            prepared
                .sources
                .values()
                .next()
                .unwrap()
                .master_id
                .as_deref(),
            Some("better-payload")
        );
    }
    assert!(
        catalog
            .prepare_evaluation(exact - numeric)
            .unwrap()
            .expression
            .is_zero()
    );
}

#[test]
fn canceled_finite_and_arbitrary_precision_origins_are_not_reported() {
    for formula in ["PR4fin+1", "symbolica::zeta(3)"] {
        let mut catalog = empty_catalog();
        let _ = add_master(&mut catalog, 1, "primitive", "PR4", None, 20);
        let a = add_master(&mut catalog, 2, "canceled-A", "basis_A", Some(formula), 20);
        let b = add_master(&mut catalog, 3, "canceled-B", "basis_B", Some(formula), 20);
        let c = add_master(&mut catalog, 4, "survivor", "basis_C", None, 10);
        let prepared = catalog.prepare_evaluation(a - b + c).unwrap();
        assert_eq!(prepared.sources.len(), 1);
        assert_eq!(prepared.origins.len(), 1);
        assert_eq!(prepared.origins[0].master_id.as_deref(), Some("survivor"));
    }
}

#[test]
fn legacy_sources_remain_opaque_and_loaded_primitives_supersede_them() {
    let mut catalog = empty_catalog();
    let original = vk_parse!("PR4fin").unwrap();
    let legacy = precision::legacy::stage_table(
        original.as_view(),
        vk_parse!("2.5`50").unwrap().as_view(),
        vk_symbol!("ep"),
        "FMFT",
    )
    .unwrap();
    let preserved = catalog.prepare_evaluation(legacy.clone()).unwrap();
    assert_eq!(preserved.expression, legacy);
    assert_eq!(preserved.sources.len(), 1);
    let modern = add_master(&mut catalog, 1, "modern-primitive", "PR4", None, 20);
    let prepared = catalog.prepare_evaluation(legacy).unwrap();
    assert_eq!(
        prepared.expression,
        catalog.prepare_evaluation(modern).unwrap().expression
    );
    assert_eq!(
        prepared.available_precision(),
        AvailablePrecision::FixedPrecision(20)
    );
    assert_eq!(
        prepared
            .sources
            .values()
            .next()
            .unwrap()
            .master_id
            .as_deref(),
        Some("modern-primitive")
    );
}

fn legacy_lookup_leaf() -> (Atom, FiniteSource) {
    crate::Vakint::initialize_vakint_symbols();
    // A historical leaf names a table position; the value is this fixture's own.
    let staged = precision::legacy::stage_table(
        vk_parse!("miD6").unwrap().as_view(),
        vk_parse!("3.25`30*ep^-1").unwrap().as_view(),
        vk_symbol!("ep"),
        "MATAD",
    )
    .unwrap();
    let mut sources = precision::legacy::collect_sources(staged.as_view()).unwrap();
    assert_eq!(sources.len(), 1);
    sources.pop_first().unwrap()
}

#[test]
fn legacy_leaf_constant_lookup_retains_its_cap_and_source_provenance() {
    let catalog = Arc::new(empty_catalog());
    let (leaf, source) = legacy_lookup_leaf();
    let requested = source.significant_digits + 20;
    let constant = catalog
        .get_constant(ConstantSymbol::try_from(leaf.clone()).unwrap(), requested)
        .unwrap();
    assert_eq!(constant.symbol().as_atom(), &leaf);
    assert_eq!(constant.requested_decimal_precision(), requested);
    assert_eq!(
        constant.available_precision(),
        AvailablePrecision::FixedPrecision(source.significant_digits)
    );
    assert_eq!(constant.numerical_value(), &Float::with_val(128, 3.25));
    assert!(
        constant.numerical_value().prec() <= precision::decimal_bits(source.significant_digits)
    );
    assert!(constant.precision().imaginary().exact_zero());
    assert_eq!(constant.provenance().version(), source.data_version);
    assert_eq!(constant.provenance().source_sha256(), source.source_sha256);
    assert_eq!(
        constant.provenance().master_id(),
        source.master_id.as_deref()
    );
    assert_eq!(constant.provenance().epsilon_power(), Some(-1));
    assert_eq!(
        constant.definition(),
        Some(&definitions::matad_master("miD6", -1).unwrap())
    );
    assert!(constant.normalisation().starts_with("MATAD normalization"));
    assert!(constant.warnings().iter().any(|warning| {
        warning.requested_decimal_precision == requested
            && warning.sources.iter().any(|contributor| {
                contributor.source_sha256 == source.source_sha256
                    && contributor.available_decimal_precision == Some(source.significant_digits)
                    && contributor.epsilon_power == Some(-1)
            })
    }));
    assert!(constant.to_atom().is_ok());
}

#[test]
fn legacy_leaf_constant_lookup_rejects_malformed_or_corrupted_identities() {
    let catalog = Arc::new(empty_catalog());
    let (leaf, _) = legacy_lookup_leaf();
    let AtomView::Fun(function) = leaf.as_view() else {
        panic!("staged source must remain an opaque function")
    };
    let arguments = function
        .iter()
        .map(|part| part.to_owned())
        .collect::<Vec<_>>();
    let mut corrupted = arguments.clone();
    corrupted[3] += Atom::one();
    let corrupted = vk_symbol!("NumericalSource").call(corrupted.as_slice());
    let error = catalog
        .get_constant(ConstantSymbol::try_from(corrupted).unwrap(), 16)
        .unwrap_err();
    assert!(matches!(error, VakintError::NumericalPrecision(_)));
    assert!(error.to_string().contains("does not match its payload"));

    let mut wrong_namespace = arguments.clone();
    wrong_namespace[0] = Atom::var(vk_symbol!("user_space::source"));
    for malformed in [
        vk_parse!("NumericalSource()").unwrap(),
        vk_symbol!("NumericalSource").call(wrong_namespace.as_slice()),
    ] {
        assert!(matches!(
            catalog.get_constant(ConstantSymbol::try_from(malformed).unwrap(), 16),
            Err(VakintError::NumericalPrecision(_))
        ));
    }
    let unrelated = vk_symbol!("user_space::NumericalSource").call(arguments.as_slice());
    assert!(matches!(
        catalog.get_constant(ConstantSymbol::try_from(unrelated).unwrap(), 16),
        Err(VakintError::Inspection(
            InspectionError::UnknownConstant { .. }
        ))
    ));
}

#[test]
fn constant_records_expose_their_defining_integrals() {
    let mut catalog = empty_catalog();
    let moment = finite_constant(&mut catalog, "K5Moment(1)");
    let undefined = finite_constant(&mut catalog, "undefined_period");
    catalog.constant_definitions.remove(&undefined);
    let catalog = Arc::new(catalog);
    let record = catalog
        .get_constant(ConstantSymbol::try_from(moment).unwrap(), 20)
        .unwrap();
    assert_eq!(
        record.definition(),
        Some(&vk_parse!("vakint::Integral(vakint::t*symbolica::bessel_k(0,vakint::t)^5,vakint::Interval(vakint::t,0,∞))").unwrap())
    );
    // A record without a definition still returns its value, with the reason.
    let undefined = catalog
        .get_constant(ConstantSymbol::try_from(undefined).unwrap(), 20)
        .unwrap();
    assert!(undefined.definition().is_none());
    assert!(
        undefined
            .definition_unavailable_reason()
            .unwrap()
            .contains("has no defining integral")
    );
    assert!(
        undefined
            .normalisation()
            .starts_with("No defining integral")
    );
    assert!(undefined.numerical_value().as_raw().is_finite());
    let zeta = catalog
        .get_constant("z3".parse::<ConstantSymbol>().unwrap(), 20)
        .unwrap();
    assert_eq!(
        zeta.definition(),
        Some(
            &definitions::builtin(vk_parse!("z3").unwrap().as_view())
                .unwrap()
                .unwrap()
        )
    );
    assert!(
        zeta.definition()
            .unwrap()
            .contains_symbol(vk_symbol!("Integral"))
    );
    assert!(zeta.normalisation().starts_with("Plain real number"));
}

#[test]
fn recorded_period_aliases_define_the_bound_master_coefficient() {
    let mut catalog = empty_catalog();
    let _ = add_master(&mut catalog, 1, "primitive", "PR4", None, 20);
    let catalog = Arc::new(catalog);
    let record = catalog
        .get_constant("PR4fin".parse::<ConstantSymbol>().unwrap(), 16)
        .unwrap();
    let definition = record.definition().unwrap();
    assert_eq!(definition, &definitions::fmft_master("PR4", 0).unwrap());
    assert!(record.normalisation().starts_with("Public Vakint/FMFT"));
    let AtomView::Fun(series) = definition.as_view() else {
        panic!("{definition} is not a series coefficient")
    };
    let arguments = series.iter().collect::<Vec<_>>();
    assert_eq!(series.get_symbol(), vk_symbol!("SeriesCoefficient"));
    assert_eq!(arguments[1], Atom::var(vk_symbol!("ep")).as_view());
    assert_eq!(arguments[3], Atom::Zero.as_view());
    let mut loop_integral = None;
    arguments[0].to_owned().visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == vk_symbol!("LoopIntegral")
        {
            loop_integral = Some(
                function
                    .iter()
                    .map(|part| part.to_owned())
                    .collect::<Vec<_>>(),
            );
            return false;
        }
        true
    });
    let loop_integral = loop_integral.unwrap();
    // Dimension, integrand and four loop momenta.
    assert_eq!(loop_integral.len(), 6);
    let AtomView::Mul(denominators) = loop_integral[1].as_view() else {
        panic!("product of propagators expected")
    };
    // PR4 has five unit propagator powers, so the Wick sign is negative.
    assert_eq!(denominators.get_nargs(), 5);
    assert!(matches!(arguments[0], AtomView::Mul(product)
        if product.iter().any(|factor| factor == Atom::num(-1).as_view())));
}

#[test]
fn constants_without_a_definition_keep_their_values_and_reason() {
    // get_constant must not start failing for legacy leaves that have no
    // defining integral: runtime inputs, estimates and incomplete literals.
    use crate::presentation::Present;
    use symbolica::domains::float::Real;
    let catalog = Arc::new(empty_catalog());
    for (backend, original, table, expected) in [
        ("FMFT", "PR3ep2", "717.31`50", "717.31"),
        ("FMFT", "PR13ep1", "27.25`50", "27.25"),
        ("Input", "input_coefficient", "0.5`20", "0.5"),
        ("PySecDec", "pysecdec_result", "0.75`20", "0.75"),
        ("AlphaLoop", "MASTER2m22EP", "0.625`30", "0.625"),
    ] {
        let staged = precision::legacy::stage_table(
            vk_parse!(original).unwrap().as_view(),
            vk_parse!(table).unwrap().as_view(),
            vk_symbol!("ep"),
            backend,
        )
        .unwrap();
        let (leaf, _) = precision::legacy::collect_sources(staged.as_view())
            .unwrap()
            .pop_first()
            .unwrap();
        let constant = catalog
            .get_constant(ConstantSymbol::try_from(leaf).unwrap(), 12)
            .unwrap_or_else(|error| panic!("{backend} {original}: {error}"));
        assert!(constant.definition().is_none(), "{backend} {original}");
        let reason = constant.definition_unavailable_reason().unwrap();
        assert!(reason.contains("no defining integral"), "{reason}");
        assert!(constant.normalisation().contains(reason));
        let difference =
            constant.numerical_value().clone() - Float::parse(expected, Some(128)).unwrap();
        assert!(
            difference.norm() < Float::with_val(128, 1e-10),
            "{backend} {original}"
        );
        let card = format!("{:#}", constant.presentation());
        assert!(card.contains("none registered"), "{card}");
        assert!(card.contains("Normalisation"), "{card}");
    }
    let defined = catalog
        .get_constant("cl4".parse::<ConstantSymbol>().unwrap(), 30)
        .unwrap();
    let card = format!("{:#}", defined.presentation());
    assert!(
        card.contains("Definition") && card.contains("Integral"),
        "{card}"
    );
    assert!(
        card.contains("Normalisation") && card.contains("Clausen value"),
        "{card}"
    );
}

#[test]
fn capped_formulas_are_evaluated_near_their_cap() {
    use symbolica::domains::float::Real;
    // A request far above a finite cap is evaluated at the cap plus a guard;
    // the capped value, precision and warnings equal those of a request at
    // the cap plus one digit.
    let mut catalog = empty_catalog();
    let period = finite_constant(&mut catalog, "K5Moment(1)");
    let coefficient = add_master(
        &mut catalog,
        1,
        "capped-master",
        "PR4",
        Some("2*K5Moment(1)+vakint::D6bar2"),
        20,
    );
    let catalog = Arc::new(catalog);
    let prepared = catalog.prepare_evaluation(coefficient).unwrap();
    assert_eq!(
        prepared.available_precision(),
        AvailablePrecision::FixedPrecision(30)
    );
    let started = std::time::Instant::now();
    let (far, far_warnings) = MasterCatalog::evaluate_prepared(&prepared, 20_000).unwrap();
    assert!(started.elapsed().as_secs() < 120, "{:?}", started.elapsed());
    let (near, _) = MasterCatalog::evaluate_prepared(&prepared, 31).unwrap();
    let difference = far.value().re.clone() - near.value().re.clone();
    assert!(difference.norm() < near.value().re.clone().norm() * Float::with_val(128, 1e-28));
    for number in [&far, &near] {
        assert!(
            number
                .precision()
                .real()
                .significant_digits()
                .is_some_and(|digits| (29..=30).contains(&digits))
        );
    }
    assert!(!far_warnings.is_empty());
    assert!(
        far_warnings
            .iter()
            .all(|warning| warning.requested_decimal_precision == 20_000)
    );
    assert!(catalog.constants.contains_key(&period));
}
