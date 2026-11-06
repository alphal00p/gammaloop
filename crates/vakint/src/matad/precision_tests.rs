use super::*;

#[test]
fn rational_poles_are_supplied_by_the_existing_exact_master_formula() {
    Vakint::initialize_vakint_symbols();
    let matad = MATAD::with_settings(VakintSettings::default());
    let source = vk_parse!("miT111").unwrap();
    let supplied = &DIRECT_SUBSTITUTIONS[&source].0;
    let verified = matad
        .preserve_exact_rational_poles(&source, supplied)
        .unwrap();
    let series = verified
        .series(vk_symbol!("ep"), Atom::Zero.as_view(), Rational::from(0))
        .unwrap();
    assert_eq!(series.coefficient((-2).into()).unwrap(), Atom::num((-3, 2)));
    assert_eq!(series.coefficient((-1).into()).unwrap(), Atom::num((-9, 2)));
    assert!(matches!(
        series.coefficient(0.into()).unwrap().as_view(),
        AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Float(..))
    ));
    let staged = crate::precision::legacy::stage_table(
        source.as_view(),
        verified.as_view(),
        vk_symbol!("ep"),
        "MATAD",
    )
    .unwrap();
    let sources = crate::precision::legacy::collect_sources(staged.as_view()).unwrap();
    assert!(!sources.is_empty());
    assert!(
        sources
            .values()
            .all(|source| source.epsilon_power.unwrap() >= 0)
    );
}

#[test]
fn numeric_values_never_establish_an_exact_pole_definition() {
    Vakint::initialize_vakint_symbols();
    let matad = MATAD::with_settings(VakintSettings::default());
    for (source, supplied) in [
        ("unknown_master", "-1.5/ep^2"),
        ("miT111", "-1.6/ep^2"),
        ("miT111", "-1.5"),
    ] {
        let source = vk_parse!(source).unwrap();
        let supplied = vk_parse!(supplied).unwrap();
        assert_eq!(
            matad
                .preserve_exact_rational_poles(&source, &supplied)
                .unwrap(),
            supplied
        );
    }
}

#[test]
fn laurent_expansion_preserves_open_metric_and_scalar_evanescences() {
    use spenso::{
        network::library::symbolic::ETS,
        structure::representation::{BaseRepName, Minkowski, initialize},
    };

    Vakint::initialize_vakint_symbols();
    initialize();
    let matad = MATAD::with_settings(VakintSettings {
        epsilon_symbol: "eps_metric".into(),
        number_of_terms_in_epsilon_expansion: 2,
        integral_normalization_factor: crate::LoopNormalizationFactor::FMFTandMATAD,
        use_dot_product_notation: true,
        ..Default::default()
    });
    let dimension = vk_parse!("d").unwrap();
    let slot = |index| Minkowski::selfless_symbol().call((&dimension, index));
    let metric = ETS.metric(slot(vk_symbol!("open_mu")), slot(vk_symbol!("open_nu")));
    let input = &metric * vk_parse!("Gam(1,1)*(1+ep)/(ep*(1-ep)*(4-2*ep))").unwrap();
    let input = matad.process_matad_form_output(input).unwrap();
    let result = matad
        .finalize_reduced_masters(
            input,
            1,
            None,
            &Atom::num(1),
            &MATADOptions {
                direct_numerical_substition: false,
                substitute_hpls: false,
                ..Default::default()
            },
        )
        .unwrap();
    // The unit-mass tadpole is 1/ep+1+O(ep); the scalar tensor factor
    // contributes another 3/8 to the finite coefficient.
    let metric = metric
        .replace(dimension.to_pattern())
        .with(vk_parse!("4-2*eps_metric").unwrap());
    let expected = metric * vk_parse!("1/(4*eps_metric)+5/8").unwrap();
    assert!(!result.contains_symbol(Symbol::DERIVATIVE));
    assert!(!result.contains_symbol(vk_symbol!("ep")));
    assert!((result - expected).expand().is_zero());
}
