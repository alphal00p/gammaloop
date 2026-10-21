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
