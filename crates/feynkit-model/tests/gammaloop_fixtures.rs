use feynkit_model::{Model, ParameterCard};

#[test]
fn normalized_scalar_and_standard_model_fixtures_load() {
    let scalar = Model::from_json(include_str!("fixtures/scalars_2p_3p.json"))
        .expect("the normalized scalar fixture should remain compatible");
    assert_eq!(scalar.name(), "scalars");
    assert!(!scalar.vertex_rules().is_empty());

    let standard_model = Model::from_json(include_str!("fixtures/sm.json"))
        .expect("the normalized Standard Model fixture should remain compatible");
    assert_eq!(standard_model.name(), "sm");
    assert_eq!(standard_model.particles().len(), 43);
    assert_eq!(standard_model.vertex_rules().len(), 153);
}

#[test]
fn normalized_restriction_card_applies_to_the_scalar_fixture() {
    let mut model = Model::from_json(include_str!("fixtures/scalars_2p_3p.json"))
        .expect("the normalized scalar fixture should remain compatible");
    let card = ParameterCard::from_json(include_str!("fixtures/restrict_default.json"))
        .expect("the normalized restriction card should remain compatible");

    model
        .apply_parameter_card(&card)
        .expect("the restriction card should contain only known external parameters");
    assert_eq!(
        model.parameter("mass_scalar_2").unwrap().value.unwrap().re,
        2.0
    );
}
