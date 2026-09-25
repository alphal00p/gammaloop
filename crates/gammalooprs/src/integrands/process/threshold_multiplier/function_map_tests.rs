use super::*;
use crate::initialisation::test_initialise;

#[test]
fn collection_deduplication_includes_function_definitions() {
    test_initialise().unwrap();
    let layout = ThresholdMultiplierLayout::new(vec![], vec![], 0, vec![], vec![]).unwrap();
    let definitions = |value: &str| {
        BTreeMap::from([
            ("shape".into(), "nested(2)".into()),
            ("nested(x)".into(), format!("x+{value}")),
        ])
    };
    let first = layout.parse_expression("shape", &definitions("1")).unwrap();
    let second = layout.parse_expression("shape", &definitions("2")).unwrap();
    let values = ThresholdMultiplierInputValues::new(&layout, F(0.0));
    let mut collection = ThresholdMultiplierEvaluatorCollection::build(
        layout,
        vec![
            (ThresholdCountertermVariantId(0), Some(first.clone())),
            (ThresholdCountertermVariantId(1), Some(second)),
            (ThresholdCountertermVariantId(2), Some(first)),
        ],
        vec![],
        &EvaluatorSettings::default(),
    )
    .unwrap()
    .unwrap();
    assert_eq!(collection.evaluators().len(), 2);
    assert_eq!(
        collection.left_variants()[0].evaluator_id,
        collection.left_variants()[2].evaluator_id
    );
    let mut metadata = EvaluationMetaData::new_empty();
    assert_eq!(
        collection.evaluators_mut()[0]
            .evaluate(&values, &mut metadata)
            .unwrap(),
        F(3.0)
    );
    assert_eq!(
        collection.evaluators_mut()[1]
            .evaluate(&values, &mut metadata)
            .unwrap(),
        F(4.0)
    );
}

#[test]
fn definitions_cannot_shadow_registered_kinematic_parameters() {
    test_initialise().unwrap();
    let layout = ThresholdMultiplierLayout::new(
        vec![Atom::var(symbol!("coefficient"))],
        vec![],
        0,
        vec![],
        vec![],
    )
    .unwrap();
    let error = layout
        .parse_expression(
            "coefficient",
            &BTreeMap::from([("coefficient".into(), "2".into())]),
        )
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("conflicts with multiplier input"),
        "{error}"
    );
}
