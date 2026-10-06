use super::*;
use crate::ComplexValue;

fn chain_model() -> Model {
    let mut definition: serde_json::Value =
        serde_json::from_str(&Model::phi4().to_json().unwrap()).unwrap();
    for (name, expression) in [("first", "2*mass"), ("second", "first+lam")] {
        definition["parameters"]
            .as_array_mut()
            .unwrap()
            .push(serde_json::json!({
                "name": name, "nature": "internal", "parameter_type": "complex",
                "expression": expression, "value": [999.0, 0.0]
            }));
    }
    definition["couplings"][0]["expression"] = "first*second".into();
    definition["couplings"][0]["value"] = serde_json::json!([888.0, 0.0]);
    Model::from_json(&definition.to_string()).unwrap()
}

#[test]
fn exact_overrides_resolve_dependencies_without_a_card_or_stale_caches() {
    let model = chain_model();
    let before = model.to_json().unwrap();
    let values = model
        .scalar_bindings(
            None,
            &BTreeMap::from([(symbol!("UFO::mass"), Atom::num(3))]),
        )
        .unwrap();
    assert_eq!(values[&symbol!("UFO::first")], Atom::num(6));
    assert_eq!(values[&symbol!("UFO::second")], Atom::num(7));
    assert_eq!(values[&symbol!("UFO::SCALAR_COUPLING")], Atom::num(42));
    assert_eq!(model.to_json().unwrap(), before);
}

#[test]
fn internal_card_values_and_exact_overrides_precede_analytic_definitions() {
    let model = chain_model();
    let before = model.to_json().unwrap();
    let mut card = ParameterCard::new();
    card.insert("mass".into(), ComplexValue::new(3.0, 0.0));
    card.insert("first".into(), ComplexValue::new(7.0, 0.0));
    let values = model
        .scalar_bindings(
            Some(&card),
            &BTreeMap::from([(symbol!("UFO::mass"), Atom::num(5))]),
        )
        .unwrap();
    assert_eq!(values[&symbol!("UFO::mass")], Atom::num(5));
    assert_eq!(values[&symbol!("UFO::first")], Atom::num(7));
    assert_eq!(values[&symbol!("UFO::SCALAR_COUPLING")], Atom::num(56));
    assert_eq!(model.to_json().unwrap(), before);
    let mut applied = model.clone();
    applied.apply_parameter_card(&card).unwrap();
    assert_eq!(
        applied
            .scalar_bindings(Some(&card), &BTreeMap::new())
            .unwrap(),
        model
            .scalar_bindings(Some(&card), &BTreeMap::new())
            .unwrap()
    );
}

#[test]
fn card_complex_binary64_values_remain_exact_rationals() {
    let model = Model::phi4();
    let mut card = ParameterCard::new();
    card.insert("mass".into(), ComplexValue::new(0.1, 0.2));
    let values = model
        .scalar_bindings(Some(&card), &BTreeMap::new())
        .unwrap();
    assert_eq!(
        values[&symbol!("UFO::mass")],
        Atom::num((3602879701896397_i64, 36028797018963968_i64))
            + Atom::num((3602879701896397_i64, 18014398509481984_i64)) * Atom::i()
    );
    card.insert("mass".into(), ComplexValue::new(f64::INFINITY, 0.0));
    assert!(matches!(
        model.scalar_bindings(Some(&card), &BTreeMap::new()),
        Err(ModelError::NonFiniteScalarBinding { .. })
    ));
}

#[test]
fn literal_aliases_resolve_and_cycles_fail_without_mutating_inputs() {
    let (a, b, outside) = (
        symbol!("literal::a_"),
        symbol!("literal::b_"),
        symbol!("literal::outside"),
    );
    let values = resolve_scalar_bindings(BTreeMap::from([
        (a, Atom::var(b) + Atom::var(outside)),
        (b, Atom::num(2)),
    ]))
    .unwrap();
    assert_eq!(values[&a], Atom::num(2) + Atom::var(outside));
    for definitions in [
        BTreeMap::from([(a, Atom::var(a))]),
        BTreeMap::from([(a, Atom::var(b) + 1), (b, Atom::var(a) + 1)]),
    ] {
        assert!(matches!(
            resolve_scalar_bindings(definitions),
            Err(ModelError::UnresolvedScalarBindings)
        ));
    }
    let model = chain_model();
    let before = model.to_json().unwrap();
    let mut card = ParameterCard::new();
    card.insert("missing".into(), ComplexValue::new(1.0, 0.0));
    assert!(matches!(
        model.scalar_bindings(Some(&card), &BTreeMap::new()),
        Err(ModelError::UnknownCardParameter { .. })
    ));
    assert_eq!(model.to_json().unwrap(), before);
}

// Archived native CLI algorithm, intentionally kept only as a migration oracle.
// The new owner must preserve every exact binding at the previously audited point.
fn legacy_cli_bindings(model: &Model, card: &ParameterCard) -> BTreeMap<Symbol, Atom> {
    let mut values = BTreeMap::new();
    for parameter in model.parameters() {
        let analytic = parameter.expression.as_ref().filter(|_| {
            parameter.nature == ParameterNature::Internal && !card.contains_key(&parameter.name)
        });
        let value = if let Some(expression) = analytic {
            expression.clone()
        } else if let Some(value) = parameter.value {
            Atom::num(Rational::try_from(value.re).unwrap())
                + Atom::num(Rational::try_from(value.im).unwrap()) * Atom::i()
        } else if let Some(expression) = &parameter.expression {
            expression.clone()
        } else {
            continue;
        };
        values.insert(symbol!(&format!("UFO::{}", parameter.name)), value);
    }
    for coupling in model.couplings() {
        values.insert(
            symbol!(&format!("UFO::{}", coupling.name)),
            coupling.expression.clone(),
        );
    }
    for _ in 0..=values.len() {
        let next = values
            .iter()
            .map(|(key, value)| {
                (
                    *key,
                    value.replace_multiple(values.iter().map(|(key, value)| {
                        Replacement::new(
                            Pattern::Literal(Atom::var(*key)),
                            Pattern::Literal(value.clone()),
                        )
                    })),
                )
            })
            .collect();
        if next == values {
            break;
        }
        values = next;
    }
    assert!(values.values().all(|value| {
        values
            .keys()
            .all(|key| !value.contains(Atom::var(*key).as_view()))
    }));
    values
}

#[test]
fn standard_model_inline_point_matches_all_legacy_gghh_bindings() {
    let original = Model::standard_model();
    let before = original.to_json().unwrap();
    let mut card = original.default_parameter_card().unwrap();
    let mut overrides = BTreeMap::new();
    for (name, value) in [("MT", 172.5), ("ymt", 172.5), ("WT", 0.0), ("WH", 0.0)] {
        card.insert(name.into(), ComplexValue::new(value, 0.0));
        overrides.insert(
            symbol!(&format!("UFO::{name}")),
            Atom::num(Rational::try_from(value).unwrap()),
        );
    }
    let mut applied = original.clone();
    applied.apply_parameter_card(&card).unwrap();
    let legacy = legacy_cli_bindings(&applied, &card);
    assert_eq!(original.scalar_bindings(None, &overrides).unwrap(), legacy);
    assert_eq!(
        original
            .scalar_bindings(Some(&card), &BTreeMap::new())
            .unwrap(),
        legacy
    );
    assert_eq!(original.to_json().unwrap(), before);
}
