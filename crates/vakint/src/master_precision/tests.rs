use super::*;
use crate::utils::vakint_macros::vk_parse;

#[test_log::test]
fn finite_source_precision_warns_once_without_changing_values_or_settings() {
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        ..Default::default()
    };
    let mut warnings = MasterPrecisionWarnings::new(&settings);
    let source = vk_parse!("test_low_precision_master").unwrap();
    let value = vk_parse!("1.23456789").unwrap();
    let saved = value.clone();
    assert!(warnings.check(&source, &value));
    assert!(!warnings.check(&source, &value));
    assert_eq!(value, saved);
    assert_eq!(settings.run_time_decimal_precision, 80);
    assert!(!warnings.check(
        &vk_parse!("exact_rational").unwrap(),
        &vk_parse!("1/3").unwrap()
    ));
    assert!(!warnings.check(&vk_parse!("exact_zero").unwrap(), &Atom::Zero));
}

#[test_log::test]
fn dependency_walk_ignores_unused_constants_and_terminates_on_cycles() {
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        ..Default::default()
    };
    let mut warnings = MasterPrecisionWarnings::new(&settings);
    let a = vk_parse!("precision_a").unwrap();
    let b = vk_parse!("precision_b").unwrap();
    let unused = vk_parse!("precision_unused").unwrap();
    let first = vk_parse!("precision_b+1.23456789").unwrap();
    let second = vk_parse!("precision_a+9.87654321").unwrap();
    let unused_value = vk_parse!("0.12345").unwrap();
    let tables = [
        (&a, &first, None),
        (&b, &second, None),
        (&unused, &unused_value, None),
    ];
    warnings.check_substitutions(a.as_view(), tables);
    assert_eq!(warnings.seen.len(), 2);
    assert!(warnings.seen.contains(&a) && warnings.seen.contains(&b));
    assert!(!warnings.seen.contains(&unused));
    warnings.check_substitutions(a.as_view(), tables);
    assert_eq!(warnings.seen.len(), 2);
}

#[test_log::test]
fn fmft_known_source_metadata_is_checked_before_runtime_resize() {
    crate::Vakint::initialize_vakint_symbols();
    let settings = VakintSettings::default();
    let mut warnings = MasterPrecisionWarnings::new(&settings);
    let source = vk_parse!("PR11dep0").unwrap();
    let target = &crate::fmft_numerics::MASTERS_NUMERIC_SUBSTITUTIONS
        .get(&source)
        .unwrap()
        .0;
    assert!(warnings.check(&source, target));
    let settings = VakintSettings {
        run_time_decimal_precision: 25,
        ..Default::default()
    };
    assert!(!MasterPrecisionWarnings::new(&settings).check(&source, target));
}
