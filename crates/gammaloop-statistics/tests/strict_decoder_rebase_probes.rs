//! Adversarial checks for the represented-sample invariant of checked moment JSON.

use gammaloop_statistics::{ExactF64Moments, MomentOrder};
use serde_json::{Value, json};

fn rational(numerator: u128, denominator: u128) -> Value {
    json!({
        "numerator": numerator.to_string(),
        "denominator": denominator.to_string(),
    })
}

fn empty_wire(order: MomentOrder) -> Value {
    serde_json::to_value(ExactF64Moments::empty(order)).unwrap()
}

fn assert_rejected(payload: Value) {
    let bytes = serde_json::to_vec(&payload).unwrap();
    assert!(
        ExactF64Moments::from_json(&bytes).is_err(),
        "checked decoder accepted impossible represented-sample moments: {payload}",
    );
}

#[test]
fn checked_json_rejects_singleton_that_is_not_a_represented_binary64_sample() {
    // 1 + 2^-53 is halfway between adjacent binary64 values, so it cannot
    // be the exact sum of one represented binary64 sample.
    let numerator = (1_u128 << 53) + 1;
    let denominator = 1_u128 << 53;
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("1");
    payload["sum"] = rational(numerator, denominator);
    payload["sum_abs"] = rational(numerator, denominator);
    payload["sum_squared"] = rational(numerator.pow(2), denominator.pow(2));
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_second_moments_inconsistent_with_two_signed_samples() {
    // With two samples, S1=0 and Sabs=2 force {-1, 1}, hence S2=2.
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("2");
    payload["sum_abs"] = rational(2, 1);
    payload["sum_squared"] = rational(3, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_fourth_moment_when_every_sample_must_be_zero() {
    // S2=0 forces every real sample to be zero, hence S4=0 as well.
    let mut payload = empty_wire(MomentOrder::Fourth);
    payload["sample_count"] = json!("2");
    payload["sum_fourth"] = rational(1, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_a_pair_with_a_rational_root_but_wrong_absolute_sum() {
    // S1=1 and S2=5 force the pair {2, -1}, hence Sabs=3.
    // 25/8 passes the general moment bounds but disagrees with that pair.
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("2");
    payload["sum"] = rational(1, 1);
    payload["sum_abs"] = rational(25, 8);
    payload["sum_squared"] = rational(5, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_a_distinct_pair_containing_a_binary64_midpoint() {
    // These are exact moments of {1+2^-53, 2}; the first rational is not binary64.
    let denominator = 1_u128 << 53;
    let numerator = denominator + 1;
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("2");
    payload["sum"] = rational(3 * denominator + 1, denominator);
    payload["sum_abs"] = payload["sum"].clone();
    payload["sum_squared"] = rational(
        numerator.pow(2) + 4 * denominator.pow(2),
        denominator.pow(2),
    );
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_pair_higher_moments_even_when_hankel_bounds_hold() {
    let mut payload = empty_wire(MomentOrder::Fourth);
    payload["sample_count"] = json!("2");
    payload["sum_abs"] = rational(2, 1);
    payload["sum_squared"] = rational(2, 1);
    payload["sum_fourth"] = rational(3, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_identical_midpoints_in_a_two_sample_state() {
    let numerator = (1_u128 << 53) + 1;
    let denominator = 1_u128 << 53;
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("2");
    payload["sum"] = rational(numerator, denominator / 2);
    payload["sum_abs"] = rational(numerator, denominator / 2);
    payload["sum_squared"] = rational(numerator.pow(2), denominator.pow(2) / 2);
    assert_rejected(payload);
}

#[test]
fn checked_json_enforces_the_fourth_moment_bound_beyond_two_samples() {
    let mut payload = empty_wire(MomentOrder::Fourth);
    payload["sample_count"] = json!("3");
    payload["sum_fourth"] = rational(1, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_inconsistent_zero_variance_higher_moments() {
    // Zero variance requires four samples equal to 1, hence S4=4.
    let mut payload = empty_wire(MomentOrder::Fourth);
    payload["sample_count"] = json!("4");
    for name in ["sum", "sum_abs", "sum_squared", "sum_cubed"] {
        payload[name] = rational(4, 1);
    }
    payload["sum_fourth"] = rational(5, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_inconsistent_constant_magnitude_higher_moments() {
    // Equality Sabs^2=N*S2 fixes every |x|=1, hence S4=4.
    let mut payload = empty_wire(MomentOrder::Fourth);
    payload["sample_count"] = json!("4");
    payload["sum_abs"] = rational(4, 1);
    payload["sum_squared"] = rational(4, 1);
    payload["sum_fourth"] = rational(8, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_second_moments_above_the_signed_absolute_bound() {
    // Positive and negative totals are both 1; their squared sums total at most 2.
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("3");
    payload["sum_abs"] = rational(2, 1);
    payload["sum_squared"] = rational(3, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_fractional_constant_magnitude_signed_counts() {
    // Three samples in {-1, 1} cannot have sum zero.
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("3");
    payload["sum_abs"] = rational(3, 1);
    payload["sum_squared"] = rational(3, 1);
    assert_rejected(payload);
}

#[test]
fn checked_json_rejects_constant_magnitude_binary64_midpoints() {
    let numerator = (1_u128 << 53) + 1;
    let denominator = 1_u128 << 53;
    let mut payload = empty_wire(MomentOrder::Second);
    payload["sample_count"] = json!("4");
    payload["sum"] = rational(numerator, denominator / 4);
    payload["sum_abs"] = rational(numerator, denominator / 4);
    payload["sum_squared"] = rational(numerator.pow(2), denominator.pow(2) / 4);
    assert_rejected(payload);
}

#[test]
fn represented_samples_round_trip_across_binary64_range_and_equality_cases() {
    let tiny = f64::from_bits(1);
    let large = 134_217_728.0_f64;
    let cases: &[&[f64]] = &[
        &[0.0],
        &[-0.0],
        &[tiny],
        &[-tiny],
        &[f64::MAX],
        &[-f64::MAX],
        &[-f64::MAX, f64::MAX],
        &[f64::MAX, f64::MAX],
        &[tiny, 1.0],
        &[large, large.next_up()],
        &[0.5, 1.25],
        &[0.5, -1.25],
        &[1.0; 4],
        &[0.0; 4],
        &[-1.0, -1.0, 1.0, 1.0],
        &[-1.0, -1.0, 1.0],
        &[-1.0, 0.0, 1.0],
        &[-tiny, -tiny, tiny, tiny],
        &[-f64::MAX, -f64::MAX, f64::MAX, f64::MAX],
    ];
    for order in [MomentOrder::Second, MomentOrder::Fourth] {
        for samples in cases {
            let mut moments = ExactF64Moments::empty(order);
            for &sample in *samples {
                moments.push(sample).unwrap();
            }
            let bytes = moments.to_canonical_json().unwrap();
            assert_eq!(
                ExactF64Moments::from_json(&bytes).unwrap(),
                moments,
                "{order:?}: {samples:?}",
            );
        }
    }
}
