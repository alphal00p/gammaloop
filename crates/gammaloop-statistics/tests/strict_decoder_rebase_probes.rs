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
