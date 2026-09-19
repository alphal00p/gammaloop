//! Public four-loop RustRed/FMFT comparisons, available without experimental features.
//!
//! The numerical oracle gates are ignored by default and require an explicit
//! FORM path. Native evaluation uses the ordinary public RustRed backend and
//! the shared harness supplies an invalid FORM path to that lane.

#[path = "experimental_rustred_4l/acceptance_inputs.rs"]
mod acceptance_inputs;
#[path = "experimental_rustred_4l/input.rs"]
mod input;
#[path = "experimental_rustred_4l/propagator_pinches.rs"]
mod propagator_pinches;
#[path = "experimental_rustred_4l/public_acceptance.rs"]
mod public_acceptance;
#[path = "experimental_rustred_4l/public_timing.rs"]
mod public_timing;
mod test_utils;

use input::retained_parent_descriptor;
use std::collections::BTreeSet;

/// Full existing numerical inventory, not only parent/dotted/pinch examples.
const FOUR_LOOP_ANALYTIC_FMFT_CASES: [&str; 15] = [
    "test_integrate_4l_h",
    "test_integrate_4l_h_squared_mass",
    "test_integrate_4l_h_rank_4",
    "test_integrate_4l_h_rank_4_additional_symbols_numerator",
    "test_integrate_4l_PR9d_from_H",
    "test_integrate_4l_PR9d_from_X",
    "test_integrate_4l_PR9d_from_H_pinch",
    "test_integrate_4l_PR9d_from_FG",
    "test_integrate_4l_PR9d_from_FG_pinch",
    "test_integrate_4l_PR11d",
    "test_integrate_4l_clover",
    "test_integrate_4l_clover_with_non_unit_scales",
    "test_integrate_4l_dotted_clover",
    "test_integrate_4l_clover_with_numerator",
    // Historical name is misleading: this input is a four-tadpole.
    "test_integrate_1l_decorated_indices_fmft",
];

/// Fixture coverage is distinct from arbitrary-index family closure.
#[test]
fn four_loop_rustred_fmft_inventory_tracks_upstream_tests() {
    let source = [
        include_str!("integral_evaluation_analytic_tests.rs"),
        include_str!("integral_evaluation_freeform_tests.rs"),
    ]
    .concat();
    for case in FOUR_LOOP_ANALYTIC_FMFT_CASES {
        assert!(
            source.contains(&format!("fn {case}()")),
            "upstream four-loop FMFT case disappeared: {case}"
        );
    }
    assert_eq!(
        FOUR_LOOP_ANALYTIC_FMFT_CASES
            .iter()
            .copied()
            .collect::<BTreeSet<_>>(),
        acceptance_inputs::cases()
            .iter()
            .map(|case| case.name)
            .collect::<BTreeSet<_>>()
    );
}

#[test]
fn all_literal_inputs_and_numerical_settings_match_existing_acceptance() {
    acceptance_inputs::fixture_checks::all_literal_inputs_and_numerical_settings_match_existing_acceptance();
}
