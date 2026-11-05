//! The frozen Forcer grid fixtures `tests/fixtures/forcer/p{1,2,3}l_grid.json`
//! and `p_anchors.json`. They total about 0.22 MB and these checks take about
//! a second, so they all run in the fast gate (`just test_vakint`).
//!
//! Each file must hold exactly the cases of its definition
//! (`vakint::p_integral::forcer_fixtures::{grid_cases, anchor_cases}`), and
//! every frozen output must parse and finalize. With FORM 5,
//!
//! ```text
//! VAKINT_FORCER_FIXTURES=write VAKINT_FORCER_FIXTURES_JOBS=4 \
//!     cargo nextest run -p vakint --test p_integral_forcer_fixtures
//! ```
//!
//! regenerates the files (the anchor values are checked by the unit tests in
//! `src/p_integral/forcer/tests.rs`).

#[path = "support/p_fixtures.rs"]
mod p_fixtures;

use p_fixtures::*;
use vakint::p_integral::forcer_fixtures::ForcerFixture;
use vakint::{LoopNormalizationFactor, RustRedMasterMode, VakintSettings};

fn retain_settings() -> VakintSettings {
    VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

/// Every case parses and finalizes exactly; scaleless cases have no output.
/// A numerator can still cancel a line into a scaleless integral, which
/// Forcer reduces to an exact `0`.
fn assert_complete(file: &str, fixture: &ForcerFixture) {
    let settings = retain_settings();
    for case in &fixture.cases {
        let values = frozen_values(case);
        assert_eq!(
            values.is_empty(),
            case.outputs.is_empty(),
            "{file}: {}",
            case.name
        );
        let value = finalized(case, &settings, RustRedMasterMode::Retain);
        let zero = values.iter().all(|(_, reduced)| reduced.is_zero());
        assert_eq!(value.is_zero(), zero, "{file}: {}", case.name);
    }
}

#[test]
fn the_anchor_fixture_is_complete() {
    let fixture = anchor_fixture();
    assert_complete("p_anchors.json", &fixture);
}

#[test]
fn the_one_loop_grid_fixture_is_complete() {
    assert_complete("p1l_grid.json", &grid_fixture(1));
}

#[test]
fn the_two_loop_grid_fixture_is_complete() {
    assert_complete("p2l_grid.json", &grid_fixture(2));
}

#[test]
fn the_three_loop_grid_fixture_is_complete() {
    assert_complete("p3l_grid.json", &grid_fixture(3));
}
