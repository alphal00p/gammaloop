//! `IDENSO_EDGE_OUTPUT=<list>`: opt-in switches that make an edge-basis
//! contraction produce less, with the same result. Comma-separated; `all`
//! enables every switch; unset (or a value naming none) leaves everything as
//! before.
//!
//! - `flatten`: the caller's rename of momentum combinations to independent
//!   vectors expands a nested linear combination of momenta first, so that
//!   the whole combination becomes one vector (gammalooprs
//!   `simplify_edge_basis` / `contract_keeping_edges`).
//! - `dots`: the caller asks the frontier for dot notation
//!   ([`crate::tensor::AlgebraContraction::Dots`]) and `to_dots` returns a
//!   terminal polynomial that its producer certified in dot form as it is,
//!   instead of rewriting the whole expression once more.
use std::sync::LazyLock;

static SWITCHES: LazyLock<Vec<String>> = LazyLock::new(|| {
    std::env::var("IDENSO_EDGE_OUTPUT")
        .map(|value| value.split(',').map(|v| v.trim().to_owned()).collect())
        .unwrap_or_default()
});

/// Whether `IDENSO_EDGE_OUTPUT` names `switch` (or `all`).
pub fn enabled(switch: &str) -> bool {
    SWITCHES.iter().any(|entry| entry == switch || entry == "all")
}
