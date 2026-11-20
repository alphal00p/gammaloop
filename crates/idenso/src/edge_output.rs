//! `IDENSO_EDGE_OUTPUT=<list>`: switches that make an edge-basis contraction
//! produce less, with the same result. Comma-separated; `all` enables every
//! switch. Unset is `flatten` (the production setting since r5); a set value
//! replaces it, so an empty value or one naming no switch (e.g. `off`) leaves
//! everything as before. The value is read once per process.
//!
//! **Recommended setting: `flatten`** (the default). It is the switch the
//! edge-basis contractions need, and its reach is limited to them.
//!
//! - `flatten`: the caller's rename of momentum combinations to independent
//!   vectors expands a nested linear combination of momenta first, so that
//!   the whole combination becomes one vector (gammalooprs
//!   `simplify_edge_basis` / `contract_keeping_edges`, which read the switch
//!   and pass it to `flat_combination`). Only these two renames change.
//! - `dots`: not limited to edge-basis contractions. It reaches every
//!   gammalooprs `simplify_view_certified` call, which then asks the frontier
//!   for dot notation ([`crate::tensor::AlgebraContraction::Dots`]): the R*
//!   contractions of every remainder, hard Taylor term and UV child, in the
//!   loop basis too. It also reaches every `to_dots` call in the process
//!   (idenso and all its callers), which then returns a terminal polynomial
//!   that its producer certified in dot form as it is, instead of rewriting
//!   the whole expression once more. The results are the same, but the
//!   output form of every one of those contractions is the frontier's; use
//!   it (or `all`) only where that whole reach has been checked.
use std::sync::LazyLock;

static SWITCHES: LazyLock<Vec<String>> = LazyLock::new(|| {
    std::env::var("IDENSO_EDGE_OUTPUT")
        .map(|value| value.split(',').map(|v| v.trim().to_owned()).collect())
        .unwrap_or_else(|_| vec!["flatten".to_owned()])
});

/// Whether `IDENSO_EDGE_OUTPUT` names `switch` (or `all`).
pub fn enabled(switch: &str) -> bool {
    SWITCHES
        .iter()
        .any(|entry| entry == switch || entry == "all")
}
