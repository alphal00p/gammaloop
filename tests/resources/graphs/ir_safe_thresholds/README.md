# IR-safe threshold fixtures

`GL297.dot` and `GL638.dot` are directive-free topology fixtures. Their focused
[cross-section tests](../../../../crates/gammalooprs/src/processes/cross_section/ir_safe_threshold_fixture_tests.rs)
own the cut selection, orientation signatures and test-specific metadata;
signatures, rather than orientation ordinals, identify the intended cases.

`GL297_socp_failures.json` preserves numerical overlap-solver failures.
`GL638_legacy_multipliers.toml` is the expanded-expression oracle for the compact
function-map regression. `triple_dotted_bubble.dot` exercises raised residues.
Keep these independent fixtures when changing the production examples.

See the canonical [threshold metadata guide](../../../../docs/products/gammaloop/content/threshold-subtraction.typ)
for solve groups, full parent bases, multipliers and validation boundaries.
