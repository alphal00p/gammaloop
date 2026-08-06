# IR-safe threshold fixtures

These are directive-free saved-DOT fixtures for the two `e+ e- > t t~ H`
NNLO graph topologies used to test IR-safe threshold subtraction.  Add
test-specific `threshold_counterterms` directives after loading or in dedicated
derived fixtures; keep these files as the normalized topology baseline.

A bare saved DOT does not retain the process-valid physical-cut selection.
Every import helper must set `global.generation.force_cuts` with `e<ID>` edge
names before `generate existing`; otherwise GammaLoop reconstructs all generic
s-t cuts and the original threshold associations are lost.  The full forced-cut
lists and the smaller target-cut lists are in the experiment README linked
below.

The files were revalidated on 2026-08-06 against graph-only generation with
the current 166-graph process definition.  The selection explicitly excluded
`ANY_RAISING` propagator and cut signatures, so every represented cut and
threshold group has multiplicity one.

| fixture | SHA-256 | generation LMB |
| --- | --- | --- |
| `GL297.dot` | `818cfcdc61868fd172d8b5bd52efe7f4801a7cc5eac4ded41b48c138e59c80b2` | `(4,5,6,11)` |
| `GL638.dot` | `c6736f0c902066f3b11d3bb7ee9a0c7a4b178bef8967518c138e3f2dbbebe0ff` | `(3,4,7,10)` |

Use exactly one generation-time orientation in the graph-specific tests, and
assert the signature itself rather than its generated ordinal:

- GL297: `(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)`
- GL638: `(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)`

For GL297 the selected orientation contains physical cut `(2,11,14)`, left
thresholds `(5,11,13)` and `(5,7)`, and right thresholds `(3,11,12)` and
`(3,9)`.  For GL638 it contains physical cut `(2,4,12)`, left thresholds
`(8,12,14)` and `(7,8)`, and right thresholds `(5,12,13)`, `(5,10)`, and the
physical-cut geometry `(2,6,10)`.  The focused Cartesian regression explicitly
disables unrelated left `(8,12,14)` and right `(2,6,10)`, duplicates `(7,8)`
into the verified `[7]` and `[3,7]` variants, and keeps `(5,12,13)` and `(5,10)`
as its two distinct right partners.  Its `[3,7]` variant requires explicit
parent LMB `(3,4,7,10)`: omitting the parent produces genuinely different
common-parent embeddings and is rejected as ambiguous.

At the audit baseline, forcing a cut together with a generation-time
orientation pattern exposed a missing-CFF-surface error for topology-wide
thresholds excluded by that pattern.  The regression now distinguishes the
topology identity from its optional selected-orientation instance: excluded
directives remain matched but dormant, while a selected instance with genuinely
missing CFF data still errors.  Tests must not conceal failures by silently
dropping the forced cut, threshold CTs, or the full-UV requirement of the final
GL297 cure.

The complete provenance, cut/E-surface map, edge map, reproduction commands,
and differences from the historical candidates are recorded in
`examples/cli/epem_a_ttxh/NNLO_ir_safe_thresholds/README.md`.
