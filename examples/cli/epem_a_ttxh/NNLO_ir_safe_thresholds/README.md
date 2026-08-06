# GL297 and GL638 threshold-subtraction revalidation

This directory records the 2026-08-06 graph-specific audit for the IR-safe
threshold-subtraction implementation.  It is an experiment area only; tests
must import the self-contained directive-free fixtures from
`tests/resources/graphs/ir_safe_thresholds` and must not depend on generated
state in this directory or on Feyngen.

The current compact directive snippets are
`GL297.threshold_counterterms.toml.inc` and
`GL638.threshold_counterterms.toml.inc`.  The `.inc` suffix distinguishes these
embedded graph documents from standalone GammaLoop CLI cards.  They mirror the
compact documents embedded in the clean-checkout graphs at
`../NNLO/graphs/GL297.dot` and `../NNLO/graphs/GL638.dot`.  The GL297 document
overrides only the four associations whose legacy CTs are genuinely two-loop.
The GL638 document covers every current association of `(7,8)` with
independently named `[7]` and `[3,7]` variants and the same illustrative
complementary squared-E-surface ratios as the dressed graph.  Those ratios are
finite smoke-test examples only.  Other thresholds remain implicit legacy
defaults, so GammaLoop—not either TOML document—constructs all left/right
Cartesian products.

For a clean checkout, use the dressed graphs and runnable cards under
`../NNLO/analyses/gl297_correlated_soft_threshold` and
`../NNLO/analyses/gl638_multiplier_components`.  The absolute paths and local
state directories retained later in this document are historical audit
provenance, not prerequisites or current runnable instructions.

No command in this audit used privilege escalation or network access.  The
graph selection explicitly removed `ANY_RAISING` propagator and cut
signatures.  Consequently the mappings below concern singleton, multiplicity-
one raised groups; raised behavior needs separate generic fixtures.

## Sources and identity

| graph | verified saved DOT | SHA-256 |
| --- | --- | --- |
| GL297 | `tests/resources/graphs/ir_safe_thresholds/GL297.dot` | `818cfcdc61868fd172d8b5bd52efe7f4801a7cc5eac4ded41b48c138e59c80b2` |
| GL638 | `tests/resources/graphs/ir_safe_thresholds/GL638.dot` | `c6736f0c902066f3b11d3bb7ee9a0c7a4b178bef8967518c138e3f2dbbebe0ff` |

As historical audit provenance, GL297 was compared byte-for-byte with the
then-local preserved state at
`/Users/vjhirsch/Documents/Work/gammaloop_single_parametric_orientation_mc/examples/cli/epem_a_ttxh/NNLO_experiment/gammaloop_state/processes/cross_sections/epem_a_tth/NNLO/GL297.dot`.
GL638 was regenerated from the then-current 166-graph process and saved in the
local, ignored audit state at
`GL638_graph_state/processes/cross_sections/epem_a_tth/NNLO/GL638.dot`.
That regenerated file also matched the previously inspected state hash.

## Historical exact graph-only regeneration

The exact command used during the audit from the repository root was the block
below.  It deliberately preserves the original absolute executable path and
local output-state name as provenance; it is not expected to run in a clean
checkout.  Current users should instead import the tracked dressed graphs or
run the analysis cards linked above.

```bash
/Users/vjhirsch/Documents/Work/gammaloop_single_parametric_orientation_mc/target/release/gammaloop \
  --clean-state -o -l info -L off -p none \
  -s examples/cli/epem_a_ttxh/NNLO_ir_safe_thresholds/GL638_graph_state \
  run -c 'import model sm-default.json; set global kv global.n_cores.feyngen=10; remove processes; generate xs e+ e- > t t~ h | e+ e- g t t~ h d d~ ghG ghG~ a QCD^2==4 QED^2==6 [{{4}} QCD=2] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --symmetrize-left-right-states true -p epem_a_tth -i NNLO --global-prefactor-num "1𝑖" --only-diagrams; select -p epem_a_tth -i NNLO --with-only-graph-names GL638 --without-raised-propagator-signatures ANY_RAISING --without-raised-cuts-signatures ANY_RAISING; save dot; quit -n'
```

It generated 166 grouped graphs, selected only GL638, and saved the hash above.
Replacing both occurrences of `GL638` with `GL297` reproduces the corresponding
graph-only selection.  A saved DOT can then be imported directly, avoiding
Feyngen in subsequent experiments.

### Direct-import cut contract

The DOT topology is self-contained, but a saved DOT does not encode the
process-valid physical-cut selection.  Importing it with the default empty
`global.generation.force_cuts` reconstructed 19 generic s-t cuts and did not
reproduce the threshold associations below.  A direct-import helper must set
the exact process cuts before `generate existing`.

The complete GL297 cut list, in the original channel order, is:

```toml
[cli_settings.global.generation]
force_cuts = [
  ["e2", "e11", "e14"],
  ["e2", "e9", "e12", "e14"],
  ["e2", "e7", "e13", "e14"],
  ["e2", "e6", "e11", "e13"],
  ["e2", "e6", "e9", "e12", "e13"],
  ["e2", "e6", "e7"],
  ["e2", "e4", "e11", "e12"],
  ["e2", "e4", "e9"],
  ["e2", "e4", "e7", "e12", "e13"],
]
```

The complete GL638 cut list is:

```toml
[cli_settings.global.generation]
force_cuts = [
  ["e2", "e6", "e12", "e13"],
  ["e2", "e6", "e10"],
  ["e2", "e6", "e7", "e13", "e14"],
  ["e2", "e4", "e12"],
  ["e2", "e4", "e10", "e13"],
  ["e2", "e4", "e7", "e14"],
]
```

Preprocessing-only structural probes may isolate the problematic target cut:

```toml
# GL297
force_cuts = [["e2", "e11", "e14"]]

# GL638
force_cuts = [["e2", "e4", "e12"]]
```

These one-cut settings must not be used for numerical LU evaluations, which
require the complete cut lists above for IR cancellation.  The lists are
fixture-loading metadata, not threshold-counterterm directives; the baseline
DOTs deliberately remain directive-free.

## Edge order, directions, and generation LMBs

The arrows below are the saved-DOT source and sink.  `back` and `none` record
an explicit DOT `dir`; an empty entry uses the default direction.  Edges 15 and
16 are the outgoing partners sewn to incoming edges 0 and 1.  Generated CFF
orientation signatures therefore have entries for edge IDs 0 through 14.

### GL297

| ID | endpoints | particle | `dir` | LMB representation |
| ---: | --- | --- | --- | --- |
| 0 | `exte0 -> 9` | `e+` | `back` | `P(0)` |
| 1 | `exte1 -> 9` | `e-` | | `P(1)` |
| 2 | `0 -> 1` | `H` | `none` | `-K(1)+K(2)` |
| 3 | `0 -> 3` | `t` | | `-K(2)+K(0)+K(1)` |
| 4 | `7 -> 0` | `t` | | `K(0)` (`lmb_id=0`) |
| 5 | `2 -> 1` | `t` | | `K(1)` (`lmb_id=1`) |
| 6 | `1 -> 6` | `t` | | `K(2)` (`lmb_id=2`) |
| 7 | `5 -> 2` | `t` | | `-P(0)-P(1)+K(1)` |
| 8 | `2 -> 9` | `a` | `none` | `-P(0)-P(1)` |
| 9 | `3 -> 4` | `t` | | `-K(2)-P(0)-P(1)+K(0)+K(1)` |
| 10 | `3 -> 8` | `a` | `none` | `P(0)+P(1)` |
| 11 | `4 -> 5` | `t` | | `K(3)` (`lmb_id=3`) |
| 12 | `4 -> 7` | `g` | `none` | `-K(2)-K(3)-P(0)-P(1)+K(0)+K(1)` |
| 13 | `5 -> 6` | `g` | `none` | `-K(1)+K(3)+P(0)+P(1)` |
| 14 | `6 -> 7` | `t` | | `-K(1)+K(2)+K(3)+P(0)+P(1)` |
| 15 | `8 -> exte15` | `e+` | `back` | sewn partner of 0 |
| 16 | `8 -> exte16` | `e-` | | sewn partner of 1 |

The generation LMB is `(4,5,6,11)`.  GammaLoop generated 160 full bases.  The
nine cut-channel bases were, in channel order,
`(11,12,13,14)`, `(6,11,12,13)`, `(7,12,13,14)`, `(6,7,12,13)`,
`(9,12,13,14)`, `(6,9,12,13)`, `(4,11,12,13)`, `(4,7,12,13)`, and
`(4,9,12,13)`.

### GL638

| ID | endpoints | particle | `dir` | LMB representation |
| ---: | --- | --- | --- | --- |
| 0 | `exte0 -> 9` | `e+` | `back` | `P(0)` |
| 1 | `exte1 -> 9` | `e-` | | `P(1)` |
| 2 | `0 -> 1` | `H` | `none` | `-K(0)+K(1)` |
| 3 | `0 -> 6` | `t` | | `K(0)` (`lmb_id=0`) |
| 4 | `7 -> 0` | `t` | | `K(1)` (`lmb_id=1`) |
| 5 | `3 -> 1` | `t` | | `-P(0)-P(1)+K(3)` |
| 6 | `1 -> 7` | `t` | | `-K(0)-P(0)-P(1)+K(1)+K(3)` |
| 7 | `2 -> 5` | `t` | | `K(2)` (`lmb_id=2`) |
| 8 | `6 -> 2` | `t` | | `-P(0)-P(1)+K(2)` |
| 9 | `2 -> 9` | `a` | `none` | `-P(0)-P(1)` |
| 10 | `4 -> 3` | `t` | | `K(3)` (`lmb_id=3`) |
| 11 | `3 -> 8` | `a` | `none` | `P(0)+P(1)` |
| 12 | `5 -> 4` | `t` | | `K(0)+P(0)+P(1)` |
| 13 | `4 -> 7` | `g` | `none` | `-K(3)+K(0)+P(0)+P(1)` |
| 14 | `5 -> 6` | `g` | `none` | `-K(0)-P(0)-P(1)+K(2)` |
| 15 | `8 -> exte15` | `e+` | `back` | sewn partner of 0 |
| 16 | `8 -> exte16` | `e-` | | sewn partner of 1 |

The generation LMB is `(3,4,7,10)`.  GammaLoop generated 127 full bases.  The
six cut-channel bases were, in channel order, `(6,12,13,14)`,
`(4,12,13,14)`, `(6,7,13,14)`, `(4,7,13,14)`, `(6,10,13,14)`, and
`(4,10,13,14)`.

## Cuts and threshold associations

An `E#` before the physical cut is that cut's E-surface.  Left and right lists
show `E-surface: edge tuple` in GammaLoop's current topology-inferred order.

### GL297

| cut | physical cut | left thresholds | right thresholds |
| ---: | --- | --- | --- |
| 0 | `E38: (2,11,14)` | `E40: (5,11,13)`; `E23: (5,7)` | `E47: (3,11,12)`; `E29: (3,9)` |
| 1 | `E31: (2,9,12,14)` | `E33: (5,9,12,13)`; `E23: (5,7)` | none |
| 2 | `E57: (2,7,13,14)` | none | `E29: (3,9)`; `E62: (3,7,12,13)` |
| 3 | `E39: (2,6,11,13)` | none | `E47: (3,11,12)`; `E29: (3,9)` |
| 4 | `E32: (2,6,9,12,13)` | none | none |
| 5 | `E22: (2,6,7)` | none | `E47: (3,11,12)`; `E29: (3,9)`; `E62: (3,7,12,13)` |
| 6 | `E48: (2,4,11,12)` | `E40: (5,11,13)`; `E23: (5,7)` | none |
| 7 | `E30: (2,4,9)` | `E40: (5,11,13)`; `E33: (5,9,12,13)`; `E23: (5,7)` | none |
| 8 | `E63: (2,4,7,12,13)` | none | none |

Thus the problematic `(5,7)` threshold is left under cuts 0, 1, 6, and 7;
`(3,9)` is right under cuts 0, 2, 3, and 5.  Cut 0 contains both at once.

### GL638

| cut | physical cut | left thresholds | right thresholds |
| ---: | --- | --- | --- |
| 0 | `E44: (2,6,12,13)` | `E29: (8,12,14)`; `E30: (7,8)` | none |
| 1 | `E22: (2,6,10)` | `E29: (8,12,14)`; `E28: (8,10,13,14)`; `E30: (7,8)`; `E11: (3,12)`; `E10: (3,10,13)`; `E47: (3,7,14)` | none |
| 2 | `E53: (2,6,7,13,14)` | none | none |
| 3 | `E15: (2,4,12)` | `E29: (8,12,14)`; `E30: (7,8)` | `E43: (5,12,13)`; `E17: (5,10)` |
| 4 | `E14: (2,4,10,13)` | `E28: (8,10,13,14)`; `E30: (7,8)` | none |
| 5 | `E48: (2,4,7,14)` | none | `E17: (5,10)`; `E52: (5,7,13,14)` |

This exactly confirms the earlier GL638 association map: `(7,8)` is left
under cuts 0, 1, 3, and 4, while cut 3 has two independent right partners.

For the selected cut-3 regression, the generalized `(3,7)` variant cannot be
constructed in any cut-compatible parent: edge 3 lies outside every such
basis.  It is valid as a graph-global two-loop subspace in the generation LMB
`(3,4,7,10)`.  Enumerating every generated common-parent solution with
`parent_lmb` omitted finds genuinely different active-topology/fixed-complement
embeddings, so inference is intentionally rejected as ambiguous.  The explicit
`parent_lmb = [3,4,7,10]` in the experiment and fixture is therefore required,
not a numerical tie-breaker.  In that parent the right-side legacy subspace is
represented by basis edge 10 and is disjoint from both left `[7]` and `[3,7]`
variants.

The current focused-test settings also expose physical-cut geometry
`(2,6,10)` as a third right-side threshold candidate under cut 3.  It is not
one of the two historical partner thresholds `(5,10)` and `(5,12,13)`, so the
fast Cartesian-product fixture explicitly disables it (and the unrelated left
threshold `(8,12,14)`) to assert exactly two left variants, two right partners,
and four iterated pairs.

## Selected orientations

Tests must select and assert these signatures, not an orientation ordinal.
Entries correspond to edge IDs 0 through 14.

### GL297

```toml
[cli_settings.global.generation.orientation_pattern]
pat = "(+,+,-,+,+,+,+,-,0,-,0,-,-,+,+)"
```

For cut 0, the generated counterterm activation gates require:

- physical cut `(2,11,14)`: `e2- e11- e14+`;
- left `(5,11,13)`: `e5+ e11- e13+`;
- left `(5,7)`: `e5+ e7-`;
- right `(3,11,12)`: `e3+ e11- e12-`;
- right `(3,9)`: `e3+ e9-`.

The selected signature satisfies the complete set simultaneously.  Three
other acyclic signatures differ only in the unconstrained `e4`/`e6` signs;
they are deliberately not part of the fixture contract.

### GL638

```toml
[cli_settings.global.generation.orientation_pattern]
pat = "(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)"
```

For cut 3, topology fixes the simultaneous gates as:

- physical cut `(2,4,12)`: `e2+ e4- e12+`;
- left `(8,12,14)`: `e8- e12+ e14+`;
- left `(7,8)`: `e7+ e8-`;
- right `(5,12,13)`: `e5- e12+ e13-`;
- right `(5,10)`: `e5- e10+`.

The selected signature topologically contains both left thresholds and both
right partners for the target cut.  Each generated Cartesian-product term
therefore identifies its left and right components without a user-visible pair
expansion.  The activation gates above were checked against the original
generated counterterms.  At the audit baseline, the then-released generation
path could not complete this same check with both a forced cut and a selected
orientation.  The implementation developed here resolves that historical
limitation, as documented below.

The historical first full-UV direct-import attempt used this exact command.
Like the graph-only command above, it records an absolute executable and an
ignored local state path and is not current clean-checkout guidance:

```bash
/Users/vjhirsch/Documents/Work/gammaloop_single_parametric_orientation_mc/target/release/gammaloop \
  --clean-state -o -l info -L off -p none \
  -s examples/cli/epem_a_ttxh/NNLO_ir_safe_thresholds/GL638_orientation_verified_state \
  run -c 'import model sm-default.json; import graphs examples/cli/epem_a_ttxh/NNLO_ir_safe_thresholds/GL638_graph_state/processes/cross_sections/epem_a_tth/NNLO/GL638.dot -p epem_a_tth -i NNLO; set global kv global.n_cores.generate=10 global.generation.orientation_pattern.pat="(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)" global.generation.evaluator.compile=false global.generation.evaluator.do_algebra=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true global.generation.evaluator.store_atom=true; generate existing -p epem_a_tth -i NNLO; save state -o; quit -n'
```

It reached integrated-UV generation, then failed in the existing Vakint path
with `could not solve momentum system for Vakint integrand` / `Could not solve
Inconsistent` at `crates/gammalooprs/src/uv/approx/integrated.rs:1062`.  This is
an independent pre-existing direct-import limitation, not evidence against the
orientation signature.

For historical metadata-only verification, the command was repeated with local
output state
`GL638_orientation_metadata_state` and the additional setting

```text
global.generation.uv.generate_integrated=false
```

It completed with exactly one stored orientation whose serialized signature is
`(+,+,+,+,-,-,-,+,-,0,+,0,+,-,+)`.  This switch is acceptable only for the
orientation audit.  The final GL297 cure regression must retain integrated UV
physics.

### Selected-orientation topology handling

At the audit baseline, combining either the full or target-only forced-cut list
with the selected GL638 orientation failed during topology threshold
preprocessing.  The failure persisted even with
`global.generation.threshold_subtraction.enable_thresholds=false`:

```text
Topology-discovered threshold surface [EdgeIndex(8), EdgeIndex(10),
EdgeIndex(13), EdgeIndex(14)] is missing from graph 'GL638' CFF surface cache
```

The missing surface is `E28: (8,10,13,14)`.  Its activation gate conflicts
with the chosen orientation, so it correctly has no selected-orientation
instance and no CFF E-surface ID.

The implementation keeps the topology identity independently of that optional
CFF identity.  This permits an explicit directive to match the topology and
remain dormant when orientation selection excludes it.  Only candidates with
a selected-orientation instance enter association classification, raised-data
lookup, and CT generation.  Conversely, if a candidate matches a selected
orientation, its E-surface must still exist in both the graph cache and that
orientation's CFF expression; a genuinely missing surface remains an error.
The GL638 fixture regressions exercise the dormant `E28` case, the genuinely
missing-cache error, and a dormant explicit `(5,10)` directive under cut
`(2,4,12)`.

## Subspace conclusions

- GL638's historical one-loop `[7]` and two-loop `[3,7]` subspaces are direct
  ordered subsets of the verified generation LMB `(3,4,7,10)`.  The completed
  resolver audit found genuinely different common-parent embeddings when the
  parent is omitted, so `parent_lmb = [3,4,7,10]` is required.
- GL297's `[5]` is a direct subset of generation LMB `(4,5,6,11)`, but `[3]`
  is not.  There are 35 generated full bases containing edge 3 (basis IDs
  69 through 103 in this state), and none is a selected cut-channel basis.
  The cut-compatible search resolves both target-cut requests in parent
  `(3,5,11,14)`: `[5]` is left with active topology `(5,6,7,13)`, while `[3]`
  is right with active topology `(3,4,9,12)`.  This confirms that `[3]` must
  not be interpreted in the generation basis and that no explicit parent is
  needed for the verified side-projected embeddings.

## Reconciliation with the research PDF

[`ttH_defo.pdf`](../../../../IMPLEMENTATION_DOCS/ttH_defo.pdf) is conceptual:
it contains no GL labels, edge IDs, cut IDs, LMBs, or orientation signatures.
Its four pages establish that ordinary one-loop-subspace cuts are already
NLO-like; the exceptional topologies require separating a soft region from a
two-loop deformation threshold; and the hard case uses the partition

```text
I_cut * eta_yellow^2 / (eta_yellow^2 + eta_blue^2)
+ I_cut * eta_blue^2 / (eta_yellow^2 + eta_blue^2).
```

The generated metadata supplies the numerical identification absent from the
PDF.  GL297 exercises forced one-loop solving for two problematic thresholds.
GL638 exercises two variants of `(7,8)` in subspaces `[7]` and `[3,7]`, with
the multiplier carrying the multichannel partition.  The exact final physics
multiplier remains follow-up research; it is not a fixture acceptance value.

## Differences from historical hypotheses

- All recorded GL297 and GL638 edge tuples, cut associations, and GL638
  generation LMB were confirmed against the current saved DOT and generated
  metadata.
- The historical GL297 orientation
  `(+,+,-,+,+,+,-,+,0,-,0,-,-,+,+)` does not activate `(5,7)` under cut 0:
  it has `e7+`.  The verified all-partner signature flips `e6/e7` to `+/−`.
- The historical GL638 orientation
  `(+,+,+,+,-,+,-,+,-,0,+,0,+,-,+)` contains the physical cut and left
  thresholds but has `e5+`, so it does not contain both right partners.  The
  verified all-partner signature changes `e5` to `-`.
- Historical orientation ordinals are intentionally discarded.  Filtering
  changes ordinals, whereas the complete edge signature is stable and is what
  tests must assert.
