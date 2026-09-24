= Unrestricted RustRed package refresh, September 2026

All seven equation-to-package producer outputs pass their data audit,
installation, combined rebuild and fresh reference-enabled reproduction. The
unchanged 83-test selection, 15 four-loop references and 16 pinch cases pass on
the final packages. All fixed timing matrices pass their numerical comparisons;
observed performance is mixed, not a strict equal-or-better result. Earlier
writer-only results are retained as `WRITER_ONLY_RESULTS.typ` in the evidence
directory and are not substituted for final-package results.

This refresh updates the RustRed dependency from
`8ad62b964de6f3a508fd165dc6ac2509f25839fe` to
`7b22f5bda441588f2e437b15fa0872461c483708`. It retains Vakint's opt-in
RustRed backend, default evaluation order, numerical conventions and original
intended unrestricted rule applicability. The four-loop packages remain
uncertified candidates; neither this refresh nor the numerical tests
assert arbitrary-index closure.

The finite four-loop walker controls are separate input experiments. Their
starting-domain bounds, extra anchor obligations and successful descendant
walks are not substituted for packaged rules or imported as package authority.
No bounded-certified package, rank cut, positive-power cut or dimension cut is
introduced here. All 74 terminal values and the normalization sidecars remain
unchanged.

== Fresh producer and installed payloads

RustRed's `examples/python/generate_vakint_artifacts.py` steers actual native
IBP generation, not loading or re-encoding cached rules. K1/K3 use the public
closing-family generator; K6 uses generic `family-close` from its explicit
equations. H/FG/BMW/X use generic `family-candidates` with the exact Vakint family
inputs, physical roots, auxiliary nonpositive coordinates, sparse backend,
numerical depth 2, SearchFinite and default ordering. No rank, positive-power,
dimension, elapsed-time or cumulative-work cutoff is added. Finite search
strategy settings and explicit input/output safety limits are not a closure
certificate or an index-domain restriction.

Each certified program is inspected and reduced in separate cold processes.
Four-loop normalization is freshly derived from native exact symmetries, then
independently decoded and replayed. The precomputed master catalogs are supplied
separately: their values are not generated, fitted or weakened. All 900 saved
sectors and 1,155 raw terminal declarations still normalize to the same 74
catalog outputs. K1/K3, all four freshly generated normalization sidecars and
all four catalog files are byte-identical to the previous shipped versions.

The new K6 has 5,640 rather than 5,639 rule cells. Independent cold decoding and
cell comparison identify a split into `n0 = 0` and `n0 <= -1` with the same three
right-hand-side terms, plus redundant/proportional guard representations; no
recurrence coefficient change was found. The four candidate programs also
differ from their historical rule payloads. Strict historical-reference
comparisons correctly fail; they are not relaxed or relabelled as equal.
Fresh reproduction against the matching installed producer output is the
separate reproducibility check.

#table(
  columns: 4,
  [Package], [Current rules], [Current native bytes], [Packaged bytes, historical → current],
  [K1], [Unchanged certified program], [2,790], [2,790 → 2,790],
  [K3], [Unchanged certified program], [36,692], [36,692 → 36,692],
  [K6], [5,640], [3,731,581], [3,725,739 → 3,731,581],
  [H], [21,318], [58,025,776], [4,678,555 → 4,657,299],
  [FG], [9,266], [21,996,303], [1,953,460 → 1,939,050],
  [BMW], [9,018], [50,065,170], [3,624,123 → 3,607,218],
  [X], [19,907], [149,404,677], [10,255,210 → 10,177,721],
)

The four gzip programs total 20,381,288 bytes, versus 20,511,348 historically.
Their 59,509 rules retain unrestricted intended applicability; the rule count
alone is neither stronger closure nor a performance measure. The envelope is
still `RRPBIN` version 1, with Certified/V6 for the lower loops and Candidates
for four loops. There is no new envelope schema or bounded-certified substitute.

The producer recipe is published at RustRed revision
`010466a05c86b25d939724a9689e8baab856e7e3`; Vakint retains the tested native
dependency pin `7b22f5bda441588f2e437b15fa0872461c483708`. From the RustRed
checkout containing the producer example:

First follow the
#link("https://github.com/alphal00p/rustred/blob/010466a05c86b25d939724a9689e8baab856e7e3/README.md#symbolica-30-development-checkout")[RustRed Symbolica checkout instructions],
including its existing `patches/symbolica/heap-pow-wide-radix.patch`. The tested
generator used that patched vendored Symbolica source; the GammaLoop consumer
used the unpatched same base revision and passed the numerical gates below.
These are distinct build environments, not a claim of bit-identical libraries.

```sh
mkdir -p TMP
export TMPDIR="$PWD/TMP"
export TMP="$TMPDIR" TEMP="$TMPDIR"
cargo build --release --locked -p rustred-app --bin rustred \
  --example normalize_candidate_terminals -j6
python examples/python/generate_vakint_artifacts.py \
  --loops 1 2 3 4 --workers 6 --output-directory "$TMPDIR/vakint-new" \
  --catalog-directory /path/to/vakint/data/rustred/four_loop \
  --reference-directory /path/to/vakint/data/rustred --execute
```

Use a new output directory under an existing workspace parent. Omit `--execute`
for a plan only. The optional reference must be the matching producer version;
an older different rule program is expected to fail its strict comparison.
The producer never installs artifacts. Publishing requires the separate cold,
numerical and performance gates described below.

== Earlier writer-only checkpoint (superseded payloads)

The trusted packaged input is inspected with `inspect_program`. K1, K3 and K6
must retain `BinaryProgramKind::Certified`; H, FG, BMW and X must retain
`BinaryProgramKind::Candidates`. Other kinds, including bounded-certified,
are refused. The procedure then:

+ Imports the original state and coefficient sections with
  `DecodedCoefficientTable::import_generated`.
+ Interns every coefficient in its original ordinal order into
  `CoefficientTableBuilder`, requiring the returned `CoefficientId` to equal
  the original ordinal. Unused coefficients are retained too.
+ Calls the current native writer's `finish`, replacing only the state and
  coefficient sections, and calls `encode_program` with the original kind,
  other section bytes and section order.
+ Requires `equivalent_generated_programs` to confirm every exact coefficient
  and variable map, original ordering, kind and structural bytes before
  publishing the output.

That initial checkpoint involved no IBP search, source replay, guard-solving
campaign or certificate promotion. A byte-identical result left the file and modification
time untouched, including the existing gzip wrapper. Thus a validated refresh
can legitimately have no binary-file diff; artificial byte churn is not a
format migration.

All seven writer exports passed exact equivalence. K1/K3/K6 were byte-identical.
At that intermediate stage, the four candidate programs differed only in 23
bytes of each 1,017-byte native symbol-state section. Their coefficient-atom,
family and program sections were byte-identical, as were native lengths and
kind. They retained the current
`RRPBIN` version-1 envelope: this is a current-writer re-export, not a new
envelope-schema migration. Its four gzip files totalled 20,467,096 bytes versus
20,511,348 originally. That transport-only checkpoint and its successful gates
remain in the evidence directory, but its K6/four-loop packages and timing
tables are superseded by actual equation-to-package generation above.

The final three-/four-loop gates use a fresh rebuild embedding the generated
packages. One-/two-loop packages remain byte-identical, so their matched
measurements remain applicable. Frozen baseline executables retain their
original embedded programs; later source-file changes cannot silently switch
their numerical or timing inputs.

== Reproducing the validation

Build the same tests on both dependency revisions with Rust 1.97.1, release
optimization, LTO disabled and six build jobs. The recorded runs used CPUs
56–61, single-threaded inner libraries, a sampled 100 GB RSS hard guard with
95 GB cooperative stop and 20 GB host headroom, no elapsed timeout and no
additional address-space limit. No running walker was modified.

```sh
export TMPDIR=/common/dev/rustred/TMP
export TMP="$TMPDIR" TEMP="$TMPDIR"
export CARGO_PROFILE_RELEASE_LTO=off CARGO_BUILD_JOBS=6
cargo test --release --locked --offline -p vakint --no-run \
  --test rustred_refresh_timing_tests --test rustred_four_loop_tests
```

The lower-loop ignored tests are `one_loop_public_scalar_timings`,
`two_loop_public_scalar_timings` and `three_loop_public_scalar_timings` in
`rustred_refresh_timing_tests`. The four-loop timing test is
`public_timing::representative_public_scalar_timings`. Invoke each frozen test
executable with its exact test name and
`--exact --ignored --nocapture --test-threads=1`.

The oracle uses the vendored FORM 5.0.1 executable through explicit
`VAKINT_REFRESH_ORACLE_FORM_PATH` and
`VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH`; its FLINT shared-library directory must
be on `LD_LIBRARY_PATH`. Native lanes deliberately retain forbidden FORM paths.
The license is inherited, never stored in receipts. The first attempted
baseline gate failed because FORM's FLINT library was unavailable; that
environment failure is retained, and the unchanged selection passed with the
explicit runtime library path.

The unchanged four-loop correctness tests are
`public_acceptance::all_fifteen_numerical_references_match_fmft_and_public_rustred`
and `propagator_pinches::sixteen_numerator_propagator_pinches_match_fmft_and_rustred`.
Unset `VAKINT_4L_CANDIDATE_FAMILY_FILTER`, and require all 15 and 16 unique case
markers, not merely two passing test functions.

== Timing interpretation

Each loop uses six fresh-process baseline/candidate pairs in fixed order
BC, CB, BC, CB, BC, CB. Every pair is retained; no reruns are selected for
favorable results. The scalar timer excludes preparation, numerical conversion
and test-process startup. First-family/parent use includes lazy native loading
and application, not a cold-filesystem guarantee. Later inputs in the same
process are labelled separately.
The host also runs the separately pinned five-loop campaign and may run
coordinated work on other cores; disjoint affinity does not eliminate shared
memory/cache contention. These are observed paired timings, not a strict
speed guarantee.

Lower-loop cases retain 31 repeated public calls and compare each result with
MATAD at 32-digit working precision, through the finite part, using the same
normalization and a symmetric relative threshold of 1e-20. Outputs must be
nonempty and finite. Public caches are intentionally retained. The warm
experimental unit is each fresh process's median, not 31 independent samples.
The unchanged four-loop harness compares its native and FMFT results at the
existing precision and tolerance with five repeated pairs per case.

Lower-loop medians across six fresh processes per version are below, in
milliseconds. B is the old dependency; C is the refreshed dependency. The last
column is the median of the six paired warm C/B ratios, which is not generally
the ratio of the two marginal medians. Each loop passed all 768 numerical
comparisons across its 12 processes.

#table(
  columns: 5,
  [Case], [First-call scope], [First B → C, ms], [Warm B → C, ms], [Paired warm C/B],
  [1L/D4], [Family load + apply], [115.369 → 114.473], [46.938 → 47.381], [1.0072],
  [1L/D6], [Next input], [46.944 → 47.379], [46.933 → 47.242], [1.0032],
  [2L/D1 cubed], [Family load + apply], [125.474 → 122.090], [53.339 → 52.406], [0.9788],
  [2L/numerator], [Next input], [49.200 → 49.620], [49.312 → 50.351], [0.9924],
  [3L/D1 squared], [Family load + apply], [1268.583 → 1200.186], [111.779 → 114.517], [1.0115],
  [3L/pinch6], [Next input], [64.165 → 63.620], [64.563 → 64.663], [1.0001],
)

The small mixed lower-loop differences do not establish a uniform speedup or
formal statistical non-inferiority. The first one-loop warm controls are
slightly slower by their paired medians, with faster and slower pairs present.
The final K6 squared-input warm paired median is 1.15% slower, with five of six
pairs slower; its pinch is mixed and essentially unchanged by this summary.
This is not a strict equal-or-better performance pass. No favorable reruns were
selected.

The final fresh-producer four-loop packages passed the six-pair matrix (648
numerical comparisons). First cubed-parent calls include lazy loading; other
rows are subsequent inputs. Units and paired-ratio convention are the same as
above.

#table(
  columns: 4,
  [Case], [First B → C, ms], [Warm B → C, ms], [Paired warm C/B],
  [H/D1 cubed], [5700.938 → 5716.695], [82.743 → 82.149], [0.9825],
  [H/expanded D7], [20.906 → 21.375], [20.899 → 21.227], [1.0120],
  [FG/D1 cubed], [978.607 → 939.737], [78.794 → 78.654], [0.9964],
  [FG/expanded D7], [664.009 → 662.445], [30.192 → 31.085], [1.0318],
  [BMW/D1 cubed], [2901.887 → 2835.692], [64.664 → 64.863], [1.0009],
  [BMW/expanded D7], [24.743 → 24.780], [22.536 → 22.617], [1.0035],
  [X/D1 cubed], [21575.094 → 21608.305], [113.830 → 112.841], [1.0041],
  [X/expanded D7], [329.032 → 320.105], [48.059 → 48.655], [1.0191],
  [Factorized clover], [15.103 → 15.457], [15.000 → 15.013], [0.9995],
)

First-call paired medians range from 0.9670 to 1.0242; warm paired medians range
from 0.9825 to 1.0318. The expanded FG input is 3.18% slower by paired warm
median, with five of six pairs slower. These mixed small differences are not a
uniform speedup, formal equality guarantee or strict equal-or-better pass.
Both baseline and refreshed candidate pass the unchanged 83-test selection
and all 15 four-loop references plus 16 pinch cases. The fixed timing matrix
has no retries or discarded observations.

== Repository checks

The final package rebuild also passes seven focused loader/catalog tests,
three four-loop fixture tests and formatting of the changed Rust files.
`cargo hakari generate` and `cargo hakari verify` pass without generated
changes. The exact evaluated CI graph generator ran locally with the existing
workspace Cargo cache and workspace-only temporary files; its graph and
regenerated `nix-ci.nix` are byte-identical. Top-level CI remains disabled.
Full CI, uploads and unrelated test suites were not run. This local metadata
regeneration is not a claim that the external Nix CI build ran.

The raw evidence, frozen binaries, writer source and commands are in the
RustRed workspace's `TMP/vakint-generic-refresh.uBYwrQ/`.
