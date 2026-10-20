= External-data migration validation

This records completed migration checks as of 20 October 2026. The focused
suites overlap; their counts are not a combined whole-repository test result.
Completed scientific checks and repository validation/delivery records appear
separately below.

#let evidence(path, body) = link(
  "https://github.com/alphal00p/vakint-external-data/blob/fced6061d2037272bc76c15bd0396b9ac7287bf7/" + path,
  body,
)
#let producer_evidence(path, body) = link(
  "https://github.com/alphal00p/vakint-external-data/blob/9fbcd9a6bcc63ac39c342b623b394c118ac318b1/validation/" + path,
  body,
)
#let cache_evidence(path, body) = link(
  "https://github.com/alphal00p/vakint-external-data/blob/15270b5bfddf535afc740de87a6d89d597930b00/validation/" + path,
  body,
)

== Releases and implementation identities

- #link("https://github.com/alphal00p/vakint-external-data/releases/tag/v0.1.0")[Immutable release v0.1.0]
  and #link("https://github.com/alphal00p/vakint-external-data/releases/tag/data-latest")[the data-latest alias]
  were published from
  #link("https://github.com/alphal00p/vakint-external-data/commit/fced6061d2037272bc76c15bd0396b9ac7287bf7")[fced6061d2037272bc76c15bd0396b9ac7287bf7].
  Each contains five individually downloadable loop archives, a manifest and
  `SHA256SUMS`. Two complete archive/manifest builds were byte-identical.
- Vakint pins Symbolica, Numerica and Graphica together at
  #link("https://github.com/ValentinHirschi/symbolica/commit/b6cca776f9dce9128fcf10a724228c9445d34070")[b6cca776f9dce9128fcf10a724228c9445d34070].
  This backports the isolated sixth-root HPL implementation while retaining
  Vakint's existing dependency lineage. The implementation is submitted in
  #link("https://github.com/symbolica-dev/symbolica/pull/95")[Symbolica PR 95],
  whose PR branch head is
  #link("https://github.com/ValentinHirschi/symbolica/commit/2bc005909a3be71b16a3aa3e1075073be17d10c9")[2bc005909a3be71b16a3aa3e1075073be17d10c9].
- Vakint pins RustRed at
  #link("https://github.com/alphal00p/rustred/commit/a666e705a92b06ae640a0d9e69b2dfe5bbe748ef")[a666e705a92b06ae640a0d9e69b2dfe5bbe748ef].
  Owned programs, guarded residual preservation and typed factor-coverage errors
  are also submitted against current upstream main in
  #link("https://github.com/alphal00p/rustred/pull/4")[RustRed PR 4], head
  #link("https://github.com/alphal00p/rustred/commit/62886f036eef01cba7d4236b6575faf50417e200")[62886f036eef01cba7d4236b6575faf50417e200].
  Saved-cache codecs and authenticated program-binding schemas are unchanged.

The licensed Symbolica library run passed 834 tests, with one ignored. A fork
coverage run lacked the license secret and aborted; local validation does not
turn that hosted check into a pass. Symbolica's CLA was handled independently
and requires no further action for this delivery.

== Data, provenance and limits

#table(
  columns: (0.6fr, 1fr, 1fr, 1fr, 1.4fr),
  [Loops], [Declared keys], [Numeric coefficients], [Symbolic coefficients], [Completed cache roots],
  [1], [1], [6], [6], [5],
  [2], [2], [13], [13], [77],
  [3], [38], [252], [252], [2,320],
  [4], [74], [468], [364], [164,352],
  [5], [905], [0], [0], [No prepared cache],
)

The release contains 739 numerical and 635 symbolic coefficients. Its 104
numerical-only four-loop coefficients have validated caps of 50 digits (100
coefficients) or 200 digits (four coefficients). Exact rational/HPL/zeta formulas
are evaluated at requested precision; working precision cannot enlarge a finite
source cap or supply missing Laurent orders. Five-loop terminal declarations
provide definitions and campaign provenance, not computed coefficients or a
complete independent master basis. See #link("masters.typ")[master capabilities].

#producer_evidence("producer-provenance-summary.json.gz", [Producer provenance])
was checked for all 52 masters containing the 104
numerical-only coefficients: acceptance hashes match the original export and
the accepted computed files match the runtime numerical payloads byte-for-byte.
These archived proofs are an evidence-only follow-up; release payloads are unchanged.
For P11 (`PR4d`), the
#producer_evidence("producer-4L-P11-validation.json.gz", [accepted producer receipt])
(SHA256
`2809d9a3621549ae19a01ae3f4265adc33ced281f1e8d843b6835e37bc42ab73`)
records 50 stable significant digits, independent refinement and a 70-digit
internal comparison. This validates the extra digits beyond the older 29-digit
FORM literal. A separate
#producer_evidence("opaque-four-loop-audit.json.gz", [legacy-table audit])
found 89 compatible opaque rows and
15 orders absent from that table; it does not use the older table to certify
the producer's additional digits.

The #evidence("provenance/migration.json", [migration receipt]) preserves checksums
for all 10,324 original artifact files. Replaced coefficient payloads remain in
immutable provenance copies. The
#evidence("provenance/coefficient-supplement.json", [coefficient supplement])
adds 49 previously omitted finite coefficients across 42 keys: 25 at 2L/3L and
24 at 4L. No new connected five-loop numerical calculation was performed.

The #evidence("validation/manifest.json", [validation manifest]) hashes archived
source formulas, original decimal references, exporters and their reports:

- All 464 historical HPL projections, representing 324 words, agree with native
  evaluation at both
  #evidence("validation/hpl-parity-200.json.gz", [200 digits]) and
  #evidence("validation/hpl-parity-20000.json.gz", [20,000 digits]).
  #evidence("validation/hpl-ginac-reference.json.gz", [Independent GiNaC checks])
  cover three weight-six values at 220 digits. Supported arguments are the two
  sixth roots `exp(±i*pi/3)`, with alphabet `{-1,0,1}` through weight six.
- The #evidence("validation/coefficient-supplement-validation.json.gz", [2L/3L proof])
  and #evidence("validation/coefficient-supplement-four-loop-validation.json.gz", [4L proof])
  compare exact formulas, native/legacy HPL values and independent retained
  MATAD numerical coefficients. The 4L audit includes every stored FMFT numeric
  entry, all 74 mappings, 22 graph-factorization certificates and all 132 existing
  factorized symbolic/numerical coefficients. All 49 added formulas
  #evidence("validation/coefficient-supplement-serialization-validation.json.gz", [parse and round-trip]).
- #evidence("validation/shared-constant-validation.json.gz", [Shared-constant checks])
  preserve the finite literal precision and rounding evidence. Removed HPL data
  and comparison values remain archived; FORM master formulas remain in Vakint.

Cache preparation reports for
#evidence("1L/cache.json", [1L]), #evidence("2L/cache.json", [2L]),
#evidence("3L/cache.json", [3L]) and #evidence("4L/cache.json", [4L])
all report their selected finite workloads complete. Guarded intermediate capture
published 735 3L and 195,449 4L records; 703 and 73,274 of those later served
selected roots. The complete packs therefore contain 2,352 and 286,527 distinct
records, including 32 and 122,175 outside their selected primers. One 4L root reached
the optional capture-metadata budget. This limits intermediate capture, not root
completion; neither these workloads nor their caches certify arbitrary-index
IBP coverage. Cached/uncached 1L–4L parity passed with the final native bindings.

== Published cache workload and size

The offline producer was run separately for 1L, 2L, 3L and 4L. The public
reproduction driver is #link("examples/prepare_bundled_caches.rs")[`prepare_bundled_caches`]:
it calls `prepare_cache` with the loop count and default `Primer` options.
The input is the complete finite primer, spanning all physical parent pinches,
at most two total extra positive powers above the surviving unit lines, and
numerator rank at most one. At 4L the parents are H, FG, BMW and X. This covers
common unit, dotted, pinched and rank-one requests while harvesting complete
guarded intermediate reductions along their paths. It is independent of the
later illustrative requests in `native_cache.rs` and of certification coverage.

Each container is `nL/cache/`, containing `bundle.rrindex`, one
`bindings/<binding-id>/cache.rrindex` per program, and hash-named
`chunks/*.rrchunk.gz`. These store exact unit-mass, unnormalized reductions and
nonzero guards; numerical master values are separate data. Writer locks and
optional preparation checkpoints are excluded from distribution.

#table(columns: 6,
  [Loops], [Selected roots], [Distinct records], [Native file bytes], [Metadata bytes], [Preparation seconds],
  [1], [5], [5], [1,235], [5,964], [0.00473],
  [2], [77], [77], [3,854], [6,445], [0.00747],
  [3], [2,320], [2,352], [160,651], [22,749], [1.76790],
  [4], [164,352], [286,527], [90,887,721], [4,564,269], [653.07154],
)

Sizes are logical file-content lengths: 9,863 native files total 91,053,461
bytes, plus 4,599,427 bytes of `cache.json` generation metadata. Native chunk
payloads are already gzip-compressed; filesystem allocation is a separate
quantity. The 4L cache comprises 9,812 chunks and five indexes. Its native
preparation took 10 min 53 s, or 11 min 9 s including the producer's startup
and publication work. These generation times come from the published producer
receipts at Vakint `e98b51c489ab9a35d01f30aa865ad79c48372624` and RustRed
`74a83e94f9f841ae113f3d151cd80567dbc9805d`; final-runtime compatibility was
validated separately through codecs, bindings and exact cached/uncached parity.

The complete compressed loop assets are 48,045, 107,633, 370,372 and
111,922,455 bytes respectively. They include IBPs, master coefficients and graphs
alongside these caches. No five-loop prepared cache is published.

The #cache_evidence("prepared-cache-inventory.json.gz", [independent inventory]),
#cache_evidence("prepared-cache-native-file-inventory.json.gz", [complete file list])
and #cache_evidence("prepared-cache-metadata-authentication.json.gz", [manifest authentication])
verify all 9,863 native file hashes and sizes. This evidence-only follow-up
does not change the immutable runtime release.

=== Measured reduction benefit

The #cache_evidence("prepared-cache-reduction-timings.json.gz", [complete timing receipt])
and #cache_evidence("prepared-cache-reduction-protocol.json.gz", [protocol/profile])
record one fresh cached/uncached process pair per request on two fixed cores,
using the unoptimized Cargo test profile with debug information. Public exact
RustRed Retain calls use unit mass and FMFT/MATAD normalization. Engine/data/program
initialization, cache attachment, three RAM repeats and exact-reference checks
are outside the first-call timings; numerical master and HPL evaluation are
excluded. All 22 processes exited zero, and all eleven cached expressions
agreed exactly with uncached references. These are shared-host samples, with
no OS page-cache flush or statistical release-build speedup claim.

#table(columns: 5,
  [Family], [First power], [Uncached seconds], [Prepared seconds], [Uncached / prepared],
  [1L], [3], [0.010444], [0.010304], [1.01],
  [2L], [2], [0.012492], [0.012369], [1.01],
  [3L], [2], [0.254533], [0.094634], [2.69],
  [4L/H], [2], [34.154138], [0.123183], [277.26],
  [4L/FG], [2], [2.537141], [0.135310], [18.75],
  [4L/BMW], [2], [0.056794], [0.064122], [0.89],
  [4L/X], [2], [157.714177], [0.318779], [494.74],
  [4L/H], [4], [47.442782], [24.538869], [1.93],
  [4L/FG], [4], [4.582351], [0.919546], [4.98],
  [4L/BMW], [4], [26.793927], [6.579722], [4.07],
  [4L/X], [4], [249.320210], [81.082004], [3.07],
)

Every other physical denominator has power one and auxiliary powers are zero.
The seven primer requests have backing entry hits and zero cached native rule
applications. The four power-four requests have no exact entry hit: H, FG, BMW
and X instead use 2,983, 52, 888 and 8,460 descendant hits respectively. Thus
intermediate records help inputs outside the selected primer.

The BMW primer already needs zero ordinary native rules; its small slowdown
shows no cache benefit for that request. The sampled 1L and 2L calls are also
effectively unchanged. Common first-use four-loop data/program setup took
74.91–98.93 seconds per fresh test process; attaching its cache took
0.140–0.168 seconds. Those separate costs matter for a one-shot calculation
and are amortized when the loaded engine handles multiple requests.

The #cache_evidence("prepared-cache-reduction-independent-review.json.gz", [independent review])
verifies the raw logs, profile, source identity, expression assertions and all
ratios. Its scope is separate from the earlier exact native guard checks.
Benchmark sources and raw process logs accompany the archived receipt.

== Completed runtime checks

- The #link("examples/official_external_data.rs")[published-release example]
  downloaded and activated all 1L–5L packs and the 1L–4L prepared caches from
  release v0.1.0. A separate process under
  `unshare --user --map-root-user --net` reused the published installation with
  no external network interface and returned identical 1L numerical output.
  Version, loop selection, cache availability and downloaded/reused status
  were machine-checked. These are functional checks, not timing benchmarks.
- The four-loop scientific test helper now shares one official-data engine and
  its resolved path. It requests maximum loop count four, version `latest` and
  `refresh=false`. With `VAKINT_TEST_EXTERNAL_DATA` absent, the path option
  remains `None`: the loader uses the actual test process working directory,
  rather than a repository or build-directory default. An explicit override
  still accepts a checkout or managed installation base. Ordinary vendored
  construction remains offline unless an external context is explicitly chosen.
  Native Unix is the supported official-download platform; the three scientific
  entry points retain their original assertions, tolerances, ignores and FORM
  requirements. See #link("external_data.typ")[test loading and execution].
- The #link("src/rustred_evaluation/masters_acceptance_tests.rs")[four new-order acceptance tests]
  passed for 2L sunset epsilon^4, 3L `miBN` epsilon^4, the 3L five-line integral
  epsilon^3 and external 4L `PR2` epsilon^3. They check loaded metadata,
  symbolic/direct numerical agreement, 70-digit parity at 80-digit working
  precision and errors at the next unavailable order.
- Ten #link("src/rustred_evaluation/official_tests.rs")[official-loader tests]
  passed: offline reuse, redirects, refresh, interrupted-download recovery,
  path/resource rejection and transactional rejection of a checksum-valid pack
  with the wrong native cache binding. Additional assertions passed for vendored
  2L/3L fallback after smaller external loads. The complete inflated-stream
  bound, including internal GNU/PAX extension metadata, also passed.
  Incomplete local checkouts fail without network access or filesystem mutation;
  unrelated valid manifests and empty loop-named directories allow managed installs.
- Focused #link("src/rustred_evaluation/masters.rs")[master-registry checks]
  passed for exact cancellation, coefficient products, missing orders, finite
  precision, oversized Float annotations, missing validated caps, tagged HPLs,
  nonfinite data and snapshot-bound coefficients. Fifteen combined numerical,
  #link("src/rustred_evaluation/experimental.rs")[four-loop finalizer] and loader
  checks passed; a subsequent capability selection passed 14 tests; the final native matrix
  passed all 16 master-capability checks.
  These include exact alias handling without masking malformed data or discarding
  a richer available representation. Scalar MATAD comparison also passed after
  independently checking exact zero poles before numerical comparison.
  A subsequent 24-test HPL/catalog/new-order selection passed after preserving
  real HR/HI projection nodes through symbolic expansion; nonlinear projection
  polynomials have exactly zero imaginary output, without a numerical zero clamp.
- The complete fifteen-case four-loop acceptance passed both the unchanged
  historical FORM/frozen-reference checks and independent FORM-retained/
  external-coefficient comparisons against RustRed at `1e-30`. The FORM-free
  H rank-four test also passed its separately produced, source-hashed
  #link("tests/inputs/h_rank_four_external_reference.json")[external-data reference].
  All sixteen expanded-numerator/pinch pairs passed for H, FG, BMW and X with
  the original `1e-20` threshold and unchanged algebraic/routing assertions.
- Eight #link("src/rustred_evaluation/cache_transfer_tests.rs")[cache-transfer tests]
  and the #link("src/rustred_evaluation/cache_tests.rs")[fresh-process 1L replay]
  passed. Fifteen materialization/numerator/routing checks passed, including
  exact cancellation and partial-coverage behavior. Cache reads remain optional
  and evaluation does not modify prepared packs.
- RustRed's pinned integration passed 559 candidate/routed/guarded-cache tests,
  53 saved-application/cache tests (one ignored), and 24 product tests. Its
  main-based PR passed 596 candidate/routed/cache tests (one ignored), the same
  53 saved-application/cache tests, and CLI compilation. Tests distinguish
  preserved residuals from hard guard/resource failures and keep partial results
  out of strict completed caches.

The recorded local all-targets Clippy pass completed with `-D warnings` on the
Float, alias and loader-bound changes. The complete repository Nix/upload result
is recorded below.
The initial offline locked-resolution proof matched all 709 Git dependency/
feature nodes against the Nix vendor output. The later payload-adapter migration
removed the obsolete separate payload package and reused the already-pinned
unified plugin; the coordinated Symbolica/Numerica/Graphica revision is unchanged.
Its 14 display tests and three integration tests passed, including native Atom
round-trip and rejection of another embedded export version. The standalone
Tydenso Wasm-target source check and native-feature exclusion check also passed.
Existing compressed Wasm assets were not regenerated, and no rebuilt-asset
interoperability result is claimed.

The final Cargo package inventory contains 316 files, with no removed
`artifacts/` tree, prepared caches or 4L/5L runtime packs. Building the crate archive
with `cargo package --allow-dirty --no-verify` stops at the existing Git-only
`rustred` dependency without a crates.io version requirement. No registry version
was invented. Publishing Vakint itself is outside this migration's scope;
workspace compilation and the published external-data packs are the delivery
targets.

== Consolidated workspace dependencies

Vakint now uses the workspace FeynKit, Spenso and Idenso packages. This removes
separate Git-sourced copies whose Symbolica state initializers collided in a
combined GammaLoop process. Symbolica/Numerica/Graphica and RustRed retain the
pinned revisions above. Required integral-family and symbolic-kinematics APIs
were incorporated from the current FeynKit dependency branch, including exact
Gaussian-rational coefficients and tiny-coefficient regressions. Stored graph
routing and physical propagator order are preserved at the local graph boundary.

The consolidated implementation passed 70 graph and 30 kinematics tests, six
combined GammaLoop precision/tensor/settings regressions, and all 52 active
FeynKit tensor tests (one existing ignored test). Twenty Python subprocess
checks passed: four graph-ingress/serialization, five external-data and eleven
native-cache checks. A fresh-process tensor JSON/DOT regression verifies exact
registered representation metadata and rejects altered tags; the previous
extension reproduces its failure. The importer reuses exact registered symbol
declarations, following the existing dependency implementation.

Actual offline Nix-vendor resolution agrees with Cargo for all 695 selected-CI
nodes and 822 all-feature nodes, including dependency edges and enabled features.
Hakari and generated CI graph checks pass. The obsolete Git-Spenso vendor
augmentation was removed after that proof. All required new fixture and stub
inputs were checked in the Nix source. The rebuilt Python extension downloaded, validated and activated all five
published packs in a fresh installation. A separate
`unshare --user --map-root-user --net` process reused that same installation and
produced exactly the same 80-digit one-loop output as a new baseline made with
the same extension. Download, baseline and isolated reuse all passed every
assertion and exited zero without recorder failures. The module SHA256 was
`ac081c32cfdbc4cda27b8e262ca86d874cd311fde8eff2e1a13cc02cb7654abb` throughout.

The consolidated scalar/tensor integration and collected-power regressions
passed, including all affected tensor frontend cases. Five old/new tensor
prepass comparisons have exactly zero algebraic differences. The revision-specific
four-loop numerical results below retain the original scientific assertions and
tolerances.

The consolidated algebra adapters required two corrections already represented
in the original dependency behavior. Collected positive powers must reach
callbacks as complete monomials; otherwise a repeated vector can be dropped.
The local QCD fixture adapter also expands metric sums around gamma/color
contraction. Independent Clifford pairings agree through ten gamma matrices
(945 pairings), and the targeted six-gamma regression passes. Neither correction
changes FORM formulas or acceptance tolerances.

Closed ordinary gamma traces also reuse the existing dimension-generic repeated
word rules before pairing recursion, following the upstream word-builder design.
All 235 selected Idenso package tests pass, including the three new repeated-slash,
formal-dimension and refusal tests; the existing ignores and default filters are
unchanged. The captured physical twenty-gamma expression agrees exactly with the
original dependency result after explicitly mapping only its old scalar DOT head
declaration to the local DOT declaration. Both preserve the same symmetric,
linear dot product and every argument; production registration checks are
unchanged. The local probe completes in about 0.26 seconds; this is a contraction
measurement, not an end-to-end QCD timing claim. The complete QCD scalar comparison passed with exact symbolic equality, as
recorded below.

The QCD fixture contracts vector coefficients exposed by expansion before gamma
collection. The independent four-gamma control and all five fixture controls
pass. Complete original and candidate preparations each produced 8,306 scalar
terms, 21 contributions and five recorded guards. Their exact irreducible
nonzero predicates agree: `D`, `D-1` and `D+2`. The complete strict Rust comparison
passed after both inputs were admitted and routed to
`I5L_32745_pinch_8_9_12` with unknown-integral fallback disabled. Exact symbolic
subtraction and rational combination returned zero. This used no IBP reduction,
numerical tolerance or DOT declaration rewrites in the complete inputs; the
explicit compatibility mapping used by isolated trace probes had zero matches
here. Source, library, binary and input hashes bind the original preparation to
`93ce035278ab1feb93e0e6215ad85eb79dbc5e75` and the candidate/comparison to the
`c398e309de1767fcecbe841f6555080e10700fb6` runtime source cohort.

A call-path audit against `d5383061f075c3762bfb2ecdc7affb03f50e6e3b` confirms
fourteen called or related method bodies and twenty-one transitive source/pin
scopes are unchanged at that revision. The sealed preparation and strict
canonicalization do not enter the indexed-scope path; the later UV-restoration
and FORM-I/O repairs are also outside those paths. This is source-path impact
evidence, not a rerun with final binaries. Numerical RustRed evaluation enters
the changed tensor prepass and is covered by the separate numerical gates below.

The complete native registry/cache/capability matrix passed all 64 selected tests,
including actual five-loop witness admission/reuse. The fresh 38-root four-loop
campaign completed with 30 newly stored roots and eight strictly reused captured
descendants. Every selected root has a validated durable record, all four cold
term/ordered-guard controls agree, and fresh-process value/mass/hit/limit checks
passed. Its two optional metadata capture exhaustions are recorded separately
from successful root completion.

A separate authentic five-loop catalog test passed after validating all 905
terminal records against the immutable v0.1.0 manifest. Retain preserves a real
family-qualified terminal and its ordered powers. Direct numerical finalization
and retained-catalog numerical resolution both reject it with its identity,
package version, declared/uncomputed status and existing parent-topology DOT
path. Both routes return exact zero when rational epsilon coefficients cancel
that terminal. This check loads no IBP application and makes no claim of
successful numerical five-loop evaluation. It is an additional check beyond the
64-test native matrix. Formatting and all-target Clippy also passed.

Fifteen inherited Idenso snapshots were refreshed after exact validation of
seventeen old/new expression pairs and a bijection preserving all six external
DOT slots. Three fresh processes reproduced the same new port order. All 235
selected Idenso package tests and the 578-test Nix Spenso group pass, with
existing skips unchanged. The refresh changed no production algebra or CI
selection.

Generated Python signatures now follow the tested
`symbolica.community.hepkit.vakint` module and include cumulative loaders, master
modes and the existing diagram-ingress function. Declared, stub and runtime
surfaces agree on all five exports. Stub regeneration/check, all 23 catalog
tests, formatting and relevant all-target Clippy passed. The runtime topology
reference still contains 27 entries; native graph isomorphism and both original
and final Python hosts confirm that its refreshed FG contraction labels describe
the same existing graph classes and supported short forms.

A further native/compact DOT discriminator defect was reproduced and submitted
against the current FeynKit branch in
#link("https://github.com/alphal00p/gammaloop/pull/156")[PR 156]. Its 123 graph
tests and all-target Clippy pass. Its exact final-source checks/upload command
completed successfully at `d39d9bf206d8f858efedc9281f034cddb50c51fa`; the PR was pushed,
marked non-draft and labeled `final-review`. Its
#link("https://nix-ci.com/gh:alphal00p:gammaloop/reject-native-dot-metadata-fallback/d39d9bf206d8f858efedc9281f034cddb50c51fa")[hosted NixCI check also passed]. Already-fixed dependency issues were
incorporated without duplicate upstream submissions.

Native multiattachment tensor-port preservation was submitted against the current
FeynKit branch in #link("https://github.com/alphal00p/gammaloop/pull/157")[PR 157],
head `070556bcab0eff348188f27b0dfc44c173ffee46`. All 50 focused unit tests,
Clippy, formatting and 74 architecture-scope checks passed. Its exact
final-source static/producer checks and upload exited zero with
`CI upload completed`; the PR is non-draft and labeled `final-review`.

FORM saved-index epsilon restoration was submitted in
#link("https://github.com/alphal00p/gammaloop/pull/158")[PR 158], head
`2e231c0e2ee6967c2ab61b070ecbbe654c5391bf`. All 48 focused unit tests, Clippy,
formatting and 74 architecture checks passed. Its exact final-source
static/producer checks and upload exited zero with `CI upload completed`;
the PR is non-draft and labeled `final-review`. Its tested base is `8cb80dfd`;
the target subsequently included the independent PR 156 graph fix at `d39d9bf2`,
with a clean merge verified separately. No merged-base rerun is claimed.
The latest dependency already contains the registered symbol-data parser fix;
this submission does not duplicate that repair. No hosted-CI result is claimed
for PRs 157 or 158.

Python build producers now use the complete manifest normal-dependency closure,
including optional UFO and FeynKit Python sources. This follows the existing
upstream source-closure fix: the earlier native-only closure supplied a dummy UFO
crate even though Python ABI features activated its unconditional `UfoLoader`
export. The evaluated corrected source contains all twenty dependency packages
with matching source hashes. The corrected dependency producer compiled the
real UFO loader and GammaLoop API successfully. No Rust API or Cargo manifest
was changed.

Two final serialization corrections follow existing upstream fixes. Hedge metadata is
passed unquoted to Linnet, which owns DOT quoting; the final split-forward
round-trip assertions remain unchanged. FORM input uses the upstream exact-tree
serializer, including explicit complex components, to preserve precedence,
symbol identities and unit imaginary coefficients. FORM scripts and master
formulas are unchanged. Two formatter regressions pass across 34 exact
round-trip cases. Four existing display expectations within the two-loop ghost
test omit the unit numerator before the imaginary unit, matching the current
printer. Each complete old/new expression is an identical exact Atom; all
other bytes of that test source, including numerical and forest-size assertions,
are preserved. These serialization paths do not alter the completed QCD
contraction or RustRed coefficient evaluation.

A subsequent UV integration regression exposed bare source-index labels beside
already typed Lorentz slots after Vakint evaluation. The existing upstream
restoration boundary supplies the missing inverse to numerator projection:
returned momenta and metrics regain their source Lorentz representation,
existing typed dimensions and labels are preserved, and numerical/dot-generated
indices keep their topology-qualified chart. Two focused restoration controls
pass, including an independent scalar identity for the previously incompatible
mixed tensor/dot sum. The complete two-loop ghost regression also passes,
including all six graphs, unchanged scientific expression checks, exact
legacy/new forest comparisons and original forest-size assertions. Matching
all-target Vakint/GammaLoop Clippy passes with `-D warnings`. No tensor parser
or safety check is bypassed.

Closed indexed scalar powers also require independent dummy scopes. The existing
upstream composite-power freshening pass preserves original indices and creates
new labels only for dots opened by that conversion. A factorized Lorentz owner
validates original multiplicities and closes raw scalar bases before collection;
validated pure-loop invariants then remain exact scalar coefficients through the
unchanged checked tensor reducer. Native tensor components, expanded source
ports and mixed-interface sums retain their provenance and checked reduction.

All 43 focused/frontend checks pass: eight owner controls, six original numerator
scope tests, three conversion/legacy/alias controls, two FORM serializer controls
and all 24 existing frontend cases. This includes denominator power eleven,
exact parent-key transport, alias collision and failure rollback, the existing
rank-20 acceptance/rank-22 refusal, and genuine ambiguity errors. Relevant
Clippy, formatting and independent owner review pass. No original numerator or
frontend assertion, tensor rank/expansion limit, FORM source, master table,
FeynKit/Idenso/Spenso implementation or dependency pin changed in this repair.
All three numerical gates also pass on the corrected `00ed8877` source: the
sixteen numerator/pinch pairs, complete fifteen-case four-loop inventory and
independent FORM-free H rank-four reference. Their reference data, scientific
assertions and tolerances are unchanged.

== Final metadata and UV integration checks

The native UV boundary now restores the known regulated projector dimension to
the retained source Lorentz representation before the backend's Laurent
expansion. Scalar epsilon coefficients, including evanescent factors multiplying
poles, still participate in that expansion; tensor dimension metadata does not
become a formal derivative. Source labels, unrelated represented dimensions and
Vakint's momentum/metric ABI remain intact. The existing strict tensor-network
slot checks are unchanged.

FORM's saved index table now returns in the caller's epsilon convention, which
also preserves vector-mapping source identity. Its parsing boundary reuses the
existing registered symbol declarations so attributed native indices round-trip
without suppressing incompatible-declaration errors. Thirteen focused checks
pass: eight FORM/restoration controls, the two raised-cut selection integrations,
the threshold-cut association integration, `addbar`, and the complete six-graph
two-loop ghost regression. The controls include both index-substitution modes,
metric/vector mappings, altered-tag refusal, derivative absence and the exact
finite evanescent contribution. Existing scientific/forest assertions remain
unchanged. Relevant all-target Clippy, formatting and diff checks pass; all
fourteen checked FORM scripts, templates, master/HPL sources are byte-identical.
The retained receipt is `uv-lorentz-metadata-final-validation.json` (SHA256
`c88b8d2fa37963215a7fb5fd186933fbc558798c19b1ed7968a11d661a6d3427`).

Eight tracked DOT fixtures also required eleven model-fingerprint metadata
updates. The two old fingerprints are reproduced exactly by restoring only the
former explicit unit numerator (`1𝑖` versus `𝑖`) in the serialized models. Reparsing those models
returns the same canonical data, and their original source assets are unchanged.
Every other fixture byte, including topology, algebra, cuts and indices, is
preserved. All eleven affected diagrams pass strict import and exact round-trip;
the old fingerprints remain rejected. These are metadata refreshes, not relaxed
model validation.

The GL20 local-event snapshot changes only 53 numerical leaves; 694 leaves,
all keys/array order, metadata, sample point and cut/event identity are unchanged.
The public f64 API reproduces three fresh archived candidate outputs exactly.
A forced arbitrary-precision run uses the same algorithm at 1,000 *bits*, with
an f64 result boundary. At this single ill-conditioned point, old and new
snapshots have the same worst complex-relative difference from that reference,
about `2.62237e-12`; their largest mutual change is about `59.4068` double
machine epsilons relative to complex scale. This is a precision-escalation
comparison, not an independent physics oracle, a general forward-error bound or
per-component ULP validation. No small imaginary component is clamped or
filtered. After independent review, only expected numerical data was refreshed;
the original test, strict full-precision assertion and tolerances are unchanged.
Three fresh strict snapshot replays pass with `INSTA_UPDATE=no`.

All 53 declared architecture scopes match their exact registered hashes. The
follow-up metadata/UV repairs add no runtime files and leave the recorded
316-file package inventory unchanged.

The current CWD-default four-loop rerun also passes all three scientific gates
in fresh network-isolated processes, with `VAKINT_TEST_EXTERNAL_DATA` unset.
Maximum four, `latest` and `refresh=false` resolve v0.1.0 beneath the actual
working directory; each process reports validated reuse and available caches
for all four packs. The numerical references, assertions and tolerances are
unchanged. This separately validates the new official-loading test path after
the earlier `00ed8877` numerical run. Its retained result receipt is
`official-cwd-four-loop-results.json` (SHA256
`7ea7dac9e6570d96a112536969e64fd52770658e01f30199e2cf0e5fa7f9eac1`).

== Final repository-CI and delivery gates

The exact `just ci-checks-and-upload` run at
`00ed8877f629f8bb42578813c8a286104cb5e3fa` passed static checks and its complete
Vakint runtime selection: 343 passed, 24 skipped. The core selections passed
341, 19 and 436 tests respectively, and the Python selection passed all five.
The integration group ran 108 tests: 98 passed and ten failed, with 67 skipped.
Those failures exposed the fingerprint metadata, tensor-interface and numerical
snapshot issues addressed above. The aggregate run exited one before upload;
319 valid producer outputs were retained. Successful groups and partial
documentation test output do not make the failed aggregate or its upload complete.

The corrected runtime source at
`d5383061f075c3762bfb2ecdc7affb03f50e6e3b` passed the exact
`just ci-checks-and-upload` command with exit zero and `CI upload completed`.
Clippy, formatting, workspace graph checks and doctests passed. All nine runtime
groups have valid successful results: 2,723 selected tests passed, comprising
344 Vakint, 108 integration, 799 core, 256 documentation, five Python and 1,211
other native package tests. Existing skips were preserved. Checks took 2,784
seconds, publication preparation 363 seconds and upload 475 seconds.

The complete source-bound receipt is
`coherent-main-ci-native-metadata-cwd-upload.receipt.json`, SHA256
`0a4ce20114955def3dd8e910bba7c321f5189f270fa630d7b2ae8107f13d68fc`.
Its nine runtime logs and valid producer outputs are retained locally. The
separate four-loop scientific checks and published-download/offline-reuse checks
above remain part of acceptance; they are not replaced by this CI selection.

The final report-only commit uses the same exact checks/upload gate before a
normal push to `origin/vakint_rustred`, without a main-repository PR. Separate
upstream submissions retain their own checks/upload, PR and review-label
requirements. Upstream merge timing is not a delivery requirement.
