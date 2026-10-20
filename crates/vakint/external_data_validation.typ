= External-data migration validation

This records completed migration checks as of 20 October 2026. The focused
suites overlap; their counts are not a combined whole-repository test result.
Outstanding end-to-end and repository-CI checks are listed separately below.

#let evidence(path, body) = link(
  "https://github.com/alphal00p/vakint-external-data/blob/fced6061d2037272bc76c15bd0396b9ac7287bf7/" + path,
  body,
)
#let producer_evidence(path, body) = link(
  "https://github.com/alphal00p/vakint-external-data/blob/9fbcd9a6bcc63ac39c342b623b394c118ac318b1/validation/" + path,
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

Both upstream PRs were open at this snapshot. Symbolica's fork coverage job
lacked the license secret and aborted; the
licensed local library run passed 834 tests, with one ignored. This does not
turn the failed hosted check into a pass. RustRed's hosted NixCI check was pending.

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
has SHA256
`2809d9a3621549ae19a01ae3f4265adc33ced281f1e8d843b6835e37bc42ab73`
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
stored 735 additional 3L and 195,449 additional 4L records. One 4L root reached
the optional capture-metadata budget. This limits intermediate capture, not root
completion; neither these workloads nor their caches certify arbitrary-index
IBP coverage. Cached/uncached 1L–4L parity passed with the final native bindings.

== Completed runtime checks

- The #link("examples/official_external_data.rs")[published-release example]
  downloaded and activated all 1L–5L packs and the 1L–4L prepared caches from
  release v0.1.0. A separate process under
  `unshare --user --map-root-user --net` reused the published installation with
  no external network interface and returned identical 1L numerical output.
  Version, loop selection, cache availability and downloaded/reused status
  were machine-checked. These are functional checks, not timing benchmarks.
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
  checks passed; the subsequent final capability selection passed 14 tests.
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
Float, alias and loader-bound changes. The repository Nix/upload gate remains
pending below.
The initial offline locked-resolution proof matched all 709 Git dependency/
feature nodes against the Nix vendor output. The later payload-adapter migration
removed the obsolete separate payload package and reused the already-pinned
unified plugin; the coordinated Symbolica/Numerica/Graphica revision is unchanged.
Its 14 display tests and three integration tests passed, including native Atom
round-trip and rejection of another embedded export version. The standalone
Tydenso Wasm-target source check and native-feature exclusion check also passed.
Existing compressed Wasm assets were not regenerated, and no rebuilt-asset
interoperability result is claimed.

The final Cargo package inventory contains 315 files, with no removed
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
inputs were checked in the Nix source. The final coherent Python extension downloaded and activated all five published
packs. A clean `unshare --user --map-root-user --net` process reused that same
installation and produced identical 80-digit one-loop numerical output, with no
network interface. Development recorder path-serialization failures were retained
and corrected; the final isolated run passed every assertion and exited zero.

The consolidated numerical rerun passed the complete fifteen-case four-loop
inventory, all sixteen numerator/pinch pairs and the independent FORM-free
H rank-four reference with unchanged assertions and tolerances. A 48-test
scalar/tensor integration selection passed; after the collected-power correction,
all 24 affected tensor frontend tests and all 52 active tensor unit tests passed
again. Both focused collected-power regressions passed, and five old/new tensor
prepass comparisons have exactly zero algebraic differences. The final Python
extension passed all twenty subprocess checks. The suites overlap.

The consolidated algebra adapters required two corrections already represented
in the original dependency behavior. Collected positive powers must reach
callbacks as complete monomials; otherwise a repeated vector can be dropped.
The local QCD fixture adapter also expands metric sums around gamma/color
contraction. Independent Clifford pairings agree through ten gamma matrices
(945 pairings), and the targeted six-gamma regression passes. Neither correction
changes FORM formulas or acceptance tolerances. The complete native registry/cache/capability matrix passed all 64 selected tests,
including actual five-loop witness admission/reuse. The fresh 38-root four-loop
campaign completed with 30 newly stored roots and eight strictly reused captured
descendants. Every selected root has a validated durable record, all four cold
term/ordered-guard controls agree, and fresh-process value/mass/hit/limit checks
passed. Its two optional metadata capture exhaustions are recorded separately
from successful root completion. Full five-loop QCD scalar input parity is still
pending its complete run.

Fifteen inherited Idenso snapshots were refreshed after exact validation of
seventeen old/new expression pairs and a bijection preserving all six external
DOT slots. Three fresh processes reproduced the same new port order. All 232
Idenso tests and the complete 575-test Nix Spenso group passed, with the existing
skips unchanged; production algebra and CI selection were not modified.

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
tests and all-target Clippy pass; the PR's separate Nix upload/readiness gate is
in progress. Already-fixed dependency issues were incorporated without duplicate
upstream submissions.

== Final repository-CI result — pending completion

To be filled after the final Nix build, required checks, cache upload and push
complete. Record the tested Vakint commit, Nix result paths, upload result,
hosted check URLs and any remaining failures. No whole-repository CI pass or
final-review readiness is claimed by the focused evidence above.
