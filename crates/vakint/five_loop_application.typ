= Symbolic five-loop saved applications

Successful external-package loading prepares the 67 TIDE A5 representatives
(48 nonfactorized and 19 factorized), using four twelve-line parents and explicit
physical-to-family slot maps. Ordinary `Vakint::new()` constructs only the
one-through-four-loop catalogue. The registry is topology input data, not a claim that every power
or numerator is covered by an application package. The 15 common coordinates
include auxiliary numerator coordinates. All numerator transport uses the
same retained simultaneous momentum routing as the propagators.

At the RustRed bridge, closed indexed scalar contractions are restored with
Vakint's existing FeynKit tensor service before any numerator expansion. The
bridge permits direct contractions but rejects any remaining projector rank;
it does not perform an additional tensor reduction. Parent momentum maps act
simultaneously on native scalar-product heads, preserving the scope of powers
such as `(k1.k2)^2` and `(q^2-m^2)^n`.
Bare loop-loop scalar products use FeynKit's native compact invariant syntax;
they are not opened into artificial free tensor indices. Mixed loop-spectator
products and genuinely open tensors keep the existing native projector path
and its rank allowance. The scalar numerator's degree is not a projector rank.

Five-loop long forms use native momentum-family equivalence, including the full
denominator permutation and a Jacobian of plus or minus one. Thus different
attachments of the same vacuum factors can share an integral class. Parametric
polynomial equality only proposes a line permutation. The existing Symbolica
basis solver tests the remaining independent-edge orientations; FeynKit must
verify the complete momentum map for numerator transport. No polynomial-only
identification is accepted. The orientation allowance is an operational limit.

Contractions use one graph quotient before rewriting propagator incidence.
This matters for shared endpoints: for census sector 31746, contracting the
parent's q9 and q10 requires the identifications `(8,3)` and `(8,4)`. Sequential
renames lose the second identification after vertex 8 disappears. The complete
equivalence class must merge vertices 3, 4 and 8. The resulting six-edge graph
has two vertices and five loops. Loop bases are complements of deterministic
native spanning forests, including all disconnected components.

== Selecting an external package

No five-loop application artifact is embedded or downloaded. Set
`VakintSettings::rustred_artifact_path` to a locally hosted exported package,
or preload it with `Vakint::with_rustred_application(&mut settings, path)`,
which also sets that settings path. An instance authenticates the selected package
once and reuses its owned memoization; clones share the synchronized cache.
Changing the selected path loads another immutable application. Removing the
settings path disables five-loop application even if a package is cached.
One-through-four-loop behavior does not require this setting.
Matching an existing lower-loop or custom topology does not load a configured
artifact merely because its path is present. A cold five-loop matching request
loads the selected artifact before attaching the immutable five-loop catalogue.
`vakint.topologies.five_loop_catalogue_statistics()` is `None` before successful
publication and reports its topology count and preparation time afterward.

`VakintSettings::rustred_application_options` controls the native application's
operational allowances. The default coefficient cache retains RustRed's 1 GiB
and 16-million aggregate term ceilings. Its default policy evicts inactive,
completed results in least-recently-used order. Entries needed by pending
combinations or the current result stay protected. A protected working set or
single entry that cannot fit still produces a resource error, not a missing-rule
or closure failure; large admitted workloads can require a larger allowance.
A different path or different options reloads
the application, so an older cache cannot bypass tightened allowances.

The loaded package also owns a separate bounded frontend witness memo. It
retains the exact witness selected by the unchanged minimum-key matcher, keyed
by the full integral and scalar numerator, relevant lowering/matching settings,
and the actual immutable registry/catalogue incarnation. Package selection and
settings validation precede lookup; a hit reruns native admission and compares
the accepted key before use. Unsupported fallbacks and errors are not cached.
The coefficient cache's LRU/LFU policy and numerical results are unchanged.

This metadata memo has at most 32 entries and 16MiB of tracked encoded Atom
payload. The byte counter is not a total RSS estimate: graph, pattern, container,
index-vector and shared-symbol overhead are excluded. Oversized entries are
skipped, and old entries are evicted; optional memo pressure does not fail a
reduction. `rustred_application_witness_cache_statistics()` reports its distinct
hit/miss, invalidation, eviction and payload counters. The existing
`clear_rustred_application_cache()` also clears retained witnesses. A poisoned
optional memo is bypassed, while failures from actual native admission propagate.

For a byte-governed cache in Rust, set `max_cached_coefficient_bytes` on the
options' `reduction_limits`, and set `max_cached_coefficient_terms` to
`usize::MAX`. In Python, `rustred_cache_bytes` selects that behavior;
`rustred_cache_terms` can explicitly add a separate aggregate term ceiling.
For example `rustred_cache_bytes=32*1024**3` requests a 32 GiB coefficient cache.
None of these settings enlarge the saved mathematical entry union.

Rust callers select the policy through
`settings.rustred_application_options.cache_policy`, using
`vakint::RustRedCachePolicy::EvictLeastRecentlyUsed` or `FailOnLimit`. Python's
`rustred_cache_policy` accepts `"evict-lru"` (default) or `"fail-on-limit"`.
This changes reuse only. The legacy raw core reducer's default remains
fail-on-limit, and the separate lower-factor cache retains its existing policy.

`settings.rustred_scalar_numerator_limits` controls caller-numerator lowering
in the saved five-loop frontend, both when choosing a verified witness and when
applying it. Its native defaults are unchanged. The charged work allowance is
`exact_algebra.max_term_operations` (16 million by default); Python exposes the
optional `rustred_scalar_max_work` keyword. This is a per-lowering-call allowance,
not a process-memory budget or a total across all matching candidates. Raising
it never changes the entry union. Tightening it is checked before cached rule
results can be used, without reloading the immutable package. The legacy
one-through-four-loop programs and the separate factor-product policy are not
configured by this setting.

`rustred_product_limits` separately bounds lower-factor sum-of-products work
and its coefficient cache (default 1 GiB and 16 million terms). It is not an
aggregate budget with the saved application's cache. Python exposes the
corresponding `rustred_product_cache_bytes` and `rustred_product_cache_terms`
options, with the same byte-only override convention. Changing either typed
limits setting replaces both owned sessions; no old cache bypasses new limits.
The existing four-loop native program also retains its own operational policy;
raising the outer product allowance does not relax that inner session's cache
or per-root limits.

Select `EvaluationMethod::RustRed` with `substitute_masters: false`. Numerical
five-loop master substitution is unsupported and returns an error. There is
no FORM fallback. In Python the constructor takes `rustred_artifact_path`,
and `load_rustred_application(path)` explicitly preloads the same package.

The standard expression contains exact coefficients of
`RustRedMaster(rustred_family_HASH, powers...)`. HASH is SHA-256 of the full
authenticated semantic family fingerprint. These masters denote unit-mass
Minkowski vacuum integrals with unnormalized loop measures. The expression
includes the common-mass factor and the selected Vakint loop normalization.
No epsilon expansion or numerical master values are inserted. Coefficients
hold in generic dimension `d = 4 - 2 epsilon`, subject to the saved nonzero
conditions. Unexpected parameter-dependent saved-program conditions are rejected;
conditions inherited from caller-owned scalar spectators and the common mass
remain explicit in the detailed symbolic result.
These are conservative execution-provenance conditions, not a minimal pole
set of the final rational expression. Poles inherited before a coefficient
cancellation remain recorded; different valid reduction paths need not have
identical condition lists.

RustRed's exact vacuum-product service prepares proven separations once after
loading. Existing one-through-four-loop native programs reduce each lower
factor under its original guards. FeynKit verifies the full momentum map to
the selected program, including a unit Jacobian. RustRed owns coefficient
products, coalescing, recursive separation and qualified-factor memoization;
the adapter does not implement an additional reducer. Numerator-aware charts
retain the complete verified momentum map, including inactive denominator rows.
Local numerators pass through the same native transport into the selected lower
program. Cross-component scalar products use the existing FeynKit tensor engine,
projecting one independent component at a time with other components as tensor
spectators. This uses native `d`; epsilon substitution happens only in the final
Vakint expression. The adapter retains the denominators of FeynKit's public
tensor-term coefficients before summing those terms. It does not claim access to
undisclosed intermediate pivot histories inside that engine.

The result contains sorted products of qualified declared masters, not a claim
of independent master products. Connected or unproved decompositions retain
their original qualified master labels. The common-mass factor and global
five-loop normalization are applied once, never once per factor.

`evaluate_rustred_symbolic` exposes the same value together with the exact
nonzero conditions and the map from master family labels to full semantic
fingerprints. Its Python counterpart returns `(value, nonzero_conditions,
families)`, with label/fingerprint pairs in `families`.
`rustred_application_cache_statistics` and
`clear_rustred_application_cache` inspect and clear the owned cache without
changing the accepted entry domain. `rustred_application_product_statistics`
reports the one-time factor preparation and explicitly counts retained terms.
`rustred_application_factor_cache_statistics` separately reports the lower
factor cache. Clearing application caches clears both mutable sessions while
retaining the immutable structural product plan.
`rustred_application_product_chart_statistics` separately counts native verified
numerator charts; these are geometry witnesses, not completed tensor reductions.
Projection occurs when such a terminal is used, so no additional unbounded
adapter coefficient cache is created. Owner decoding and validation are reported
by `rustred_application_owner_preparation_statistics`.

Catalogue preparation, package import, owner preparation and reduction are
separate stages. Owner programs are prepared lazily on first use; rule results
and lower-factor products have separate bounded mutable caches. Input matching,
verified witness selection and numerator lowering still run for each request;
a warm saved-rule cache does not imply that the whole frontend is cached.
The family's native inverse scalar-product basis is reused during lowering.
`rustred_application_cache_work_statistics` reports evicted entries, charged
coefficient terms/bytes, and recomputed nodes. Python includes these fields in
`rustred_application_cache_statistics()`. Work counters are cumulative, include
failed requests, and survive cache clearing. Recomputed nodes count repeated
cache misses within a single request; no unbounded cross-request history is
maintained. Cache byte allowances describe retained coefficient payloads, not
total process memory or an RSS limit.

RustRed combines scalar coefficients, spectators and common-mass factors before
requesting admission for lowered keys. Thus an expanded denominator-cancelling
numerator is admitted on its surviving endpoints; intermediate terms which
cancel exactly do not need independent coverage. Every remaining nonzero key
follows the package's authenticated route to its declared owner and must belong
to the recorded entry union. A numerator's total momentum degree does not
establish this. A squared momentum can leave both a denominator term and a mass
term, which may occupy different entry domains.
Missing rules, unsupported entries and unresolved terminals remain errors.
There is no adapter-wide cutoff at a particular numerator degree or total
denominator power. Native scalar-algebra, tensor-projector and verified-matching
work limits are operational bounds, separate from the authenticated entry union.
Cancellation is performed within each supplied numerator before admission;
separately supplied integral roots are admitted individually, not accepted on
the assumption that different unsupported roots will cancel later.
The five-loop handoff retains products and powers after restoring closed scalar
contractions. RustRed changes the complete polynomial to denominator coordinates
before expansion and native coefficient collection. Legacy unweighted lower-loop
handoffs retain their original expanded-input behavior.

== Native diagram ingress

`VakintExpression::from_diagram` consumes FeynKit's routed vacuum diagram and
family. `DiagramIntegralOptions::physical_slots` lists zero-based family
indices in sorted graph-edge order; the default is the physical prefix.
All other slots must have nonpositive powers. Zero physical powers contract
edges. Negative physical powers also multiply their original denominators
into the numerator. These operations preserve the original stored routing
simultaneously. A contraction removing an independent loop is rejected because
the remaining graph cannot represent that scaleless integration.
Authenticated zero-flow massive propagators are extracted as exact native
scalar powers before the graph contraction, preserving the five-loop measure.
A positive-power massless zero-flow propagator is rejected. Symbolic scalar
poles remain in the detailed result's nonzero conditions. Exponents outside the
native integer-power representation range are rejected before materialization.
Large supported numeric powers can still require substantial coefficient
materialization; that representation check is not a memory-cost guarantee.

The native Python function `integral_from_diagram` accepts the same
`physical_slots` and `powers` options. Compact DOT and exported DOT are parsed
by FeynKit; Vakint adds no DOT parser. The `five_loop_scalar_diagrams` example
uses `Process::generate_diagrams` and the scalar model with self loops retained
to create exported DOT fixtures. Its graph count is not an integral-class count.

The source routing data is transcribed from RustRed's
`tide_five_loop_parent_vertices.json`, `tide_five_loop_manifest.json`, and
`tide_five_loop_routing_witnesses.json`, based on Thomas Luthe's 2015 thesis,
tables A.1 and B.2. Unit-mass denominators are `q^2 - 1`; no numerical equality
to the thesis's Euclidean convention is implied.

The native generated-graph check maps all 166 connected phi-cubed/phi-quartic
vacuum graphs to 59 registered integral classes. The 67-class registry also
covers eight classes absent as direct representatives of that scalar generation
problem; the separate native DOT pinch round trips cover all 67. Thus neither
166 nor 59 is the number of masters or a replacement for the TIDE census.
The separate contraction-closure check enumerates 7,548 source contraction
representatives, deduplicates 101 native incidence classes, and verifies all
67 integral classes, including all 19 factorized classes, with no unmatched
class. Generation permits tadpoles and forbids zero-flow bridges; only
loop-preserving pinches enter this census. The massive-bridge ingress check
separately covers connected products excluded by the generation filter.

== Reproducible entry points

Build and run the native examples with the `feynkit-ingress` feature:

```sh
cargo run --release -p vakint --features feynkit-ingress \
  --example five_loop_scalar_diagrams -- /path/to/new/dot-directory
cargo run --release -p vakint --features feynkit-ingress \
  --example five_loop_saved_application -- /path/to/application --admitted-suite
```

The admitted suite uses three generation-54 witnesses: scalar point 162,
inverse-denominator numerator point 60, and scalar rank-eleven point 0.
It checks exact warm reuse, expanded numerator equality, simultaneous unimodular
basis changes with vertex/edge relabeling and orientation reversal, nonunit mass,
independent one-loop factor reductions, cache clearing, and settings opt-out.
Six additional positive-line controls compare an expanded cancelling numerator
against a separately constructed and parsed pinched DOT: single, doubled,
multiple and mixed dotted cancellations, including nonunit mass. Native
all-row momentum maps and weighted lowering first authenticate their common
admitted key without loading a package. Loaded comparisons check both the
qualified value and nonconstant guard obligations. Massive-bridge comparisons
check unit, nonunit and symbolic masses, retaining the symbolic mass pole.
The point-60 numerator `(q^2-1)^5` is inhomogeneous, with maximum momentum degree
ten and scalar degree five; point 0 has maximum momentum degree twenty-two.
Neither is a claim about arbitrary homogeneous polynomials. These are individual
admitted witnesses, not rectangular rank/power coverage claims.

The same example accepts a generated DOT path in place of `--admitted-suite`.
That mode builds the homogeneous degree-ten numerator `(D1-D2)^5` (GeV^10)
with two dotted propagators, and succeeds
only if all surviving endpoints lie in the selected package's actual union.
Failure of admission is a reported unsupported input, not a successful reduction.
For the pinned generation-54 package, the full parent-32745 two-dot control is
outside the recorded union: one surviving key is
`[-3,1,1,1,2,1,1,1,1,1,0,1,0,0,1]`. The separate homogeneous degree-ten control
on the admitted five-tadpole support succeeds. These are deliberately distinct
tests; the easier factorized case is not evidence for full-parent coverage.

For a scalar input with explicitly selected powers, use the separate generic
DOT mode instead of that degree-ten regression fixture:

```sh
cargo run --release -p vakint --features feynkit-ingress \
  --example five_loop_saved_application -- /path/to/application \
  --dot /path/to/scalar-vacuum.dot --powers 1,1,1,1,1,4
```

This example mode reads either supported FeynKit DOT dialect with the built-in
scalar model, takes numerator1, and normalizes its sole scalar mass to1 through
native simultaneous substitutions. The integer power list is in native family
denominator order; its length must match exactly. There is no topology-name
dispatch, new DOT parser or implicit rule generation. Existing named example
modes retain their behavior. Native routing and the loaded saved-domain union,
not an edge count or drawing symmetry, decide whether the input is supported.

The new ignored `rustred_five_loop_dot_tests` regression constructs a six-line
banana from native effective vertices, exports and reparses its stable DOT,
and checks powers `(1,1,1,1,1,4)` against the known admitted R0/D9 common-family
key. It requires an actual simultaneous routing witness and successful saved
admission before comparing exact qualified-master values and guards, including
a warm repeat with an invalid FORM path. Fast companion tests cover option
parsing, native routing preservation and power-count mismatches.

On 15 October 2026, all three fast tests and the actual ignored external test
passed against the separate completed generation1 R0/D10 package (329 terminal
maps), using RustRed `78cbc5c1`, Symbolica `400a53c4` and FeynKit `024f96d6` optimized
libraries. The external test took 231.84s, including 75.014s package preparation.
Whole-frontend first and warm applications took 51.866s and 50.630s, respectively:
these figures demonstrate correctness, not an end-to-end warm-cache speedup.
The terminal application cache recorded 47,804 cumulative hits and 3,753 rules,
retaining 7,627 integrals and 7.807MB of coefficients. The test did not separately
time each frontend phase, so a routing bottleneck remains a hypothesis rather
than a measured attribution in that test. Exact DOT/reference values, qualified
families and normalized guard sets agree. No numerical masters or FORM were used.

A separate run of the public `--dot` example then passed its standard/detailed,
simultaneous-basis/relabel/orientation, indexed-input, guard-domain and cache
option controls in 353.045s supervised wall, with 2.884GB sampled peak RSS. Its
existing phase logs resolve the bottleneck: the 44.808s warm total includes
44.789s of admission-aware witness search but only 0.000572s of saved-rule
application and 0.001525s of terminal handling. First application was 53.832s,
including 41.045s witness search and 12.708s saved-rule application. Package
preparation was 73.366s, reported separately. Thus intermediate coefficient
memoization works; repeated frontend witness selection dominates this DOT
control's warm time. These are separate runs, not a controlled speed comparison.
This one R0/D9 connected control does not cover the full R6 scope, all graph
classes, the complete D10 boundary or the physical numerator ladder.

The subsequent prepared-witness-memo gate passed five focused tests, an external
exact minimum-versus-uncached test (including clone/registry/policy mutation,
stale-key revalidation and optional-mutex poisoning), the existing 4MiB-LRU
regression, and this same public DOT example. Its observed whole-frontend warm
time was 0.032923s instead of 44.808307s; hit witness selection took 0.003011s and
cached IBP application 0.001266s. Four witness hits and two misses retained two
entries, with 4,767 tracked Atom bytes. Exact results, guards and all example
invariance checks passed; the IBP work remained 3,753 rules.

This is a repeated-input benefit, not a first-use or general throughput claim.
The first application was 64.196s versus 53.832s before, and an uncached rerouted
witness search was 109.639s versus 46.879s, despite identical printed inputs.
The source of that cold-time variation is not established. The whole public
example took 250.845s versus 353.045s on the same CPU set, in single shared-host
runs; these are observations rather than a confidence-bounded speedup estimate.
The native package and its coverage were unchanged.

A subsequent frozen new-then-old repeat passed all comparisons again. First
application was 78.078s new versus 78.374s old; warm application was 0.022868s
versus 76.409812s. The uncached rerouted search was 83.088s versus 113.597s, so a
consistent new-only cold regression was not reproduced. Every full search still
visited 721 graph mappings, 23,040 orientations and 1,440 verified witnesses; memo
reuse reduced the number of full searches from six to two. The complete example
took 238.774s new versus 625.756s old, with byte-identical canonical final output.
Shared-host load and cold-time variability remain qualifications; neither this
control nor its warm ratio establishes general five-loop throughput or closure.

The Python `examples/five_loop_saved_application.py` uses the installed native
HepKit DOT reader, family constructor and Vakint ingress. Its optional `--powers`
argument follows the native family order. It checks an exact warm repeat of the
value, guards and family identities, then prints a bounded preview and counts.
Use `--verbose` to print the complete guarded result and family identities.

The ignored external-package tests require `VAKINT_SAVED_APPLICATION`; the
813-witness admission test additionally requires
`VAKINT_SAVED_APPLICATION_INVENTORY`. That inventory is a supplied finite
witness set, not an exhaustive enumeration of every point in the domain union.
The admission test reconstructs each common-family integrand, proves its
weighted-lowering identity, and runs ordinary canonicalization, tensor
preprocessing and verified witness selection. It asserts that no owner program
was prepared and no IBP rule was applied. The bounded-cache test separately
requires actual eviction and compares exact repeated values, guards and family
identities while checking the retained byte allowance.

The 813-witness diagnostic can explicitly select
`VAKINT_SAVED_APPLICATION_TRANSPORT_PROFILE=materialize-512m`. This test-only
profile is shared with the two-case transport diagnostic. It allows 512 million
native and aggregate transport operations; 128 million exponent and endpoint
power entries; eight million operational nodes, projected endpoints and native
polynomial terms; and two billion bytes each for retained endpoint keys and
coefficient clones. Native exact-coefficient, scalar-lowering, matching and
reduction policies remain unchanged. The complete selected options are printed
before loading. The profile cannot be combined with the older single
`VAKINT_SAVED_APPLICATION_NATIVE_POLYNOMIAL_OPERATIONS` override.

The default-policy census retained 811 passes and two resource failures; raising
only the native-operation allowance did not address the independent aggregate
operation and exponent-entry caps. The complete profile admitted both remaining
witnesses through the unchanged production matcher without owner preparation
or IBP application, peaking at 3.72 GB in that two-case run. This is an explicit
operational choice, not a larger mathematical entry union or a production
default.

On 13 October 2026, the unchanged full census then passed all 813 witnesses
under this same complete profile, with one preparation worker. The successful
IDs exactly matched the inventory, all reconstruction identities passed, and
owner preparation, rule-cache entries and IBP applications remained zero.
The admission loop took 971.550 seconds; the whole test including loading took
1046.92 seconds, with sampled peak RSS 3.847 GB on shared-host CPUs 8–15.
The frozen actual-Git diagnostic used RustRed `1fca5bd0` and executable SHA256
`eaa4dbf2f0e85caa26c4c6f863266cf3c9d425714c5dd25e3edeb65e83e45cdb`.
These are finite admission witnesses, not 813 completed reductions or an
exhaustive proof over the entry union. The earlier default-policy failures are
not erased or reclassified as successful default-limit runs.

A separate ordering hypothesis was tested and rejected: the problematic
candidate already had just one coalesced rank-ten key and no rank-zero key to
check first. The manual diagnostic retains that negative result; it is not an
acceptance test. No production predicate ordering or fast-rejection API was
introduced to obtain the successful full-census result.

== Application validation and remaining performance gap

The post-weighted-lowering, actual-Git `1fca5bd0` application suite passed all
three primary inputs, expanded and simultaneous-basis-change comparisons,
all six positive-line cancellation controls, all three massive-bridge controls,
the homogeneous degree-ten control, and every guard-domain comparison.
The frozen executable SHA256 was
`dd5307d2c1f9a72ede886450f310a7a261e88a4e677f0c4b95683c7f4a39afdb`.
It completed in 3203.631 seconds (53 minutes 24 seconds), with sampled peak RSS
2.884 GB. The current frozen library also passed the 52-terminal factor-product
gate in 202.90 seconds, plus six numerator, four matcher and one five-loop
basis-invariance tests.

This successful application run was materially slower than the earlier
`bf5b2882` run (1365.807 seconds); the runs were not a controlled benchmark
pair. Four frontend admitted-witness searches account for almost all of the
observed increase: their combined time rose from 896.439 to 2718.036 seconds.
Each search retained exactly 73 graph maps, 2304 orientation candidates and
144 verified witnesses in both runs. The evidence therefore points to higher
per-candidate frontend cost, not an expanded witness search space. Each
candidate repeats scalar sealing, parent transport, weighted lowering and
entry admission; existing timers do not separate those inner costs. A live
profile was requested but the process had already finished, so no finer
attribution is claimed and no run was restarted. This is an open performance
issue, separate from correctness and from saved-rule cache effectiveness.
