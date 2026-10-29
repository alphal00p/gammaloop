= First-principles mark geometry

== Source and review boundary

This separate WP11 drawing lane starts at approved source36 integration
`100da20603a083d031c31bd4e8c6970fc44dda2f`. Its prerequisite review approved
asset diagnostics and cleared the drawing gate, not all WP3 or final packaging.
The original drawing bookmark `generation-drawing-ports` remains at `98cbf5b8`.
The new bookmark is `generation-mark-geometry`; main integration is untouched.
Source38 `15dd7e1b188dad112d7887bdc6d769125b2697a1` is reconciled here; no
source39–46 change or historical binary is imported. This implementation awaits
NEW independent Astra review. WP11/13 final artifact gates remain open.

The #link("evidence/generation-integration/mark-geometry/receipt.json")[receipt]
records exact source, tool, package, notice and output identities, changed-file
and asset inventories, fresh versus reused validation and reproducibility limits.
The #link("evidence/generation-integration/mark-geometry/audit.json")[audit]
records provenance and notice-owner decisions.

== Geometry and public contract

Kurvst's existing MarkSpec owner constructs thirteen leaf shapes or a composite
in a local tip/shaft frame. Public shape names, accepted fields, Typst data and
Python options remain unchanged. Composite children still omit fit/shorten;
finite checks, selected versus candidate batching, carrier clamping, native SVG,
MathML, static notebook output and optional-widget contracts remain owned by
their existing implementations. This drawing geometry intentionally uses f64;
no physics, numerator factorization or general generation code changes.

A polygon's vertices move along intersections of its offset edges. Its sharpest
corner determines its miter limit, with numerical headroom. A head that cannot
hold its stroke grows before an inset edge collapses or reverses, including deep
stealth notches. Paint keeps the nominal box when that box can hold the stroke.
The drawn rear intersection stops the shaft and places the next composite part.
The straight chevron uses tangents from the tip to its round caps; TikZ uses
quarter-circle arms; round darts use the inner box. Chord-fitted heads rotate
with the changed contact. The default diamond is exactly the default square
rotated 45 degrees. Default ray spoke sets retain their previous visual phase
except counts divisible by four from eight onward; explicit phases are retained.
Leaf/composite enums remove unreachable string-dispatch fallbacks.

Tests exercise mathematical bounds and contacts, sharp miters, scaling and
stroke-dominated growth, deep notches, reversal, gaps/composites, ray symmetry,
diamond rotation, curve fitting, painted coverage and packed/native transport.
Intended visible changes are those geometry/contact rules, not new label or
physics conventions. Updated expectations express the changed geometric formula;
no tolerance or golden snapshot is loosened.

== Round 1 geometry corrections

Scope B approved the original source/notice/assets audit at `cc7bc03a`; scope A
required corrections for numerical polygon inset failure, full-turn barb contacts
and one test import. The old candidate, binaries, WASM and reviews remain retained
under `generation-mark-geometry-round1` and the original review directories.
The #link("evidence/generation-integration/mark-geometry/round1-correction/receipt.json")[correction receipt]
records the old-fail/new-pass native and public Typst/WASM producer identities.

Polygons now retain normalized coordinates until growth and inset are resolved.
Power-of-two vector normalization preserves component bits near opposite edges.
A compensated determinant and stable half-angle identities compute corner slides
without the former angular cutoff. Each edge's linear contraction supplies its
analytic collapse limit, replacing an absolute-resolution bisection. Uniform
growth is computed in normalized units and then converted back to drawing units.
Tiny positive dimensions, including subnormal uniform sizes and the closest
representable inset below 100%, therefore cannot
masquerade as zero-area polygons. The as-given exception is restricted to parameter-derived zero area:
zero dimensions or the doubled-back 100% notch. Numerical normalization cannot
turn a positive-area input into this exception. Nonrepresentable numerical results propagate an error rather than a
successful uninset polygon; the reviewed valid positive-area cases stay accepted.
No clamping, substitute shape, unbounded growth loop or narrowed public option
range is used.

Barbs enumerate axis crossings in their signed angular interval. A complete
positive or negative turn stops at the circle's rear; a neighboring open turn
retains its open gap and front crossing. Reversal transforms the computed front
and rear through the existing local frame; composites inherit the actual rear.
New native and public Typst fixtures cover the two old numerical mechanisms,
true degeneracy, signed full/neighboring turns, curve and composite contacts.
The PathSeg test type is imported at test-module scope. Notice/package ownership
and the scope-B legal/source-audit boundaries are unchanged. The corrected code
and changed artifact bindings require NEW independent Astra reviews.

== Source and notice audit

The pre-port Rust mark owner equals source38's parent. It explicitly credits
Tiptoe 0.4.0 and reproduces its per-shape geometry. The pinned Nix Tiptoe source
and MIT permission text were inspected, including its marks, line, stroke/dimension
utility and curve fitting owners. The replacement comes from source38's local
polygon-offset, tangent and Kurbo-arc constructions, not a new transcription of
Tiptoe or LGPL CeTZ code. It removes the old trigonometric dart formula, fixed
miter limit, TikZ magic control points, manual arc builder and old branch fallback.

Retained shape names, public field vocabulary, dimensional defaults and simple
bar/arc/rectangle meanings preserve the earlier API. Those compatibility data
and elementary constructions are documented provenance; no Tiptoe package source,
renderer, resource or material per-shape implementation is retained in the new
runtime geometry. Existing Rust serialization, batching and carrier/path owners
are GammaLoop/Kurvst implementation, not the Tiptoe renderer. Historical notes,
archived evidence, old revisions and the read-only audit source keep their history.

Consequently the obsolete Kurvst and notebook-host Tiptoe license copies are
removed with their embedding assertions and Nix input references. Standalone
Linnet packaging is absent and is not revived. Common root-license selectors
are unchanged. Clarabel and NetworkX notices remain byte-exact and delivered.
Both sync recipes remove only the obsolete generated Tiptoe notice, copy every
remaining package notice and preserve unrelated files. Upgrade tests require
idempotence, restored missing root licenses and preserved user notices; sdist,
API extraction, Nix payload and fixture-wheel checks reject obsolete attribution.

== Validation and remaining gates

The original round 1 results below are source-bound to `cc7bc03a`; corrected
producer/check results and refreshed artifacts are recorded in the correction receipt.

Original focused results: 344 Kurvst/Linnest tests with native SVG enabled pass,
with only the fixture-rewriting test skipped; all 46 Clinnet tests, five packaging
regressions and three actual Typst/WASM layout tests pass. Scoped Clippy with
denied warnings, workspace formatting and Python formatting pass. The unchanged
canonical registry oracle and independent Git hashing validate 104 scopes,
66 titles and 16 negative controls; semantic manual table headers are preserved.
Five affected notes render to PDF and HTML using the canonical document wrapper.

The bundle contains 29 files matching the checkout, built from 1402 matching
source files. The independent final Nix producer rebuild compares byte-identical
for both WASM modules and its whole package. It normally reuses source-qualified
dependency artifacts: no cold dependency rebuild or general AHash determinism
claim is made. All 34 website SVGs reproduce; 24 intentionally change and ten
remain unchanged. The manual and three example PDFs reproduce at timestamp zero.
The two-theme mark gallery reproduces in SVG, PDF and PNG, and was visually
inspected with the affected two-theme diagram. The new note's PDF/HTML also
reproduces. Earlier bare-note bundle exports required the existing document
wrapper; that tooling correction required no authored-source fallback.

Changed WASM and the manual/example/website outputs use the canonical current
Nix and timestamp-zero Typst owners. Visual checks cover representative marks,
curves and stroked composites plus affected diagrams in both themes. Native
checks compile current Rust owners; an old notebook-host extension is never used
as evidence for changed marks. Public API/stub metadata is unchanged.

Source-qualified build provenance and independent byte reproduction are distinct.
Nix output/derivation and installed asset identities are recorded. Any remaining
Linnest data/name-section byte mismatch is explicitly qualified even when its
executable sections agree. Cached reuse alone is not reproducibility evidence.
The canonical site qualification remains frozen at `a372`; this lane supplies
bounded evidence for its changed documentation and figures.

Fresh independent Astra review is required before parent integration. Complete
repaired/compressed native/Pyodide and externally locked Community wheels still
need METADATA/RECORD/notices and installed-behavior inspection. Fixture wheels
are packaging tests, not those release artifacts. No solver replacement, callback
fix, frame/weight diagnostics, full-site or physics qualification is claimed.

== Round 2 painted-miter correction

Scope B approved the round 2 source/notice/artifact bindings at `e426ff5f`;
that candidate is retained under `generation-mark-geometry-round2`. Scope A
found that finite thin-triangle tips and deep-stealth wings were still dropped
by the downstream Kurbo stroke evaluator. The frozen native and public WASM
producers reproduce both failures before rebuilding. The prior barb/contact,
import and inset corrections remain in place.

The authoritative Kurbo 0.13.1 stroke owner is now a local patch, selected by
both the root and notebook-host manifests. Within the accepted correctness task,
the root coordinator authorized the narrow exact Kurbo 0.13.1 vendor patch and
necessary root/notebook-host patch-selection entries; all unrelated dependency
versions, pins and locks stay fixed. The existing Nix local-patch discovery includes this owner in filtered
build contexts. Its published Apache-2.0 and MIT license texts are retained;
the existing Clarabel/NetworkX notices and common selector are unchanged.

The stroke owner rescales tangents by powers of two, compensates cross products,
and rationalizes the miter ratio near a half-turn. Offset intersections are
computed in the vertex's local frame before translation. This repairs actual
paint coverage, rather than changing limits or bounds metadata. Real zero-turn
and miter-limit clipping retain their existing behavior. Public mark names,
options, wire format, stubs and the native/Typst painting contracts are unchanged.

Complete-bound and actual painted-interior tests use independent lobe inradii
across three scales, skinny triangles and deep notches. Their new analytic
expectations account explicitly for coordinate rounding amplified by reciprocal
notch depth; no existing tolerance is relaxed. A separate 100-digit oracle
intersects the exact emitted f64 offset lines to distinguish construction
conditioning from stroke-evaluation error. Prior assertions and true-degeneracy
controls remain active. The correction receipt records fresh checks, rebuilt
assets, exact native/WASM producers and remaining qualification limits.

See the #link("evidence/generation-integration/mark-geometry/round2-correction/receipt.json")[round 2 correction receipt].
Fresh independent geometry and artifact-binding reviews remain mandatory.
Final native/Pyodide/Community wheels and full-site qualification remain open;
normal dependency-cache reuse does not establish cold-build determinism.

The API's embed-only feature selection initially omitted Kurbo from Maturin's
source archive although Cargo still parses the root path patch. The existing
sdist include owner now explicitly ships the complete published Kurbo directory,
including both license texts and provenance. The packaging regression compares
every vendored file with the actual source archive. This source-closure repair
adds no wheel or external Community artifact qualification.

The pinned raster renderer still leaves the extreme skinny-triangle tip crop
blank; full extreme extents are below display resolution. Those images are
paired with numeric bounds and are not paint-coverage acceptance evidence.
Native winding/complete-bound controls and public WASM bounds, shaft and
composite controls provide the corresponding numeric checks. Ordinary galleries,
barb contacts, deep-notch tip crops and diagrams are inspected in both themes.
