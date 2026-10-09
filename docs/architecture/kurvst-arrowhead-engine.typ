= Kurvst Arrowhead Engine

#quote(block: true)[
*Status: Accepted design, recorded 2026-10-04; implementation in progress.*
P0 and P1 are verified. P2–P7 are implemented with focused validation; final
performance and image acceptance remain pending. This is not
a claim that the complete stock-CeTZ renderer or performance targets are ready.
]

== Decision and delivery boundary

Move numeric arrowhead geometry, label-candidate footprints, hit-area placement,
and collision sampling into Kurvst. One Rust geometry engine serves native
Python SVG and Typst Wasm; stock CeTZ paints finished geometry. Kurvst owns
curve and mark geometry; Linnest owns placement decisions, label search, and
link semantics; CeTZ owns painting and its existing Boolean operations.
Eventually delete the vendored CeTZ tree: 47 files, five core patches, and its
pinned rebuild recipe. Do not port LGPL CeTZ patches into the new MIT engine.

The prerequisite consolidation was finalized by the parent, with 58 native
tests reported green, as `uyrstsvl` /
`da456712676c16db00da377b2447008dc28c73e7`. This separate documentation change
is `rnzmqrrv`. These are provenance records, not checks reproduced for this note.
Delivery order is consolidation, private reconciliation with feynkit, then
regeneration of P0 against that base. Reconciliation is recorded as `lxptuzuk` /
`cf565f42b9c36ee20ac171a84f80c8347c889de4`; the shared feynkit branch was not
rewritten. The phase records below distinguish implementation evidence from
remaining end-to-end acceptance.

== Mark catalogue and data contract

`kurvst.mark` constructors produce data only. Python `linnet.Mark(name, ...)`
mirrors the same names, fields, and defaults and serializes the same dictionary;
there is no Python `kurvst` module. Names, parameters, and defaults follow
#link("https://typst.app/universe/package/tiptoe/0.4.0/")[tiptoe 0.4.0 (MIT)],
with attribution retained in the fresh implementation.

The catalogue has thirteen shapes plus `combine`. Percent lengths below are
ratios of line stroke thickness, not path length. Mixed dimensions add fixed
length and stroke-relative length.

#table(
  columns: (auto, 1fr, 1fr),
  [*Mark*], [*Accepted size / shape defaults*], [*Placement / prior use*],
  [`triangle`], [Stealth geometry with inset 0%], [Native flow today],
  [`straight`], [Length 3pt + 450%], [Native momentum today],
  [`stealth`], [Length 3pt + 450%, inset 40%], [CeTZ today],
  [`round`], [Stealth with round joins], [New],
  [`tikz`], [Width 3pt + 450%], [New],
  [`barb`], [Width 3pt + 450%, arc 180deg], [New],
  [`hooks`], [Width 3pt + 450%, arc 180deg], [New],
  [`bar`], [Width 2.4pt + 360%], [Centred],
  [`bracket`], [Width 2.4pt + 360%], [End aligned],
  [`circle`], [Length 400%], [Centred],
  [`square`], [Length 400%], [Centred],
  [`diamond`], [Length 565.69%], [Centred],
  [`rays`], [n = 4, length 280%], [Centred],
  [`combine`], [Marks and gaps in sequence], [One composite],
)

Parameters include `length`, `width`, `inset`, `fill`, `stroke`, `rev`,
`align`, `arc`, `n`, and `phase`. Automatic width follows the per-shape ratios.
Lengths may be fixed, stroke-width ratios, or mixed, such as `3pt + 450%`.
Fill defaults to the line stroke paint; sizing uses the line stroke thickness.
The tiptoe end distance defines the shaft stop, including the stealth inset
notch. Arrow tips land on the path end. Bar, circle, square, diamond, and rays
are centred unless `align: "end"`.

`combine` is one composite with aggregate contacts, end distance, and footprint.
Apply the end fit once to the composite, not independently to its children.

== Carrier, placement, and curved-end fits

Preserve `mark-position`, `mark-orientation`, `mark-direction`, and `mark-shift`,
including centre placement, arc-length ratios, signed shifts, inward endpoint
clamping, and label footprints. Stations and shifts use the derived visible
carrier *after offsets and node outsets, before arrow shortening or bending*.
Footprints and collision geometry instead consume the painted shaft.
Kurvst arc length is the chosen source of truth: accept its tiny numerical drift
while keeping Python and Typst geometry identical.

#table(
  columns: (auto, 1fr),
  [*Fit*], [*Accepted rule*],
  [`chord` (default)],
  [Put tip and back on the original carrier; trim the shaft at the back without
   changing the carrier. This preserves the current arrow behavior.],
  [`bend`],
  [Follow the end tangent and retract the endpoint by mark end distance times
   `shorten`, retaining the control points as tiptoe does. Default `shorten`
   is 100%; 0% is unbent.],
)

Fits apply at ends. Interior marks straddle their station on a chord.
The painted shaft owns footprint and collision geometry; path patterns follow
the bent base when `bend` is selected, including phase behavior.

The curved-end comparison has no image asset in this workspace. Its intended
reference uses the same cubic and a 22 × 16pt triangle: `bend` with `shorten: 0%`
leaves the line exiting the side; `bend` at 100% reshapes the endpoint to the
head centre; `chord` turns the head to the original carrier while leaving the
curve unchanged. Large heads exaggerate the differences. P0 will generate real
references; no missing `curved-ends.png` is linked here.

== Execution boundary: two bounded batches

Each drawing uses *two bounded batches*, both sharing templates and the engine:

+ Prepare candidate geometry and footprints for Linnest's label search.
+ Paint the selected placements using finished geometry.

This is not a one-call design. Candidate preparation avoids repeated Typst
geometry evaluation; numeric hitboxes and prepared content drawables avoid
per-target painting work. Link target selection and semantics remain Linnest
responsibilities, separate from the mark catalogue.

== API and migration

Mark constructors and the native Python configuration are implemented. The
Typst renderer migration below is still in progress.

These sketches describe the accepted interface direction. They do not establish
currently available constructors or native configuration options.

```typst
#import "../../crates/kurvst/typst/lib.typ": kurvst
#let flow = kurvst.mark.triangle(length: 0.21cm, width: 0.1575cm)
#let momentum = kurvst.mark.straight(length: 0.16cm, width: 0.12cm)
#let cut = kurvst.mark.combine(kurvst.mark.bar(), 2pt, kurvst.mark.bar())
#let bent = kurvst.mark.stealth(fit: "bend", shorten: 100%)
#let source-style = (
  stroke: ink + 0.7pt,
  mark: flow,
  mark-position: "center-if-dangling",
  mark-orientation: "edge",
)
// source-style is supplied to the drawing style configuration.
```

```python
import linnet

flow = linnet.Mark("triangle", length="0.21cm", width="0.1575cm")
momentum = linnet.Mark("stealth", inset="25%", fit="bend")
config = {
    "style": {
        "edge-style": {
            "flow-arrow": flow,
            "momentum-arrow": momentum,
        }
    }
}
diagram.render(config=config, momenta=True)
```

Native configuration gains `flow-arrow` and `momentum-arrow`. Reject CeTZ mark
dictionaries with an error showing the replacement; remove `linnet.Mark` CeTZ
fields `scale`, `anchor`, and `shorten_to`, with useful migration errors.
Switch physics `fermion-arrow-mark` and `momentum-arrow-mark` together with
documentation and examples. Preserve the notebook look using explicit existing
sizes: tiptoe defaults are larger (roughly 7.5pt for straight at a 1pt line,
versus roughly 4.5pt today).

On the measured feynkit revision, native marks are hardcoded: momentum chevron
0.16 × 0.12cm and flow triangle 0.21 × 0.1575cm. There are no native arrow options
there. Current Typst vendored-only calls are `mark.geometry`, `mark.footprints`,
and `draw.content-many`; their replacements are planned, not present here.

== Reported measurements and rationale

These figures were supplied with the accepted plan and were *not reproduced*
for this documentation change. Context: 24 drawing-corpus diagrams, Typst
0.15.0, feynkit `78c33d71`, measured 2026-10-04; medians of three interleaved runs.

#table(
  columns: (1fr, auto, auto, auto),
  [*Renderer*], [*Total*], [*Momentum*], [*Plain*],
  [Vendored CeTZ], [26.7s], [~1.3s], [~1.0s],
  [Vendored, three Linnest batched calls replaced], [42.3s], [~2.4s], [~1.2s],
  [Stock CeTZ 0.5.1], [70.0s], [~5.5s], [~1.5s],
)

Stock rasterized all 24 diagrams to the same pixels in that reported comparison;
SVG differences were roughly 300 invisible hit rectangles. This observation is
not a blanket pixel-identity requirement for the migration.

Of stock's extra ~4.75s momentum time, candidate footprints account for 79%,
hit areas 16%, collision lines 2%, and generic drawing 3%. Of the extra ~0.6s
plain time, hit areas account for ~75%, collision lines ~14%, and generic
drawing ~9%. Stock momentum used 2.7 million Typst calls versus 186 thousand;
the Wasm interpreter was reported to cost roughly 20 times native execution.
These observations motivate numeric preparation rather than merely replacing
the painter. P0 must record host, toolchain, cache state, timing boundaries, and
paired baseline comparisons before using these figures as performance gates.

== Delivery phases and evidence

Each phase is its own change after the reconciled feynkit base, with corpus
timings recorded. S/M/L are planning sizes, not completion indicators.

=== P0 — S: reproducible corpus and visual policy

Check in the 24-diagram corpus, regenerated templates, median timing script,
and PNG references. Establish explicit lenient tolerances and image-diff review
for arrows and labels, *not blanket pixel identity*. Require cross-renderer
outline agreement and correct link targets and hit areas. Record provenance
(host, toolchain, cache, timing boundaries) and paired baseline comparison.
Regenerate this baseline only after consolidation and feynkit reconciliation.

The P0 harness, immutable input corpus, and reconciled reference images are
now recorded in `tests/drawing-bench/`. The captured three-round baseline is
15.010 s on M4 Pro / Typst 0.15.1, with 24/24 reference PNGs byte-identical
to the supplied snapshot at 150 ppi. Source, package, host, and raw timing
identities accompany the images. Numerical tolerances for future arrow/label
changes still require explicit review; renderer migration is not claimed.

=== P1 — L: shared geometry engine, no renderer migration

Implement fresh MIT Kurvst templates for all thirteen marks and combine,
sizing, end distances, placement, both fits, inward clamps, contacts, painted
shaft, and footprints. Add batched native/plugin APIs and data constructors.
No user-visible renderer migration yet. Compare straight shapes and end
distances with tiptoe; test bend/shorten including short last segments,
painted-footprint coverage, shaft/back contact, zero-length paths, and near-end
placement.

P1 now provides the shared native/Wasm engine and strict data-only Typst
constructors for thirteen shapes and composites. Its 110 native tests and
native/Wasm parity checks pass. The approved parity tolerance is
`1e-12 * max(1, abs(actual), abs(expected))` for geometry coordinates and
lengths only; structure, styles, indices, and other data remain exact.
Candidate results remain independent, while selected results share one
authoritative painted shaft per carrier.

=== P2 — M: native and Python

Replace hardcoded `linnest::svg::marks` with the engine, add native
`flow-arrow` / `momentum-arrow` and the Python mirror. Preserve the default
notebook look with explicit sizes.

The native renderer now uses one drawing-wide candidate call and one selected
call, with deduplicated templates and shared shafts. Fourteen focused SVG tests
pass. A fresh CPython 3.13 ABI3 wheel passes all 100 API/mark/streaming tests, including
actual Typst rendering of every shape and composites with both fits. Named paints
resolve through Typst's own color constants, not CSS names. FeynKit and Spynso3
caller checks and both mark-configuration integration tests pass. Native timing
and full end-to-end acceptance remain pending.

The maintainer approved the shared engine's default appearance on 2026-10-09:
retain the explicit physical sizes and inset stroke geometry rather than
recreating legacy protruding miters. The isolated raw-paint comparison at
150ppi changed 67/2814 pixels for the flow triangle and 74/2814 for momentum;
maximum bounding-extent changes were 0.4272pt and 1.6512pt respectively.
This approval does not waive corpus image review or native/Typst parity.

=== P3 — M: Typst painting

Replace `mark.geometry` with engine geometry painted through CeTZ; accept only
Kurvst marks. Migrate physics, manual, and map-style uses. Remove one of the
three vendored-only calls.

The Typst route now consumes the shared candidate and selected batches.
Physics templates, the map-style example and Python selectors use canonical
mark data; CeTZ edge-mark dictionaries are rejected. All five Clinnet public
behavior tests pass against the stock CeTZ 0.5.1 cache. Named preparation
anchors remain on the unshortened reference carrier; final painted heads
contribute their engine outlines to conservative canvas bounds.

=== P4 — M: candidate footprints

Replace `mark.footprints` and the slow paint fallback with shared candidate-batch
engine footprints. This addresses the main 79% momentum penalty.

Candidate templates and placement indices are collected drawing-wide before
label search. Outlines determine collision boxes; retained sizing thickness
cannot turn a shaft with no paint into an obstacle. This is implemented and
covered by focused tests, not yet a corpus performance claim.

=== P5 — S: hit areas

Use numeric hitboxes and prepared content drawables instead of per-target CeTZ
`content-many`. Address roughly 75% of the plain penalty and 16% of momentum's;
keep link semantics separate from the catalogue.

The authored fixed-size identity-target helper passes focused stock CeTZ 0.5.1
and 0.5.2 checks, preserving both expected SVG links. Integrated public
rendering tests also pass; corpus acceptance remains pending.

=== P6 — S: collision sampling

Land the collision-sampling prototype JJ change `tvxztqnq` with one plugin batch.
Its reported vendored-corpus improvement was 1.8%; collision accounts for
roughly 2–14% of the stock penalty. These figures are not new measurements.

The prototype revision was unavailable. The implementation instead reuses the
existing Linnest `label_stroke_lines` owner, with strict numeric validation and
one batched traversal of painted strokes. Twenty-three label-placement tests
and a focused stock-CeTZ fixture pass. The drawing owner now gathers all
painted strokes for one traversal; timing acceptance remains pending.

=== P7 — M: stock CeTZ and vendor removal

Target *stock CeTZ 0.5.2* (the measurements above used 0.5.1). Delete vendor,
patches, and `REBUILD`; migrate behavior fixtures and Nix packaging; rerun the
corpus. Generic bounds/sample speedups may optionally be upstreamed.

The patched 47-file package tree, five patches and rebuild recipe are deleted.
Nix selects stock CeTZ 0.5.2 directly. Offline Python rendering embeds one
unmodified upstream release archive, with its LGPL license and SHA256
`77cf8490114ae04c6e665a11efa691d284a0cadb9719771b5708c1197292f23f`;
there is no editable CeTZ fork or runtime download. Bundled package identities
replace complete staged packages, including conflicting file/directory entries,
without changing external stores. Stock rendering tests pass.

== Validation and remaining acceptance

The first integrated corpus completed all 144 samples but regressed to 48.784s
on stock CeTZ 0.5.1 against the same-host 15.010s P0 baseline. Profiles excluded
stock painting: current and stock package variants differed by only 0.38%.
Growing drawing-wide values crossed Typst function boundaries once per
candidate, creating quadratic preparation and projection overhead.
Registration, bounds projection and preview materialization are now bulk
operations; the shared engine supplies ordered conservative bounds directly.
Exact context/paint/geometry equality and unchanged PNGs cover these
optimizations. The problem-case stock time fell from approximately 2.99s to
0.982s before the final corpus rerun. This is improvement evidence, not a
claim that the same-host no-regression gate is met.

The first complete stock 0.5.2 corpus measured 18.127s, 20.76% slower than the
same-host P0 baseline. It passes the separate 26.7s total ceiling, but fails the
no-regression requirement and the 1.3s worst-case momentum ceiling. The 1.0s
plain ceiling also fails on the largest diagram; that diagram already took
1.077s at P0, so this is not a newly introduced plain-rendering slowdown.
All 24 staged-archive/Nix-stock PNG pairs match exactly. Eighteen differ from
P0, with unchanged dimensions and at most 1.52% changed pixels in the reviewed
high-change cases. Image review remains required. Later parser and packed-path
transport improvements and removal of unused painted-length metadata are now
covered by a second complete 144-sample corpus: *17.456s*, or *16.29% slower*
than P0. All 24 package-paired PNGs match, and all 24 outputs match the preceding
corpus exactly. The worst momentum case takes 1.583s; the worst plain case
takes 1.079s. Performance acceptance is therefore still open, despite passing
the separate 26.7s total ceiling.

`tests/drawing-bench/native.py` provides a real installed-Python, fixed
four-loop native-render check, without generating a diagram population or
including Typst label typesetting. After renewing the local license and
preserving sub-resolution positive trim connectors, the real installed run
passed: 10 warmups, 21 timed samples, median *12.186ms* against the 100ms
absolute ceiling. All 28 catalogue/composite and fit combinations also render
with their requested native paints. No genuine old same-host native timing
baseline exists, so native non-regression is not independently established.
Credentials are supplied only through the process environment.

The installed-extension test exposed distinct Python type identities across
the host and standalone Linnet binaries. `Mark.from_dict` / `to_native` now
provide an explicit portable data boundary through the existing codec, rather
than duck typing or exception-driven fallback. The actual separately loaded
extension regression passes.

== Acceptance and risks

On the recorded reference setup, with paired comparisons, the stock corpus
must reach at most 26.7s total, 1.3s momentum, and 1.0s plain. Native rendering
must not regress from the reported 60–100ms four-loop baseline. P0 must make
the reference setup reproducible rather than treating these numbers as portable.

Both renderers must use the same geometry for thirteen marks plus combine.
Bend must match tiptoe within the recorded tolerance; chord must match current
arrows under the P0 visual policy. Cross-renderer outlines, correct links and
hit areas, painted footprint coverage, and shaft contacts remain correctness
requirements. Completion also requires no vendored-only calls, vendor removal,
and passing stock-CeTZ tests.

- Interpreter fuel and instruction count: keep candidate loops free of
  trigonometry and powers; cache templates and arc-length tables.
- Label drift: review image differences with explicit lenient tolerances;
  tiny Kurvst arc-length drift is chosen, not hidden by separate renderer rules.
- Breaking API: migrate both arrow options, examples, and physics styles
  together; reject obsolete dictionaries and fields with actionable errors.
- Bent path patterns: validate phase against the painted bent base.
- Licensing: fresh MIT geometry with tiptoe attribution, not copied LGPL
  CeTZ patch code.

This note records acceptance of the design only. None of the P0–P7 acceptance
checks, migration outcomes, or performance targets is claimed complete.
