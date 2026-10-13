= Kurvst Arrowhead Engine

#quote(block: true)[
*Status: Accepted design, recorded 2026-10-04; P0–P7 implemented.*
Linnest draws with stock CeTZ 0.5.2 and the shared Kurvst mark engine; no CeTZ
is vendored. Performance and image acceptance remain open decisions.
]

== Decision and delivery boundary

Move numeric arrowhead geometry, label-candidate footprints, hit-area placement,
and collision sampling into Kurvst. One Rust geometry engine serves native
Python SVG and Typst Wasm; stock CeTZ paints finished geometry. Kurvst owns
curve and mark geometry; Linnest owns placement decisions, label search, and
link semantics; CeTZ owns painting and its existing Boolean operations.
Eventually delete the vendored CeTZ tree: 47 files, five core patches, and its
pinned rebuild recipe. Do not port LGPL CeTZ patches into the new MIT engine.

The prerequisite Kurvst consolidation lands with this note: the shared curve
API, stroke outlines, arcs, pattern splitting, and removal of the `optimize`
and close-mode options. Delivery order is: consolidate, capture the P0 corpus
baseline against the consolidated base, then implement the phases below.

== Mark catalogue and data contract

`kurvst.mark` constructors produce data only. The Python `Mark(name, ...)` of
`linnet_py` and `symbolica.community.graph` mirrors the same names, fields, and
defaults and serializes the same dictionary; there is no Python `kurvst` module. Names, parameters, and defaults follow
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

Mark constructors, the native Python configuration and the Typst renderer
migration are implemented as sketched below.

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
from symbolica.community.graph import GraphStyleOptions, Mark, RenderSettings

flow = Mark("triangle", length="0.21cm", width="0.1575cm")
momentum = Mark("stealth", inset="25%", fit="bend")
settings = RenderSettings(
    style=GraphStyleOptions(
        edge_style={"flow-arrow": flow, "momentum-arrow": momentum}
    )
)
diagram.render(config=settings, momenta=True)
```

Native configuration gains `flow-arrow` and `momentum-arrow`. Reject CeTZ mark
dictionaries with an error showing the replacement; remove the Python `Mark` CeTZ
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

Each phase is its own change after the consolidated base, with corpus
timings recorded. S/M/L are planning sizes, not completion indicators.

=== P0 — S: reproducible corpus and visual policy

Check in the 24-diagram corpus, regenerated templates, median timing script,
and PNG references. Establish explicit lenient tolerances and image-diff review
for arrows and labels, *not blanket pixel identity*. Require cross-renderer
outline agreement and correct link targets and hit areas. Record provenance
(host, toolchain, cache, timing boundaries) and paired baseline comparison.
Regenerate this baseline only after consolidation.

The P0 harness, immutable input corpus, and reference images are recorded in
`tests/drawing-bench/`. The captured three-round baseline is 15.010 s on M4 Pro
/ Typst 0.15.1, with 24/24 reference PNGs byte-identical to the supplied
snapshot at 150 ppi; the consolidated sources reproduce all 24 PNGs
byte-for-byte on Linux. Source, package, host, and raw timing identities
accompany the images. Numerical tolerances for future arrow/label changes
still require explicit review; renderer migration is not claimed.

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

The native renderer uses one drawing-wide candidate call and one selected
call, with deduplicated templates and one painted shaft per carrier. Linnest's
SVG tests cover the call count, every catalogue shape and composites with both
fits, independent leaf paints and hit areas, and painted candidate bounds.
Native configuration reads `flow-arrow` and `momentum-arrow` from the graph
edge style (`GraphStyleOptions.edge_style`); linnet-py turns `Mark` values and their public
dictionaries into the numeric `MarkSpec` and its leaf paints on every native
path, for ordinary graphs and domain scenes alike. It checks them, and
dictionaries under the Typst `mark` key, when a style is configured, so a CeTZ
dictionary or an invalid mark is an error before either renderer is chosen. A
domain scene whose other settings need Typst draws the same configured heads
there. Named paints resolve through
Typst's own color constants, not CSS names. Painted head outlines count in the
native label-overlap score, so label feedback moves labels off heads as well as
shafts.

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

The Typst route consumes the shared candidate and selected batches. Physics
templates, the map-style example, linnet-py's scene template and Python
selectors use canonical mark data; CeTZ edge-mark dictionaries are rejected.
All Clinnet tests, including the native/Wasm mark parity checks, pass with
stock CeTZ 0.5.2 as the only package source. Named preparation
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

The authored fixed-size identity-target helper passes focused stock CeTZ checks,
preserving both expected SVG links. A document whose own canvas uses stock CeTZ
0.5.1 or 0.5.2 can embed Linnest edge halves, Kurvst paths and a full drawing
with Kurvst heads; both produce the same SVG.

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

The patched 47-file package tree, five patches and rebuild recipe are deleted,
and no CeTZ copy remains in the repository or the wheels. Nix selects stock
CeTZ 0.5.2 for the documentation and the development shell. linnet-py's Typst
route reads `@preview/cetz:0.5.2` from the local package path and cache, in the
Typst CLI's order, and never downloads; a missing package is reported with the
directories searched. Native rendering needs no Typst package.

== Validation and remaining acceptance

Measured on itphlies (AMD EPYC 9754, Linux, Typst 0.15.1) at a shared-host load
average near 70. The 24-case corpus compiles with stock CeTZ 0.5.2 as the only
package source. At 150 ppi, six PNGs equal P0 byte for byte and eighteen
differ. Every momentum case redraws its momentum heads, and particle-flow
heads are filled with the line paint instead of ink, with the engine's explicit
sizes and inset stroke geometry. Label search places labels differently in
some momentum cases: nonplanar K3,3 no longer overlaps two labels, the Higgs
scattering figure is 8 px shorter, and in the six-leg Higgs theta figure one
label now crosses a dotted edge. Image review remains required.

Paired, interleaved fresh-process compiles over five rounds compare these
sources with the P0 sources and the formerly vendored, patched CeTZ 0.5.1:
stock CeTZ is 1.167 times slower in the median round (range 1.164–1.224), with
sums of case medians of 32.6 s against 27.7 s. Momentum cases take 1.23–1.40
times as long; particle and notebook cases 0.97–1.02 times. Performance
acceptance is therefore still open; the M4 Pro ceilings above are not
comparable on this host.

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

P0–P7 are implemented. The performance targets and the image review above
remain open acceptance decisions.
