= Kurvst Arrowhead Engine

#quote(block: true)[
*Status: Accepted design, recorded 2026-10-04; P0–P7 NOT IMPLEMENTED.*
This is a planned architecture note, not an executable API reference or a claim
that renderer migration has landed. API examples below are proposed sketches.
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
Delivery order is: finish consolidation, rebase feynkit, then regenerate the P0
corpus baseline against the reconciled base. This documentation task performs
no rebase, implementation, CI change, or generated-asset update.

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

== Proposed API and migration (not executable today)

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

== Delivery phases: P0 baseline captured, P1 implemented, P2–P7 NOT IMPLEMENTED

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
constructors for thirteen shapes and composites. Its 93 native tests and
native/Wasm parity checks pass. The approved parity tolerance is
`1e-12 * max(1, abs(actual), abs(expected))` for geometry coordinates and
lengths only; structure, styles, indices, and other data remain exact.
Candidate results remain independent, while selected results share one
authoritative painted shaft per carrier. No renderer has migrated yet.

=== P2 — M: native and Python

Replace hardcoded `linnest::svg::marks` with the engine, add native
`flow-arrow` / `momentum-arrow` and the Python mirror. Preserve the default
notebook look with explicit sizes.

=== P3 — M: Typst painting

Replace `mark.geometry` with engine geometry painted through CeTZ; accept only
Kurvst marks. Migrate physics, manual, and map-style uses. Remove one of the
three vendored-only calls.

=== P4 — M: candidate footprints

Replace `mark.footprints` and the slow paint fallback with shared candidate-batch
engine footprints. This addresses the main 79% momentum penalty.

=== P5 — S: hit areas

Use numeric hitboxes and prepared content drawables instead of per-target CeTZ
`content-many`. Address roughly 75% of the plain penalty and 16% of momentum's;
keep link semantics separate from the catalogue.

=== P6 — S: collision sampling

Land the collision-sampling prototype JJ change `tvxztqnq` with one plugin batch.
Its reported vendored-corpus improvement was 1.8%; collision accounts for
roughly 2–14% of the stock penalty. These figures are not new measurements.

=== P7 — M: stock CeTZ and vendor removal

Target *stock CeTZ 0.5.2* (the measurements above used 0.5.1). Delete vendor,
patches, and `REBUILD`; migrate behavior fixtures and Nix packaging; rerun the
corpus. Generic bounds/sample speedups may optionally be upstreamed.

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
