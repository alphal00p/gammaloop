// Render from the checked measurement record; never edit timing values here.
// Visual reference: https://github.com/alphal00p/ninjars/blob/main/docs/benchmark/benchmark.pdf
#let report = json("data/results.json")
#let selected = report.at("selected_reduction_run", default: none)
#let data = if selected == none { report } else { report.reduction_runs.at(selected) }
#assert(data.validation_passed, message: "Only validated benchmark results may be published")
#let navy = rgb("#1f3b63")
#let pale = rgb("#eef3f9")

#set document(title: "Vakint reduction benchmark")
#set page(paper: "a4", margin: (x: 1.6cm, y: 1.7cm), numbering: "1 / 1")
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em)
#show raw: set text(font: "DejaVu Sans Mono", size: 8pt)
#show link: set text(fill: rgb("#1f5fa8"))

#let fixed(x, d) = {
  let parts = str(calc.round(x, digits: d)).split(".")
  if d == 0 { return parts.first() }
  let fraction = if parts.len() > 1 { parts.at(1) } else { "" }
  parts.first() + "." + fraction + "0" * (d - fraction.len())
}
#let timing(x) = fixed(x, if x >= 10 { 1 } else if x >= 1 { 2 } else if x >= 0.1 { 3 } else { 4 })
#let multiplier(x) = fixed(x, if x < 0.1 { 3 } else { 2 })

// Same logarithmic green–white–red scale as the NinjaRS reference report.
#let palette = gradient.linear(
  rgb("#5aae61"), rgb("#d9f0d3"), rgb("#ffffff"), rgb("#fddbc7"), rgb("#d6604d"),
)
#let ratio-fill(r) = {
  let t = calc.clamp(0.5 + calc.log(r, base: 2) / 4, 0, 1)
  palette.sample(t * 100%)
}
#let comparison(family, kind) = {
  let row = data.summary.find(row => row.family == family and row.kind == kind)
  let native = row.cached.median_seconds
  let reference = row.fmft.median_seconds
  let ratio = native / reference
  table.cell(fill: ratio-fill(ratio), {
    set par(leading: 0.35em)
    text(size: 17pt, weight: "bold")[#timing(native) #text(size: 10pt)[s]]
    linebreak()
    text(size: 9pt, fill: luma(65))[FMFT #timing(reference) s]
    h(7pt)
    text(size: 10pt, weight: "semibold")[×#multiplier(ratio)]
  })
}

#align(center)[
  #v(0.4em)
  #text(23pt, weight: "bold", fill: navy)[Vakint reduction benchmark]
  #v(-0.2em)
  #text(12pt, fill: luma(80))[Cached RustRed versus the FORM / FMFT reduction route]
  #v(0.2em)
  #text(9.5pt, fill: luma(110))[
    #data.completed_utc.split("T").first() · #data.protocol.cpu_model · CPU affinity #str(data.protocol.cpu_affinity)
  ]
]

#v(0.6em)
#block(fill: pale, inset: 10pt, radius: 4pt, width: 100%)[
  #set text(9.5pt)
  *Reading guide.* Bold values are RustRed times; FMFT times appear below.
  Multipliers are _RustRed ÷ FMFT_: green means faster, white parity, red slower.
  Both receive identical four-loop scalar inputs. The first physical denominator
  has power 2 (base) or 4 (held out); all others have power 1, auxiliaries 0.
]

#v(0.5em)
#text(14pt, weight: "bold", fill: navy)[Reduction time]
#v(0.3em)
Median seconds over #data.protocol.repetitions fresh processes per backend and case.

#v(0.35em)
#{
  set text(font: "DejaVu Sans", size: 9pt)
  set par(justify: false)
  show table.cell.where(y: 0): set text(fill: white, weight: "bold")
  table(
    columns: (0.7fr, 1.5fr, 1.5fr),
    align: (x, y) => if x == 0 { left + horizon } else { center + horizon },
    inset: (x: 12pt, y: 12pt),
    stroke: (x, y) => (
      top: if y <= 1 { none } else { 0.4pt + luma(225) },
      left: if x >= 2 and y > 0 { 0.4pt + luma(225) } else { none },
    ),
    fill: (x, y) => if y == 0 { navy } else { none },
    table.header([Topology], [Base · power 2], [Held out · power 4]),
    ..for family in ("H", "FG", "BMW", "X") {
      (
        text(size: 13pt, weight: "semibold", family),
        comparison(family, "base"),
        comparison(family, "held-out"),
      )
    },
  )
}
#v(0.45em)
#align(center)[
  #set text(8pt)
  #box(width: 9cm, height: 0.3cm, fill: palette, radius: 2pt)
  #linebreak()
  #box(width: 9cm)[≤ ×1/4 #h(1fr) ×1/2 #h(1fr) ×1 #h(1fr) ×2 #h(1fr) ≥ ×4]
]

#v(0.6em)
#block(stroke: 0.6pt + luma(200), inset: 10pt, radius: 4pt, width: 100%)[
  #set text(9.5pt)
  *Method and correctness.*
  #if data.protocol.at("configured_profile", default: none) != none [
    Release build: optimization level #data.protocol.configured_profile.opt_level,
    #data.protocol.configured_profile.lto LTO;
  ] else [
    Build profile: #data.protocol.build_profile (caller-declared);
    optimization/LTO settings are not inferred.
  ]
  #(data.protocol.rustc.split(" (").first() + "; " +
    data.protocol.form_version.split("\n").first().split(" (").first() + ".")
  Both retain masters in FMFT/MATAD normalization, with identical symbolic mass.
  Exact rational parity passed all eight cases at unit mass and renormalization
  scale, without truncating epsilon or using numerical tolerances.

  Base cases hit prepared roots without applying rules. Held-out cases reuse
  cached descendants after root misses. Timing covers first public reduction
  calls, including FORM version checks, launch, I/O and parsing. Data loading,
  cache attachment, master evaluation and validation are excluded. Fresh
  processes run sequentially, alternating backends, with the recorded CPU
  affinity on a shared host.
  Filesystem pages are warm; this does not measure cold startup.
]

#v(0.5em)
#set text(8.5pt, fill: luma(100))
#if selected == none [Historical measurements, unchanged.] else [Selected validated run: #raw(selected).]
Source: #raw(data.protocol.source_commit.slice(0, 12)).
#if data.protocol.source_dirty [Modified source tree.]
Data: #raw(data.protocol.data_version).
Samples and provenance: #raw("data/results.json").
Refresh instructions: #raw("../../AGENTS.md").

#pagebreak()
#let ingress = report.five_loop_dot_ingress_acceptance
#let pending = report.five_loop_performance_status
#assert(ingress.validation_passed, message: "Five-loop ingress checks must pass")
#set text(font: "Libertinus Serif", size: 10pt, fill: black)

#align(center)[
  #text(22pt, weight: "bold", fill: navy)[Five-loop validation and benchmarks]
  #v(0.2em)
  #text(11.5pt, fill: luma(80))[640 declared, uncomputed masters · performance measurements pending]
]

#v(0.55em)
#text(13pt, weight: "bold", fill: navy)[Historical DOT ingress acceptance]
#v(0.2em)
#ingress.completed_utc.split("T").first() · one sequential debug-build pass.
This checks graph definitions; it does not time IBP reduction, master evaluation,
package loading or cache preparation.

#v(0.35em)
#{
  set text(font: "DejaVu Sans", size: 9pt)
  set par(justify: false)
  table(
    columns: (2.1fr, 1fr),
    inset: (x: 11pt, y: 8pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else if calc.rem(y, 2) == 0 { pale } else { none },
    table.header(
      text(fill: white, weight: "bold")[Validation workload],
      text(fill: white, weight: "bold")[Recorded result],
    ),
    [Master graphs / dotted / with numerator],
    [#ingress.graphs.masters_checked / #ingress.graphs.dotted_masters / #ingress.graphs.numerator_masters],
    [Exact routing, full numerators, native ingress and DOT round trips], [Passed],
    [Complete graph-validation pass],
    text(size: 16pt, weight: "bold", fill: navy)[#timing(ingress.graphs.elapsed_seconds) s],
    [Per graph: minimum / median / maximum],
    [#timing(ingress.graphs.minimum_seconds) / #timing(ingress.graphs.median_seconds) / #timing(ingress.graphs.maximum_seconds) s],
    [Numerical master values / prepared 5L cache in this run], [Absent / absent],
  )
}

#v(0.3em)
#set text(9pt)
Physical powers, signed indices, mass terms, full numerators, loop-basis slots
(including zero) and integration Jacobians were checked exactly. These checks
establish the stored graph definitions, without claiming that the catalogue is
an independent master basis. One observed run supports no statistical speedup
claim. Earlier interrupted checks are not a timing baseline.

#v(0.5em)
#text(13pt, weight: "bold", fill: navy)[Five-loop reduction comparison]
#v(0.25em)
#{
  set text(font: "DejaVu Sans", size: 9pt)
  set par(justify: false)
  table(
    columns: (1.8fr, 1fr, 1fr),
    inset: (x: 10pt, y: 11pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else if x > 0 { rgb("#fff3dc") } else { pale },
    table.header(
      text(fill: white, weight: "bold")[Retained-master workload],
      text(fill: white, weight: "bold")[No persistent cache],
      text(fill: white, weight: "bold")[With persistent cache],
    ),
    ..for case in pending.cases {
      (case.label, [Pending], [Pending])
    },
  )
}

#v(0.35em)
#block(fill: pale, inset: 9pt, radius: 4pt, width: 100%)[
  #set text(9pt)
  *Required before filling these cells.* Complete reductions and exact result/guard
  agreement; a completed, strict cache workload; three fresh-process runs per
  mode with fixed affinity and worker counts. Report graph preparation, loading,
  lowering, IBP, cache loading and end-to-end time separately. Cache size,
  preparation cost, intermediate reuse and held-out regressions accompany the
  comparison. No cache speedup is claimed yet.
]

#v(0.35em)
#set text(8pt, fill: luma(100))
Historical ingress: source #raw(ingress.protocol.source_commit.slice(0, 12))
(modified); RustRed #raw(ingress.protocol.rustred_commit.slice(0, 12));
Symbolica #raw(ingress.protocol.symbolica_commit.slice(0, 12)).
Exact patch, executable, data and per-graph digests remain in JSON.

#v(0.2em)
Latest validated integration cohort, separate from these historical timings:
RustRed #raw(pending.validated_compatibility_cohort.rustred_revision.slice(0, 12));
Symbolica #raw(pending.validated_compatibility_cohort.symbolica_numerica_graphica_revision.slice(0, 12)).
Published upgrade: Symbolica
#raw(pending.pending_upgrade.symbolica_numerica_graphica_revision.slice(0, 12)),
FeynKit #raw(pending.pending_upgrade.feynkit_owner_source_revision.slice(0, 12));
final coordinated acceptance pending. Full revisions: #raw("data/results.json").
