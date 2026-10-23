// Render from the checked measurement record; never edit timing values here.
// Visual reference: https://github.com/alphal00p/ninjars/blob/main/docs/benchmark/benchmark.pdf
#let data = json("data/results.json")
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
  *Method and correctness.* Release build: optimization level 3, fat LTO;
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
Historical measurements, unchanged. Source: #raw(data.protocol.source_commit.slice(0, 12)).
#if data.protocol.source_dirty [Modified source tree.]
Data: #raw(data.protocol.data_version).
Samples and provenance: #raw("data/results.json").
Refresh instructions: #raw("../../AGENTS.md").

#pagebreak()
#let ingress = data.five_loop_dot_ingress_acceptance
#assert(ingress.validation_passed, message: "Five-loop ingress checks must pass")
#set text(font: "Libertinus Serif", size: 10.5pt, fill: black)

#align(center)[
  #text(22pt, weight: "bold", fill: navy)[Five-loop DOT ingress acceptance]
  #v(0.2em)
  #text(12pt, fill: luma(80))[640 unevaluated masters · debug verification]
  #v(0.2em)
  #text(9.5pt, fill: luma(110))[
    #ingress.completed_utc.split("T").first() · #ingress.protocol.cpu_model
  ]
]

#v(0.65em)
#block(fill: pale, inset: 11pt, radius: 4pt, width: 100%)[
  *Measurement boundary.* This is one sequential verification pass through all
  published five-loop master graphs, using the native FeynKit parser and
  Vakint input conversion. It is a debug build with exact assertions enabled.
  It does not time package loading, IBP reduction, master evaluation or cache
  preparation. The release/FMFT comparison on page 1 retains its original
  date, source revision and measurements.
]

#v(0.65em)
#{
  set text(font: "DejaVu Sans", size: 9.5pt)
  set par(justify: false)
  table(
    columns: (2.1fr, 1fr),
    inset: (x: 11pt, y: 9pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy } else if calc.rem(y, 2) == 0 { pale } else { none },
    table.header(
      text(fill: white, weight: "bold")[Validation workload],
      text(fill: white, weight: "bold")[Result],
    ),
    [Master graphs checked], [#ingress.graphs.masters_checked / 640],
    [Graphs with dotted propagators], [#ingress.graphs.dotted_masters],
    [Graphs with numerator factors], [#ingress.graphs.numerator_masters],
    [Exact routing and complete numerator checks], [Passed],
    [Native DOT round trips and Vakint ingress], [Passed],
    [Five-loop numerical values / prepared cache], [Absent / absent],
  )
}

#v(0.65em)
#block(fill: pale, inset: 12pt, radius: 4pt, width: 100%)[
  #text(17pt, weight: "bold", fill: navy)[
    #timing(ingress.graphs.elapsed_seconds) s
  ]
  #h(8pt) complete graph-validation pass
  #linebreak()
  #text(9.5pt)[
    Per graph: #timing(ingress.graphs.minimum_seconds)–#timing(ingress.graphs.maximum_seconds) s;
    median #timing(ingress.graphs.median_seconds) s.
    One observed run; no statistical performance claim.
  ]
]

#v(0.55em)
*What was checked.* Each stored family and signed index tuple is mapped exactly
to the graph's physical propagators, powers, mass terms and complete numerator.
Loop-basis annotations, including slot zero, and the integration-variable
Jacobian are verified. The same integral is then reconstructed through native
Vakint ingress and a FeynKit DOT round trip. These checks establish the graph
definitions; they do not claim linear independence of the 640 retained masters.

#v(0.4em)
*Implementation relevant to this pass.* Vakint now substitutes contracted
scalar products directly inside factored numerators. FeynKit accepts identical
stored numerator expressions before attempting expanded comparison. Both
changes avoid unnecessary expansion while retaining exact fallback checks.
Earlier interrupted checks are not a completed timing baseline, so no speedup
factor is inferred from them.

#v(0.5em)
#set text(8.5pt, fill: luma(100))
Consumer base: #raw(ingress.protocol.source_commit.slice(0, 12)) (modified source;
exact patch and executable digests recorded in JSON).
RustRed: #raw(ingress.protocol.rustred_commit.slice(0, 12));
Symbolica: #raw(ingress.protocol.symbolica_commit.slice(0, 12)).
Build: #ingress.protocol.build_profile.
Measured receipt and data manifest digests:
#raw("data/results.json → five_loop_dot_ingress_acceptance").
