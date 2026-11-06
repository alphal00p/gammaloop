// vakint p-integrals: RustRed vs Forcer (benchmark v1). Every number is read
// from data/p_results.json (examples/benchmark_p_integrals.rs) and
// data/p_exact.json (the exact-agreement tier of the Forcer battery).
// Build: `just vakint-p-benchmark-pdf` in the dev shell (Typst 0.15.1).

#let results = json("data/p_results.json")
#let exact = json("data/p_exact.json")
#assert(results.schema == "vakint.p-benchmark.v1", message: "data/p_results.json is not vakint.p-benchmark.v1")
#assert(exact.schema == "vakint.p-exact.v1", message: "data/p_exact.json is not vakint.p-exact.v1")

// The selected run and every run of the same source and build: a resumed run
// continues the interrupted one.
#let selected = results.runs.find(run => run.id == results.selected_run)
#assert(
  selected != none or results.selected_run == none,
  message: "data/p_results.json has no run " + repr(results.selected_run),
)
#let protocol = if selected == none { (:) } else { selected.protocol }
#let same-source(run) = (
  selected != none
    and ("gammaloop_sha", "rustred_sha", "build_profile", "source_dirty").all(
      name => run.protocol.at(name) == protocol.at(name),
    )
)
#let measurements = results.runs.filter(same-source).map(run => run.measurements).flatten()

// A filled exact file marks only the source it was computed on: gammaloop_sha
// is the last commit outside docs/benchmark, rustred_sha the Cargo.lock pin.
#let exact-marks = selected != none and exact.at("stub", default: false) != true
#if exact-marks {
  for name in ("gammaloop_sha", "rustred_sha") {
    let (found, measured) = (exact.at(name, default: none), protocol.at(name))
    assert(
      found == measured,
      message: "p_exact.json has " + name + " " + repr(found) + ", the selected run " + measured,
    )
  }
}

// The deterministic date of the document: the start of the selected run.
#let stamp = protocol.at("started_utc", default: none)
#let date = if stamp == none { none } else {
  datetime(
    year: int(stamp.slice(0, 4)),
    month: int(stamp.slice(5, 7)),
    day: int(stamp.slice(8, 10)),
  )
}

#set document(title: "vakint p-integrals: RustRed vs Forcer", author: "Valentin Hirschi", date: date)
#set page(paper: "a4", margin: (x: 1.3cm, y: 1.4cm), numbering: "1 / 1")
#set text(font: "Libertinus Serif", size: 10pt, lang: "en")
#set par(justify: true, leading: 0.55em)
#show raw: set text(font: "DejaVu Sans Mono", size: 7.5pt)
#show figure.where(kind: table): set figure.caption(position: top)
#show figure.caption: set text(size: 9pt)

#let navy = rgb("#1f3b63")
#let pale = rgb("#eef3f9")
#let grey = luma(232)

// ---------------------------------------------------------------- helpers --

#let is-num(x) = type(x) == int or type(x) == float

#let fixed(x, d) = {
  let s = str(calc.round(x, digits: d))
  if d == 0 { return s }
  let parts = s.split(".")
  let frac = if parts.len() > 1 { parts.at(1) } else { "" }
  parts.at(0) + "." + frac + "0" * (d - frac.len())
}

// Significant figures suited to timings in seconds; never a printed zero.
#let fmt(x) = {
  if not is-num(x) { return "—" }
  let a = calc.abs(x)
  if a < 0.001 { return "< 0.001" }
  fixed(x, if a >= 100 { 0 } else if a >= 10 { 1 } else if a >= 1 { 2 } else { 3 })
}
#let fmt-ratio(r) = if r >= 10 { fixed(r, 0) } else if r >= 0.01 { fixed(r, 2) } else {
  "1/" + str(int(calc.round(1 / r)))
}

#let median(values) = {
  let sorted = values.sorted()
  let n = sorted.len()
  if n == 0 { none } else { (sorted.at(int((n - 1) / 2)) + sorted.at(int(n / 2))) / 2 }
}

// The measurements of one (cell, backend, threads, state).
#let samples(cell, backend, threads, state) = measurements.filter(
  m => m.cell == cell and m.backend == backend and m.threads == threads and m.state == state,
)
// Median over the repetitions; none without a complete measurement.
#let seconds(cell, backend, threads, state, field: "evaluate_seconds") = median(
  samples(cell, backend, threads, state)
    .filter(m => m.status == "ok")
    .map(m => m.at(field, default: none))
    .filter(is-num),
)
// The cap a state went over, the longest when several ran; none if it did not.
#let cap-of(cell, backend, threads, state) = {
  let caps = samples(cell, backend, threads, state)
    .filter(m => m.status == "capped")
    .map(m => m.cap_seconds)
  if caps.len() == 0 { none } else { calc.max(..caps) }
}
#let over-cap(seconds) = if seconds >= 60 {
  "> " + str(int(calc.round(seconds / 60))) + " min"
} else {
  "> " + str(seconds) + " s"
}
// One repetition at four loops (Forcer only), the protocol's below.
#let repetitions(loops) = if loops >= 4 { 1 } else { protocol.at("repetitions", default: 3) }
// A dagger on a median over fewer repetitions than the protocol's.
#let partial(cell, backend, threads, state, loops) = {
  let done = samples(cell, backend, threads, state).filter(m => m.status == "ok").len()
  if done > 0 and done < repetitions(loops) { "†" } else { "" }
}

#let ratio(a, b) = if is-num(a) and is-num(b) and a > 0 and b > 0 { a / b } else { none }

// Diverging log2 palette: green (RustRed faster), white (parity), red (slower).
#let palette = gradient.linear(
  rgb("#5aae61"), rgb("#d9f0d3"), rgb("#ffffff"), rgb("#fddbc7"), rgb("#d6604d"),
)
#let ratio-fill(r) = {
  if r == none { return none }
  // log2 scale: x1/4 -> 0 %, x1 -> 50 %, x4 -> 100 %
  palette.sample(calc.clamp(0.5 + calc.log(r, base: 2) / 4, 0, 1) * 100%)
}

// The mark of a grid cell, or of the physics cell of `loops` when `dots` is none.
#let exact-mark(loops, dots, rank) = {
  let found = if exact-marks { exact.cells.find(c => if dots == none {
    c.at("physics", default: none) not in (none, false) and c.L == loops
  } else {
    c.L == loops and c.at("dots", default: none) == dots and c.at("rank", default: none) == rank
  }) }
  if found == none { [] } else if found.agree {
    text(fill: rgb("#1b7a3e"), size: 7pt)[#h(2pt)✓]
  } else {
    text(fill: rgb("#b2182b"), size: 7pt, weight: "bold")[#h(2pt)✗]
  }
}

// The Forcer state of the table: cold, against RustRed cold; warm at four
// loops, the only state measured there.
#let forcer-state(loops) = if loops >= 4 { "warm" } else { "cold" }
#let forcer-time(cell, threads, loops) = seconds(cell, "forcer", threads, forcer-state(loops))
#let timing-text(cell, backend, threads, state, loops) = {
  let value = seconds(cell, backend, threads, state)
  let cap = cap-of(cell, backend, threads, state)
  if value != none {
    fmt(value) + partial(cell, backend, threads, state, loops)
  } else if cap != none { over-cap(cap) } else { "—" }
}
#let cell-count(cell, threads) = {
  let found = measurements.find(m => m.cell == cell)
  if found == none { none } else { found.n }
}

// One table cell: RustRed cold, Forcer, the ratio, the two warm states, the
// exact mark and the integral count. Four loops are Forcer-only in v1.
#let cmp-cell(cell, loops, threads, dots: none, rank: none) = {
  let own = seconds(cell, "rustred", threads, "cold")
  let r = ratio(own, forcer-time(cell, threads, loops))
  let n = cell-count(cell, threads)
  table.cell(fill: if loops >= 4 { grey } else { ratio-fill(r) }, {
    set par(leading: 0.3em, justify: false)
    if dots == none { text(size: 6.5pt, weight: "bold")[#loops L #h(2pt)] }
    if loops >= 4 {
      text(size: 7pt, fill: luma(90))[n/a (M3)]
    } else if own == none and cap-of(cell, "rustred", threads, "cold") == none {
      text(size: 7pt, fill: luma(90))[pending]
    } else {
      text(weight: "bold", size: 8.2pt, timing-text(cell, "rustred", threads, "cold", loops))
    }
    exact-mark(loops, dots, rank)
    linebreak()
    text(size: 6.5pt, fill: luma(40))[F #timing-text(cell, "forcer", threads, forcer-state(loops), loops)]
    if r != none { text(size: 6.8pt, weight: "semibold")[#h(2pt)×#fmt-ratio(r)] }
    if own != none {
      linebreak()
      text(size: 5.6pt, fill: luma(45))[
        ws #timing-text(cell, "rustred", threads, "warm-same", loops)
        #h(1pt) wc #timing-text(cell, "rustred", threads, "warm-shared", loops)
      ]
    }
    if n != none {
      linebreak()
      text(size: 5.6pt, fill: luma(60))[n=#n]
    }
  })
}

#let header(body) = text(font: "Roboto", size: 7.5pt, weight: "bold", fill: white, body)

#let timing-table(threads) = {
  let ranks = (0, 2, 4, 6, 8)
  let rows = ()
  for loops in (1, 2, 3, 4) {
    for dots in (0, 1, 2) {
      if dots == 0 {
        rows.push(table.cell(rowspan: 3, align: horizon, text(weight: "bold")[#loops L]))
      }
      rows.push(text(size: 8pt)[#dots])
      for rank in ranks {
        rows.push(cmp-cell(
          "L" + str(loops) + "_d" + str(dots) + "_r" + str(rank), loops, threads, dots: dots, rank: rank,
        ))
      }
    }
  }
  rows.push(table.cell(colspan: 2, align: horizon, text(size: 8pt, weight: "bold")[phys.]))
  for loops in (2, 3, 4) {
    rows.push(cmp-cell("physics_L" + str(loops), loops, threads))
  }
  rows.push(table.cell(colspan: 2, align: horizon + center, fill: grey, text(size: 7pt, fill: luma(90))[5L canary: n/a]))
  table(
    columns: (1.05cm, 0.75cm, ..((1fr,) * 5)),
    align: (x, y) => if x < 2 { center + horizon } else { left + horizon },
    inset: (x: 3pt, y: 2.2pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy },
    header[L], header[dots],
    ..ranks.map(rank => header[rank #rank]),
    ..rows,
  )
}

// Ratios of the 3L cells of rank >= 4 at 16 threads (G10b).
#let g10b = {
  let values = ()
  for dots in (0, 1, 2) {
    for rank in (4, 6, 8) {
      let cell = "L3_d" + str(dots) + "_r" + str(rank)
      let r = ratio(seconds(cell, "rustred", 16, "cold"), forcer-time(cell, 16, 3))
      if r != none { values.push(r) }
    }
  }
  values
}
// A capped or missing cell has no ratio, so the range says how many it covers.
#let range-of(values, cells) = if values.len() == 0 { [pending] } else {
  [×#fmt-ratio(calc.min(..values)) to ×#fmt-ratio(calc.max(..values)), #values.len() of #cells cells]
}
#let load-seconds = (
  measurements
    .filter(m => m.backend == "rustred" and m.status == "ok")
    .map(m => m.at("load_seconds", default: none))
    .filter(is-num)
)

// ------------------------------------------------------------------ title --

#align(center)[
  #text(19pt, weight: "bold", fill: navy)[vakint p-integrals: RustRed vs Forcer]
  #v(-0.5em)
  #text(9pt, fill: luma(100))[
    #if date == none [no selected run] else { date.display("[year]-[month]-[day]") }
    #if selected != none [· #protocol.host · gammaloop #protocol.gammaloop_sha.slice(0, 9)#if protocol.source_dirty [ (dirty)] · rustred #protocol.rustred_sha.slice(0, 8)]
  ]
]

#block(fill: pale, inset: 7pt, radius: 3pt, width: 100%)[
  #set text(8.6pt)
  *Reading guide.* Each cell is the sum over all tops of that loop count at the given dots, times a
  numerator of that rank. Bold: RustRed cold (seconds, median). F: Forcer. ×: the ratio RustRed cold
  over Forcer, cold against cold because neither replays results; green is RustRed faster.
  ws, wc: RustRed warm on its own cache, and warm on one cache shared along the cells of that loop
  count. ✓: result exactly equal to Forcer's. n: integrals in the input. †: fewer repetitions than
  the protocol. Grey: Forcer only in v1.
]

#block(stroke: 0.6pt + luma(200), inset: 6pt, radius: 3pt, width: 100%)[
  #set text(8.4pt)
  *Performance gates.* #h(1em)
  4L on par with Forcer at 1 thread, better at 16: n/a in v1 (RustRed 4L: M3). #h(1em)
  Cold start dominated by loading: load #if load-seconds.len() == 0 [pending] else [#fmt(median(load-seconds)) s (median)]. #h(1em)
  G10b, 3L rank ≥ 4, 16 threads, RustRed cold ≤ Forcer (recorded, not blocking): #range-of(g10b, 9).
]

#v(-0.3em)
#align(center, stack(
  dir: ltr,
  spacing: 4pt,
  text(size: 7pt)[RustRed faster ×1/4],
  box(width: 5cm, height: 7pt, fill: gradient.linear(
    rgb("#5aae61"), rgb("#d9f0d3"), rgb("#ffffff"), rgb("#fddbc7"), rgb("#d6604d"),
  )),
  text(size: 7pt)[×4 Forcer faster],
))
#figure(timing-table(1), caption: [1 thread (`n_cores` 1, `form`).], supplement: [Table])
#figure(timing-table(16), caption: [16 threads (`n_cores` 16, `tform -w16`).], supplement: [Table])

// ------------------------------------------------------------ methodology --

// A sum over a part of the measured cells would read as the whole: a cell
// without the sub-line (no math span measured) leaves a dash.
#let per-loops(backend, state, field) = (1, 2, 3, 4).map(loops => {
  let cells = (
    measurements
      .filter(m => (
        m.loops == loops and m.backend == backend and m.threads == 1 and m.state == state and m.status == "ok"
      ))
      .map(m => m.cell)
      .dedup()
  )
  let values = cells
    .map(cell => seconds(cell, backend, 1, state, field: field))
    .filter(is-num)
  if values.len() == 0 or values.len() < cells.len() { [—] } else { [#fmt(values.sum())] }
})

#block(breakable: false)[
  #set text(8.6pt)
  *Sub-lines* (seconds, sum over the cells measured at 1 thread).
  #table(
    columns: (auto, 1fr, 1fr, 1fr, 1fr),
    stroke: 0.4pt + luma(215),
    inset: (x: 4pt, y: 2pt),
    fill: (x, y) => if y == 0 { navy },
    header[], ..(1, 2, 3, 4).map(loops => header[#loops L]),
    [RustRed reduction only], ..per-loops("rustred", "cold", "reduction_seconds"),
    [RustRed load (excluded)], ..per-loops("rustred", "cold", "load_seconds"),
    [Forcer FORM subprocess], ..per-loops("forcer", "cold", "reduction_seconds"),
    [Raw FORM reference], ..per-loops("forcer", "raw", "evaluate_seconds"),
  )

  *Method.* Timed: the wall time of `Vakint::evaluate` in Retain mode (exact reduction onto masters,
  same frontend and finalizer for both backends) in one fresh process per measurement, with the
  input parsed, the engine and library loaded and the 1L bubble evaluated once beforehand, by either
  backend, on directories of its own; repetitions:
  #protocol.at("repetitions", default: 3), median. RustRed cold: empty ordered-cache directory, shipped caches off. ws: a fresh
  process on a read-only copy of the cache of its cold run. wc: per loop count, one chain over the
  cells in table order sharing one cache from empty, each cell measured once in the chain. Forcer
  cold: private empty materialisation directory; Forcer warm: materialised. Forcer has no result
  cache. The raw reference runs the generated FORM programs directly, one process per top. OS pages
  are not flushed. 4L and the 4L physics cell: Forcer only, 16 threads, warm, one repetition, capped
  per cell. Physics row: massless-quark (L-1)-rung gluon self-energy, Feynman gauge, transverse
  projector, external momentum kept.

  *Provenance.* #if selected == none [No selected run.] else [
    Run #raw(selected.id), profile #protocol.build_profile#if protocol.source_dirty [, dirty tree],
    #raw(protocol.form_version.split("\n").first()), Forcer #protocol.forcer_commit.slice(0, 7),
    #protocol.cpu_model on #protocol.host,
    CPUs #protocol.cpus.pairs().map(pair => pair.at(0) + " thread: " + pair.at(1)).join("; "),
    binding #if protocol.numactl == none [none] else { raw(protocol.numactl) },
    RustRed arithmetic #protocol.ordered_weight_arithmetic.pairs().map(pair => pair.at(0) + " thread: " + pair.at(1).split("(").first()).join("; "),
    #raw(str(measurements.len())) measurements in `data/p_results.json`.
  ]
]
