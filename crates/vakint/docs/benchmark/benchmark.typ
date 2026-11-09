// vakint p-integrals: RustRed vs Forcer (benchmark v1). Every number is read
// from data/p_results.json (examples/benchmark_p_integrals.rs) and
// data/p_exact.json (the exact-agreement tier of the Forcer battery).
// Build: `just vakint-p-benchmark-pdf` in the dev shell (Typst 0.15.1).

#let results = json("data/p_results.json")
#let exact = json("data/p_exact.json")
#assert(results.schema == "vakint.p-benchmark.v1", message: "data/p_results.json is not vakint.p-benchmark.v1")
#assert(exact.schema == "vakint.p-exact.v1", message: "data/p_exact.json is not vakint.p-exact.v1")

// The selected run and every run of the same source and build: a resumed run
// continues the interrupted one. The source includes the identity of the
// RustRed p library (the built-in one, or the one external data list).
#let selected = results.runs.find(run => run.id == results.selected_run)
#assert(
  selected != none or results.selected_run == none,
  message: "data/p_results.json has no run " + repr(results.selected_run),
)
#let protocol = if selected == none { (:) } else { selected.protocol }
#let library-of(p) = p.at("p_library", default: none)
#let identity-of(p) = { let library = library-of(p); if library == none { none } else { library.identity } }
#let same-source-as(base) = run => (
  ("gammaloop_sha", "rustred_sha", "build_profile", "source_dirty").all(
    name => run.protocol.at(name) == base.at(name),
  )
    and identity-of(run.protocol) == identity-of(base)
)
// Each measurement carries the repetitions its run planned (the dagger below).
#let gather(base) = results.runs.filter(same-source-as(base)).map(run => run.measurements.map(
  m => m + (planned: run.protocol.repetitions),
)).flatten()
// The provenance prints a run's CPUs for every run merged with it: each of
// them must have pinned every thread count it measured.
#let check-cpus(base) = for run in results.runs.filter(same-source-as(base)) {
  assert(
    run.protocol.threads.all(t => (
      str(t) in run.protocol.cpus and run.protocol.cpus.at(str(t)) == base.cpus.at(str(t), default: none)
    )),
    message: "run " + run.id + " has a thread count without --cpus, or CPUs other than those of " + repr(base.at("started_utc", default: none)),
  )
}
#let measured = if selected == none { () } else { gather(protocol) }
#if selected != none { check-cpus(protocol) }

// The cited run (benchmark v1, M2-T10c): its source fills what the selected
// source does not measure: every Forcer row and RustRed at 1L-3L.
#let cited-run = results.runs.find(run => run.id == results.at("cited_run", default: none))
#assert(
  cited-run != none or results.at("cited_run", default: none) == none,
  message: "data/p_results.json cites no recorded run " + repr(results.at("cited_run", default: none)),
)
#let cited-protocol = if cited-run == none { none } else { cited-run.protocol }
#let cited = if cited-run == none { () } else { gather(cited-protocol) }
#if cited-run != none { check-cpus(cited-protocol) }
#let measurements = measured + cited

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

// RustRed's Ordered weight arithmetic of the tables; the exact series has a
// table of its own. Forcer measurements carry "none".
#let headline = "modular_polynomial"
#let arithmetic-of(m) = { assert(m.backend != "rustred" or "arithmetic" in m, message: "RustRed record without arithmetic"); m.at("arithmetic", default: "none") }

#let loops-of(cell) = if cell.starts-with("physics_L") { int(cell.slice(9)) } else { int(cell.slice(1, 2)) }
// Where a cell of the tables comes from: RustRed at 1L-3L from the cited
// source (the v1 grid, built-in 1-3 library), at 4L from the measured one;
// Forcer from the measured source where it measured the cell, else the cited.
// Without a cited run everything is the selected source's.
#let rustred-rows(loops) = if loops < 4 and cited-run != none { cited } else { measured }
#let pick(list, cell, backend, threads, state, arithmetic) = list.filter(
  m => (
    m.cell == cell and m.backend == backend and m.threads == threads and m.state == state
      and (backend != "rustred" or arithmetic-of(m) == arithmetic)
  ),
)
#let samples(cell, backend, threads, state, arithmetic: headline) = if backend == "rustred" {
  pick(rustred-rows(loops-of(cell)), cell, backend, threads, state, arithmetic)
} else {
  let own = pick(measured, cell, backend, threads, state, arithmetic)
  if own.len() > 0 or cited-run == none { own } else { pick(cited, cell, backend, threads, state, arithmetic) }
}
// Median over the repetitions; none without a complete measurement.
#let median-of(list, field) = median(list.filter(m => m.status == "ok").map(m => m.at(field, default: none)).filter(is-num))
#let seconds(cell, backend, threads, state, field: "evaluate_seconds", arithmetic: headline) = median-of(
  samples(cell, backend, threads, state, arithmetic: arithmetic),
  field,
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
// A dagger on a median over fewer repetitions than its run planned.
#let partial(cell, backend, threads, state) = {
  let all = samples(cell, backend, threads, state)
  let done = all.filter(m => m.status == "ok").len()
  if done > 0 and done < calc.max(..all.map(m => m.planned)) { "†" } else { "" }
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

// The mark of a grid cell, or of the physics cell of `loops` when `dots` is
// none; n.c. for a cell the exact tier lists as not exact-checked (its cost).
#let same-cell(loops, dots, rank) = c => if dots == none {
  c.at("physics", default: none) not in (none, false) and c.L == loops
} else {
  c.L == loops and c.at("dots", default: none) == dots and c.at("rank", default: none) == rank
}
#let exact-mark(loops, dots, rank) = {
  let found = if exact-marks { exact.cells.find(same-cell(loops, dots, rank)) }
  let unchecked = exact-marks and exact.at("not_exact_checked", default: ()).any(same-cell(loops, dots, rank))
  if found == none and unchecked { text(fill: luma(100), size: 6pt)[#h(2pt)n.c.] } else if found == none { [] } else if found.agree {
    text(fill: rgb("#1b7a3e"), size: 7pt)[#h(2pt)✓]
  } else {
    text(fill: rgb("#b2182b"), size: 7pt, weight: "bold")[#h(2pt)✗]
  }
}

// The Forcer state of a cell: cold, against RustRed cold; warm where only
// that was measured (v1 measured 4L warm only).
#let forcer-state(cell, threads) = if seconds(cell, "forcer", threads, "cold") == none and seconds(cell, "forcer", threads, "warm") != none { "warm" } else { "cold" }
#let forcer-time(cell, threads) = seconds(cell, "forcer", threads, forcer-state(cell, threads))
#let timing-text(cell, backend, threads, state) = {
  let value = seconds(cell, backend, threads, state)
  let cap = cap-of(cell, backend, threads, state)
  if value != none {
    fmt(value) + partial(cell, backend, threads, state)
  } else if cap != none { over-cap(cap) } else { "—" }
}
#let cell-count(cell) = {
  let found = measurements.find(m => m.cell == cell)
  if found == none { none } else { found.n }
}
#let cold-ratio(cell, threads) = ratio(seconds(cell, "rustred", threads, "cold"), forcer-time(cell, threads))

// One table cell: RustRed cold, Forcer, the ratio, the two warm states, the
// exact mark and the integral count.
#let cmp-cell(cell, loops, threads, dots: none, rank: none) = {
  let own = seconds(cell, "rustred", threads, "cold")
  let r = cold-ratio(cell, threads)
  let n = cell-count(cell)
  let forcer = forcer-state(cell, threads)
  table.cell(fill: ratio-fill(r), {
    set par(leading: 0.3em, justify: false)
    if dots == none { text(size: 6.5pt, weight: "bold")[#loops L #h(2pt)] }
    if own == none and cap-of(cell, "rustred", threads, "cold") == none {
      text(size: 7pt, fill: luma(90))[pending]
    } else {
      text(weight: "bold", size: 8.2pt, timing-text(cell, "rustred", threads, "cold"))
    }
    exact-mark(loops, dots, rank)
    linebreak()
    text(size: 6.5pt, fill: luma(40))[F#if forcer == "warm" [w] #timing-text(cell, "forcer", threads, forcer)]
    if r != none { text(size: 6.8pt, weight: "semibold")[#h(2pt)×#fmt-ratio(r)] }
    if own != none {
      linebreak()
      text(size: 5.6pt, fill: luma(45))[
        ws #timing-text(cell, "rustred", threads, "warm-same")
        #h(1pt) wc #timing-text(cell, "rustred", threads, "warm-shared")
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

// Ratios of the 3L cells of rank >= 4 at 16 threads (G10b of M2).
#let g10b = {
  let values = ()
  for dots in (0, 1, 2) {
    for rank in (4, 6, 8) {
      let r = cold-ratio("L3_d" + str(dots) + "_r" + str(rank), 16)
      if r != none { values.push(r) }
    }
  }
  values
}
// A capped or missing cell has no ratio, so the range says how many it covers.
#let range-of(values, cells) = if values.len() == 0 { [pending] } else {
  [×#fmt-ratio(calc.min(..values)) to ×#fmt-ratio(calc.max(..values)), #values.len() of #cells cells]
}
#let load-seconds(list) = (
  list
    .filter(m => m.backend == "rustred" and m.status == "ok" and arithmetic-of(m) == headline)
    .map(m => m.at("load_seconds", default: none))
    .filter(is-num)
)

// The 4L speed target (G10b of M3; recorded, blocking only in M5), per cell:
// on par with Forcer at 1 thread, better at 16 threads, and a cold start
// dominated by the artefact load (the load of a fresh process at least the
// one-off cost of a cold evaluation over a warm one, cold - ws).
#let par-tolerance = 1.1
#let four-loop-cells = (0, 1, 2).map(dots => (0, 2, 4, 6, 8).map(rank => "L4_d" + str(dots) + "_r" + str(rank))).flatten() + ("physics_L4",)
#let verdict(ok) = if ok == none { text(fill: luma(90))[—] } else if ok {
  text(fill: rgb("#1b7a3e"), weight: "bold")[met]
} else { text(fill: rgb("#b2182b"), weight: "bold")[missed] }
#let target-1(cell) = { let r = cold-ratio(cell, 1); if r == none { none } else { r <= par-tolerance } }
#let target-16(cell) = { let r = cold-ratio(cell, 16); if r == none { none } else { r < 1 } }
#let one-off(cell, threads) = {
  let (cold, warm) = (seconds(cell, "rustred", threads, "cold"), seconds(cell, "rustred", threads, "warm-same"))
  if is-num(cold) and is-num(warm) { calc.max(cold - warm, 0) } else { none }
}
#let cell-load(cell, threads) = seconds(cell, "rustred", threads, "cold", field: "load_seconds")
#let target-load(cell) = {
  let (load, extra) = (cell-load(cell, 1), one-off(cell, 1))
  if is-num(load) and is-num(extra) { load >= extra } else { none }
}
#let tally(check) = {
  let found = four-loop-cells.map(check).filter(ok => ok != none)
  if found.len() == 0 [pending] else [#found.filter(ok => ok).len() of #found.len() met]
}
// The breakdown of a missed thread count: load, reduction (RustRed's math
// span), transport and the rest of evaluate (frontend, finalizer), and the
// cache save made after the clock.
#let breakdown(cell, threads) = {
  let field(name) = seconds(cell, "rustred", threads, "cold", field: name)
  let (evaluate, reduction) = (field("evaluate_seconds"), field("reduction_seconds"))
  let rest = if is-num(evaluate) and is-num(reduction) { calc.max(evaluate - reduction, 0) } else { none }
  [#threads: #fmt(field("load_seconds")) / #fmt(reduction) / #fmt(rest) / #fmt(field("save_seconds"))]
}
#let target-table = {
  let rows = ()
  for cell in four-loop-cells {
    let misses = ((1, target-1(cell)), (16, target-16(cell))).filter(pair => pair.at(1) == false).map(pair => pair.at(0))
    let r1 = cold-ratio(cell, 1)
    let r16 = cold-ratio(cell, 16)
    rows.push(text(size: 7pt)[#cell.replace("_", " ")])
    rows.push(text(size: 7pt)[#if r1 == none [—] else [×#fmt-ratio(r1)] #verdict(target-1(cell))])
    rows.push(text(size: 7pt)[#if r16 == none [—] else [×#fmt-ratio(r16)] #verdict(target-16(cell))])
    rows.push(text(size: 7pt)[#fmt(cell-load(cell, 1)) / #fmt(one-off(cell, 1)) #verdict(target-load(cell))])
    rows.push(text(size: 6.5pt)[#misses.map(threads => breakdown(cell, threads)).join(linebreak())])
  }
  table(
    columns: (auto, auto, auto, auto, 1fr),
    align: (x, y) => if x == 0 { left } else { right },
    inset: (x: 3pt, y: 1.6pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy },
    header[cell], header[1 thread (≤ ×#par-tolerance)], header[16 threads (< ×1)], header[load / cold − ws, 1 thread],
    header[misses: threads: load / reduction / rest / save (s, RustRed cold)],
    ..rows,
  )
}

// 1L-3L under the scope 1-4 program (one cell per loop count, 1 thread): the
// measured RustRed cells against the cited ones of the built-in 1-3 program.
#let cross-cells = measured.filter(m => m.backend == "rustred" and m.loops < 4 and m.threads == 1).map(m => m.cell).dedup()
#let cross-table = {
  let rows = ()
  for cell in cross-cells {
    let own(state, field) = median-of(pick(measured, cell, "rustred", 1, state, headline), field)
    let old(state, field) = median-of(pick(cited, cell, "rustred", 1, state, headline), field)
    rows.push(text(size: 7pt)[#cell.replace("_", " ")])
    for value in (
      own("cold", "evaluate_seconds"), own("warm-same", "evaluate_seconds"), old("cold", "evaluate_seconds"),
      own("cold", "reduction_seconds"), old("cold", "reduction_seconds"),
      own("cold", "load_seconds"), old("cold", "load_seconds"),
    ) { rows.push(text(size: 7pt)[#fmt(value)]) }
  }
  table(
    columns: (auto, ..((1fr,) * 7)),
    align: (x, y) => if x == 0 { left } else { right },
    inset: (x: 3pt, y: 1.6pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y < 2 { navy },
    table.cell(rowspan: 2, header[cell]), table.cell(colspan: 3, header[evaluate]),
    table.cell(colspan: 2, header[reduction only]), table.cell(colspan: 2, header[load]),
    header[1-4 cold], header[1-4 ws], header[1-3 cold], header[1-4], header[1-3], header[1-4], header[1-3],
    ..rows,
  )
}

// ------------------------------------------------------------------ title --

#let library = library-of(protocol)
#let library-name = if library == none { "p-library-1-3" } else { library.directory }
#let library-origin = if library == none { "built-in" } else { library.origin }

#align(center)[
  #text(19pt, weight: "bold", fill: navy)[vakint p-integrals: RustRed vs Forcer]
  #v(-0.5em)
  #text(9pt, fill: luma(100))[
    #if date == none [no selected run] else { date.display("[year]-[month]-[day]") }
    #if selected != none [· #protocol.host · gammaloop #protocol.gammaloop_sha.slice(0, 9)#if protocol.source_dirty [ (dirty)] · rustred #protocol.rustred_sha.slice(0, 8) · #library-name (#library-origin)]
  ]
]

#block(fill: pale, inset: 7pt, radius: 3pt, width: 100%)[
  #set text(8.6pt)
  *Reading guide.* Each cell is the sum over all tops of that loop count at the given dots, times a
  numerator of that rank. Bold: RustRed cold (seconds, median). F: Forcer cold (Fw: warm, where
  only that was measured). ×: the ratio RustRed cold over Forcer; green is RustRed faster.
  ws, wc: RustRed warm on its own cache, and warm on one cache shared along the cells of that loop
  count. ✓: result exactly equal to Forcer's; n.c.: not exact-checked (cost; `p_exact.json`). n: integrals in the input. †: fewer repetitions than
  its run planned.
  #if cited-run != none [*Sources:* 4L RustRed measured on the scope 1-4 library (#library-name, #library-origin);
    every Forcer row and RustRed at 1L-3L cited from benchmark v1 (run #raw(cited-run.id), gammaloop
    #cited-protocol.gammaloop_sha.slice(0, 9), built-in 1-3 library). Identity caveat: Forcer, FORM and vakint's
    Forcer path are unchanged since v1, and the built-in 1-3 route is unchanged; Table 6 measures 1L-3L under the
    1-4 program (one cell per loop count; the cross-program cost in full: M3-T19 S2).]
  *Caveat:* indicative timings on a shared, loaded host.
  #if cited-run != none [Benchmark v1 (cited rows):] else [Benchmark v1:] median foreign share of the measuring CPUs 1L n/a (before the sampler), 2L 0.82, 3L 0.12, 4L 0.52; v1's RustRed cold and wc include the after-evaluation cache save (fsync on the ZFS work dir): 122 of 624 saving children stalled 2-18 s, none of the non-saving children nor Forcer at 1L-3L.
  #if measured.any(m => "save_seconds" in m) [Measured rows save the cache after the clock (save seconds in Table 5 and the sub-lines).]
]

#block(stroke: 0.6pt + luma(200), inset: 6pt, radius: 3pt, width: 100%)[
  #set text(8.4pt)
  *Performance gates* (recorded, not blocking; the 4L target blocks in M5). #h(1em)
  4L per cell (Table 5): on par with Forcer at 1 thread (≤ ×#par-tolerance): #tally(target-1)\; better at 16 threads:
  #tally(target-16)\; cold start dominated by loading: #tally(target-load). #h(1em)
  Load of a fresh process: #let loads = load-seconds(measured); #if loads.len() == 0 [pending] else [#fmt(median(loads)) s (median, #library-name)]. #h(1em)
  G10b of v1, 3L rank ≥ 4, 16 threads, RustRed cold ≤ Forcer: #range-of(g10b, 9).
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
#figure(timing-table(1), caption: [1 thread (`n_cores` 1, `form`); RustRed arithmetic #raw(headline).], supplement: [Table])
#figure(timing-table(16), caption: [16 threads (`n_cores` 16, `tform -w16`); RustRed arithmetic #raw(headline).], supplement: [Table])

// ------------------------------------------------- arithmetic and reuse --

// Exact against modular_polynomial, RustRed cold, on three grid cells and the
// physics cell per loop count.
#let exact-cells = (1, 2, 3).map(loops => {
  let names = ("d0_r0", "d1_r4", "d2_r8").map(cell => "L" + str(loops) + "_" + cell)
  if loops >= 2 { names.push("physics_L" + str(loops)) }
  names.map(name => (loops, name))
}).flatten().chunks(2)
#let gain(a, b) = { let r = ratio(a, b); if r == none { [—] } else { [×#fmt-ratio(r)] } }
#let exact-table = {
  let rows = ()
  for (loops, cell) in exact-cells {
    rows.push(text(size: 7pt)[#cell.replace("_", " ")])
    for threads in (1, 16) {
      let ex = seconds(cell, "rustred", threads, "cold", arithmetic: "exact", field: "reduction_seconds")
      let mp = seconds(cell, "rustred", threads, "cold", field: "reduction_seconds")
      rows.push(text(size: 7pt)[#fmt(ex)])
      rows.push(text(size: 7pt)[#fmt(mp)])
      rows.push(text(size: 7pt, weight: "semibold")[#gain(ex, mp)])
    }
  }
  table(
    columns: (auto, ..((1fr,) * 6)),
    align: (x, y) => if x == 0 { left } else { right },
    inset: (x: 3pt, y: 1.6pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y < 2 { navy },
    table.cell(rowspan: 2, header[cell]), table.cell(colspan: 3, header[1 thread]), table.cell(colspan: 3, header[16 threads]),
    ..((header[exact], header[mod.], header[gain]) * 2),
    ..rows,
  )
}

// Measured vakint-level reuse per loop count: the sums of the medians of a
// warm state over those of RustRed cold, on the cells where both exist.
#let reuse-ratio(loops, threads, state) = {
  let cells = rustred-rows(loops).filter(m => m.loops == loops and m.backend == "rustred").map(m => m.cell).dedup()
  let pairs = cells
    .map(cell => (seconds(cell, "rustred", threads, state), seconds(cell, "rustred", threads, "cold")))
    .filter(pair => pair.all(is-num))
  if pairs.len() == 0 { [—] } else {
    let r = pairs.map(pair => pair.at(0)).sum() / pairs.map(pair => pair.at(1)).sum()
    [×#fmt-ratio(r)]
  }
}
// Cache counters of the RustRed children (headline arithmetic) of one loop count and state.
#let children(loops, state) = rustred-rows(loops).filter(m => (
  m.loops == loops and m.backend == "rustred" and m.state == state and m.status == "ok"
    and arithmetic-of(m) == headline and m.at("cache", default: none) != none
))
// Warm-same replays its cold cache: no miss and no insert in any child.
#let clean-replays(loops) = {
  let all = children(loops, "warm-same")
  if all.len() == 0 { [—] } else {
    [#all.filter(m => m.cache.misses == 0 and m.cache.inserts == 0).len() / #all.len()]
  }
}
// Warm-shared hits and misses summed over the chain of one loop count, the
// median over repetitions, at 16 threads.
#let chain-sum(loops, field) = {
  let all = children(loops, "warm-shared").filter(m => m.threads == 16)
  let sums = all.map(m => m.repetition).dedup().map(rep => all.filter(m => m.repetition == rep).map(m => m.cache.at(field)).sum())
  if sums.len() == 0 { [—] } else { [#calc.round(median(sums))] }
}
#let reuse-table = table(
  columns: (auto, ..((1fr,) * 6)),
  align: (x, y) => if x == 0 { left } else { right },
  inset: (x: 3pt, y: 1.6pt),
  stroke: 0.4pt + luma(215),
  fill: (x, y) => if y < 2 { navy },
  table.cell(rowspan: 2, header[L]), table.cell(colspan: 2, header[1 thread]), table.cell(colspan: 2, header[16 threads]),
  table.cell(rowspan: 2, header[ws clean]), table.cell(rowspan: 2, header[wc hits / misses]),
  ..((header[ws/cold], header[wc/cold]) * 2),
  ..(1, 2, 3, 4).map(loops => (
    text(size: 7pt)[#loops L],
    ..(1, 16).map(threads => ("warm-same", "warm-shared").map(state => text(size: 7pt)[#reuse-ratio(loops, threads, state)])).flatten(),
    text(size: 7pt)[#clean-replays(loops)],
    text(size: 7pt)[#chain-sum(loops, "hits") / #chain-sum(loops, "misses")],
  )).flatten(),
)

#let v1 = if cited-run == none [] else [ (1L-3L: benchmark v1)]
#block(breakable: false, grid(
  columns: (1fr, 1fr),
  column-gutter: 8pt,
  figure(exact-table, supplement: [Table], caption: [RustRed cold, arithmetic `exact` against `modular_polynomial` (mod.), reduction-only seconds (median); gain = exact / mod.#v1]),
  figure(reuse-table, supplement: [Table], caption: [Reuse, `modular_polynomial`: warm time over cold time per loop count (sums over the cells). ws clean: warm-same children with 0 misses and 0 inserts. wc: the chain's cache counters at 16 threads.#v1 #if cited-run != none [The 1L-3L ratios are dominated by the cache-save stalls of v1's cold runs (caveat).]]),
))

#figure(target-table, supplement: [Table], caption: [4L speed target per cell (G10b; recorded, not blocking): RustRed cold over Forcer#if cited-run != none [ (cited)] at 1 and 16 threads; load of a fresh process against the one-off cost of a cold evaluation (cold − ws); for a missed thread count the breakdown of the cold run. Inside the reduction about a third is step generation and two thirds the ordered weight elimination (M3-T19 S2, 4 threads).])

#figure(cross-table, supplement: [Table], caption: [1L-3L under the scope 1-4 program (#library-name), 1 thread, against the built-in 1-3 program of the cited rows; seconds (median). The 1-4 load is the price of one program for every loop count (M3 decision D-M3-1).])

// ------------------------------------------------------------ methodology --

// A sum over a part of the measured cells would read as the whole: a cell
// without the sub-line (no math span measured) leaves a dash.
#let per-loops(backend, state, field) = (1, 2, 3, 4).map(loops => {
  let source = if backend == "rustred" { rustred-rows(loops) } else { measurements }
  let cells = (
    source
      .filter(m => (
        m.loops == loops and m.backend == backend and m.threads == 1 and m.state == state and m.status == "ok"
          and (backend != "rustred" or arithmetic-of(m) == headline)
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
    [RustRed cache save (after the clock)], ..per-loops("rustred", "cold", "save_seconds"),
    [Forcer FORM subprocess], ..per-loops("forcer", "cold", "reduction_seconds"),
    [Raw FORM reference], ..per-loops("forcer", "raw", "evaluate_seconds"),
  )

  *Method.* Timed: the wall time of `Vakint::evaluate` in Retain mode (exact reduction onto masters,
  same frontend and finalizer for both backends) in one fresh process per measurement, with the
  input parsed, the engine and library loaded and the 1L bubble with four dots (in no cell: vakint
  memoises Forcer's FORM runs per process) evaluated once beforehand, by either backend, on directories of its own; repetitions:
  #protocol.at("repetitions", default: 3), median. RustRed cold: empty ordered-cache directory, shipped caches off#if measured.any(m => "save_seconds" in m) [; the cache is saved after the clock]. ws: a fresh
  process on a read-only copy of the cache of its cold run. wc: per loop count, one chain over the
  cells in table order sharing one cache from empty, each cell measured once in the chain. Forcer
  cold: private empty materialisation directory; Forcer warm: materialised. Forcer has no result
  cache. The raw reference runs the generated FORM programs directly, one process per top. OS pages
  are not flushed. #if protocol.at("cap_seconds", default: none) != none [Cap per measurement: #over-cap(protocol.cap_seconds).]
  #if cited-run != none [Cited 4L Forcer rows: 16 threads, warm, one repetition (benchmark v1).]
  Physics row: massless-quark (L-1)-rung gluon self-energy, Feynman gauge, transverse
  projector, external momentum kept.

  *Provenance.* #if selected == none [No selected run.] else [
    Run #raw(selected.id), profile #protocol.build_profile#if protocol.source_dirty [, dirty tree],
    #raw(protocol.form_version.split("\n").first()), Forcer #protocol.forcer_commit.slice(0, 7),
    #protocol.cpu_model on #protocol.host,
    CPUs #protocol.cpus.pairs().map(pair => pair.at(0) + " thread: " + pair.at(1)).join("; "),
    binding #if protocol.numactl == none [none] else { raw(protocol.numactl) },
    RustRed arithmetic #raw(headline) in Tables 1, 2, 4-6, `exact` in Table 3 (each child records the arithmetic it resolved),
    `RUSTRED_ORDERED_ADJOINT_TAIL` #if protocol.at("adjoint_tail_env", default: none) == none [unset] else [set],
    #if library == none [no official external data loaded] else if library.origin == "built-in" [built-in p library, no official external data loaded] else [
      p library #raw(library.directory) listed by the external data at #raw(str(library.root))#if library.at("root_head", default: none) != none [ (commit #raw(library.root_head.slice(0, 8)))], shipped caches off],
    #raw(str(measured.len())) measurements in `data/p_results.json`#if cited-run != none [; cited: run #raw(cited-run.id), gammaloop #raw(cited-protocol.gammaloop_sha.slice(0, 9)), rustred #raw(cited-protocol.rustred_sha.slice(0, 8)), #raw(str(cited.len())) measurements].
  ]
]
