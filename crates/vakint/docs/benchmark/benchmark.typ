// vakint p-integrals: RustRed vs Forcer, the clean benchmark (bench-clean v1).
// Every number is read from data/p_clean.json (one harness, one build, 1L-4L,
// 1 and 16 threads; artefact bench-clean.v1 copied verbatim). data/p_results.json
// (benchmark M3, gammaloop fd8b44641) is read only for its pinned outputs, and
// data/p_exact.json (M3 exact-agreement tier) for the agreement marks.
// Build: `just vakint-p-benchmark-pdf` in the dev shell (Typst 0.15.1).

#let data = json("data/p_clean.json")
#let pinned = json("data/p_results.json")
#let exact = json("data/p_exact.json")
#assert(data.schema == "bench-clean.v1", message: "data/p_clean.json is not bench-clean.v1")
#assert(pinned.schema == "vakint.p-benchmark.v1", message: "data/p_results.json is not vakint.p-benchmark.v1")
#assert(exact.schema == "vakint.p-exact.v1", message: "data/p_exact.json is not vakint.p-exact.v1")

// The measured build and the pushed tip it stands for (runtime equivalence:
// section "Measured build and pushed tip"; the evidence is listed there).
#let measured-gl = data.build.gammaloop_shas
#let measured-rr = data.build.rustred_shas
#assert(measured-gl.len() == 1 and measured-rr.len() == 1, message: "p_clean.json mixes builds")
#let measured-gl = measured-gl.first()
#let measured-rr = measured-rr.first()
// data/p_receipt.json records both builds and the evidence of their equivalence.
#let receipt = json("data/p_receipt.json")
#assert(receipt.schema == "vakint.p-benchmark-receipt.v1", message: "data/p_receipt.json is not vakint.p-benchmark-receipt.v1")
#assert(
  receipt.measured.gammaloop == measured-gl and receipt.measured.rustred == measured-rr
    and (receipt.measured.binary_sha256,) == data.build.binary_sha256
    and receipt.measured.p_library_1_4_sha256sums_sha256 == data.packages.at("p-library-1-4").sha256sums_sha256,
  message: "data/p_receipt.json names another measured build than data/p_clean.json",
)
#assert(receipt.relation in ("runtime-equivalent", "identical", "measured-on-parent"),
  message: "the receipt attests neither identity, runtime equivalence nor a stated parent build")
#let same-build = receipt.relation == "identical"
#let parent-build = receipt.relation == "measured-on-parent"
#assert(
  not same-build or (receipt.published_with.gammaloop == measured-gl and receipt.published_with.rustred == measured-rr),
  message: "receipt says identical but names another pushed build",
)
#let tip-gl = receipt.published_with.gammaloop
#let tip-rr = receipt.published_with.rustred
#let tip-xd = receipt.published_with.vakint_external_data
#for row in data.rows {
  assert(
    row.gammaloop_sha == (measured-gl,) and row.rustred_sha == (measured-rr,),
    message: "row " + row.cell + " was not measured on the one build",
  )
}

#let stamp = data.date_range.last()
#let date = datetime(year: int(stamp.slice(0, 4)), month: int(stamp.slice(5, 7)), day: int(stamp.slice(8, 10)))

#set document(title: "vakint p-integrals: RustRed vs Forcer (clean benchmark)", author: "Valentin Hirschi", date: date)
#set page(paper: "a4", margin: (x: 1.3cm, y: 1.4cm), numbering: "1 / 1")
#set text(font: "Libertinus Serif", size: 10pt, lang: "en")
#set par(justify: true, leading: 0.55em)
#show raw: set text(font: "DejaVu Sans Mono", size: 7.5pt)
#show figure.where(kind: table): set figure.caption(position: top)
#show figure.caption: set text(size: 9pt)

#let navy = rgb("#1f3b63")
#let pale = rgb("#eef3f9")
#let grey = luma(232)
#let red = rgb("#b2182b")
#let green = rgb("#1b7a3e")

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
#let fmt-ratio(r) = if not is-num(r) { "—" } else if r >= 10 { fixed(r, 1) } else if r >= 0.01 { fixed(r, 2) } else {
  "1/" + str(int(calc.round(1 / r)))
}
// Peak RSS in MiB (ru_maxrss); printed in MiB or GiB.
#let fmt-mb(mb) = if not is-num(mb) { "—" } else if mb >= 1000 { fixed(mb / 1024, 1) + " GiB" } else {
  str(int(calc.round(mb))) + " MiB"
}
#let median(values) = {
  let sorted = values.sorted()
  let n = sorted.len()
  if n == 0 { none } else { (sorted.at(int((n - 1) / 2)) + sorted.at(int(n / 2))) / 2 }
}
#let ratio(a, b) = if is-num(a) and is-num(b) and a > 0 and b > 0 { a / b } else { none }

// ------------------------------------------------------------- the rows --

#let key(cell, backend, threads, state) = cell + "|" + backend + "|" + str(threads) + "|" + state
#let by-key = {
  let d = (:)
  for row in data.rows { d.insert(key(row.cell, row.backend, row.threads, row.state), row) }
  d
}
#assert(by-key.len() == data.rows.len(), message: "p_clean.json has duplicate rows")
#let row(cell, backend, threads, state) = by-key.at(key(cell, backend, threads, state), default: none)
#let val(cell, backend, threads, state, field) = {
  let r = row(cell, backend, threads, state)
  if r == none { none } else { r.at(field, default: none) }
}
#let ranks = (0, 2, 4, 6, 8)
#let grid-cells(loops) = (0, 1, 2).map(dots => ranks.map(rank => "L" + str(loops) + "_d" + str(dots) + "_r" + str(rank))).flatten()
#let cells-of(loops) = grid-cells(loops) + if loops >= 2 { ("physics_L" + str(loops),) } else { () }
#let all-cells = (1, 2, 3, 4).map(cells-of).flatten()

// The ratios of the summary (RustRed over Forcer, medians over repetitions):
// cold evaluate; whole cold process without RustRed's cache save over Forcer's
// whole cold process; peak RSS: RustRed's measuring process against Forcer's
// process plus its FORM subprocess (sum of the two peaks).
#let ev-ratio(cell, t) = ratio(val(cell, "rustred", t, "cold", "evaluate"), val(cell, "forcer", t, "cold", "evaluate"))
#let proc-ratio(cell, t) = ratio(
  val(cell, "rustred", t, "cold", "process_excl_save"),
  val(cell, "forcer", t, "cold", "process"),
)
#let rss-rr(cell, t) = val(cell, "rustred", t, "cold", "peak_rss_mb")
#let rss-f(cell, t) = {
  let (own, kids) = (val(cell, "forcer", t, "cold", "peak_rss_mb"), val(cell, "forcer", t, "cold", "children_peak_rss_mb"))
  if is-num(own) { own + if is-num(kids) { kids } else { 0 } } else { none }
}
#let rss-ratio(cell, t) = ratio(rss-rr(cell, t), rss-f(cell, t))
#let slower(t) = all-cells.filter(cell => { let r = ev-ratio(cell, t); r != none and r > 1 })
#let worst(cells, f) = {
  let found = cells.filter(c => f(c) != none)
  if found.len() == 0 { none } else { found.sorted(key: c => f(c)).last() }
}

// --------------------------------------------- outputs against the pins --

// The pinned outputs of benchmark M3 (and v1), per cell and backend.
#let pins = {
  let d = (:)
  for run in pinned.runs {
    for m in run.measurements {
      let k = m.cell + "|" + m.backend
      let h = m.at("output_sha256", default: none)
      if h != none and m.status == "ok" { d.insert(k, d.at(k, default: ()) + (h,)) }
    }
  }
  d
}
// Re-checked here: a row the data set calls identical has its output among the
// committed pins, a row without a pin has none committed.
#let identical-rows = data.rows.filter(r => r.output_check == "identical")
#let unpinned-rows = data.rows.filter(r => r.output_check == "no pinned reference")
#assert(identical-rows.len() + unpinned-rows.len() == data.rows.len(), message: "an output differs from its pin")
#assert(data.output_failures.len() == 0, message: "p_clean.json records output failures")
#for r in identical-rows {
  assert(r.output_sha256 in pins.at(r.cell + "|" + r.backend, default: ()), message: "pin not in p_results.json: " + r.cell)
}
#for r in unpinned-rows {
  assert(pins.at(r.cell + "|" + r.backend, default: ()).len() == 0, message: "unpinned row has a pin: " + r.cell)
}
#let unpinned-cells = unpinned-rows.map(r => r.cell).dedup()
// Within one cell and backend every row (state, thread count, repetition median) gave one output.
#let one-output = all-cells.all(cell => ("rustred", "forcer").all(b => (
  data.rows.filter(r => r.cell == cell and r.backend == b).map(r => r.output_sha256).dedup().len() == 1
)))

// Exact-agreement marks of the M3 tier (RustRed = Forcer, exactly), carried
// by the byte-identical outputs; n.c.: the tier did not check the cell (cost).
#let parse-cell(cell) = if cell.starts-with("physics_L") { (int(cell.slice(9)), none, none) } else {
  (int(cell.slice(1, 2)), int(cell.slice(4, 5)), int(cell.slice(7)))
}
#let same-cell(loops, dots, rank) = c => if dots == none {
  c.at("physics", default: none) not in (none, false) and c.L == loops
} else {
  c.L == loops and c.at("dots", default: none) == dots and c.at("rank", default: none) == rank
}
#let exact-mark(cell) = {
  let (loops, dots, rank) = parse-cell(cell)
  let found = exact.cells.find(same-cell(loops, dots, rank))
  let unchecked = exact.at("not_exact_checked", default: ()).any(same-cell(loops, dots, rank))
  if found == none and unchecked { text(fill: luma(100), size: 6pt)[#h(2pt)n.c.] } else if found == none { [] } else if found.agree {
    text(fill: green, size: 7pt)[#h(2pt)✓]
  } else { text(fill: red, size: 7pt, weight: "bold")[#h(2pt)✗] }
}

// ----------------------------------------------------- foreign load (16t) --

#let job-hygiene(r) = {
  let found = r.jobs.map(j => data.hygiene.at(j, default: none)).filter(h => h != none)
  if found.len() == 0 { none } else { found.first() }
}
#let hygiene-values(prefix, field) = data.hygiene.pairs().filter(p => p.at(0).starts-with(prefix)).map(p => p.at(1).at(field)).filter(is-num)
#let foreign16 = hygiene-values("t16_", "median_total_foreign_cores")
#let foreign16-max = hygiene-values("t16_", "max_total_foreign_cores")
#let foreign1 = hygiene-values("t1_", "median_total_foreign_cores")
#let foreign-of(cell, backend, t) = {
  let r = row(cell, backend, t, "cold")
  let h = if r == none { none } else { job-hygiene(r) }
  if h == none { none } else { h.median_total_foreign_cores }
}

// --------------------------------------------------------------- colours --

// Diverging log2 palette: green (RustRed faster), white (parity), red (slower).
#let palette = gradient.linear(rgb("#5aae61"), rgb("#d9f0d3"), rgb("#ffffff"), rgb("#fddbc7"), rgb("#d6604d"))
#let ratio-fill(r) = if r == none { none } else { palette.sample(calc.clamp(0.5 + calc.log(r, base: 2) / 4, 0, 1) * 100%) }
#let header(body) = text(font: "Roboto", size: 7.5pt, weight: "bold", fill: white, body)

// A dagger on a median over fewer repetitions than the three planned.
#let partial(cell, backend, t, state) = {
  let r = row(cell, backend, t, state)
  if r != none and r.n_reps < 3 { "†" } else { "" }
}
#let secs(cell, backend, t, state) = fmt(val(cell, backend, t, state, "evaluate")) + partial(cell, backend, t, state)

// One table cell: RustRed cold, Forcer cold, the ratio, the warm states, the
// peak-RSS ratio and the exact mark.
#let cmp-cell(cell, t, physics: false) = {
  let r = ev-ratio(cell, t)
  let (loops, _, _) = parse-cell(cell)
  table.cell(fill: ratio-fill(r), {
    set par(leading: 0.3em, justify: false)
    if physics { text(size: 6.5pt, weight: "bold")[#loops L #h(2pt)] }
    text(weight: "bold", size: 8.2pt, secs(cell, "rustred", t, "cold"))
    exact-mark(cell)
    linebreak()
    text(size: 6.5pt, fill: luma(40))[F #secs(cell, "forcer", t, "cold")]
    text(size: 6.8pt, weight: "semibold")[#h(2pt)×#fmt-ratio(r)]
    linebreak()
    text(size: 5.6pt, fill: luma(45))[ws #secs(cell, "rustred", t, "warm-same") #h(1pt) wp #secs(cell, "rustred", t, "warm-process") #h(1pt) Fw #secs(cell, "forcer", t, "warm")]
    linebreak()
    let m = rss-ratio(cell, t)
    text(size: 5.6pt, fill: if m != none and m > 1 { red } else { luma(45) })[RSS #fmt-mb(rss-rr(cell, t)) ×#fmt-ratio(m)]
  })
}

#let timing-table(t) = {
  let rows = ()
  for loops in (1, 2, 3, 4) {
    for dots in (0, 1, 2) {
      if dots == 0 { rows.push(table.cell(rowspan: 3, align: horizon, text(weight: "bold")[#loops L])) }
      rows.push(text(size: 8pt)[#dots])
      for rank in ranks { rows.push(cmp-cell("L" + str(loops) + "_d" + str(dots) + "_r" + str(rank), t)) }
    }
  }
  rows.push(table.cell(colspan: 2, align: horizon, text(size: 8pt, weight: "bold")[phys.]))
  for loops in (2, 3, 4) { rows.push(cmp-cell("physics_L" + str(loops), t, physics: true)) }
  rows.push(table.cell(colspan: 2, align: horizon + center, fill: grey, text(size: 7pt, fill: luma(90))[5L: not part of this benchmark]))
  table(
    columns: (1.05cm, 0.75cm, ..((1fr,) * 5)),
    align: (x, y) => if x < 2 { center + horizon } else { left + horizon },
    inset: (x: 3pt, y: 2.2pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y == 0 { navy },
    header[L], header[dots], ..ranks.map(rank => header[rank #rank]),
    ..rows,
  )
}

// ---------------------------------------------------------------- summary --

#let summary-table = {
  let rows = ()
  for loops in (1, 2, 3, 4) {
    for t in (1, 16) {
      let cells = cells-of(loops)
      let stat(f) = {
        let values = cells.map(c => f(c, t)).filter(is-num)
        let w = worst(cells, c => f(c, t))
        (median(values), if w == none { none } else { f(w, t) }, w)
      }
      let (ev, pr, rs) = (stat(ev-ratio), stat(proc-ratio), stat(rss-ratio))
      let late = cells.filter(c => { let r = ev-ratio(c, t); r != none and r > 1 })
      rows += (
        text(size: 7pt)[#loops L], text(size: 7pt)[#t], text(size: 7pt)[#cells.len()],
        text(size: 7pt)[×#fmt-ratio(ev.at(0))], text(size: 7pt)[×#fmt-ratio(ev.at(1)) #text(size: 6pt)[(#ev.at(2).replace("_", " "))]],
        text(size: 7pt)[×#fmt-ratio(pr.at(0))], text(size: 7pt)[×#fmt-ratio(pr.at(1))],
        text(size: 7pt)[×#fmt-ratio(rs.at(0))], text(size: 7pt)[×#fmt-ratio(rs.at(1)) #text(size: 6pt)[(#rs.at(2).replace("_", " "))]],
        table.cell(fill: if late.len() > 0 { rgb("#fddbc7") }, text(size: 7pt, weight: if late.len() > 0 { "bold" } else { "regular" })[#if late.len() == 0 [none] else { late.map(c => c.replace("_", " ")).join(", ") }]),
      )
    }
  }
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, 1fr),
    align: (x, y) => if x < 3 { center } else { left },
    inset: (x: 3pt, y: 1.8pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y < 2 { navy },
    table.cell(rowspan: 2, header[L]), table.cell(rowspan: 2, header[thr.]), table.cell(rowspan: 2, header[cells]),
    table.cell(colspan: 2, header[cold evaluate]), table.cell(colspan: 2, header[cold process]),
    table.cell(colspan: 2, header[peak RSS]), table.cell(rowspan: 2, header[slower than Forcer]),
    header[median], header[worst], header[median], header[worst], header[median], header[worst],
    ..rows,
  )
}

// 4L in detail: the cold runs of both backends, per thread count.
#let detail-table = {
  let rows = ()
  for cell in cells-of(4) {
    for t in (1, 16) {
      let rr(field) = val(cell, "rustred", t, "cold", field)
      let ff(field) = val(cell, "forcer", t, "cold", field)
      let r = ev-ratio(cell, t)
      let fc = foreign-of(cell, "rustred", t)
      let ffc = foreign-of(cell, "forcer", t)
      rows += (
        text(size: 6.6pt)[#cell.replace("_", " ")], text(size: 6.6pt)[#t],
        text(size: 6.6pt, weight: "bold")[#fmt(rr("evaluate"))#partial(cell, "rustred", t, "cold")], text(size: 6.6pt)[#fmt(rr("reduction"))],
        text(size: 6.6pt)[#fmt(rr("load"))], text(size: 6.6pt)[#fmt(rr("save"))], text(size: 6.6pt)[#fmt(rr("process_excl_save"))],
        text(size: 6.6pt, weight: "bold")[#fmt(ff("evaluate"))], text(size: 6.6pt)[#fmt(ff("process"))],
        table.cell(fill: ratio-fill(r), text(size: 6.6pt, weight: "bold")[×#fmt-ratio(r)]),
        table.cell(fill: ratio-fill(proc-ratio(cell, t)), text(size: 6.6pt)[×#fmt-ratio(proc-ratio(cell, t))]),
        text(size: 6.6pt)[#fmt-mb(rss-rr(cell, t))], text(size: 6.6pt)[#fmt-mb(rss-f(cell, t))],
        text(size: 6.6pt, fill: if rss-ratio(cell, t) > 1 { red } else { black })[×#fmt-ratio(rss-ratio(cell, t))],
        text(size: 6.6pt)[#if t == 1 [—] else [#if is-num(fc) { fixed(fc, 1) } else { "<30 s" } / #if is-num(ffc) { fixed(ffc, 1) } else { "<30 s" }]],
      )
    }
  }
  table(
    columns: (auto, auto, ..((1fr,) * 13)),
    align: (x, y) => if x < 2 { left } else { right },
    inset: (x: 2.5pt, y: 1.5pt),
    stroke: 0.4pt + luma(215),
    fill: (x, y) => if y < 2 { navy },
    table.cell(rowspan: 2, header[cell]), table.cell(rowspan: 2, header[thr.]),
    table.cell(colspan: 5, header[RustRed cold (s)]), table.cell(colspan: 2, header[Forcer cold (s)]),
    table.cell(colspan: 2, header[RustRed / Forcer]), table.cell(colspan: 3, header[peak RSS]),
    table.cell(rowspan: 2, header[foreign cores RR / F]),
    header[evaluate], header[reduction], header[load], header[save], header[process − save],
    header[evaluate], header[process], header[evaluate], header[process],
    header[RustRed], header[Forcer], header[ratio],
    ..rows,
  )
}

// ------------------------------------------------------------------ title --

#let n-measurements = data.raw_measurements
#let pkg4 = data.packages.at("p-library-1-4")
#let pkg3 = data.packages.at("p-library-1-3")

#align(center)[
  #text(19pt, weight: "bold", fill: navy)[vakint p-integrals: RustRed vs Forcer]
  #v(-0.5em)
  #text(9pt, fill: luma(100))[
    clean benchmark #raw(data.id) (#raw(data.schema)) · #date.display("[year]-[month]-[day]") · #data.host.join(", ") ·
    document commit base gammaloop #raw(tip-gl.slice(0, 9)), rustred #raw(tip-rr.slice(0, 8)) ·
    measured build gammaloop #raw(measured-gl.slice(0, 9)), rustred #raw(measured-rr.slice(0, 8)) (#receipt.relation.replace("-", " "): last section)
  ]
]

#let worst4(t) = worst(cells-of(4), c => rss-ratio(c, t))
#let pairs3 = (1, 2, 3).map(cells-of).flatten().map(c => ((c, 1), (c, 16))).flatten().chunks(2).filter(p => rss-ratio(..p) != none)
#let low3 = pairs3.filter(p => rss-ratio(..p) < 1)
#let high3 = pairs3.filter(p => rss-ratio(..p) >= 1)
#let below4 = cells-of(4).filter(c => (1, 16).any(t => rss-ratio(c, t) <= 1))
#block(fill: pale, inset: 7pt, radius: 3pt, width: 100%)[
  #set text(8.6pt)
  *In short.* One harness, one build and one data set for every number: #data.rows.len() rows (cell × backend ×
  threads × state, medians of #n-measurements measurements in #data.runs.len() runs, #data.date_range.first().slice(11, 16)–#data.date_range.last().slice(11, 16)Z),
  1L–4L, RustRed (#receipt.at("package_label", default: "M5 combined package")) against Forcer, at 1 and 16 threads, cold and warm, time and peak RSS. Older
  benchmark rows are not mixed in.
  - *Time.* RustRed is slower than Forcer (cold evaluate) in #slower(1).len() of #all-cells.len() cells at 1 thread#if slower(1).len() > 0 [
    (#slower(1).map(c => [#raw(c) ×#fmt-ratio(ev-ratio(c, 1))]).join(", "))] and in #slower(16).len() at 16 threads#if slower(16).len() > 0 [
    (#slower(16).map(c => [#raw(c) ×#fmt-ratio(ev-ratio(c, 16))]).join(", "))]; in every other cell it is faster.
    From 1 to 16 threads the 4L cold evaluate changes by a median ×#fmt-ratio(median(cells-of(4).map(c => ratio(val(c, "rustred", 16, "cold", "evaluate"), val(c, "rustred", 1, "cold", "evaluate"))).filter(is-num)))
    for RustRed and ×#fmt-ratio(median(cells-of(4).map(c => ratio(val(c, "forcer", 16, "cold", "evaluate"), val(c, "forcer", 1, "cold", "evaluate"))).filter(is-num)))
    for Forcer (`tform -w16`).
  - *16-thread numbers are upper bounds.* The host was shared and loaded: on the 16 CPUs of a 16-thread job
    unpinned single-thread processes of other users took a median of #fmt(median(foreign16)) cores per job (median over
    #foreign16.len() jobs; up to #fmt(calc.max(..foreign16-max)) cores in single samples). Both backends were exposed alike;
    1-thread jobs ran with their SMT sibling held idle; other users' load on that sibling: median #fmt(median(foreign1)) cores.
  - *Memory.* At 1L–3L RustRed's peak RSS is below Forcer's (process plus FORM) in #low3.len() of #pairs3.len() cells × thread
    counts (median ×#fmt-ratio(median(pairs3.map(p => rss-ratio(..p))))#if high3.len() > 0 [; above: #high3.map(p => [#raw(p.at(0)) at #p.at(1) thread(s) ×#fmt-ratio(rss-ratio(..p))]).join(", ")]).
    At 4L it is above Forcer's in #(cells-of(4).len() - below4.len()) of #cells-of(4).len() cells#if below4.len() > 0 [ (below: #below4.map(c => [#raw(c) ×#fmt-ratio(rss-ratio(c, 1)) / ×#fmt-ratio(rss-ratio(c, 16))]).join(", "))]: median ×#fmt-ratio(median(cells-of(4).map(c => rss-ratio(c, 1)).filter(is-num))) at 1 thread and
    ×#fmt-ratio(median(cells-of(4).map(c => rss-ratio(c, 16)).filter(is-num))) at 16 threads, worst #raw(worst4(1))
    (#fmt-mb(rss-rr(worst4(1), 1)) against #fmt-mb(rss-f(worst4(1), 1)), ×#fmt-ratio(rss-ratio(worst4(1), 1))) and
    #raw(worst4(16)) at 16 threads (#fmt-mb(rss-rr(worst4(16), 16)) against #fmt-mb(rss-f(worst4(16), 16)), ×#fmt-ratio(rss-ratio(worst4(16), 16))).
  - *Results.* All #identical-rows.len() rows with a pinned reference are byte-identical to it (the outputs pinned by
    benchmark M3 in `data/p_results.json`, re-checked when this document is built); #unpinned-rows.len() RustRed rows
    (#unpinned-cells.map(raw).join(", ")) have no pinned reference (M3 capped them)#if one-output [; every cell gave one
    output per backend across states, thread counts and repetitions].
]

#figure(summary-table, supplement: [Table], caption: [RustRed over Forcer per loop count and thread count, medians over
  the cells (grid and physics cell) of the per-cell medians over repetitions. Cold evaluate: wall time of
  `Vakint::evaluate`. Cold process: the whole RustRed process (load, warm-up, evaluate, exit) without its cache save, over
  Forcer's whole process. Peak RSS: RustRed's process against Forcer's process plus its FORM subprocess. Slower: cold
  evaluate ratio > 1.])

#v(-0.3em)
#align(center, stack(
  dir: ltr,
  spacing: 4pt,
  text(size: 7pt)[RustRed faster ×1/4],
  box(width: 5cm, height: 7pt, fill: gradient.linear(rgb("#5aae61"), rgb("#d9f0d3"), rgb("#ffffff"), rgb("#fddbc7"), rgb("#d6604d"))),
  text(size: 7pt)[×4 Forcer faster],
))

#block(fill: pale, inset: 6pt, radius: 3pt, width: 100%)[
  #set text(8.2pt)
  *Reading guide* (Tables 2 and 3). Each cell sums all tops of that loop count at the given dots, times a numerator of
  that rank. Bold: RustRed cold (seconds, median of 3). F: Forcer cold; ×: RustRed cold over Forcer cold (colour); ws:
  RustRed warm, a fresh process on its own cold cache; wp: the warm repeat inside the cold process; Fw: Forcer warm
  (materialised). RSS: RustRed's cold peak RSS and its ratio to Forcer's (red: RustRed higher). ✓: result exactly equal to
  Forcer's (M3 exact tier, `p_exact.json`, carried by the byte-identical outputs); n.c.: not exact-checked (cost).
  †: two repetitions instead of three.
]

#figure(timing-table(1), caption: [1 thread (`n_cores` 1, `form`); RustRed arithmetic #raw(data.harness.arithmetic).], supplement: [Table])
#figure(timing-table(16), caption: [16 threads (`n_cores` 16, `tform -w16`) on 16 CPUs of one NUMA node; upper bounds (foreign load, Table 1 note); RustRed arithmetic #raw(data.harness.arithmetic).], supplement: [Table])

#figure(detail-table, supplement: [Table], caption: [4L in detail, cold runs, medians. RustRed: evaluate (the timed
  call), its reduction span, the load of engine and library before the clock, the cache save after the clock, and the
  whole process without that save. Forcer: evaluate and its whole process. Foreign cores: median load of other users'
  processes on the job's CPUs (16 threads only; 1-thread jobs had none to speak of).])

// ------------------------------------------------------------- equivalence --

#heading(level: 2, numbering: none)[Measured build and pushed tip]
#set text(8.6pt)
#if parent-build [
The rows were measured on the pushed tip of #receipt.measured.at("label", default: "the measured build"): one frozen harness binary
(sha256 #raw(data.build.binary_sha256.first().slice(0, 16))…) built with `just vakint-p-benchmark-build` from gammaloop
#raw(measured-gl.slice(0, 9)) (stamp = that commit, clean tree) pinning rustred #raw(measured-rr.slice(0, 8)); the 4L library
is `p_integral/p-library-1-4` of vakint-external-data #raw(tip-xd.slice(0, 8)) (#receipt.at("ext_data_tag", default: ""))
read through `--external-data` (#raw(pkg4.sha256sums_sha256.slice(0, 16))…), the 1L–3L library the built-in one (identity
ends #raw(pkg3.identity.split(",").last())). This document is committed on gammaloop #raw(tip-gl.slice(0, 9)) (rustred
#raw(tip-rr.slice(0, 8))). #receipt.evidence.pin_change *These timings are therefore not those of the commit that carries
this document*; #receipt.evidence.at("next_build", default: "")
] else if same-build [
The rows were measured with the pushed tip itself: one frozen harness binary (sha256 #raw(data.build.binary_sha256.first().slice(0, 16))…)
built with `just vakint-p-benchmark-build` from gammaloop #raw(measured-gl.slice(0, 9)) (stamp = that commit, clean tree)
pinning rustred #raw(measured-rr.slice(0, 8)); the 4L library is `p_integral/p-library-1-4` of vakint-external-data
#raw(tip-xd.slice(0, 8)) (#receipt.at("ext_data_tag", default: "")) read through `--external-data` (#raw(pkg4.sha256sums_sha256.slice(0, 16))…,
#receipt.evidence.at("p_library_1_4", default: "")), the 1L–3L library the built-in one (identity ends
#raw(pkg3.identity.split(",").last())). No equivalence argument is needed. #receipt.evidence.at("verification_of_the_tip", default: "")
] else [
The rows were measured with one frozen harness binary (sha256 #raw(data.build.binary_sha256.first().slice(0, 16))…) built
from gammaloop #raw(measured-gl.slice(0, 9)) (local branch, never pushed) pinning rustred #raw(measured-rr.slice(0, 8)).
They stand for the pushed tip — gammaloop #raw(tip-gl.slice(0, 9)), rustred #raw(tip-rr.slice(0, 8)), vakint-external-data
#raw(tip-xd.slice(0, 8)) — because the two builds are runtime-equivalent for every timed call (path by path in
`data/p_receipt.json`: #receipt.evidence.gammaloop_changed_paths.len() gammaloop paths,
#receipt.evidence.rustred_changed_source_paths.len() rustred source files):
+ *gammaloop* #raw(measured-gl.slice(0, 9))..#raw(tip-gl.slice(0, 9)) (4 commits): outside the rustred pin in
  `Cargo.toml`/`Cargo.lock`, only test code changes (`p_data.rs` inside `mod tests`, `tests/p_integral_rustred_tests.rs`,
  `tests/p_library_shipping_tests.rs`); the built-in library `crates/vakint/data/p` (the 1L–3L program) is the same git tree.
+ *rustred* #raw(measured-rr.slice(0, 8))..#raw(tip-rr.slice(0, 8)) (5 commits): test files and goldens; the package
  comparison `p_library_differences` (used by `rustred p package` and tests, never on an evaluation path); `rustfmt`
  of the scheme lane (token-identical up to added braces); `clippy` attributes and two needless borrows. One change is
  reachable at run time: the rational-constant fast path of `invert_symbolic_matrix` (`family/exact.rs`) is now also
  conditional on diagnostics being off and on every entry being a coefficient of the context. Ordinary execution has
  diagnostics off (the default) and only such coefficients, so it takes the same path; the added cost is one membership
  test per matrix entry.
+ *4L library*: `p_integral/p-library-1-4` of vakint-external-data #raw(tip-xd.slice(0, 8)) is the measured package
  (#raw(pkg4.sha256sums_sha256.slice(0, 16))…): all 472 files byte-identical by sha256. The 1L–3L library is the
  built-in one, identical in both builds (identity ends #raw(pkg3.identity.split(",").last())).
+ The pushed tip passed the vakint gate (1097/1097 tests, R\* 60/60) and the Forcer cross-check of all 22 tops of
  L4 d0_r8/d1_r8 through vakint; these verify results, not timings.
+ *Spot re-time*: #receipt.evidence.spot_retime.what. #receipt.evidence.spot_retime.cells.map(c => [#raw(c.cell) #c.state
  #fmt(c.tip) s against #fmt(c.measured) s (×#fmt-ratio(c.tip / c.measured))]).join("; "). Outputs
  #receipt.evidence.spot_retime.outputs.
]

// ------------------------------------------------------------ methodology --

*Method.* #data.harness.clock. Work directory #data.harness.work_dir; #data.harness.placement. States: RustRed cold
(empty ordered-cache directory, shipped caches off), warm-same (a fresh process on a read-only copy of its cold cache),
warm-process (in the cold process, after the save); Forcer cold (private empty materialisation directory) and warm
(materialised); Forcer has no result cache. Three repetitions per state (†: two), median; cap 3600 s per child (#{
  let capped = data.rows.filter(r => r.n_capped > 0 or r.n_skipped > 0).len()
  if capped == 0 [no row was capped or skipped] else [#capped rows capped or skipped]
}). OS pages are not flushed. Physics row: massless-quark (L−1)-rung gluon self-energy, Feynman gauge,
transverse projector, external momentum kept.

*Provenance.* `data/p_clean.json` = artefact #raw(data.id), written #raw(data.written_utc), copied verbatim; host #data.host.join(", ") (#data.cpu_model.join(", ")),
#raw(data.form.first().split(" (").first()), Forcer #raw(data.forcer_commit.first().slice(0, 7))\; RustRed at 1L–3L with the
built-in `p-library-1-3`, at 4L with `p-library-1-4` listed by external data. Outputs compared with the #data.reference.
Exact tier: `data/p_exact.json` (gammaloop #raw(exact.gammaloop_sha.slice(0, 9)), rustred #raw(exact.rustred_sha.slice(0, 8)),
#exact.cells.len() cells agree#if exact.cells.any(c => not c.agree) [, some do not]); of the 4L cells it marks n.c., the gate
cross-check above covers the 22 tops of L4 d0_r8/d1_r8. The M3 benchmark rows remain in `data/p_results.json` for their
pinned outputs only; none of their timings is used here.
