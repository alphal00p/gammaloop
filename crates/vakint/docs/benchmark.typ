// Vakint reduction benchmark. Generated data only: regenerate with
//   crates/vakint/docs/benchmark/make-benchmark.sh DATA_ROOT BENCHMARK_BINARY OUT_DIR [--install]
// Inputs (typst --input): results=<results.json>, graphs=<graphs.json>, tables=<tables.json>, draft=true|false.
// tables.json: benchmark/tools/precision_tables.py on the same data root (precision of masters and results).

#import "precision_tables.typ": precision-section
#let results = json(sys.inputs.at("results", default: "benchmark/data/suite.json"))
#let graphs = json(sys.inputs.at("graphs", default: "benchmark/data/graphs.json"))
#let tables = json(sys.inputs.at("tables", default: "benchmark/data/tables.json"))
#let draft = sys.inputs.at("draft", default: if results.at("draft", default: false) { "true" } else { "false" }) == "true"

#let navy = rgb("173f65")
#let ink = rgb("1d2733")
#let pale = rgb("eaf2fa")
#let muted = rgb("5b6876")
#let sans = ("DejaVu Sans", "Roboto", "Noto Sans", "Libertinus Serif")
#let ident = results.identity
#let host = ident.at("host", default: (:))
#let day = results.assembled_utc.slice(0, 10)
#let timing = results.at("timing", default: (:))

#set document(title: "Vakint reduction benchmark", author: "Valentin Hirschi")
#set page(
  paper: "a4",
  margin: (x: 16mm, top: 16mm, bottom: 16mm),
  footer: context align(right, text(size: 8pt, fill: rgb("687584"))[
    Vakint reduction benchmark · #counter(page).display("1")]),
  background: if draft {
    rotate(-32deg, text(size: 110pt, weight: "bold", fill: rgb(200, 40, 40, 28), font: sans)[DRAFT])
  },
)
#set text(font: "Libertinus Serif", size: 10.4pt, fill: ink, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 0.8em)
#show heading.where(level: 1): it => block(above: 1.1em, below: 0.6em, text(14pt, weight: "bold", fill: navy, it.body))
#show raw: set text(font: "DejaVu Sans Mono", size: 7.6pt)
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)
#show figure.where(kind: table): set block(breakable: true)

// ---------------------------------------------------------------- helpers --

#let is-num(x) = type(x) == int or type(x) == float
#let pl(n, one, many: none) = [#n #if n == 1 { one } else if many != none { many } else { one + "s" }]
#let fixed(x, d) = {
  let s = str(calc.round(x, digits: d))
  if d == 0 { return s }
  let parts = s.split(".")
  let frac = if parts.len() > 1 { parts.at(1) } else { "" }
  parts.at(0) + "." + frac + "0" * (d - frac.len())
}
#let sig(x) = {
  let a = calc.abs(x)
  if a >= 100 { fixed(x, 0) } else if a >= 10 { fixed(x, 1) } else if a >= 1 { fixed(x, 2) } else { fixed(x, 3) }
}
// Unit chosen for a duration, and the duration in that unit.
#let unit-of(s) = if s < 1e-3 { (1e6, "µs") } else if s < 1 { (1e3, "ms") } else if s < 600 { (1, "s") } else if s < 7200 { (1 / 60, "min") } else { (1 / 3600, "h") }
#let dur(s) = {
  if not is-num(s) { return [—] }
  let (f, u) = unit-of(s)
  [#sig(s * f)#h(1pt)#u]
}
// `x` in the unit of `ref` (the minimum next to its median), no unit.
#let in-unit(x, ref) = { let (f, _) = unit-of(ref); sig(x * f) }
#let ratio-text(r) = {
  if r < 0.095 { [×1/#str(int(calc.round(1 / r)))] }
  else if r < 10 { [×#fixed(r, if r < 1 { 2 } else { 1 })] }
  else { [×#str(int(calc.round(r)))] }
}
// Diverging palette of the ninjars benchmark, log2 scale: ×1/64 (green, RustRed faster) to ×64 (red).
#let span2 = 6
#let palette = gradient.linear(rgb("#5aae61"), rgb("#d9f0d3"), rgb("#ffffff"), rgb("#fddbc7"), rgb("#d6604d"))
#let ratio-fill(r) = palette.sample(calc.clamp(0.5 + calc.log(r, base: 2) / (2 * span2), 0, 1) * 100%)
#let na-fill = rgb("f3f5f8")
// Secondary text: grey on light fills, ink on the saturated ends of the palette.
#let lum(c) = { let (r, g, b, ..) = c.rgb().components(); 0.2126 * r / 100% + 0.7152 * g / 100% + 0.0722 * b / 100% }
#let second-ink(fill) = if fill != none and lum(fill) < 0.72 { ink } else { luma(85) }

#let cell-of(case, key) = case.cells.at(key, default: (status: "pending"))
#let med(c) = if c.status == "ok" { c.median } else { none }
// Small note next to a value: the minimum of the runs, or "1 run" when it is not a repetition.
#let spread(c, fill) = {
  if c.status != "ok" { return none }
  let note = if c.n == 1 [1 run] else [(#in-unit(c.min, c.median))]
  text(size: 5.6pt, weight: "regular", fill: second-ink(fill))[#h(2pt)#note]
}
// Multi-core cells (heavy five-loop reductions) name their core count on a line of their own.
#let cores-of(c) = if c.status == "ok" { c.at("n_cores", default: 1) } else { 1 }
#let cores-line(c, fill) = text(size: 6.1pt, fill: second-ink(fill))[on #cores-of(c) cores]
#let blank = text(size: 6.9pt)[#sym.space.nobreak]

// One comparison cell, three lines: RustRed median (min) / reference median (min) / ratio.
// `show-alone`: when the reference does not apply, still show the RustRed value (else only "n/a").
#let ibp-line(c) = {
  let m = c.at("reduction_math_median", default: none)
  if c.status == "ok" and m != none [IBP #dur(m)]
}
#let cmp(rr, ref, label, show-alone: true, alone-note: [no FORM backend], extra: none, ref-spread: true) = {
  let a = med(rr)
  let b = med(ref)
  if ref.status == "n/a" and not show-alone {
    return table.cell(fill: na-fill, text(size: 7pt, fill: luma(150))[n/a])
  }
  let fill = if a != none and b != none { ratio-fill(a / b) } else if ref.status == "n/a" { na-fill } else { none }
  let ink2 = second-ink(fill)
  table.cell(fill: fill, {
    set par(leading: 0.34em, justify: false)
    if a != none { text(weight: "bold", size: 8.4pt, fill: ink, dur(a)); spread(rr, fill) } else {
      text(size: 7.6pt, fill: if rr.status == "failed" { rgb("b42318") } else { muted }, rr.status)
    }
    linebreak()
    if ref.status == "n/a" { text(size: 6.3pt, fill: luma(120), alone-note) }
    else if b != none { text(size: 6.6pt, fill: ink2)[#label #dur(b)]; if ref-spread { spread(ref, fill) } }
    else { text(size: 6.6pt, fill: if ref.status == "failed" { rgb("b42318") } else { luma(110) })[#label #ref.status] }
    linebreak()
    if a != none and b != none { text(size: 7pt, weight: "semibold", fill: ink, ratio-text(a / b)) }
    else if cores-of(rr) > 1 { cores-line(rr, fill) } else { blank }
    if a != none and b != none and cores-of(rr) > 1 { linebreak(); cores-line(rr, fill) }
    if extra != none { linebreak(); text(size: 6.1pt, fill: ink2, extra) }
  })
}
// Warm RustRed call against the cold one.
#let warm(case) = cmp(cell-of(case, "rustred.reduction.warm"), cell-of(case, "rustred.reduction.cold"), "cold", ref-spread: false)
#let agree(case) = {
  let a = case.agreement
  table.cell({
    set par(leading: 0.32em, justify: false)
    if a.status == "n/a" or a.min_digits == none { text(size: 7.4pt, fill: luma(130))[—] }
    else {
      let ok = a.status == "agree"
      text(size: 8.2pt, weight: "bold", fill: if ok { rgb("1b7a3e") } else { rgb("b42318") })[#fixed(a.min_digits, 0)]
      linebreak()
      text(size: 6.4pt, fill: luma(90))[#if ok [digits] else [DISAGREE]]
    }
  })
}
#let ref-pair(case) = {
  let l = case.loops
  if l <= 3 { ("matad", "M") } else if l == 4 { ("fmft", "F") } else { ("matad", "M") }
}
#let cell-label = (
  "rustred.reduction.cold": "reduction", "rustred.reduction.warm": "warm reduction",
  "rustred.full.default": "full call", "matad.reduction": "MATAD reduction", "matad.full": "MATAD full call",
  "alphaloop.reduction": "alphaLoop reduction", "alphaloop.full": "alphaLoop full call",
  "fmft.reduction": "FMFT reduction", "fmft.full": "FMFT full call",
)
#let case-label(c) = [#c.loops\L #c.title#if c.at("input", default: (:)).at("master", default: none) != none [ (#c.input.master)], #c.subtitle]

// ------------------------------------------------------------------ title --

#let core-claim = if timing.at("all_single_core", default: true) [every timed process on one core] else [
  one core per timed process except #pl(timing.multi_core_cells.len(), "marked cell")]
#align(center)[
  #v(0.2em)
  #text(22pt, weight: "bold", fill: navy)[Vakint reduction benchmark]
  #v(-0.35em)
  #text(12pt, fill: luma(80))[RustRed versus alphaLoop, MATAD and FMFT · reduction graphs from one to five loops · precision of the results]
  #v(0.15em)
  #text(9pt, fill: luma(110))[#day · #host.at("cpu_model", default: "") · #core-claim · vakint #ident.vakint_version]
]

#block(fill: pale, inset: 10pt, radius: 4pt, width: 100%)[
  #set text(9.3pt)
  *Reading guide.* Bold: RustRed median, then the reference's (M = MATAD, F = FMFT, αL = alphaLoop) and the ratio
  RustRed ÷ reference; green = RustRed faster, red = slower (log scale); (·) = minimum over the runs.
  _Reduction_: masters retained; _full_: complete call with master values; _warm_: the identical call repeated
  in the process (step cache filled). Every timed call follows an untimed warm-up on another integral, so one-time
  setup is excluded (listed below the table). _IBP_: RustRed's time in the reduction itself. Five loops (no FORM
  backend): a master is recognised in any form and returned without a reduction walk; other inputs run RustRed's
  modular ordered walk (`modular_sampled`, on the cores a cell names).
]

= Timings <timings>

#let header-style(body) = {
  set text(7.6pt, font: sans)
  set par(justify: false)
  set table.cell(breakable: false)
  show table.cell.where(y: 0): set text(fill: white, weight: "bold")
  body
}
#let groups = results.groups
#figure(
  kind: table,
  caption: [*RustRed versus the FORM backends* (MATAD and alphaLoop: one to three loops; FMFT: four loops; none
    reaches five). Terminals: distinct RustRed master symbols in the retained result (one to three loops: terminals of
    the closed program, symmetric copies included). Agreement: digits shared with RustRed at $m = mu = 1$, worst
    order of $epsilon$.],
  header-style(table(
    columns: (1fr, 1.1cm, 2.32cm, 2.32cm, 2.32cm, 2.32cm, 2.32cm, 1.25cm),
    align: (x, y) => if x == 0 { left + horizon } else { center + horizon },
    inset: (x: 4pt, y: 4pt),
    stroke: (x, y) => (top: if y <= 1 { none } else { 0.4pt + luma(225) },
                       left: if x >= 2 and y > 0 { 0.4pt + luma(225) } else { none }),
    fill: (x, y) => if y == 0 { navy } else { none },
    table.header([Test case], [Termi-\ nals], [Reduction\ vs M / F], [Reduction\ vs αL], [Warm\ reduction],
      [Full call\ vs M / F], [Full call\ vs αL], [Agree-\ ment]),
    ..for g in groups {
      let rows = results.cases.filter(c => c.loops == g.loops)
      if rows.len() == 0 { () } else {
        (table.cell(colspan: 8, fill: pale, text(size: 7.6pt, weight: "bold", fill: navy)[#g.label]),)
        for c in rows {
          let (rk, rl) = ref-pair(c)
          let cold = cell-of(c, "rustred.reduction.cold")
          let full = cell-of(c, "rustred.full.default")
          let nm = cold.at("masters_retained", default: none)
          (
            table.cell(text(size: 8.1pt)[#c.title#{let m = c.at("input", default: (:)).at("master", default: none); if m != none [ #text(fill: luma(95))[(#m)]]} \ #text(fill: luma(95), size: 7.4pt, c.subtitle)]),
            table.cell(text(size: 7.8pt)[#if nm != none { str(nm) } else { [—] }]),
            cmp(cold, cell-of(c, rk + ".reduction"), rl, extra: ibp-line(cold)),
            cmp(cold, cell-of(c, "alphaloop.reduction"), "αL", show-alone: false),
            warm(c),
            cmp(full, cell-of(c, rk + ".full"), rl),
            cmp(full, cell-of(c, "alphaloop.full"), "αL", show-alone: false),
            agree(c),
          )
        }
      }
    },
  )),
) <tab:timings>

#align(center)[
  #set text(7.5pt, font: sans, fill: luma(60))
  #box(width: 9cm, height: 0.32cm, fill: palette, radius: 2pt)
  #linebreak()
  #box(width: 9cm)[×1/64 #h(1fr) ×1/8 #h(1fr) ×1 #h(1fr) ×8 #h(1fr) ×64]
]

// One-time cost per process (excluded from every cell).
#let once = results.at("one_time_per_process", default: (:))
#if once.len() > 0 {
  let loops = range(1, 6).map(str).filter(l => l in once)
  let get(l, cell, key) = once.at(l).at(cell, default: (:)).at(key, default: none)
  let loader(l) = { let v = get(l, "rustred.reduction.cold", "loader_median"); if v == none { get(l, "rustred.full.default", "loader_median") } else { v } }
  let rows = (
    ([RustRed data load], l => loader(l)),
    ([RustRed full call], l => get(l, "rustred.full.default", "setup_median")),
    ([MATAD full call], l => get(l, "matad.full", "setup_median")),
    ([alphaLoop full call], l => get(l, "alphaloop.full", "setup_median")),
    ([FMFT full call], l => get(l, "fmft.full", "setup_median")),
  )
  block(width: 100%, above: 0.9em)[
    #set text(7.6pt, font: sans)
    #text(weight: "bold", fill: navy)[Once per process] #text(fill: luma(90))[(excluded from every cell; first-call setup = first minus second call on the warm-up integral)]
    #v(-0.3em)
    #table(columns: (3.4cm,) + loops.map(_ => 1.6cm), stroke: none, inset: (x: 3pt, y: 2.2pt),
      align: (x, y) => if x == 0 { left } else { right },
      table.hline(stroke: 0.4pt + luma(200)),
      [], ..loops.map(l => text(fill: navy, weight: "bold")[#l\L]),
      table.hline(stroke: 0.4pt + luma(220)),
      ..for (label, f) in rows {
        let vals = loops.map(l => f(l))
        if vals.all(v => v == none) { () } else {
          (label, ..vals.map(v => if v == none { text(fill: luma(160))[—] } else { dur(calc.max(v, 0)) }))
        }
      },
      table.hline(stroke: 0.4pt + luma(200)),
    )
  ]
}

#let failed = results.cases.map(c => c.cells.pairs().filter(p => p.at(1).status == "failed").map(p => (c, p.at(0), p.at(1)))).flatten().chunks(3)
#let pending = results.cases.map(c => (c, c.cells.pairs().filter(p => p.at(1).status == "pending").map(p => cell-label.at(p.at(0), default: p.at(0))))).filter(p => p.at(1).len() > 0)
#if failed.len() > 0 {
  block(fill: rgb("fdecea"), inset: 6pt, radius: 3pt, width: 100%, text(size: 7.3pt)[
    *Failed calls.*
    #for (c, k, v) in failed {
      let s = v.at("error_summary", default: (text: v.at("error", default: "")))
      [#linebreak()#case-label(c) · #cell-label.at(k, default: k) (#v.at("n", default: 1)/#v.at("n", default: 1) runs): #text(fill: luma(60), s.text)]
    }])
}
#if pending.len() > 0 {
  block(fill: rgb("fff3df"), inset: 7pt, radius: 3pt, width: 100%, text(size: 7.6pt)[
    *Not measured in this build* (#pl(pending.map(p => p.at(1).len()).sum(), "cell")):
    #pending.map(((c, ks)) => [#case-label(c): #ks.join(", ")]).join("; ").])
}

#block(width: 100%, inset: (x: 0pt, y: 4pt))[
  #text(9pt, weight: "bold", fill: navy)[Provenance]
#{
  set text(7.6pt)
  set par(justify: false)
  let src = ident.at("source", default: (:))
  let dg = ident.at("data_git", default: none)
  let cn = results.at("contention", default: (:))
  let pct(x) = if x == none { [—] } else [#calc.round(x * 100, digits: 1)%]
  let form = host.at("form_version", default: "—").replace(regex("\(\w{3}\s+\d+\s+\d{4},?\s*\)?"), "").replace(regex("\s+\(\s*\)"), "").trim()
  table(columns: (auto, 1fr), stroke: none, inset: (x: 3pt, y: 2pt),
    [Data], if dg != none and dg.at("commit", default: none) != none [
      vakint-external-data #raw(str(dg.at("describe", default: dg.commit))) · commit #raw(dg.commit.slice(0, 12))#if dg.at("dirty_files", default: 0) > 0 [ (+#pl(dg.dirty_files, "modified file"))]
    ] else [snapshot (no git checkout): #raw(str(ident.data_root))],
    [Data manifest], raw(str(ident.at("root_manifest_sha256", default: "—"))),
    [Vakint], [#ident.vakint_version · built from #raw(str(src.at("commit", default: "unknown")).slice(0, calc.min(12, str(src.at("commit", default: "unknown")).len())))#{
      let d = src.at("dirty_paths", default: none)
      if d == none [ (build identity unknown)] else if d > 0 [ (+#pl(d, "modified path"))] else [ (clean)]
    } · profile #src.at("profile", default: "—")],
    [Harness binary], raw(str(ident.executable_sha256)),
    [FORM], [#form],
    [Host], [#host.at("cpu_model", default: "—"), #host.at("logical_cpus", default: "—") logical CPUs, shared; timed on #pl(host.at("cpus_used", default: ()).len(), "pinned CPU") in #pl(cn.at("l3_domains", default: 0), "L3 domain") (listed in suite.json)],
    [Contention], [#cn.at("samples", default: "—") samples, #cn.at("contended_samples", default: "—") still contended after re-measurement (multi-core samples are measured once; thresholds: SMT sibling #pct(cn.at("thresholds", default: (:)).at("sibling_busy", default: none)), other tasks #pct(cn.at("thresholds", default: (:)).at("foreign_cpu", default: none)), call run-queue wait #pct(cn.at("thresholds", default: (:)).at("run_delay", default: none)) when other tasks ran on the timed CPU) · SMT sibling busy #pct(cn.at("sibling_busy_median", default: none)) median (max #pct(cn.at("sibling_busy_max", default: none))) · timed call waiting in the run queue #pct(cn.at("call_run_delay_median", default: none)) (max #pct(cn.at("call_run_delay_max", default: none))) · other tasks on the timed CPU #pct(cn.at("foreign_cpu_median", default: none)) (max #pct(cn.at("foreign_cpu_max", default: none))) · 1-min load ≤ #cn.at("load1_max", default: "—")],
    [Regenerate], [#raw("crates/vakint/docs/benchmark/make-benchmark.sh DATA_ROOT BINARY OUT --install"); after a data release #raw("DATA_ROOT - OUT --refresh-tables --install") (timings kept)],
    [Built], [#results.assembled_utc #if draft [· *DRAFT*: not the published data or build]],
  )
}
]

// ---------------------------------------------------------------- graphs --

#pagebreak(weak: true)
= Reduction graphs <graphs>

Nodes are the vacuum topologies of the shipped reduction, on the row of their number of lines; an edge removes one
line (thicker: several lines), as in Fig. 21 of B. Ruijl's thesis (Leiden, 2017). Products (one-vertex-reducible
graphs) are merged by their factors; numbers count masters (one to three loops: RustRed terminals).

#let c-master = rgb("2a9d8f")
#let c-master-line = rgb("1d6f65")
#let c-product = rgb("f4a259")
#let c-product-line = rgb("c9762a")
#let c-away-line = luma(140)
#let node-class(n) = if n.product { "product" } else if n.masters.len() > 0 { "master" } else { "away" }
#let node-fill(n) = {
  let k = node-class(n)
  if k == "master" { c-master } else if k == "product" { if n.masters.len() > 0 { c-product } else { white } } else { white }
}
#let node-stroke(n) = {
  let k = node-class(n)
  0.7pt + if k == "master" { c-master-line } else if k == "product" { c-product-line } else { c-away-line }
}
#let edge-paint(n, alpha) = {
  let k = node-class(n)
  let (r, g, b) = if k == "master" { (42, 157, 143) } else if k == "product" { (214, 130, 60) } else { (110, 120, 135) }
  rgb(r, g, b, alpha)
}
#let legend(columns: 5) = {
  set text(7.4pt, font: sans, fill: luma(50))
  set par(justify: false)
  let sw(f, s) = box(baseline: 1.5pt, ellipse(width: 11pt, height: 7.5pt, fill: f, stroke: s))
  align(center, grid(columns: columns, column-gutter: 12pt, row-gutter: 5pt, align: left + horizon,
    [#sw(c-master, 0.7pt + c-master-line) topology with masters],
    [#sw(c-product, 0.7pt + c-product-line) product of lower-loop masters],
    [#sw(white, 0.7pt + c-product-line) product, no master of its own],
    [#sw(white, 0.7pt + c-away-line) reduced away],
    [#box(baseline: 1.5pt, ellipse(width: 11pt, height: 7.5pt, fill: white, stroke: 1.6pt + navy)) owns a rule program],
  ))
}
#legend(columns: 3)

#let draw-graph(g, width: 17cm, row: 1.25cm, labels: true, rx: 0.62cm, ry: 0.27cm, label-size: 6.8pt, edge-alpha: 150, edge-base: 0.35pt, max-unit: 3.2cm) = {
  let nodes = g.nodes
  let layers = g.layers
  let span = calc.max(g.span, 1)
  let maxb = calc.max(1, ..nodes.map(n => n.owner_bytes))
  let left = 1.3cm
  let w = width - left - 2 * rx
  let pos = (:)
  let by-id = (:)
  let wide = nodes.any(n => n.x > 0)
  let unit = if wide { calc.min(w / span, max-unit) } else { 0cm }
  let offset = if wide { (w - unit * span) / 2 } else { w / 2 }
  for n in nodes {
    pos.insert(n.id, (left + rx + offset + n.x * unit, ry + 0.15cm + n.layer * row))
    by-id.insert(n.id, n)
  }
  let height = (layers - 1) * row + 2 * ry + 0.5cm
  let lines-of = (:)
  for n in nodes { lines-of.insert(str(n.layer), n.lines) }
  box(width: width, height: height, {
    // layer labels
    for (layer, lines) in lines-of {
      place(dx: 0pt, dy: ry + 0.15cm + int(layer) * row - 0.18cm,
        text(size: 6.8pt, font: sans, fill: luma(120))[#pl(lines, "line")])
    }
    // edges first, coloured by the class of the topology they reach; long edges fainter
    let edges = g.arrows.map(e => {
      let (x1, _) = pos.at(e.from)
      let (x2, _) = pos.at(e.to)
      (e, calc.abs((x2 - x1) / w))
    }).sorted(key: p => -p.at(1))
    for (e, d) in edges {
      let (x1, y1) = pos.at(e.from)
      let (x2, y2) = pos.at(e.to)
      let y1b = y1 + ry
      let y2t = y2 - ry
      let dy = (y2t - y1b) * 0.5
      let wdt = edge-base + 0.45pt * calc.log(e.lines, base: 2)
      let a = int(edge-alpha * calc.clamp(1.15 - 2.2 * d, 0.3, 1))
      place(curve(
        stroke: (paint: edge-paint(by-id.at(e.to), a), thickness: wdt, cap: "round"),
        curve.move((x1, y1b)),
        curve.cubic((x1, y1b + dy), (x2, y2t - dy), (x2, y2t)),
      ))
    }
    // nodes
    for n in nodes {
      let (x, y) = pos.at(n.id)
      let ring = if n.owner_bytes > 0 { 1.1pt + 1.1pt * calc.log(1 + n.owner_bytes) / calc.log(1 + maxb) } else { 0pt }
      if ring > 0pt {
        place(dx: x - rx - ring, dy: y - ry - ring,
          ellipse(width: 2 * rx + 2 * ring, height: 2 * ry + 2 * ring, fill: navy, stroke: none))
      }
      let fill = node-fill(n)
      place(dx: x - rx, dy: y - ry, ellipse(width: 2 * rx, height: 2 * ry, fill: fill, stroke: node-stroke(n)))
      let m = n.masters.len()
      let on-dark = fill == c-master
      let label = if labels { n.name } else if m > 0 { str(m) } else { "" }
      place(dx: x - rx, dy: y - ry, box(width: 2 * rx, height: 2 * ry, align(center + horizon,
        text(size: label-size, font: sans, fill: if on-dark { white } else { ink }, weight: if labels { "regular" } else { "bold" })[#label])))
      if labels and m > 0 {
        place(dx: x + rx * 0.62, dy: y - ry - 0.22cm, text(size: 6pt, font: sans, weight: "bold",
          fill: if n.product { c-product-line } else { c-master-line })[#m])
      }
    }
  })
}

#let gstats(L) = {
  let s = graphs.loops.at(str(L)).stats
  let what = if L <= 3 { "terminal" } else { "master" }
  [#pl(s.topologies, "topology", many: "topologies") · #pl(s.transitions, "edge") · #pl(s.masters, what)]
}
#let fig-graph(L, body, cap) = figure(body, caption: [*#(L)-loop reduction graph.* #cap #gstats(L).], kind: image)

// equal heights (the two-loop graph's), so that both captions start on the same line
#grid(columns: (1fr, 1fr), column-gutter: 8mm, align: top,
  fig-graph(1, block(height: 1.1cm + 2 * 0.27cm + 0.5cm, align(horizon, draw-graph(graphs.loops.at("1"), width: 8.2cm, row: 1cm))),
    [T1 carries the one-loop master; dots and numerators reduce inside it.]),
  fig-graph(2, draw-graph(graphs.loops.at("2"), width: 8.2cm, row: 1.1cm),
    [Removing any sunset line gives T1².]),
)
#fig-graph(3, draw-graph(graphs.loops.at("3"), width: 17cm, row: 1.35cm, rx: 0.7cm, ry: 0.3cm, label-size: 7.2pt),
  [MATAD names; symmetric sectors are separate terminals of the closed program.])
#fig-graph(4, draw-graph(graphs.loops.at("4"), width: 17.8cm, row: 2.1cm, rx: 0.82cm, ry: 0.34cm, label-size: 7.4pt),
  [Rings: the FMFT family programs H, X, BMW, FG; irreducible topologies named by their FMFT master PR$k$.])

#page(flipped: true, margin: (x: 12mm, y: 12mm))[
  #legend()
  #fig-graph(5, draw-graph(graphs.loops.at("5"), width: 27cm, row: 2.3cm, labels: false, rx: 0.36cm, ry: 0.27cm, label-size: 6.6pt, edge-alpha: 120, edge-base: 0.25pt),
    [Rings: the #graphs.loops.at("5").stats.owner_programs owner sectors of the saved application (width grows with their rule tables).])
]

// ------------------------------------------------------------- precision --

#precision-section(tables)
