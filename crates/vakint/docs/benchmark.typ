// Vakint reduction benchmark. Generated data only: regenerate with
//   crates/vakint/docs/benchmark/make-benchmark.sh DATA_ROOT BENCHMARK_BINARY OUT_DIR [--install]
// Inputs (typst --input): results=<results.json>, graphs=<graphs.json>, draft=true|false.

#let results = json(sys.inputs.at("results", default: "benchmark/data/suite.json"))
#let graphs = json(sys.inputs.at("graphs", default: "benchmark/data/graphs.json"))
#let draft = sys.inputs.at("draft", default: if results.at("draft", default: false) { "true" } else { "false" }) == "true"

#let navy = rgb("173f65")
#let ink = rgb("243343")
#let pale = rgb("eaf2fa")
#let muted = rgb("5b6876")
#let sans = ("DejaVu Sans", "Roboto", "Noto Sans", "Libertinus Serif")
#let ident = results.identity
#let host = ident.at("host", default: (:))
#let day = results.assembled_utc.slice(0, 10)

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
#show raw: set text(font: "DejaVu Sans Mono", size: 8pt)
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)
#show figure.where(kind: table): set block(breakable: true)

// ---------------------------------------------------------------- helpers --

#let is-num(x) = type(x) == int or type(x) == float
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
// Seconds with an adapted unit.
#let dur(s) = {
  if not is-num(s) { return [—] }
  if s < 1e-3 { [#sig(s * 1e6)#h(1pt)µs] }
  else if s < 1 { [#sig(s * 1e3)#h(1pt)ms] }
  else if s < 600 { [#sig(s)#h(1pt)s] }
  else if s < 36000 { [#sig(s / 60)#h(1pt)min] }
  else { [#sig(s / 3600)#h(1pt)h] }
}
#let ratio-text(r) = {
  if r < 0.095 { [×1/#str(int(calc.round(1 / r)))] }
  else if r < 10 { [×#fixed(r, if r < 1 { 2 } else { 1 })] }
  else { [×#str(int(calc.round(r)))] }
}
// Diverging palette, log2 scale: ×1/64 (deep green, RustRed faster) to ×64 (deep red).
#let span2 = 6
#let palette = gradient.linear(rgb("#3f9a4c"), rgb("#a6d79b"), rgb("#ffffff"), rgb("#f9c4a8"), rgb("#cf4f3d"))
#let ratio-fill(r) = palette.sample(calc.clamp(0.5 + calc.log(r, base: 2) / (2 * span2), 0, 1) * 100%)
#let na-fill = rgb("f3f5f8")

#let cell-of(case, key) = case.cells.at(key, default: (status: "pending"))
#let med(c) = if c.status == "ok" { c.median } else { none }

// One comparison cell: RustRed value in bold, the reference value and RustRed/reference below.
// `show-alone`: when the reference does not apply, still show the RustRed value (else only "n/a").
#let ibp-line(c) = {
  let m = c.at("reduction_math_median", default: none)
  if c.status == "ok" and m != none { text(size: 6.1pt, fill: luma(105))[IBP #dur(m)] }
}
#let cmp(rr, ref, label, show-alone: true, alone-note: [—], extra: none) = {
  let a = med(rr)
  let b = med(ref)
  if ref.status == "n/a" and not show-alone {
    return table.cell(fill: na-fill, text(size: 7pt, fill: luma(150))[n/a])
  }
  let fill = if a != none and b != none { ratio-fill(a / b) } else if ref.status == "n/a" { na-fill } else { none }
  table.cell(fill: fill, {
    set par(leading: 0.32em, justify: false)
    if a != none { text(weight: "bold", size: 8.4pt, dur(a)) } else {
      text(size: 7.6pt, fill: if rr.status == "failed" { rgb("b42318") } else { muted }, rr.status)
    }
    linebreak()
    if ref.status == "n/a" { text(size: 6.4pt, fill: luma(130), alone-note) }
    else if b != none {
      text(size: 6.6pt, fill: luma(80))[#label #dur(b)]
      if a != none { text(size: 6.9pt, weight: "semibold")[#h(3pt)#ratio-text(a / b)] }
    } else { text(size: 6.6pt, fill: if ref.status == "failed" { rgb("b42318") } else { luma(110) })[#label #ref.status] }
    if extra != none { linebreak(); extra }
  })
}
// Warm RustRed call against the cold one.
#let warm(case) = {
  let w = cell-of(case, "rustred.reduction.warm")
  let c = cell-of(case, "rustred.reduction.cold")
  cmp(w, c, "cold", extra: ibp-line(w))
}
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

// ------------------------------------------------------------------ title --

#align(center)[
  #v(0.2em)
  #text(22pt, weight: "bold", fill: navy)[Vakint reduction benchmark]
  #v(-0.35em)
  #text(12pt, fill: luma(80))[RustRed versus alphaLoop, MATAD and FMFT · reduction graphs from one to five loops]
  #v(0.15em)
  #text(9pt, fill: luma(110))[#day · #host.at("cpu_model", default: "") · every timed process on one core · vakint #ident.vakint_version]
]

#block(fill: pale, inset: 10pt, radius: 4pt, width: 100%)[
  #set text(9.3pt)
  *Reading guide.* Each coloured cell shows the RustRed time in bold, the reference below it and the
  ratio _RustRed ÷ reference_; green means RustRed is faster, red slower, white is parity
  (log scale, ×1/64 to ×64). _Reduction_: the public `evaluate` call with masters retained
  (RustRed `Retain`; FORM backends without master expansion or substitution). _Full_: the complete
  public call with master values (RustRed `Numerical`; FORM backends with their defaults).
  _Cold_: fresh process, no step cache; _warm_: the same call repeated in that process (Ordered step
  cache populated); _full_ uses the engine defaults (shipped step caches on). Engine creation and data
  loading are excluded. Median of fresh processes on one CPU, normalization $overline("MS")$, $L+1$
  orders in $epsilon$. _IBP_: the part of a RustRed call spent in the reduction itself (rule
  application and ordered scheduling); the rest is input processing, topology matching, numerator
  lowering and output. _Agreement_: decimal digits shared by every backend's numerical result with
  RustRed's at $m = mu = 1$, worst order.
]

= Timings <timings>

#let header-style(body) = {
  set text(7.6pt, font: sans)
  set par(justify: false)
  show table.cell.where(y: 0): set text(fill: white, weight: "bold")
  body
}
#let groups = results.groups
#figure(
  kind: table,
  caption: [*RustRed versus the FORM backends.* Reference letters: M = MATAD, F = FMFT, αL = alphaLoop.
    MATAD and alphaLoop apply to one to three loops, FMFT to four loops; no FORM backend reaches five loops
    (n/a). Masters: distinct masters in RustRed's retained result.],
  header-style(table(
    columns: (1fr, 1.15cm, 2.3cm, 2.3cm, 2.3cm, 2.3cm, 2.3cm, 1.3cm),
    align: (x, y) => if x == 0 { left + horizon } else { center + horizon },
    inset: (x: 4pt, y: 4pt),
    stroke: (x, y) => (top: if y <= 1 { none } else { 0.4pt + luma(225) },
                       left: if x >= 2 and y > 0 { 0.4pt + luma(225) } else { none }),
    fill: (x, y) => if y == 0 { navy } else { none },
    table.header([Test case], [Masters], [Reduction\ vs M / F], [Reduction\ vs αL], [Warm\ reduction],
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
            table.cell(text(size: 8.1pt)[#c.title \ #text(fill: luma(95), size: 7.4pt, c.subtitle)]),
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

#let pend = results.cases.map(c => c.cells.pairs().filter(p => p.at(1).status == "pending").map(p => c.id + " " + p.at(0))).flatten()
#let failed = results.cases.map(c => c.cells.pairs().filter(p => p.at(1).status == "failed").map(p => (c.id + " " + p.at(0), p.at(1).at("error", default: "")))).flatten().chunks(2)
#if failed.len() > 0 {
  block(fill: rgb("fdecea"), inset: 7pt, radius: 3pt, width: 100%, text(size: 7.6pt)[
    *Failed calls* (the error is the backend's own):
    #for (cell, err) in failed [#linebreak()#raw(cell): #text(fill: luma(70), err.slice(0, calc.min(err.len(), 240)))]])
}
#if pend.len() > 0 {
  block(fill: rgb("fff3df"), inset: 7pt, radius: 3pt, width: 100%, text(size: 7.6pt)[
    *Not measured in this build* (#pend.len() cells): #pend.join(", ").])
}

#block(width: 100%, inset: (x: 0pt, y: 4pt))[
  #text(9pt, weight: "bold", fill: navy)[Provenance]
#{
  set text(8pt)
  let src = ident.at("source", default: (:))
  table(columns: (auto, 1fr), stroke: none, inset: (x: 3pt, y: 2pt),
    [Data root], raw(str(ident.data_root)),
    [Root manifest], raw(str(ident.at("root_manifest_sha256", default: "—"))),
    [Vakint], [#ident.vakint_version · source #raw(str(src.at("commit", default: "—")).slice(0, 12)) #if src.at("dirty_files", default: 0) > 0 [(+#src.dirty_files modified files)]],
    [Harness binary], raw(str(ident.executable_sha256)),
    [FORM], [#host.at("form_version", default: "—")],
    [Host], [#host.at("cpu_model", default: "—"), #host.at("logical_cpus", default: "—") logical CPUs; timed on CPUs #host.at("cpus_used", default: ()).map(str).join(", ")],
    [Contention], [#{let c = results.at("contention", default: (:)); [SMT sibling of the timed CPU busy #c.at("sibling_busy_median", default: "—") (median), #c.at("sibling_busy_max", default: "—") (max); 1-min load up to #c.at("load1_max", default: "—")]}],
    [Regenerate], raw("crates/vakint/docs/benchmark/make-benchmark.sh DATA_ROOT BENCHMARK_BINARY OUT_DIR --install"),
    [Built], [#results.assembled_utc #if draft [· *DRAFT*: not the published data]],
  )
}
]

// ---------------------------------------------------------------- graphs --

= Reduction graphs <graphs>

Each node is a vacuum topology of the shipped reduction (graphs up to isomorphism), on the row of its
number of lines; an edge joins a topology to the one obtained by removing a line (thicker: several lines
give the same topology), following Fig. 21 of B. Ruijl's thesis (Leiden, 2017). Removing a line that is a loop of its
own gives zero and is not drawn. One-vertex-reducible graphs (vacuum products) are merged by their
factors. Colour gives the role in RustRed's reduction:

#let teal(t) = color.mix((rgb("ffffff"), (1 - t) * 100%), (rgb("0f766e"), t * 100%))
#let green(t) = color.mix((rgb("ffffff"), (1 - t) * 100%), (rgb("3f9a4c"), t * 100%))
#let node-fill(n, maxm) = {
  let m = n.masters.len()
  if m == 0 { white } else {
    let t = 0.28 + 0.62 * calc.log(1 + m) / calc.log(1 + maxm)
    if n.product { green(t) } else { teal(t) }
  }
}
#let node-stroke(n) = if n.product { 0.7pt + rgb("3f9a4c") } else if n.masters.len() > 0 { 0.7pt + rgb("0f766e") } else { 0.7pt + luma(120) }
#let legend = {
  set text(7.6pt, font: sans, fill: luma(50))
  let sw(f, s) = box(baseline: 1.5pt, ellipse(width: 10pt, height: 7pt, fill: f, stroke: s))
  grid(columns: 3, column-gutter: 14pt, row-gutter: 5pt, align: left + horizon,
    [#sw(teal(0.75), 0.7pt + rgb("0f766e")) topology with masters (darker: more)],
    [#sw(green(0.75), 0.7pt + rgb("3f9a4c")) product of lower-loop masters],
    [#sw(white, 0.7pt + luma(120)) reduced away by the rule tables],
    [#sw(white, 0.7pt + rgb("3f9a4c")) product without own master],
    [#box(baseline: 1.5pt, ellipse(width: 10pt, height: 7pt, fill: white, stroke: 1.6pt + navy)) owns a rule program],
    [#box(baseline: 1.5pt, line(length: 14pt, stroke: 1.2pt + luma(140))) line removal (width: lines)],
  )
}
#legend

#let draw-graph(g, width: 17cm, row: 1.25cm, labels: true, rx: 0.62cm, ry: 0.27cm, label-size: 6.8pt, edge-alpha: 150, edge-base: 0.35pt, max-unit: 3.2cm) = {
  let nodes = g.nodes
  let layers = g.layers
  let span = calc.max(g.span, 1)
  let maxm = calc.max(1, ..nodes.map(n => n.masters.len()))
  let maxb = calc.max(1, ..nodes.map(n => n.owner_bytes))
  let left = 1.3cm
  let w = width - left - 2 * rx
  let pos = (:)
  let wide = nodes.any(n => n.x > 0)
  let unit = if wide { calc.min(w / span, max-unit) } else { 0cm }
  let offset = if wide { (w - unit * span) / 2 } else { w / 2 }
  for n in nodes {
    pos.insert(n.id, (left + rx + offset + n.x * unit, ry + 0.15cm + n.layer * row))
  }
  let height = (layers - 1) * row + 2 * ry + 0.5cm
  let lines-of = (:)
  for n in nodes { lines-of.insert(str(n.layer), n.lines) }
  box(width: width, height: height, {
    // layer labels
    for (layer, lines) in lines-of {
      place(dx: 0pt, dy: ry + 0.15cm + int(layer) * row - 0.18cm,
        text(size: 6.8pt, font: sans, fill: luma(120))[#lines #if lines == 1 [line] else [lines]])
    }
    // edges first
    for e in g.arrows {
      let (x1, y1) = pos.at(e.from)
      let (x2, y2) = pos.at(e.to)
      let y1b = y1 + ry
      let y2t = y2 - ry
      let dy = (y2t - y1b) * 0.5
      let wdt = edge-base + 0.45pt * calc.log(e.lines, base: 2)
      place(curve(
        stroke: (paint: rgb(110, 122, 140, edge-alpha), thickness: wdt, cap: "round"),
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
      place(dx: x - rx, dy: y - ry, ellipse(width: 2 * rx, height: 2 * ry, fill: node-fill(n, maxm), stroke: node-stroke(n)))
      let m = n.masters.len()
      let dark = m > 0 and (0.28 + 0.62 * calc.log(1 + m) / calc.log(1 + maxm)) > 0.6
      let label = if labels { n.name } else if m > 0 { str(m) } else { "" }
      place(dx: x - rx, dy: y - ry, box(width: 2 * rx, height: 2 * ry, align(center + horizon,
        text(size: label-size, font: sans, fill: if dark { white } else { ink }, weight: if labels { "regular" } else { "bold" })[#label])))
      if labels and m > 0 {
        place(dx: x + rx * 0.62, dy: y - ry - 0.21cm, text(size: 5.8pt, font: sans, fill: rgb("0f766e"), weight: "bold")[#m])
      }
    }
  })
}

#let gstats(L) = {
  let g = graphs.loops.at(str(L))
  let s = g.stats
  [#s.topologies topologies (#s.irreducible irreducible, #s.products products), #s.removable_lines removable lines
   in #s.transitions edges, #s.masters masters on #s.master_topologies topologies]
}
#let fig-graph(L, body, cap) = figure(body, caption: [*#(L)-loop reduction graph.* #cap #gstats(L).], kind: image)

#grid(columns: (1fr, 1fr), column-gutter: 8mm,
  fig-graph(1, draw-graph(graphs.loops.at("1"), width: 8.2cm, row: 1cm),
    [The tadpole T1 carries the one-loop master; dots and numerators reduce inside it.]),
  fig-graph(2, draw-graph(graphs.loops.at("2"), width: 8.2cm, row: 1.1cm),
    [Removing any of the three sunset lines gives the product T1².]),
)
#fig-graph(3, draw-graph(graphs.loops.at("3"), width: 17cm, row: 1.1cm),
  [MATAD names; small numbers count the RustRed terminals of a topology (symmetric sectors stay separate terminals of the closed three-loop program).])
#fig-graph(4, draw-graph(graphs.loops.at("4"), width: 17cm, row: 1.15cm, rx: 0.7cm),
  [Rings: the four FMFT family programs H, X, BMW, FG; irreducible topologies are named by their FMFT master PR$k$.])

#page(flipped: true, margin: (x: 12mm, y: 12mm))[
  #fig-graph(5, draw-graph(graphs.loops.at("5"), width: 27cm, row: 1.95cm, labels: false, rx: 0.36cm, ry: 0.27cm, label-size: 6.6pt, edge-alpha: 105, edge-base: 0.25pt),
    [Numbers: masters of the shipped basis in a topology. Rings: topologies holding the #graphs.loops.at("5").stats.owner_programs owner sectors of the saved application (ring width grows with the size of their rule tables).])
  #legend
]
