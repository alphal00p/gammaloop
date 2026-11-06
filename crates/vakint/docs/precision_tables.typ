// Precision tables of the benchmark document (generated data only):
//   crates/vakint/docs/benchmark/tools/precision_tables.py --data-root DATA_ROOT --json tables.json
// `precision-section(t)` renders table 1 (every master x every Laurent order the data hold) and table 2 (the L-loop
// result of any input, worst case over inputs, spurious poles included) from that JSON.

#let pt-sans = ("DejaVu Sans", "Roboto", "Noto Sans", "Libertinus Serif")
#let pt-navy = rgb("173f65")
#let pt-style = (
  "S": (fill: rgb("cdebc5"), ink: rgb("1b5e20")),
  "S*": (fill: rgb("fde0b4"), ink: rgb("8a4b00")),
  "N": (fill: rgb("f8c9c4"), ink: rgb("8e1b12")),
)
#let pt-na = rgb("f3f5f8")
#let pt-key(o) = if o < 0 { "-" + str(-o) } else { str(o) }
#let pt-eps(o) = [$epsilon^(#o)$]
#let pt-text(c) = if c.class == "S" [S] else [#c.class\(#str(c.X)\)]
#let pt-marks(c, size) = {
  let m = c.at("marks", default: ())
  if m.len() > 0 { text(size: size, weight: "regular", super(m.join(","))) }
}

// One table-1 cell: coloured class, '—' where the order is not held.
#let pt-cell1(c) = {
  if c == none { return table.cell(fill: pt-na, text(fill: luma(170))[—]) }
  let s = pt-style.at(c.class)
  table.cell(fill: s.fill, text(fill: s.ink, weight: "semibold")[#pt-text(c)#pt-marks(c, 4.6pt)])
}

// One table-2 cell: large coloured cell; partial = grey '—' with its note letter.
#let pt-cell2(c) = {
  if c.class == none {
    if c.at("na", default: none) == "below" { return table.cell(fill: white, []) }
    if c.at("na", default: none) == "partial" {
      return table.cell(fill: pt-na, {
        text(size: 9pt, fill: luma(120))[—#pt-marks(c, 6pt)]
        linebreak()
        text(size: 5.6pt, fill: luma(110))[some inputs: #c.returned]
      })
    }
    return table.cell(fill: pt-na, text(size: 9pt, fill: luma(170))[—])
  }
  let s = pt-style.at(c.class)
  table.cell(fill: s.fill, {
    text(size: 9.2pt, weight: "bold", fill: s.ink, if c.class == "S" [S] else [#c.class])
    if c.class != "S" { linebreak(); text(size: 7.4pt, weight: "semibold", fill: s.ink)[#str(c.X) digits] }
    pt-marks(c, 6pt)
  })
}

#let pt-legend = {
  set text(7.6pt, font: pt-sans)
  let sw(k, body) = box(baseline: 2pt, rect(width: auto, height: 9pt, inset: (x: 3pt, y: 1.5pt), radius: 1.5pt,
    fill: pt-style.at(k).fill, text(size: 6.8pt, weight: "bold", fill: pt-style.at(k).ink, body)))
  grid(columns: 4, column-gutter: 12pt, align: left + horizon,
    [#sw("S")[S] exact, periods at any precision],
    [#sw("S*")[S\*(X)] exact, periods to at most X digits],
    [#sw("N")[N(X)] numerical only, at most X digits],
    [#box(baseline: 2pt, rect(height: 9pt, inset: (x: 3pt, y: 1.5pt), radius: 1.5pt, fill: pt-na, text(size: 6.8pt, fill: luma(130))[—])) not available],
  )
}

#let pt-table2(t) = {
  let o2 = t.orders_table2
  set text(font: pt-sans)
  table(
    columns: (1.2cm,) + o2.map(_ => 1fr),
    align: center + horizon,
    inset: (x: 2pt, y: 5pt),
    stroke: 0.6pt + white,
    table.header(
      table.cell(fill: pt-navy, text(size: 7.6pt, fill: white, weight: "bold")[L]),
      ..o2.map(o => table.cell(fill: pt-navy, text(size: 7.6pt, fill: white, weight: "bold", pt-eps(o)))),
    ),
    ..for r in t.table2 {
      (table.cell(fill: rgb("eaf2fa"), text(size: 8.6pt, weight: "bold", fill: pt-navy)[#r.loop\L]),
       ..o2.map(o => pt-cell2(r.cells.at(pt-key(o)))))
    },
  )
}

#let pt-table1(t) = {
  let o1 = t.orders_table1
  let counts = t.counts_table1
  set text(font: pt-sans, size: 5.5pt)
  table(
    // wider columns where S*(X) labels occur
    columns: (1.9cm,) + o1.map(o => if t.table1.any(r => {
      let c = r.cells.at(pt-key(o))
      c != none and c.class == "S*"
    }) { 1.75fr } else { 1fr }),
    align: (x, y) => if x == 0 { left + horizon } else { center + horizon },
    inset: (x: 1.4pt, y: 1.22pt),
    stroke: 0.4pt + white,
    table.header(
      table.cell(fill: pt-navy, text(fill: white, weight: "bold")[master]),
      ..o1.map(o => table.cell(fill: pt-navy, text(fill: white, weight: "bold", pt-eps(o)))),
    ),
    ..{
      let out = ()
      let last = none
      for r in t.table1 {
        if r.loop != last {
          let c = counts.at(str(r.loop) + "L")
          let n = t.table1.filter(x => x.loop == r.loop).len()
          out.push(table.cell(colspan: o1.len() + 1, fill: rgb("eaf2fa"),
            text(weight: "bold", fill: pt-navy)[#r.loop\L · #n master#if n != 1 [s] · #c.at("S") S, #c.at("S*") S\*, #c.at("N") N]))
          last = r.loop
        }
        let alias = r.at("name_640", default: none)
        out.push(table.cell(text(fill: rgb("1d2733"))[#r.loop\L/#r.name#if alias != none [ #text(fill: luma(120))[(#alias)]]]))
        for o in o1 { out.push(pt-cell1(r.cells.at(pt-key(o)))) }
      }
      out
    },
  )
}

// The whole section, on landscape pages: provenance line, legend, table 1 (two columns), table 2 with its notes, method.
#let precision-section(t) = page(flipped: true, margin: (x: 11mm, top: 11mm, bottom: 13mm), {
  let src = t.sources
  let f5 = src.five_loop_values
  let dg = t.at("data_git", default: none)
  [= Precision of the masters and of the results <precision>]
  v(-0.4em)
  block(fill: if t.status == "final" { rgb("eaf2fa") } else { rgb("fff3df") }, inset: 6pt, radius: 3pt, width: 100%,
    text(size: 7.2pt, font: pt-sans)[
    *#upper(t.status)* · data tree #if t.at("label", default: none) != none [#t.label] else [#raw(t.data_root)]
    #if dg != none [(git #raw(dg.commit.slice(0, 12)), tree #raw(dg.tree.slice(0, 12)))] ·
    #src.version root manifest #raw(src.data_root_manifest_sha256.slice(0, 16)) ·
    4L constants #raw(src.constants_sha256.at("4L", default: "").slice(0, 16)) ·
    5L values build #f5.build_utc (producer-update #raw(f5.producer_update_sha256.slice(0, 16))) ·
    generated #t.generated_utc · tables #raw(t.content_sha256.slice(0, 16))
  ])
  grid(columns: (1fr, auto), column-gutter: 10pt, align: horizon,
    text(size: 9pt)[Class of every Laurent coefficient, worst case X rounded down. *Precision table 1*: every master the data tree
      ships, at every order it holds. *Precision table 2*: the result of any single-scale vacuum input with $L$ loops (RustRed
      backend), worst case over the inputs; an order that needs master coefficients the data do not hold (spurious poles of
      the reduction) is not available.],
    pt-legend)
  v(0.2em)
  text(9pt, weight: "bold", fill: pt-navy)[Precision table 1. Master coefficients]
  text(7.4pt, fill: luma(80))[ — every master of the data tree (rows, by loop count; 5L: the 640-pack name in brackets) at every
    Laurent order held (columns); — = not held. #super[p] X rests in part on a constant computed with two profiles of one code
    beyond its independent confirmation; #super[n] Numerical mode returns more digits.]
  columns(2, gutter: 6mm, {
    pt-table1(t)
    set text(size: 6pt, font: pt-sans)
    set par(justify: false, leading: 0.4em)
    if t.two_profile_constants.len() > 0 [
      #super[p] Two-profile constants (digits; independent confirmation): #t.two_profile_constants.join("; ").
    ]
    if t.numerical_mode_more_digits.len() > 0 [
      #super[n] Numerical mode uses the best payload of a master's coefficient group (same exact basis expression), here a
      validated numerical value: #t.numerical_mode_more_digits.map(x => [#x.at(0) S\*(#x.at(1)) → #x.at(2) d]).join("; ").
    ]
  })
  v(0.6em)
  text(9pt, weight: "bold", fill: pt-navy)[Precision table 2. The $L$-loop result]
  text(7.4pt, fill: luma(80))[ — class and digits reachable for every admitted input, per Laurent order of the result (letters:
    notes below).]
  v(-0.3em)
  pt-table2(t)
  block(width: 100%, text(size: 6.9pt, font: pt-sans)[
    #set par(justify: false, leading: 0.42em, spacing: 0.45em)
    #for (i, n) in t.notes_table2.enumerate() [
      #text(weight: "bold")[#str.from-unicode(97 + i)] #n \
    ]
  ])
  block(width: 100%, above: 0.8em, text(size: 7pt)[
    #text(9pt, weight: "bold", fill: pt-navy)[How the precision tables are made]
    #set par(justify: true, leading: 0.46em, spacing: 0.48em)
    #for n in t.method [- #n]
  ])
})
